from __future__ import annotations

import json
import logging
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Optional, Sequence

from .context import ValidationContext
from .registry import REGISTERED_CHECKS
from .results import Level, ValidationResult


def _discover_files(root: Path) -> List[Path]:
    files = [p for p in root.rglob("*") if p.is_file()]
    files.sort(key=lambda p: str(p.relative_to(root)))
    return files


def _files_for_check(check: str, files: Sequence[Path]) -> List[Path]:
    lower = check.lower()
    if lower in {"naming", "relationships", "sequences", "metadata"}:
        return list(files)
    if lower in {"plddt", "pae"}:
        return [p for p in files if p.suffix.lower() == ".json"]
    if lower == "pdb":
        return [p for p in files if p.suffix.lower() == ".pdb"]
    return []


def _sorted_unique(items: Iterable[str]) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            ordered.append(item)
    return ordered


def _resolve_checks(requested: Optional[Sequence[str]]) -> List[str]:
    if not requested or requested == ["all"]:
        return sorted(REGISTERED_CHECKS)
    resolved = _sorted_unique(requested)
    unknown = [name for name in resolved if name not in REGISTERED_CHECKS]
    if unknown:
        available = ", ".join(sorted(REGISTERED_CHECKS))
        raise ValueError(f"Unknown validation checks: {', '.join(unknown)}. Available: {available}")
    return resolved


def _run_single_check(check: str, files: List[Path], ctx: ValidationContext) -> List[ValidationResult]:
    try:
        return REGISTERED_CHECKS[check](files, ctx)
    except Exception as exc:  # pragma: no cover - defensive
        ctx.logger.exception("Validation check '%s' failed", check)
        return [
            ValidationResult(
                check=check,
                file=ctx.root,
                level=Level.ERROR,
                code="internal_error",
                message=f"{type(exc).__name__}: {exc}",
                suggested_fix="Inspect logs for stack trace and fix validator implementation.",
            )
        ]


def run_validations(
    root: Path,
    checks: Optional[Sequence[str]] = None,
    config: Optional[Dict] = None,
    *,
    parallel: bool = True,
) -> List[ValidationResult]:
    if not root.exists():
        raise FileNotFoundError(root)
    root = root.resolve()
    config = dict(config or {})

    logger = logging.getLogger("afdb_integration_kit.validation")
    selected = _resolve_checks(checks)
    ctx = ValidationContext(root=root, config=config, selected_checks=set(selected), logger=logger)

    files = _discover_files(root)
    files_by_check = {name: _files_for_check(name, files) for name in selected}

    results: Dict[str, List[ValidationResult]] = {name: [] for name in selected}

    if parallel and len(selected) > 1:
        with ThreadPoolExecutor(max_workers=min(len(selected), 8)) as executor:
            futures = {
                executor.submit(_run_single_check, name, files_by_check[name], ctx): name for name in selected
            }
            for future in as_completed(futures):
                name = futures[future]
                results[name] = future.result()
    else:
        for name in selected:
            results[name] = _run_single_check(name, files_by_check[name], ctx)

    ordered: List[ValidationResult] = []
    for name in selected:
        ordered.extend(sorted(results[name], key=lambda res: (res.file, res.code)))
    return ordered


def run_validation_check(
    check: str,
    files: Sequence[Path],
    config: Optional[Dict] = None,
) -> List[ValidationResult]:
    if check not in REGISTERED_CHECKS:
        available = ", ".join(sorted(REGISTERED_CHECKS)) or "<none>"
        raise ValueError(f"Unknown validation check '{check}'. Available: {available}")
    resolved_files = [Path(path).resolve() for path in files]
    if not resolved_files:
        raise ValueError("files must not be empty")
    logger = logging.getLogger("afdb_integration_kit.validation")
    ctx = ValidationContext(
        root=resolved_files[0].parent,
        config=dict(config or {}),
        selected_checks={check},
        logger=logger,
    )
    return _run_single_check(check, resolved_files, ctx)


def summarise(results: Sequence[ValidationResult]) -> Dict[str, Dict[str, int]]:
    counts_by_level: Dict[str, int] = defaultdict(int)
    counts_by_check: Dict[str, Dict[str, int]] = {}

    for res in results:
        counts_by_level[res.level.value] += 1
        bucket = counts_by_check.setdefault(res.check, defaultdict(int))  # type: ignore[assignment]
        bucket[res.level.value] += 1  # type: ignore[index]

    # materialise nested defaultdicts
    checks_summary = {
        check: {level: count for level, count in sorted(level_counts.items())}
        for check, level_counts in counts_by_check.items()
    }

    return {
        "levels": {level: counts_by_level[level] for level in sorted(counts_by_level)},
        "checks": checks_summary,
    }


def write_results(results: Sequence[ValidationResult], path: Path, fmt: Literal["json", "jsonl"]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "json":
        payload = [res.to_dict() for res in results]
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        return
    if fmt == "jsonl":
        with path.open("w", encoding="utf-8") as handle:
            for res in results:
                handle.write(json.dumps(res.to_dict()) + "\n")
        return
    raise ValueError(f"Unsupported format: {fmt}")


__all__ = ["run_validation_check", "run_validations", "summarise", "write_results"]
