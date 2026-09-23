#!/usr/bin/env python3
"""Run and compare frozen homodimer reconciliation outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence


IGNORED_COMPONENTS = frozenset({"input", "logs", ".cache"})
IGNORED_FILENAMES = frozenset(
    {"pipeline_results.json", "report.html", ".DS_Store", "Thumbs.db"}
)
TOOL = "ColabFold v1.6.0 / AlphaFold-Multimer"
ORCHESTRATION_SOURCE_COMMIT = "08e0b5f"
PARITY_POLICY = {
    "name": "ipsae-operational-columns",
    "version": 1,
    "path": "ipsae/ipsae_summary.csv",
    "canonicalized_columns": ["pdb_path", "processing_time_ms"],
}


@dataclass(frozen=True)
class ArtifactDigest:
    relative_path: str
    size: int
    sha256: str


@dataclass(frozen=True)
class ParityResult:
    identical: bool
    added: tuple[str, ...]
    removed: tuple[str, ...]
    changed: tuple[str, ...]
    raw_changed: tuple[str, ...] = ()
    accepted_operational_drift: tuple[str, ...] = ()
    policy_errors: tuple[str, ...] = ()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_retained(relative_path: Path) -> bool:
    return (
        not IGNORED_COMPONENTS.intersection(relative_path.parts)
        and relative_path.name not in IGNORED_FILENAMES
    )


def inventory_output(output_dir: Path) -> tuple[ArtifactDigest, ...]:
    """Return exact digests for every retained regular file below output_dir."""
    root = Path(output_dir)
    if not root.exists():
        return ()
    artifacts = []
    for path in root.rglob("*"):
        relative_path = path.relative_to(root)
        if path.is_file() and _is_retained(relative_path):
            artifacts.append(
                ArtifactDigest(
                    relative_path=relative_path.as_posix(),
                    size=path.stat().st_size,
                    sha256=_sha256(path),
                )
            )
    return tuple(sorted(artifacts, key=lambda item: item.relative_path))


def inventory_fixture(fixture_dir: Path) -> tuple[ArtifactDigest, ...]:
    """Return exact digests for every regular file in the frozen fixture."""
    root = Path(fixture_dir)
    artifacts = []
    for path in root.rglob("*"):
        if path.is_file():
            artifacts.append(
                ArtifactDigest(
                    relative_path=path.relative_to(root).as_posix(),
                    size=path.stat().st_size,
                    sha256=_sha256(path),
                )
            )
    return tuple(sorted(artifacts, key=lambda item: item.relative_path))


def counts_by_top_level(
    artifacts: Sequence[ArtifactDigest],
) -> dict[str, int]:
    counts = Counter(item.relative_path.split("/", 1)[0] for item in artifacts)
    return dict(sorted(counts.items()))


def compare_outputs(
    baseline: Sequence[ArtifactDigest],
    candidate: Sequence[ArtifactDigest],
    baseline_root: Path | None = None,
    candidate_root: Path | None = None,
) -> ParityResult:
    baseline_by_path = {item.relative_path: item for item in baseline}
    candidate_by_path = {item.relative_path: item for item in candidate}
    baseline_paths = set(baseline_by_path)
    candidate_paths = set(candidate_by_path)
    added = tuple(sorted(candidate_paths - baseline_paths))
    removed = tuple(sorted(baseline_paths - candidate_paths))
    raw_changed = tuple(
        sorted(
            path
            for path in baseline_paths & candidate_paths
            if baseline_by_path[path] != candidate_by_path[path]
        )
    )
    changed = set(raw_changed)
    accepted = []
    policy_errors = []
    policy_path = str(PARITY_POLICY["path"])
    if (
        policy_path in baseline_paths & candidate_paths
        and baseline_root is not None
        and candidate_root is not None
    ):
        baseline_csv, baseline_error = _canonical_ipsae_csv(
            Path(baseline_root) / policy_path, "baseline"
        )
        candidate_csv, candidate_error = _canonical_ipsae_csv(
            Path(candidate_root) / policy_path, "candidate"
        )
        policy_errors.extend(error for error in (baseline_error, candidate_error) if error)
        if not policy_errors:
            if baseline_csv[0] != candidate_csv[0]:
                policy_errors.append(
                    f"{PARITY_POLICY['path']}: header bytes differ between outputs"
                )
                changed.add(policy_path)
            elif baseline_csv == candidate_csv:
                if policy_path in changed:
                    changed.remove(policy_path)
                    accepted.append(policy_path)
            else:
                changed.add(policy_path)
    return ParityResult(
        identical=not (added or removed or changed or policy_errors),
        added=added,
        removed=removed,
        changed=tuple(sorted(changed)),
        raw_changed=raw_changed,
        accepted_operational_drift=tuple(accepted),
        policy_errors=tuple(policy_errors),
    )


def _canonical_ipsae_csv(
    path: Path, label: str
) -> tuple[tuple[bytes, ...] | None, str | None]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        return None, f"{label} {PARITY_POLICY['path']}: malformed CSV: {exc}"
    if not raw:
        return None, f"{label} {PARITY_POLICY['path']}: malformed CSV: empty file"

    lines = raw.splitlines(keepends=True)
    if b"".join(lines) != raw:
        return None, f"{label} {PARITY_POLICY['path']}: malformed CSV: line split"
    for line_number, line in enumerate(lines, start=1):
        if line.endswith(b"\r") or (
            b"\r" in line and not line.endswith(b"\r\n")
        ):
            return (
                None,
                f"{label} {PARITY_POLICY['path']}: malformed CSV: "
                f"unsupported line terminator on line {line_number}",
            )

    header_line = lines[0]
    header_body, _ = _separate_line_ending(header_line)
    if b'"' in header_body:
        if header_body.count(b'"') % 2:
            return None, f"{label} {PARITY_POLICY['path']}: malformed CSV header"
        return (
            None,
            f"{label} {PARITY_POLICY['path']}: ambiguous header dialect",
        )
    header = header_body.split(b",")
    indexes = {}
    for column in PARITY_POLICY["canonicalized_columns"]:
        column_bytes = str(column).encode("ascii")
        count = header.count(column_bytes)
        if count == 0:
            return (
                None,
                f"{label} {PARITY_POLICY['path']}: missing canonicalized column {column}",
            )
        if count != 1:
            return (
                None,
                f"{label} {PARITY_POLICY['path']}: duplicate canonicalized column {column}",
            )
        indexes[column_bytes] = header.index(column_bytes)
    canonical_rows = [header_line]
    for row_number, line in enumerate(lines[1:], start=2):
        body, ending = _separate_line_ending(line)
        if b'"' in body:
            if body.count(b'"') % 2:
                return (
                    None,
                    f"{label} {PARITY_POLICY['path']}: malformed CSV: "
                    f"unbalanced quote on row {row_number}",
                )
            return (
                None,
                f"{label} {PARITY_POLICY['path']}: ambiguous data-row dialect "
                f"on row {row_number}",
            )
        row = body.split(b",")
        if len(row) != len(header):
            return (
                None,
                f"{label} {PARITY_POLICY['path']}: malformed CSV: "
                f"row {row_number} has {len(row)} cells, expected {len(header)}",
            )
        canonical = list(row)
        for column, index in indexes.items():
            canonical[index] = b"<canonical:" + column + b">"
        canonical_rows.append(b",".join(canonical) + ending)
    return tuple(canonical_rows), None


def _separate_line_ending(line: bytes) -> tuple[bytes, bytes]:
    if line.endswith(b"\r\n"):
        return line[:-2], b"\r\n"
    if line.endswith(b"\n"):
        return line[:-1], b"\n"
    return line, b""


def pipeline_command(repo: Path, fixture: Path, output_dir: Path) -> list[str]:
    config = fixture / "config"
    return [
        sys.executable,
        str(repo / "scripts" / "production_pipeline.py"),
        "--repo-dir",
        str(repo),
        "--input-dir",
        str(fixture / "input"),
        "--output-dir",
        str(output_dir),
        "--mapping-file",
        str(config / "model_ids.txt"),
        "--chain-mapping",
        str(config / "chain_manifest.csv"),
        "--dataset-config",
        str(config / "dataset_config.json"),
        "--provider-json",
        str(config / "provider.json"),
        "--uniprot-db",
        str(config / "uniprot_example_subset.duckdb"),
        "--tool-used",
        TOOL,
        "--homodimer-tool-used",
        TOOL,
        "--workers",
        "1",
        "--batch-size",
        "1",
        "--analysis-batch-size",
        "1",
        "--dssp-algorithm",
        "pydssp",
        "--clash-device",
        "cpu",
        "--python-cmd",
        sys.executable,
        "--no-cache",
    ]


def run_pipeline(repo: Path, fixture: Path, output_dir: Path) -> None:
    """Run production_pipeline.py from repo against a fresh output tree."""
    output_dir = Path(output_dir)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    subprocess.run(
        pipeline_command(Path(repo), Path(fixture), output_dir),
        cwd=repo,
        check=True,
    )


def source_commit(repo: Path) -> str:
    marker = repo / ".source_commit"
    if marker.is_file():
        return marker.read_text(encoding="utf-8").strip()
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def runtime_image_record(path: Path, expected_sha256: str) -> dict[str, object]:
    valid_expected = re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha256) is not None
    readable = path.is_file() and os.access(path, os.R_OK)
    actual_sha256 = _sha256(path) if readable else None
    record: dict[str, object] = {
        "path": str(path),
        "expected_sha256": expected_sha256,
        "actual_sha256": actual_sha256,
        "readable_regular_file": readable,
        "size": path.stat().st_size if readable else None,
        "match": bool(
            valid_expected
            and actual_sha256 is not None
            and actual_sha256.lower() == expected_sha256.lower()
        ),
    }
    return record


def runtime_image_integrity_failure(record: dict[str, object]) -> str | None:
    if re.fullmatch(r"[0-9a-fA-F]{64}", str(record["expected_sha256"])) is None:
        return "invalid expected SHA-256 format"
    if not record["readable_regular_file"]:
        return "runtime image is not a readable regular file"
    if not record["match"]:
        return "runtime image SHA-256 mismatch"
    return None


def _write_report(path: Path, report: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-repo", type=Path, required=True)
    parser.add_argument("--candidate-repo", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--report-json", type=Path, required=True)
    parser.add_argument("--runtime-image-path", type=Path, required=True)
    parser.add_argument("--runtime-image-sha256", required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    baseline_output = args.work_dir / "baseline"
    candidate_output = args.work_dir / "candidate"
    image_record = runtime_image_record(
        args.runtime_image_path, args.runtime_image_sha256
    )
    report: dict[str, object] = {
        "baseline_commit": source_commit(args.baseline_repo),
        "candidate_commit": source_commit(args.candidate_repo),
        "orchestration_source_commit": ORCHESTRATION_SOURCE_COMMIT,
        "comparison_policy": PARITY_POLICY,
        "runtime_image": image_record,
        "fixture": {
            "path": str(args.fixture),
            "inventory": [asdict(item) for item in inventory_fixture(args.fixture)],
        },
        "runs": [],
    }
    integrity_failure = runtime_image_integrity_failure(image_record)
    if integrity_failure is not None:
        report["integrity_failure"] = integrity_failure
        _write_report(args.report_json, report)
        return 1

    exit_code = 0
    try:
        for label, repo, output in (
            ("baseline", args.baseline_repo, baseline_output),
            ("candidate", args.candidate_repo, candidate_output),
        ):
            command = pipeline_command(repo, args.fixture, output)
            run_pipeline(repo, args.fixture, output)
            artifacts = inventory_output(output)
            report["runs"].append(
                {
                    "label": label,
                    "command": command,
                    "inventory": [asdict(item) for item in artifacts],
                    "counts_by_top_level": counts_by_top_level(artifacts),
                }
            )
        baseline_inventory = inventory_output(baseline_output)
        candidate_inventory = inventory_output(candidate_output)
        result = compare_outputs(
            baseline_inventory,
            candidate_inventory,
            baseline_output,
            candidate_output,
        )
        report["parity"] = asdict(result)
        report["raw_differences"] = {
            "added": list(result.added),
            "removed": list(result.removed),
            "changed": list(result.raw_changed),
        }
        report["canonical_differences"] = {
            "added": list(result.added),
            "removed": list(result.removed),
            "changed": list(result.changed),
            "policy_errors": list(result.policy_errors),
        }
        baseline_by_path = {
            item.relative_path: item for item in baseline_inventory
        }
        candidate_by_path = {
            item.relative_path: item for item in candidate_inventory
        }
        report["accepted_operational_drift"] = [
            {
                "path": path,
                "policy": PARITY_POLICY,
                "baseline_raw": asdict(baseline_by_path[path]),
                "candidate_raw": asdict(candidate_by_path[path]),
            }
            for path in result.accepted_operational_drift
        ]
        if not result.identical:
            exit_code = 1
    except subprocess.CalledProcessError as exc:
        report["pipeline_failure"] = {
            "returncode": exc.returncode,
            "command": [str(value) for value in exc.cmd],
        }
        exit_code = 1
    finally:
        _write_report(args.report_json, report)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
