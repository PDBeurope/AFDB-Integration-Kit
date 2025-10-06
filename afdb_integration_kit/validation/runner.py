from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from .defaults import ensure_default_validators
from .registry import (
    FormatCallable,
    ValidationHook,
    get_validator,
    list_validators,
)


@dataclass
class ValidationResult:
    name: str
    ok: bool
    report: Dict[str, Any]
    formatter: Optional[FormatCallable]
    description: Optional[str]
    options: Dict[str, Any]


def _resolve_hooks(checks: Optional[Sequence[str]]) -> Sequence[ValidationHook]:
    if checks is None:
        return list(list_validators())
    return [get_validator(name) for name in checks]


def run_validations(
    root: Path | str,
    *,
    checks: Optional[Sequence[str]] = None,
    overrides: Optional[Mapping[str, Dict[str, Any]]] = None,
) -> Tuple[bool, Sequence[ValidationResult]]:
    """Run one or more registered validations for ``root``.

    Args:
        root: Dataset directory to validate.
        checks: Optional explicit sequence of validator names. Defaults to all registered.
        overrides: Optional mapping of validator name to keyword overrides passed into the
            underlying validation callable.

    Returns:
        (overall_ok, results)
    """
    if overrides is None:
        overrides = {}

    ensure_default_validators()
    path = Path(root).expanduser().resolve()
    hooks = _resolve_hooks(checks)

    results: list[ValidationResult] = []
    overall_ok = True

    for hook in hooks:
        kwargs = hook.build_kwargs(overrides.get(hook.name))
        ok, report = hook.run(path, **kwargs)
        results.append(
            ValidationResult(
                name=hook.name,
                ok=ok,
                report=report,
                formatter=hook.formatter,
                description=hook.description,
                options=kwargs,
            )
        )
        if not ok:
            overall_ok = False

    return overall_ok, results


__all__ = ["ValidationResult", "run_validations"]
