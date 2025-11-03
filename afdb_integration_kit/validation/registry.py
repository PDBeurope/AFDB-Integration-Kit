from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List

from .context import ValidationContext
from .results import ValidationResult

CheckCallable = Callable[[list[Path], ValidationContext], list[ValidationResult]]

REGISTERED_CHECKS: Dict[str, CheckCallable] = {}


def register_check(name: str) -> Callable[[CheckCallable], CheckCallable]:
    """Decorator to register a validation callable under ``name``."""

    def decorator(func: CheckCallable) -> CheckCallable:
        if name in REGISTERED_CHECKS:
            raise ValueError(f"Validation check '{name}' already registered")
        REGISTERED_CHECKS[name] = func
        return func

    return decorator


def list_registered_checks() -> List[str]:
    return sorted(REGISTERED_CHECKS.keys())


__all__ = ["CheckCallable", "REGISTERED_CHECKS", "list_registered_checks", "register_check"]
