from __future__ import annotations

from functools import lru_cache
from importlib import import_module
from typing import Iterable

_DEFAULT_MODULES: tuple[str, ...] = (
    "afdb_integration_kit.quality_assessment.naming",
    "afdb_integration_kit.quality_assessment.pLDDT",
)


@lru_cache(maxsize=None)
def ensure_default_validators(modules: Iterable[str] | None = None) -> None:
    """Import built-in validator modules so they self-register."""
    to_import = modules if modules is not None else _DEFAULT_MODULES
    for dotted in to_import:
        import_module(dotted)


__all__ = ["ensure_default_validators"]
