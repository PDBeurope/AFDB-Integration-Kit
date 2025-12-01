"""ColabFold support utilities for converting outputs into AFDB ingestion format."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .converter import convert_file, main

__all__ = ["convert_file", "main"]


def __getattr__(name: str):
    if name in __all__:
        from .converter import convert_file, main  # lazy import to avoid runpy warning
        return {"convert_file": convert_file, "main": main}[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
