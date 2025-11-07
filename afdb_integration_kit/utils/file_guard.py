from __future__ import annotations

import logging
from pathlib import Path

import typer

logger = logging.getLogger("afdb_integration_kit")


def require_non_empty_file(path: Path, *, description: str | None = None) -> None:
    """
    Ensure the given path exists, is a file, and is not zero bytes.
    """
    label = description or str(path)
    if not path.exists():
        logger.error("%s not found: %s", label, path)
        raise typer.Exit(code=1)
    if not path.is_file():
        logger.error("%s is not a file: %s", label, path)
        raise typer.Exit(code=1)
    try:
        size = path.stat().st_size
    except OSError as exc:  # pragma: no cover - filesystem specific
        logger.error("Unable to stat %s (%s): %s", label, path, exc)
        raise typer.Exit(code=1)
    if size == 0:
        logger.error("%s is empty: %s", label, path)
        raise typer.Exit(code=1)
