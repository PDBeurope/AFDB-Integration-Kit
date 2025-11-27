"""
Sharding helpers that wrap the streaming sharder script.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

from uniprot.scripts.shard_uniprot import shard_file, stable_shard_for_accession

logger = logging.getLogger(__name__)


def shard_release(
    inputs: Sequence[Path],
    outdir: Path,
    release: str,
    shard_count: int = 8,
    gzip_level: int = 5,
) -> None:
    """Shard one or more UniProt flat files by accession."""
    for path in inputs:
        shard_file(
            input_path=path,
            outdir=outdir,
            release=release,
            shard_count=shard_count,
            gzip_level=gzip_level,
        )


__all__ = ["shard_release", "stable_shard_for_accession"]
