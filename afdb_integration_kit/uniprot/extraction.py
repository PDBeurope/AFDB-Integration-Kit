"""
Programmatic subset extraction mirroring the CLI script.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional, Sequence

from uniprot.scripts.extract_subset import (
    ParquetBufferWriter,
    ensure_output_dir,
    get_entry_schema,
    load_targets,
    load_targets_from_mapping,
    process_inputs,
    stable_shard_for_accession,
)

logger = logging.getLogger(__name__)


def _load_targets(target_files: Iterable[Path], mapping_files: Iterable[Path]) -> set[str]:
    targets: set[str] = set()
    for path in target_files:
        targets |= load_targets(path)
    for path in mapping_files:
        targets |= load_targets_from_mapping(path)
    return targets


def extract_subset_to_parquet(
    inputs: Sequence[Path],
    outdir: Path,
    release: str,
    targets: Optional[Sequence[Path]] = None,
    mappings: Optional[Sequence[Path]] = None,
    shard_count: Optional[int] = None,
    shard_index: Optional[int] = None,
    batch_size: int = 1000,
) -> None:
    """
    Stream UniProt flat files and write the matching subset to a Parquet file.

    This mirrors the behaviour of ``uniprot/scripts/extract_subset.py`` but as a
    callable function.
    """
    target_files = targets or []
    mapping_files = mappings or []
    target_set = _load_targets(target_files, mapping_files)
    if not target_set:
        raise ValueError("No accessions provided (targets or mappings required).")

    if (shard_count is None) ^ (shard_index is None):
        raise ValueError("Provide both shard_count and shard_index when filtering by shard.")
    if shard_count is not None:
        if shard_index is None or shard_index < 0 or shard_index >= shard_count:
            raise ValueError("shard_index must be between 0 and shard_count-1.")
        target_set = {
            ac for ac in target_set if stable_shard_for_accession(ac, shard_count) == shard_index
        }
        logger.info(
            "Filtered targets to shard %d of %d (%d remaining).",
            shard_index,
            shard_count,
            len(target_set),
        )
        if not target_set:
            ensure_output_dir(outdir)
            ParquetBufferWriter(outdir / "entry.parquet", get_entry_schema(), batch_size).close()
            logger.info("No targets in shard %d; wrote empty parquet.", shard_index)
            return

    ensure_output_dir(outdir)
    writer = ParquetBufferWriter(outdir / "entry.parquet", get_entry_schema(), batch_size)
    try:
        process_inputs(inputs, target_set, release, writer)
    finally:
        writer.close()


__all__ = ["extract_subset_to_parquet"]
