"""
Merge shard-level Parquet files into a single Parquet output.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


def merge_parquet_shards(inputs: Iterable[str], output: Path, compression: str = "zstd") -> None:
    dataset = ds.dataset(list(inputs), format="parquet")
    writer: pq.ParquetWriter | None = None
    try:
        for batch in dataset.to_batches():
            if writer is None:
                writer = pq.ParquetWriter(output, dataset.schema, compression=compression)
            writer.write_table(pa.Table.from_batches([batch]))
    finally:
        if writer is not None:
            writer.close()


__all__ = ["merge_parquet_shards"]
