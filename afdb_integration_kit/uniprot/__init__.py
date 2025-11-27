"""
Programmatic UniProt utilities wrapping the streaming scripts.

These helpers let you shard releases, extract subsets, and merge shard outputs
without shelling out to the CLI entrypoints in ``uniprot/scripts``.
"""

from .sharding import shard_release
from .extraction import extract_subset_to_parquet
from .merge import merge_parquet_shards

__all__ = ["shard_release", "extract_subset_to_parquet", "merge_parquet_shards"]
