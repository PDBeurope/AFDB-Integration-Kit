"""Parallel processing utilities for validation."""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import orjson


def load_json_files_parallel(
    paths: List[Path],
    max_workers: Optional[int] = None,
) -> Dict[Path, Union[object, Exception]]:
    """Load multiple JSON files in parallel using threads (I/O bound).

    Args:
        paths: List of paths to JSON files.
        max_workers: Maximum number of threads to use (default: all available CPUs).

    Returns:
        Dictionary mapping paths to loaded JSON data or Exception if loading failed.
    """
    if not paths:
        return {}

    if max_workers is None:
        max_workers = os.cpu_count() or 8

    def load_one(p: Path) -> Tuple[Path, Union[object, Exception]]:
        try:
            return p, orjson.loads(p.read_bytes())
        except Exception as exc:
            return p, exc

    effective_workers = min(len(paths), max_workers)

    if effective_workers <= 1:
        # Sequential for single file
        return dict(load_one(p) for p in paths)

    with ThreadPoolExecutor(max_workers=effective_workers) as executor:
        results = list(executor.map(load_one, paths))

    return dict(results)


def load_files_parallel(
    paths: List[Path],
    max_workers: Optional[int] = None,
) -> Dict[Path, Union[bytes, Exception]]:
    """Load multiple files as bytes in parallel using threads (I/O bound).

    Args:
        paths: List of paths to files.
        max_workers: Maximum number of threads to use (default: all available CPUs).

    Returns:
        Dictionary mapping paths to file bytes or Exception if loading failed.
    """
    if not paths:
        return {}

    if max_workers is None:
        max_workers = os.cpu_count() or 8

    def load_one(p: Path) -> Tuple[Path, Union[bytes, Exception]]:
        try:
            return p, p.read_bytes()
        except Exception as exc:
            return p, exc

    effective_workers = min(len(paths), max_workers)

    if effective_workers <= 1:
        return dict(load_one(p) for p in paths)

    with ThreadPoolExecutor(max_workers=effective_workers) as executor:
        results = list(executor.map(load_one, paths))

    return dict(results)


__all__ = ["load_json_files_parallel", "load_files_parallel"]
