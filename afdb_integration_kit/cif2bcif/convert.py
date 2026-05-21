"""
CIF to BinaryCIF conversion.

Prefers Biotite (when version is supported) for in-process table-copy conversion;
falls back to Mol* cif2bcif CLI for Mol*/gemmi compatibility when Biotite is
unavailable, too old, or conversion fails.
Uses ProcessPoolExecutor for CPU-bound parallel execution.
"""
import gzip
import logging
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Tuple

import numpy as np

MOLSTAR_CIF2BCIF_CMD = "cif2bcif"

# BinaryCIF mask values: 0=present, 1="." (inapplicable), 2="?" (missing)
_MASK_PRESENT = np.uint8(0)
_MASK_DOT = np.uint8(1)
_MASK_QMARK = np.uint8(2)
_MISSING_TOKENS = frozenset((".", "?"))


def _build_bcif_column(raw_strings, binary_cif_column: type[Any]):
    """Build a BinaryCIFColumn with proper mask and type detection.

    Detects missing value tokens ('.' and '?'), builds a BinaryCIF mask,
    and promotes data to the tightest numeric type (int32 -> float64 -> string).
    """
    values = list(raw_strings)
    if not values:
        return binary_cif_column(np.array([], dtype="U1"))

    mask = np.array(
        [_MASK_DOT if v == "." else (_MASK_QMARK if v == "?" else _MASK_PRESENT)
         for v in values],
        dtype=np.uint8,
    )
    has_missing = bool(np.any(mask != _MASK_PRESENT))
    present_idx = [i for i, v in enumerate(values) if v not in _MISSING_TOKENS]

    # Try int32
    try:
        data = np.zeros(len(values), dtype=np.int32)
        for i in present_idx:
            data[i] = int(values[i])
        return binary_cif_column(data, mask if has_missing else None)
    except (ValueError, OverflowError):
        pass

    # Try float64
    try:
        data = np.zeros(len(values), dtype=np.float64)
        for i in present_idx:
            data[i] = float(values[i])
        return binary_cif_column(data, mask if has_missing else None)
    except (ValueError, OverflowError):
        pass

    # String fallback
    data = np.array(values, dtype="U")
    return binary_cif_column(data, mask if has_missing else None)


BIOTITE_MIN_VERSION = (0, 40, 0)

_biotite_version_ok_cache: bool | None = None


def _biotite_version_ok() -> bool:
    """Return True if Biotite is installed at a supported version."""
    global _biotite_version_ok_cache
    if _biotite_version_ok_cache is not None:
        return _biotite_version_ok_cache
    try:
        import biotite

        raw = getattr(biotite, "__version__", "0")
        parts = raw.split(".")[:3]
        ver = []
        for p in parts:
            digits = "".join(c for c in str(p) if c.isdigit()) or "0"
            ver.append(int(digits))
        while len(ver) < 3:
            ver.append(0)
        _biotite_version_ok_cache = tuple(ver) >= BIOTITE_MIN_VERSION
    except Exception:
        _biotite_version_ok_cache = False
    return _biotite_version_ok_cache


logger = logging.getLogger("cif2bcif")
logger.setLevel(logging.INFO)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def _process_single_cif(args: Tuple[Path, Path, str]) -> Tuple[str, bool]:
    """
    Process a single CIF file. Module-level function for ProcessPoolExecutor pickling.

    Args:
        args: Tuple of (input_file, output_dir, extension)

    Returns:
        Tuple of (filename, success)
    """
    input_file, output_dir, ext = args
    output_file = output_dir / (input_file.stem + ext)
    success = run_cif2bcif(input_file, output_file)
    return input_file.name, success


def _run_molstar_cif2bcif(input_file: Path, output_file: Path) -> bool:
    """
    Convert CIF to BCIF using Mol* (cif2bcif). Tries PATH first, then npx.
    Produces BCIF compatible with Mol* and gemmi. Gemmi cannot write BCIF.
    """
    for cmd in (_molstar_cmd_path(), _molstar_cmd_npx()):
        if not cmd:
            continue
        result = subprocess.run(
            cmd + [str(input_file), str(output_file)],
            capture_output=True,
            text=True,
            timeout=900,
        )
        if result.returncode == 0:
            return True
        logger.debug(
            f"Mol* cif2bcif failed: {result.stderr or result.stdout}"
        )
    return False


def _molstar_cmd_path() -> list[str] | None:
    """Command list if cif2bcif is on PATH."""
    cmd = shutil.which(MOLSTAR_CIF2BCIF_CMD)
    return [cmd] if cmd else None


def _molstar_cmd_npx() -> list[str] | None:
    """Command list for npx -p molstar cif2bcif (no global install)."""
    if not shutil.which("npx"):
        return None
    return ["npx", "--yes", "-p", "molstar", MOLSTAR_CIF2BCIF_CMD]


def _run_biotite_cif2bcif(
    input_file: Path, output_file: Path, tmpdir: str | None = None
) -> bool:
    """
    Convert CIF to BinaryCIF using Biotite (table-copy, no AtomArray).
    Preserves all categories and metadata via a temporary file.
    """
    output_str = str(output_file)
    base_tmpdir = tmpdir or os.environ.get("TMPDIR") or tempfile.gettempdir()
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        from biotite.structure.io.pdbx import CIFFile, BinaryCIFFile, BinaryCIFBlock
        from biotite.structure.io.pdbx.bcif import BinaryCIFCategory, BinaryCIFColumn

        cif = CIFFile.read(str(input_file))

        bcif = BinaryCIFFile()
        for block_name in cif.keys():
            cif_block = cif[block_name]
            bcif_block = BinaryCIFBlock()
            for cat_name in cif_block.keys():
                cif_cat = cif_block[cat_name]
                columns = {}
                for col_name in cif_cat.keys():
                    col = cif_cat[col_name]
                    # Get raw strings so we can detect '.' and '?' for masks
                    raw = col.as_array(str)
                    columns[col_name] = _build_bcif_column(raw, BinaryCIFColumn)
                bcif_cat = BinaryCIFCategory(columns)
                bcif_block[cat_name] = bcif_cat
            bcif[block_name] = bcif_block

        if output_str.endswith(".bcif.gz"):
            tmp_path = Path(base_tmpdir) / (output_file.name + ".tmp.bcif")
            bcif.write(str(tmp_path))
            try:
                with open(tmp_path, "rb") as f_in:
                    with gzip.open(output_file, "wb") as f_out:
                        f_out.write(f_in.read())
            finally:
                if tmp_path.exists():
                    tmp_path.unlink(missing_ok=True)
        else:
            tmp_path = Path(base_tmpdir) / (output_file.name + ".tmp")
            bcif.write(str(tmp_path))
            try:
                tmp_path.replace(output_file)
            except OSError:
                # Cross-device rename (e.g. /tmp -> /lustre); fall back to copy+delete
                shutil.move(str(tmp_path), str(output_file))
        return True

    except Exception as e:
        logger.debug(f"Biotite conversion failed for {input_file}: {e}")
        import traceback
        traceback.print_exc()
        return False


def run_cif2bcif(
    input_file: Path, output_file: Path, tmpdir: str | None = None
) -> bool:
    """
    Convert CIF to BinaryCIF. Prefers Biotite (when version supported) for in-process
    table-copy; falls back to Mol* cif2bcif CLI for Mol*/gemmi compatibility.

    Args:
        input_file: Input CIF file path
        output_file: Output BCIF file path (.bcif or .bcif.gz)
        tmpdir: Optional directory for temporary writes. Defaults to TMPDIR or
                the system temp directory.

    Returns:
        True on success, False on failure
    """
    output_str = str(output_file)
    if not (output_str.endswith(".bcif") or output_str.endswith(".bcif.gz")):
        logger.warning(
            f"Output file extension '{output_file.suffix}' is not '.bcif' or '.bcif.gz'"
        )

    logger.info(f"Converting {input_file} to {output_file}")

    if _biotite_version_ok():
        if _run_biotite_cif2bcif(input_file, output_file, tmpdir=tmpdir):
            logger.info(f"Conversion complete (Biotite): {output_file}")
            return True
        logger.warning("Biotite conversion failed; trying Mol* cif2bcif.")
    else:
        logger.warning(
            "Biotite missing or version < %s; using Mol* cif2bcif. "
            "Install/upgrade biotite for faster in-process conversion.",
            ".".join(str(x) for x in BIOTITE_MIN_VERSION),
        )
    if _run_molstar_cif2bcif(input_file, output_file):
        logger.info(f"Conversion complete (Mol*): {output_file}")
        return True
    logger.error(f"Conversion failed for {input_file} (Biotite and Mol* both failed)")
    return False


def run_batch_cif2bcif(
    input_dir: Path,
    output_dir: Path,
    workers: int = 8,
    gzip: bool = False,
    pattern: str = "*.cif"
) -> Tuple[int, int]:
    """
    Batch CIF to BinaryCIF conversion. Uses Biotite first (when version supported),
    then Mol* CLI on failure; ProcessPoolExecutor for parallel execution.

    Args:
        input_dir: Directory containing input CIF files
        output_dir: Output directory for BCIF files
        workers: Number of parallel workers
        gzip: Whether to gzip output files (.bcif.gz)
        pattern: Glob pattern for input files

    Returns:
        Tuple of (success_count, error_count)
    """
    input_files = list(input_dir.glob(pattern))
    if not input_files:
        logger.warning(f"No files matching '{pattern}' in {input_dir}")
        return 0, 0

    output_dir.mkdir(parents=True, exist_ok=True)

    ext = ".bcif.gz" if gzip else ".bcif"

    # Prepare arguments for each file (module-level function for pickling)
    work_items = [(f, output_dir, ext) for f in input_files]

    logger.info(
        f"Processing {len(input_files)} files with {workers} workers (ProcessPool)"
    )

    success_count = 0
    error_count = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = executor.map(_process_single_cif, work_items)
        for fname, ok in results:
            if ok:
                success_count += 1
                logger.info(f"[OK] {fname}")
            else:
                error_count += 1
                logger.error(f"[FAILED] {fname}")

    logger.info(
        f"Batch conversion complete: {success_count} success, {error_count} errors"
    )
    return success_count, error_count


# Legacy function for backwards compatibility with molstar cif2bcif subprocess
def process_file(input_path, output_path):
    """Legacy subprocess wrapper - kept for compatibility but not recommended."""
    import subprocess
    command = ["cif2bcif", str(input_path), str(output_path)]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode, input_path, result.stdout, result.stderr
