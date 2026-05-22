"""
Secondary structure calculation with multiple algorithm support.

Supports four algorithms for secondary structure assignment:
1. mkdssp: External DSSP binary, preserving the historical default behavior
2. P-SEA (Biotite's annotate_sse): Geometry-based, ~95% agreement with DSSP
3. PyDSSP: Simplified H-bond based DSSP, ~97% agreement with DSSP
4. TM-align: CA-CA distance based (from Foldseek/TM-align), very fast

The Python algorithms produce 3-state output: helix, strand, coil (turn mapped
to coil for 3-state). They use gemmi for mmCIF I/O to produce complete
_struct_conf annotations with all standard fields (label/auth identifiers,
insertion codes).
"""
import logging
import subprocess
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple, Literal, Optional

import gemmi
import numpy as np

logger = logging.getLogger("dssp")
logger.setLevel(logging.INFO)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

Algorithm = Literal["mkdssp", "psea", "pydssp", "tmalign"]
DEFAULT_ALGORITHM: Algorithm = "mkdssp"

PSEA_CODE_MAP = {
    'a': 'HELX_P',
    'b': 'STRN',
    'c': None,
}

PYDSSP_CODE_MAP = {
    'H': 'HELX_P',
    'E': 'STRN',
    '-': None,
}

TMALIGN_CODE_MAP = {
    1: None,
    2: 'HELX_P',
    3: None,
    4: 'STRN',
}

ResidueInfo = namedtuple("ResidueInfo", [
    "label_asym_id",
    "auth_asym_id",
    "label_seq_id",
    "auth_seq_id",
    "comp_id",
    "ins_code",
])

_CURRENT_ALGORITHM: Algorithm = DEFAULT_ALGORITHM
_CURRENT_DEVICE: str = "cpu"

_BACKBONE_ATOMS = ("N", "CA", "C", "O")

_STRUCT_CONF_FIELDS = [
    "id", "conf_type_id",
    "beg_label_comp_id", "beg_label_asym_id", "beg_label_seq_id",
    "pdbx_beg_PDB_ins_code",
    "end_label_comp_id", "end_label_asym_id", "end_label_seq_id",
    "pdbx_end_PDB_ins_code",
    "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
    "end_auth_comp_id", "end_auth_asym_id", "end_auth_seq_id",
]


def _make_residue_info(chain: gemmi.Chain, residue: gemmi.Residue) -> ResidueInfo:
    icode = residue.seqid.icode
    has_icode = icode not in ('\0', ' ', '')
    label_seq = str(residue.label_seq) if residue.label_seq else str(residue.seqid.num)
    return ResidueInfo(
        label_asym_id=residue.subchain or chain.name,
        auth_asym_id=chain.name,
        label_seq_id=label_seq,
        auth_seq_id=str(residue.seqid.num),
        comp_id=residue.name,
        ins_code=str(icode) if has_icode else "?",
    )


def _extract_polymer_residue_info(structure: gemmi.Structure) -> List[ResidueInfo]:
    """Extract per-residue metadata for all polymer residues in the first model."""
    info: List[ResidueInfo] = []
    for chain in structure[0]:
        polymer = chain.get_polymer()
        if not polymer:
            continue
        for residue in polymer:
            info.append(_make_residue_info(chain, residue))
    return info


def _extract_backbone_coords(
    structure: gemmi.Structure,
) -> Tuple[Optional[np.ndarray], List[ResidueInfo]]:
    """
    Extract N, CA, C, O backbone coordinates and residue metadata.

    Returns (coords shaped (n_residues, 4, 3), residue_info) or (None, [])
    if any polymer residue is missing a backbone atom.
    """
    coords_list: list = []
    residue_info: List[ResidueInfo] = []

    for chain in structure[0]:
        polymer = chain.get_polymer()
        if not polymer:
            continue
        for residue in polymer:
            atom_coords = []
            for atom_name in _BACKBONE_ATOMS:
                atom = residue.find_atom(atom_name, '\0')
                if atom is None:
                    return None, []
                atom_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
            coords_list.append(atom_coords)
            residue_info.append(_make_residue_info(chain, residue))

    if not coords_list:
        return None, []
    return np.array(coords_list, dtype=np.float32), residue_info


def _extract_ca_coords(
    structure: gemmi.Structure,
) -> Tuple[Optional[np.ndarray], List[ResidueInfo]]:
    """
    Extract CA coordinates and residue metadata.

    Returns (coords shaped (n_residues, 3), residue_info) or (None, [])
    if any polymer residue is missing CA.
    """
    coords_list: list = []
    residue_info: List[ResidueInfo] = []

    for chain in structure[0]:
        polymer = chain.get_polymer()
        if not polymer:
            continue
        for residue in polymer:
            ca = residue.find_atom("CA", '\0')
            if ca is None:
                return None, []
            coords_list.append([ca.pos.x, ca.pos.y, ca.pos.z])
            residue_info.append(_make_residue_info(chain, residue))

    if not coords_list:
        return None, []
    return np.array(coords_list, dtype=np.float32), residue_info


def _tmalign_sec_str(
    d13: float, d14: float, d15: float, d24: float, d25: float, d35: float
) -> int:
    """
    Classify secondary structure based on CA-CA distances.
    Implements the sec_str function from TM-align.

    Returns: 1 = coil, 2 = helix, 3 = turn, 4 = strand
    """
    delta_helix = 2.1
    if (
        abs(d15 - 6.37) < delta_helix
        and abs(d14 - 5.18) < delta_helix
        and abs(d25 - 5.18) < delta_helix
        and abs(d13 - 5.45) < delta_helix
        and abs(d24 - 5.45) < delta_helix
        and abs(d35 - 5.45) < delta_helix
    ):
        return 2

    delta_strand = 1.42
    if (
        abs(d15 - 13.0) < delta_strand
        and abs(d14 - 10.4) < delta_strand
        and abs(d25 - 10.4) < delta_strand
        and abs(d13 - 6.1) < delta_strand
        and abs(d24 - 6.1) < delta_strand
        and abs(d35 - 6.1) < delta_strand
    ):
        return 4

    if d15 < 8:
        return 3
    return 1


def _compute_sse_tmalign(ca_coords: np.ndarray) -> np.ndarray:
    """Compute SSE using TM-align's make_sec algorithm from CA coordinates."""
    n_residues = len(ca_coords)
    sse = np.ones(n_residues, dtype=np.int32)

    for i in range(2, n_residues - 2):
        d13 = np.linalg.norm(ca_coords[i - 2] - ca_coords[i])
        d14 = np.linalg.norm(ca_coords[i - 2] - ca_coords[i + 1])
        d15 = np.linalg.norm(ca_coords[i - 2] - ca_coords[i + 2])
        d24 = np.linalg.norm(ca_coords[i - 1] - ca_coords[i + 1])
        d25 = np.linalg.norm(ca_coords[i - 1] - ca_coords[i + 2])
        d35 = np.linalg.norm(ca_coords[i] - ca_coords[i + 2])
        sse[i] = _tmalign_sec_str(d13, d14, d15, d24, d25, d35)

    return sse


def _compute_sse_pydssp(coords: np.ndarray, device: str = "cpu") -> np.ndarray:
    """Compute SSE using PyDSSP from pre-extracted backbone coordinates."""
    try:
        import pydssp
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "The 'pydssp' package is required for --algorithm pydssp. "
            "Install production dependencies with `uv pip install '.[production]'`."
        ) from exc

    if device != "cpu" and _torch_cuda_available():
        import torch

        coords_tensor = torch.tensor(coords, dtype=torch.float32, device=device)
        return pydssp.assign(coords_tensor, out_type='c3')
    return pydssp.assign(coords, out_type='c3')


def _compute_sse_psea(input_file: Path) -> Optional[np.ndarray]:
    """
    Compute SSE using P-SEA (Biotite). Requires a secondary Biotite read
    since annotate_sse operates on Biotite AtomArrays.

    Returns array of SSE codes ('a', 'b', 'c') per residue, or None on failure.
    """
    try:
        from biotite.structure import annotate_sse
        from biotite.structure.io.pdbx import CIFFile as BiotiteCIFFile, get_structure

        cif = BiotiteCIFFile.read(str(input_file))
        atom_array = get_structure(cif, model=1)
        return annotate_sse(atom_array)
    except Exception as e:
        logger.error(f"P-SEA computation failed: {e}")
        return None


def _make_range_dict(conf_type_id: str, beg: ResidueInfo, end: ResidueInfo) -> dict:
    return {
        'conf_type_id': conf_type_id,
        'beg_label_comp_id': beg.comp_id,
        'beg_label_asym_id': beg.label_asym_id,
        'beg_label_seq_id': beg.label_seq_id,
        'pdbx_beg_PDB_ins_code': beg.ins_code,
        'end_label_comp_id': end.comp_id,
        'end_label_asym_id': end.label_asym_id,
        'end_label_seq_id': end.label_seq_id,
        'pdbx_end_PDB_ins_code': end.ins_code,
        'beg_auth_comp_id': beg.comp_id,
        'beg_auth_asym_id': beg.auth_asym_id,
        'beg_auth_seq_id': beg.auth_seq_id,
        'end_auth_comp_id': end.comp_id,
        'end_auth_asym_id': end.auth_asym_id,
        'end_auth_seq_id': end.auth_seq_id,
    }


def _find_secondary_structure_ranges(
    sse_codes: np.ndarray,
    residue_info: List[ResidueInfo],
    code_map: dict,
) -> List[dict]:
    """
    Convert per-residue SSE codes into contiguous ranges for _struct_conf.

    Each range contains all fields needed for a complete _struct_conf entry
    including label/auth identifiers and insertion codes.
    """
    ranges: List[dict] = []
    if len(sse_codes) == 0:
        return ranges

    current_type = None
    start_idx = 0

    for i, code in enumerate(sse_codes):
        ss_type = code_map.get(code)

        if ss_type != current_type:
            if current_type is not None and i > start_idx:
                ranges.append(_make_range_dict(
                    current_type, residue_info[start_idx], residue_info[i - 1],
                ))
            current_type = ss_type
            start_idx = i

    if current_type is not None and len(sse_codes) > start_idx:
        ranges.append(_make_range_dict(
            current_type, residue_info[start_idx], residue_info[-1],
        ))

    return ranges


def _get_algorithm_label(algorithm: str) -> str:
    labels = {
        "mkdssp": "mkdssp",
        "psea": "P-SEA",
        "pydssp": "PyDSSP",
        "tmalign": "TM-align",
    }
    return labels.get(algorithm, algorithm)


def _torch_cuda_available() -> bool:
    """Return whether torch CUDA is available without making torch a core import."""
    try:
        import torch
    except ModuleNotFoundError:
        return False
    return bool(torch.cuda.is_available())


def _add_struct_conf_to_cif(
    block: gemmi.cif.Block,
    ranges: List[dict],
    algorithm: str,
) -> None:
    """Add _struct_conf and _struct_conf_type categories to a gemmi CIF block."""
    if not ranges:
        return

    criteria_label = _get_algorithm_label(algorithm)

    conf_types = sorted(set(r['conf_type_id'] for r in ranges))
    type_loop = block.init_loop("_struct_conf_type.", ["id", "criteria"])
    for ct in conf_types:
        type_loop.add_row([ct, criteria_label])

    conf_loop = block.init_loop("_struct_conf.", _STRUCT_CONF_FIELDS)
    type_counters: dict = {}
    for r in ranges:
        ct = r['conf_type_id']
        type_counters[ct] = type_counters.get(ct, 0) + 1
        conf_id = f"{ct}{type_counters[ct]}"

        conf_loop.add_row([
            conf_id,
            r['conf_type_id'],
            r['beg_label_comp_id'],
            r['beg_label_asym_id'],
            r['beg_label_seq_id'],
            r['pdbx_beg_PDB_ins_code'],
            r['end_label_comp_id'],
            r['end_label_asym_id'],
            r['end_label_seq_id'],
            r['pdbx_end_PDB_ins_code'],
            r['beg_auth_comp_id'],
            r['beg_auth_asym_id'],
            r['beg_auth_seq_id'],
            r['end_auth_comp_id'],
            r['end_auth_asym_id'],
            r['end_auth_seq_id'],
        ])


def _run_mkdssp(input_file: Path, output_file: Path) -> bool:
    """Run the external mkdssp binary, preserving the legacy subprocess path."""
    try:
        result = subprocess.run(
            ["mkdssp", str(input_file), str(output_file)],
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        logger.error(
            "mkdssp executable not found. Install DSSP or choose "
            "--algorithm psea, --algorithm pydssp, or --algorithm tmalign."
        )
        return False

    if result.returncode != 0:
        logger.error(f"mkdssp failed for {input_file}: {result.stderr}")
        return False
    return True


def run_dssp(
    input_file: Path,
    output_file: Path,
    algorithm: Algorithm = DEFAULT_ALGORITHM,
    device: str = "cpu",
) -> bool:
    """
    Compute secondary structure and write annotated CIF.

    Reads the input CIF with gemmi (preserving all existing categories),
    computes SSE with the chosen algorithm, adds _struct_conf / _struct_conf_type,
    and writes the annotated CIF.
    """
    if not str(output_file).endswith(".cif"):
        raise ValueError("Output file must have a .cif extension")

    if algorithm == "mkdssp":
        return _run_mkdssp(input_file, output_file)

    try:
        doc = gemmi.cif.read(str(input_file))
        block = doc.sole_block()
        structure = gemmi.make_structure_from_block(block)

        if algorithm == "psea":
            sse = _compute_sse_psea(input_file)
            if sse is None:
                logger.error(f"P-SEA failed for {input_file}")
                return False
            residue_info = _extract_polymer_residue_info(structure)
            if len(sse) != len(residue_info):
                logger.error(
                    f"Residue count mismatch for P-SEA in {input_file}: "
                    f"biotite={len(sse)}, gemmi={len(residue_info)}"
                )
                return False
            code_map = PSEA_CODE_MAP

        elif algorithm == "pydssp":
            coords, residue_info = _extract_backbone_coords(structure)
            if coords is None:
                logger.error(f"Backbone extraction failed for {input_file}")
                return False
            sse = _compute_sse_pydssp(coords, device=device)
            code_map = PYDSSP_CODE_MAP

        elif algorithm == "tmalign":
            ca_coords, residue_info = _extract_ca_coords(structure)
            if ca_coords is None:
                logger.error(f"CA extraction failed for {input_file}")
                return False
            sse = _compute_sse_tmalign(ca_coords)
            code_map = TMALIGN_CODE_MAP

        else:
            logger.error(f"Unknown algorithm: {algorithm}")
            return False

        ranges = _find_secondary_structure_ranges(sse, residue_info, code_map)
        _add_struct_conf_to_cif(block, ranges, algorithm)

        doc.write_file(str(output_file))
        return True

    except ModuleNotFoundError as e:
        logger.error(f"Secondary structure failed for {input_file}: {e}")
        return False
    except Exception as e:
        logger.error(f"Secondary structure failed for {input_file}: {e}")
        import traceback
        traceback.print_exc()
        return False


def _process_single_dssp(args: Tuple[Path, Path]) -> Tuple[str, bool]:
    """Worker function for batch DSSP. Works with both ProcessPool and ThreadPool."""
    input_file, output_dir = args
    output_file = output_dir / input_file.name
    success = run_dssp(
        input_file, output_file,
        algorithm=_CURRENT_ALGORITHM, device=_CURRENT_DEVICE,
    )
    return input_file.name, success


def run_batch_dssp(
    input_dir: Path,
    output_dir: Path,
    workers: int = 8,
    pattern: str = "*.cif",
    algorithm: Algorithm = DEFAULT_ALGORITHM,
    device: str = "cpu",
) -> Tuple[int, int]:
    """
    Batch secondary structure processing with parallel execution.

    When device is 'cuda', uses ThreadPoolExecutor (CUDA contexts are not
    fork-safe) and runs PyDSSP H-bond / helix / strand math on GPU.
    Otherwise uses ProcessPoolExecutor for CPU parallelism.
    """
    global _CURRENT_ALGORITHM, _CURRENT_DEVICE
    _CURRENT_ALGORITHM = algorithm
    _CURRENT_DEVICE = device

    input_files = list(input_dir.glob(pattern))
    if not input_files:
        logger.warning(f"No files matching '{pattern}' in {input_dir}")
        return 0, 0

    output_dir.mkdir(parents=True, exist_ok=True)
    work_items = [(f, output_dir) for f in input_files]

    algo_name = _get_algorithm_label(algorithm)
    device_label = f"GPU ({device})" if device != "cpu" else "CPU"
    logger.info(
        f"Processing {len(input_files)} files with {workers} workers "
        f"using {algo_name} on {device_label}"
    )

    success_count = 0
    error_count = 0

    use_gpu = device != "cpu" and _torch_cuda_available()
    PoolClass = ThreadPoolExecutor if use_gpu else ProcessPoolExecutor

    with PoolClass(max_workers=workers) as executor:
        futures = {
            executor.submit(_process_single_dssp, item): item
            for item in work_items
        }
        for future in as_completed(futures):
            fname, ok = future.result()
            if ok:
                success_count += 1
                if success_count % 100 == 0:
                    logger.info(f"Processed {success_count}/{len(input_files)} files")
            else:
                error_count += 1
                logger.error(f"[FAILED] {fname}")

    logger.info(f"Batch complete: {success_count} success, {error_count} errors")
    return success_count, error_count


def run_command(input_path, output_path):
    """Legacy subprocess wrapper - kept for compatibility but not recommended."""
    command = ["mkdssp", str(input_path), str(output_path)]
    result = subprocess.run(command, capture_output=True, text=True)
    return result.returncode, input_path, result.stdout, result.stderr
