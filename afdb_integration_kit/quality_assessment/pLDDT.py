from __future__ import annotations
import math

import orjson
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

# Reuse naming patterns so we discover sibling files by AFID/version
from .naming import PATTERNS, REQUIRED_TYPES, RE_AFID_ANYWHERE
try:  # Backwards compatibility with legacy validation registry
    from afdb_integration_kit.validation.registry import ValidationHook, register_validator  # type: ignore
except ImportError:  # pragma: no cover
    class ValidationHook:  # type: ignore[no-redef]
        def __init__(self, name, run, formatter=None, description=None, default_kwargs=None):
            self.name = name
            self.run = run
            self.formatter = formatter
            self.description = description
            self.default_kwargs = default_kwargs or {}

        def build_kwargs(self, overrides=None):
            data = dict(self.default_kwargs)
            if overrides:
                data.update(overrides)
            return data

    def register_validator(hook):  # type: ignore[no-redef]
        return None

# Try to import gemmi for structure residue counts and B-factor checks.
# If not available, we degrade gracefully.
try:
    import gemmi  # type: ignore
    _HAVE_GEMMI = True
except Exception:
    _HAVE_GEMMI = False


def _find_all_files(root: Path) -> List[Path]:
    return [p for p in root.rglob("*") if p.is_file()]


def _group_entries(files: List[Path]) -> Dict[Tuple[str, str], Dict[str, Path]]:
    entries: Dict[Tuple[str, str], Dict[str, Path]] = {}
    watch = ("plddt", "cif", "pae")
    for p in files:
        name = p.name
        for kind, pat in PATTERNS.items():
            if kind not in watch:
                continue
            m = pat.match(name)
            if not m:
                continue
            afid, ver = m.group(1), m.group(2)
            key = (afid, ver)
            bucket = entries.setdefault(key, {})
            if kind == "cif":
                bucket["cif"] = p
            else:
                bucket[kind] = p
            break
    return entries



def _filter_entries(
    entries: Dict[Tuple[str, str], Dict[str, Path]],
    ids_pairs: Optional[Set[Tuple[str, str]]] = None,
    ids_afids: Optional[Set[str]] = None,
) -> Dict[Tuple[str, str], Dict[str, Path]]:
    if not ids_pairs and not ids_afids:
        return entries
    out: Dict[Tuple[str, str], Dict[str, Path]] = {}
    for (afid, ver), kinds in entries.items():
        if ids_pairs and (afid, ver) in ids_pairs:
            out[(afid, ver)] = kinds
        elif ids_afids and afid in ids_afids:
            out[(afid, ver)] = kinds
    return out


def _load_json(path: Path) -> Tuple[bool, Any, Optional[str]]:
    try:
        return True, orjson.loads(path.read_bytes()), None
    except Exception as e:
        return False, None, str(e)


def _extract_plddt_arrays(obj: Any) -> Tuple[Optional[List[int]], Optional[List[float]], Optional[List[str]], str]:
    """
    Accepts either:
      - AFDB-style object with keys: residueNumber, confidenceScore, confidenceCategory (optional)
      - Legacy bare array [float, ...] (we synthesize residueNumber 1..N, no categories)
    Returns (residue_numbers, scores, categories, mode)
    """
    if isinstance(obj, dict):
        rn = obj.get("residueNumber")
        cs = obj.get("confidenceScore")
        cc = obj.get("confidenceCategory")
        if isinstance(rn, list) and isinstance(cs, list) and (cc is None or isinstance(cc, list)):
            return rn, cs, cc, "object"
    if isinstance(obj, list):
        # legacy bare array of numbers
        if all(isinstance(x, (int, float)) for x in obj) and len(obj) > 0:
            n = len(obj)
            return list(range(1, n + 1)), [float(x) for x in obj], None, "array"
    return None, None, None, "unknown"


def _check_categories(cc: Optional[List[str]]) -> Tuple[bool, Optional[Set[str]]]:
    if cc is None:
        return True, None
    allowed = {"V", "H", "M", "L", "D"}
    vals = set(cc)
    return vals.issubset(allowed), vals - allowed


def _check_scores_range(scores: List[float]) -> Tuple[bool, List[int]]:
    bad: List[int] = []
    for i, v in enumerate(scores):
        if not (isinstance(v, (int, float)) and math.isfinite(v) and 0.0 <= float(v) <= 100.0):
            bad.append(i)
    return (len(bad) == 0), bad


def _check_residue_numbers(rn: List[int]) -> Tuple[bool, Optional[str]]:
    if len(rn) == 0:
        return False, "empty residueNumber"
    # They specify 1-based sequential list; be forgiving but still strict by default
    for i, val in enumerate(rn, start=1):
        if not isinstance(val, int):
            return False, f"non-integer residueNumber at index {i-1}"
        if val != i:
            return False, f"non-sequential residueNumber at index {i-1} (got {val}, expected {i})"
    return True, None


def _structure_residue_count(struct_path: Path) -> Tuple[Optional[int], Optional[str]]:
    if not _HAVE_GEMMI:
        return None, "gemmi_not_available"
    try:
        st = gemmi.read_structure(str(struct_path))

        def count_if(pred) -> int:
            c = 0
            for model in st:
                for chain in model:
                    for res in chain:
                        try:
                            if pred(res):
                                c += 1
                        except Exception:
                            pass
            return c

        # 1) Preferred: amino-acid residues (works for AF mmCIFs)
        n = count_if(lambda r: r.is_amino_acid())
        if n > 0:
            return n, None

        # 2) Fallback: residues that contain a CA atom (case-insensitive)
        n = count_if(lambda r: any(a.name.upper() == "CA" for a in r))
        if n > 0:
            return n, None

        # 3) Fallback: residues with a defined numeric seqid
        n = count_if(lambda r: getattr(r, "seqid", None) and getattr(r.seqid, "num", None) is not None)
        if n > 0:
            return n, None

        return None, "no_residues_matched"
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"


def _bfactor_min_max(struct_path: Path) -> Optional[Tuple[float, float]]:
    if not _HAVE_GEMMI:
        return None
    try:
        st = gemmi.read_structure(str(struct_path))
        vals: List[float] = []
        has_ca = False

        # Prefer CA atoms from amino-acid residues
        for model in st:
            for chain in model:
                for res in chain:
                    if not res.is_amino_acid():
                        continue
                    for atom in res:
                        if atom.name.upper() == "CA":
                            has_ca = True
                            vals.append(float(atom.b_iso))

        # If no CA collected, fall back to all atoms in amino-acid residues
        if not has_ca:
            for model in st:
                for chain in model:
                    for res in chain:
                        if not res.is_amino_acid():
                            continue
                        for atom in res:
                            vals.append(float(atom.b_iso))

        if not vals:
            return None
        return (min(vals), max(vals))
    except Exception:
        return None


def _pae_dim(pae_path: Path) -> Optional[int]:
    ok, obj, _ = _load_json(pae_path)
    if not ok:
        return None
    # Accept either top-level matrix or a dict with a matrix-like value
    cand = None
    if isinstance(obj, list) and obj and isinstance(obj[0], list):
        cand = obj
    elif isinstance(obj, dict):
        for k in ("predicted_aligned_error", "pae", "data", "matrix"):
            v = obj.get(k)
            if isinstance(v, list) and v and isinstance(v[0], list):
                cand = v
                break
    if cand is None:
        return None
    n = len(cand)
    # quick sanity: first few rows are length n
    for row in cand[:5]:
        if not isinstance(row, list) or len(row) != n:
            return None
    return n


def validate_plddt_file(
    plddt_path: Path,
    *,
    struct_path: Optional[Path] = None,
    pae_path: Optional[Path] = None,
    bfactor_tolerance: float = 1.0,
) -> Dict[str, Any]:
    """
    Validate a single pLDDT JSON.
    """
    res: Dict[str, Any] = {
        "path": str(plddt_path),
        "parse_ok": False,
        "mode": None,  # 'object' or 'array'
        "length": None,
        "range_ok": None,
        "out_of_range_indices": [],
        "resnums_ok": None,
        "categories_ok": None,
        "invalid_categories": None,
        "structure_residue_count": None,
        "len_matches_structure": None,
        "pae_dim": None,
        "len_matches_pae": None,
        "bfactor_min_max": None,
        "bfactor_matches_json_minmax": None,
        "status": "FAIL",
        "notes": [],
    }

    ok, obj, err = _load_json(plddt_path)
    if not ok:
        res["notes"].append(f"JSON parse error: {err}")
        return res
    rn, scores, cats, mode = _extract_plddt_arrays(obj)
    if rn is None or scores is None:
        res["notes"].append("Unrecognized structure: expected AFDB object or bare array")
        return res

    res["parse_ok"] = True
    res["mode"] = mode
    res["length"] = len(scores)

    # number range
    rng_ok, bad_idx = _check_scores_range(scores)
    res["range_ok"] = rng_ok
    res["out_of_range_indices"] = bad_idx

    # residue numbers
    rn_ok, rn_msg = _check_residue_numbers(rn)
    res["resnums_ok"] = rn_ok
    if not rn_ok and rn_msg:
        res["notes"].append(rn_msg)

    # categories
    cat_ok, invalid = _check_categories(cats)
    res["categories_ok"] = cat_ok
    res["invalid_categories"] = list(invalid) if invalid else None

    # structure residue count
    if struct_path is not None:
        nres, serr = _structure_residue_count(struct_path)
        res["had_structure_file"] = True
        res["structure_residue_count"] = nres
        res["structure_error"] = serr
    else:
        res["had_structure_file"] = False

    # PAE dimension
    if pae_path is not None:
        n = _pae_dim(pae_path)
        res["pae_dim"] = n
        if n is not None:
            res["len_matches_pae"] = (n == len(scores))

    # B-factor min/max spot-check
    if struct_path is not None:
        bmm = _bfactor_min_max(struct_path)
        res["bfactor_min_max"] = bmm
        if bmm is not None:
            js_min = min(scores)
            js_max = max(scores)
            ok_min = abs(bmm[0] - js_min) <= bfactor_tolerance
            ok_max = abs(bmm[1] - js_max) <= bfactor_tolerance
            res["bfactor_matches_json_minmax"] = (ok_min and ok_max)

    # status decision
    hard_checks = [
        res["parse_ok"] is True,
        res["range_ok"] is True,
        res["resnums_ok"] is True,
    ]
    # If we had a structure, require length match
    if res["structure_residue_count"] is not None:
        hard_checks.append(res["len_matches_structure"] is True)
    # If we had a PAE, require length match
    if res["pae_dim"] is not None:
        hard_checks.append(res["len_matches_pae"] is True)
    # categories are recommended, not strictly required
    res["status"] = "PASS" if all(hard_checks) else "FAIL"
    return res


def validate_plddt_dataset(
    root: Path,
    *,
    ids_pairs: Optional[Set[Tuple[str, str]]] = None,
    ids_afids: Optional[Set[str]] = None,
    skip_pae: bool = False,
    bfactor_tolerance: float = 1.0,
    with_structure: bool = False,
) -> Tuple[bool, Dict[str, Any]]:
    root = root.resolve()
    files = _find_all_files(root)
    entries_all = _group_entries(files)
    entries = _filter_entries(entries_all, ids_pairs, ids_afids)

    results: List[Dict[str, Any]] = []
    for (afid, ver), kinds in sorted(entries.items()):
        if "plddt" not in kinds:
            continue
        plddt = kinds["plddt"]
        struct = kinds.get("cif") if with_structure else None
        pae = None if skip_pae else kinds.get("pae")

        r = validate_plddt_file(
            plddt,
            struct_path=struct,
            pae_path=pae,
            bfactor_tolerance=bfactor_tolerance,
        )
        # Ensure the flag is always present for the summary
        r.update({
            "afid": afid,
            "version": ver,
            "had_structure_file": struct is not None,
        })
        results.append(r)

    # summary buckets
    total = len(results)
    passed = sum(1 for r in results if r["status"] == "PASS")
    failed = total - passed

    buckets = {
        "parse_error": [r for r in results if not r["parse_ok"]],
        "bad_range": [r for r in results if r["parse_ok"] and r["range_ok"] is False],
        "bad_resnums": [r for r in results if r["parse_ok"] and r["resnums_ok"] is False],
        "len_mismatch_struct": [r for r in results if r.get("structure_residue_count") is not None and r.get("len_matches_structure") is False],
        "len_mismatch_pae": [r for r in results if r.get("pae_dim") is not None and r.get("len_matches_pae") is False],
        "bfactor_mismatch": [r for r in results if r.get("bfactor_min_max") is not None and r.get("bfactor_matches_json_minmax") is False],
        "structure_parse_error": [r for r in results if
                                  r.get("had_structure_file") and r.get("structure_residue_count") is None and r.get(
                                      "structure_error")],
        "structure_zero_count": [r for r in results if
                                 r.get("had_structure_file") and r.get("structure_residue_count") is None and not r.get(
                                     "structure_error")],
    }

    report = {
        "dataset_root": str(root),
        "summary": {
            "files_checked": total,
            "passed": passed,
            "failed": failed,
            "have_structure_files": sum(1 for r in results if r.get("had_structure_file")),
            "have_structure": sum(1 for r in results if r.get("structure_residue_count") is not None),
            "have_pae": sum(1 for r in results if r.get("pae_dim") is not None),
        },
        "buckets": {k: len(v) for k, v in buckets.items()},
        "results": results,
    }
    return (failed == 0), report


def _join_examples(ids: List[str], limit: int) -> str:
    if not ids:
        return ""
    head = ", ".join(ids[:limit])
    tail = len(ids) - limit
    return (head + (f" ... (+{tail} more)" if tail > 0 else ""))


def format_human_plddt(
    report: Dict[str, Any],
    *,
    errors_only: bool = False,
    verbose: bool = False,
    limit: int = 5,
) -> str:
    """
    Human-readable summary for pLDDT checks.
    """
    lines: List[str] = []
    lines.append(f"# pLDDT check: {report['dataset_root']}")
    s = report["summary"]
    lines.append(f"Files: {s['files_checked']}  Passed: {s['passed']}  Failed: {s['failed']}")
    if 'have_structure_files' in s:
        lines.append(
            f"With structure files: {s['have_structure_files']}  Parsed: {s['have_structure']}  With PAE: {s['have_pae']}")
    else:
        lines.append(f"With structure: {s['have_structure']}  With PAE: {s['have_pae']}")

    b = report["buckets"]
    if any(v > 0 for v in b.values()):
        lines.append("\nIssues by type:")
        for key, label in [
            ("parse_error", "JSON parse error"),
            ("bad_range", "scores out of 0–100"),
            ("bad_resnums", "residueNumber not 1..N"),
            ("len_mismatch_struct", "length ≠ structure residues"),
            ("len_mismatch_pae", "length ≠ PAE dimension"),
            ("bfactor_mismatch", "B-factor min/max ≠ JSON min/max"),
            ("structure_parse_error", "structure parse failed"),
            ("structure_zero_count", "structure had zero residue count"),
        ]:
            if b.get(key, 0):
                lines.append(f"  - {label}: {b[key]}")

    results = report["results"]
    if not results:
        return "\n".join(lines)

    if errors_only or verbose:
        lines.append("")
        to_print = results
        if errors_only:
            to_print = [r for r in results if r["status"] != "PASS"]
        count = 0
        cap = None if verbose else limit
        for r in to_print:
            afid = r["afid"]; ver = r["version"]
            bits = []
            if not r["parse_ok"]:
                bits.append("parse_error")
            if r["range_ok"] is False:
                idxs = r["out_of_range_indices"][:3]
                more = len(r["out_of_range_indices"]) - len(idxs)
                bits.append(f"bad_range idx={idxs}{'..' if more>0 else ''}")
            if r["resnums_ok"] is False:
                bits.append("bad_resnums")
            if r.get("structure_residue_count") is not None and r.get("len_matches_structure") is False:
                bits.append(f"len≠struct ({r['length']} vs {r['structure_residue_count']})")
            if r.get("pae_dim") is not None and r.get("len_matches_pae") is False:
                bits.append(f"len≠PAE ({r['length']} vs {r['pae_dim']})")
            if r.get("bfactor_min_max") is not None and r.get("bfactor_matches_json_minmax") is False:
                bits.append(f"bfactor!=json_minmax (struct {r['bfactor_min_max']})")
            status = r["status"]
            lines.append(f"{afid} {ver}: {status}  {'; '.join(bits) if bits else 'ok'}")
            count += 1
            if cap is not None and count >= cap:
                remaining = len(to_print) - count
                if remaining > 0:
                    lines.append(f"... ({remaining} more)")
                break

    return "\n".join(lines)


register_validator(
    ValidationHook(
        name="plddt",
        run=validate_plddt_dataset,
        formatter=format_human_plddt,
        description="Validate pLDDT JSON files and optional cross-checks against structure/PAE data.",
    )
)
