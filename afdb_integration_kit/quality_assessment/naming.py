# afdb_integration_kit/quality_assessment/naming.py
from __future__ import annotations
import re
from pathlib import Path
from typing import Dict, List, Tuple, Any, Iterable, Set

# Canonical per-entry patterns
PATTERNS = {
    "pdb":   re.compile(r"^(AF-\d{16})-model-(v\d+)\.pdb$"),
    "cif":   re.compile(r"^(AF-\d{16})-model-(v\d+)\.cif$"),
    "bcif":  re.compile(r"^(AF-\d{16})-model-(v\d+)\.bcif$"),
    "plddt": re.compile(r"^(AF-\d{16})-confidence-(v\d+)\.json$"),
    "pae":   re.compile(r"^(AF-\d{16})-predicted_aligned_error-(v\d+)\.json$"),
    "msa":   re.compile(r"^(AF-\d{16})-msa-(v\d+)\.a3m$"),
    # Dataset-level model metadata batches
    "model_batch": re.compile(r"^AF-metadata-(\d+)-of-(\d+)\.json$"),
}

REQUIRED_TYPES = ["pdb", "cif", "bcif", "plddt", "pae", "msa"]
RE_AFID_ANYWHERE = re.compile(r"(AF-\d{16})")
RE_VERSION_SUFFIX = re.compile(r"-v\d+\b")
RE_PROVIDER_JSON_NAME = re.compile(r"^[A-Za-z0-9._-]+\.json$")
FASTA_NAME = "sequences.fasta"

# Accept lines like: AF-0000000000000001-v2  or AF-0000000000000001 v2  or AF-0000000000000001
RE_AFID_VER_LINE = re.compile(r"^\s*(AF-\d{16})(?:[ \t_-]*(v\d+))?\s*$")


def get_required_types() -> List[str]:
    return list(REQUIRED_TYPES)


def parse_ids_file(path: Path) -> tuple[Set[tuple[str, str]], Set[str]]:
    """
    Returns:
      - pairs: set of (AFID, vN) explicitly requested
      - afids: set of AFIDs with no version specified
    Blank lines and lines starting with '#' are ignored.
    """
    pairs: Set[tuple[str, str]] = set()
    afids: Set[str] = set()
    for line in path.read_text().splitlines():
        if not line or line.lstrip().startswith("#"):
            continue
        m = RE_AFID_VER_LINE.match(line)
        if not m:
            continue
        afid, ver = m.group(1), m.group(2)
        if ver:
            pairs.add((afid, ver))
        else:
            afids.add(afid)
    return pairs, afids


def _find_all_files(root: Path) -> List[Path]:
    return [p for p in root.rglob("*") if p.is_file()]


def _group_entries(files: List[Path]) -> Dict[Tuple[str, str], Dict[str, Path]]:
    """
    Group by (AFID, vN) using canonical per-entry patterns only.
    Returns: {(AF-<16>, vN): {kind: Path, ...}, ...}
    """
    entries: Dict[Tuple[str, str], Dict[str, Path]] = {}
    for p in files:
        name = p.name
        for kind, pat in PATTERNS.items():
            if kind == "model_batch":
                continue
            m = pat.match(name)
            if m:
                afid, ver = m.group(1), m.group(2)
                key = (afid, ver)
                entries.setdefault(key, {})
                entries[key][kind] = p
                break
    return entries


def _find_non_compliant(files: List[Path]) -> List[Path]:
    """
    Files that contain an AF-<16> token but do NOT match any canonical per-entry pattern.
    """
    canonical = [pat for k, pat in PATTERNS.items() if k != "model_batch"]
    non_compliant = []
    for p in files:
        name = p.name
        if not RE_AFID_ANYWHERE.search(name):
            continue
        if any(pat.match(name) for pat in canonical):
            continue
        non_compliant.append(p)
    return non_compliant


def _filter_entries(
    entries: Dict[Tuple[str, str], Dict[str, Path]],
    wanted_pairs: Set[tuple[str, str]] | None,
    wanted_afids: Set[str] | None,
) -> Dict[Tuple[str, str], Dict[str, Path]]:
    if not wanted_pairs and not wanted_afids:
        return entries
    out: Dict[Tuple[str, str], Dict[str, Path]] = {}
    for (afid, ver), kinds in entries.items():
        if wanted_pairs and (afid, ver) in wanted_pairs:
            out[(afid, ver)] = kinds
        elif wanted_afids and afid in wanted_afids:
            out[(afid, ver)] = kinds
    return out


def validate_dataset_naming(
    root: Path,
    *,
    ids_pairs: Set[tuple[str, str]] | None = None,
    ids_afids: Set[str] | None = None,
) -> tuple[bool, dict]:
    """
    Perform naming + presence checks only. Returns (ok, report_dict).
    If ids_pairs/ids_afids are provided, only those entries are counted and printed.
    """
    root = root.resolve()
    if not root.is_dir():
        return False, {"error": f"{root} is not a directory"}

    files = _find_all_files(root)
    all_entries = _group_entries(files)
    entries = _filter_entries(all_entries, ids_pairs, ids_afids)

    # Dataset-level
    fasta_paths = [p for p in files if p.name == FASTA_NAME]
    batch_files = [p for p in files if PATTERNS["model_batch"].match(p.name)]

    # Summarise model-metadata batches in a tolerant way
    totals_map: Dict[int, List[int]] = {}
    for p in batch_files:
        m = PATTERNS["model_batch"].match(p.name)
        if not m:
            continue
        start = int(m.group(1))
        total = int(m.group(2))
        totals_map.setdefault(total, []).append(start)

    batch_summary: Dict[str, Any] = {
        "count": len(batch_files),
        # map of total -> how many batch files with that total
        "totals": {str(t): len(starts) for t, starts in sorted(totals_map.items())},
        "single_total_stats": None,  # filled only if exactly one total is present
        "notes": [],
    }

    # If exactly one 'total' value exists, add light stats
    if len(totals_map) == 1:
        (only_total, starts) = next(iter(totals_map.items()))
        starts_sorted = sorted(starts)
        single_stats = {
            "total": only_total,
            "starts_min": int(starts_sorted[0]) if starts_sorted else None,
            "starts_max": int(starts_sorted[-1]) if starts_sorted else None,
            "unique_starts": len(set(starts)),
        }
        # Optional heuristic: estimate typical step (batch size) if we have ≥2 files
        if len(starts_sorted) >= 2:
            diffs = [b - a for a, b in zip(starts_sorted, starts_sorted[1:]) if b > a]
            if diffs:
                diffs_sorted = sorted(diffs)
                mid = len(diffs_sorted) // 2
                typical_step = diffs_sorted[mid] if len(diffs_sorted) % 2 == 1 else (diffs_sorted[mid - 1] +
                                                                                     diffs_sorted[mid]) // 2
                if typical_step > 0:
                    import math
                    expected = math.ceil(only_total / typical_step)
                    batch_summary["notes"].append(
                        f"heuristic step ≈ {typical_step}; observed files {len(starts)}; rough expected ≈ {expected}"
                    )
                single_stats["typical_step"] = typical_step if diffs else None
        batch_summary["single_total_stats"] = single_stats

    # Provider metadata candidates:
    per_entry_jsons = set()
    for (_afid, _ver), kinds in all_entries.items():
        for k, path in kinds.items():
            if k in ("plddt", "pae"):
                per_entry_jsons.add(path)

    provider_candidates: List[Path] = []
    for p in files:
        if p.suffix != ".json":
            continue
        if p in per_entry_jsons:
            continue
        if PATTERNS["model_batch"].match(p.name):
            continue
        if not RE_PROVIDER_JSON_NAME.match(p.name):
            continue
        if RE_AFID_ANYWHERE.search(p.name):
            continue
        if RE_VERSION_SUFFIX.search(p.stem):
            continue
        provider_candidates.append(p)

    # Foldseek index presence by names only
    has_ffindex = any(p.suffix == ".ffindex" for p in files)
    has_ffdata = any(p.suffix == ".ffdata" for p in files)

    # Per-entry presence matrix
    per_entry_reports: List[dict] = []
    any_issue = False

    for (afid, ver) in sorted(entries.keys()):
        kinds = entries[(afid, ver)]
        missing = [k for k in REQUIRED_TYPES if k not in kinds]
        per_entry_reports.append({
            "afid": afid,
            "version": ver,
            "present": {k: (k in kinds) for k in REQUIRED_TYPES},
            "filenames": {k: kinds[k].name for k in kinds},
            "missing": missing,
            "status": "PASS" if not missing else "FAIL",
        })
        if missing:
            any_issue = True

    # If we filtered and ended up with zero selected, that is still a signal
    if not per_entry_reports:
        any_issue = True

    # Non-compliant AF-* filenames (scanned on all files)
    odd = _find_non_compliant(files)

    # Dataset-level issues
    provider_status = "OK" if len(provider_candidates) == 1 else ("MISSING" if len(provider_candidates) == 0 else "MULTIPLE")
    fasta_status = "OK" if len(fasta_paths) == 1 else ("MISSING" if len(fasta_paths) == 0 else "MULTIPLE")

    if provider_status != "OK":
        any_issue = True
    if fasta_status != "OK":
        any_issue = True
    # Foldseek not required

    report: dict[str, Any] = {
        "dataset_root": str(root),
        "sequences_fasta": {
            "status": fasta_status,
            "paths": [str(p.relative_to(root)) for p in fasta_paths],
        },
        "provider_metadata": {
            "status": provider_status,
            "candidates": [str(p.relative_to(root)) for p in provider_candidates],
        },
        "model_metadata_batches": batch_summary,
        "foldseek_index": {
            "status": "PRESENT" if (has_ffindex and has_ffdata) else ("PARTIAL" if (has_ffindex or has_ffdata) else "MISSING"),
        },
        "entries": per_entry_reports,
        "non_compliant": [str(p.relative_to(root)) for p in sorted(odd)],
        "summary": {
            "entry_count_total": len(all_entries),
            "entry_count_selected": len(per_entry_reports),
            "issues": any_issue,
        },
        "status": "PASS" if not any_issue else "FAIL",
    }

    return (not any_issue), report


def _join_examples(ids: List[str], limit: int) -> str:
    if not ids:
        return ""
    show = ids[:limit]
    suffix = f" ... (+{len(ids) - limit} more)" if len(ids) > limit else ""
    return ", ".join(show) + suffix


def format_human(
    report: dict,
    *,
    errors_only: bool = False,
    verbose: bool = False,
    limit: int = 5,
) -> str:
    """
    Human-readable output.
      - Default: compact summary with issue buckets and a few examples
      - errors_only: one-line per failing entry
      - verbose: one-line per entry, including passing ones
    """
    lines: List[str] = []
    lines.append(f"# Dataset: {report['dataset_root']}")
    # FASTA
    fasta = report["sequences_fasta"]
    if fasta["status"] == "OK":
        lines.append(f"- sequences.fasta: PRESENT at {fasta['paths'][0]}")
    elif fasta["status"] == "MISSING":
        lines.append("- sequences.fasta: NOT FOUND")
    else:
        examples = ", ".join(fasta["paths"][:limit])
        more = f" (and {len(fasta['paths'])-limit} more)" if len(fasta["paths"]) > limit else ""
        lines.append(f"- sequences.fasta: MULTIPLE FOUND: {examples}{more}")
    # Provider
    prov = report["provider_metadata"]
    if prov["status"] == "OK":
        lines.append(f"- provider metadata: PRESENT as {prov['candidates'][0]}")
    elif prov["status"] == "MISSING":
        lines.append("- provider metadata: NOT FOUND")
    else:
        examples = ", ".join(prov["candidates"][:limit])
        more = f" (and {len(prov['candidates'])-limit} more)" if len(prov["candidates"]) > limit else ""
        lines.append(f"- provider metadata: MULTIPLE CANDIDATES: {examples}{more}")
    # Batches
    batches = report["model_metadata_batches"]
    if batches["count"] > 0:
        lines.append(f"- model metadata batches: PRESENT ({batches['count']} file(s))")
        if batches["totals"]:
            # e.g. totals seen: 6000×12, 120000×315
            totals_bits = [f"{t}×{c}" for t, c in sorted(batches["totals"].items(), key=lambda x: int(x[0]))]
            lines.append(f"  totals seen: {', '.join(totals_bits)}")
        if batches["single_total_stats"]:
            sts = batches["single_total_stats"]
            lines.append(
                f"  starts (total={sts['total']}): "
                f"min={sts['starts_min']}, max={sts['starts_max']}, unique={sts['unique_starts']}"
            )
        for note in batches.get("notes", []):
            lines.append(f"  note: {note}")
    else:
        lines.append("- model metadata batches (AF-metadata-<start>-of-<total>.json): NOT FOUND")

    entries = report["entries"]
    n_sel = report["summary"]["entry_count_selected"]
    if not entries:
        lines.append("\nNo entries selected. Nothing to report.")
        return "\n".join(lines)

    # Summary
    complete = sum(1 for e in entries if e["status"] == "PASS")
    issues = n_sel - complete
    pct_ok = 0.0 if n_sel == 0 else (100.0 * complete / n_sel)
    lines.append(f"\nEntries: {n_sel}")
    lines.append(f"Complete: {complete} ({pct_ok:.1f}%)")
    lines.append(f"With issues: {issues} ({100.0 - pct_ok:.1f}%)")

    # Bucket by missing type
    missing_buckets: Dict[str, List[str]] = {k: [] for k in REQUIRED_TYPES}
    for e in entries:
        if e["missing"]:
            for t in e["missing"]:
                missing_buckets[t].append(f"{e['afid']}")

    any_bucket = False
    for t in REQUIRED_TYPES:
        ids = missing_buckets[t]
        if not ids:
            continue
        any_bucket = True
        lines.append(f"\nMissing by type:")
        break
    for t in REQUIRED_TYPES:
        ids = sorted(set(missing_buckets[t]))
        if not ids:
            continue
        lines.append(f"  - {t}: {len(ids)} / {n_sel} ({(100.0*len(ids)/n_sel):.1f}%)   e.g. " +
                     _join_examples(ids, limit))

    # Non-compliant filenames
    if report["non_compliant"]:
        lines.append(f"\nNon-compliant filenames: {len(report['non_compliant'])}")
        for name in report["non_compliant"][:limit]:
            lines.append(f"  - {name}")

    # Drill-down modes
    if errors_only or verbose:
        lines.append("")
        to_print = entries
        if errors_only:
            to_print = [e for e in entries if e["missing"]]
        cap = limit if not verbose else None
        count = 0
        for e in to_print:
            miss = ", ".join(e["missing"]) if e["missing"] else "OK"
            lines.append(f"{e['afid']} {e['version']}: {miss}")
            count += 1
            if cap is not None and count >= cap:
                remaining = len(to_print) - count
                if remaining > 0:
                    lines.append(f"... ({remaining} more)")
                break

    return "\n".join(lines)
