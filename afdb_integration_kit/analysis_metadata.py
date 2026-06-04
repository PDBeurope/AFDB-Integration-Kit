from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence


ANALYSIS_METADATA_FIELDS = [
    "timestamp",
    "run_name",
    "model_id",
    "original_id",
    "upload_status",
    "passes_quality_threshold",
    "ipsae_max",
    "pdockq2_max",
    "quality_ipsae_threshold",
    "quality_pdockq2_threshold",
    "source_archive",
    "shard_id",
    "task_id",
    "batch_id",
    "batch_started_at",
    "batch_finished_at",
    "failure_reason",
    "expected_output_files_json",
    "scores_json",
]


def _build_analysis_metadata_rows(
    *,
    work_dir: Path,
    run_name: str,
    original_ids: Sequence[str],
    model_ids: Sequence[str],
    status_by_id: Dict[str, str],
    failure_reasons: Dict[str, str],
    source_archive: str,
    shard_id: int,
    task_id: int,
    batch_id: int,
    batch_started_at: str,
    batch_finished_at: str,
    ipsae_threshold: float,
    pdockq2_threshold: float,
    include_scores: bool,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    ipsae_by_model = _read_ipsae_summary(work_dir / "ipsae" / "ipsae_summary.csv")

    for index, model_id in enumerate(model_ids):
        scores = ipsae_by_model.get(model_id, {}).copy()
        scores.update(_read_clash_scores(work_dir, model_id))
        scores.update(_read_interface_scores(work_dir, model_id))

        ipsae_max = _max_score(scores, "ipsae_")
        pdockq2_max = _max_score(scores, "pDockQ2_")
        passes = (
            ipsae_max is not None
            and pdockq2_max is not None
            and ipsae_max >= ipsae_threshold
            and pdockq2_max >= pdockq2_threshold
        )

        row = {
            "timestamp": batch_finished_at,
            "run_name": run_name,
            "model_id": model_id,
            "original_id": original_ids[index] if index < len(original_ids) else "",
            "upload_status": status_by_id.get(model_id, ""),
            "passes_quality_threshold": "true" if passes else "false",
            "ipsae_max": ipsae_max,
            "pdockq2_max": pdockq2_max,
            "quality_ipsae_threshold": ipsae_threshold,
            "quality_pdockq2_threshold": pdockq2_threshold,
            "source_archive": source_archive,
            "shard_id": shard_id,
            "task_id": task_id,
            "batch_id": batch_id,
            "batch_started_at": batch_started_at,
            "batch_finished_at": batch_finished_at,
            "failure_reason": failure_reasons.get(model_id, ""),
            "expected_output_files_json": json.dumps(_expected_output_files(model_id)),
            "scores_json": json.dumps(scores if include_scores else {}, sort_keys=True),
        }
        rows.append(row)

    return rows


def _score_keys(csv_path: Path) -> List[str]:
    keys = set()
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                scores = json.loads(row.get("scores_json") or "{}")
            except json.JSONDecodeError:
                continue
            keys.update(scores.keys())
    return sorted(keys)


def _coerce_row(row: Dict[str, str], score_keys: Iterable[str]) -> Dict[str, Any]:
    parsed: Dict[str, Any] = dict(row)
    for key in score_keys:
        scores = _loads_scores(row)
        parsed[f"score_{key}"] = _coerce_float(scores.get(key))
    return parsed


def _read_ipsae_summary(path: Path) -> Dict[str, Dict[str, float]]:
    if not path.exists():
        return {}

    by_model: Dict[str, Dict[str, float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            model_id = _model_id_from_path(row.get("pdb_path", ""))
            if not model_id:
                continue
            scores = {
                key: value
                for key, value in row.items()
                if key != "pdb_path" and _coerce_float(value) is not None
            }
            by_model[model_id] = {
                key: _coerce_float(value) for key, value in scores.items()
            }
    return by_model


def _read_clash_scores(work_dir: Path, model_id: str) -> Dict[str, int]:
    path = work_dir / "clash_interface_analysis" / f"{model_id}-model_v1_clashes.json"
    data = _read_json(path)
    scores: Dict[str, int] = {}
    for site in data.get("sites", []):
        label = site.get("label", "")
        annotations = site.get("additional_site_annotations", {})
        n_clashes = annotations.get("n_clashes")
        if label == "backbone_clashes":
            scores["N_clash_backbone"] = n_clashes
        elif label == "heavy_atom_clashes":
            scores["N_clash_heavyAtom"] = n_clashes
    return {key: value for key, value in scores.items() if value is not None}


def _read_interface_scores(work_dir: Path, model_id: str) -> Dict[str, int]:
    path = work_dir / "clash_interface_analysis" / f"{model_id}-model_v1_interface.json"
    data = _read_json(path)
    total = 0
    for site in data.get("sites", []):
        annotations = site.get("additional_site_annotations", {})
        total += len(annotations.get("interactions", []))
    return {"N_interface_interactions": total} if total else {}


def _read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}


def _model_id_from_path(path: str) -> str:
    stem = Path(path).name
    suffix = "-model_v1.pdb"
    return stem[: -len(suffix)] if stem.endswith(suffix) else Path(stem).stem


def _max_score(scores: Dict[str, Any], prefix: str) -> float | None:
    values = [
        value
        for key, raw_value in scores.items()
        if key.startswith(prefix)
        for value in [_coerce_float(raw_value)]
        if value is not None
    ]
    return max(values) if values else None


def _coerce_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _loads_scores(row: Dict[str, str]) -> Dict[str, Any]:
    try:
        return json.loads(row.get("scores_json") or "{}")
    except json.JSONDecodeError:
        return {}


def _expected_output_files(model_id: str) -> List[str]:
    return [
        f"{model_id}-model_v1.cif",
        f"{model_id}-model_v1.bcif",
        f"{model_id}-model_v1.pdb",
        f"{model_id}-confidence_v1.json",
        f"{model_id}-predicted_aligned_error_v1.json",
        f"{model_id}-model_v1_clashes.json",
        f"{model_id}-model_v1_interface.json",
    ]
