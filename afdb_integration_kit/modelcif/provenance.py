from __future__ import annotations

from typing import Any


DEFAULT_ALPHAFOLD2_VERSION = "2.3.2"
DEFAULT_IPSAE_PAE_CUTOFF = "10.0"
DEFAULT_IPSAE_DIST_CUTOFF = "15.0"
SECONDARY_STRUCTURE_DESCRIPTION = (
    "Secondary-structure assignment and annotation extraction from predicted coordinates"
)
IPSAE_DESCRIPTION = "Interface scoring / QA metrics (ipSAE, pDockQ, pDockQ2, LIS)"
IPSAE_PROTOCOL_DETAILS = "Post-processing interface QA metrics computed with ipSAE"


def _as_category(payload: dict[str, Any], category_name: str) -> dict[str, list[Any]]:
    categories = payload.setdefault("categories", {})
    if not isinstance(categories, dict):
        raise TypeError("ModelCIF metadata payload must contain a 'categories' dictionary.")
    category = categories.setdefault(category_name, {})
    if not isinstance(category, dict):
        raise TypeError(f"Category {category_name} must be a dictionary.")
    return category


def _get_list(category: dict[str, Any], key: str) -> list[Any]:
    value = category.get(key, [])
    return value if isinstance(value, list) else []


def _collect_row_dicts(category: dict[str, Any]) -> list[dict[str, Any]]:
    keys = list(category.keys())
    if not keys:
        return []
    lengths = [len(v) for v in category.values() if isinstance(v, list)]
    row_count = max(lengths, default=0)
    rows: list[dict[str, Any]] = []
    for index in range(row_count):
        row: dict[str, Any] = {}
        for key in keys:
            values = category.get(key, [])
            if isinstance(values, list) and index < len(values):
                row[key] = values[index]
        rows.append(row)
    return rows


def _set_rows(category: dict[str, list[Any]], rows: list[dict[str, Any]], fields: list[str]) -> None:
    category.clear()
    for field in fields:
        category[field] = [row.get(field, "?") for row in rows]


def _canonicalize_dssp_algorithm(dssp_algorithm: str | None, software_names: list[str]) -> str:
    if dssp_algorithm:
        lowered = dssp_algorithm.strip().lower()
        if lowered == "pydssp":
            return "pydssp"
        if lowered == "mkdssp":
            return "mkdssp"
    for name in software_names:
        lowered = str(name).strip().lower()
        if lowered == "pydssp":
            return "pydssp"
        if lowered in {"dssp", "mkdssp"}:
            return "mkdssp"
    return "mkdssp"


def _resolve_alphafold_version(
    software_rows: list[dict[str, Any]], *, allow_default: bool
) -> str:
    for row in software_rows:
        raw_name = str(row.get("name", "") or "").strip()
        if raw_name.startswith("AlphaFold"):
            raw_version = str(row.get("version", "") or "").strip()
            if raw_version:
                return raw_version
            break
    if allow_default:
        return DEFAULT_ALPHAFOLD2_VERSION
    raise ValueError(
        "ModelCIF metadata input is missing _software.version for the AlphaFold provenance row."
    )


def _infer_parameter_value(category: dict[str, Any], parameter_name: str, default: str) -> str:
    names = _get_list(category, "name")
    values = _get_list(category, "value")
    for idx, name in enumerate(names):
        if str(name) == parameter_name and idx < len(values):
            value = str(values[idx])
            if value and value not in {"?", "."}:
                return value
    return default


def _detect_is_complex(payload: dict[str, Any], explicit: bool | None) -> bool:
    if explicit is not None:
        return explicit
    chains = payload.get("chains", [])
    if isinstance(chains, list):
        chain_ids = {
            str(chain.get("chain_id"))
            for chain in chains
            if isinstance(chain, dict) and chain.get("chain_id")
        }
        if len(chain_ids) > 1:
            return True
    return False


def normalize_modelcif_provenance(
    payload: dict[str, Any],
    *,
    dssp_algorithm: str | None = None,
    is_complex: bool | None = None,
    allow_default_alphafold_version: bool = False,
) -> None:
    categories = payload.setdefault("categories", {})
    if not isinstance(categories, dict):
        raise TypeError("ModelCIF metadata payload must contain a 'categories' dictionary.")

    software = _as_category(payload, "_software")
    software_rows = _collect_row_dicts(software)
    software_names = [str(row.get("name", "") or "") for row in software_rows]

    detected_complex = _detect_is_complex(payload, is_complex)
    resolved_dssp = _canonicalize_dssp_algorithm(dssp_algorithm, software_names)
    alphafold_version = _resolve_alphafold_version(
        software_rows, allow_default=allow_default_alphafold_version
    )

    parameter_category = _as_category(payload, "_ma_software_parameter")
    pae_cutoff = _infer_parameter_value(
        parameter_category, "pae_cutoff", DEFAULT_IPSAE_PAE_CUTOFF
    )
    dist_cutoff = _infer_parameter_value(
        parameter_category, "dist_cutoff", DEFAULT_IPSAE_DIST_CUTOFF
    )

    dssp_name = "PyDSSP" if resolved_dssp == "pydssp" else "DSSP"
    dssp_type = "library" if resolved_dssp == "pydssp" else "package"
    alphafold_name = "AlphaFold-Multimer" if detected_complex else "AlphaFold"

    software_rows_out: list[dict[str, Any]] = [
        {
            "pdbx_ordinal": "1",
            "name": alphafold_name,
            "version": alphafold_version,
            "type": "package",
            "description": "Structure prediction",
            "classification": "model building",
        }
    ]
    if detected_complex:
        software_rows_out.append(
            {
                "pdbx_ordinal": "2",
                "name": "ipSAE",
                "version": "?",
                "type": "package",
                "description": IPSAE_DESCRIPTION,
                "classification": "data processing",
            }
        )
        software_rows_out.append(
            {
                "pdbx_ordinal": "3",
                "name": dssp_name,
                "version": "?",
                "type": dssp_type,
                "description": SECONDARY_STRUCTURE_DESCRIPTION,
                "classification": "data extraction",
            }
        )
    else:
        software_rows_out.append(
            {
                "pdbx_ordinal": "2",
                "name": dssp_name,
                "version": "?",
                "type": dssp_type,
                "description": SECONDARY_STRUCTURE_DESCRIPTION,
                "classification": "data extraction",
            }
        )

    _set_rows(
        software,
        software_rows_out,
        [
            "pdbx_ordinal",
            "name",
            "version",
            "type",
            "description",
            "classification",
        ],
    )

    software_group = _as_category(payload, "_ma_software_group")
    if detected_complex:
        software_group_rows = [
            {"ordinal_id": "1", "group_id": "1", "software_id": "1"},
            {
                "ordinal_id": "2",
                "group_id": "2",
                "software_id": "2",
                "parameter_group_id": "2",
            },
            {"ordinal_id": "3", "group_id": "3", "software_id": "3"},
        ]
    else:
        software_group_rows = [
            {"ordinal_id": "1", "group_id": "1", "software_id": "1"},
            {"ordinal_id": "2", "group_id": "2", "software_id": "2"},
        ]
    _set_rows(
        software_group,
        software_group_rows,
        ["ordinal_id", "group_id", "software_id", "parameter_group_id"],
    )

    if detected_complex:
        parameter_rows = [
            {
                "parameter_id": "1",
                "group_id": "2",
                "data_type": "float",
                "name": "pae_cutoff",
                "value": pae_cutoff,
                "description": "PAE cutoff used by ipSAE",
            },
            {
                "parameter_id": "2",
                "group_id": "2",
                "data_type": "float",
                "name": "dist_cutoff",
                "value": dist_cutoff,
                "description": "Distance cutoff used by ipSAE",
            },
        ]
        _set_rows(
            parameter_category,
            parameter_rows,
            ["parameter_id", "group_id", "data_type", "name", "value", "description"],
        )
    else:
        categories.pop("_ma_software_parameter", None)

    protocol = _as_category(payload, "_ma_protocol_step")
    if detected_complex:
        protocol_rows = [
            {
                "ordinal_id": "1",
                "protocol_id": "1",
                "step_id": "1",
                "method_type": "modeling",
                "step_name": "model inference",
                "details": "Predicted structure generated with AlphaFold-Multimer",
                "software_group_id": "1",
            },
            {
                "ordinal_id": "2",
                "protocol_id": "1",
                "step_id": "2",
                "method_type": "other",
                "step_name": "interface scoring",
                "details": IPSAE_PROTOCOL_DETAILS,
                "software_group_id": "2",
            },
            {
                "ordinal_id": "3",
                "protocol_id": "1",
                "step_id": "3",
                "method_type": "other",
                "step_name": "secondary structure assignment",
                "details": "Post-processing secondary-structure assignment using DSSP",
                "software_group_id": "3",
            },
        ]
    else:
        protocol_rows = [
            {
                "ordinal_id": "1",
                "protocol_id": "1",
                "step_id": "1",
                "method_type": "modeling",
                "step_name": "model inference",
                "details": "Predicted structure generated with AlphaFold",
                "software_group_id": "1",
            },
            {
                "ordinal_id": "2",
                "protocol_id": "1",
                "step_id": "2",
                "method_type": "other",
                "step_name": "secondary structure assignment",
                "details": "Post-processing secondary-structure assignment using DSSP",
                "software_group_id": "2",
            },
        ]
    _set_rows(
        protocol,
        protocol_rows,
        [
            "ordinal_id",
            "protocol_id",
            "step_id",
            "method_type",
            "step_name",
            "details",
            "software_group_id",
        ],
    )
