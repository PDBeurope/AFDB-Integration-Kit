"""Shared protein-description fallback rules for manifest consumers."""

from __future__ import annotations

from typing import Any, Mapping

import orjson


def _string_values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        values: list[str] = []
        for item in value:
            values.extend(_string_values(item))
        return values
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("[") and text.endswith("]"):
            try:
                return _string_values(orjson.loads(text))
            except orjson.JSONDecodeError:
                pass
        return [text] if text else []
    text = str(value).strip()
    return [text] if text else []


def protein_description(
    manifest_name: str | None,
    entry: Mapping[str, Any],
    accession: str,
) -> str:
    """Return the canonical manifest/UniProt protein description."""
    candidates = (
        _string_values(manifest_name),
        _string_values(entry.get("protein_full_names")),
        _string_values(entry.get("protein_short_names")),
        _string_values(entry.get("entry_name")),
        _string_values(accession),
    )
    for values in candidates:
        if values:
            return values[0]
    return accession
