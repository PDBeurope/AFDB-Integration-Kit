from collections import defaultdict
from typing import Any, Dict, List
import re
from afdb_integration_kit.utils.constant import (
    CAT_ATOM_SITE,
    ITEM_LABEL_ASYM_ID,
    ITEM_AUTH_ASYM_ID,
    ITEM_LABEL_COMP_ID,
    ITEM_AUTH_COMP_ID,
    ITEM_LABEL_SEQ_ID,
    ITEM_AUTH_SEQ_ID,
    CAT_CELL,
    CAT_SYMMETRY,
)
import gemmi
import logging

logger = logging.getLogger("afdb_integration_kit")

# Pre-configured WriteOptions for reuse (avoid recreating per file)
_WRITE_OPTIONS_ALIGNED = gemmi.cif.WriteOptions()
_WRITE_OPTIONS_ALIGNED.align_loops = 50
_WRITE_OPTIONS_ALIGNED.align_pairs = 50
_WRITE_OPTIONS_ALIGNED.prefer_pairs = True

_WRITE_OPTIONS_MINIMAL = gemmi.cif.WriteOptions()
_WRITE_OPTIONS_MINIMAL.align_loops = 0
_WRITE_OPTIONS_MINIMAL.align_pairs = 0
_WRITE_OPTIONS_MINIMAL.prefer_pairs = False


_CIF_RESERVED_PREFIXES = ("data_", "loop_", "save_", "stop_", "global_")
_CIF_SPECIAL_FIRST_CHARS = frozenset("_#$'\"[];")
_QUOTED_VALUE_RE = re.compile(r"'([^']*)'")


def _cif_value_needs_quoting(value: str) -> bool:
    """Return True if the CIF 1.1 spec requires this value to be quoted."""
    if not value:
        return True
    if value[0] in _CIF_SPECIAL_FIRST_CHARS:
        return True
    if any(ch <= " " for ch in value):
        return True
    if value.lower().startswith(_CIF_RESERVED_PREFIXES):
        return True
    return False


def _strip_unnecessary_cif_quotes(filepath: str) -> None:
    """Remove single quotes that gemmi adds but the CIF spec does not require.

    Reads the file, replaces every ``'value'`` token with ``value`` when the
    inner string is a valid unquoted CIF 1.1 data value, then writes back.
    """
    with open(filepath, "r", encoding="utf-8") as fh:
        content = fh.read()

    def _unquote_if_safe(match: re.Match) -> str:
        inner = match.group(1)
        if _cif_value_needs_quoting(inner):
            return match.group(0)
        return inner

    cleaned = _QUOTED_VALUE_RE.sub(_unquote_if_safe, content)

    if cleaned != content:
        with open(filepath, "w", encoding="utf-8") as fh:
            fh.write(cleaned)


def _parse_seq_id(raw: str):
    """Return (True, int_value) if numeric, else (False, 0)."""
    try:
        return True, int(raw)
    except (ValueError, TypeError):
        return False, 0


def _renumber_seq_id_per_chain(
    chain_ids: List[str], auth_seq_ids: List[str]
) -> List[str]:
    """Produce 1-based label_seq_id values per chain.

    ColabFold multimer PDB files may use continuous or offset auth_seq_id
    numbering across chains (e.g. chain B: -89..320).  mmCIF parent-child
    constraints require _entity_poly_seq.num (derived from label_seq_id)
    to be a positive 1-based index per entity.  This function builds a
    mapping from each chain's sorted unique auth_seq_id values to a
    1-based sequence, then applies that mapping to every atom row.
    Non-numeric values (e.g. "?", ".") are treated as distinct and
    assigned the next ordinal so validators see only positive integers.
    """
    # Per chain: collect unique (raw_seq_id, sort_key) for stable ordering
    # Numeric: sort by int value. Non-numeric: sort after numerics by first-seen index.
    unique_per_chain: Dict[str, Dict[str, int]] = defaultdict(dict)
    first_seen: Dict[str, int] = defaultdict(lambda: 0)
    for chain_id, raw in zip(chain_ids, auth_seq_ids):
        if raw not in unique_per_chain[chain_id]:
            unique_per_chain[chain_id][raw] = first_seen[chain_id]
            first_seen[chain_id] += 1

    renumber_map: Dict[str, Dict[str, str]] = {}
    for chain_id, raw_to_index in unique_per_chain.items():
        numeric_items = []
        non_numeric_items = []
        for raw, first_idx in raw_to_index.items():
            is_num, val = _parse_seq_id(raw)
            if is_num:
                numeric_items.append((val, raw))
            else:
                non_numeric_items.append((first_idx, raw))
        numeric_items.sort(key=lambda x: x[0])
        non_numeric_items.sort(key=lambda x: x[0])
        ordered_raw = [r for (_, r) in numeric_items] + [r for (_, r) in non_numeric_items]
        renumber_map[chain_id] = {
            raw: str(new_idx)
            for new_idx, raw in enumerate(ordered_raw, start=1)
        }

    return [
        renumber_map[chain_id][raw_seq_id]
        for chain_id, raw_seq_id in zip(chain_ids, auth_seq_ids)
    ]


class CifDataStorage:
    """A container for holding and writing mmCIF data."""

    def __init__(self):
        self.data: Dict[str, Dict[str, List[Any]]] = {}

    def set_items(self, category_name: str, items_dict: Dict[str, List[Any]]):
        if category_name not in self.data:
            self.data[category_name] = {}
        for item, values in items_dict.items():
            self.data[category_name][item] = values

    def set_item(self, category_name: str, item_name: str, item_value: Any):
        if category_name not in self.data:
            self.data[category_name] = {}
        self.data[category_name][item_name] = item_value

    def get_data(self) -> Dict[str, Dict[str, List[Any]]]:
        return self.data

    def populate_from_cif_block(self, cif_block: gemmi.cif.Block):
        """Initializes the storage from a gemmi cif.Block."""
        for category in cif_block.get_mmcif_category_names():
            items = cif_block.get_mmcif_category(category)
            self.set_items(category, items)
        # Perform initial data mappings required for consistency
        self.data[CAT_ATOM_SITE][ITEM_LABEL_ASYM_ID] = self.data[CAT_ATOM_SITE][
            ITEM_AUTH_ASYM_ID
        ]
        self.data[CAT_ATOM_SITE][ITEM_AUTH_COMP_ID] = self.data[CAT_ATOM_SITE][
            ITEM_LABEL_COMP_ID
        ]
        renumbered = _renumber_seq_id_per_chain(
            self.data[CAT_ATOM_SITE][ITEM_AUTH_ASYM_ID],
            self.data[CAT_ATOM_SITE][ITEM_AUTH_SEQ_ID],
        )
        self.data[CAT_ATOM_SITE][ITEM_LABEL_SEQ_ID] = renumbered
        self.data[CAT_ATOM_SITE][ITEM_AUTH_SEQ_ID] = renumbered
        if CAT_SYMMETRY in self.data:
            del(self.data[CAT_SYMMETRY])
        if CAT_CELL in self.data:
            del(self.data[CAT_CELL])
    
    def convert_values_to_none(self):
        """Converts all '?' values to None using list comprehension for speed."""
        for category, items in self.data.items():
            for item, values in items.items():
                # List comprehension is ~2-3x faster than index-based in-place modification
                self.data[category][item] = [None if v == "?" else v for v in values]

    def write_to_cif(self, output_file: str, block_name: str = "model", skip_alignment: bool = False):
        """Writes the stored data to an mmCIF file.

        Args:
            output_file: Path to write the CIF file.
            block_name: Name for the CIF block.
            skip_alignment: If True, skip column alignment for faster writes (useful for intermediate files).
        """
        logger.info("Writing CIF file...")
        self.convert_values_to_none()
        doc = gemmi.cif.Document()
        block = doc.add_new_block(block_name)
        for category, items in self.data.items():
            block.set_mmcif_category(category, items)

        # Use pre-cached WriteOptions for efficiency
        write_options = _WRITE_OPTIONS_MINIMAL if skip_alignment else _WRITE_OPTIONS_ALIGNED
        doc.write_file(output_file, write_options)
        _strip_unnecessary_cif_quotes(output_file)
        logger.info(f"mmCIF file written to: {output_file}")