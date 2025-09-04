from typing import Any, Dict, List
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
        self.data[CAT_ATOM_SITE][ITEM_LABEL_SEQ_ID] = self.data[CAT_ATOM_SITE][
            ITEM_AUTH_SEQ_ID
        ]
        del(self.data[CAT_SYMMETRY])
        del(self.data[CAT_CELL])

    def write_to_cif(self, output_file: str, block_name: str = "model"):
        """Writes the stored data to an mmCIF file."""
        logger.info("Writing CIF file...")
        doc = gemmi.cif.Document()
        block = doc.add_new_block(block_name)
        for category, items in self.data.items():
            block.set_mmcif_category(category, items)

        write_options = gemmi.cif.WriteOptions()
        write_options.align_loops = 50
        write_options.align_pairs = 50
        write_options.prefer_pairs = True
        doc.write_file(output_file, write_options)
        logger.info(f"mmCIF file written to: {output_file}")