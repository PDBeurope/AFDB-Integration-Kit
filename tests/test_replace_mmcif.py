import json
import pytest
from unittest.mock import MagicMock, patch
from afdb_integration_kit.modelcif_replace.replace import replace_mmcif_with_json


@patch("gemmi.cif")
def test_replace_success(mock_cif, tmp_path):
    # Mock mmCIF document and block
    mock_doc = MagicMock()
    mock_block = MagicMock()
    mock_cif.read_file.return_value = mock_doc
    mock_doc.sole_block.return_value = mock_block

    # Mock category object
    mock_category = MagicMock()
    mock_block.get_mmcif_category.return_value = mock_category

    # Create a fake JSON file
    json_data = {"categories": {"test_cat": {"key1": "val1"}}}
    json_file = tmp_path / "input.json"
    mmcif_file = tmp_path / "input.cif"
    output_file = tmp_path / "output.cif"

    json_file.write_text(json.dumps(json_data))
    mmcif_file.write_text("dummy content")

    replace_mmcif_with_json(str(mmcif_file), str(json_file), str(output_file))

    mock_block.get_mmcif_category.assert_called_once_with("test_cat")
    mock_category.update.assert_called_once_with({"key1": "val1"})
    mock_block.set_mmcif_category.assert_called_once()
    mock_doc.write_file.assert_called_once()


@patch("gemmi.cif")
def test_invalid_json_structure(mock_cif, tmp_path):
    # Mock mmCIF document and block so we get past file parsing
    mock_doc = MagicMock()
    mock_block = MagicMock()
    mock_cif.read_file.return_value = mock_doc
    mock_doc.sole_block.return_value = mock_block

    # JSON missing 'categories'
    json_data = {"invalid": {}}
    json_file = tmp_path / "input.json"
    mmcif_file = tmp_path / "input.cif"
    output_file = tmp_path / "output.cif"

    json_file.write_text(json.dumps(json_data))
    mmcif_file.write_text("dummy content")

    with pytest.raises(ValueError, match="categories"):
        replace_mmcif_with_json(str(mmcif_file), str(json_file), str(output_file))


@patch("gemmi.cif")
def test_category_update_failure(mock_cif, tmp_path, caplog):
    # Mock mmCIF document and block
    mock_doc = MagicMock()
    mock_block = MagicMock()
    mock_cif.read_file.return_value = mock_doc
    mock_doc.sole_block.return_value = mock_block

    # Raise error when accessing category
    mock_block.get_mmcif_category.side_effect = Exception("Category not found")

    # Create a fake JSON file
    json_data = {"categories": {"missing_cat": {"key1": "val1"}}}
    json_file = tmp_path / "input.json"
    mmcif_file = tmp_path / "input.cif"
    output_file = tmp_path / "output.cif"

    json_file.write_text(json.dumps(json_data))
    mmcif_file.write_text("dummy content")

    replace_mmcif_with_json(str(mmcif_file), str(json_file), str(output_file))

    assert "Failed to replace category" in caplog.text
