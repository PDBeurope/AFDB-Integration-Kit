from __future__ import annotations

from pathlib import Path

from afdb_integration_kit.modelpdb.generate import generate_pdb_headers


def test_generate_pdb_headers_preserves_heterodimer_entities(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    cif_file = repo_root / "examples/colabfold_complex_e2e/dssp/AF-0000000300000101-model_v1.cif"
    pdb_file = repo_root / "examples/colabfold_complex_e2e/input/AF-0000000300000101-model_v1.pdb"
    provider_json = repo_root / "examples/colabfold_complex_e2e/config/provider.json"
    output_pdb = tmp_path / "heterodimer-model_v1.pdb"

    generate_pdb_headers(str(cif_file), str(pdb_file), str(output_pdb), str(provider_json))
    contents = output_pdb.read_text(encoding="utf-8")

    assert "TITLE     COMPLEX OF CYTOCHROME P450/LARGE SUBUNIT RIBOSOMAL PROTEIN L3" in contents
    assert "COMPND    MOL_ID: 1" in contents
    assert "COMPND   3 CHAIN: A" in contents
    assert "COMPND   2 MOLECULE: LARGE SUBUNIT RIBOSOMAL PROTEIN L3;" in contents
    assert "COMPND   3 CHAIN: B" in contents
