import json
from pathlib import Path

import pytest

from afdb_integration_kit.colabfold.converter import convert_file


@pytest.mark.skip(reason="Test fixture files not present in repository")
def test_convert_file_adds_chain_metadata(tmp_path: Path) -> None:
    scores_json = Path("examples/multimer_examples/test_fffe7/test_fffe7_scores_rank_001_alphafold2_multimer_v3_model_2_seed_000.json")
    pdb_file = Path("examples/multimer_examples/test_fffe7/test_fffe7_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb")

    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "model_entity_id,chain_id,uniprot_ac\n"
        "TEST_MODEL,A,P35354\n"
        "TEST_MODEL,B,P35354\n"
    )

    output_paths = convert_file(
        str(scores_json),
        str(pdb_file),
        outdir=tmp_path,
        manifest_path=str(manifest),
        model_entity_id="TEST_MODEL",
    )

    with open(output_paths["plddt"]) as handle:
        plddt_payload = json.load(handle)
    with open(output_paths["pae"]) as handle:
        pae_payload = json.load(handle)[0]

    expected_chains = [
        {"name": "P35354", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 604},
        {"name": "P35354", "label_asym_id": "B", "sequenceStart": 605, "sequenceEnd": 1208},
    ]

    assert plddt_payload["chains"] == expected_chains
    assert pae_payload["chains"] == expected_chains
    assert plddt_payload["residueNumber"][0] == 1
    assert plddt_payload["residueNumber"][-1] == 1208
