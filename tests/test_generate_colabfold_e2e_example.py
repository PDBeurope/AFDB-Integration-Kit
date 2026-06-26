from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_helper_module():
    repo_root = Path(__file__).resolve().parent.parent
    module_path = repo_root / "scripts" / "generate_colabfold_e2e_example.py"
    spec = importlib.util.spec_from_file_location("generate_colabfold_e2e_example", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_chain_manifest_rows_use_local_chain_ranges() -> None:
    module = _load_helper_module()

    models = [
        module.ExampleModel(
            category="homodimer",
            example_id="AF-0000000066074510",
            directory=Path("/tmp/homodimer"),
            chain_spans=[
                {
                    "chain_id": "A",
                    "sequence_start": 1,
                    "sequence_end": 461,
                    "residue_count": 461,
                    "uniprot_ac": "Q46806",
                },
                {
                    "chain_id": "B",
                    "sequence_start": 462,
                    "sequence_end": 922,
                    "residue_count": 461,
                    "uniprot_ac": "Q46806",
                },
            ],
        ),
        module.ExampleModel(
            category="heterodimer",
            example_id="AF-0000000300000101",
            directory=Path("/tmp/heterodimer"),
            chain_spans=[
                {
                    "chain_id": "A",
                    "sequence_start": 1,
                    "sequence_end": 30,
                    "residue_count": 30,
                    "uniprot_ac": "A0ABS2QMZ4",
                },
                {
                    "chain_id": "B",
                    "sequence_start": 31,
                    "sequence_end": 63,
                    "residue_count": 33,
                    "uniprot_ac": "A0ABS2QMF5",
                },
            ],
        ),
    ]

    assert module.chain_manifest_rows(models) == [
        {
            "model_entity_id": "AF-0000000066074510",
            "entity_id": "1",
            "chain_id": "A",
            "uniprot_ac": "Q46806",
            "sequence_start": 1,
            "sequence_end": 461,
        },
        {
            "model_entity_id": "AF-0000000066074510",
            "entity_id": "1",
            "chain_id": "B",
            "uniprot_ac": "Q46806",
            "sequence_start": 1,
            "sequence_end": 461,
        },
        {
            "model_entity_id": "AF-0000000300000101",
            "entity_id": "1",
            "chain_id": "A",
            "uniprot_ac": "A0ABS2QMZ4",
            "sequence_start": 1,
            "sequence_end": 30,
        },
        {
            "model_entity_id": "AF-0000000300000101",
            "entity_id": "2",
            "chain_id": "B",
            "uniprot_ac": "A0ABS2QMF5",
            "sequence_start": 1,
            "sequence_end": 33,
        },
    ]


def test_ipsae_enrichment_maps_model_and_chain_metrics(tmp_path) -> None:
    module = _load_helper_module()

    ipsae_csv = tmp_path / "ipsae_summary.csv"
    ipsae_csv.write_text(
        "pdb_path,pae_cutoff,dist_cutoff,iptm_af,ipsae_AB,ipsae_BA,"
        "pDockQ2_AB,pDockQ2_BA,LIS_AB,LIS_BA,pDockQ,n0chn,processing_time_ms\n"
        "/tmp/AF-TEST-model_v1.pdb,10,15,0.82,0.69,0.58,0.33,0.22,0.41,0.31,0.51,734,7\n",
        encoding="utf-8",
    )

    model_dir = tmp_path / "model_jsons"
    chain_dir = tmp_path / "chain_jsons"
    model_dir.mkdir()
    chain_dir.mkdir()
    module.json_dump(model_dir / "AF-TEST.json", {"entryId": "AF-TEST"})
    module.json_dump(
        chain_dir / "AF-TEST.json",
        [
            {"uniqueId": "AF-TEST_A"},
            {"uniqueId": "AF-TEST_B"},
        ],
    )

    assert module.enrich_model_jsons(model_dir, ipsae_csv) == 1
    assert module.enrich_chain_jsons(chain_dir, ipsae_csv) == 1

    model_payload = module.json.loads((model_dir / "AF-TEST.json").read_text())
    assert model_payload["complexPredictionAccuracy_iptm_af"] == 0.82
    assert model_payload["complexPredictionAccuracy_ipsae_AB"] == 0.69
    assert model_payload["complexPredictionAccuracy_ipsae_BA"] == 0.58
    assert model_payload["complexPredictionAccuracy_pDockQ2_AB"] == 0.33
    assert model_payload["complexPredictionAccuracy_pDockQ"] == 0.51
    assert model_payload["complexPredictionAccuracy_ipsae_n0chn"] == 734.0

    chain_payload = module.json.loads((chain_dir / "AF-TEST.json").read_text())
    chain_a_keys = set(chain_payload[0])
    chain_b_keys = set(chain_payload[1])
    assert "complexPredictionAccuracy_ipsae_AB" in chain_a_keys
    assert "complexPredictionAccuracy_ipsae_BA" not in chain_a_keys
    assert "complexPredictionAccuracy_pDockQ2_AB" in chain_a_keys
    assert "complexPredictionAccuracy_pDockQ2_BA" not in chain_a_keys
    assert "complexPredictionAccuracy_ipsae_BA" in chain_b_keys
    assert "complexPredictionAccuracy_ipsae_AB" not in chain_b_keys
    assert "complexPredictionAccuracy_pDockQ2_BA" in chain_b_keys
    assert "complexPredictionAccuracy_pDockQ2_AB" not in chain_b_keys
    assert chain_payload[0]["complexPredictionAccuracy_pDockQ"] == 0.51
    assert chain_payload[1]["complexPredictionAccuracy_pDockQ"] == 0.51
