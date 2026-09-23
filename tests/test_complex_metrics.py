from pathlib import Path

from afdb_integration_kit.complex_metrics import (
    build_chain_enrichment,
    build_model_enrichment,
    parse_ipsae_csv,
)


def test_parse_ipsae_csv_supports_homodimer_names_and_coerces_numbers(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "ipsae_summary.csv"
    csv_path.write_text(
        "pdb_path,processing_time_ms,ipsae_AB,comment\n"
        "/work/AF-0000000066074510-model_v1.pdb,12,0.625,kept\n",
        encoding="utf-8",
    )

    assert parse_ipsae_csv(csv_path) == {
        "AF-0000000066074510": {
            "ipsae_AB": 0.625,
            "comment": "kept",
        }
    }


def test_parse_ipsae_csv_returns_empty_mapping_when_file_is_missing(
    tmp_path: Path,
) -> None:
    assert parse_ipsae_csv(tmp_path / "missing.csv") == {}


def test_build_model_enrichment_preserves_json_keys_and_clash_metrics() -> None:
    result = build_model_enrichment(
        {
            "ipsae_AB": 0.61,
            "iptm_af": 0.72,
            "pae_cutoff": 10.0,
            "ignored": 99.0,
        },
        {"N_clash_backbone": 2, "N_clash_heavyAtom": 7},
        [
            "ipsae",
            "iptm_af",
            "pae_cutoff",
            "N_clash_backbone",
            "N_clash_heavyAtom",
        ],
    )

    assert result == {
        "complexPredictionAccuracy_ipsae_AB": 0.61,
        "complexPredictionAccuracy_iptm_af": 0.72,
        "complexPredictionAccuracy_ipsae_pae_cutoff": 10.0,
        "complexPredictionAccuracy_N_clash_backbone": 2,
        "complexPredictionAccuracy_N_clash_heavyAtom": 7,
    }


def test_build_chain_enrichment_routes_directional_metrics_by_source_chain() -> None:
    row = {
        "ipsae_AB": 0.61,
        "ipsae_BA": 0.42,
        "iptm_af": 0.72,
    }
    metrics = ["ipsae", "iptm_af"]

    assert build_chain_enrichment(row, "A", metrics) == {
        "complexPredictionAccuracy_ipsae_AB": 0.61,
        "complexPredictionAccuracy_iptm_af": 0.72,
    }
    assert build_chain_enrichment(row, "B", metrics) == {
        "complexPredictionAccuracy_ipsae_BA": 0.42,
        "complexPredictionAccuracy_iptm_af": 0.72,
    }
