import json
from pathlib import Path

import pytest

from afdb_integration_kit.colabfold.converter import (
    _chain_spans_from_pdb_gemmi,
    _chain_spans_from_pdb_legacy,
    _load_chain_metadata_from_duckdb,
    cleanup_caches,
    convert_file,
    prefetch_duckdb_metadata,
)


def _write_scores_json(path: Path, *, plddt: list[float], pae: list[list[float]], max_pae: float) -> None:
    path.write_text(
        json.dumps({"plddt": plddt, "pae": pae, "max_pae": max_pae}),
        encoding="utf-8",
    )


def _write_test_pdb(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "ATOM      1  N   ALA A   1      11.104  13.207  14.399  1.00 20.00           N",
                "ATOM      2  CA  ALA A   1      11.500  12.000  15.000  1.00 20.00           C",
                "ATOM      3  N   GLY A   2      12.000  11.000  15.500  1.00 20.00           N",
                "ATOM      4  N   SER B  10      13.000  10.000  16.000  1.00 20.00           N",
                "ATOM      5  N   THR B  10A     14.000   9.000  16.500  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _create_duckdb(path: Path, rows: list[tuple[str, str, str]]) -> None:
    duckdb = pytest.importorskip("duckdb")
    con = duckdb.connect(str(path))
    try:
        con.execute(
            """
            CREATE TABLE entry (
                primary_ac VARCHAR,
                protein_full_names VARCHAR,
                sequence VARCHAR
            )
            """
        )
        con.executemany("INSERT INTO entry VALUES (?, ?, ?)", rows)
    finally:
        con.close()


def test_gemmi_chain_spans_match_legacy_parser(tmp_path: Path) -> None:
    pdb_path = tmp_path / "test.pdb"
    _write_test_pdb(pdb_path)

    legacy, legacy_total = _chain_spans_from_pdb_legacy(pdb_path)
    gemmi, gemmi_total = _chain_spans_from_pdb_gemmi(pdb_path)

    expected = [
        {"name": "A", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 2},
        {"name": "B", "label_asym_id": "B", "sequenceStart": 3, "sequenceEnd": 4},
    ]

    assert legacy == expected
    assert gemmi == expected
    assert legacy_total == gemmi_total == 4


def test_convert_file_rounds_outputs_and_preserves_global_chain_spans(tmp_path: Path) -> None:
    scores_json = tmp_path / "scores.json"
    pdb_file = tmp_path / "test.pdb"
    _write_scores_json(
        scores_json,
        plddt=[95.123, 89.999, 49.994, 29.995],
        pae=[
            [0.111, 1.9994, 2.555, 3.3333],
            [1.4444, 0.222, 4.444, 5.555],
            [2.777, 4.666, 0.3333, 6.666],
            [3.888, 5.999, 6.111, 0.4444],
        ],
        max_pae=6.666,
    )
    _write_test_pdb(pdb_file)

    output_paths = convert_file(str(scores_json), str(pdb_file), outdir=str(tmp_path))

    with open(output_paths["plddt"], encoding="utf-8") as handle:
        plddt_payload = json.load(handle)
    with open(output_paths["pae"], encoding="utf-8") as handle:
        pae_payload = json.load(handle)[0]

    assert plddt_payload["confidenceScore"] == [95.12, 90.0, 49.99, 30.0]
    assert plddt_payload["confidenceCategory"] == ["V", "H", "L", "L"]
    assert plddt_payload["residueNumber"] == [1, 2, 3, 4]
    assert plddt_payload["chains"] == [
        {"name": "A", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 2},
        {"name": "B", "label_asym_id": "B", "sequenceStart": 3, "sequenceEnd": 4},
    ]
    assert pae_payload["predicted_aligned_error"] == [
        [0.11, 2.0, 2.56, 3.33],
        [1.44, 0.22, 4.44, 5.56],
        [2.78, 4.67, 0.33, 6.67],
        [3.89, 6.0, 6.11, 0.44],
    ]
    assert pae_payload["max_predicted_aligned_error"] == 6.67


def test_convert_file_auto_expands_homomultimer_manifest_with_duckdb(tmp_path: Path) -> None:
    duckdb_path = tmp_path / "entries.duckdb"
    manifest_path = tmp_path / "manifest.csv"
    scores_json = tmp_path / "scores.json"
    pdb_file = tmp_path / "test.pdb"

    _create_duckdb(duckdb_path, [("P11111", "Example protein", "AC")])
    manifest_path.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac\n"
        "AF-0000000000000001,1,A,P11111\n",
        encoding="utf-8",
    )
    _write_scores_json(
        scores_json,
        plddt=[90.0, 80.0, 70.0, 60.0],
        pae=[
            [0.0, 1.0, 2.0, 3.0],
            [1.0, 0.0, 2.0, 3.0],
            [2.0, 2.0, 0.0, 1.0],
            [3.0, 3.0, 1.0, 0.0],
        ],
        max_pae=3.0,
    )
    _write_test_pdb(pdb_file)

    output_paths = convert_file(
        str(scores_json),
        str(pdb_file),
        outdir=str(tmp_path),
        manifest_path=str(manifest_path),
        model_entity_id="AF-0000000000000001",
        duckdb_path=str(duckdb_path),
    )

    with open(output_paths["plddt"], encoding="utf-8") as handle:
        plddt_payload = json.load(handle)

    assert plddt_payload["chains"] == [
        {
            "name": "Example protein",
            "label_asym_id": "A",
            "sequenceStart": 1,
            "sequenceEnd": 2,
        },
        {
            "name": "Example protein",
            "label_asym_id": "B",
            "sequenceStart": 3,
            "sequenceEnd": 4,
        },
    ]
    assert plddt_payload["residueNumber"] == [1, 2, 3, 4]


def test_partial_prefetch_does_not_break_later_duckdb_lookups(tmp_path: Path) -> None:
    duckdb_path = tmp_path / "entries.duckdb"
    _create_duckdb(
        duckdb_path,
        [
            ("P11111", "Protein one", "AC"),
            ("Q22222", "Protein two", "GT"),
        ],
    )

    cleanup_caches()
    prefetch_duckdb_metadata(str(duckdb_path), ["P11111"])

    chains, residue_numbers, effective_manifest = _load_chain_metadata_from_duckdb(
        duckdb_path,
        manifest_chains=[{"chain_id": "B", "uniprot_ac": "Q22222", "entity_id": "2"}],
    )

    assert chains == [
        {
            "name": "Protein two",
            "label_asym_id": "B",
            "sequenceStart": 1,
            "sequenceEnd": 2,
        }
    ]
    assert residue_numbers == [1, 2]
    assert effective_manifest == [{"chain_id": "B", "uniprot_ac": "Q22222", "entity_id": "2"}]

    cleanup_caches()
