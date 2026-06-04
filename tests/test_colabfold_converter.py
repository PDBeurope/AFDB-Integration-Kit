import json
import re
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


def test_convert_file_uses_manifest_sequence_overrides_with_duckdb(tmp_path: Path) -> None:
    duckdb_path = tmp_path / "entries.duckdb"
    manifest_path = tmp_path / "manifest.csv"
    scores_json = tmp_path / "scores.json"
    pdb_file = tmp_path / "test.pdb"

    _create_duckdb(
        duckdb_path,
        [
            ("P11111", "Protein one", "ABCDEFGH"),
            ("Q22222", "Protein two", "JKLMNOPQR"),
        ],
    )
    manifest_path.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,sequence_start,sequence_end,is_fragment\n"
        "AF-0000000000000002,1,A,P11111,1,2,true\n"
        "AF-0000000000000002,2,B,Q22222,1,3,true\n",
        encoding="utf-8",
    )
    _write_scores_json(
        scores_json,
        plddt=[90.0, 80.0, 70.0, 60.0, 50.0],
        pae=[
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [1.0, 0.0, 2.0, 3.0, 4.0],
            [2.0, 2.0, 0.0, 1.0, 2.0],
            [3.0, 3.0, 1.0, 0.0, 1.0],
            [4.0, 4.0, 2.0, 1.0, 0.0],
        ],
        max_pae=4.0,
    )
    _write_test_pdb(pdb_file)

    output_paths = convert_file(
        str(scores_json),
        str(pdb_file),
        outdir=str(tmp_path),
        manifest_path=str(manifest_path),
        model_entity_id="AF-0000000000000002",
        duckdb_path=str(duckdb_path),
        out_chain_manifest=str(tmp_path / "chain_manifest.csv"),
        out_model_manifest=str(tmp_path / "model_manifest.csv"),
    )

    with open(output_paths["plddt"], encoding="utf-8") as handle:
        plddt_payload = json.load(handle)

    assert plddt_payload["chains"] == [
        {
            "name": "Protein one",
            "label_asym_id": "A",
            "sequenceStart": 1,
            "sequenceEnd": 2,
        },
        {
            "name": "Protein two",
            "label_asym_id": "B",
            "sequenceStart": 3,
            "sequenceEnd": 5,
        },
    ]
    assert plddt_payload["residueNumber"] == [1, 2, 3, 4, 5]

    chain_manifest_lines = (tmp_path / "chain_manifest.csv").read_text(encoding="utf-8").splitlines()
    assert chain_manifest_lines == [
        "model_entity_id,entity_id,chain_id,uniprot_ac,is_fragment,is_isoform,entity_type,sequence_start,sequence_end,average_plddt,fraction_plddt_very_low,fraction_plddt_low,fraction_plddt_confident,fraction_plddt_very_high",
        "AF-0000000000000002,1,A,P11111,true,,protein,1,2,85.0,0.0,0.0,1.0,0.0",
        "AF-0000000000000002,2,B,Q22222,true,,protein,1,3,60.0,0.0,0.667,0.333,0.0",
    ]


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


REAL_EXAMPLES_DIR = Path(__file__).parent / "fixtures" / "colabfold_real_examples"
AF_ID_RE = re.compile(r"^AF-\d{16}$")


def _load_real_examples() -> list[dict]:
    with open(REAL_EXAMPLES_DIR / "manifest.json", encoding="utf-8") as handle:
        return json.load(handle)["examples"]


def _fixture_path(example: dict, suffix: str) -> Path:
    category_dir = REAL_EXAMPLES_DIR / f"{example['category']}s" / example["example_id"]
    matches = [f["name"] for f in example["files"] if f["name"].endswith(suffix)]
    assert len(matches) == 1
    return category_dir / matches[0]


def test_real_colabfold_fixture_names_are_single_af_ids() -> None:
    for example in _load_real_examples():
        example_id = example["example_id"]
        assert AF_ID_RE.match(example_id)
        for file_info in example["files"]:
            assert file_info["name"].startswith(f"{example_id}-")
        for chain in example["chain_spans"]:
            assert chain["uniprot_ac"]


@pytest.mark.parametrize("example", _load_real_examples(), ids=lambda ex: ex["example_id"])
def test_convert_file_handles_curated_real_colabfold_examples(
    tmp_path: Path,
    example: dict,
) -> None:
    scores_json = _fixture_path(example, "-scores_v1.json")
    pdb_file = _fixture_path(example, "-model_v1.pdb")

    output_paths = convert_file(
        str(scores_json),
        str(pdb_file),
        outdir=str(tmp_path / example["example_id"]),
        model_entity_id=example["example_id"],
    )

    with open(scores_json, encoding="utf-8") as handle:
        raw_scores = json.load(handle)
    with open(output_paths["plddt"], encoding="utf-8") as handle:
        plddt_payload = json.load(handle)
    with open(output_paths["pae"], encoding="utf-8") as handle:
        pae_payload = json.load(handle)[0]

    expected_chains = [
        {
            "name": chain["chain_id"],
            "label_asym_id": chain["chain_id"],
            "sequenceStart": chain["sequence_start"],
            "sequenceEnd": chain["sequence_end"],
        }
        for chain in example["chain_spans"]
    ]

    assert plddt_payload["chains"] == expected_chains
    assert pae_payload["chains"] == expected_chains
    assert len(plddt_payload["confidenceScore"]) == example["plddt_length"]
    assert plddt_payload["residueNumber"] == list(range(1, example["plddt_length"] + 1))
    assert len(pae_payload["predicted_aligned_error"]) == example["pae_dimension"]
    assert all(
        len(row) == example["pae_dimension"]
        for row in pae_payload["predicted_aligned_error"]
    )
    assert pae_payload["max_predicted_aligned_error"] == round(raw_scores["max_pae"], 2)
