"""Focused tests for the streaming fragment-manifest enricher."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import duckdb

from uniprot.scripts.enrich_fragment_manifest import enrich_manifest, main


FIXTURE_DIRECTORY = Path(__file__).parent / "fixtures" / "fragment_manifest"


def _database(path: Path, rows: list[tuple[str, int, int, str]]) -> None:
    con = duckdb.connect(str(path))
    con.execute(
        """
        CREATE TABLE fragment_metadata (
            uniprot_ac VARCHAR,
            sequence_start INTEGER,
            sequence_end INTEGER,
            protein_name VARCHAR
        )
        """
    )
    if rows:
        con.executemany(
            "INSERT INTO fragment_metadata VALUES (?, ?, ?, ?)", rows
        )
    con.close()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_wide_manifest_normalises_entities_and_enriches_names(
    tmp_path: Path,
) -> None:
    db = tmp_path / "metadata.duckdb"
    _database(
        db,
        [
            ("P27409", 1, 46, "Fragment one"),
            ("P27409", 47, 92, "Fragment two"),
        ],
    )
    source = tmp_path / "wide.csv"
    source.write_text(
        "afdb_id,chain_a,chain_a_uniprot,chain_a_is_polyprotein,"
        "chain_a_start,chain_a_end,chain_b,chain_b_uniprot,"
        "chain_b_is_polyprotein,chain_b_start,chain_b_end,iptm\n"
        "AF-MONOMER,A,P11111,false,1,100,,,,,,0.1\n"
        "AF-HOMO,A,P27409,true,1,46,B,P27409,true,1,46,0.2\n"
        "AF-HETERO,A,P27409,true,1,46,B,P27409,true,47,92,0.3\n",
        encoding="utf-8",
    )
    output = tmp_path / "canonical.csv"

    report = enrich_manifest(source, output, db)

    assert report["schema"] == "wide"
    assert report["counts"]["input_rows"] == 3
    assert report["counts"]["output_rows"] == 5
    rows = _read_rows(output)
    assert [(row["model_entity_id"], row["entity_id"]) for row in rows] == [
        ("AF-MONOMER", "1"),
        ("AF-HOMO", "1"),
        ("AF-HOMO", "1"),
        ("AF-HETERO", "1"),
        ("AF-HETERO", "2"),
    ]
    assert [row["protein_name"] for row in rows] == [
        "",
        "Fragment one",
        "Fragment one",
        "Fragment one",
        "Fragment two",
    ]
    assert [row["iptm"] for row in rows] == [
        "0.1", "0.2", "0.2", "0.3", "0.3"
    ]


def test_actual_virus_wide_header_uses_literal_chain_labels(
    tmp_path: Path,
) -> None:
    db = tmp_path / "metadata.duckdb"
    _database(
        db,
        [
            ("P27409", 1, 46, "NS1"),
            ("P27409", 961, 1071, "Viral genome-linked protein"),
        ],
    )
    output = tmp_path / "canonical.csv"

    report = enrich_manifest(
        FIXTURE_DIRECTORY / "virus_wide_manifest.csv", output, db
    )

    assert report["schema"] == "wide"
    rows = _read_rows(output)
    assert [(row["model_entity_id"], row["chain_id"]) for row in rows] == [
        ("AF-0000000211971324", "A"),
        ("AF-0000000212005744", "A"),
        ("AF-0000000212039401", "A"),
        ("AF-0000000212039401", "B"),
    ]
    assert [row["protein_name"] for row in rows] == [
        "", "NS1", "NS1", "Viral genome-linked protein"
    ]
    assert [row["virus_name"] for row in rows] == [
        "Canis familiaris oral papillomavirus 1",
        "feline calicivirus",
        "feline calicivirus",
        "feline calicivirus",
    ]
    assert "chain_a_id" not in rows[0]
    assert "chain_b_length" not in rows[0]


def test_canonical_manifest_preserves_extra_and_existing_name(
    tmp_path: Path,
) -> None:
    db = tmp_path / "metadata.duckdb"
    _database(db, [("P27409", 1, 46, "Database name")])
    source = tmp_path / "canonical.csv"
    source.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,is_fragment,"
        "sequence_start,sequence_end,protein_name,average_plddt\n"
        "AF-1,99,A,P27409,true,1,46,Provided name,92.3\n"
        "AF-1,98,B,P27409,true,1,46,,91.8\n"
        "AF-1,42,C,P11111,false,1,100,,90.1\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.csv"

    report = enrich_manifest(source, output, db)

    assert report["schema"] == "canonical"
    assert report["counts"]["fragment_names_preserved"] == 1
    rows = _read_rows(output)
    assert [row["entity_id"] for row in rows] == ["1", "1", "2"]
    assert [row["protein_name"] for row in rows] == [
        "Provided name", "Database name", ""
    ]
    assert [row["average_plddt"] for row in rows] == ["92.3", "91.8", "90.1"]


def test_missing_name_reports_and_strict_mode_keeps_destination(
    tmp_path: Path,
) -> None:
    db = tmp_path / "metadata.duckdb"
    _database(db, [])
    source = tmp_path / "canonical.csv"
    source.write_text(
        "model_entity_id,chain_id,uniprot_ac,is_fragment,sequence_start,"
        "sequence_end\n"
        "AF-1,A,P27409,true,1,46\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.csv"
    output.write_text("old output\n", encoding="utf-8")
    report_path = tmp_path / "report.json"

    exit_code = main(
        [
            "--input", str(source), "--output", str(output), "--db", str(db),
            "--strict", "--report", str(report_path),
        ]
    )

    assert exit_code == 1
    assert output.read_text(encoding="utf-8") == "old output\n"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["status"] == "strict_failure"
    assert report["counts"]["unmatched_fragment_lookups"] == 1
    assert report["unmatched_fragment_lookups"] == [
        {"uniprot_ac": "P27409", "sequence_start": 1, "sequence_end": 46}
    ]


def test_ambiguous_fragment_lookup_is_reported_without_blocking(
    tmp_path: Path,
) -> None:
    db = tmp_path / "metadata.duckdb"
    _database(
        db,
        [("P27409", 1, 46, "First name"), ("P27409", 1, 46, "Second name")],
    )
    source = tmp_path / "canonical.csv"
    source.write_text(
        "model_entity_id,chain_id,uniprot_ac,is_fragment,sequence_start,"
        "sequence_end\n"
        "AF-1,A,P27409,true,1,46\n",
        encoding="utf-8",
    )
    output = tmp_path / "output.csv"

    report = enrich_manifest(source, output, db)

    assert report["counts"]["ambiguous_fragment_lookups"] == 1
    assert _read_rows(output)[0]["protein_name"] == ""


def test_streams_representative_wide_rows(tmp_path: Path) -> None:
    db = tmp_path / "metadata.duckdb"
    _database(db, [])
    source = tmp_path / "wide.csv"
    header = (
        "afdb_id,chain_a,chain_a_uniprot,chain_a_is_polyprotein,"
        "chain_a_start,chain_a_end\n"
    )
    with source.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for number in range(5_000):
            handle.write(f"AF-{number},A,P{number:05d},false,1,100\n")
    output = tmp_path / "output.csv"

    report = enrich_manifest(source, output, db)

    assert report["counts"]["input_rows"] == 5_000
    assert report["counts"]["output_rows"] == 5_000
    with output.open(newline="") as handle:
        assert sum(1 for _ in handle) == 5_001
