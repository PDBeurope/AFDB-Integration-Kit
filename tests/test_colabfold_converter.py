import csv
import json
from pathlib import Path

import duckdb
import pytest

from afdb_integration_kit.colabfold import converter
from afdb_integration_kit.colabfold.converter import (
    _bind_chains_to_structure,
    _compute_plddt_metrics,
    _load_manifest_chains,
    convert_file,
)


_AA1TO3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def _cand(acc, name, seq, entity_id="1"):
    return {"uniprot_ac": acc, "entity_id": entity_id, "name": name, "seq": seq}


def _phys(label, seq):
    return {"label": label, "length": len(seq), "seq": seq}


def _write_scores_json(path, *, plddt, pae, max_pae):
    path.write_text(
        json.dumps({"plddt": plddt, "pae": pae, "max_pae": max_pae}),
        encoding="utf-8",
    )


def _create_duckdb(path, entries):
    """Create a minimal DuckDB ``entry`` table for the converter slow path."""
    con = duckdb.connect(str(path))
    try:
        con.execute(
            "CREATE TABLE entry(primary_ac VARCHAR, protein_full_names VARCHAR[], sequence VARCHAR)"
        )
        for primary_ac, protein_name, sequence in entries:
            con.execute(
                "INSERT INTO entry VALUES (?, ?, ?)",
                [primary_ac, [protein_name], sequence],
            )
    finally:
        con.close()


def _atom_line(serial, resname, chain, resseq, x):
    # Minimal CA-only ATOM record with standard PDB column positions.
    return (
        f"ATOM  {serial:>5} "      # cols 1-12
        f" CA "                      # 13-16 atom name
        f" "                         # 17 altLoc
        f"{resname:>3} "             # 18-21 resName + space
        f"{chain}"                   # 22 chainID
        f"{resseq:>4} "              # 23-27 resSeq + iCode
        f"   "                       # 28-30
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}"  # 31-54 x,y,z
        f"  1.00 90.00           C  "
    )


def _write_pdb(path, chain_specs):
    """chain_specs: list of (chain_id, resname, n_residues)."""
    lines, serial = [], 1
    for chain_id, resname, n in chain_specs:
        for r in range(1, n + 1):
            lines.append(_atom_line(serial, resname, chain_id, r, float(serial)))
            serial += 1
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def _write_test_pdb(path, chains):
    """chains: list of (chain_id, one-letter sequence)."""
    lines, serial = [], 1
    for chain_id, sequence in chains:
        for r, code in enumerate(sequence, start=1):
            resname = _AA1TO3.get(code, "ALA")
            lines.append(_atom_line(serial, resname, chain_id, r, float(serial)))
            serial += 1
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize(
    ("second_range", "provided_ids", "expected_ids"),
    [
        ((961, 1071), ("1", "2"), ["1", "2"]),
        ((961, 1071), ("", ""), ["1", "2"]),
        ((1, 46), ("", ""), ["1", "1"]),
    ],
)
def test_manifest_entity_assignment_uses_fragment_component_identity(
    tmp_path: Path,
    second_range: tuple[int, int],
    provided_ids: tuple[str, str],
    expected_ids: list[str],
) -> None:
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,is_fragment,"
        "sequence_start,sequence_end\n"
        f"AF-TEST,{provided_ids[0]},A,P27409,true,1,46\n"
        f"AF-TEST,{provided_ids[1]},B,P27409,true,"
        f"{second_range[0]},{second_range[1]}\n",
        encoding="utf-8",
    )

    _, rows = _load_manifest_chains(manifest, "AF-TEST")

    assert [row["entity_id"] for row in rows] == expected_ids


def test_bind_chains_corrects_length_swap():
    # Structure order: A=short (94), B=long (117). Manifest/candidate order is
    # REVERSED (the bug condition). Result must follow the structure.
    physical = [_phys("A", "A" * 94), _phys("B", "G" * 117)]
    candidates = [_cand("ACC_LONG", "Long protein", "G" * 117, "1"),
                  _cand("ACC_SHORT", "Short protein", "A" * 94, "2")]
    chains, residue_numbers, effective = _bind_chains_to_structure(physical, candidates)

    assert [(c["label_asym_id"], c["name"], c["sequenceEnd"]) for c in chains] == [
        ("A", "Short protein", 94),
        ("B", "Long protein", 117),
    ]
    assert all(c["sequenceStart"] == 1 for c in chains)
    assert len(residue_numbers) == 94 + 117
    assert residue_numbers[0] == 1 and residue_numbers[-1] == 211
    # effective_chains must map each structure label to the correct accession.
    assert {(e["chain_id"], e["uniprot_ac"]) for e in effective} == {
        ("A", "ACC_SHORT"), ("B", "ACC_LONG")
    }


def test_bind_chains_disambiguates_equal_length_by_sequence():
    # Both chains 100 aa: length cannot disambiguate, sequence must.
    seq_a = "A" * 100
    seq_b = "G" * 100
    physical = [_phys("A", seq_a), _phys("B", seq_b)]
    candidates = [_cand("ACC_B", "B protein", seq_b), _cand("ACC_A", "A protein", seq_a)]
    chains, _, effective = _bind_chains_to_structure(physical, candidates)
    assert [(c["label_asym_id"], c["name"]) for c in chains] == [
        ("A", "A protein"), ("B", "B protein")
    ]
    assert {(e["chain_id"], e["uniprot_ac"]) for e in effective} == {
        ("A", "ACC_A"), ("B", "ACC_B")
    }


def test_bind_chains_homodimer_keeps_structure_order():
    seq = "M" * 50
    physical = [_phys("A", seq), _phys("B", seq)]
    candidates = [_cand("ACC", "Same protein", seq), _cand("ACC", "Same protein", seq)]
    chains, _, effective = _bind_chains_to_structure(physical, candidates)
    assert [c["label_asym_id"] for c in chains] == ["A", "B"]
    assert all(c["name"] == "Same protein" for c in chains)
    assert all(e["uniprot_ac"] == "ACC" for e in effective)


def test_bind_chains_preserves_unresolved_uniprot_residues():
    """ColabFold can score an ``X`` residue that has no coordinate record.

    Chain matching must tolerate that absent coordinate while retaining the
    full UniProt/pLDDT span in the confidence metadata.
    """
    physical = [_phys("A", "AAA"), _phys("B", "AAA")]
    candidates = [
        _cand("ACC", "Same protein", "XAAA"),
        _cand("ACC", "Same protein", "XAAA"),
    ]

    chains, residue_numbers, effective = _bind_chains_to_structure(physical, candidates)

    assert [c["sequenceEnd"] for c in chains] == [4, 4]
    assert len(residue_numbers) == 8
    assert all(e["uniprot_ac"] == "ACC" for e in effective)


def test_bind_chains_propagates_manifest_provenance():
    """The bound effective rows must retain per-chain provenance fields."""
    physical = [_phys("A", "AAA"), _phys("B", "GGG")]
    candidates = [
        {
            "uniprot_ac": "ACC_A",
            "entity_id": "1",
            "name": "A protein",
            "seq": "AAA",
            "is_fragment": "true",
            "is_isoform": "false",
            "entity_type": "polypeptide",
            "sequence_start": "10",
            "sequence_end": "12",
            "protein_name": "Named A",
        },
        {
            "uniprot_ac": "ACC_B",
            "entity_id": "2",
            "name": "B protein",
            "seq": "GGG",
            "is_fragment": "false",
            "is_isoform": "true",
            "entity_type": "protein",
            "sequence_start": "",
            "sequence_end": "",
        },
    ]

    _, _, effective = _bind_chains_to_structure(physical, candidates)

    by_chain = {e["chain_id"]: e for e in effective}
    assert by_chain["A"]["is_fragment"] == "true"
    assert by_chain["A"]["sequence_start"] == "10"
    assert by_chain["A"]["sequence_end"] == "12"
    assert by_chain["A"]["protein_name"] == "Named A"
    assert by_chain["B"]["is_isoform"] == "true"
    assert "protein_name" not in by_chain["B"]


@pytest.mark.skip(reason="Test fixture files not present in repository")
def test_convert_file_adds_chain_metadata(tmp_path: Path) -> None:
    # Fixture files (not present in the repository):
    # examples/multimer_examples/test_fffe7/test_fffe7_scores_rank_001_alphafold2_multimer_v3_model_2_seed_000.json
    # examples/multimer_examples/test_fffe7/test_fffe7_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb
    pass


def test_convert_file_preserves_manifest_protein_name_in_generated_manifest(
    tmp_path: Path,
) -> None:
    duckdb_path = tmp_path / "entries.duckdb"
    manifest_path = tmp_path / "manifest.csv"
    scores_json = tmp_path / "scores.json"
    pdb_file = tmp_path / "test.pdb"
    chain_manifest = tmp_path / "chain_manifest.csv"

    _create_duckdb(duckdb_path, [("P11111", "Whole protein", "AC")])
    manifest_path.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,protein_name\n"
        "AF-0000000000000003,1,A,P11111,Named fragment\n",
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
    _write_test_pdb(pdb_file, [("A", "AC"), ("B", "AC")])

    output_paths = convert_file(
        str(scores_json),
        str(pdb_file),
        outdir=str(tmp_path),
        manifest_path=str(manifest_path),
        model_entity_id="AF-0000000000000003",
        duckdb_path=str(duckdb_path),
        out_chain_manifest=str(chain_manifest),
    )

    plddt_payload = json.loads(
        Path(output_paths["plddt"]).read_text(encoding="utf-8")
    )
    assert [chain["name"] for chain in plddt_payload["chains"]] == [
        "Named fragment",
        "Named fragment",
    ]
    header = chain_manifest.read_text(encoding="utf-8").splitlines()[0]
    assert "protein_name" in header
    rows = list(
        csv.DictReader(chain_manifest.read_text(encoding="utf-8").splitlines())
    )
    assert all(row["protein_name"] == "Named fragment" for row in rows)


def test_load_chain_metadata_uses_protein_description_for_chain_name(
    tmp_path: Path,
) -> None:
    """The resolved ``protein_description`` must be the emitted chain name."""
    pdb_file = tmp_path / "model.pdb"
    _write_test_pdb(pdb_file, [("A", "AC")])

    duckdb_path = tmp_path / "fake.duckdb"
    duckdb_path.write_text("")
    cache_key = str(duckdb_path.resolve())
    converter._DUCKDB_METADATA_CACHE[cache_key] = {
        "P11111": {
            "primary_ac": "P11111",
            "protein_full_names": ["DuckDB name"],
            "protein_short_names": ["Short name"],
            "sequence": "AC",
        },
    }
    try:
        # A manifest protein_name takes precedence through protein_description.
        chains, _, _ = converter._load_chain_metadata_from_duckdb(
            duckdb_path,
            manifest_chains=[
                {
                    "chain_id": "A",
                    "entity_id": "1",
                    "uniprot_ac": "P11111",
                    "protein_name": "Manifest name",
                }
            ],
            pdb_path=pdb_file,
        )
        assert chains[0]["name"] == "Manifest name"

        # Without a manifest name, fall back to DuckDB protein_full_names.
        chains, _, _ = converter._load_chain_metadata_from_duckdb(
            duckdb_path,
            manifest_chains=[
                {"chain_id": "A", "entity_id": "1", "uniprot_ac": "P11111"}
            ],
            pdb_path=pdb_file,
        )
        assert chains[0]["name"] == "DuckDB name"
    finally:
        converter._DUCKDB_METADATA_CACHE.pop(cache_key, None)


def test_compute_plddt_metrics_indexes_manifest_by_chain_id() -> None:
    """Manifest rows must be looked up by chain ID, not manifest order."""
    chains = [
        {"name": "Alanine", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 3},
        {"name": "Glycine", "label_asym_id": "B", "sequenceStart": 1, "sequenceEnd": 5},
    ]
    # Manifest order is reversed relative to the emitted chain order.
    manifest = [
        {
            "chain_id": "B",
            "uniprot_ac": "ACC_B",
            "entity_id": "2",
            "is_fragment": "true",
            "sequence_start": "10",
            "sequence_end": "14",
        },
        {
            "chain_id": "A",
            "uniprot_ac": "ACC_A",
            "entity_id": "1",
            "is_fragment": "false",
        },
    ]

    rows, _ = _compute_plddt_metrics([90.0] * 8, chains, manifest)

    by_chain = {r["chain_id"]: r for r in rows}
    assert by_chain["A"]["uniprot_ac"] == "ACC_A"
    assert by_chain["A"]["entity_id"] == "1"
    assert by_chain["B"]["uniprot_ac"] == "ACC_B"
    assert by_chain["B"]["sequence_start"] == 10
    assert by_chain["B"]["sequence_end"] == 14


def test_compute_plddt_metrics_falls_back_when_ranges_are_missing() -> None:
    chains = [
        {"name": "Alanine", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 3},
    ]
    manifest = [{"chain_id": "A", "uniprot_ac": "ACC_A", "entity_id": "1"}]

    rows, _ = _compute_plddt_metrics([90.0, 80.0, 70.0], chains, manifest)

    assert rows[0]["sequence_start"] == 1
    assert rows[0]["sequence_end"] == 3


def test_compute_plddt_metrics_uses_fragment_range_overrides() -> None:
    chains = [
        {"name": "Fragment", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 3},
    ]
    manifest = [
        {
            "chain_id": "A",
            "uniprot_ac": "ACC_A",
            "entity_id": "1",
            "is_fragment": "true",
            "sequence_start": "100",
            "sequence_end": "102",
        }
    ]

    rows, _ = _compute_plddt_metrics([90.0, 80.0, 70.0], chains, manifest)

    assert rows[0]["sequence_start"] == 100
    assert rows[0]["sequence_end"] == 102
    assert rows[0]["is_fragment"] == "true"


def test_compute_plddt_metrics_propagates_provenance_fields() -> None:
    chains = [
        {"name": "Named", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 2},
    ]
    manifest = [
        {
            "chain_id": "A",
            "uniprot_ac": "ACC_A",
            "entity_id": "7",
            "is_fragment": "true",
            "is_isoform": "true",
            "entity_type": "polypeptide",
            "protein_name": "Custom name",
        }
    ]

    rows, _ = _compute_plddt_metrics([90.0, 80.0], chains, manifest)

    row = rows[0]
    assert row["is_fragment"] == "true"
    assert row["is_isoform"] == "true"
    assert row["entity_type"] == "polypeptide"
    assert row["protein_name"] == "Custom name"


def test_compute_plddt_metrics_defaults_entity_type_and_omits_protein_name() -> None:
    chains = [
        {"name": "Plain", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 2},
    ]
    manifest = [{"chain_id": "A", "uniprot_ac": "ACC_A", "entity_id": "1"}]

    rows, _ = _compute_plddt_metrics([90.0, 80.0], chains, manifest)

    assert rows[0]["entity_type"] == "protein"
    assert "protein_name" not in rows[0]


def test_compute_plddt_metrics_rounds_fractions_deterministically() -> None:
    """``round_float`` must replace builtin ``round`` for confidence fractions.

    1351/2000 == 0.6755 and 649/2000 == 0.3245 sit on binary-float boundaries
    where builtin ``round`` disagrees with Decimal half-even.
    """
    n = 2000
    values = [40.0] * 1351 + [80.0] * 649
    chains = [
        {"name": "Chain A", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": n},
    ]
    manifest = [{"chain_id": "A", "uniprot_ac": "P00000", "entity_id": "1"}]

    rows, model_avg = _compute_plddt_metrics(values, chains, manifest)

    row = rows[0]
    assert row["average_plddt"] == 52.98
    assert row["fraction_plddt_very_low"] == 0.676  # builtin round would give 0.675
    assert row["fraction_plddt_low"] == 0.0
    assert row["fraction_plddt_confident"] == 0.324  # builtin round would give 0.325
    assert row["fraction_plddt_very_high"] == 0.0
    assert model_avg == 52.98


def test_compute_plddt_metrics_rounds_means_deterministically() -> None:
    """Chain and model means must use ``round_float`` on half-cent boundaries."""
    chains = [
        {"name": "Chain A", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 1},
    ]
    manifest = [{"chain_id": "A", "uniprot_ac": "P00000", "entity_id": "1"}]

    rows, model_avg = _compute_plddt_metrics([2.675], chains, manifest)

    # builtin round(2.675, 2) == 2.67; Decimal half-even gives 2.68.
    assert rows[0]["average_plddt"] == 2.68
    assert model_avg == 2.68


def test_convert_file_uses_manifest_sequence_overrides_with_duckdb(tmp_path: Path) -> None:
    """Fragment provenance ranges must be emitted from the manifest, while the
    pLDDT/PAE chain boundaries keep following the physical structure."""
    duckdb_path = tmp_path / "entries.duckdb"
    manifest_path = tmp_path / "manifest.csv"
    scores_json = tmp_path / "scores.json"
    pdb_file = tmp_path / "test.pdb"
    chain_manifest = tmp_path / "chain_manifest.csv"

    _create_duckdb(
        duckdb_path,
        [
            ("P11111", "Protein one", "ABCDEFGHIJKLMNOP"),
            ("Q22222", "Protein two", "QRSTUVWXYZABCDEFG"),
        ],
    )
    manifest_path.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac,sequence_start,sequence_end,is_fragment\n"
        "TEST_MODEL,1,A,P11111,10,12,true\n"
        "TEST_MODEL,2,B,Q22222,5,9,true\n",
        encoding="utf-8",
    )
    _write_test_pdb(pdb_file, [("A", "AAA"), ("B", "GGGGG")])
    _write_scores_json(
        scores_json,
        plddt=[90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0],
        pae=[[0.0] * 8 for _ in range(8)],
        max_pae=3.0,
    )

    output_paths = convert_file(
        str(scores_json),
        str(pdb_file),
        outdir=str(tmp_path),
        manifest_path=str(manifest_path),
        model_entity_id="TEST_MODEL",
        duckdb_path=str(duckdb_path),
        out_chain_manifest=str(chain_manifest),
    )

    plddt_payload = json.loads(
        Path(output_paths["plddt"]).read_text(encoding="utf-8")
    )
    assert [
        (c["label_asym_id"], c["sequenceStart"], c["sequenceEnd"])
        for c in plddt_payload["chains"]
    ] == [("A", 1, 3), ("B", 1, 5)]

    rows = list(
        csv.DictReader(chain_manifest.read_text(encoding="utf-8").splitlines())
    )
    assert [
        (r["chain_id"], r["sequence_start"], r["sequence_end"], r["is_fragment"])
        for r in rows
    ] == [("A", "10", "12", "true"), ("B", "5", "9", "true")]


def test_convert_file_orders_chains_by_structure(tmp_path: Path):
    """End-to-end: manifest lists the accessions in the WRONG (reversed) order;
    the converter must still emit chain names/boundaries that match the structure."""
    # Structure: chain A = 3 ALA residues, chain B = 5 GLY residues.
    pdb_file = tmp_path / "model.pdb"
    _write_pdb(pdb_file, [("A", "ALA", 3), ("B", "GLY", 5)])
    n = 8

    scores = {
        "plddt": [90.0] * n,
        "pae": [[0.0] * n for _ in range(n)],
        "max_pae": 5.0,
    }
    scores_json = tmp_path / "scores.json"
    scores_json.write_text(json.dumps(scores))

    # Manifest order is reversed vs the structure (chain A -> the 5-res accession).
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "model_entity_id,chain_id,uniprot_ac\n"
        "TEST,A,ACC_GLY5\n"
        "TEST,B,ACC_ALA3\n"
    )

    duckdb_path = tmp_path / "fake.duckdb"
    duckdb_path.write_text("")  # not read; prefetch cache is injected below
    cache_key = str(duckdb_path.resolve())
    converter._DUCKDB_METADATA_CACHE[cache_key] = {
        "ACC_ALA3": {"primary_ac": "ACC_ALA3", "protein_full_names": ["Alanine protein"], "sequence": "AAA"},
        "ACC_GLY5": {"primary_ac": "ACC_GLY5", "protein_full_names": ["Glycine protein"], "sequence": "GGGGG"},
    }
    try:
        out = convert_file(
            str(scores_json),
            str(pdb_file),
            outdir=str(tmp_path),
            manifest_path=str(manifest),
            model_entity_id="TEST",
            duckdb_path=str(duckdb_path),
        )
        pae_payload = json.loads(Path(out["pae"]).read_text())[0]
        plddt_payload = json.loads(Path(out["plddt"]).read_text())
    finally:
        converter._DUCKDB_METADATA_CACHE.pop(cache_key, None)

    expected = [
        {"name": "Alanine protein", "label_asym_id": "A", "sequenceStart": 1, "sequenceEnd": 3},
        {"name": "Glycine protein", "label_asym_id": "B", "sequenceStart": 1, "sequenceEnd": 5},
    ]
    assert pae_payload["chains"] == expected
    assert plddt_payload["chains"] == expected


# ---------------------------------------------------------------------------
# Duplicate-chain-row dedup (stage_03 residue-count 2x double-count fix).
# A renamed heterodimer can receive duplicate chain rows in the work manifest
# (forward + swapped compound IDs renaming to one AF-ID), so the converter sees
# 4 chains for a 2-chain structure -> fallback path -> residue count doubles.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "dup_rows",
    [
        # verbatim duplicate rows (the identical-row case)
        "TEST,1,A,ACC_ALA3\nTEST,2,B,ACC_GLY5\n",
        # swapped duplicate rows (forward + entity/chain-flipped swapped compound)
        "TEST,1,A,ACC_GLY5\nTEST,2,B,ACC_ALA3\n",
    ],
)
def test_convert_file_dedups_duplicate_manifest_chain_rows(tmp_path: Path, dup_rows: str):
    """The work manifest lists each chain twice (4 rows for a 2-chain model).
    The converter must dedup to 2 chains so the residue count matches the
    structure (8, not 16) and chains bind to the structure order."""
    pdb_file = tmp_path / "model.pdb"
    _write_pdb(pdb_file, [("A", "ALA", 3), ("B", "GLY", 5)])
    n = 8
    scores = {"plddt": [90.0] * n, "pae": [[0.0] * n for _ in range(n)], "max_pae": 5.0}
    scores_json = tmp_path / "scores.json"
    scores_json.write_text(json.dumps(scores))

    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "model_entity_id,entity_id,chain_id,uniprot_ac\n"
        "TEST,1,A,ACC_ALA3\n"
        "TEST,2,B,ACC_GLY5\n" + dup_rows
    )

    duckdb_path = tmp_path / "fake.duckdb"
    duckdb_path.write_text("")
    cache_key = str(duckdb_path.resolve())
    converter._DUCKDB_METADATA_CACHE[cache_key] = {
        "ACC_ALA3": {"primary_ac": "ACC_ALA3", "protein_full_names": ["Alanine protein"], "sequence": "AAA"},
        "ACC_GLY5": {"primary_ac": "ACC_GLY5", "protein_full_names": ["Glycine protein"], "sequence": "GGGGG"},
    }
    try:
        out = convert_file(
            str(scores_json),
            str(pdb_file),
            outdir=str(tmp_path),
            manifest_path=str(manifest),
            model_entity_id="TEST",
            duckdb_path=str(duckdb_path),
        )
        plddt_payload = json.loads(Path(out["plddt"]).read_text())
    finally:
        converter._DUCKDB_METADATA_CACHE.pop(cache_key, None)

    # residueNumber must follow the structure (8), not the doubled count (16).
    assert len(plddt_payload["residueNumber"]) == 8
    # deduped to 2 chains, bound to the structure order.
    assert [(c["label_asym_id"], c["name"], c["sequenceEnd"]) for c in plddt_payload["chains"]] == [
        ("A", "Alanine protein", 3),
        ("B", "Glycine protein", 5),
    ]


def test_load_chain_metadata_dedups_duplicate_rows(tmp_path: Path):
    """Unit: 4 duplicate manifest rows + a 2-chain structure -> 2 chains, residue
    count == structure length, and the structure-binding path is taken."""
    pdb_file = tmp_path / "model.pdb"
    _write_pdb(pdb_file, [("A", "ALA", 3), ("B", "GLY", 5)])
    manifest_chains = [
        {"chain_id": "A", "entity_id": "1", "uniprot_ac": "ACC_ALA3"},
        {"chain_id": "B", "entity_id": "2", "uniprot_ac": "ACC_GLY5"},
        {"chain_id": "A", "entity_id": "1", "uniprot_ac": "ACC_ALA3"},
        {"chain_id": "B", "entity_id": "2", "uniprot_ac": "ACC_GLY5"},
    ]
    duckdb_path = tmp_path / "fake.duckdb"
    duckdb_path.write_text("")
    cache_key = str(duckdb_path.resolve())
    converter._DUCKDB_METADATA_CACHE[cache_key] = {
        "ACC_ALA3": {"primary_ac": "ACC_ALA3", "protein_full_names": ["Alanine protein"], "sequence": "AAA"},
        "ACC_GLY5": {"primary_ac": "ACC_GLY5", "protein_full_names": ["Glycine protein"], "sequence": "GGGGG"},
    }
    try:
        chains, residue_numbers, effective = converter._load_chain_metadata_from_duckdb(
            duckdb_path, manifest_chains=manifest_chains, pdb_path=pdb_file,
        )
    finally:
        converter._DUCKDB_METADATA_CACHE.pop(cache_key, None)

    assert len(chains) == 2
    assert len(residue_numbers) == 8
    assert len(effective) == 2
    assert [c["label_asym_id"] for c in chains] == ["A", "B"]


def test_load_chain_metadata_homomultimer_still_expands(tmp_path: Path):
    """Guard: dedup must NOT break genuine homomultimer auto-expansion — a 1-row
    manifest + a 2-chain structure (same accession) must still expand to 2."""
    pdb_file = tmp_path / "homo.pdb"
    _write_pdb(pdb_file, [("A", "MET", 50), ("B", "MET", 50)])
    manifest_chains = [{"chain_id": "A", "entity_id": "1", "uniprot_ac": "ACC_HOMO"}]
    duckdb_path = tmp_path / "fake.duckdb"
    duckdb_path.write_text("")
    cache_key = str(duckdb_path.resolve())
    converter._DUCKDB_METADATA_CACHE[cache_key] = {
        "ACC_HOMO": {"primary_ac": "ACC_HOMO", "protein_full_names": ["Homo protein"], "sequence": "M" * 50},
    }
    try:
        _chains, residue_numbers, effective = converter._load_chain_metadata_from_duckdb(
            duckdb_path, manifest_chains=manifest_chains, pdb_path=pdb_file,
        )
    finally:
        converter._DUCKDB_METADATA_CACHE.pop(cache_key, None)

    assert len(effective) == 2
    assert len(residue_numbers) == 100
