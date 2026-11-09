"""Regression tests for ALTERNATIVE PRODUCTS isoform parsing.

A3R4N1 was merged into P0DOJ1/2/3 in UniProt, so ``A3R4N1-1`` is a legacy
isoform ID listed as a synonym (``IsoId=P0DOJ3-1, A3R4N1-1``) in all three
entries. It must not become a row, otherwise ``primary_ac`` is duplicated.
"""

from pathlib import Path

from uniprot.scripts.extract_subset import (
    build_entry_payload,
    extract_accessions,
    parse_alt_product_groups,
    parse_alt_products,
)

FIXTURES = Path(__file__).parent / "fixtures" / "uniprot"


def _load(name):
    return (FIXTURES / name).read_text().splitlines()


def test_synonym_isoform_ids_are_not_returned():
    lines = _load("P0DOJ3_VP2_POVK3.dat")
    isoforms = parse_alt_products(lines, "P0DOJ3")
    assert isoforms == {"P0DOJ3-1": ["Displayed"], "P0DOJ3-2": ["VSP_059273"]}


def test_synonyms_are_kept_on_the_group():
    groups = parse_alt_product_groups(_load("P0DOJ3_VP2_POVK3.dat"), "P0DOJ3")
    assert [(g.isoid, g.synonyms) for g in groups] == [
        ("P0DOJ3-1", ("A3R4N1-1",)),
        ("P0DOJ3-2", ("A3R4N1-2",)),
    ]


def test_external_isoforms_are_skipped():
    lines = _load("P0DOI4_VP1_POVK3.dat")
    assert list(parse_alt_products(lines, "P0DOI4")) == ["P0DOI4-1"]


def test_entry_rows_are_unique_and_have_correct_sequences():
    lines = _load("P0DOJ3_VP2_POVK3.dat")
    rows = build_entry_payload(lines, extract_accessions(lines), "2025_04", True)
    acs = [r["primary_ac"] for r in rows]
    assert acs == ["P0DOJ3", "P0DOJ3-1", "P0DOJ3-2"]
    assert len(set(acs)) == len(acs)
    by_ac = {r["primary_ac"]: r for r in rows}
    assert by_ac["P0DOJ3-1"]["length"] == 400
    assert by_ac["P0DOJ3-1"]["md5"] == "a1fba4cc534e575606d4cc60164a32ea"
    assert by_ac["P0DOJ3-2"]["length"] == 257
    assert by_ac["P0DOJ3-2"]["md5"] == "06de9d1563101397e05e8674de365957"
    assert not any(a.startswith("A3R4N1") for a in acs)


def test_single_position_var_seq_isoform_is_built():
    # CDC42_MOUSE isoform 1 = VSP_040585 (single residue 163) + VSP_040586 (182..191)
    lines = _load("P60766_CDC42_MOUSE.dat")
    assert parse_alt_products(lines, "P60766") == {
        "P60766-2": ["Displayed"],
        "P60766-1": ["VSP_040585", "VSP_040586"],
    }
    rows = build_entry_payload(lines, extract_accessions(lines), "2025_04", True)
    by_ac = {r["primary_ac"]: r for r in rows}
    assert set(by_ac) == {"P60766", "P60766-1", "P60766-2"}
    canonical = by_ac["P60766"]["sequence"]
    iso1 = by_ac["P60766-1"]["sequence"]
    assert len(iso1) == 191 and iso1[162] == "R" and iso1[181:] == "TQPKRKCCIF"
    assert iso1[:162] == canonical[:162]
    assert not any(a.startswith("P21181") for a in by_ac)


def test_wrapped_var_seq_note_gives_full_replacement():
    # SNAP25a/b differ by a 32-residue swap whose /note wraps over two FT lines;
    # both isoforms are 206 aa.
    lines = _load("P60879_SNP25_MOUSE.dat")
    rows = build_entry_payload(lines, extract_accessions(lines), "2025_04", True)
    by_ac = {r["primary_ac"]: r for r in rows}
    assert by_ac["P60879-2"]["length"] == 206
    assert by_ac["P60879-2"]["sequence"][57:89] == "DRVEEGMNHINQDMKEAEKNLKDLGKCCGLFI"


def test_wrapped_var_seq_note_with_long_replacement():
    # THA_RAT isoform Alpha-1 replaces 371..492 (122 aa) with a 41-aa tail.
    lines = _load("P63059_THA_RAT.dat")
    rows = build_entry_payload(lines, extract_accessions(lines), "2025_04", True)
    by_ac = {r["primary_ac"]: r for r in rows}
    canonical_len = by_ac["P63059"]["length"]
    assert by_ac["P63059-1"]["length"] == canonical_len
    assert by_ac["P63059-2"]["length"] == 410  # 492 - 122 + 40
