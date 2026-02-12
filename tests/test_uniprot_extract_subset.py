from uniprot.scripts.extract_subset import (
    build_entry_payload,
    parse_alt_products,
    parse_de_sections,
    parse_var_seq,
)


def test_parse_de_sections_includes_subname_full() -> None:
    lines = [
        "ID   Q64890_ADEM1            Unreviewed;        97 AA.",
        "AC   Q64890;",
        "DE   SubName: Full=Adenovirus type 1 early regions 1A and 1B DNA {ECO:0000313|EMBL:AAA42426.1};",
    ]

    full_names, short_names = parse_de_sections(lines)

    assert full_names == ["Adenovirus type 1 early regions 1A and 1B DNA"]
    assert short_names == []


def test_parse_de_sections_prefers_recname_before_altname() -> None:
    lines = [
        "DE   RecName: Full=Genome polyprotein;",
        "DE   AltName: Full=Core protein;",
        "DE   RecName: Short=NS1;",
    ]

    full_names, short_names = parse_de_sections(lines)

    assert full_names[0] == "Genome polyprotein"
    assert "Core protein" in full_names
    assert short_names == ["NS1"]


def test_parse_de_sections_includes_alt_special_name_fields() -> None:
    lines = [
        "DE   AltName: Allergen=Api m 1;",
        "DE   AltName: Biotech=Etanercept;",
        "DE   AltName: CD_antigen=CD177;",
        "DE   AltName: INN=Adalimumab;",
    ]

    full_names, short_names = parse_de_sections(lines)

    assert full_names == ["Api m 1", "Etanercept", "CD177", "Adalimumab"]
    assert short_names == []


def test_parse_alt_products_splits_multiple_isoids() -> None:
    lines = [
        "CC   -!- ALTERNATIVE PRODUCTS:",
        "CC       Event=Alternative splicing; Named isoforms=3;",
        "CC         Name=1;",
        "CC           IsoId=P12345-1, P12345-2; Sequence=Displayed;",
        "CC         Name=2;",
        "CC           IsoId=P12345-3; Sequence=VSP_000001, VSP_000002;",
        "CC   -!- SUBCELLULAR LOCATION: Cytoplasm.",
    ]

    isoforms = parse_alt_products(lines)

    assert isoforms["P12345-1"] == ["Displayed"]
    assert isoforms["P12345-2"] == ["Displayed"]
    assert isoforms["P12345-3"] == ["VSP_000001", "VSP_000002"]


def test_build_entry_payload_skips_external_isoforms() -> None:
    lines = [
        "ID   TEST_HUMAN              Reviewed;         5 AA.",
        "AC   P12345;",
        "DE   RecName: Full=Test protein;",
        "SQ   SEQUENCE   5 AA;  500 MW;  ABCDEF1234567890 CRC64;",
        "     MTEST",
        "CC   -!- ALTERNATIVE PRODUCTS:",
        "CC       Event=Alternative splicing; Named isoforms=2;",
        "CC         Name=1;",
        "CC           IsoId=P12345-1; Sequence=Displayed;",
        "CC         Name=2;",
        "CC           IsoId=P12345-2; Sequence=External;",
    ]

    rows = build_entry_payload(lines, ["P12345"], "2025_04", reviewed=True)
    accessions = [row["primary_ac"] for row in rows]

    assert "P12345-1" in accessions
    assert "P12345-2" not in accessions


def test_parse_var_seq_multiline_note_preserves_full_replacement() -> None:
    lines = [
        "FT   VAR_SEQ         1056..1097",
        "FT                   /note=\"SPLRHDGTPVPARRRPLGHGFGLAHPGMMQELQARLGRPKPQ -> RWEDRL",
        "FT                   RPGVRDQPGQHSKIPIF (in isoform 3)\"",
        "FT                   /id=\"VSP_060194\"",
    ]

    varseqs = parse_var_seq(lines)

    assert varseqs["VSP_060194"] == (1056, 1097, "RWEDRLRPGVRDQPGQHSKIPIF")
