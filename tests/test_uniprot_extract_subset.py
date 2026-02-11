from uniprot.scripts.extract_subset import parse_de_sections


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
