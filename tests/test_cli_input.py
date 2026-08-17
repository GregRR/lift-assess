import pytest

from liftassess.cli_input import parse_ucsc_locus, ucsc_assembly_identifier
from liftassess.models import AssemblyIdentifier, GenomicInterval


def test_ucsc_assembly_identifier_records_exact_db_without_inferred_aliases() -> None:
    assembly = ucsc_assembly_identifier("canFam3")

    assert assembly == AssemblyIdentifier(name="canFam3", provider="UCSC")
    assert assembly.accession is None
    assert assembly.aliases == ()


@pytest.mark.parametrize("db", ["", "can Fam3", "canFam3/other"])
def test_ucsc_assembly_identifier_rejects_invalid_database_names(db: str) -> None:
    with pytest.raises(ValueError):
        ucsc_assembly_identifier(db)


def test_parse_ucsc_locus_converts_display_coordinates_to_half_open() -> None:
    assembly = ucsc_assembly_identifier("canFam3")

    interval = parse_ucsc_locus("chr16:12345-12400", assembly=assembly)

    assert interval == GenomicInterval(
        assembly=assembly,
        sequence_name="chr16",
        start=12_344,
        end=12_400,
    )


def test_parse_ucsc_locus_accepts_browser_comma_grouping() -> None:
    assembly = ucsc_assembly_identifier("canFam3")

    interval = parse_ucsc_locus("chrUn_AAEX03020568:12,345-12,400", assembly=assembly)

    assert interval.start == 12_344
    assert interval.end == 12_400


def test_parse_ucsc_locus_single_base_preserves_one_base_length() -> None:
    assembly = ucsc_assembly_identifier("canFam3")

    interval = parse_ucsc_locus("chr1:1-1", assembly=assembly)

    assert interval.start == 0
    assert interval.end == 1
    assert interval.length == 1


@pytest.mark.parametrize(
    ("locus", "message"),
    [
        ("", "sequence:start-end"),
        ("chr1", "sequence:start-end"),
        ("chr1:1", "sequence:start-end"),
        ("chr1:0-1", "at least 1"),
        ("chr1:10-9", "greater than or equal"),
        ("chr1:1,23-200", "invalid comma grouping"),
        ("chr1:1 0-20", "positive integer"),
        ("chr 1:10-20", "sequence name"),
    ],
)
def test_parse_ucsc_locus_rejects_malformed_display_coordinates(
    locus: str,
    message: str,
) -> None:
    assembly = ucsc_assembly_identifier("canFam3")

    with pytest.raises(ValueError, match=message):
        parse_ucsc_locus(locus, assembly=assembly)
