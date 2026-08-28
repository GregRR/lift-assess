import pytest

from liftassess.batch_input import BatchInputError, parse_bed_batch
from liftassess.models import AssemblyIdentifier, GenomicInterval

ASSEMBLY = AssemblyIdentifier(name="hg38", provider="UCSC")


def test_parse_bed_batch_preserves_native_half_open_coordinates_and_label() -> None:
    records = parse_bed_batch(
        ["chr1\t10\t20\tfirst\n", "chr2\t30\t31\n"],
        assembly=ASSEMBLY,
    )

    assert records[0].record_id == "row-1"
    assert records[0].source_line_number == 1
    assert records[0].label == "first"
    assert records[0].source_interval == GenomicInterval(
        assembly=ASSEMBLY,
        sequence_name="chr1",
        start=10,
        end=20,
    )
    assert records[1].record_id == "row-2"
    assert records[1].label is None
    assert records[1].source_interval.length == 1


def test_parse_bed_batch_ignores_blank_comment_track_and_browser_lines() -> None:
    records = parse_bed_batch(
        [
            "# comment\n",
            "track name=test\n",
            "\n",
            "browser position chr1:1-10\n",
            "chr1\t0\t1\tpoint\n",
        ],
        assembly=ASSEMBLY,
    )

    assert len(records) == 1
    assert records[0].record_id == "row-1"
    assert records[0].source_line_number == 5


@pytest.mark.parametrize(
    ("line", "message"),
    [
        ("chr1 0 1\n", "tab-delimited"),
        ("chr 1\t0\t1\n", "sequence name"),
        ("chr1\tone\t2\n", "start"),
        ("chr1\t0\ttwo\n", "end"),
        ("chr1\t10\t10\n", "span at least one base"),
        ("chr1\t10\t9\n", "span at least one base"),
    ],
)
def test_parse_bed_batch_rejects_invalid_rows(line: str, message: str) -> None:
    with pytest.raises(BatchInputError, match=message):
        parse_bed_batch([line], assembly=ASSEMBLY)


def test_parse_bed_batch_rejects_input_without_data_rows() -> None:
    with pytest.raises(BatchInputError, match="at least one"):
        parse_bed_batch(["# comment\n", "\n"], assembly=ASSEMBLY)
