"""Pure parsers for first-class batch interval input."""

from collections.abc import Iterable

from .batch import BatchInputRecord
from .models import AssemblyIdentifier, GenomicInterval


class BatchInputError(ValueError):
    """Raised when one batch input row has invalid interval semantics."""


def parse_bed_batch(
    lines: Iterable[str],
    *,
    assembly: AssemblyIdentifier,
) -> tuple[BatchInputRecord, ...]:
    """Parse BED3-or-later rows into canonical non-empty source intervals.

    Blank lines, ``#`` comments, and UCSC ``track``/``browser`` header lines are
    ignored.  Data rows must be tab-delimited and provide at least BED columns
    ``chrom``, ``chromStart``, and ``chromEnd``.  The optional fourth BED column is
    preserved as a display label.  BED coordinates remain 0-based, half-open.
    """

    records: list[BatchInputRecord] = []
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "track ", "browser ")):
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            raise BatchInputError(
                f"BED line {line_number} must contain at least 3 tab-delimited columns"
            )

        sequence_name = fields[0]
        _validate_sequence_name(
            sequence_name,
            line_number=line_number,
            input_name="BED",
        )

        start = _parse_non_negative_coordinate(
            fields[1], line_number=line_number, field="start", input_name="BED"
        )
        end = _parse_non_negative_coordinate(
            fields[2], line_number=line_number, field="end", input_name="BED"
        )
        if end <= start:
            raise BatchInputError(
                f"BED line {line_number} must span at least one base; "
                "chromEnd must be greater than chromStart"
            )

        label = fields[3] if len(fields) >= 4 and fields[3] else None
        records.append(
            BatchInputRecord(
                record_id=f"row-{len(records) + 1}",
                source_interval=GenomicInterval(
                    assembly=assembly,
                    sequence_name=sequence_name,
                    start=start,
                    end=end,
                ),
                source_line_number=line_number,
                label=label,
            )
        )

    if not records:
        raise BatchInputError("BED input must contain at least one non-empty data row")
    return tuple(records)


def parse_interval_table_batch(
    lines: Iterable[str],
    *,
    assembly: AssemblyIdentifier,
) -> tuple[BatchInputRecord, ...]:
    """Parse a simple tab-delimited 1-based inclusive interval table.

    The first non-empty, non-comment line must be exactly
    ``sequence<TAB>start<TAB>end`` or
    ``sequence<TAB>start<TAB>end<TAB>label``. Data rows use the same 1-based,
    inclusive coordinate convention as the single-locus CLI. They are normalized to
    the canonical 0-based, half-open intervals used by the batch engine.
    """

    records: list[BatchInputRecord] = []
    expected_field_count: int | None = None
    header_line_number: int | None = None

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields = line.split("\t")
        if expected_field_count is None:
            header = tuple(fields)
            if header not in {
                ("sequence", "start", "end"),
                ("sequence", "start", "end", "label"),
            }:
                raise BatchInputError(
                    "interval-table header must be tab-delimited and exactly "
                    "'sequence\\tstart\\tend' or "
                    "'sequence\\tstart\\tend\\tlabel'"
                )
            expected_field_count = len(header)
            header_line_number = line_number
            continue

        if len(fields) != expected_field_count:
            raise BatchInputError(
                f"interval-table line {line_number} must contain exactly "
                f"{expected_field_count} tab-delimited columns to match the header"
            )

        sequence_name = fields[0]
        _validate_sequence_name(
            sequence_name,
            line_number=line_number,
            input_name="interval-table",
        )
        start_1based = _parse_positive_coordinate(
            fields[1],
            line_number=line_number,
            field="start",
            input_name="interval-table",
        )
        end_1based = _parse_positive_coordinate(
            fields[2],
            line_number=line_number,
            field="end",
            input_name="interval-table",
        )
        if end_1based < start_1based:
            raise BatchInputError(
                f"interval-table line {line_number} end must be greater than or "
                "equal to start for 1-based inclusive coordinates"
            )

        label = fields[3] if expected_field_count == 4 and fields[3] else None
        records.append(
            BatchInputRecord(
                record_id=f"row-{len(records) + 1}",
                source_interval=GenomicInterval(
                    assembly=assembly,
                    sequence_name=sequence_name,
                    start=start_1based - 1,
                    end=end_1based,
                ),
                source_line_number=line_number,
                label=label,
            )
        )

    if expected_field_count is None:
        raise BatchInputError(
            "interval-table input must contain a header and at least one data row"
        )
    if not records:
        assert header_line_number is not None
        raise BatchInputError(
            "interval-table input must contain at least one data row after the header "
            f"on line {header_line_number}"
        )
    return tuple(records)


def _validate_sequence_name(
    sequence_name: str,
    *,
    line_number: int,
    input_name: str,
) -> None:
    if not sequence_name or any(character.isspace() for character in sequence_name):
        raise BatchInputError(
            f"{input_name} line {line_number} sequence name must not be empty "
            "or contain whitespace"
        )


def _parse_non_negative_coordinate(
    text: str,
    *,
    line_number: int,
    field: str,
    input_name: str,
) -> int:
    if not text or not text.isdigit():
        raise BatchInputError(
            f"{input_name} line {line_number} {field} must be a non-negative integer"
        )
    return int(text)


def _parse_positive_coordinate(
    text: str,
    *,
    line_number: int,
    field: str,
    input_name: str,
) -> int:
    value = _parse_non_negative_coordinate(
        text,
        line_number=line_number,
        field=field,
        input_name=input_name,
    )
    if value < 1:
        raise BatchInputError(
            f"{input_name} line {line_number} {field} must be at least 1 for "
            "1-based inclusive coordinates"
        )
    return value
