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
        if not sequence_name or any(character.isspace() for character in sequence_name):
            raise BatchInputError(
                f"BED line {line_number} sequence name must not be empty "
                "or contain whitespace"
            )

        start = _parse_bed_coordinate(fields[1], line_number=line_number, field="start")
        end = _parse_bed_coordinate(fields[2], line_number=line_number, field="end")
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


def _parse_bed_coordinate(text: str, *, line_number: int, field: str) -> int:
    if not text or not text.isdigit():
        raise BatchInputError(
            f"BED line {line_number} {field} must be a non-negative integer"
        )
    value = int(text)
    if value < 0:
        raise BatchInputError(
            f"BED line {line_number} {field} must be a non-negative integer"
        )
    return value
