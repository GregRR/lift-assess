"""Pure input-boundary helpers for the planned liftAssess CLI.

The CLI accepts UCSC database identifiers and UCSC-style display loci, while the
scientific core uses structured assembly identities and 0-based, half-open intervals.
Keeping these conversions here makes the boundary explicit and independently testable
before resource discovery, downloading, or report rendering is involved.

For v1, a UCSC database name is canonicalized conservatively to exactly that database
identifier plus the ``UCSC`` provider. The CLI does not infer NCBI accessions or
biological assembly aliases, and it does not redefine ``AssemblyIdentifier`` equality.
General cross-provider assembly alias resolution remains outside v1 scope.
"""

from .models import AssemblyIdentifier, GenomicInterval
from .resources import _validate_ucsc_db


def ucsc_assembly_identifier(db: str) -> AssemblyIdentifier:
    """Return the minimal structured identity for one explicit UCSC database name."""

    _validate_ucsc_db(db)
    return AssemblyIdentifier(name=db, provider="UCSC")


def parse_ucsc_locus(
    locus: str,
    *,
    assembly: AssemblyIdentifier,
) -> GenomicInterval:
    """Parse a UCSC-style 1-based inclusive locus into canonical coordinates.

    ``chr16:12345-12400`` becomes the internal half-open interval
    ``[12344, 12400)``. Comma-grouped display coordinates copied from a browser,
    such as ``chr16:12,345-12,400``, are accepted as well. Only surrounding
    whitespace is ignored; whitespace inside the locus is rejected rather than
    guessed through.
    """

    text = locus.strip()
    if not text or text.count(":") != 1:
        raise ValueError("locus must use UCSC-style 'sequence:start-end' coordinates")

    sequence_name, coordinate_text = text.split(":", 1)
    if not sequence_name or any(character.isspace() for character in sequence_name):
        raise ValueError("locus sequence name must not be empty or contain whitespace")
    if coordinate_text.count("-") != 1:
        raise ValueError("locus must use UCSC-style 'sequence:start-end' coordinates")

    start_text, end_text = coordinate_text.split("-", 1)
    start = _parse_display_coordinate(start_text, field="start")
    end = _parse_display_coordinate(end_text, field="end")
    if end < start:
        raise ValueError("locus end must be greater than or equal to start")

    return GenomicInterval(
        assembly=assembly,
        sequence_name=sequence_name,
        start=start - 1,
        end=end,
    )


def _parse_display_coordinate(text: str, *, field: str) -> int:
    if not text or any(character.isspace() for character in text):
        raise ValueError(f"locus {field} must be a positive integer")

    if "," in text:
        groups = text.split(",")
        if (
            not 1 <= len(groups[0]) <= 3
            or not groups[0].isdigit()
            or any(len(group) != 3 or not group.isdigit() for group in groups[1:])
        ):
            raise ValueError(f"locus {field} has invalid comma grouping")
    elif not text.isdigit():
        raise ValueError(f"locus {field} must be a positive integer")

    value = int(text.replace(",", ""))
    if value < 1:
        raise ValueError(f"locus {field} must be at least 1 in display coordinates")
    return value
