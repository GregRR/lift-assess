"""Minimal streaming reader for the UCSC net text format.

The parser preserves the target-section identity, record indentation depth, and
v1-relevant optional fields without assigning evidence or matching net records
to normalized candidates. UCSC net indentation is semantic: one leading space
represents one hierarchy level.
"""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from enum import Enum

from .models import MappingOrientation


class NetFormatError(ValueError):
    """Raised when UCSC net text is structurally invalid."""


class NetRecordKind(str, Enum):
    """The two record classes defined by the UCSC net format."""

    FILL = "fill"
    GAP = "gap"


class NetClassification(str, Enum):
    """UCSC net classifications produced by netSyntenic/netClass."""

    TOP = "top"
    SYNTENIC = "syn"
    INVERSION = "inv"
    NON_SYNTENIC = "nonSyn"


@dataclass(frozen=True)
class NetRecord:
    """One ``fill`` or ``gap`` record from a UCSC net section.

    ``depth`` is the semantic indentation level: top-level records have depth 1,
    their children depth 2, and so on. ``attributes`` preserves every optional
    name/value pair exactly as text, including fields liftAssess does not yet
    interpret. The typed convenience fields expose only v1-relevant attributes.
    """

    target_name: str
    target_sequence_size: int
    depth: int
    kind: NetRecordKind
    target_start: int
    target_span_size: int
    query_name: str
    orientation: MappingOrientation
    query_start: int
    query_span_size: int
    chain_id: int | None = None
    score: int | None = None
    aligned_bases: int | None = None
    duplicated_query_bases: int | None = None
    classification: NetClassification | None = None
    attributes: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if not self.target_name:
            raise ValueError("net target name must not be empty")
        if self.target_sequence_size < 0:
            raise ValueError("net target sequence size must be non-negative")
        if self.depth < 1:
            raise ValueError("net record depth must be at least 1")
        if self.target_start < 0 or self.target_span_size < 0:
            raise ValueError("net target coordinates must be non-negative")
        if self.target_start + self.target_span_size > self.target_sequence_size:
            raise ValueError("net target span exceeds target sequence bounds")
        if not self.query_name:
            raise ValueError("net query name must not be empty")
        if self.query_start < 0 or self.query_span_size < 0:
            raise ValueError("net query coordinates must be non-negative")
        if self.chain_id is not None and self.chain_id < 0:
            raise ValueError("net chain ID must be non-negative")
        if self.score is not None and self.score < 0:
            raise ValueError("net score must be non-negative")
        if self.aligned_bases is not None and self.aligned_bases < 0:
            raise ValueError("net aligned bases must be non-negative")
        if self.duplicated_query_bases is not None and self.duplicated_query_bases < 0:
            raise ValueError("net duplicated query bases must be non-negative")

    @property
    def target_end(self) -> int:
        """Return the 0-based, half-open target end coordinate."""

        return self.target_start + self.target_span_size

    @property
    def query_end(self) -> int:
        """Return the stored query end coordinate from start + span size."""

        return self.query_start + self.query_span_size


@dataclass(frozen=True)
class _NetSection:
    target_name: str
    target_sequence_size: int


def iter_net_records(lines: Iterable[str]) -> Iterator[NetRecord]:
    """Yield ``fill``/``gap`` records from UCSC net text in file order.

    The function streams individual records and does not materialize a whole
    chromosome net tree. Hierarchy is preserved as the integer ``depth`` on
    each record. Blank lines and ``#`` comment/metadata lines are ignored.
    """

    section: _NetSection | None = None
    previous_depth: int | None = None

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if "\t" in line[: len(line) - len(line.lstrip(" \t"))]:
            raise NetFormatError(
                f"line {line_number}: net hierarchy indentation must use spaces"
            )

        indent = len(line) - len(line.lstrip(" "))
        content = line[indent:]
        fields = content.split()

        if fields[0] == "net":
            if indent != 0:
                raise NetFormatError(
                    f"line {line_number}: net section header must not be indented"
                )
            section = _parse_section_header(fields, line_number)
            previous_depth = None
            continue

        if section is None:
            raise NetFormatError(
                f"line {line_number}: expected net section header before records"
            )
        if indent < 1:
            raise NetFormatError(
                f"line {line_number}: net fill/gap record must be indented"
            )
        if previous_depth is None:
            if indent != 1:
                raise NetFormatError(
                    f"line {line_number}: first record in a net section must have depth 1"
                )
        elif indent > previous_depth + 1:
            raise NetFormatError(
                f"line {line_number}: net hierarchy depth may increase by only one level"
            )

        try:
            record = _parse_record(fields, section, indent)
        except ValueError as exc:
            raise NetFormatError(f"line {line_number}: {exc}") from exc

        previous_depth = indent
        yield record


def _parse_section_header(fields: list[str], line_number: int) -> _NetSection:
    if len(fields) != 3:
        raise NetFormatError(
            f"line {line_number}: net section header must contain exactly 3 fields"
        )
    _, target_name, target_size_text = fields
    if not target_name:
        raise NetFormatError(f"line {line_number}: net target name must not be empty")
    try:
        target_size = int(target_size_text)
    except ValueError as exc:
        raise NetFormatError(
            f"line {line_number}: net target sequence size must be an integer"
        ) from exc
    if target_size < 0:
        raise NetFormatError(
            f"line {line_number}: net target sequence size must be non-negative"
        )
    return _NetSection(target_name=target_name, target_sequence_size=target_size)


def _parse_record(fields: list[str], section: _NetSection, depth: int) -> NetRecord:
    if len(fields) < 7:
        raise ValueError("net record must contain at least 7 fixed fields")
    if (len(fields) - 7) % 2:
        raise ValueError("net optional fields must be name/value pairs")

    kind = _parse_kind(fields[0])
    target_start = _parse_nonnegative_int(fields[1], "net target start")
    target_span_size = _parse_nonnegative_int(fields[2], "net target span size")
    query_name = fields[3]
    if not query_name:
        raise ValueError("net query name must not be empty")
    orientation = _parse_orientation(fields[4])
    query_start = _parse_nonnegative_int(fields[5], "net query start")
    query_span_size = _parse_nonnegative_int(fields[6], "net query span size")

    attributes = tuple(zip(fields[7::2], fields[8::2], strict=True))
    attribute_by_name: dict[str, str] = {}
    for name, value in attributes:
        if not name:
            raise ValueError("net optional field name must not be empty")
        if name in attribute_by_name:
            raise ValueError(f"duplicate net optional field: {name}")
        attribute_by_name[name] = value

    chain_id = _optional_nonnegative_int(attribute_by_name, "id", "net chain ID")
    score = _optional_nonnegative_int(attribute_by_name, "score", "net score")
    aligned_bases = _optional_nonnegative_int(
        attribute_by_name, "ali", "net aligned bases"
    )
    duplicated_query_bases = _optional_nonnegative_int(
        attribute_by_name, "qDup", "net duplicated query bases"
    )
    classification = _optional_classification(attribute_by_name.get("type"))

    return NetRecord(
        target_name=section.target_name,
        target_sequence_size=section.target_sequence_size,
        depth=depth,
        kind=kind,
        target_start=target_start,
        target_span_size=target_span_size,
        query_name=query_name,
        orientation=orientation,
        query_start=query_start,
        query_span_size=query_span_size,
        chain_id=chain_id,
        score=score,
        aligned_bases=aligned_bases,
        duplicated_query_bases=duplicated_query_bases,
        classification=classification,
        attributes=attributes,
    )


def _parse_kind(value: str) -> NetRecordKind:
    try:
        return NetRecordKind(value)
    except ValueError as exc:
        raise ValueError("net record class must be 'fill' or 'gap'") from exc


def _parse_orientation(value: str) -> MappingOrientation:
    if value == "+":
        return MappingOrientation.SAME
    if value == "-":
        return MappingOrientation.REVERSE
    raise ValueError("net relative orientation must be '+' or '-'")


def _parse_nonnegative_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{label} must be an integer") from exc
    if parsed < 0:
        raise ValueError(f"{label} must be non-negative")
    return parsed


def _optional_nonnegative_int(
    attributes: dict[str, str], name: str, label: str
) -> int | None:
    value = attributes.get(name)
    if value is None:
        return None
    return _parse_nonnegative_int(value, label)


def _optional_classification(value: str | None) -> NetClassification | None:
    if value is None:
        return None
    try:
        return NetClassification(value)
    except ValueError as exc:
        raise ValueError(
            "net classification must be one of top, syn, inv, or nonSyn"
        ) from exc
