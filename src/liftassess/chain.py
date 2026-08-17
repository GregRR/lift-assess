"""Minimal UCSC chain-format reader.

This module parses chain structure only. It does not generate normalized
candidates, annotate chains with nets, score evidence, or make assessments.
UCSC chain coordinates are kept in their native 0-based, half-open convention;
query coordinates can be exposed in forward-reference coordinates explicitly.
"""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from enum import Enum

from .models import MappingOrientation


class ChainFormatError(ValueError):
    """Raised when input does not conform to the chain structure we consume."""


_CHAIN_CANDIDATE_SEPARATOR = ":chain:"


class ChainStrand(str, Enum):
    """Strand token used by the UCSC chain format."""

    PLUS = "+"
    MINUS = "-"


def chain_candidate_id(mapping_source_id: str, chain_id: int) -> str:
    """Return the canonical identity for a candidate generated from one chain."""

    if not mapping_source_id:
        raise ValueError("mapping source ID must not be empty")
    if chain_id < 0:
        raise ValueError("chain ID must be non-negative")
    return f"{mapping_source_id}{_CHAIN_CANDIDATE_SEPARATOR}{chain_id}"


def chain_id_from_candidate_id(candidate_id: str) -> int | None:
    """Return the UCSC chain ID from a canonical chain candidate ID, if present."""

    mapping_source_id, separator, chain_text = candidate_id.rpartition(
        _CHAIN_CANDIDATE_SEPARATOR
    )
    if not separator or not mapping_source_id:
        return None
    try:
        chain_id = int(chain_text)
    except ValueError:
        return None
    if chain_id < 0 or chain_text != str(chain_id):
        return None
    return chain_id


@dataclass(frozen=True)
class ChainBlock:
    """One ungapped chain block and the gaps that follow it, if any."""

    size: int
    target_gap: int | None = None
    query_gap: int | None = None

    def __post_init__(self) -> None:
        if self.size <= 0:
            raise ValueError("chain block size must be positive")
        if (self.target_gap is None) != (self.query_gap is None):
            raise ValueError("chain block gaps must both be present or both be absent")
        if self.target_gap is not None and self.target_gap < 0:
            raise ValueError("chain target gap must be non-negative")
        if self.query_gap is not None and self.query_gap < 0:
            raise ValueError("chain query gap must be non-negative")

    @property
    def is_terminal(self) -> bool:
        """Return whether this is the final block in a chain."""

        return self.target_gap is None

    def gaps_after(self) -> tuple[int, int]:
        """Return the target/query gaps after a non-terminal block."""

        target_gap = self.target_gap
        query_gap = self.query_gap
        if target_gap is None or query_gap is None:
            raise ValueError("terminal chain block has no following gaps")
        return target_gap, query_gap


@dataclass(frozen=True)
class ChainAlignedBlock:
    """One ungapped chain block in forward-reference coordinates.

    The chain file stores negative-strand query coordinates in reverse-complement
    space. Exposing each aligned block through this type keeps that conversion in
    one place so projection and evidence code cannot drift into different
    off-by-one conventions.
    """

    target_start: int
    target_end: int
    query_forward_start: int
    query_forward_end: int
    orientation: MappingOrientation

    def __post_init__(self) -> None:
        target_length = self.target_end - self.target_start
        query_length = self.query_forward_end - self.query_forward_start
        if target_length <= 0 or query_length <= 0:
            raise ValueError("aligned chain block must span at least one base")
        if target_length != query_length:
            raise ValueError("aligned chain block target/query lengths must match")

    def query_interval_for_target_interval(
        self, target_start: int, target_end: int
    ) -> tuple[int, int]:
        """Map a contained target-side subinterval onto forward query coordinates."""

        if target_start < self.target_start or target_end > self.target_end:
            raise ValueError("target subinterval lies outside aligned chain block")
        if target_end <= target_start:
            raise ValueError("target subinterval must span at least one base")

        offset_start = target_start - self.target_start
        offset_end = target_end - self.target_start
        if self.orientation is MappingOrientation.SAME:
            return (
                self.query_forward_start + offset_start,
                self.query_forward_start + offset_end,
            )
        return (
            self.query_forward_end - offset_end,
            self.query_forward_end - offset_start,
        )


@dataclass(frozen=True)
class ChainRecord:
    """One UCSC chain record in native chain coordinates."""

    score: float
    target_name: str
    target_size: int
    target_strand: ChainStrand
    target_start: int
    target_end: int
    query_name: str
    query_size: int
    query_strand: ChainStrand
    query_start: int
    query_end: int
    chain_id: int
    blocks: tuple[ChainBlock, ...]

    def __post_init__(self) -> None:
        if not self.target_name or not self.query_name:
            raise ValueError("chain sequence names must not be empty")
        if self.target_strand is not ChainStrand.PLUS:
            raise ValueError("UCSC chain target strand must be '+'")
        _validate_span("target", self.target_size, self.target_start, self.target_end)
        _validate_span("query", self.query_size, self.query_start, self.query_end)
        if not self.blocks:
            raise ValueError("chain must contain at least one alignment block")
        if not self.blocks[-1].is_terminal:
            raise ValueError("last chain block must be terminal")
        if any(block.is_terminal for block in self.blocks[:-1]):
            raise ValueError("only the last chain block may be terminal")

        target_consumed = sum(block.size for block in self.blocks) + sum(
            block.target_gap or 0 for block in self.blocks[:-1]
        )
        query_consumed = sum(block.size for block in self.blocks) + sum(
            block.query_gap or 0 for block in self.blocks[:-1]
        )
        if target_consumed != self.target_end - self.target_start:
            raise ValueError("chain blocks do not span the target header interval")
        if query_consumed != self.query_end - self.query_start:
            raise ValueError("chain blocks do not span the query header interval")

    @property
    def orientation(self) -> MappingOrientation:
        """Return query orientation relative to the chain target."""

        if self.query_strand is ChainStrand.PLUS:
            return MappingOrientation.SAME
        return MappingOrientation.REVERSE

    def query_interval_to_forward(
        self, native_start: int, native_end: int
    ) -> tuple[int, int]:
        """Convert native chain-query coordinates to forward-reference coordinates."""

        if (
            native_start < 0
            or native_end < native_start
            or native_end > self.query_size
        ):
            raise ValueError("query interval is outside sequence bounds")
        if self.query_strand is ChainStrand.PLUS:
            return native_start, native_end
        return self.query_size - native_end, self.query_size - native_start

    @property
    def query_forward_start(self) -> int:
        """Return query start in forward-reference coordinates."""

        start, _ = self.query_interval_to_forward(self.query_start, self.query_end)
        return start

    @property
    def query_forward_end(self) -> int:
        """Return query end in forward-reference coordinates."""

        _, end = self.query_interval_to_forward(self.query_start, self.query_end)
        return end

    def iter_aligned_blocks(self) -> Iterator[ChainAlignedBlock]:
        """Yield exact ungapped blocks with query coordinates normalized forward.

        This iterator is the canonical block-coordinate traversal for consumers
        that need aligned geometry. Gap evidence still uses the raw block gaps,
        but projection and reciprocal-best matching should use these normalized
        blocks instead of reimplementing chain cursor arithmetic.
        """

        target_cursor = self.target_start
        query_cursor = self.query_start

        for block in self.blocks:
            target_end = target_cursor + block.size
            query_native_end = query_cursor + block.size
            query_forward_start, query_forward_end = self.query_interval_to_forward(
                query_cursor, query_native_end
            )
            yield ChainAlignedBlock(
                target_start=target_cursor,
                target_end=target_end,
                query_forward_start=query_forward_start,
                query_forward_end=query_forward_end,
                orientation=self.orientation,
            )

            if block.is_terminal:
                return
            target_gap, query_gap = block.gaps_after()
            target_cursor = target_end + target_gap
            query_cursor = query_native_end + query_gap


@dataclass(frozen=True)
class _ChainHeader:
    score: float
    target_name: str
    target_size: int
    target_strand: ChainStrand
    target_start: int
    target_end: int
    query_name: str
    query_size: int
    query_strand: ChainStrand
    query_start: int
    query_end: int
    chain_id: int


def iter_chain_records(lines: Iterable[str]) -> Iterator[ChainRecord]:
    """Yield UCSC chain records from an iterable of text lines.

    Records are parsed incrementally so callers do not need to load a complete
    chain file into memory. To accept UCSC metadata prefixes without weakening
    record parsing, ``#`` metadata/comment lines are ignored only while looking for
    the next chain header. A final record may end at EOF immediately after its
    terminal block; otherwise records are separated by a blank line as in the
    documented format.
    """

    iterator = iter(enumerate(lines, start=1))

    while True:
        header_item = _next_chain_header_candidate(iterator)
        if header_item is None:
            return
        header_line_number, header_line = header_item
        header = _parse_header(header_line, header_line_number)
        blocks: list[ChainBlock] = []

        while True:
            try:
                line_number, raw_line = next(iterator)
            except StopIteration as exc:
                raise ChainFormatError(
                    f"line {header_line_number}: chain ended before terminal block"
                ) from exc

            stripped = raw_line.strip()
            if not stripped:
                raise ChainFormatError(
                    f"line {line_number}: blank line before terminal chain block"
                )
            if stripped.startswith("#"):
                raise ChainFormatError(
                    f"line {line_number}: metadata/comment line inside chain record"
                )

            parts = stripped.split()
            if len(parts) == 3:
                blocks.append(
                    _make_block(
                        line_number=line_number,
                        size=_parse_int(parts[0], line_number, "block size"),
                        target_gap=_parse_int(parts[1], line_number, "target gap"),
                        query_gap=_parse_int(parts[2], line_number, "query gap"),
                    )
                )
                continue

            if len(parts) != 1:
                raise ChainFormatError(
                    f"line {line_number}: chain block must contain 1 or 3 fields"
                )

            blocks.append(
                _make_block(
                    line_number=line_number,
                    size=_parse_int(parts[0], line_number, "terminal block size"),
                )
            )
            record = _build_record(header, tuple(blocks), header_line_number)

            try:
                separator_line_number, separator = next(iterator)
            except StopIteration:
                yield record
                return

            if separator.strip():
                raise ChainFormatError(
                    f"line {separator_line_number}: expected blank line after chain record"
                )

            yield record
            break


def _next_chain_header_candidate(
    iterator: Iterator[tuple[int, str]],
) -> tuple[int, str] | None:
    """Return the next nonblank, non-metadata line outside a chain record."""

    for line_number, line in iterator:
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            return line_number, line
    return None


def _parse_header(line: str, line_number: int) -> _ChainHeader:
    parts = line.strip().split()
    if not parts or parts[0] != "chain":
        raise ChainFormatError(f"line {line_number}: expected chain header")
    if len(parts) != 13:
        raise ChainFormatError(
            f"line {line_number}: chain header must contain 13 fields including 'chain'"
        )

    try:
        target_strand = ChainStrand(parts[4])
        query_strand = ChainStrand(parts[9])
    except ValueError as exc:
        raise ChainFormatError(
            f"line {line_number}: chain strand must be '+' or '-'"
        ) from exc

    return _ChainHeader(
        score=_parse_float(parts[1], line_number, "score"),
        target_name=parts[2],
        target_size=_parse_int(parts[3], line_number, "target size"),
        target_strand=target_strand,
        target_start=_parse_int(parts[5], line_number, "target start"),
        target_end=_parse_int(parts[6], line_number, "target end"),
        query_name=parts[7],
        query_size=_parse_int(parts[8], line_number, "query size"),
        query_strand=query_strand,
        query_start=_parse_int(parts[10], line_number, "query start"),
        query_end=_parse_int(parts[11], line_number, "query end"),
        chain_id=_parse_int(parts[12], line_number, "chain id"),
    )


def _build_record(
    header: _ChainHeader,
    blocks: tuple[ChainBlock, ...],
    header_line_number: int,
) -> ChainRecord:
    try:
        return ChainRecord(
            score=header.score,
            target_name=header.target_name,
            target_size=header.target_size,
            target_strand=header.target_strand,
            target_start=header.target_start,
            target_end=header.target_end,
            query_name=header.query_name,
            query_size=header.query_size,
            query_strand=header.query_strand,
            query_start=header.query_start,
            query_end=header.query_end,
            chain_id=header.chain_id,
            blocks=blocks,
        )
    except ValueError as exc:
        raise ChainFormatError(f"line {header_line_number}: {exc}") from exc


def _make_block(
    *,
    line_number: int,
    size: int,
    target_gap: int | None = None,
    query_gap: int | None = None,
) -> ChainBlock:
    try:
        return ChainBlock(size=size, target_gap=target_gap, query_gap=query_gap)
    except ValueError as exc:
        raise ChainFormatError(f"line {line_number}: {exc}") from exc


def _parse_float(value: str, line_number: int, field_name: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ChainFormatError(
            f"line {line_number}: {field_name} must be numeric"
        ) from exc


def _parse_int(value: str, line_number: int, field_name: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ChainFormatError(
            f"line {line_number}: {field_name} must be an integer"
        ) from exc


def _validate_span(name: str, size: int, start: int, end: int) -> None:
    if size < 0:
        raise ValueError(f"chain {name} size must be non-negative")
    if start < 0 or end <= start or end > size:
        raise ValueError(f"chain {name} interval is outside sequence bounds")
