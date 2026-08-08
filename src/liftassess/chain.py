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


class ChainStrand(str, Enum):
    """Strand token used by the UCSC chain format."""

    PLUS = "+"
    MINUS = "-"


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

    @property
    def query_forward_start(self) -> int:
        """Return query start in forward-reference coordinates."""

        if self.query_strand is ChainStrand.PLUS:
            return self.query_start
        return self.query_size - self.query_end

    @property
    def query_forward_end(self) -> int:
        """Return query end in forward-reference coordinates."""

        if self.query_strand is ChainStrand.PLUS:
            return self.query_end
        return self.query_size - self.query_start


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
    chain file into memory. A final record may end at EOF immediately after its
    terminal block; otherwise records are separated by a blank line as in the
    documented format.
    """

    iterator = iter(enumerate(lines, start=1))

    while True:
        header_item = _next_nonblank(iterator)
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


def _next_nonblank(
    iterator: Iterator[tuple[int, str]],
) -> tuple[int, str] | None:
    for line_number, line in iterator:
        if line.strip():
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
