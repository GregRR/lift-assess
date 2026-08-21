"""Benchmark-only selective chain traversal retained from Milestone 18 prototyping.

Production liftAssess uses the ordinary full traversal as its unindexed fallback and the
reusable chain index for region-addressable access.  This helper remains under
``scripts/`` only so the historical selective-materialization benchmark can be rerun
without maintaining a second traversal implementation in the installed package.
"""

from __future__ import annotations

from collections.abc import Iterator

from liftassess.chain import (
    ChainBlock,
    ChainFormatError,
    ChainRecord,
    ChainStrand,
    _build_record,
    _ChainHeader,
    _make_block,
    _next_chain_header_candidate,
    _parse_header,
    _parse_int,
    _validate_span,
)
from liftassess.resource_stream import ResourcePath, open_text_resource


def iter_chain_file_overlapping_target_interval(
    path: ResourcePath,
    *,
    target_name: str,
    target_start: int,
    target_end: int,
) -> Iterator[ChainRecord]:
    """Yield selectively materialized chains for benchmark comparison only."""

    with open_text_resource(path) as lines:
        yield from _iter_chain_records_overlapping_target_interval(
            lines,
            target_name=target_name,
            target_start=target_start,
            target_end=target_end,
        )


def _iter_chain_records_overlapping_target_interval(
    lines: Iterator[str],
    *,
    target_name: str,
    target_start: int,
    target_end: int,
) -> Iterator[ChainRecord]:
    if not target_name:
        raise ValueError("target name must not be empty")
    if target_start < 0 or target_end <= target_start:
        raise ValueError("target interval must be non-empty and non-negative")

    iterator = iter(enumerate(lines, start=1))
    while True:
        header_item = _next_chain_header_candidate(iterator)
        if header_item is None:
            return
        header_line_number, header_line = header_item
        header = _parse_header(header_line, header_line_number)
        materialize = (
            header.target_name == target_name
            and target_start < header.target_end
            and target_end > header.target_start
        )

        if materialize:
            record, reached_eof = _parse_selected_chain_record(
                iterator,
                header=header,
                header_line_number=header_line_number,
            )
        else:
            record = None
            reached_eof = _validate_and_skip_chain_record(
                iterator,
                header=header,
                header_line_number=header_line_number,
            )

        if header.target_name == target_name and target_end > header.target_size:
            raise ValueError("source interval exceeds chain target sequence bounds")
        if record is not None:
            yield record
        if reached_eof:
            return


def _parse_selected_chain_record(
    iterator: Iterator[tuple[int, str]],
    *,
    header: _ChainHeader,
    header_line_number: int,
) -> tuple[ChainRecord, bool]:
    blocks: list[ChainBlock] = []
    while True:
        line_number, stripped = _next_selective_block_line(
            iterator, header_line_number=header_line_number
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
        return record, _consume_selective_separator(iterator)


def _validate_and_skip_chain_record(
    iterator: Iterator[tuple[int, str]],
    *,
    header: _ChainHeader,
    header_line_number: int,
) -> bool:
    target_consumed = 0
    query_consumed = 0
    while True:
        line_number, stripped = _next_selective_block_line(
            iterator, header_line_number=header_line_number
        )
        parts = stripped.split()
        if len(parts) == 3:
            size = _parse_int(parts[0], line_number, "block size")
            target_gap = _parse_int(parts[1], line_number, "target gap")
            query_gap = _parse_int(parts[2], line_number, "query gap")
            _validate_skipped_block_values(
                line_number=line_number,
                size=size,
                target_gap=target_gap,
                query_gap=query_gap,
            )
            target_consumed += size + target_gap
            query_consumed += size + query_gap
            continue
        if len(parts) != 1:
            raise ChainFormatError(
                f"line {line_number}: chain block must contain 1 or 3 fields"
            )

        size = _parse_int(parts[0], line_number, "terminal block size")
        _validate_skipped_block_values(line_number=line_number, size=size)
        target_consumed += size
        query_consumed += size
        _validate_skipped_record_summary(
            header,
            header_line_number=header_line_number,
            target_consumed=target_consumed,
            query_consumed=query_consumed,
        )
        return _consume_selective_separator(iterator)


def _next_selective_block_line(
    iterator: Iterator[tuple[int, str]],
    *,
    header_line_number: int,
) -> tuple[int, str]:
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
    return line_number, stripped


def _consume_selective_separator(iterator: Iterator[tuple[int, str]]) -> bool:
    try:
        separator_line_number, separator = next(iterator)
    except StopIteration:
        return True
    if separator.strip():
        raise ChainFormatError(
            f"line {separator_line_number}: expected blank line after chain record"
        )
    return False


def _validate_skipped_block_values(
    *,
    line_number: int,
    size: int,
    target_gap: int | None = None,
    query_gap: int | None = None,
) -> None:
    if size <= 0:
        raise ChainFormatError(f"line {line_number}: chain block size must be positive")
    if (target_gap is None) != (query_gap is None):
        raise ChainFormatError(
            f"line {line_number}: chain block gaps must both be present or both be absent"
        )
    if target_gap is not None and target_gap < 0:
        raise ChainFormatError(
            f"line {line_number}: chain target gap must be non-negative"
        )
    if query_gap is not None and query_gap < 0:
        raise ChainFormatError(
            f"line {line_number}: chain query gap must be non-negative"
        )


def _validate_skipped_record_summary(
    header: _ChainHeader,
    *,
    header_line_number: int,
    target_consumed: int,
    query_consumed: int,
) -> None:
    try:
        if not header.target_name or not header.query_name:
            raise ValueError("chain sequence names must not be empty")
        if header.target_strand is not ChainStrand.PLUS:
            raise ValueError("UCSC chain target strand must be '+'")
        _validate_span(
            "target", header.target_size, header.target_start, header.target_end
        )
        _validate_span("query", header.query_size, header.query_start, header.query_end)
        if target_consumed != header.target_end - header.target_start:
            raise ValueError("chain blocks do not span the target header interval")
        if query_consumed != header.query_end - header.query_start:
            raise ValueError("chain blocks do not span the query header interval")
    except ValueError as exc:
        raise ChainFormatError(f"line {header_line_number}: {exc}") from exc
