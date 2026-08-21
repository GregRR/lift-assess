#!/usr/bin/env python3
"""Benchmark a full resource-level 64-KiB blocked chain index in one pass.

This experiment is the full-scale counterpart to the targeted Milestone-18
sequence/bin/block-store prototypes. It scans the original cached UCSC chain exactly
once, simultaneously records baseline candidates for requested loci, and builds one
resource-level index containing:

- a source-sequence table;
- 64-KiB genomic-bin memberships;
- exact chain target intervals and encounter-order record IDs; and
- every serialized chain record exactly once in independently zlib-compressed blocks.

After the build, each requested locus is queried through the index and its candidate
tuple must exactly match the candidates recorded during the original full pass.
"""

from __future__ import annotations

import argparse
import io
import json
import shutil
import sqlite3
import time
import zlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from liftassess.chain import ChainRecord, iter_chain_records
from liftassess.models import (
    AssemblyIdentifier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.projection import project_interval_through_chain
from liftassess.resource_cache import load_cached_ucsc_resource_bundle
from liftassess.resource_files import iter_chain_file

CandidateTuple = tuple[NormalizedCandidate, ...]
_INDEX_FORMAT = "liftassess-chain-full-block-index-v1"
_DEFAULT_BIN_WIDTH = 65_536
_DEFAULT_BLOCK_SIZE = 1024 * 1024
_CHAIN_BATCH_SIZE = 10_000
_MEMBERSHIP_BATCH_SIZE = 100_000


@dataclass(frozen=True)
class _RequestedLocus:
    text: str
    interval: GenomicInterval


@dataclass(frozen=True)
class _IndexPaths:
    database: Path
    blocks: Path
    manifest: Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bin-width", type=int, default=_DEFAULT_BIN_WIDTH)
    parser.add_argument("--block-size", type=int, default=_DEFAULT_BLOCK_SIZE)
    parser.add_argument("source_db")
    parser.add_argument("target_db")
    parser.add_argument(
        "locus",
        nargs="+",
        help="one or more 1-based inclusive source loci, chr:start-end",
    )
    return parser.parse_args()


def _parse_locus(text: str, assembly: AssemblyIdentifier) -> _RequestedLocus:
    try:
        sequence_name, coordinates = text.rsplit(":", 1)
        start_text, end_text = coordinates.split("-", 1)
        start = int(start_text)
        end = int(end_text)
    except ValueError as exc:
        raise ValueError(
            "locus must use chr:start-end with integer coordinates"
        ) from exc
    if start <= 0 or end < start:
        raise ValueError("locus must be non-empty and 1-based inclusive")
    return _RequestedLocus(
        text=text,
        interval=GenomicInterval(
            assembly=assembly,
            sequence_name=sequence_name,
            start=start - 1,
            end=end,
        ),
    )


def _format_chain_record(record: ChainRecord) -> bytes:
    header = " ".join(
        (
            "chain",
            repr(record.score),
            record.target_name,
            str(record.target_size),
            record.target_strand.value,
            str(record.target_start),
            str(record.target_end),
            record.query_name,
            str(record.query_size),
            record.query_strand.value,
            str(record.query_start),
            str(record.query_end),
            str(record.chain_id),
        )
    )
    lines = [header]
    for block in record.blocks:
        if block.is_terminal:
            lines.append(str(block.size))
        else:
            target_gap, query_gap = block.gaps_after()
            lines.append(f"{block.size}\t{target_gap}\t{query_gap}")
    return ("\n".join(lines) + "\n\n").encode("ascii")


def _paths(output_dir: Path) -> _IndexPaths:
    return _IndexPaths(
        database=output_dir / "index.sqlite3",
        blocks=output_dir / "records.blocks",
        manifest=output_dir / "manifest.json",
    )


def _configure_database(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("PRAGMA temp_store = MEMORY")
    connection.executescript(
        """
        CREATE TABLE sequences (
            sequence_id INTEGER PRIMARY KEY,
            name TEXT NOT NULL UNIQUE,
            minimum_target_size INTEGER NOT NULL
        );
        CREATE TABLE chains (
            record_id INTEGER PRIMARY KEY,
            sequence_id INTEGER NOT NULL,
            target_start INTEGER NOT NULL,
            target_end INTEGER NOT NULL,
            block_id INTEGER NOT NULL,
            block_offset INTEGER NOT NULL,
            record_length INTEGER NOT NULL
        );
        CREATE TABLE blocks (
            block_id INTEGER PRIMARY KEY,
            file_offset INTEGER NOT NULL,
            compressed_length INTEGER NOT NULL,
            uncompressed_length INTEGER NOT NULL
        );
        CREATE TABLE bin_memberships (
            sequence_id INTEGER NOT NULL,
            bin_id INTEGER NOT NULL,
            record_id INTEGER NOT NULL,
            PRIMARY KEY (sequence_id, bin_id, record_id)
        ) WITHOUT ROWID;
        """
    )


def _write_compressed_block(
    block_store: io.BufferedWriter,
    connection: sqlite3.Connection,
    *,
    block_id: int,
    payload: bytes,
) -> None:
    compressed = zlib.compress(payload)
    file_offset = block_store.tell()
    block_store.write(compressed)
    connection.execute(
        "INSERT INTO blocks VALUES (?, ?, ?, ?)",
        (block_id, file_offset, len(compressed), len(payload)),
    )


def _flush_rows(
    connection: sqlite3.Connection,
    chain_rows: list[tuple[int, int, int, int, int, int, int]],
    membership_rows: list[tuple[int, int, int]],
) -> None:
    if chain_rows:
        connection.executemany(
            "INSERT INTO chains VALUES (?, ?, ?, ?, ?, ?, ?)",
            chain_rows,
        )
        chain_rows.clear()
    if membership_rows:
        connection.executemany(
            "INSERT INTO bin_memberships VALUES (?, ?, ?)",
            membership_rows,
        )
        membership_rows.clear()


def _build_index(
    *,
    chain_path: Path,
    output_dir: Path,
    bin_width: int,
    block_size: int,
    loci: tuple[_RequestedLocus, ...],
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[
    _IndexPaths,
    dict[str, CandidateTuple],
    int,
    int,
    int,
    float,
]:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    paths = _paths(output_dir)

    loci_by_sequence: dict[str, list[_RequestedLocus]] = defaultdict(list)
    for locus in loci:
        loci_by_sequence[locus.interval.sequence_name].append(locus)
    baseline_lists: dict[str, list[NormalizedCandidate]] = defaultdict(list)

    sequence_ids: dict[str, int] = {}
    sequence_min_sizes: dict[str, int] = {}
    chain_rows: list[tuple[int, int, int, int, int, int, int]] = []
    membership_rows: list[tuple[int, int, int]] = []
    block_payload = bytearray()
    block_id = 0
    block_count = 0
    record_count = 0
    membership_count = 0

    started = time.perf_counter()
    with (
        sqlite3.connect(paths.database) as connection,
        paths.blocks.open("wb") as block_store,
    ):
        _configure_database(connection)
        for record_id, chain in enumerate(iter_chain_file(chain_path)):
            sequence_id = sequence_ids.get(chain.target_name)
            if sequence_id is None:
                sequence_id = len(sequence_ids)
                sequence_ids[chain.target_name] = sequence_id
                sequence_min_sizes[chain.target_name] = chain.target_size
                connection.execute(
                    "INSERT INTO sequences VALUES (?, ?, ?)",
                    (sequence_id, chain.target_name, chain.target_size),
                )
            elif chain.target_size < sequence_min_sizes[chain.target_name]:
                sequence_min_sizes[chain.target_name] = chain.target_size
                connection.execute(
                    "UPDATE sequences SET minimum_target_size = ? "
                    "WHERE sequence_id = ?",
                    (chain.target_size, sequence_id),
                )

            payload = _format_chain_record(chain)
            if block_payload and len(block_payload) + len(payload) > block_size:
                _write_compressed_block(
                    block_store,
                    connection,
                    block_id=block_id,
                    payload=bytes(block_payload),
                )
                block_count += 1
                block_id += 1
                block_payload.clear()

            block_offset = len(block_payload)
            block_payload.extend(payload)
            chain_rows.append(
                (
                    record_id,
                    sequence_id,
                    chain.target_start,
                    chain.target_end,
                    block_id,
                    block_offset,
                    len(payload),
                )
            )

            first_bin = chain.target_start // bin_width
            last_bin = (chain.target_end - 1) // bin_width
            for bin_id in range(first_bin, last_bin + 1):
                membership_rows.append((sequence_id, bin_id, record_id))
                membership_count += 1

            for locus in loci_by_sequence.get(chain.target_name, ()):
                candidate = project_interval_through_chain(
                    locus.interval,
                    chain,
                    target_assembly=target_assembly,
                    mapping_provenance=provenance,
                )
                if candidate is not None:
                    baseline_lists[locus.text].append(candidate)

            record_count += 1
            if (
                len(chain_rows) >= _CHAIN_BATCH_SIZE
                or len(membership_rows) >= _MEMBERSHIP_BATCH_SIZE
            ):
                _flush_rows(connection, chain_rows, membership_rows)

        if block_payload:
            _write_compressed_block(
                block_store,
                connection,
                block_id=block_id,
                payload=bytes(block_payload),
            )
            block_count += 1
        _flush_rows(connection, chain_rows, membership_rows)
        connection.commit()
        connection.execute("ANALYZE")
        connection.commit()

    baselines = {locus.text: tuple(baseline_lists[locus.text]) for locus in loci}
    return (
        paths,
        baselines,
        record_count,
        membership_count,
        block_count,
        time.perf_counter() - started,
    )


def _indexed_records(
    paths: _IndexPaths,
    interval: GenomicInterval,
    *,
    bin_width: int,
) -> tuple[tuple[ChainRecord, ...], int, int, int]:
    first_bin = interval.start // bin_width
    last_bin = (interval.end - 1) // bin_width
    with sqlite3.connect(paths.database) as connection:
        sequence_row = connection.execute(
            "SELECT sequence_id, minimum_target_size FROM sequences WHERE name = ?",
            (interval.sequence_name,),
        ).fetchone()
        if sequence_row is None:
            return (), 0, last_bin - first_bin + 1, 0
        sequence_id, minimum_target_size = sequence_row
        if interval.end > minimum_target_size:
            raise ValueError("source interval exceeds chain target sequence bounds")

        rows = tuple(
            connection.execute(
                """
                SELECT c.record_id, c.block_id, c.block_offset, c.record_length,
                       b.file_offset, b.compressed_length, b.uncompressed_length
                FROM chains AS c
                JOIN blocks AS b ON b.block_id = c.block_id
                WHERE c.record_id IN (
                    SELECT m.record_id
                    FROM bin_memberships AS m
                    WHERE m.sequence_id = ? AND m.bin_id BETWEEN ? AND ?
                )
                  AND c.sequence_id = ?
                  AND c.target_start < ?
                  AND c.target_end > ?
                ORDER BY c.record_id
                """,
                (
                    sequence_id,
                    first_bin,
                    last_bin,
                    sequence_id,
                    interval.end,
                    interval.start,
                ),
            )
        )

    payload = bytearray()
    blocks_read = 0
    current_block_id: int | None = None
    current_block_payload = b""
    with paths.blocks.open("rb") as block_store:
        for (
            _record_id,
            selected_block_id,
            block_offset,
            record_length,
            file_offset,
            compressed_length,
            uncompressed_length,
        ) in rows:
            if selected_block_id != current_block_id:
                block_store.seek(file_offset)
                compressed = block_store.read(compressed_length)
                if len(compressed) != compressed_length:
                    raise SystemExit("compressed block store ended unexpectedly")
                current_block_payload = zlib.decompress(compressed)
                if len(current_block_payload) != uncompressed_length:
                    raise SystemExit("compressed block length changed during decoding")
                current_block_id = selected_block_id
                blocks_read += 1
            record = current_block_payload[block_offset : block_offset + record_length]
            if len(record) != record_length:
                raise SystemExit("record slice exceeds decoded compression block")
            payload.extend(record)

    records = tuple(iter_chain_records(io.StringIO(payload.decode("ascii"))))
    if len(records) != len(rows):
        raise SystemExit("indexed record count changed during parsing")
    return records, len(rows), last_bin - first_bin + 1, blocks_read


def _indexed_candidates(
    paths: _IndexPaths,
    locus: _RequestedLocus,
    *,
    bin_width: int,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, int, int, int, float]:
    started = time.perf_counter()
    records, selected_records, bins_touched, blocks_read = _indexed_records(
        paths,
        locus.interval,
        bin_width=bin_width,
    )
    candidates = tuple(
        candidate
        for chain in records
        if (
            candidate := project_interval_through_chain(
                locus.interval,
                chain,
                target_assembly=target_assembly,
                mapping_provenance=provenance,
            )
        )
        is not None
    )
    return (
        candidates,
        selected_records,
        bins_touched,
        blocks_read,
        time.perf_counter() - started,
    )


def _write_manifest(
    paths: _IndexPaths,
    *,
    source_db: str,
    target_db: str,
    source_chain_sha256_identifier: str,
    source_chain_size_bytes: int,
    bin_width: int,
    block_size: int,
    record_count: int,
    membership_count: int,
    block_count: int,
    build_seconds: float,
) -> None:
    payload = {
        "format": _INDEX_FORMAT,
        "source_db": source_db,
        "target_db": target_db,
        "source_chain_sha256_identifier": source_chain_sha256_identifier,
        "source_chain_size_bytes": source_chain_size_bytes,
        "bin_width": bin_width,
        "block_size": block_size,
        "record_count": record_count,
        "membership_count": membership_count,
        "block_count": block_count,
        "build_seconds": build_seconds,
        "database_size_bytes": paths.database.stat().st_size,
        "block_store_size_bytes": paths.blocks.stat().st_size,
    }
    paths.manifest.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    args = _parse_args()
    if args.bin_width <= 0:
        raise SystemExit("bin width must be positive")
    if args.block_size <= 0:
        raise SystemExit("block size must be positive")

    bundle = load_cached_ucsc_resource_bundle(
        args.cache_root,
        args.source_db,
        args.target_db,
    )
    if bundle is None:
        raise SystemExit("no complete verified cached UCSC bundle was found")

    source_assembly = AssemblyIdentifier(args.source_db, "UCSC")
    target_assembly = AssemblyIdentifier(args.target_db, "UCSC")
    loci = tuple(_parse_locus(text, source_assembly) for text in args.locus)
    provenance = ProvenanceSource("benchmark-chain", "benchmark chain resource")

    print(f"chain: {bundle.chain.path}")
    print(f"chain SHA-256: {bundle.chain.sha256}")
    print(f"bin width: {args.bin_width} bases")
    print(f"compression block target: {args.block_size} bytes")
    print("building one full resource-level index and recording baseline candidates...")

    (
        paths,
        baselines,
        record_count,
        membership_count,
        block_count,
        build_seconds,
    ) = _build_index(
        chain_path=bundle.chain.path,
        output_dir=args.output_dir,
        bin_width=args.bin_width,
        block_size=args.block_size,
        loci=loci,
        target_assembly=target_assembly,
        provenance=provenance,
    )
    _write_manifest(
        paths,
        source_db=args.source_db,
        target_db=args.target_db,
        source_chain_sha256_identifier=bundle.chain.sha256,
        source_chain_size_bytes=bundle.chain.size_bytes,
        bin_width=args.bin_width,
        block_size=args.block_size,
        record_count=record_count,
        membership_count=membership_count,
        block_count=block_count,
        build_seconds=build_seconds,
    )

    database_size = paths.database.stat().st_size
    block_store_size = paths.blocks.stat().st_size
    print(
        f"full index build: records={record_count} memberships={membership_count} "
        f"blocks={block_count} elapsed={build_seconds:.3f}s"
    )
    print(
        f"full index size: database={database_size / (1024**3):.3f} GiB "
        f"blocks={block_store_size / (1024**3):.3f} GiB "
        f"total={(database_size + block_store_size) / (1024**3):.3f} GiB"
    )

    for locus in loci:
        baseline = baselines[locus.text]
        (
            indexed,
            selected_records,
            bins_touched,
            blocks_read,
            indexed_elapsed,
        ) = _indexed_candidates(
            paths,
            locus,
            bin_width=args.bin_width,
            target_assembly=target_assembly,
            provenance=provenance,
        )
        print()
        print(f"source locus: {args.source_db} {locus.text}")
        print(f"full-pass baseline: candidates={len(baseline)}")
        print(
            f"full-index query: candidates={len(indexed)} "
            f"elapsed={indexed_elapsed:.3f}s bins={bins_touched} "
            f"selected_records={selected_records} blocks_read={blocks_read}"
        )
        if indexed != baseline:
            raise SystemExit(
                f"candidate mismatch for {locus.text}: "
                f"baseline={len(baseline)} indexed={len(indexed)}"
            )
        print(f"equivalence: PASS ({len(indexed)} candidates)")


if __name__ == "__main__":
    main()
