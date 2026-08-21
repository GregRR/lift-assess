#!/usr/bin/env python3
"""Prototype a compact blocked-compression record store for the 64-KiB chain index.

This benchmark reuses the exact SQLite bin indexes and uncompressed record stores
created by ``benchmark_chain_bin_index.py``. It does not rescan the original UCSC
chain or the sequence shards. Each chain record remains stored exactly once, but the
record stream is repacked into independently compressed encounter-order blocks.

The compact SQLite database preserves the existing 64-KiB bin memberships and exact
target intervals while replacing raw record offsets with block + within-block offsets.
A query reads only the compressed blocks containing exact interval-overlapping records,
restores encounter order, parses with the production chain parser, and projects with
the production projection function. Candidate tuples are compared exactly against the
existing SQLite + uncompressed-record-store prototype.
"""

from __future__ import annotations

import argparse
import io
import json
import shutil
import sqlite3
import time
import zlib
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

CandidateTuple = tuple[NormalizedCandidate, ...]
_INDEX_FORMAT = "liftassess-chain-bin-index-v1"
_COMPACT_FORMAT = "liftassess-chain-block-store-v1"
_DEFAULT_BLOCK_SIZE = 1024 * 1024
_INSERT_BATCH_SIZE = 10_000


@dataclass(frozen=True)
class _RequestedLocus:
    text: str
    interval: GenomicInterval


@dataclass(frozen=True)
class _IndexPaths:
    database: Path
    records: Path


@dataclass(frozen=True)
class _CompactPaths:
    database: Path
    blocks: Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
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


def _sequence_key(sequence_name: str) -> str:
    return sequence_name.encode("utf-8").hex()


def _index_paths(index_dir: Path, sequence_name: str) -> _IndexPaths:
    key = _sequence_key(sequence_name)
    return _IndexPaths(
        database=index_dir / f"seq-{key}.sqlite3",
        records=index_dir / f"seq-{key}.records",
    )


def _compact_paths(output_dir: Path, sequence_name: str) -> _CompactPaths:
    key = _sequence_key(sequence_name)
    return _CompactPaths(
        database=output_dir / f"seq-{key}.sqlite3",
        blocks=output_dir / f"seq-{key}.blocks",
    )


def _load_index_manifest(index_dir: Path) -> dict[str, object]:
    path = index_dir / "manifest.json"
    if not path.is_file():
        raise SystemExit(f"bin-index manifest is missing: {path}")
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("format") != _INDEX_FORMAT:
        raise SystemExit("bin-index manifest has an unsupported format")
    return raw


def _validate_index_paths(paths: _IndexPaths, sequence_name: str) -> None:
    if not paths.database.is_file() or not paths.records.is_file():
        raise SystemExit(f"bin index is incomplete for sequence {sequence_name}")


def _configure_compact_database(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("PRAGMA temp_store = MEMORY")
    connection.executescript(
        """
        CREATE TABLE chains (
            record_id INTEGER PRIMARY KEY,
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
            bin_id INTEGER NOT NULL,
            record_id INTEGER NOT NULL,
            PRIMARY KEY (bin_id, record_id)
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


def _build_compact_sequence_store(
    *,
    source_paths: _IndexPaths,
    output_dir: Path,
    sequence_name: str,
    block_size: int,
) -> tuple[_CompactPaths, int, int, float]:
    compact_paths = _compact_paths(output_dir, sequence_name)
    compact_paths.database.unlink(missing_ok=True)
    compact_paths.blocks.unlink(missing_ok=True)

    started = time.perf_counter()
    block_id = 0
    block_payload = bytearray()
    chain_rows: list[tuple[int, int, int, int, int, int]] = []
    record_count = 0
    block_count = 0

    with (
        sqlite3.connect(source_paths.database) as source_connection,
        sqlite3.connect(compact_paths.database) as compact_connection,
        source_paths.records.open("rb") as source_records,
        compact_paths.blocks.open("wb") as block_store,
    ):
        _configure_compact_database(compact_connection)
        source_rows = source_connection.execute(
            """
            SELECT record_id, target_start, target_end, record_offset, record_length
            FROM chains
            ORDER BY record_id
            """
        )
        for (
            record_id,
            target_start,
            target_end,
            record_offset,
            record_length,
        ) in source_rows:
            source_records.seek(record_offset)
            payload = source_records.read(record_length)
            if len(payload) != record_length:
                raise SystemExit("source record store ended unexpectedly")

            if block_payload and len(block_payload) + record_length > block_size:
                _write_compressed_block(
                    block_store,
                    compact_connection,
                    block_id=block_id,
                    payload=bytes(block_payload),
                )
                block_count += 1
                block_id += 1
                block_payload.clear()

            inner_offset = len(block_payload)
            block_payload.extend(payload)
            chain_rows.append(
                (
                    record_id,
                    target_start,
                    target_end,
                    block_id,
                    inner_offset,
                    record_length,
                )
            )
            record_count += 1

            if len(chain_rows) >= _INSERT_BATCH_SIZE:
                compact_connection.executemany(
                    "INSERT INTO chains VALUES (?, ?, ?, ?, ?, ?)",
                    chain_rows,
                )
                chain_rows.clear()

        if block_payload:
            _write_compressed_block(
                block_store,
                compact_connection,
                block_id=block_id,
                payload=bytes(block_payload),
            )
            block_count += 1
        if chain_rows:
            compact_connection.executemany(
                "INSERT INTO chains VALUES (?, ?, ?, ?, ?, ?)",
                chain_rows,
            )

        compact_connection.execute(
            "ATTACH DATABASE ? AS source_index", (str(source_paths.database),)
        )
        compact_connection.execute(
            "INSERT INTO bin_memberships "
            "SELECT bin_id, record_id FROM source_index.bin_memberships"
        )
        compact_connection.commit()
        compact_connection.execute("DETACH DATABASE source_index")
        compact_connection.execute("ANALYZE")
        compact_connection.commit()

    return (
        compact_paths,
        record_count,
        block_count,
        time.perf_counter() - started,
    )


def _sqlite_records(
    paths: _IndexPaths,
    *,
    target_start: int,
    target_end: int,
    bin_width: int,
) -> tuple[tuple[ChainRecord, ...], int]:
    first_bin = target_start // bin_width
    last_bin = (target_end - 1) // bin_width
    with sqlite3.connect(paths.database) as connection:
        rows = tuple(
            connection.execute(
                """
                SELECT c.record_id, c.record_offset, c.record_length
                FROM chains AS c
                WHERE c.record_id IN (
                    SELECT b.record_id
                    FROM bin_memberships AS b
                    WHERE b.bin_id BETWEEN ? AND ?
                )
                  AND c.target_start < ?
                  AND c.target_end > ?
                ORDER BY c.record_id
                """,
                (first_bin, last_bin, target_end, target_start),
            )
        )

    payload = bytearray()
    with paths.records.open("rb") as record_store:
        for _record_id, offset, length in rows:
            record_store.seek(offset)
            chunk = record_store.read(length)
            if len(chunk) != length:
                raise SystemExit("indexed chain record store ended unexpectedly")
            payload.extend(chunk)

    records = tuple(iter_chain_records(io.StringIO(payload.decode("ascii"))))
    if len(records) != len(rows):
        raise SystemExit("SQLite baseline record count changed during parsing")
    return records, len(rows)


def _compact_records(
    paths: _CompactPaths,
    *,
    target_start: int,
    target_end: int,
    bin_width: int,
) -> tuple[tuple[ChainRecord, ...], int, int, int]:
    first_bin = target_start // bin_width
    last_bin = (target_end - 1) // bin_width
    with sqlite3.connect(paths.database) as connection:
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
                    WHERE m.bin_id BETWEEN ? AND ?
                )
                  AND c.target_start < ?
                  AND c.target_end > ?
                ORDER BY c.record_id
                """,
                (first_bin, last_bin, target_end, target_start),
            )
        )

    payload = bytearray()
    blocks_read = 0
    compressed_bytes_read = 0
    current_block_id: int | None = None
    current_block_payload = b""
    with paths.blocks.open("rb") as block_store:
        for (
            _record_id,
            block_id,
            block_offset,
            record_length,
            file_offset,
            compressed_length,
            uncompressed_length,
        ) in rows:
            if block_id != current_block_id:
                block_store.seek(file_offset)
                compressed = block_store.read(compressed_length)
                if len(compressed) != compressed_length:
                    raise SystemExit("compressed block store ended unexpectedly")
                current_block_payload = zlib.decompress(compressed)
                if len(current_block_payload) != uncompressed_length:
                    raise SystemExit("compressed block length metadata is inconsistent")
                current_block_id = block_id
                blocks_read += 1
                compressed_bytes_read += compressed_length

            end = block_offset + record_length
            chunk = current_block_payload[block_offset:end]
            if len(chunk) != record_length:
                raise SystemExit("record location exceeds its compressed block")
            payload.extend(chunk)

    records = tuple(iter_chain_records(io.StringIO(payload.decode("ascii"))))
    if len(records) != len(rows):
        raise SystemExit("compact record count changed during parsing")
    return records, len(rows), blocks_read, compressed_bytes_read


def _project_records(
    records: tuple[ChainRecord, ...],
    locus: _RequestedLocus,
    *,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> CandidateTuple:
    return tuple(
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


def _sqlite_candidates(
    paths: _IndexPaths,
    locus: _RequestedLocus,
    *,
    bin_width: int,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, int, float]:
    started = time.perf_counter()
    records, selected_records = _sqlite_records(
        paths,
        target_start=locus.interval.start,
        target_end=locus.interval.end,
        bin_width=bin_width,
    )
    candidates = _project_records(
        records,
        locus,
        target_assembly=target_assembly,
        provenance=provenance,
    )
    return candidates, selected_records, time.perf_counter() - started


def _compact_candidates(
    paths: _CompactPaths,
    locus: _RequestedLocus,
    *,
    bin_width: int,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, int, int, int, float]:
    started = time.perf_counter()
    records, selected_records, blocks_read, compressed_bytes_read = _compact_records(
        paths,
        target_start=locus.interval.start,
        target_end=locus.interval.end,
        bin_width=bin_width,
    )
    candidates = _project_records(
        records,
        locus,
        target_assembly=target_assembly,
        provenance=provenance,
    )
    return (
        candidates,
        selected_records,
        blocks_read,
        compressed_bytes_read,
        time.perf_counter() - started,
    )


def _write_manifest(
    output_dir: Path,
    *,
    source_manifest: dict[str, object],
    block_size: int,
    sequence_stats: dict[str, dict[str, object]],
) -> None:
    manifest = {
        "format": _COMPACT_FORMAT,
        "source_chain_sha256": source_manifest.get("source_chain_sha256"),
        "source_chain_size_bytes": source_manifest.get("source_chain_size_bytes"),
        "bin_width": source_manifest.get("bin_width"),
        "block_size": block_size,
        "sequences": sequence_stats,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    args = _parse_args()
    if args.block_size <= 0:
        raise SystemExit("block size must be positive")

    source_manifest = _load_index_manifest(args.index_dir)
    raw_bin_width = source_manifest.get("bin_width")
    if not isinstance(raw_bin_width, int) or raw_bin_width <= 0:
        raise SystemExit("bin-index manifest has invalid bin width")
    bin_width = raw_bin_width

    source_assembly = AssemblyIdentifier(args.source_db, "UCSC")
    target_assembly = AssemblyIdentifier(args.target_db, "UCSC")
    loci = tuple(_parse_locus(text, source_assembly) for text in args.locus)
    sequence_names = tuple(sorted({locus.interval.sequence_name for locus in loci}))
    provenance = ProvenanceSource("benchmark-chain", "benchmark chain resource")

    if args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True)

    print(f"source chain SHA-256: {source_manifest.get('source_chain_sha256')}")
    print(f"bin width: {bin_width} bases")
    print(f"compression block target: {args.block_size} bytes")
    print(f"targeted source sequences: {', '.join(sequence_names)}")
    print("repacking each chain record once into independently compressed blocks...")

    compact_paths: dict[str, _CompactPaths] = {}
    source_paths_by_sequence: dict[str, _IndexPaths] = {}
    sequence_stats: dict[str, dict[str, object]] = {}
    total_started = time.perf_counter()
    for sequence_name in sequence_names:
        source_paths = _index_paths(args.index_dir, sequence_name)
        _validate_index_paths(source_paths, sequence_name)
        source_paths_by_sequence[sequence_name] = source_paths
        paths, record_count, block_count, elapsed = _build_compact_sequence_store(
            source_paths=source_paths,
            output_dir=args.output_dir,
            sequence_name=sequence_name,
            block_size=args.block_size,
        )
        compact_paths[sequence_name] = paths
        compact_size = paths.database.stat().st_size + paths.blocks.stat().st_size
        source_size = (
            source_paths.database.stat().st_size + source_paths.records.stat().st_size
        )
        sequence_stats[sequence_name] = {
            "database": paths.database.name,
            "blocks": paths.blocks.name,
            "record_count": record_count,
            "block_count": block_count,
            "database_size_bytes": paths.database.stat().st_size,
            "block_store_size_bytes": paths.blocks.stat().st_size,
            "build_seconds": elapsed,
        }
        print(
            f"compact {sequence_name}: records={record_count} blocks={block_count} "
            f"build={elapsed:.3f}s compact={compact_size / (1024**2):.2f} MiB "
            f"sqlite+records={source_size / (1024**2):.2f} MiB "
            f"ratio={compact_size / source_size:.3f}x"
        )
    print(f"total compact build: {time.perf_counter() - total_started:.3f}s")

    _write_manifest(
        args.output_dir,
        source_manifest=source_manifest,
        block_size=args.block_size,
        sequence_stats=sequence_stats,
    )

    for locus in loci:
        source_paths = source_paths_by_sequence[locus.interval.sequence_name]
        baseline, baseline_selected, baseline_elapsed = _sqlite_candidates(
            source_paths,
            locus,
            bin_width=bin_width,
            target_assembly=target_assembly,
            provenance=provenance,
        )
        (
            compact,
            selected,
            blocks_read,
            bytes_read,
            compact_elapsed,
        ) = _compact_candidates(
            compact_paths[locus.interval.sequence_name],
            locus,
            bin_width=bin_width,
            target_assembly=target_assembly,
            provenance=provenance,
        )

        print()
        print(f"source locus: {args.source_db} {locus.text}")
        print(
            f"sqlite baseline: candidates={len(baseline)} "
            f"elapsed={baseline_elapsed:.3f}s "
            f"selected_records={baseline_selected}"
        )
        print(
            f"blocked query: candidates={len(compact)} elapsed={compact_elapsed:.3f}s "
            f"selected_records={selected} blocks_read={blocks_read} "
            f"compressed_read={bytes_read / 1024:.1f} KiB"
        )
        if compact != baseline:
            raise SystemExit(
                f"candidate mismatch for {locus.text}: "
                f"baseline={len(baseline)} compact={len(compact)}"
            )
        print(f"equivalence: PASS ({len(compact)} candidates)")
        if compact_elapsed > 0:
            print(
                f"relative speed vs sqlite: {baseline_elapsed / compact_elapsed:.2f}x"
            )


if __name__ == "__main__":
    main()
