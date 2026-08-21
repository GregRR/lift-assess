#!/usr/bin/env python3
"""Prototype a 64-KiB exact interval index over Milestone-18 sequence shards.

This benchmark reuses sequence shards produced by ``benchmark_chain_sequence_shards.py``.
It does not rescan the original multi-gigabyte provider chain. For each selected source
sequence, it writes every chain record exactly once to an uncompressed local record
store and records exact 0-based target intervals in a SQLite bin-membership index.
Chains spanning bin boundaries receive multiple lightweight membership rows, not
multiple payload copies.

Queries use the bins only to find a small superset of potentially overlapping records,
then apply the exact half-open interval predicate before parsing records with the
existing chain parser and projecting them with the production projection function.
The script compares those candidates against a scan of the corresponding sequence shard
and requires exact tuple equality, including encounter order.
"""

from __future__ import annotations

import argparse
import io
import json
import shutil
import sqlite3
import time
from dataclasses import dataclass
from pathlib import Path

from _selective_chain_traversal import iter_chain_file_overlapping_target_interval

from liftassess.chain import ChainRecord, iter_chain_records
from liftassess.models import (
    AssemblyIdentifier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.projection import project_interval_through_chain
from liftassess.resource_files import iter_chain_file

CandidateTuple = tuple[NormalizedCandidate, ...]
_SEQUENCE_SHARD_FORMAT = "liftassess-chain-sequence-shards-v1"
_INDEX_FORMAT = "liftassess-chain-bin-index-v1"
_DEFAULT_BIN_WIDTH = 65_536
_INSERT_BATCH_SIZE = 10_000


@dataclass(frozen=True)
class _RequestedLocus:
    text: str
    interval: GenomicInterval


@dataclass(frozen=True)
class _IndexPaths:
    database: Path
    records: Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bin-width", type=int, default=_DEFAULT_BIN_WIDTH)
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


def _sequence_key(sequence_name: str) -> str:
    return sequence_name.encode("utf-8").hex()


def _index_paths(output_dir: Path, sequence_name: str) -> _IndexPaths:
    key = _sequence_key(sequence_name)
    return _IndexPaths(
        database=output_dir / f"seq-{key}.sqlite3",
        records=output_dir / f"seq-{key}.records",
    )


def _load_shard_manifest(shard_dir: Path) -> dict[str, object]:
    path = shard_dir / "manifest.json"
    if not path.is_file():
        raise SystemExit(f"sequence-shard manifest is missing: {path}")
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("format") != _SEQUENCE_SHARD_FORMAT:
        raise SystemExit("sequence-shard manifest has an unsupported format")
    return raw


def _shard_path(
    shard_dir: Path,
    manifest: dict[str, object],
    sequence_name: str,
) -> Path:
    raw_shards = manifest.get("shards")
    if not isinstance(raw_shards, dict):
        raise SystemExit("sequence-shard manifest has invalid shard metadata")
    item = raw_shards.get(sequence_name)
    if not isinstance(item, dict) or not isinstance(item.get("path"), str):
        raise SystemExit(f"sequence-shard manifest lacks sequence {sequence_name}")
    path = shard_dir / item["path"]
    if not path.is_file():
        raise SystemExit(f"sequence shard is missing: {path}")
    return path


def _configure_database(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("PRAGMA temp_store = MEMORY")
    connection.executescript(
        """
        CREATE TABLE chains (
            record_id INTEGER PRIMARY KEY,
            target_start INTEGER NOT NULL,
            target_end INTEGER NOT NULL,
            record_offset INTEGER NOT NULL,
            record_length INTEGER NOT NULL
        );
        CREATE TABLE bin_memberships (
            bin_id INTEGER NOT NULL,
            record_id INTEGER NOT NULL,
            PRIMARY KEY (bin_id, record_id)
        ) WITHOUT ROWID;
        """
    )


def _flush_rows(
    connection: sqlite3.Connection,
    chain_rows: list[tuple[int, int, int, int, int]],
    membership_rows: list[tuple[int, int]],
) -> None:
    connection.executemany(
        "INSERT INTO chains VALUES (?, ?, ?, ?, ?)",
        chain_rows,
    )
    connection.executemany(
        "INSERT INTO bin_memberships VALUES (?, ?)",
        membership_rows,
    )
    chain_rows.clear()
    membership_rows.clear()


def _build_sequence_index(
    *,
    shard_path: Path,
    output_dir: Path,
    sequence_name: str,
    bin_width: int,
) -> tuple[_IndexPaths, int, int, float]:
    paths = _index_paths(output_dir, sequence_name)
    paths.database.unlink(missing_ok=True)
    paths.records.unlink(missing_ok=True)

    chain_rows: list[tuple[int, int, int, int, int]] = []
    membership_rows: list[tuple[int, int]] = []
    record_count = 0
    membership_count = 0
    started = time.perf_counter()

    with (
        sqlite3.connect(paths.database) as connection,
        paths.records.open("wb") as record_store,
    ):
        _configure_database(connection)
        for record_id, chain in enumerate(iter_chain_file(shard_path)):
            if chain.target_name != sequence_name:
                raise SystemExit(
                    f"shard {shard_path} contains unexpected target sequence "
                    f"{chain.target_name}"
                )
            payload = _format_chain_record(chain)
            offset = record_store.tell()
            record_store.write(payload)
            chain_rows.append(
                (
                    record_id,
                    chain.target_start,
                    chain.target_end,
                    offset,
                    len(payload),
                )
            )

            first_bin = chain.target_start // bin_width
            last_bin = (chain.target_end - 1) // bin_width
            for bin_id in range(first_bin, last_bin + 1):
                membership_rows.append((bin_id, record_id))
                membership_count += 1
            record_count += 1

            if len(chain_rows) >= _INSERT_BATCH_SIZE:
                _flush_rows(connection, chain_rows, membership_rows)

        if chain_rows:
            _flush_rows(connection, chain_rows, membership_rows)
        connection.commit()
        connection.execute("ANALYZE")
        connection.commit()

    return paths, record_count, membership_count, time.perf_counter() - started


def _baseline_candidates(
    shard_path: Path,
    locus: _RequestedLocus,
    *,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, float]:
    started = time.perf_counter()
    candidates = tuple(
        candidate
        for chain in iter_chain_file_overlapping_target_interval(
            shard_path,
            target_name=locus.interval.sequence_name,
            target_start=locus.interval.start,
            target_end=locus.interval.end,
        )
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
    return candidates, time.perf_counter() - started


def _indexed_records(
    paths: _IndexPaths,
    *,
    target_start: int,
    target_end: int,
    bin_width: int,
) -> tuple[tuple[ChainRecord, ...], int, int]:
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
        raise SystemExit("indexed record count changed during parsing")
    return records, len(rows), last_bin - first_bin + 1


def _indexed_candidates(
    paths: _IndexPaths,
    locus: _RequestedLocus,
    *,
    bin_width: int,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, int, int, float]:
    started = time.perf_counter()
    records, selected_records, bins_touched = _indexed_records(
        paths,
        target_start=locus.interval.start,
        target_end=locus.interval.end,
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
        time.perf_counter() - started,
    )


def _write_manifest(
    output_dir: Path,
    *,
    shard_manifest: dict[str, object],
    bin_width: int,
    sequence_stats: dict[str, dict[str, object]],
) -> None:
    manifest = {
        "format": _INDEX_FORMAT,
        "bin_width": bin_width,
        "source_chain_sha256": shard_manifest.get("source_chain_sha256"),
        "source_chain_size_bytes": shard_manifest.get("source_chain_size_bytes"),
        "sequences": sequence_stats,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    args = _parse_args()
    if args.bin_width <= 0:
        raise SystemExit("bin width must be positive")

    shard_manifest = _load_shard_manifest(args.shard_dir)
    source_assembly = AssemblyIdentifier(args.source_db, "UCSC")
    target_assembly = AssemblyIdentifier(args.target_db, "UCSC")
    loci = tuple(_parse_locus(text, source_assembly) for text in args.locus)
    sequence_names = tuple(sorted({locus.interval.sequence_name for locus in loci}))
    provenance = ProvenanceSource("benchmark-chain", "benchmark chain resource")

    if args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True)

    print(f"source chain SHA-256: {shard_manifest.get('source_chain_sha256')}")
    print(f"bin width: {args.bin_width} bases")
    print(f"targeted source sequences: {', '.join(sequence_names)}")
    print("building exact bin-membership indexes from existing sequence shards...")

    index_paths: dict[str, _IndexPaths] = {}
    sequence_stats: dict[str, dict[str, object]] = {}
    total_started = time.perf_counter()
    for sequence_name in sequence_names:
        shard_path = _shard_path(args.shard_dir, shard_manifest, sequence_name)
        paths, record_count, membership_count, elapsed = _build_sequence_index(
            shard_path=shard_path,
            output_dir=args.output_dir,
            sequence_name=sequence_name,
            bin_width=args.bin_width,
        )
        index_paths[sequence_name] = paths
        sequence_stats[sequence_name] = {
            "database": paths.database.name,
            "records": paths.records.name,
            "record_count": record_count,
            "membership_count": membership_count,
            "database_size_bytes": paths.database.stat().st_size,
            "record_store_size_bytes": paths.records.stat().st_size,
            "build_seconds": elapsed,
        }
        print(
            f"index {sequence_name}: records={record_count} "
            f"memberships={membership_count} build={elapsed:.3f}s "
            f"db={paths.database.stat().st_size / (1024**2):.2f} MiB "
            f"records_file={paths.records.stat().st_size / (1024**2):.2f} MiB"
        )
    print(f"total index build: {time.perf_counter() - total_started:.3f}s")

    _write_manifest(
        args.output_dir,
        shard_manifest=shard_manifest,
        bin_width=args.bin_width,
        sequence_stats=sequence_stats,
    )

    for locus in loci:
        shard_path = _shard_path(
            args.shard_dir, shard_manifest, locus.interval.sequence_name
        )
        baseline, baseline_elapsed = _baseline_candidates(
            shard_path,
            locus,
            target_assembly=target_assembly,
            provenance=provenance,
        )
        indexed, selected_records, bins_touched, indexed_elapsed = _indexed_candidates(
            index_paths[locus.interval.sequence_name],
            locus,
            bin_width=args.bin_width,
            target_assembly=target_assembly,
            provenance=provenance,
        )

        print()
        print(f"source locus: {args.source_db} {locus.text}")
        print(
            f"sequence-shard baseline: candidates={len(baseline)} "
            f"elapsed={baseline_elapsed:.3f}s"
        )
        print(
            f"bin-index query: candidates={len(indexed)} elapsed={indexed_elapsed:.3f}s "
            f"bins={bins_touched} selected_records={selected_records}"
        )
        if indexed != baseline:
            raise SystemExit(
                f"candidate mismatch for {locus.text}: "
                f"baseline={len(baseline)} indexed={len(indexed)}"
            )
        print(f"equivalence: PASS ({len(indexed)} candidates)")
        if indexed_elapsed > 0:
            print(
                f"speedup vs sequence shard: {baseline_elapsed / indexed_elapsed:.2f}x"
            )


if __name__ == "__main__":
    main()
