#!/usr/bin/env python3
"""Prototype compact independently compressed 64-KiB chain-bin packs.

This benchmark consumes the SQLite + uncompressed record-store artifacts produced by
``benchmark_chain_bin_index.py``. It does not read the original multi-gigabyte provider
chain or the intermediate sequence shards while building the pack.

Each genomic bin becomes one independently gzip-compressed frame in a single pack file.
A chain that spans bin boundaries is referenced in each overlapping frame; the measured
64-KiB membership duplication is small. Every packed record carries its original
sequence-shard encounter ID plus exact target bounds so a query can de-duplicate records
from multiple bins, apply the exact half-open overlap predicate before chain parsing,
and restore reproducible encounter order.

The script compares packed-bin candidates against the already-validated SQLite bin index
and requires exact tuple equality.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import sqlite3
import struct
import time
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
_PACK_FORMAT = "liftassess-chain-bin-pack-v1"
_RECORD_HEADER = struct.Struct("<QQQI")


@dataclass(frozen=True)
class _RequestedLocus:
    text: str
    interval: GenomicInterval


@dataclass(frozen=True)
class _IndexPaths:
    database: Path
    records: Path


@dataclass(frozen=True)
class _PackPaths:
    pack: Path
    index: Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
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


def _load_index_manifest(index_dir: Path) -> dict[str, object]:
    path = index_dir / "manifest.json"
    if not path.is_file():
        raise SystemExit(f"bin-index manifest is missing: {path}")
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("format") != _INDEX_FORMAT:
        raise SystemExit("bin-index manifest has an unsupported format")
    return raw


def _index_paths(
    index_dir: Path,
    manifest: dict[str, object],
    sequence_name: str,
) -> _IndexPaths:
    raw_sequences = manifest.get("sequences")
    if not isinstance(raw_sequences, dict):
        raise SystemExit("bin-index manifest has invalid sequence metadata")
    item = raw_sequences.get(sequence_name)
    if not isinstance(item, dict):
        raise SystemExit(f"bin-index manifest lacks sequence {sequence_name}")
    database_name = item.get("database")
    records_name = item.get("records")
    if not isinstance(database_name, str) or not isinstance(records_name, str):
        raise SystemExit("bin-index manifest has invalid artifact paths")
    paths = _IndexPaths(
        database=index_dir / database_name,
        records=index_dir / records_name,
    )
    if not paths.database.is_file() or not paths.records.is_file():
        raise SystemExit(f"bin-index artifacts are missing for {sequence_name}")
    return paths


def _pack_paths(output_dir: Path, sequence_name: str) -> _PackPaths:
    key = _sequence_key(sequence_name)
    return _PackPaths(
        pack=output_dir / f"seq-{key}.binpack",
        index=output_dir / f"seq-{key}.json",
    )


def _read_record(record_store: io.BufferedReader, offset: int, length: int) -> bytes:
    record_store.seek(offset)
    payload = record_store.read(length)
    if len(payload) != length:
        raise SystemExit("indexed chain record store ended unexpectedly")
    return payload


def _build_pack(
    *,
    index_paths: _IndexPaths,
    pack_paths: _PackPaths,
) -> tuple[int, int, int, float]:
    started = time.perf_counter()
    frame_index: dict[str, dict[str, int]] = {}
    total_memberships = 0
    uncompressed_bytes = 0

    with (
        sqlite3.connect(index_paths.database) as connection,
        index_paths.records.open("rb") as record_store,
        pack_paths.pack.open("wb") as pack_file,
    ):
        bin_ids = tuple(
            row[0]
            for row in connection.execute(
                "SELECT DISTINCT bin_id FROM bin_memberships ORDER BY bin_id"
            )
        )
        for bin_id in bin_ids:
            rows = tuple(
                connection.execute(
                    """
                    SELECT c.record_id, c.target_start, c.target_end,
                           c.record_offset, c.record_length
                    FROM bin_memberships AS b
                    JOIN chains AS c ON c.record_id = b.record_id
                    WHERE b.bin_id = ?
                    ORDER BY c.record_id
                    """,
                    (bin_id,),
                )
            )
            payload = bytearray()
            for record_id, target_start, target_end, offset, length in rows:
                record = _read_record(record_store, offset, length)
                payload.extend(
                    _RECORD_HEADER.pack(
                        record_id,
                        target_start,
                        target_end,
                        len(record),
                    )
                )
                payload.extend(record)
            compressed = gzip.compress(bytes(payload), compresslevel=6, mtime=0)
            frame_offset = pack_file.tell()
            pack_file.write(compressed)
            frame_index[str(bin_id)] = {
                "offset": frame_offset,
                "length": len(compressed),
                "memberships": len(rows),
            }
            total_memberships += len(rows)
            uncompressed_bytes += len(payload)

    pack_paths.index.write_text(
        json.dumps(
            {
                "format": _PACK_FORMAT,
                "frames": frame_index,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return (
        len(frame_index),
        total_memberships,
        uncompressed_bytes,
        time.perf_counter() - started,
    )


def _record_from_payload(payload: bytes) -> ChainRecord:
    records = tuple(iter_chain_records(io.StringIO(payload.decode("ascii"))))
    if len(records) != 1:
        raise SystemExit("packed record payload did not contain exactly one chain")
    return records[0]


def _sqlite_candidates(
    index_paths: _IndexPaths,
    locus: _RequestedLocus,
    *,
    bin_width: int,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, int, float]:
    started = time.perf_counter()
    first_bin = locus.interval.start // bin_width
    last_bin = (locus.interval.end - 1) // bin_width
    with sqlite3.connect(index_paths.database) as connection:
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
                (first_bin, last_bin, locus.interval.end, locus.interval.start),
            )
        )

    records: list[ChainRecord] = []
    with index_paths.records.open("rb") as record_store:
        for _record_id, offset, length in rows:
            records.append(
                _record_from_payload(_read_record(record_store, offset, length))
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
    return candidates, len(rows), time.perf_counter() - started


def _load_pack_index(path: Path) -> dict[int, tuple[int, int]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("format") != _PACK_FORMAT:
        raise SystemExit("bin-pack index has an unsupported format")
    frames = raw.get("frames")
    if not isinstance(frames, dict):
        raise SystemExit("bin-pack index has invalid frame metadata")
    result: dict[int, tuple[int, int]] = {}
    for bin_text, value in frames.items():
        if not isinstance(value, dict):
            raise SystemExit("bin-pack frame metadata is invalid")
        offset = value.get("offset")
        length = value.get("length")
        if not isinstance(offset, int) or not isinstance(length, int):
            raise SystemExit("bin-pack frame offset/length is invalid")
        result[int(bin_text)] = (offset, length)
    return result


def _packed_candidates(
    pack_paths: _PackPaths,
    locus: _RequestedLocus,
    *,
    bin_width: int,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, int, int, float]:
    started = time.perf_counter()
    frame_index = _load_pack_index(pack_paths.index)
    first_bin = locus.interval.start // bin_width
    last_bin = (locus.interval.end - 1) // bin_width

    records_by_id: dict[int, bytes] = {}
    memberships_read = 0
    with pack_paths.pack.open("rb") as pack_file:
        for bin_id in range(first_bin, last_bin + 1):
            frame = frame_index.get(bin_id)
            if frame is None:
                continue
            offset, length = frame
            pack_file.seek(offset)
            compressed = pack_file.read(length)
            if len(compressed) != length:
                raise SystemExit("bin-pack file ended unexpectedly")
            payload = memoryview(gzip.decompress(compressed))
            cursor = 0
            while cursor < len(payload):
                header_end = cursor + _RECORD_HEADER.size
                if header_end > len(payload):
                    raise SystemExit("truncated packed-record header")
                (
                    record_id,
                    target_start,
                    target_end,
                    record_length,
                ) = _RECORD_HEADER.unpack(payload[cursor:header_end])
                cursor = header_end
                record_end = cursor + record_length
                if record_end > len(payload):
                    raise SystemExit("truncated packed chain record")
                memberships_read += 1
                if (
                    target_start < locus.interval.end
                    and target_end > locus.interval.start
                    and record_id not in records_by_id
                ):
                    records_by_id[record_id] = bytes(payload[cursor:record_end])
                cursor = record_end

    records = tuple(
        _record_from_payload(records_by_id[record_id])
        for record_id in sorted(records_by_id)
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
        len(records),
        memberships_read,
        time.perf_counter() - started,
    )


def main() -> None:
    args = _parse_args()
    index_manifest = _load_index_manifest(args.index_dir)
    bin_width = index_manifest.get("bin_width")
    if not isinstance(bin_width, int) or bin_width <= 0:
        raise SystemExit("bin-index manifest has invalid bin width")

    source_assembly = AssemblyIdentifier(args.source_db, "UCSC")
    target_assembly = AssemblyIdentifier(args.target_db, "UCSC")
    loci = tuple(_parse_locus(text, source_assembly) for text in args.locus)
    sequence_names = tuple(sorted({locus.interval.sequence_name for locus in loci}))
    provenance = ProvenanceSource("benchmark-chain", "benchmark chain resource")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"source chain SHA-256: {index_manifest.get('source_chain_sha256')}")
    print(f"bin width: {bin_width} bases")
    print(f"targeted source sequences: {', '.join(sequence_names)}")
    print("packing independently compressed genomic bins from the existing index...")

    index_paths_by_sequence: dict[str, _IndexPaths] = {}
    pack_paths_by_sequence: dict[str, _PackPaths] = {}
    total_started = time.perf_counter()
    for sequence_name in sequence_names:
        index_paths = _index_paths(args.index_dir, index_manifest, sequence_name)
        pack_paths = _pack_paths(args.output_dir, sequence_name)
        pack_paths.pack.unlink(missing_ok=True)
        pack_paths.index.unlink(missing_ok=True)
        frames, memberships, uncompressed_bytes, elapsed = _build_pack(
            index_paths=index_paths,
            pack_paths=pack_paths,
        )
        index_paths_by_sequence[sequence_name] = index_paths
        pack_paths_by_sequence[sequence_name] = pack_paths
        original_bytes = (
            index_paths.database.stat().st_size + index_paths.records.stat().st_size
        )
        packed_bytes = pack_paths.pack.stat().st_size + pack_paths.index.stat().st_size
        print(
            f"pack {sequence_name}: frames={frames} memberships={memberships} "
            f"build={elapsed:.3f}s packed={packed_bytes / (1024**2):.2f} MiB "
            f"sqlite+records={original_bytes / (1024**2):.2f} MiB "
            f"payload_before_compression={uncompressed_bytes / (1024**2):.2f} MiB"
        )
    print(f"total pack build: {time.perf_counter() - total_started:.3f}s")

    for locus in loci:
        sequence_name = locus.interval.sequence_name
        baseline, baseline_records, baseline_elapsed = _sqlite_candidates(
            index_paths_by_sequence[sequence_name],
            locus,
            bin_width=bin_width,
            target_assembly=target_assembly,
            provenance=provenance,
        )
        packed, selected_records, memberships_read, packed_elapsed = _packed_candidates(
            pack_paths_by_sequence[sequence_name],
            locus,
            bin_width=bin_width,
            target_assembly=target_assembly,
            provenance=provenance,
        )

        print()
        print(f"source locus: {args.source_db} {locus.text}")
        print(
            f"sqlite baseline: candidates={len(baseline)} "
            f"elapsed={baseline_elapsed:.3f}s selected_records={baseline_records}"
        )
        print(
            f"compressed-bin query: candidates={len(packed)} "
            f"elapsed={packed_elapsed:.3f}s selected_records={selected_records} "
            f"memberships_read={memberships_read}"
        )
        if packed != baseline:
            raise SystemExit(
                f"candidate mismatch for {locus.text}: "
                f"baseline={len(baseline)} packed={len(packed)}"
            )
        print(f"equivalence: PASS ({len(packed)} candidates)")
        if packed_elapsed > 0:
            print(f"relative speed vs sqlite: {baseline_elapsed / packed_elapsed:.2f}x")


if __name__ == "__main__":
    main()
