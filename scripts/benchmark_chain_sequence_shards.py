#!/usr/bin/env python3
"""Prototype sequence-sharded chain access for Milestone 18.

This benchmark builds temporary gzip shards for only the source sequence names used by
one or more requested loci. The build performs exactly one full parse of the original
chain resource, preserving the existing parser's validation semantics. Each selected
``ChainRecord`` is serialized into a derived per-sequence shard; the original provider
resource remains the scientific/provenance source.

During that same full pass, the script derives baseline candidates for every requested
locus with the production ``project_interval_through_chain`` function. It then queries
the derived shards through the selective-materialization parser and requires exact
candidate-tuple equality. This tests whether a persistent sequence-addressable local
representation is a worthwhile next step without yet changing production assessment.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import time
from collections import defaultdict
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from _selective_chain_traversal import iter_chain_file_overlapping_target_interval

from liftassess.chain import ChainRecord
from liftassess.models import (
    AssemblyIdentifier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.projection import project_interval_through_chain
from liftassess.resource_cache import CachedResource, load_cached_ucsc_resource_bundle
from liftassess.resource_files import iter_chain_file

CandidateTuple = tuple[NormalizedCandidate, ...]
_INDEX_FORMAT = "liftassess-chain-sequence-shards-v1"


@dataclass(frozen=True)
class _RequestedLocus:
    text: str
    interval: GenomicInterval


@dataclass
class _ShardWriter:
    path: Path
    handle: TextIO
    record_count: int = 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("source_db")
    parser.add_argument("target_db")
    parser.add_argument(
        "locus",
        nargs="+",
        help="one or more 1-based inclusive source loci, chr:start-end",
    )
    parser.add_argument(
        "--reuse",
        action="store_true",
        help="reuse an existing matching shard set instead of rebuilding it",
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


def _format_chain_record(record: ChainRecord) -> str:
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
    return "\n".join(lines) + "\n\n"


def _shard_filename(sequence_name: str) -> str:
    # Sequence names can contain punctuation. Hex encoding keeps the path reversible
    # without assuming a naming convention or risking path separators.
    return f"seq-{sequence_name.encode('utf-8').hex()}.chain.gz"


def _manifest_path(output_dir: Path) -> Path:
    return output_dir / "manifest.json"


def _write_manifest(
    output_dir: Path,
    *,
    chain_resource: CachedResource,
    source_db: str,
    target_db: str,
    writers: dict[str, _ShardWriter],
    build_seconds: float,
) -> None:
    manifest = {
        "format": _INDEX_FORMAT,
        "source_db": source_db,
        "target_db": target_db,
        "source_chain_sha256": chain_resource.sha256,
        "source_chain_size_bytes": chain_resource.size_bytes,
        "source_chain_url": chain_resource.source_url,
        "build_seconds": build_seconds,
        "shards": {
            sequence_name: {
                "path": writer.path.name,
                "record_count": writer.record_count,
                "size_bytes": writer.path.stat().st_size,
            }
            for sequence_name, writer in sorted(writers.items())
        },
    }
    _manifest_path(output_dir).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _load_reusable_shards(
    output_dir: Path,
    *,
    chain_resource: CachedResource,
    source_db: str,
    target_db: str,
    sequence_names: tuple[str, ...],
) -> dict[str, Path]:
    manifest_path = _manifest_path(output_dir)
    if not manifest_path.is_file():
        raise SystemExit(f"--reuse requested but manifest is missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = (
        manifest.get("format") == _INDEX_FORMAT
        and manifest.get("source_db") == source_db
        and manifest.get("target_db") == target_db
        and manifest.get("source_chain_sha256") == chain_resource.sha256
        and manifest.get("source_chain_size_bytes") == chain_resource.size_bytes
    )
    if not expected:
        raise SystemExit(
            "existing shard manifest does not match the cached chain resource"
        )

    raw_shards = manifest.get("shards")
    if not isinstance(raw_shards, dict):
        raise SystemExit("existing shard manifest has invalid shard metadata")
    paths: dict[str, Path] = {}
    for sequence_name in sequence_names:
        item = raw_shards.get(sequence_name)
        if not isinstance(item, dict) or not isinstance(item.get("path"), str):
            raise SystemExit(f"existing shard set lacks sequence {sequence_name}")
        path = output_dir / item["path"]
        if not path.is_file():
            raise SystemExit(f"existing shard file is missing: {path}")
        paths[sequence_name] = path
    return paths


def _build_targeted_shards(
    *,
    chain_resource: CachedResource,
    output_dir: Path,
    source_db: str,
    target_db: str,
    loci: tuple[_RequestedLocus, ...],
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[dict[str, Path], dict[str, CandidateTuple], float]:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    sequence_names = tuple(sorted({locus.interval.sequence_name for locus in loci}))
    baseline_lists: dict[str, list[NormalizedCandidate]] = defaultdict(list)
    loci_by_sequence: dict[str, list[_RequestedLocus]] = defaultdict(list)
    for locus in loci:
        loci_by_sequence[locus.interval.sequence_name].append(locus)

    started = time.perf_counter()
    with ExitStack() as stack:
        writers: dict[str, _ShardWriter] = {}
        for sequence_name in sequence_names:
            path = output_dir / _shard_filename(sequence_name)
            handle = stack.enter_context(
                gzip.open(path, "wt", encoding="ascii", newline="\n")
            )
            writers[sequence_name] = _ShardWriter(path=path, handle=handle)

        for chain in iter_chain_file(chain_resource.path):
            writer = writers.get(chain.target_name)
            if writer is None:
                continue
            writer.handle.write(_format_chain_record(chain))
            writer.record_count += 1
            for locus in loci_by_sequence[chain.target_name]:
                candidate = project_interval_through_chain(
                    locus.interval,
                    chain,
                    target_assembly=target_assembly,
                    mapping_provenance=provenance,
                )
                if candidate is not None:
                    baseline_lists[locus.text].append(candidate)
    build_seconds = time.perf_counter() - started

    _write_manifest(
        output_dir,
        chain_resource=chain_resource,
        source_db=source_db,
        target_db=target_db,
        writers=writers,
        build_seconds=build_seconds,
    )
    paths = {sequence_name: writer.path for sequence_name, writer in writers.items()}
    baselines = {locus.text: tuple(baseline_lists[locus.text]) for locus in loci}
    return paths, baselines, build_seconds


def _query_shard(
    path: Path,
    locus: _RequestedLocus,
    *,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> tuple[CandidateTuple, float]:
    started = time.perf_counter()
    candidates = tuple(
        candidate
        for chain in iter_chain_file_overlapping_target_interval(
            path,
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


def main() -> None:
    args = _parse_args()
    bundle = load_cached_ucsc_resource_bundle(
        args.cache_root,
        args.source_db,
        args.target_db,
    )
    if bundle is None:
        raise SystemExit(
            "no complete verified cached bundle found for "
            f"{args.source_db}→{args.target_db} under {args.cache_root}"
        )

    source_assembly = AssemblyIdentifier(args.source_db, "UCSC")
    target_assembly = AssemblyIdentifier(args.target_db, "UCSC")
    loci = tuple(_parse_locus(text, source_assembly) for text in args.locus)
    sequence_names = tuple(sorted({locus.interval.sequence_name for locus in loci}))
    provenance = ProvenanceSource("benchmark-chain", "benchmark chain resource")

    print(f"chain: {bundle.chain.path}")
    print(f"chain SHA-256: {bundle.chain.sha256}")
    print(f"targeted source sequences: {', '.join(sequence_names)}")

    baselines: dict[str, CandidateTuple] | None
    if args.reuse:
        shard_paths = _load_reusable_shards(
            args.output_dir,
            chain_resource=bundle.chain,
            source_db=args.source_db,
            target_db=args.target_db,
            sequence_names=sequence_names,
        )
        baselines = None
        print("shard build: REUSED existing matching shard set")
    else:
        print("building targeted sequence shards with one full validated chain pass...")
        shard_paths, baselines, build_seconds = _build_targeted_shards(
            chain_resource=bundle.chain,
            output_dir=args.output_dir,
            source_db=args.source_db,
            target_db=args.target_db,
            loci=loci,
            target_assembly=target_assembly,
            provenance=provenance,
        )
        print(f"shard build: {build_seconds:.3f}s")

    manifest = json.loads(_manifest_path(args.output_dir).read_text(encoding="utf-8"))
    for sequence_name in sequence_names:
        item = manifest["shards"][sequence_name]
        print(
            f"shard {sequence_name}: records={item['record_count']} "
            f"size={item['size_bytes'] / (1024**2):.2f} MiB"
        )

    for locus in loci:
        candidates, elapsed = _query_shard(
            shard_paths[locus.interval.sequence_name],
            locus,
            target_assembly=target_assembly,
            provenance=provenance,
        )
        print()
        print(f"source locus: {args.source_db} {locus.text}")
        print(f"shard query: candidates={len(candidates)} elapsed={elapsed:.3f}s")
        if baselines is None:
            print(
                "equivalence: NOT CHECKED (--reuse does not repeat the full baseline pass)"
            )
            continue
        expected = baselines[locus.text]
        if candidates != expected:
            raise SystemExit(
                f"candidate mismatch for {locus.text}: "
                f"baseline={len(expected)} shard={len(candidates)}"
            )
        print(f"equivalence: PASS ({len(candidates)} candidates)")


if __name__ == "__main__":
    main()
