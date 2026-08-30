#!/usr/bin/env python3
"""Benchmark full versus selective chain traversal for one source interval.

This Milestone-18 prototype isolates the chain candidate-generation hot path. It
compares the existing full parser with query-aware selective materialization while
feeding both record streams through the same projection implementation. The two
candidate tuples must be exactly equal before the speedup is reported.

The selective parser preserves the current parser's strict header, block, and span
validation for skipped records. The performance hypothesis is therefore specifically
whether avoiding irrelevant ``ChainBlock``/``ChainRecord`` materialization is enough to
justify this traversal boundary.
"""

from __future__ import annotations

import argparse
import time
from collections.abc import Callable, Iterable
from pathlib import Path

from _selective_chain_traversal import iter_chain_file_overlapping_target_interval

from liftassess.chain import ChainRecord
from liftassess.models import (
    AssemblyIdentifier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.projection import iter_candidates_from_chains
from liftassess.resource_cache import load_cached_ucsc_resource_bundle
from liftassess.resource_files import iter_chain_file

CandidateTuple = tuple[NormalizedCandidate, ...]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--chain", type=Path, help="local chain or chain.gz resource")
    source.add_argument(
        "--cache-root",
        type=Path,
        help="liftAssess cache root; the verified cached chain is selected automatically",
    )
    parser.add_argument("source_db", help="source assembly/database label")
    parser.add_argument("target_db", help="target assembly/database label")
    parser.add_argument(
        "locus",
        nargs="+",
        help="one or more 1-based inclusive source loci, chr:start-end",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="number of timed full/selective traversals (default: 1)",
    )
    return parser.parse_args()


def _parse_locus(text: str, assembly: AssemblyIdentifier) -> GenomicInterval:
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
    return GenomicInterval(
        assembly=assembly,
        sequence_name=sequence_name,
        start=start - 1,
        end=end,
    )


def _candidate_tuple(
    source_interval: GenomicInterval,
    chains: Iterable[ChainRecord],
    *,
    target_assembly: AssemblyIdentifier,
    provenance: ProvenanceSource,
) -> CandidateTuple:
    return tuple(
        iter_candidates_from_chains(
            source_interval,
            chains,
            target_assembly=target_assembly,
            mapping_provenance=provenance,
        )
    )


def _time_call(
    label: str,
    call: Callable[[], CandidateTuple],
    *,
    repeat: int,
) -> tuple[CandidateTuple, list[float]]:
    samples: list[float] = []
    first_result: CandidateTuple | None = None
    for _ in range(repeat):
        started = time.perf_counter()
        result = call()
        samples.append(time.perf_counter() - started)
        if first_result is None:
            first_result = result
        elif result != first_result:
            raise SystemExit(
                f"{label} traversal was not deterministic across repetitions"
            )
    assert first_result is not None
    print(
        f"{label}: candidates={len(first_result)} "
        f"min={min(samples):.3f}s mean={sum(samples) / len(samples):.3f}s "
        f"max={max(samples):.3f}s"
    )
    return first_result, samples


def _resolve_chain_path(args: argparse.Namespace) -> Path:
    chain_path = args.chain
    if chain_path is not None:
        if not isinstance(chain_path, Path):
            raise SystemExit("--chain must resolve to a filesystem path")
        return chain_path
    assert args.cache_root is not None
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
    return bundle.chain.path


def _benchmark_locus(
    *,
    chain_path: Path,
    source_db: str,
    target_db: str,
    locus: str,
    repeat: int,
) -> None:
    source_assembly = AssemblyIdentifier(source_db, "UCSC")
    target_assembly = AssemblyIdentifier(target_db, "UCSC")
    source_interval = _parse_locus(locus, source_assembly)
    provenance = ProvenanceSource("benchmark-chain", "benchmark chain resource")

    def full_call() -> CandidateTuple:
        return _candidate_tuple(
            source_interval,
            iter_chain_file(chain_path),
            target_assembly=target_assembly,
            provenance=provenance,
        )

    def selective_call() -> CandidateTuple:
        return _candidate_tuple(
            source_interval,
            iter_chain_file_overlapping_target_interval(
                chain_path,
                target_name=source_interval.sequence_name,
                target_start=source_interval.start,
                target_end=source_interval.end,
            ),
            target_assembly=target_assembly,
            provenance=provenance,
        )

    print(f"source locus: {source_db} {locus}")
    print(f"target assembly: {target_db}")
    print("timing both traversals and verifying exact candidate equivalence...")
    full, full_samples = _time_call("full", full_call, repeat=repeat)
    selective, selective_samples = _time_call(
        "selective", selective_call, repeat=repeat
    )
    if full != selective:
        raise SystemExit("candidate mismatch between full and selective traversal")
    print(f"equivalence: PASS ({len(full)} candidates)")
    full_mean = sum(full_samples) / len(full_samples)
    selective_mean = sum(selective_samples) / len(selective_samples)
    print(f"speedup: {full_mean / selective_mean:.2f}x")


def main() -> None:
    args = _parse_args()
    if args.repeat <= 0:
        raise SystemExit("--repeat must be positive")

    chain_path = _resolve_chain_path(args)
    print(f"chain: {chain_path}")
    for index, locus in enumerate(args.locus):
        if index:
            print()
        _benchmark_locus(
            chain_path=chain_path,
            source_db=args.source_db,
            target_db=args.target_db,
            locus=locus,
            repeat=args.repeat,
        )


if __name__ == "__main__":
    main()
