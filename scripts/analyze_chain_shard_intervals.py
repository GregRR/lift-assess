#!/usr/bin/env python3
"""Measure chain-span/binning behavior in existing Milestone-18 sequence shards.

This script is intentionally analysis-only. It reuses sequence shards already built by
``benchmark_chain_sequence_shards.py`` so that candidate interval-index designs can be
compared without another multi-gigabyte source-chain traversal.

For each sequence shard and candidate bin width, it reports how many bin memberships
would be created by a simple overlap-bin representation and the largest bin occupancy.
No index files are written.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from liftassess.resource_files import iter_chain_file

_SEQUENCE_SHARD_FORMAT = "liftassess-chain-sequence-shards-v1"
_DEFAULT_BIN_WIDTHS = (65_536, 262_144, 1_048_576, 4_194_304)


@dataclass
class _BinStats:
    width: int
    memberships: int = 0
    long_records: int = 0
    bin_counts: dict[int, int] | None = None

    def __post_init__(self) -> None:
        self.bin_counts = defaultdict(int)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-dir", type=Path, required=True)
    parser.add_argument(
        "--bin-width",
        type=int,
        action="append",
        dest="bin_widths",
        help="bin width in source bases; repeat to compare widths",
    )
    parser.add_argument(
        "sequence",
        nargs="*",
        help="optional sequence names; default is every shard in the manifest",
    )
    return parser.parse_args()


def _load_manifest(shard_dir: Path) -> dict[str, object]:
    path = shard_dir / "manifest.json"
    if not path.is_file():
        raise SystemExit(f"sequence-shard manifest is missing: {path}")
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("format") != _SEQUENCE_SHARD_FORMAT:
        raise SystemExit("sequence-shard manifest has an unsupported format")
    return raw


def _selected_shards(
    shard_dir: Path,
    manifest: dict[str, object],
    requested: tuple[str, ...],
) -> tuple[tuple[str, Path], ...]:
    raw_shards = manifest.get("shards")
    if not isinstance(raw_shards, dict):
        raise SystemExit("sequence-shard manifest has invalid shard metadata")

    names = requested or tuple(sorted(str(name) for name in raw_shards))
    result: list[tuple[str, Path]] = []
    for name in names:
        item = raw_shards.get(name)
        if not isinstance(item, dict) or not isinstance(item.get("path"), str):
            raise SystemExit(f"sequence-shard manifest lacks sequence {name}")
        path = shard_dir / item["path"]
        if not path.is_file():
            raise SystemExit(f"sequence shard is missing: {path}")
        result.append((name, path))
    return tuple(result)


def _format_ratio(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "n/a"
    return f"{numerator / denominator:.3f}x"


def _analyze_sequence(
    sequence_name: str,
    path: Path,
    widths: tuple[int, ...],
) -> None:
    stats = {width: _BinStats(width) for width in widths}
    record_count = 0
    total_span = 0
    maximum_span = 0
    target_size: int | None = None

    for chain in iter_chain_file(path):
        if chain.target_name != sequence_name:
            raise SystemExit(
                f"shard {path} contains unexpected target sequence {chain.target_name}"
            )
        if target_size is None:
            target_size = chain.target_size
        elif target_size != chain.target_size:
            raise SystemExit(
                f"shard {path} contains inconsistent target sizes for {sequence_name}"
            )

        span = chain.target_end - chain.target_start
        record_count += 1
        total_span += span
        maximum_span = max(maximum_span, span)

        for width, item in stats.items():
            first_bin = chain.target_start // width
            last_bin = (chain.target_end - 1) // width
            memberships = last_bin - first_bin + 1
            item.memberships += memberships
            if memberships > 1:
                item.long_records += 1
            assert item.bin_counts is not None
            for bin_id in range(first_bin, last_bin + 1):
                item.bin_counts[bin_id] += 1

    print()
    print(f"sequence: {sequence_name}")
    print(f"shard: {path}")
    print(f"records: {record_count}")
    if target_size is not None:
        print(f"target size: {target_size} bases")
    if record_count:
        print(f"mean chain target span: {total_span / record_count:.1f} bases")
        print(f"maximum chain target span: {maximum_span} bases")

    for width in widths:
        item = stats[width]
        assert item.bin_counts is not None
        maximum_occupancy = max(item.bin_counts.values(), default=0)
        populated_bins = len(item.bin_counts)
        print(
            f"bin {width}: memberships={item.memberships} "
            f"duplication={_format_ratio(item.memberships, record_count)} "
            f"cross_bin_records={item.long_records} "
            f"populated_bins={populated_bins} "
            f"max_bin_records={maximum_occupancy}"
        )


def main() -> None:
    args = _parse_args()
    widths = tuple(sorted(set(args.bin_widths or _DEFAULT_BIN_WIDTHS)))
    if any(width <= 0 for width in widths):
        raise SystemExit("bin widths must be positive")

    manifest = _load_manifest(args.shard_dir)
    print(f"source chain SHA-256: {manifest.get('source_chain_sha256')}")
    print("candidate bin widths: " + ", ".join(str(width) for width in widths))

    for sequence_name, path in _selected_shards(
        args.shard_dir, manifest, tuple(args.sequence)
    ):
        _analyze_sequence(sequence_name, path, widths)


if __name__ == "__main__":
    main()
