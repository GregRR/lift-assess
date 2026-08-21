"""Command-line preparation of reusable liftAssess chain indexes.

The preparation command is deliberately cache-only. It never discovers or downloads
provider resources; it builds a derived acceleration artifact only from an already
verified cached UCSC bundle. Normal CLI assessment remains usable without an index and
falls back to the original full traversal when no matching usable index exists. Lower-level
library callers retain explicit control over query-time index-corruption recovery.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import TextIO

from .chain_index import (
    ChainIndexCorruptionError,
    ChainIndexManifest,
    build_cached_chain_index,
    chain_index_cache_path,
    load_cached_chain_index,
)
from .cli import default_user_cache_root
from .resource_cache import load_cached_ucsc_resource_bundle

_PROGRESS_BAR_WIDTH = 20


def main(argv: Sequence[str] | None = None) -> int:
    """Run ``prepare-liftassess-index`` and return its process exit code."""

    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return _run(args, stdout=sys.stdout, stderr=sys.stderr)
    except (ChainIndexCorruptionError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="prepare-liftassess-index",
        description=(
            "Build a reusable local chain index from an already verified liftAssess "
            "UCSC cache bundle. This command never contacts UCSC."
        ),
    )
    parser.add_argument(
        "source_db", help="source UCSC database identifier, e.g. canFam3"
    )
    parser.add_argument(
        "target_db", help="target UCSC database identifier, e.g. canFam4"
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="resource cache directory (default: platform user cache)",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help=(
            "discard an existing derived chain index for the exact cached chain and "
            "build it again"
        ),
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress nonessential status and progress messages",
    )
    return parser


def _run(
    args: argparse.Namespace,
    *,
    stdout: TextIO,
    stderr: TextIO,
) -> int:
    cache_root = args.cache_dir or default_user_cache_root()
    _status("Checking/verifying local UCSC cache...", quiet=args.quiet, stderr=stderr)
    bundle = load_cached_ucsc_resource_bundle(
        cache_root,
        args.source_db,
        args.target_db,
    )
    if bundle is None:
        print(
            "error: index preparation requires a complete verified cached UCSC bundle "
            f"for {args.source_db}→{args.target_db} under {cache_root}; run "
            "assess-liftover for this assembly pair first to acquire/verify resources",
            file=stderr,
        )
        return 1

    resource = bundle.chain
    index_path = chain_index_cache_path(cache_root, resource.sha256)

    if args.rebuild and index_path.exists():
        _status(
            f"Removing derived chain index before rebuild: {index_path}",
            quiet=args.quiet,
            stderr=stderr,
        )
        shutil.rmtree(index_path)
    elif not args.rebuild:
        try:
            existing = load_cached_chain_index(
                cache_root,
                resource,
                verify_database=True,
            )
        except ChainIndexCorruptionError as exc:
            print(
                "error: the cached chain index is unusable; rerun with --rebuild to "
                f"replace this derived artifact ({exc})",
                file=stderr,
            )
            return 1
        if existing is not None:
            _print_summary(
                existing.manifest,
                index_path=index_path,
                source_db=args.source_db,
                target_db=args.target_db,
                prepared=False,
                stdout=stdout,
            )
            return 0

    available_bytes = shutil.disk_usage(cache_root).free
    _status(
        f"Preparing reusable chain index for {args.source_db}→{args.target_db}...",
        quiet=args.quiet,
        stderr=stderr,
    )
    _status(
        f"Source chain: {_format_bytes(resource.size_bytes)}; "
        f"available filesystem space: {_format_bytes(available_bytes)}",
        quiet=args.quiet,
        stderr=stderr,
    )
    _status(
        "Index construction parses the complete chain once and may take many minutes "
        "for multi-gigabyte resources.",
        quiet=args.quiet,
        stderr=stderr,
    )

    progress_display = _IndexBuildProgressDisplay(
        total_bytes=resource.size_bytes,
        stderr=stderr,
    )
    progress_callback = None
    if not args.quiet and _is_interactive_terminal(stderr):
        progress_display.start()
        progress_callback = progress_display.update

    try:
        result = build_cached_chain_index(
            cache_root,
            resource,
            progress_callback=progress_callback,
        )
    except BaseException:
        if progress_callback is not None:
            progress_display.abort()
        raise
    if progress_callback is not None:
        progress_display.finish()

    _print_summary(
        result.manifest,
        index_path=result.index.root,
        source_db=args.source_db,
        target_db=args.target_db,
        prepared=True,
        stdout=stdout,
    )
    return 0


class _IndexBuildProgressDisplay:
    """Render exact raw source-byte progress while the index is built."""

    def __init__(self, *, total_bytes: int, stderr: TextIO) -> None:
        self._total_bytes = max(total_bytes, 0)
        self._stderr = stderr
        self._started = False
        self._last_percent = -1

    def start(self) -> None:
        self._render(0)
        self._started = True

    def update(self, bytes_read: int) -> None:
        bounded = min(max(bytes_read, 0), self._total_bytes)
        percent = _progress_percent(bounded, self._total_bytes)
        if percent == self._last_percent:
            return
        self._render(bounded)

    def finish(self) -> None:
        self._render(self._total_bytes)
        self._stderr.write("\n")
        self._stderr.flush()
        self._started = False

    def abort(self) -> None:
        if self._started:
            self._stderr.write("\n")
            self._stderr.flush()
            self._started = False

    def _render(self, bytes_read: int) -> None:
        percent = _progress_percent(bytes_read, self._total_bytes)
        if self._started:
            self._stderr.write("\r")
        filled = round(_PROGRESS_BAR_WIDTH * percent / 100)
        bar = "█" * filled + "—" * (_PROGRESS_BAR_WIDTH - filled)
        self._stderr.write(
            f"Indexing chain [{bar}] {percent:3d}% "
            f"({_format_bytes(bytes_read)} / {_format_bytes(self._total_bytes)})"
        )
        self._stderr.flush()
        self._last_percent = percent


def _print_summary(
    manifest: ChainIndexManifest,
    *,
    index_path: Path,
    source_db: str,
    target_db: str,
    prepared: bool,
    stdout: TextIO,
) -> None:
    state = "Prepared" if prepared else "Already prepared"
    total_size = (
        manifest.database_size_bytes
        + manifest.block_store_size_bytes
        + manifest.lookup_catalog_size_bytes
    )
    print(f"{state}: reusable chain index for {source_db}→{target_db}", file=stdout)
    print(f"Index: {index_path}", file=stdout)
    print(f"Source chain: {manifest.source_chain_sha256_identifier}", file=stdout)
    print(f"Records: {manifest.record_count}", file=stdout)
    print(f"Bin memberships: {manifest.membership_count}", file=stdout)
    print(f"Compressed blocks: {manifest.block_count}", file=stdout)
    print(f"Derived index size: {_format_bytes(total_size)}", file=stdout)


def _status(message: str, *, quiet: bool, stderr: TextIO) -> None:
    if not quiet:
        print(message, file=stderr)


def _is_interactive_terminal(stream: TextIO) -> bool:
    try:
        return stream.isatty()
    except (AttributeError, OSError):
        return False


def _progress_percent(bytes_read: int, total_bytes: int) -> int:
    if total_bytes <= 0:
        return 100 if bytes_read > 0 else 0
    return min(100, max(0, int(bytes_read * 100 / total_bytes)))


def _format_bytes(value: int) -> str:
    bounded = max(value, 0)
    gib = 1024**3
    mib = 1024**2
    kib = 1024
    if bounded >= gib:
        return f"{bounded / gib:.2f} GiB"
    if bounded >= mib:
        return f"{bounded / mib:.2f} MiB"
    if bounded >= kib:
        return f"{bounded / kib:.1f} KiB"
    return f"{bounded} B"


if __name__ == "__main__":
    raise SystemExit(main())
