"""Streaming adapters for local UCSC chain/net resource files.

The parser and engine layers operate on iterables of text lines / parsed records so
that they remain independent of storage.  This module is the thin local-file boundary:
it opens plain-text or gzip-compressed resources, keeps the file handle alive only for
the duration of iteration, and feeds the existing parsers into the UCSC engine.

It deliberately does not download resources, infer provider terms, or manufacture
provenance.  Callers must supply provenance describing the files they chose; a future
downloader/cache layer can construct that provenance from recorded URLs, checksums,
retrieval metadata, and applicable provider terms without coupling those concerns to
parsing.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from os import PathLike
from pathlib import Path
from typing import TextIO, TypeAlias

from .chain import ChainRecord, iter_chain_records
from .engine import build_ucsc_candidates
from .models import (
    AssemblyIdentifier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestResourceCompleteness,
)
from .net import NetRecord, iter_net_records

ResourcePath: TypeAlias = str | PathLike[str]


def iter_chain_file(path: ResourcePath) -> Iterator[ChainRecord]:
    """Yield chain records from a local plain-text or gzip resource.

    The file is streamed; it is not read into memory as one object. Compression is
    detected from the gzip magic bytes rather than the filename so renamed local
    resources do not depend on a ``.gz`` suffix for correct decoding.
    """

    with _open_text_resource(path) as lines:
        yield from iter_chain_records(lines)


def iter_net_file(path: ResourcePath) -> Iterator[NetRecord]:
    """Yield net records from a local plain-text or gzip resource."""

    with _open_text_resource(path) as lines:
        yield from iter_net_records(lines)


def build_ucsc_candidates_from_files(
    source_interval: GenomicInterval,
    chain_path: ResourcePath,
    *,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_path: ResourcePath | None = None,
    net_provenance: ProvenanceSource | None = None,
    reciprocal_best_chain_path: ResourcePath | None = None,
    reciprocal_best_provenance: ProvenanceSource | None = None,
    reciprocal_best_completeness: ReciprocalBestResourceCompleteness | None = None,
) -> tuple[NormalizedCandidate, ...]:
    """Build UCSC candidates directly from local resource files.

    This is a storage adapter around :func:`build_ucsc_candidates`, not a second
    candidate-generation implementation.  The underlying engine still consumes each
    parser stream once across the full candidate set, so local gzip files are not
    rescanned once per candidate and are never materialized wholesale in memory.

    Optional resource/provenance groups retain the engine's existing validation: net
    path + provenance must be supplied together, and reciprocal-best path + provenance
    + completeness must be supplied together.
    """

    net_records = iter_net_file(net_path) if net_path is not None else None
    reciprocal_best_chains = (
        iter_chain_file(reciprocal_best_chain_path)
        if reciprocal_best_chain_path is not None
        else None
    )

    return build_ucsc_candidates(
        source_interval,
        iter_chain_file(chain_path),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_records=net_records,
        net_provenance=net_provenance,
        reciprocal_best_chains=reciprocal_best_chains,
        reciprocal_best_provenance=reciprocal_best_provenance,
        reciprocal_best_completeness=reciprocal_best_completeness,
    )


def _open_text_resource(path: ResourcePath) -> TextIO:
    """Open one local UCSC text resource without loading it into memory.

    UCSC currently distributes these resources as ``.gz`` files, but local filenames
    are not evidence of their actual encoding. Inspecting the gzip magic bytes keeps
    the storage boundary content-driven: a renamed gzip resource still decompresses,
    while an uncompressed file with a ``.gz`` suffix is still read as plain text.
    Only the two-byte signature is read before the file is reopened for streaming.
    """

    resource_path = Path(path)
    with resource_path.open(mode="rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"

    if is_gzip:
        return gzip.open(resource_path, mode="rt", encoding="utf-8", newline="")
    return resource_path.open(mode="rt", encoding="utf-8", newline="")
