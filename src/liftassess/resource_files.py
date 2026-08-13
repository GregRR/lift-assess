"""Streaming adapters for local UCSC chain/net resource files.

The parser and engine layers operate on iterables of text lines / parsed records so
that they remain independent of storage.  This module is the thin local-file boundary:
it opens plain-text or gzip-compressed resources, keeps the file handle alive only for
the duration of iteration, and feeds the existing parsers into the UCSC engine.

It deliberately does not download resources or infer provider terms. Callers supply
provenance describing the files they chose; ``provenance_source_for_file`` in the
resource-identity layer constructs the content-addressed file node, while callers still
provide any upstream alignment/process provenance that file bytes cannot establish. A
future downloader/cache layer can add recorded URLs, retrieval metadata, provider
checksum verification, and applicable terms without coupling those concerns to parsing.
"""

from __future__ import annotations

import gzip
import hashlib
import io
from collections.abc import Buffer, Iterator
from contextlib import contextmanager
from os import PathLike
from pathlib import Path
from typing import Protocol, TextIO, TypeAlias

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
from .resource_identity import (
    ResourceIdentityMismatchError,
    _sha256_checksum_from_file_provenance,
)

ResourcePath: TypeAlias = str | PathLike[str]

_CHUNK_SIZE = 1024 * 1024


class _Digest(Protocol):
    def update(self, data: bytes) -> None:
        ...

    def hexdigest(self) -> str:
        ...


class _HashingRawReader(io.RawIOBase):
    """Read one binary stream while hashing exactly the bytes returned upstream."""

    def __init__(self, raw: io.RawIOBase, digest: _Digest) -> None:
        super().__init__()
        self._raw = raw
        self._digest = digest

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Buffer) -> int:
        view = memoryview(buffer)
        data = self._raw.read(len(view))
        if not data:
            return 0
        self._digest.update(data)
        view[: len(data)] = data
        return len(data)

    def close(self) -> None:
        try:
            self._raw.close()
        finally:
            super().close()


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

    Every file provenance node used for parsing must carry one canonical SHA-256
    identifier. The raw bytes streamed into each consumed parser are checked against
    that digest before results can return, so provenance cannot silently describe a
    pre-mutation version of the file.

    Optional resource/provenance groups retain the engine's existing validation: net
    path + provenance must be supplied together, and reciprocal-best path + provenance
    + completeness must be supplied together.
    """

    chains = _iter_chain_file_with_provenance(chain_path, chain_provenance)
    # Unpaired optional inputs intentionally fall through to engine validation;
    # only a complete path/provenance pair may use the verified file stream.
    net_records = (
        _iter_net_file_with_provenance(net_path, net_provenance)
        if net_path is not None and net_provenance is not None
        else iter_net_file(net_path) if net_path is not None else None
    )
    reciprocal_best_chains = (
        _iter_chain_file_with_provenance(
            reciprocal_best_chain_path, reciprocal_best_provenance
        )
        if reciprocal_best_chain_path is not None
        and reciprocal_best_provenance is not None
        else iter_chain_file(reciprocal_best_chain_path)
        if reciprocal_best_chain_path is not None
        else None
    )

    return build_ucsc_candidates(
        source_interval,
        chains,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_records=net_records,
        net_provenance=net_provenance,
        reciprocal_best_chains=reciprocal_best_chains,
        reciprocal_best_provenance=reciprocal_best_provenance,
        reciprocal_best_completeness=reciprocal_best_completeness,
    )


def _iter_chain_file_with_provenance(
    path: ResourcePath,
    provenance: ProvenanceSource,
) -> Iterator[ChainRecord]:
    expected_sha256 = _sha256_checksum_from_file_provenance(provenance)
    with _open_text_resource(path, expected_sha256=expected_sha256) as lines:
        yield from iter_chain_records(lines)


def _iter_net_file_with_provenance(
    path: ResourcePath,
    provenance: ProvenanceSource,
) -> Iterator[NetRecord]:
    expected_sha256 = _sha256_checksum_from_file_provenance(provenance)
    with _open_text_resource(path, expected_sha256=expected_sha256) as lines:
        yield from iter_net_records(lines)


@contextmanager
def _open_text_resource(
    path: ResourcePath,
    *,
    expected_sha256: str | None = None,
) -> Iterator[TextIO]:
    """Open and stream one local UCSC text resource from a single file handle.

    Compression is detected from the first two bytes of that same open handle, so
    renamed resources do not depend on filename suffixes and there is no separate
    probe/reopen window.  When ``expected_sha256`` is supplied, the exact raw bytes
    that feed the decoder/decompressor are hashed as they are read.  Successful
    exhaustion of the parser therefore verifies that the provenance digest identifies
    the bytes actually assessed; a changed file fails loudly instead of returning
    candidates with stale provenance.
    """

    resource_path = Path(path)
    raw = io.FileIO(resource_path, mode="r")
    try:
        prefix = raw.read(2)
        raw.seek(0)
        is_gzip = prefix == b"\x1f\x8b"

        digest: _Digest | None = (
            hashlib.sha256() if expected_sha256 is not None else None
        )
        binary_raw: io.RawIOBase
        if digest is None:
            binary_raw = raw
        else:
            binary_raw = _HashingRawReader(raw, digest)

        buffered = io.BufferedReader(binary_raw)
        binary_text: io.BufferedIOBase
        if is_gzip:
            binary_text = gzip.GzipFile(fileobj=buffered, mode="rb")
        else:
            binary_text = buffered

        text = io.TextIOWrapper(binary_text, encoding="utf-8", newline="")
        try:
            yield text
            if digest is not None and expected_sha256 is not None:
                # Parser exhaustion should already have reached EOF.  Draining is
                # defensive and also includes any raw trailing bytes in the exact
                # artifact identity rather than only the decoded text payload.
                while text.read(_CHUNK_SIZE):
                    pass
                if is_gzip:
                    while buffered.read(_CHUNK_SIZE):
                        pass
                actual = digest.hexdigest()
                if actual != expected_sha256:
                    raise ResourceIdentityMismatchError(
                        f"SHA256 provenance mismatch for {resource_path}: "
                        f"expected {expected_sha256}, got {actual}"
                    )
        finally:
            text.close()
            if is_gzip:
                buffered.close()
    finally:
        raw.close()
