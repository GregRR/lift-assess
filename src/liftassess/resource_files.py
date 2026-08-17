"""Streaming adapters for local UCSC chain/net resource files.

The parser and engine layers operate on iterables of text lines / parsed records so
that they remain independent of storage.  This module is the thin local-file boundary:
it opens plain-text or gzip-compressed resources, keeps the file handle alive only for
the duration of iteration, and feeds the existing parsers into the UCSC engine.

It deliberately does not download resources or infer provider terms. Callers may
supply local files plus provenance directly, or pass a fully acquired cache bundle
through the bridge below. In both cases upstream alignment/process provenance remains
explicit because file bytes and retrieval metadata cannot establish that dependency on
their own.
"""

from __future__ import annotations

import gzip
import hashlib
import io
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from os import PathLike
from pathlib import Path
from typing import Any, Protocol, TextIO, TypeAlias

from .chain import ChainRecord, iter_chain_records
from .engine import build_ucsc_candidates
from .models import (
    AssemblyIdentifier,
    EvidenceAvailabilityTier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceIdentifier,
    ProvenanceIdentifierKind,
    ProvenanceSource,
    ReciprocalBestResourceCompleteness,
)
from .net import NetRecord, iter_net_records
from .resource_cache import (
    CachedResource,
    CachedUCSCResourceBundle,
    UCSCBundleResourceRole,
)
from .resource_identity import (
    ResourceIdentityMismatchError,
    _sha256_checksum_from_file_provenance,
)

ResourcePath: TypeAlias = str | PathLike[str]
ResourceReadProgressCallback: TypeAlias = Callable[
    [UCSCBundleResourceRole, int, int], None
]
RawReadProgressCallback: TypeAlias = Callable[[int], None]

_CHUNK_SIZE = 1024 * 1024


class _Digest(Protocol):
    def update(self, data: bytes) -> None: ...

    def hexdigest(self) -> str: ...


class _HashingRawReader(io.RawIOBase):
    """Read one binary stream while hashing exactly the bytes returned upstream.

    Progress, when requested by the cached-bundle CLI path, is measured from these
    exact raw artifact bytes rather than decompressed records or a time estimate.
    Reporting is throttled by a caller-selected byte interval so terminal UI updates
    cannot dominate multi-gigabyte parsing.
    """

    def __init__(
        self,
        raw: io.RawIOBase,
        digest: _Digest,
        *,
        progress_callback: RawReadProgressCallback | None = None,
        progress_interval_bytes: int = _CHUNK_SIZE,
    ) -> None:
        super().__init__()
        if progress_interval_bytes <= 0:
            raise ValueError("progress interval must be positive")
        self._raw = raw
        self._digest = digest
        self._progress_callback = progress_callback
        self._progress_interval_bytes = progress_interval_bytes
        self._bytes_read = 0
        self._last_reported_bytes = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Any) -> int:
        # Python 3.11 has no public type for the general buffer protocol; using
        # Any here keeps this private RawIOBase adapter importable on the declared
        # Python 3.11 floor while preserving the runtime buffer check in memoryview.
        view = memoryview(buffer)
        data = self._raw.read(len(view))
        if not data:
            return 0
        self._digest.update(data)
        self._bytes_read += len(data)
        if (
            self._progress_callback is not None
            and self._bytes_read - self._last_reported_bytes
            >= self._progress_interval_bytes
        ):
            self._progress_callback(self._bytes_read)
            self._last_reported_bytes = self._bytes_read
        view[: len(data)] = data
        return len(data)

    def finish_progress(self) -> None:
        """Publish the exact final raw-byte count after identity verification."""

        if (
            self._progress_callback is not None
            and self._last_reported_bytes != self._bytes_read
        ):
            self._progress_callback(self._bytes_read)
            self._last_reported_bytes = self._bytes_read

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

    return _build_ucsc_candidates_from_files(
        source_interval,
        chain_path,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_path=net_path,
        net_provenance=net_provenance,
        reciprocal_best_chain_path=reciprocal_best_chain_path,
        reciprocal_best_provenance=reciprocal_best_provenance,
        reciprocal_best_completeness=reciprocal_best_completeness,
    )


def _build_ucsc_candidates_from_files(
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
    chain_progress: RawReadProgressCallback | None = None,
    chain_progress_interval_bytes: int = _CHUNK_SIZE,
    net_progress: RawReadProgressCallback | None = None,
    net_progress_interval_bytes: int = _CHUNK_SIZE,
    reciprocal_best_progress: RawReadProgressCallback | None = None,
    reciprocal_best_progress_interval_bytes: int = _CHUNK_SIZE,
) -> tuple[NormalizedCandidate, ...]:
    chains = _iter_chain_file_with_provenance(
        chain_path,
        chain_provenance,
        progress_callback=chain_progress,
        progress_interval_bytes=chain_progress_interval_bytes,
    )
    # Unpaired optional inputs intentionally fall through to engine validation;
    # only a complete path/provenance pair may use the verified file stream.
    net_records = (
        _iter_net_file_with_provenance(
            net_path,
            net_provenance,
            progress_callback=net_progress,
            progress_interval_bytes=net_progress_interval_bytes,
        )
        if net_path is not None and net_provenance is not None
        else iter_net_file(net_path)
        if net_path is not None
        else None
    )
    reciprocal_best_chains = (
        _iter_chain_file_with_provenance(
            reciprocal_best_chain_path,
            reciprocal_best_provenance,
            progress_callback=reciprocal_best_progress,
            progress_interval_bytes=reciprocal_best_progress_interval_bytes,
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


def build_ucsc_candidates_from_cached_bundle(
    source_interval: GenomicInterval,
    bundle: CachedUCSCResourceBundle,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    progress_callback: ResourceReadProgressCallback | None = None,
) -> tuple[NormalizedCandidate, ...]:
    """Build candidates from one fully acquired UCSC cache bundle.

    The cache already records the exact SHA-256 identity assigned at acquisition, so
    this bridge constructs content-addressed file provenance directly from those
    recorded digests instead of performing a separate pre-parse rehash. The existing
    file-backed parser still hashes every *consumed* raw file stream and checks it
    against that identity before candidates can return, preserving the mutation/TOCTOU
    protection at the scientific-use boundary.

    ``alignment_provenance`` remains caller-supplied. A cached URL plus file digest can
    identify the external artifacts liftAssess acquired, but cannot by itself establish
    the upstream alignment/process provenance needed to reason about evidence
    dependence.

    A complete COMPARATIVE cache bundle intentionally contains five provider files,
    while the current v1 candidate engine consumes three of them: the all-chain for
    candidate generation, the ordinary classified net for chain/net evidence, and the
    reciprocal-best chain for membership geometry. UCSC's current automation pipeline
    produces the ordinary net through chainNet/netSyntenic and netClass, then creates
    the optional ``*.syn.net.gz`` by filtering that ordinary net for synteny. Its
    reciprocal-best pipeline publishes both net and chain resources, while the v1
    membership implementation operates on the reciprocal-best chain geometry. The
    syntenic net and reciprocal-best net therefore remain available on ``bundle`` as
    retrieval/provenance context but are not silently substituted into parsers that do
    not use them.

    Primary implementation references checked 2026-08-14:
    https://raw.githubusercontent.com/ucscGenomeBrowser/kent/refs/heads/master/src/hg/utils/automation/doBlastzChainNet.pl
    https://raw.githubusercontent.com/ucscGenomeBrowser/kent/refs/heads/master/src/hg/utils/automation/doRecipBest.pl
    """

    _validate_cached_bundle_assemblies(
        source_interval,
        bundle,
        target_assembly=target_assembly,
    )

    chain_provenance = _cached_bundle_resource_provenance(
        bundle,
        UCSCBundleResourceRole.CHAIN,
        bundle.chain,
        alignment_provenance=alignment_provenance,
    )

    if bundle.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return _build_ucsc_candidates_from_files(
            source_interval,
            bundle.chain.path,
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
            chain_progress=_resource_progress_callback(
                progress_callback, UCSCBundleResourceRole.CHAIN, bundle.chain
            ),
            chain_progress_interval_bytes=_progress_interval_bytes(bundle.chain),
        )

    # CachedUCSCResourceBundle enforces the complete five-resource COMPARATIVE shape.
    # These assertions narrow the dataclass invariants for the type checker; only the
    # three resources below are direct inputs to the current candidate engine.
    assert bundle.net is not None
    assert bundle.reciprocal_best_chain is not None

    net_provenance = _cached_bundle_resource_provenance(
        bundle,
        UCSCBundleResourceRole.NET,
        bundle.net,
        alignment_provenance=alignment_provenance,
    )
    reciprocal_best_provenance = _cached_bundle_resource_provenance(
        bundle,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
        bundle.reciprocal_best_chain,
        alignment_provenance=alignment_provenance,
    )

    # Passing a CachedUCSCResourceBundle is the caller's claim that these cache
    # records represent the complete published bundle produced by the acquisition
    # boundary. As with the lower-level COMPLETE_RESOURCE API, liftAssess can verify
    # the bytes it consumes but cannot independently prove that a manually constructed
    # external object was not truncated before being described as complete.
    return _build_ucsc_candidates_from_files(
        source_interval,
        bundle.chain.path,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_path=bundle.net.path,
        net_provenance=net_provenance,
        reciprocal_best_chain_path=bundle.reciprocal_best_chain.path,
        reciprocal_best_provenance=reciprocal_best_provenance,
        reciprocal_best_completeness=(
            ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
        ),
        chain_progress=_resource_progress_callback(
            progress_callback, UCSCBundleResourceRole.CHAIN, bundle.chain
        ),
        chain_progress_interval_bytes=_progress_interval_bytes(bundle.chain),
        net_progress=_resource_progress_callback(
            progress_callback, UCSCBundleResourceRole.NET, bundle.net
        ),
        net_progress_interval_bytes=_progress_interval_bytes(bundle.net),
        reciprocal_best_progress=_resource_progress_callback(
            progress_callback,
            UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
            bundle.reciprocal_best_chain,
        ),
        reciprocal_best_progress_interval_bytes=_progress_interval_bytes(
            bundle.reciprocal_best_chain
        ),
    )


def _resource_progress_callback(
    callback: ResourceReadProgressCallback | None,
    role: UCSCBundleResourceRole,
    resource: CachedResource,
) -> RawReadProgressCallback | None:
    if callback is None:
        return None

    def report(bytes_read: int) -> None:
        callback(role, bytes_read, resource.size_bytes)

    return report


def _progress_interval_bytes(resource: CachedResource) -> int:
    # Aim for roughly one hundred measured updates while keeping small resources
    # responsive and avoiding a callback for every decompressor read.
    return max(64 * 1024, resource.size_bytes // 100)


def _validate_cached_bundle_assemblies(
    source_interval: GenomicInterval,
    bundle: CachedUCSCResourceBundle,
    *,
    target_assembly: AssemblyIdentifier,
) -> None:
    if not _assembly_represents_ucsc_db(source_interval.assembly, bundle.source_db):
        raise ValueError(
            "source interval assembly does not represent cached bundle source db "
            f"{bundle.source_db!r}"
        )
    if not _assembly_represents_ucsc_db(target_assembly, bundle.target_db):
        raise ValueError(
            "target assembly does not represent cached bundle target db "
            f"{bundle.target_db!r}"
        )


def _assembly_represents_ucsc_db(
    assembly: AssemblyIdentifier,
    db: str,
) -> bool:
    """Match only an explicitly recorded UCSC db name/alias; do no alias resolution."""

    return db == assembly.name or db in assembly.aliases


def _cached_bundle_resource_provenance(
    bundle: CachedUCSCResourceBundle,
    role: UCSCBundleResourceRole,
    resource: CachedResource,
    *,
    alignment_provenance: ProvenanceSource,
) -> ProvenanceSource:
    """Build the canonical file-provenance node for one consumed bundle resource.

    Candidate evidence and assessment reporting must refer to the same structural
    provenance node for the same cached bytes. Keeping the role-specific labels here
    prevents the two boundaries from silently drifting apart.
    """

    if role is UCSCBundleResourceRole.CHAIN:
        label = f"UCSC {bundle.source_db}→{bundle.target_db} chain resource"
    elif role is UCSCBundleResourceRole.NET:
        label = f"UCSC {bundle.source_db}→{bundle.target_db} net resource"
    elif role is UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN:
        label = (
            f"UCSC {bundle.source_db}→{bundle.target_db} reciprocal-best chain resource"
        )
    else:
        raise ValueError("bundle resource role is not consumed by the v1 engine")

    return _provenance_for_cached_resource(
        resource,
        label=label,
        derived_from=(alignment_provenance,),
    )


def _provenance_for_cached_resource(
    resource: CachedResource,
    *,
    label: str,
    derived_from: tuple[ProvenanceSource, ...],
) -> ProvenanceSource:
    identifier = ProvenanceIdentifier(
        kind=ProvenanceIdentifierKind.SHA256,
        value=resource.sha256,
    )
    return ProvenanceSource(
        source_id=f"file:{identifier.value}",
        label=label,
        identifiers=(identifier,),
        derived_from=derived_from,
    )


def _iter_chain_file_with_provenance(
    path: ResourcePath,
    provenance: ProvenanceSource,
    *,
    progress_callback: RawReadProgressCallback | None = None,
    progress_interval_bytes: int = _CHUNK_SIZE,
) -> Iterator[ChainRecord]:
    expected_sha256 = _sha256_checksum_from_file_provenance(provenance)
    with _open_text_resource(
        path,
        expected_sha256=expected_sha256,
        progress_callback=progress_callback,
        progress_interval_bytes=progress_interval_bytes,
    ) as lines:
        yield from iter_chain_records(lines)


def _iter_net_file_with_provenance(
    path: ResourcePath,
    provenance: ProvenanceSource,
    *,
    progress_callback: RawReadProgressCallback | None = None,
    progress_interval_bytes: int = _CHUNK_SIZE,
) -> Iterator[NetRecord]:
    expected_sha256 = _sha256_checksum_from_file_provenance(provenance)
    with _open_text_resource(
        path,
        expected_sha256=expected_sha256,
        progress_callback=progress_callback,
        progress_interval_bytes=progress_interval_bytes,
    ) as lines:
        yield from iter_net_records(lines)


@contextmanager
def _open_text_resource(
    path: ResourcePath,
    *,
    expected_sha256: str | None = None,
    progress_callback: RawReadProgressCallback | None = None,
    progress_interval_bytes: int = _CHUNK_SIZE,
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
        hashing_raw: _HashingRawReader | None = None
        if digest is None:
            if progress_callback is not None:
                raise ValueError(
                    "raw-byte progress requires an expected SHA256 identity"
                )
            binary_raw = raw
        else:
            hashing_raw = _HashingRawReader(
                raw,
                digest,
                progress_callback=progress_callback,
                progress_interval_bytes=progress_interval_bytes,
            )
            binary_raw = hashing_raw

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
                assert hashing_raw is not None
                hashing_raw.finish_progress()
        finally:
            text.close()
            if is_gzip:
                buffered.close()
    finally:
        raw.close()
