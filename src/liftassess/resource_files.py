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

from collections.abc import Callable, Iterable, Iterator
from typing import TypeAlias

from .chain import ChainRecord, iter_chain_records
from .chain_index import ChainIndex
from .engine import build_ucsc_candidates, build_ucsc_chain_candidates_for_intervals
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
    CachedUCSCChainResource,
    CachedUCSCResourceBundle,
    UCSCBundleResourceRole,
)
from .resource_identity import _sha256_checksum_from_file_provenance
from .resource_stream import (
    RawReadProgressCallback,
    ResourcePath,
    open_text_resource,
)

ResourceReadProgressCallback: TypeAlias = Callable[
    [UCSCBundleResourceRole, int, int], None
]

_CHUNK_SIZE = 1024 * 1024


def iter_chain_file(path: ResourcePath) -> Iterator[ChainRecord]:
    """Yield chain records from a local plain-text or gzip resource.

    The file is streamed; it is not read into memory as one object. Compression is
    detected from the gzip magic bytes rather than the filename so renamed local
    resources do not depend on a ``.gz`` suffix for correct decoding.
    """

    with open_text_resource(path) as lines:
        yield from iter_chain_records(lines)


def iter_net_file(path: ResourcePath) -> Iterator[NetRecord]:
    """Yield net records from a local plain-text or gzip resource."""

    with open_text_resource(path) as lines:
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
    return _build_ucsc_candidates_with_chain_records(
        source_interval,
        chains,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_path=net_path,
        net_provenance=net_provenance,
        reciprocal_best_chain_path=reciprocal_best_chain_path,
        reciprocal_best_provenance=reciprocal_best_provenance,
        reciprocal_best_completeness=reciprocal_best_completeness,
        net_progress=net_progress,
        net_progress_interval_bytes=net_progress_interval_bytes,
        reciprocal_best_progress=reciprocal_best_progress,
        reciprocal_best_progress_interval_bytes=reciprocal_best_progress_interval_bytes,
    )


def _build_ucsc_candidates_with_chain_records(
    source_interval: GenomicInterval,
    chains: Iterable[ChainRecord],
    *,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_path: ResourcePath | None = None,
    net_provenance: ProvenanceSource | None = None,
    reciprocal_best_chain_path: ResourcePath | None = None,
    reciprocal_best_provenance: ProvenanceSource | None = None,
    reciprocal_best_completeness: ReciprocalBestResourceCompleteness | None = None,
    net_progress: RawReadProgressCallback | None = None,
    net_progress_interval_bytes: int = _CHUNK_SIZE,
    reciprocal_best_progress: RawReadProgressCallback | None = None,
    reciprocal_best_progress_interval_bytes: int = _CHUNK_SIZE,
) -> tuple[NormalizedCandidate, ...]:
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
    chain_index: ChainIndex | None = None,
) -> tuple[NormalizedCandidate, ...]:
    """Build candidates from one fully acquired UCSC cache bundle.

    The cache already records the exact SHA-256 identity assigned at acquisition, so
    this bridge constructs content-addressed file provenance directly from those
    recorded digests instead of performing a separate pre-parse rehash. Without an
    index, the existing file-backed parser hashes every consumed raw file stream and
    checks it against that identity before candidates can return. With a validated
    chain index, candidate generation instead consumes exact chain records from the
    derived artifact while preserving the original chain provenance; index metadata is
    bound to the cached chain SHA-256/size and selected compressed blocks are verified
    before parsing.

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
    indexed_chains: tuple[ChainRecord, ...] | None = None
    if chain_index is not None:
        _validate_chain_index_resource(chain_index, bundle.chain)
        indexed_chains = chain_index.records_for_interval(source_interval)

    if bundle.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        if indexed_chains is not None:
            return build_ucsc_candidates(
                source_interval,
                indexed_chains,
                target_assembly=target_assembly,
                chain_provenance=chain_provenance,
            )
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
    if indexed_chains is not None:
        return _build_ucsc_candidates_with_chain_records(
            source_interval,
            indexed_chains,
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
            net_path=bundle.net.path,
            net_provenance=net_provenance,
            reciprocal_best_chain_path=bundle.reciprocal_best_chain.path,
            reciprocal_best_provenance=reciprocal_best_provenance,
            reciprocal_best_completeness=(
                ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
            ),
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


def build_ucsc_chain_candidates_for_intervals_from_cached_chain(
    source_intervals: Iterable[GenomicInterval],
    chain_context: CachedUCSCChainResource,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    progress_callback: ResourceReadProgressCallback | None = None,
    chain_index: ChainIndex | None = None,
) -> tuple[tuple[NormalizedCandidate, ...], ...]:
    """Build chain-only candidates for many intervals from one cached chain.

    With an index, each interval uses region-addressable lookup. Without one, the
    original chain is verified and parsed exactly once across the full interval set.
    """

    intervals = tuple(source_intervals)
    if not intervals:
        return ()
    for interval in intervals:
        if not _assembly_represents_ucsc_db(interval.assembly, chain_context.source_db):
            raise ValueError(
                "source interval assembly does not represent cached chain source db "
                f"{chain_context.source_db!r}"
            )
    if not _assembly_represents_ucsc_db(target_assembly, chain_context.target_db):
        raise ValueError(
            "target assembly does not represent cached chain target db "
            f"{chain_context.target_db!r}"
        )
    chain_provenance = _cached_chain_resource_provenance(
        chain_context,
        alignment_provenance=alignment_provenance,
    )
    if chain_index is not None:
        _validate_chain_index_resource(chain_index, chain_context.chain)
        return tuple(
            build_ucsc_candidates(
                interval,
                chain_index.records_for_interval(interval),
                target_assembly=target_assembly,
                chain_provenance=chain_provenance,
            )
            for interval in intervals
        )

    chains = _iter_chain_file_with_provenance(
        chain_context.chain.path,
        chain_provenance,
        progress_callback=_resource_progress_callback(
            progress_callback, UCSCBundleResourceRole.CHAIN, chain_context.chain
        ),
        progress_interval_bytes=_progress_interval_bytes(chain_context.chain),
    )
    return build_ucsc_chain_candidates_for_intervals(
        intervals,
        chains,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
    )


def build_ucsc_chain_candidates_for_intervals_from_cached_bundle(
    source_intervals: Iterable[GenomicInterval],
    bundle: CachedUCSCResourceBundle,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    progress_callback: ResourceReadProgressCallback | None = None,
    chain_index: ChainIndex | None = None,
) -> tuple[tuple[NormalizedCandidate, ...], ...]:
    """Compatibility wrapper using only ``bundle.chain`` for reverse execution."""

    return build_ucsc_chain_candidates_for_intervals_from_cached_chain(
        source_intervals,
        CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        target_assembly=target_assembly,
        alignment_provenance=alignment_provenance,
        progress_callback=progress_callback,
        chain_index=chain_index,
    )


def _validate_chain_index_resource(
    chain_index: ChainIndex, resource: CachedResource
) -> None:
    if chain_index.manifest.source_chain_sha256_identifier != resource.sha256:
        raise ValueError("chain index source identity does not match cached chain")
    if chain_index.manifest.source_chain_size_bytes != resource.size_bytes:
        raise ValueError("chain index source size does not match cached chain")


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


def _cached_chain_resource_provenance(
    chain_context: CachedUCSCChainResource,
    *,
    alignment_provenance: ProvenanceSource,
) -> ProvenanceSource:
    """Build canonical provenance for one consumed directional chain."""

    return _provenance_for_cached_resource(
        chain_context.chain,
        label=(
            f"UCSC {chain_context.source_db}→{chain_context.target_db} chain resource"
        ),
        derived_from=(alignment_provenance,),
    )


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
    with open_text_resource(
        path,
        expected_sha256_hex=expected_sha256,
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
    with open_text_resource(
        path,
        expected_sha256_hex=expected_sha256,
        progress_callback=progress_callback,
        progress_interval_bytes=progress_interval_bytes,
    ) as lines:
        yield from iter_net_records(lines)
