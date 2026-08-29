"""Acquire UCSC resources into an explicit content-addressed local cache.

Discovery answers which provider resources exist; this module handles the separate
act of planning and retrieving already-discovered UCSC URLs.  Retrieval is deliberately
explicit about provider terms, streams bytes without materializing large resources,
verifies UCSC-published MD5 metadata when available, and stores the resulting exact
artifact by liftAssess's canonical SHA-256 identity.

The cache root is always supplied by the caller.  liftAssess does not create a cache
inside the source tree; the CLI supplies the platform-specific default cache root at
its boundary.

Primary provider references checked 2026-08-13:

- UCSC data/software licensing: https://genome.ucsc.edu/license/
- restricted canFam3 liftOver-chain terms: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/
- canFam3/canFam4 comparative terms and files: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/
- comparative MD5 metadata: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt
- reciprocal-best MD5 metadata: https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/md5sum.txt
- UCSC download guidance: https://genome.ucsc.edu/goldenpath/help/ftp.html

The canFam3/canFam4 comparison demonstrates why checksum lookup is exact-filename
based: its MD5 file covers some comparison artifacts but not ``canFam3.canFam4.net.gz``.
Remote metadata inspection uses body-free HTTP HEAD requests and preserves only headers
actually advertised by the provider.  Live UCSC checks on 2026-08-14 verified byte-range,
``Content-Range``, strong-ETag, and ``If-Range`` behavior for the comparative fixture.
Acquisition can therefore retain and resume validator-bound partial HTTPS transfers when
those exact preconditions are advertised, while falling back to the original fresh streaming
path when resumable metadata is unavailable.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import Enum
from pathlib import Path, PurePosixPath
from types import TracebackType
from typing import BinaryIO, Protocol, Self, TypeAlias, cast
from urllib.error import HTTPError, URLError
from urllib.parse import SplitResult, urljoin, urlsplit, urlunsplit
from urllib.request import Request, urlopen

from ._version import __version__
from .models import EvidenceAvailabilityTier
from .resource_identity import (
    ResourceChecksumAlgorithm,
    ResourceChecksumMismatchError,
    compute_resource_checksum,
    sha256_hex_from_identifier,
)
from .resources import UCSCResourceBundle

ResourcePath: TypeAlias = str | os.PathLike[str]
CacheVerificationProgressCallback: TypeAlias = Callable[[int, int, bool], None]
ResourceTransferProgressCallback: TypeAlias = Callable[[int, int | None], None]

_UCSC_HOSTS = frozenset({"hgdownload.soe.ucsc.edu", "hgdownload.gi.ucsc.edu"})
_UCSC_LICENSE_URL = "https://genome.ucsc.edu/license/"
_USER_AGENT = f"liftAssess/{__version__} resource-acquisition"
_CHUNK_SIZE = 1024 * 1024
_MD5_RE = re.compile(r"^[0-9a-fA-F]{32}$")
_CACHE_SCHEMA_VERSION = 1
_CONTENT_RANGE_RE = re.compile(r"^bytes (\d+)-(\d+)/(\d+)$")


class UCSCResourceAcquisitionError(RuntimeError):
    """A UCSC resource could not be acquired or cached reliably."""


class UCSCResourceTermsAcknowledgementRequired(UCSCResourceAcquisitionError):
    """Retrieval was requested without explicit acknowledgement of provider terms."""


class UCSCBundleAcquisitionPlanAcknowledgementRequired(UCSCResourceAcquisitionError):
    """Bundle retrieval was requested without acknowledging its explicit transfer plan."""


class _ResumeRestartRequired(RuntimeError):
    """A partial download cannot be safely extended and must restart fresh."""


class UCSCResourceClass(str, Enum):
    """UCSC publication classes currently produced by the resource resolver."""

    ASSEMBLY_METADATA = "ASSEMBLY_METADATA"
    COMPARATIVE = "COMPARATIVE"
    LIFTOVER_CHAIN = "LIFTOVER_CHAIN"


class UCSCBundleResourceRole(str, Enum):
    """One resource role in a discovered UCSC evidence bundle."""

    CHAIN = "CHAIN"
    NET = "NET"
    SYNTENIC_NET = "SYNTENIC_NET"
    RECIPROCAL_BEST_CHAIN = "RECIPROCAL_BEST_CHAIN"
    RECIPROCAL_BEST_NET = "RECIPROCAL_BEST_NET"


UCSCBundleTransferProgressCallback: TypeAlias = Callable[
    [UCSCBundleResourceRole, int, int | None, bool], None
]


@dataclass(frozen=True)
class UCSCResourceTerms:
    """Provider terms references that must be surfaced before retrieval.

    ``restricted_liftover_chain`` distinguishes UCSC's dedicated
    ``liftOver/*.over.chain.gz`` files, for which UCSC currently states that
    downloading/using indicates EULA acceptance and that free use is limited to the
    described non-commercial/nonprofit cases.  Comparative resources are a distinct
    publication class and retain their own directory terms URL instead. UCSC database
    table dumps used for assembly metadata are separately documented as freely usable
    and therefore do not require an acknowledgement gate before retrieval.
    """

    resource_class: UCSCResourceClass
    general_terms_url: str
    directory_terms_url: str
    restricted_liftover_chain: bool


@dataclass(frozen=True)
class ProviderChecksum:
    """Provider-published integrity metadata for one downloaded file."""

    algorithm: ResourceChecksumAlgorithm
    value: str
    source_url: str


@dataclass(frozen=True)
class CachedResource:
    """One exact UCSC artifact stored in the caller-supplied cache."""

    path: Path
    source_url: str
    retrieved_at: str
    sha256: str
    size_bytes: int
    provider_checksum: ProviderChecksum | None
    terms: UCSCResourceTerms
    cache_hit: bool


@dataclass(frozen=True)
class CachedUCSCChainResource:
    """One cached directional chain with its exact publication class.

    Actual reverse mapping consumes only a chain resource. Keeping the UCSC database
    direction and evidence tier attached prevents comparative all-chain and filtered
    liftOver chains from being substituted silently for one another.
    """

    source_db: str
    target_db: str
    evidence_tier: EvidenceAvailabilityTier
    chain: CachedResource

    def __post_init__(self) -> None:
        _validate_bundle_resource_binding(
            UCSCBundleResourceRole.CHAIN,
            self.chain.source_url,
            source_db=self.source_db,
            target_db=self.target_db,
            evidence_tier=self.evidence_tier,
        )


@dataclass(frozen=True)
class _CachedResourceIndexEntry:
    """Structurally valid cache-index metadata awaiting SHA-256 verification."""

    artifact_path: Path
    retrieved_at: str
    sha256_hex: str
    size_bytes: int
    provider_checksum: ProviderChecksum | None


@dataclass(frozen=True)
class UCSCRemoteResourceMetadata:
    """Metadata returned by a body-free HTTP HEAD request for one UCSC resource.

    ``content_length_bytes`` records the provider's HTTP ``Content-Length`` when it is
    present; it is not guessed from directory-listing display text.  Inspection requests
    identity encoding.  If a provider nevertheless reports a non-identity
    ``Content-Encoding``, the raw header value is preserved but excluded from bundle
    transfer-size totals because it may not describe the cached resource bytes.
    ``accept_ranges`` is preserved as advertised and is not, by itself, treated as proof
    that resumable HTTP acquisition is safe.
    """

    url: str
    terms: UCSCResourceTerms
    content_length_bytes: int | None
    accept_ranges: str | None
    last_modified: str | None
    etag: str | None
    content_encoding: str | None

    def __post_init__(self) -> None:
        expected_terms = ucsc_resource_terms(self.url)
        if self.terms != expected_terms:
            raise ValueError(
                "remote resource metadata terms must match the resource URL classification"
            )
        if self.content_length_bytes is not None and self.content_length_bytes < 0:
            raise ValueError("remote resource Content-Length cannot be negative")


@dataclass(frozen=True)
class UCSCBundleTransferInspectionItem:
    """Remote metadata for one exact role/URL in a bundle acquisition plan."""

    role: UCSCBundleResourceRole
    metadata: UCSCRemoteResourceMetadata

    def __post_init__(self) -> None:
        _validate_resource_role_filename(self.role, self.metadata.url)

    @property
    def identity_content_length_bytes(self) -> int | None:
        """Return the usable exact-byte size advertised for this resource, if any."""

        return _identity_content_length_bytes(self.metadata)


@dataclass(frozen=True)
class UCSCBundleTransferInspection:
    """Body-free remote metadata inspection for an existing bundle plan."""

    source_db: str
    target_db: str
    evidence_tier: EvidenceAvailabilityTier
    items: tuple[UCSCBundleTransferInspectionItem, ...]

    def __post_init__(self) -> None:
        expected_roles = _bundle_roles_for_tier(self.evidence_tier)
        actual_roles = tuple(item.role for item in self.items)
        if actual_roles != expected_roles:
            raise ValueError(
                "bundle transfer inspection must contain the exact ordered resource roles "
                f"for {self.evidence_tier.value}: expected {expected_roles}, got "
                f"{actual_roles}"
            )
        for item in self.items:
            _validate_bundle_resource_binding(
                item.role,
                item.metadata.url,
                source_db=self.source_db,
                target_db=self.target_db,
                evidence_tier=self.evidence_tier,
            )

    @property
    def known_content_length_bytes(self) -> int:
        """Sum usable provider Content-Length values for identity-encoded resources."""

        return sum(
            _identity_content_length_bytes(item.metadata) or 0 for item in self.items
        )

    @property
    def total_content_length_bytes(self) -> int | None:
        """Return the complete identity-encoded bundle size only when fully known."""

        lengths = tuple(
            _identity_content_length_bytes(item.metadata) for item in self.items
        )
        if any(length is None for length in lengths):
            return None
        return sum(cast(int, length) for length in lengths)


def _identity_content_length_bytes(
    metadata: UCSCRemoteResourceMetadata,
) -> int | None:
    if (
        metadata.content_encoding is not None
        and metadata.content_encoding.casefold() != "identity"
    ):
        return None
    return metadata.content_length_bytes


@dataclass(frozen=True)
class UCSCBundleAcquisitionItem:
    """One inspectable resource entry in a bundle transfer plan."""

    role: UCSCBundleResourceRole
    url: str
    terms: UCSCResourceTerms

    def __post_init__(self) -> None:
        expected_terms = ucsc_resource_terms(self.url)
        if self.terms != expected_terms:
            raise ValueError(
                "bundle acquisition item terms must match the resource URL "
                "classification"
            )
        _validate_resource_role_filename(self.role, self.url)


@dataclass(frozen=True)
class UCSCBundleAcquisitionPlan:
    """Explicit pre-transfer plan for one discovered UCSC resource bundle."""

    source_db: str
    target_db: str
    evidence_tier: EvidenceAvailabilityTier
    items: tuple[UCSCBundleAcquisitionItem, ...]

    def __post_init__(self) -> None:
        expected_roles = _bundle_roles_for_tier(self.evidence_tier)
        actual_roles = tuple(item.role for item in self.items)
        if actual_roles != expected_roles:
            raise ValueError(
                "bundle acquisition plan must contain the exact ordered resource roles "
                f"for {self.evidence_tier.value}: expected {expected_roles}, got "
                f"{actual_roles}"
            )
        for item in self.items:
            _validate_bundle_resource_binding(
                item.role,
                item.url,
                source_db=self.source_db,
                target_db=self.target_db,
                evidence_tier=self.evidence_tier,
            )


@dataclass(frozen=True)
class CachedUCSCResourceBundle:
    """Complete local cache records for one discovered UCSC evidence bundle."""

    source_db: str
    target_db: str
    evidence_tier: EvidenceAvailabilityTier
    chain: CachedResource
    net: CachedResource | None = None
    syntenic_net: CachedResource | None = None
    reciprocal_best_chain: CachedResource | None = None
    reciprocal_best_net: CachedResource | None = None

    def __post_init__(self) -> None:
        comparative = (
            self.net,
            self.syntenic_net,
            self.reciprocal_best_chain,
            self.reciprocal_best_net,
        )
        if self.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
            if any(resource is None for resource in comparative):
                raise ValueError(
                    "COMPARATIVE cached bundle requires net, syntenic net, and "
                    "reciprocal-best chain/net resources"
                )
        elif any(resource is not None for resource in comparative):
            raise ValueError(
                "LIFTOVER_ONLY cached bundle cannot carry comparative resources"
            )

        resources_by_role = (
            (UCSCBundleResourceRole.CHAIN, self.chain),
            (UCSCBundleResourceRole.NET, self.net),
            (UCSCBundleResourceRole.SYNTENIC_NET, self.syntenic_net),
            (UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN, self.reciprocal_best_chain),
            (UCSCBundleResourceRole.RECIPROCAL_BEST_NET, self.reciprocal_best_net),
        )
        for role, resource in resources_by_role:
            if resource is None:
                continue
            _validate_bundle_resource_binding(
                role,
                resource.source_url,
                source_db=self.source_db,
                target_db=self.target_db,
                evidence_tier=self.evidence_tier,
            )


class _BinaryResponse(Protocol):
    status: int

    def read(self, size: int = -1) -> bytes: ...

    def getheader(self, name: str, default: str | None = None) -> str | None: ...

    def __enter__(self) -> Self: ...

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None: ...


URLopener = Callable[[Request], _BinaryResponse]
Clock = Callable[[], datetime]


class _ResourceAcquirer(Protocol):
    def __call__(
        self,
        url: str,
        cache_root: ResourcePath,
        *,
        terms_acknowledged: bool,
        refresh: bool = False,
        progress_callback: ResourceTransferProgressCallback | None = None,
    ) -> CachedResource: ...


def ucsc_resource_terms(url: str) -> UCSCResourceTerms:
    """Return the terms references for one resolver-produced UCSC resource URL."""

    parts = _validate_ucsc_resource_url(url)
    segments = parts.path.strip("/").split("/")
    filename = PurePosixPath(parts.path).name

    if (
        len(segments) == 4
        and segments[0] == "goldenPath"
        and segments[2] == "database"
        and filename in {"chromInfo.txt.gz", "chromAlias.txt.gz"}
    ):
        return UCSCResourceTerms(
            resource_class=UCSCResourceClass.ASSEMBLY_METADATA,
            general_terms_url=_UCSC_LICENSE_URL,
            directory_terms_url=_directory_url(parts, segments, 2),
            restricted_liftover_chain=False,
        )

    if "liftOver" in segments[:-1] and filename.endswith(".over.chain.gz"):
        directory_index = segments.index("liftOver")
        return UCSCResourceTerms(
            resource_class=UCSCResourceClass.LIFTOVER_CHAIN,
            general_terms_url=_UCSC_LICENSE_URL,
            directory_terms_url=_directory_url(parts, segments, directory_index),
            restricted_liftover_chain=True,
        )

    comparative_indexes = [
        index
        for index, segment in enumerate(segments[:-1])
        if segment.startswith("vs") and len(segment) > 2
    ]
    if comparative_indexes:
        # Terms for reciprocalBest children live on the comparison directory's
        # README/index, not necessarily on the immediate subdirectory page.
        directory_index = comparative_indexes[0]
        return UCSCResourceTerms(
            resource_class=UCSCResourceClass.COMPARATIVE,
            general_terms_url=_UCSC_LICENSE_URL,
            directory_terms_url=_directory_url(parts, segments, directory_index),
            restricted_liftover_chain=False,
        )

    raise ValueError(
        "unsupported UCSC resource URL outside assembly-metadata/comparative/liftOver "
        "paths"
    )


def plan_ucsc_bundle_acquisition(
    bundle: UCSCResourceBundle,
) -> UCSCBundleAcquisitionPlan:
    """Build an inspectable no-network transfer plan from a discovered bundle.

    Planning does not create the cache or contact UCSC.  It enumerates the exact URLs
    required by the bundle's evidence-availability tier and records the applicable
    provider terms for each resource.  After reviewing those terms, callers can pass
    the resulting plan to ``inspect_ucsc_bundle_transfer_plan`` with explicit terms
    acknowledgement for body-free provider metadata before separately acknowledging
    and executing the transfer plan.
    """

    urls: tuple[tuple[UCSCBundleResourceRole, str], ...]
    if bundle.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        net_url = bundle.net_url
        syntenic_net_url = bundle.syntenic_net_url
        reciprocal_best_chain_url = bundle.reciprocal_best_chain_url
        reciprocal_best_net_url = bundle.reciprocal_best_net_url
        if (
            net_url is None
            or syntenic_net_url is None
            or reciprocal_best_chain_url is None
            or reciprocal_best_net_url is None
        ):
            raise ValueError("COMPARATIVE discovered bundle is incomplete")
        urls = (
            (UCSCBundleResourceRole.CHAIN, bundle.chain_url),
            (UCSCBundleResourceRole.NET, net_url),
            (UCSCBundleResourceRole.SYNTENIC_NET, syntenic_net_url),
            (UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN, reciprocal_best_chain_url),
            (UCSCBundleResourceRole.RECIPROCAL_BEST_NET, reciprocal_best_net_url),
        )
    else:
        urls = ((UCSCBundleResourceRole.CHAIN, bundle.chain_url),)

    return UCSCBundleAcquisitionPlan(
        source_db=bundle.source_db,
        target_db=bundle.target_db,
        evidence_tier=bundle.evidence_tier,
        items=tuple(
            UCSCBundleAcquisitionItem(
                role=role,
                url=url,
                terms=ucsc_resource_terms(url),
            )
            for role, url in urls
        ),
    )


def inspect_ucsc_resource(
    url: str,
    *,
    terms_acknowledged: bool,
) -> UCSCRemoteResourceMetadata:
    """Inspect one UCSC resource with an HTTP HEAD request and no body transfer.

    This is metadata inspection, not acquisition: it does not create a cache or download
    resource bytes.  Because it still contacts the provider, it requires the same
    explicit terms acknowledgement as resource acquisition.  Missing HTTP metadata is
    preserved as ``None`` rather than inferred from a directory listing or filename.
    """

    return _inspect_ucsc_resource(
        url,
        terms_acknowledged=terms_acknowledged,
        open_url=_open_url,
    )


def _inspect_ucsc_resource(
    url: str,
    *,
    terms_acknowledged: bool,
    open_url: URLopener,
) -> UCSCRemoteResourceMetadata:
    terms = ucsc_resource_terms(url)
    _require_terms_acknowledgement(terms, terms_acknowledged=terms_acknowledged)
    request = Request(
        url,
        headers={
            "User-Agent": _USER_AGENT,
            "Accept-Encoding": "identity",
        },
        method="HEAD",
    )
    try:
        with open_url(request) as response:
            content_length = _response_content_length(response)
            accept_ranges = _optional_response_header(response, "Accept-Ranges")
            last_modified = _optional_response_header(response, "Last-Modified")
            etag = _optional_response_header(response, "ETag")
            content_encoding = _optional_response_header(response, "Content-Encoding")
    except HTTPError as exc:
        raise UCSCResourceAcquisitionError(
            f"failed to inspect UCSC resource metadata {url}: HTTP {exc.code}"
        ) from exc
    except (URLError, TimeoutError, OSError) as exc:
        raise UCSCResourceAcquisitionError(
            f"failed to inspect UCSC resource metadata {url}: {exc}"
        ) from exc

    return UCSCRemoteResourceMetadata(
        url=url,
        terms=terms,
        content_length_bytes=content_length,
        accept_ranges=accept_ranges,
        last_modified=last_modified,
        etag=etag,
        content_encoding=content_encoding,
    )


def inspect_ucsc_bundle_transfer_plan(
    plan: UCSCBundleAcquisitionPlan,
    *,
    terms_acknowledged: bool,
) -> UCSCBundleTransferInspection:
    """Inspect remote metadata for every item in a bundle plan without GET bodies.

    The existing no-network plan remains the source of role/URL/terms intent.  After
    explicit terms acknowledgement, this separate step contacts each exact URL with HEAD
    and returns only provider-advertised HTTP metadata.  It does not execute or
    acknowledge the transfer plan.
    """

    return _inspect_ucsc_bundle_transfer_plan(
        plan,
        terms_acknowledged=terms_acknowledged,
        inspect_resource=lambda url: inspect_ucsc_resource(
            url,
            terms_acknowledged=True,
        ),
    )


class _ResourceInspector(Protocol):
    def __call__(self, url: str) -> UCSCRemoteResourceMetadata: ...


def _inspect_ucsc_bundle_transfer_plan(
    plan: UCSCBundleAcquisitionPlan,
    *,
    terms_acknowledged: bool,
    inspect_resource: _ResourceInspector,
) -> UCSCBundleTransferInspection:
    if not terms_acknowledged:
        # Fail before inspecting the first item so bundle metadata inspection cannot make
        # partial provider requests without the caller's explicit terms acknowledgement.
        _require_terms_acknowledgement(
            plan.items[0].terms,
            terms_acknowledged=False,
        )
    inspected = tuple(
        UCSCBundleTransferInspectionItem(
            role=item.role,
            metadata=inspect_resource(item.url),
        )
        for item in plan.items
    )
    return UCSCBundleTransferInspection(
        source_db=plan.source_db,
        target_db=plan.target_db,
        evidence_tier=plan.evidence_tier,
        items=inspected,
    )


def _cached_bundle_candidates(
    root: Path,
    source_db: str,
    target_db: str,
) -> dict[
    tuple[EvidenceAvailabilityTier, UCSCBundleResourceRole],
    list[tuple[Path, str]],
]:
    index_root = root / "by-url"
    candidates: dict[
        tuple[EvidenceAvailabilityTier, UCSCBundleResourceRole],
        list[tuple[Path, str]],
    ] = {}
    if not index_root.is_dir():
        return candidates

    for index_path in sorted(index_root.glob("*.json")):
        source_url = _cache_index_source_url(index_path)
        if source_url is None:
            continue
        binding = _cached_bundle_binding_for_url(
            source_url,
            source_db=source_db,
            target_db=target_db,
        )
        if binding is None:
            continue
        candidates.setdefault(binding, []).append((index_path, source_url))
    return candidates


def resolve_cached_ucsc_chain_resource_metadata(
    cache_root: ResourcePath,
    source_db: str,
    target_db: str,
    *,
    evidence_tier: EvidenceAvailabilityTier,
) -> CachedUCSCChainResource | None:
    """Resolve one exact-class directional chain structurally without hashing."""

    root = Path(cache_root)
    candidates = _cached_bundle_candidates(root, source_db, target_db)
    resource = _first_structural_cached_candidate(
        root,
        candidates.get((evidence_tier, UCSCBundleResourceRole.CHAIN), []),
    )
    if resource is None:
        return None
    return CachedUCSCChainResource(
        source_db=source_db,
        target_db=target_db,
        evidence_tier=evidence_tier,
        chain=resource,
    )


def load_cached_ucsc_chain_resource(
    cache_root: ResourcePath,
    source_db: str,
    target_db: str,
    *,
    evidence_tier: EvidenceAvailabilityTier,
    trusted_artifact_sha256_identifiers: frozenset[str] = frozenset(),
) -> CachedUCSCChainResource | None:
    """Load and verify only one exact-class directional chain from the cache.

    Net and reciprocal-best artifacts are intentionally not read. A trusted SHA-256
    identity may be supplied only by package code after a validated derived index has
    already established the exact source-chain identity.
    """

    trusted_sha256_hexes = frozenset(
        sha256_hex_from_identifier(identifier)
        for identifier in trusted_artifact_sha256_identifiers
    )
    root = Path(cache_root)
    candidates = _cached_bundle_candidates(root, source_db, target_db)
    resource = _first_verified_cached_candidate(
        root,
        candidates.get((evidence_tier, UCSCBundleResourceRole.CHAIN), []),
        trusted_sha256_hexes=trusted_sha256_hexes,
    )
    if resource is None:
        return None
    return CachedUCSCChainResource(
        source_db=source_db,
        target_db=target_db,
        evidence_tier=evidence_tier,
        chain=resource,
    )


def resolve_cached_ucsc_resource_bundle_metadata(
    cache_root: ResourcePath,
    source_db: str,
    target_db: str,
    *,
    evidence_tier: EvidenceAvailabilityTier | None = None,
) -> CachedUCSCResourceBundle | None:
    """Resolve the best structurally complete cached bundle without hashing bytes.

    This package-internal helper exists only so callers can look for a derived chain
    index before deciding whether rereading the original chain is necessary. Cache
    index metadata, exact artifact presence, and exact artifact size are checked, but
    this function makes no integrity claim about artifact contents. ``evidence_tier``
    selects one exact publication class when supplied; otherwise COMPARATIVE remains
    preferred. Normal callers should use :func:`load_cached_ucsc_resource_bundle`.
    """

    root = Path(cache_root)
    candidates = _cached_bundle_candidates(root, source_db, target_db)
    if not candidates:
        return None

    tiers: tuple[EvidenceAvailabilityTier, ...]
    if evidence_tier is not None:
        tiers = (evidence_tier,)
    else:
        tiers = (
            EvidenceAvailabilityTier.COMPARATIVE,
            EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )
    for selected_tier in tiers:
        bundle = _load_structural_cached_bundle_tier(
            root,
            candidates,
            source_db=source_db,
            target_db=target_db,
            evidence_tier=selected_tier,
        )
        if bundle is not None:
            return bundle
    return None


def load_cached_ucsc_resource_bundle(
    cache_root: ResourcePath,
    source_db: str,
    target_db: str,
    *,
    evidence_tier: EvidenceAvailabilityTier | None = None,
    progress_callback: CacheVerificationProgressCallback | None = None,
) -> CachedUCSCResourceBundle | None:
    """Load the best complete fully verified bundle already present in a local cache.

    This public loader preserves the original integrity contract: every selected
    provider artifact is checked for exact size and SHA-256 before the bundle returns.
    It performs no network access. When ``evidence_tier`` is omitted, it prefers a
    complete ``COMPARATIVE`` bundle over a ``LIFTOVER_ONLY`` chain; when supplied, it
    returns only that exact publication class.
    """

    return load_cached_ucsc_resource_bundle_for_indexed_assessment(
        cache_root,
        source_db,
        target_db,
        evidence_tier=evidence_tier,
        progress_callback=progress_callback,
    )


def load_cached_ucsc_resource_bundle_for_indexed_assessment(
    cache_root: ResourcePath,
    source_db: str,
    target_db: str,
    *,
    evidence_tier: EvidenceAvailabilityTier | None = None,
    progress_callback: CacheVerificationProgressCallback | None = None,
    trusted_artifact_sha256_identifiers: frozenset[str] = frozenset(),
) -> CachedUCSCResourceBundle | None:
    """Internal loader supporting exact identities prevalidated by derived artifacts.

    ``trusted_artifact_sha256_identifiers`` is used only after another package boundary
    has independently established the exact content identity. The CLI currently uses
    it for an original chain whose reusable chain index has passed source-identity and
    lookup-catalog validation. Matching cache metadata and exact file size are still
    required;
    all other selected provider artifacts retain normal SHA-256 verification.

    Trusted artifacts are excluded from ``progress_callback`` byte totals because their
    original bytes are not reread. This private hook must not be exposed as a general
    caller-controlled bypass of the public cache-integrity contract.
    """

    trusted_sha256_hexes = frozenset(
        sha256_hex_from_identifier(identifier)
        for identifier in trusted_artifact_sha256_identifiers
    )
    root = Path(cache_root)
    candidates = _cached_bundle_candidates(root, source_db, target_db)
    if not candidates:
        return None

    verification = _CacheVerificationTracker(progress_callback)
    tiers: tuple[EvidenceAvailabilityTier, ...]
    if evidence_tier is not None:
        tiers = (evidence_tier,)
    else:
        tiers = (
            EvidenceAvailabilityTier.COMPARATIVE,
            EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )
    for selected_tier in tiers:
        plan = _cache_verification_plan(
            root,
            candidates,
            selected_tier,
            trusted_sha256_hexes=trusted_sha256_hexes,
        )
        tier_verification = verification if plan else None
        if tier_verification is not None:
            tier_verification.start_attempt(plan)

        acquired: dict[UCSCBundleResourceRole, CachedResource] = {}
        for role in _bundle_roles_for_tier(selected_tier):
            resource = _first_verified_cached_candidate(
                root,
                candidates.get((selected_tier, role), []),
                verification=tier_verification,
                trusted_sha256_hexes=trusted_sha256_hexes,
            )
            if resource is None:
                break
            acquired[role] = resource
        else:
            if selected_tier is EvidenceAvailabilityTier.COMPARATIVE:
                bundle = CachedUCSCResourceBundle(
                    source_db=source_db,
                    target_db=target_db,
                    evidence_tier=selected_tier,
                    chain=acquired[UCSCBundleResourceRole.CHAIN],
                    net=acquired[UCSCBundleResourceRole.NET],
                    syntenic_net=acquired[UCSCBundleResourceRole.SYNTENIC_NET],
                    reciprocal_best_chain=acquired[
                        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN
                    ],
                    reciprocal_best_net=acquired[
                        UCSCBundleResourceRole.RECIPROCAL_BEST_NET
                    ],
                )
            else:
                bundle = CachedUCSCResourceBundle(
                    source_db=source_db,
                    target_db=target_db,
                    evidence_tier=selected_tier,
                    chain=acquired[UCSCBundleResourceRole.CHAIN],
                )
            verification.complete()
            return bundle
    return None


def _load_structural_cached_bundle_tier(
    root: Path,
    candidates: dict[
        tuple[EvidenceAvailabilityTier, UCSCBundleResourceRole],
        list[tuple[Path, str]],
    ],
    *,
    source_db: str,
    target_db: str,
    evidence_tier: EvidenceAvailabilityTier,
) -> CachedUCSCResourceBundle | None:
    resources: dict[UCSCBundleResourceRole, CachedResource] = {}
    for role in _bundle_roles_for_tier(evidence_tier):
        resource = _first_structural_cached_candidate(
            root,
            candidates.get((evidence_tier, role), []),
        )
        if resource is None:
            return None
        resources[role] = resource

    if evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return CachedUCSCResourceBundle(
            source_db=source_db,
            target_db=target_db,
            evidence_tier=evidence_tier,
            chain=resources[UCSCBundleResourceRole.CHAIN],
        )

    return CachedUCSCResourceBundle(
        source_db=source_db,
        target_db=target_db,
        evidence_tier=evidence_tier,
        chain=resources[UCSCBundleResourceRole.CHAIN],
        net=resources[UCSCBundleResourceRole.NET],
        syntenic_net=resources[UCSCBundleResourceRole.SYNTENIC_NET],
        reciprocal_best_chain=resources[UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN],
        reciprocal_best_net=resources[UCSCBundleResourceRole.RECIPROCAL_BEST_NET],
    )


def _first_structural_cached_candidate(
    root: Path,
    candidates: Sequence[tuple[Path, str]],
) -> CachedResource | None:
    for index_path, source_url in candidates:
        try:
            terms = ucsc_resource_terms(source_url)
        except ValueError:
            continue
        entry = _load_cached_resource_index_entry(
            root,
            index_path,
            source_url=source_url,
        )
        if entry is None:
            continue
        return CachedResource(
            path=entry.artifact_path,
            source_url=source_url,
            retrieved_at=entry.retrieved_at,
            sha256=f"sha256:{entry.sha256_hex}",
            size_bytes=entry.size_bytes,
            provider_checksum=entry.provider_checksum,
            terms=terms,
            cache_hit=True,
        )
    return None


class _CacheVerificationTracker:
    """Aggregate exact checksum-read work without declaring integrity too early."""

    def __init__(
        self,
        callback: CacheVerificationProgressCallback | None,
    ) -> None:
        self._callback = callback
        self._planned_paths: set[Path] = set()
        self._total_bytes = 0
        self._bytes_hashed = 0
        self._candidate_base = 0

    def start_attempt(self, plan: Sequence[tuple[Path, int]]) -> None:
        self._planned_paths = {index_path for index_path, _ in plan}
        self._total_bytes = sum(size_bytes for _, size_bytes in plan)
        self._bytes_hashed = 0
        self._candidate_base = 0
        self._emit(complete=False)

    def checksum_callback(
        self,
        index_path: Path,
        *,
        size_bytes: int,
    ) -> Callable[[int, int], None]:
        if index_path not in self._planned_paths:
            # A structurally valid first candidate can still fail its digest.  If the
            # loader then checks an alternate URL for the same role, include that extra
            # real hashing work rather than pretending the original denominator still
            # describes what is being verified.
            self._planned_paths.add(index_path)
            self._total_bytes += size_bytes
            self._emit(complete=False)

        self._candidate_base = self._bytes_hashed

        def update(bytes_hashed: int, total_bytes: int) -> None:
            if total_bytes != size_bytes:
                raise ValueError(
                    "cache verification progress total does not match indexed artifact size"
                )
            bounded = min(max(bytes_hashed, 0), total_bytes)
            if self._callback is not None:
                self._callback(
                    self._candidate_base + bounded,
                    self._total_bytes,
                    False,
                )

        return update

    def candidate_hashed(self, size_bytes: int) -> None:
        self._bytes_hashed = self._candidate_base + size_bytes

    def complete(self) -> None:
        if self._callback is not None and self._total_bytes >= 0:
            self._callback(self._total_bytes, self._total_bytes, True)

    def _emit(self, *, complete: bool) -> None:
        if self._callback is not None:
            self._callback(self._bytes_hashed, self._total_bytes, complete)


def _cache_verification_plan(
    root: Path,
    candidates: dict[
        tuple[EvidenceAvailabilityTier, UCSCBundleResourceRole],
        list[tuple[Path, str]],
    ],
    evidence_tier: EvidenceAvailabilityTier,
    *,
    trusted_sha256_hexes: frozenset[str] = frozenset(),
) -> tuple[tuple[Path, int], ...]:
    """Return the structurally viable artifact prefix the loader will hash.

    Planning reads cache metadata and checks exact artifact presence/size but deliberately
    does not hash bytes. The returned prefix stops at the first required role with no
    structurally viable candidate, matching the loader's own fail-closed role order. This
    lets progress describe real SHA-256 work even for an incomplete bundle without
    inventing bytes for roles that will never be hashed. The normal loader path remains
    authoritative for SHA-256 verification and can still try an alternate candidate if a
    planned artifact fails its digest.
    """

    plan: list[tuple[Path, int]] = []
    for role in _bundle_roles_for_tier(evidence_tier):
        selected: tuple[Path, int] | None = None
        for index_path, source_url in candidates.get((evidence_tier, role), []):
            try:
                ucsc_resource_terms(source_url)
            except ValueError:
                continue
            entry = _load_cached_resource_index_entry(
                root,
                index_path,
                source_url=source_url,
            )
            if entry is not None:
                selected = (
                    index_path,
                    0 if entry.sha256_hex in trusted_sha256_hexes else entry.size_bytes,
                )
                break
        if selected is None:
            break
        plan.append(selected)
    return tuple(plan)


def _cache_index_source_url(index_path: Path) -> str | None:
    try:
        payload = json.loads(index_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    source_url = payload.get("source_url")
    return source_url if isinstance(source_url, str) and source_url else None


def _cached_bundle_binding_for_url(
    url: str,
    *,
    source_db: str,
    target_db: str,
) -> tuple[EvidenceAvailabilityTier, UCSCBundleResourceRole] | None:
    for role in _bundle_roles_for_tier(EvidenceAvailabilityTier.COMPARATIVE):
        try:
            _validate_bundle_resource_binding(
                role,
                url,
                source_db=source_db,
                target_db=target_db,
                evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            )
        except ValueError:
            continue
        return EvidenceAvailabilityTier.COMPARATIVE, role

    try:
        _validate_bundle_resource_binding(
            UCSCBundleResourceRole.CHAIN,
            url,
            source_db=source_db,
            target_db=target_db,
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )
    except ValueError:
        return None
    return EvidenceAvailabilityTier.LIFTOVER_ONLY, UCSCBundleResourceRole.CHAIN


def _first_verified_cached_candidate(
    root: Path,
    candidates: Sequence[tuple[Path, str]],
    *,
    verification: _CacheVerificationTracker | None = None,
    trusted_sha256_hexes: frozenset[str] = frozenset(),
) -> CachedResource | None:
    for index_path, source_url in candidates:
        try:
            terms = ucsc_resource_terms(source_url)
        except ValueError:
            continue
        resource = _read_verified_cache_entry(
            root,
            index_path,
            source_url=source_url,
            terms=terms,
            verification=verification,
            trusted_sha256_hexes=trusted_sha256_hexes,
        )
        if resource is not None:
            return resource
    return None


def acquire_ucsc_resource_bundle(
    plan: UCSCBundleAcquisitionPlan,
    cache_root: ResourcePath,
    *,
    transfer_plan_acknowledged: bool,
    terms_acknowledged: bool,
    refresh: bool = False,
    progress_callback: UCSCBundleTransferProgressCallback | None = None,
) -> CachedUCSCResourceBundle:
    """Acquire every resource in an explicitly acknowledged bundle plan.

    This operation is complete-or-error at the returned-object boundary: a
    ``CachedUCSCResourceBundle`` is created only after every planned resource has been
    acquired or verified in cache.  Successfully published content-addressed artifacts
    from an earlier item are intentionally retained if a later item fails; they are
    valid immutable cache entries and can be reused on retry.

    ``transfer_plan_acknowledged`` is separate from provider-terms acknowledgement.
    It exists so a caller cannot pass a discovered COMPARATIVE bundle directly into a
    function that silently starts several transfers, including potentially very large
    UCSC resources.  The single-resource layer can now resume checksum-verifiable HTTPS
    transfers when the provider advertises the required strong validator/range metadata.

    ``progress_callback`` reports resource bytes complete toward the exact representation
    size when known.  A resumable resource may therefore begin above zero when a
    validator-bound partial is retained.  The final boolean distinguishes a verified
    cache hit from body-transfer progress so callers do not mislabel reused bytes as
    downloaded during the current run.
    """

    return _acquire_ucsc_resource_bundle(
        plan,
        cache_root,
        transfer_plan_acknowledged=transfer_plan_acknowledged,
        terms_acknowledged=terms_acknowledged,
        refresh=refresh,
        progress_callback=progress_callback,
        acquire_resource=acquire_ucsc_resource,
    )


def _acquire_ucsc_resource_bundle(
    plan: UCSCBundleAcquisitionPlan,
    cache_root: ResourcePath,
    *,
    transfer_plan_acknowledged: bool,
    terms_acknowledged: bool,
    refresh: bool,
    progress_callback: UCSCBundleTransferProgressCallback | None = None,
    acquire_resource: _ResourceAcquirer,
) -> CachedUCSCResourceBundle:
    if not transfer_plan_acknowledged:
        raise UCSCBundleAcquisitionPlanAcknowledgementRequired(
            "UCSC bundle acquisition requires explicit acknowledgement of the "
            "pre-transfer plan before any resource acquisition is attempted"
        )

    acquired: dict[UCSCBundleResourceRole, CachedResource] = {}
    for item in plan.items:
        resource_progress = (
            None
            if progress_callback is None
            else _resource_transfer_progress_for_role(progress_callback, item.role)
        )
        resource = acquire_resource(
            item.url,
            cache_root,
            terms_acknowledged=terms_acknowledged,
            refresh=refresh,
            progress_callback=resource_progress,
        )
        acquired[item.role] = resource
        if progress_callback is not None:
            # Body-transfer callbacks deliberately do not fire for a cache hit.  Emit
            # a distinct terminal event so the CLI can say "cached" rather than
            # misleadingly presenting those bytes as downloaded this run.  A final
            # event also gives unknown-length fresh transfers an exact denominator
            # once the complete resource has been acquired and verified.
            progress_callback(
                item.role,
                resource.size_bytes,
                resource.size_bytes,
                resource.cache_hit,
            )

    chain = acquired[UCSCBundleResourceRole.CHAIN]
    if plan.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return CachedUCSCResourceBundle(
            source_db=plan.source_db,
            target_db=plan.target_db,
            evidence_tier=plan.evidence_tier,
            chain=chain,
        )

    return CachedUCSCResourceBundle(
        source_db=plan.source_db,
        target_db=plan.target_db,
        evidence_tier=plan.evidence_tier,
        chain=chain,
        net=acquired[UCSCBundleResourceRole.NET],
        syntenic_net=acquired[UCSCBundleResourceRole.SYNTENIC_NET],
        reciprocal_best_chain=acquired[UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN],
        reciprocal_best_net=acquired[UCSCBundleResourceRole.RECIPROCAL_BEST_NET],
    )


def _resource_transfer_progress_for_role(
    progress_callback: UCSCBundleTransferProgressCallback,
    role: UCSCBundleResourceRole,
) -> ResourceTransferProgressCallback:
    """Bind one single-resource transfer stream to its bundle resource role."""

    def update(bytes_complete: int, total_bytes: int | None) -> None:
        progress_callback(role, bytes_complete, total_bytes, False)

    return update


def _bundle_roles_for_tier(
    tier: EvidenceAvailabilityTier,
) -> tuple[UCSCBundleResourceRole, ...]:
    if tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return (UCSCBundleResourceRole.CHAIN,)
    return (
        UCSCBundleResourceRole.CHAIN,
        UCSCBundleResourceRole.NET,
        UCSCBundleResourceRole.SYNTENIC_NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
        UCSCBundleResourceRole.RECIPROCAL_BEST_NET,
    )


def _validate_resource_role_filename(role: UCSCBundleResourceRole, url: str) -> None:
    """Reject a resource whose filename type does not match its declared role.

    This is intentionally a local filename-shape check.  Directional assembly binding
    belongs to the containing plan/cached-bundle validation, which has source/target
    context.
    """

    filename = PurePosixPath(_validate_ucsc_resource_url(url).path).name
    if role is UCSCBundleResourceRole.CHAIN:
        matches = filename.endswith((".all.chain.gz", ".over.chain.gz"))
    elif role is UCSCBundleResourceRole.NET:
        matches = filename.endswith(".net.gz") and not filename.endswith(
            (".syn.net.gz", ".rbest.net.gz")
        )
    elif role is UCSCBundleResourceRole.SYNTENIC_NET:
        matches = filename.endswith(".syn.net.gz")
    elif role is UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN:
        matches = filename.endswith(".rbest.chain.gz")
    else:
        matches = filename.endswith(".rbest.net.gz")

    if not matches:
        raise ValueError(
            f"bundle resource role {role.value} does not match filename {filename!r}"
        )


def _validate_bundle_resource_binding(
    role: UCSCBundleResourceRole,
    url: str,
    *,
    source_db: str,
    target_db: str,
    evidence_tier: EvidenceAvailabilityTier,
) -> None:
    """Require the exact directional UCSC filename for one bundle role.

    Hosting directory is deliberately not part of this identity check because UCSC
    can publish directional reciprocal-best files under the sibling comparison tree.
    The requested source/target pair and filename are authoritative for bundle role
    binding; discovery separately verifies where that exact file is hosted.
    """

    _validate_resource_role_filename(role, url)
    filename = PurePosixPath(_validate_ucsc_resource_url(url).path).name
    expected = _expected_bundle_resource_filename(
        role,
        source_db=source_db,
        target_db=target_db,
        evidence_tier=evidence_tier,
    )
    if filename != expected:
        raise ValueError(
            f"bundle resource {role.value} must use directional filename "
            f"{expected!r} for {source_db}->{target_db}, got {filename!r}"
        )


def _expected_bundle_resource_filename(
    role: UCSCBundleResourceRole,
    *,
    source_db: str,
    target_db: str,
    evidence_tier: EvidenceAvailabilityTier,
) -> str:
    if not source_db or not target_db:
        raise ValueError("bundle source_db and target_db must not be empty")

    if evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        if role is not UCSCBundleResourceRole.CHAIN:
            raise ValueError("LIFTOVER_ONLY bundle supports only the CHAIN role")
        target_title = target_db[0].upper() + target_db[1:]
        return f"{source_db}To{target_title}.over.chain.gz"

    suffix_by_role = {
        UCSCBundleResourceRole.CHAIN: "all.chain.gz",
        UCSCBundleResourceRole.NET: "net.gz",
        UCSCBundleResourceRole.SYNTENIC_NET: "syn.net.gz",
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN: "rbest.chain.gz",
        UCSCBundleResourceRole.RECIPROCAL_BEST_NET: "rbest.net.gz",
    }
    return f"{source_db}.{target_db}.{suffix_by_role[role]}"


def load_cached_ucsc_resource(
    cache_root: ResourcePath,
    url: str,
) -> CachedResource | None:
    """Load and SHA-256 verify one exact UCSC URL from the local cache only.

    This performs no provider access. It is useful for small assembly-metadata tables
    whose exact URLs are already known from an earlier verified discovery/acquisition
    step, while preserving the same content-addressed cache integrity contract used by
    evidence resources.
    """

    terms = ucsc_resource_terms(url)
    root = Path(cache_root)
    return _read_verified_cache_entry(
        root,
        _url_index_path(root, url),
        source_url=url,
        terms=terms,
    )


def acquire_ucsc_resource(
    url: str,
    cache_root: ResourcePath,
    *,
    terms_acknowledged: bool,
    refresh: bool = False,
    progress_callback: ResourceTransferProgressCallback | None = None,
) -> CachedResource:
    """Acquire one UCSC resource into a caller-selected content-addressed cache.

    Evidence-resource callers must explicitly acknowledge the applicable UCSC and
    directory-specific terms. For restricted liftOver chains this acknowledgement is
    especially important because UCSC currently states that downloading or using the
    files indicates EULA acceptance. The assembly ``chromInfo``/``chromAlias`` table
    dumps are exempt because UCSC's database download directory explicitly states that
    its files and tables are freely usable for any purpose.

    UCSC ``md5sum.txt`` metadata is verified when the resource's parent directory
    publishes an exact filename entry.  MD5 remains transfer-integrity metadata only;
    the cached artifact is named and reported by SHA-256. When no exact provider
    checksum was published at retrieval time, that metadata is retained with the URL
    index. A verified cache hit is intentionally usable offline and does not claim that
    the remote URL is unchanged; callers set ``refresh=True`` to contact UCSC and
    reacquire current bytes.  When UCSC also publishes an exact provider checksum and
    HEAD supplies a strong ETag, exact identity-encoded size, and byte-range support, an
    interrupted transfer may retain a validator-bound partial and resume it safely.

    ``progress_callback`` is body-transfer progress, expressed as exact resource bytes
    complete toward the known representation size (or ``None`` when the body length is
    unavailable).  On resume, its first event can be nonzero because a retained prefix
    is already validated against the same URL/size/strong-ETag representation.  Verified
    cache hits do not emit body-transfer events.  Reaching the reported byte total means
    the response body is complete; checksum/final publication can still be in progress.
    """

    return _acquire_ucsc_resource(
        url,
        cache_root,
        terms_acknowledged=terms_acknowledged,
        refresh=refresh,
        progress_callback=progress_callback,
        open_url=_open_url,
        now=lambda: datetime.now(UTC),
    )


def _acquire_ucsc_resource(
    url: str,
    cache_root: ResourcePath,
    *,
    terms_acknowledged: bool,
    refresh: bool = False,
    progress_callback: ResourceTransferProgressCallback | None = None,
    open_url: URLopener,
    now: Clock,
) -> CachedResource:
    terms = ucsc_resource_terms(url)
    _require_terms_acknowledgement(terms, terms_acknowledged=terms_acknowledged)

    root = Path(cache_root)
    root.mkdir(parents=True, exist_ok=True)

    index_path = _url_index_path(root, url)
    if not refresh:
        cached = _read_verified_cache_entry(
            root,
            index_path,
            source_url=url,
            terms=terms,
        )
        if cached is not None:
            return cached

    provider_checksum = _read_ucsc_provider_md5(url, open_url)
    remote_metadata = (
        _inspect_for_resumable_acquisition(url, open_url)
        if provider_checksum is not None
        else None
    )
    if (
        provider_checksum is not None
        and remote_metadata is not None
        and _supports_resumable_https(remote_metadata)
    ):
        return _download_and_cache_resumable(
            url,
            root,
            index_path=index_path,
            provider_checksum=provider_checksum,
            terms=terms,
            remote_metadata=remote_metadata,
            progress_callback=progress_callback,
            open_url=open_url,
            now=now,
        )

    return _download_and_cache(
        url,
        root,
        index_path=index_path,
        provider_checksum=provider_checksum,
        terms=terms,
        progress_callback=progress_callback,
        open_url=open_url,
        now=now,
    )


def _require_terms_acknowledgement(
    terms: UCSCResourceTerms,
    *,
    terms_acknowledged: bool,
) -> None:
    if terms.resource_class is UCSCResourceClass.ASSEMBLY_METADATA:
        return
    if terms_acknowledged:
        return

    restriction = (
        " UCSC identifies dedicated liftOver chain files as restricted and states "
        "that downloading/using them indicates EULA acceptance."
        if terms.restricted_liftover_chain
        else ""
    )
    raise UCSCResourceTermsAcknowledgementRequired(
        "UCSC provider access requires explicit acknowledgement that the applicable "
        "provider/directory terms were reviewed and permit the intended "
        f"use.{restriction} Terms: {terms.general_terms_url} and "
        f"{terms.directory_terms_url}"
    )


def _inspect_for_resumable_acquisition(
    url: str,
    open_url: URLopener,
) -> UCSCRemoteResourceMetadata | None:
    """Best-effort HEAD inspection used only to decide whether resume is available.

    Acquisition itself remains usable when HEAD metadata is unavailable: a failure here
    falls back to the existing fresh streaming path.  Provider terms have already been
    acknowledged before this helper is reached.
    """

    try:
        return _inspect_ucsc_resource(
            url,
            terms_acknowledged=True,
            open_url=open_url,
        )
    except UCSCResourceAcquisitionError:
        return None


def _supports_resumable_https(metadata: UCSCRemoteResourceMetadata) -> bool:
    """Return whether the advertised metadata is sufficient for safe byte-range resume.

    v1 resume deliberately requires a strong ETag, a known identity-encoded total size,
    and explicit byte-range support.  The acquisition caller separately requires an exact
    provider checksum before it selects this path, so a retained prefix cannot become a
    cached artifact on SHA-256 self-identity alone.  When any precondition is absent, the
    caller falls back to a fresh streaming transfer.
    """

    if metadata.content_length_bytes is None or metadata.content_length_bytes <= 0:
        return False
    if _identity_content_length_bytes(metadata) is None:
        return False
    if not _accepts_byte_ranges(metadata.accept_ranges):
        return False
    return _strong_etag(metadata.etag) is not None


def _accepts_byte_ranges(value: str | None) -> bool:
    if value is None:
        return False
    return any(token.strip().casefold() == "bytes" for token in value.split(","))


def _strong_etag(value: str | None) -> str | None:
    if value is None:
        return None
    etag = value.strip()
    if (
        len(etag) < 2
        or etag.casefold().startswith("w/")
        or not etag.startswith('"')
        or not etag.endswith('"')
    ):
        return None
    return etag


def _resource_get_request(
    url: str,
    *,
    range_start: int | None = None,
    if_range: str | None = None,
) -> Request:
    headers = {
        "User-Agent": _USER_AGENT,
        "Accept-Encoding": "identity",
    }
    if range_start is not None:
        headers["Range"] = f"bytes={range_start}-"
        if if_range is not None:
            headers["If-Range"] = if_range
    return Request(url, headers=headers)


def _resumable_partial_path(
    root: Path,
    url: str,
    *,
    total_size: int,
    etag: str,
) -> Path:
    """Return a deterministic partial path bound to URL + strong representation ETag.

    The URL hash keeps unrelated resources separate.  The ETag hash plus exact total
    length prevents bytes from a previous representation from being appended after the
    provider changes the object.  No sidecar metadata is needed to make that binding.
    """

    url_key = hashlib.sha256(url.encode("utf-8")).hexdigest()
    etag_key = hashlib.sha256(etag.encode("utf-8")).hexdigest()
    return root / "partials" / url_key[:2] / url_key / f"{etag_key}-{total_size}.part"


def _download_and_cache_resumable(
    url: str,
    root: Path,
    *,
    index_path: Path,
    provider_checksum: ProviderChecksum,
    terms: UCSCResourceTerms,
    remote_metadata: UCSCRemoteResourceMetadata,
    progress_callback: ResourceTransferProgressCallback | None,
    open_url: URLopener,
    now: Clock,
) -> CachedResource:
    total_size = cast(int, remote_metadata.content_length_bytes)
    etag = cast(str, _strong_etag(remote_metadata.etag))
    partial_path = _resumable_partial_path(
        root,
        url,
        total_size=total_size,
        etag=etag,
    )
    partial_path.parent.mkdir(parents=True, exist_ok=True)
    partial_path.touch(exist_ok=True)

    try:
        size_on_disk = partial_path.stat().st_size
    except OSError as exc:
        raise UCSCResourceAcquisitionError(
            f"failed to inspect resumable UCSC partial download {partial_path}: {exc}"
        ) from exc

    if size_on_disk > total_size:
        # A file larger than the representation it claims to prefix is unusable.  It
        # cannot be repaired by appending more bytes, so reset it before any request.
        try:
            partial_path.write_bytes(b"")
        except OSError as exc:
            raise UCSCResourceAcquisitionError(
                f"failed to reset invalid UCSC partial download {partial_path}: {exc}"
            ) from exc
        size_on_disk = 0

    if progress_callback is not None:
        # A retained validator-bound prefix is already complete acquisition state for
        # this exact representation, so resumable progress may legitimately start
        # above zero.  Subsequent events are emitted only after newly received bytes
        # have been written to the partial.
        progress_callback(size_on_disk, total_size)

    snapshot_path: Path | None = None
    snapshot_sha256: str | None = None
    snapshot_md5: str | None = None

    try:
        with partial_path.open("r+b") as output:
            output.seek(size_on_disk)

            if size_on_disk < total_size:
                request = _resource_get_request(
                    url,
                    range_start=size_on_disk if size_on_disk > 0 else None,
                    if_range=etag if size_on_disk > 0 else None,
                )
                try:
                    with open_url(request) as response:
                        if size_on_disk > 0:
                            response_bytes = _validate_resume_response(
                                response,
                                url=url,
                                range_start=size_on_disk,
                                total_size=total_size,
                                expected_etag=etag,
                            )
                        else:
                            response_bytes = _validate_fresh_resumable_response(
                                response,
                                url=url,
                                total_size=total_size,
                                expected_etag=etag,
                            )

                        bytes_before_response = size_on_disk
                        received = _stream_expected_response_bytes(
                            response,
                            output,
                            expected_bytes=response_bytes,
                            progress_callback=(
                                None
                                if progress_callback is None
                                else lambda response_bytes_received: progress_callback(
                                    bytes_before_response + response_bytes_received,
                                    total_size,
                                )
                            ),
                        )
                        size_on_disk += received
                except _ResumeRestartRequired:
                    raise
                except (HTTPError, URLError, TimeoutError, OSError) as exc:
                    _sync_partial_file(output)
                    size_on_disk = output.tell()
                    if size_on_disk == 0 and partial_path.exists():
                        partial_path.unlink()
                    raise UCSCResourceAcquisitionError(
                        f"interrupted UCSC resource download {url}; verified partial "
                        f"state was retained for retry: {exc}"
                    ) from exc

            _sync_partial_file(output)
            if size_on_disk == total_size:
                (
                    snapshot_path,
                    snapshot_sha256,
                    snapshot_md5,
                ) = _snapshot_resumable_partial(
                    output,
                    root,
                    total_size=total_size,
                    url=url,
                )
    except _ResumeRestartRequired:
        # Never splice a response that fails the Range/validator contract onto the
        # existing prefix.  Discard this resumable state and restart through the fresh
        # path; the fresh downloader writes to its own unique temporary file.
        partial_path.unlink(missing_ok=True)
        return _download_and_cache(
            url,
            root,
            index_path=index_path,
            provider_checksum=provider_checksum,
            terms=terms,
            progress_callback=progress_callback,
            open_url=open_url,
            now=now,
        )

    if size_on_disk != total_size:
        raise UCSCResourceAcquisitionError(
            f"incomplete UCSC resource download {url}: expected {total_size} bytes, "
            f"received {size_on_disk}; partial state retained for retry"
        )

    # The deterministic partial path is intentionally shared resumable state.  Never
    # publish that inode directly: another process may still have it open for writing,
    # and an ``os.replace`` would let that stale descriptor mutate the supposedly
    # immutable content-addressed artifact after publication.  The snapshot above was
    # therefore copied through this process's already-open handle into a unique private
    # file and re-hashed independently before publication.
    if snapshot_path is None or snapshot_sha256 is None or snapshot_md5 is None:
        raise UCSCResourceAcquisitionError(
            f"completed UCSC resumable download {url} could not be finalized safely"
        )

    try:
        if snapshot_md5 != provider_checksum.value:
            partial_path.unlink(missing_ok=True)
            raise ResourceChecksumMismatchError(
                f"md5 checksum mismatch for downloaded UCSC resource {url}: "
                f"expected {provider_checksum.value}, got {snapshot_md5}"
            )

        result = _publish_completed_download(
            snapshot_path,
            url=url,
            root=root,
            index_path=index_path,
            sha256_hex=snapshot_sha256,
            size_bytes=size_on_disk,
            provider_checksum=provider_checksum,
            terms=terms,
            now=now,
        )
        partial_path.unlink(missing_ok=True)
        return result
    finally:
        snapshot_path.unlink(missing_ok=True)


def _snapshot_resumable_partial(
    source: BinaryIO,
    root: Path,
    *,
    total_size: int,
    url: str,
) -> tuple[Path, str, str]:
    """Copy shared resumable state into a private, independently verified snapshot."""

    temp_dir = root / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix="resume-snapshot-", suffix=".part", dir=temp_dir
    )
    snapshot_path = Path(temp_name)
    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)

    try:
        with os.fdopen(fd, "wb") as destination:
            source.seek(0)
            remaining = total_size
            while remaining:
                chunk = source.read(min(_CHUNK_SIZE, remaining))
                if not chunk:
                    raise UCSCResourceAcquisitionError(
                        f"resumable UCSC partial for {url} changed while being finalized"
                    )
                destination.write(chunk)
                sha256.update(chunk)
                md5.update(chunk)
                remaining -= len(chunk)

            if source.read(1):
                raise UCSCResourceAcquisitionError(
                    f"resumable UCSC partial for {url} exceeded its expected size "
                    "while being finalized"
                )

            destination.flush()
            os.fsync(destination.fileno())

        return snapshot_path, sha256.hexdigest(), md5.hexdigest()
    except Exception:
        snapshot_path.unlink(missing_ok=True)
        raise


def _validate_fresh_resumable_response(
    response: _BinaryResponse,
    *,
    url: str,
    total_size: int,
    expected_etag: str,
) -> int:
    if response.status != 200:
        raise _ResumeRestartRequired(
            f"fresh UCSC resumable GET returned HTTP {response.status}"
        )
    _require_identity_content_encoding(response, url=url)
    response_etag = _optional_response_header(response, "ETag")
    if response_etag != expected_etag:
        raise _ResumeRestartRequired(
            "fresh GET representation does not match HEAD ETag"
        )
    response_length = _response_content_length(response)
    if response_length is not None and response_length != total_size:
        raise _ResumeRestartRequired(
            "fresh GET length does not match HEAD Content-Length"
        )
    return total_size


def _validate_resume_response(
    response: _BinaryResponse,
    *,
    url: str,
    range_start: int,
    total_size: int,
    expected_etag: str,
) -> int:
    if response.status != 206:
        raise _ResumeRestartRequired(
            f"resume request returned HTTP {response.status} instead of 206"
        )
    _require_identity_content_encoding(response, url=url)
    response_etag = _optional_response_header(response, "ETag")
    if response_etag is not None and response_etag != expected_etag:
        raise _ResumeRestartRequired("resume response ETag changed")

    content_range = _optional_response_header(response, "Content-Range")
    if content_range is None:
        raise _ResumeRestartRequired("resume response omitted Content-Range")
    match = _CONTENT_RANGE_RE.fullmatch(content_range)
    if match is None:
        raise _ResumeRestartRequired("resume response returned malformed Content-Range")
    start, end, total = (int(value) for value in match.groups())
    if start != range_start or end < start or total != total_size or end >= total_size:
        raise _ResumeRestartRequired(
            "resume response Content-Range does not match request"
        )

    expected_bytes = end - start + 1
    response_length = _response_content_length(response)
    if response_length is not None and response_length != expected_bytes:
        raise _ResumeRestartRequired(
            "resume response Content-Length does not match Content-Range"
        )
    return expected_bytes


def _stream_expected_response_bytes(
    response: _BinaryResponse,
    output: BinaryIO,
    *,
    expected_bytes: int,
    progress_callback: Callable[[int], None] | None = None,
) -> int:
    remaining = expected_bytes
    received = 0
    while remaining:
        chunk = response.read(min(_CHUNK_SIZE, remaining))
        if not chunk:
            break
        output.write(chunk)
        received += len(chunk)
        remaining -= len(chunk)
        if progress_callback is not None:
            progress_callback(received)

    if remaining:
        return received
    if response.read(1):
        raise _ResumeRestartRequired("UCSC response exceeded its advertised byte range")
    return received


def _sync_partial_file(output: BinaryIO) -> None:
    output.flush()
    os.fsync(output.fileno())


def _require_identity_content_encoding(
    response: _BinaryResponse,
    *,
    url: str,
) -> None:
    encoding = _optional_response_header(response, "Content-Encoding")
    if encoding is not None and encoding.casefold() != "identity":
        raise UCSCResourceAcquisitionError(
            f"UCSC resource {url} returned unsupported Content-Encoding {encoding!r}; "
            "exact provider bytes require identity encoding"
        )


def _download_and_cache(
    url: str,
    root: Path,
    *,
    index_path: Path,
    provider_checksum: ProviderChecksum | None,
    terms: UCSCResourceTerms,
    progress_callback: ResourceTransferProgressCallback | None,
    open_url: URLopener,
    now: Clock,
) -> CachedResource:
    temp_dir = root / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix="download-", suffix=".part", dir=temp_dir)
    temp_path = Path(temp_name)

    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False) if provider_checksum is not None else None
    try:
        with os.fdopen(fd, "wb") as output:
            try:
                with open_url(_resource_get_request(url)) as response:
                    _require_identity_content_encoding(response, url=url)
                    expected_size = _response_content_length(response)
                    size_bytes = 0
                    if progress_callback is not None:
                        progress_callback(0, expected_size)
                    while chunk := response.read(_CHUNK_SIZE):
                        output.write(chunk)
                        size_bytes += len(chunk)
                        sha256.update(chunk)
                        if md5 is not None:
                            md5.update(chunk)
                        if progress_callback is not None:
                            progress_callback(size_bytes, expected_size)
                    if expected_size is not None and size_bytes != expected_size:
                        raise UCSCResourceAcquisitionError(
                            f"incomplete UCSC resource download {url}: expected "
                            f"{expected_size} bytes, received {size_bytes}"
                        )
            except (HTTPError, URLError, TimeoutError, OSError) as exc:
                raise UCSCResourceAcquisitionError(
                    f"failed to download UCSC resource {url}: {exc}"
                ) from exc
            output.flush()
            os.fsync(output.fileno())

        sha256_hex = sha256.hexdigest()
        if provider_checksum is not None and md5 is not None:
            actual_md5 = md5.hexdigest()
            if actual_md5 != provider_checksum.value:
                raise ResourceChecksumMismatchError(
                    f"md5 checksum mismatch for downloaded UCSC resource {url}: "
                    f"expected {provider_checksum.value}, got {actual_md5}"
                )

        return _publish_completed_download(
            temp_path,
            url=url,
            root=root,
            index_path=index_path,
            sha256_hex=sha256_hex,
            size_bytes=size_bytes,
            provider_checksum=provider_checksum,
            terms=terms,
            now=now,
        )
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _publish_completed_download(
    completed_path: Path,
    *,
    url: str,
    root: Path,
    index_path: Path,
    sha256_hex: str,
    size_bytes: int,
    provider_checksum: ProviderChecksum | None,
    terms: UCSCResourceTerms,
    now: Clock,
) -> CachedResource:
    artifact_path = _artifact_path(root, sha256_hex)
    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    if artifact_path.exists():
        existing_sha256 = compute_resource_checksum(
            artifact_path, ResourceChecksumAlgorithm.SHA256
        )
        if existing_sha256 == sha256_hex:
            completed_path.unlink(missing_ok=True)
        else:
            os.replace(completed_path, artifact_path)
    else:
        try:
            os.replace(completed_path, artifact_path)
        except FileNotFoundError:
            # Another process may have completed the same exact representation while this
            # process was streaming it.  Accept that race only if the final content-addressed
            # artifact now exists and verifies to the digest computed here.
            if not artifact_path.is_file():
                raise
            existing_sha256 = compute_resource_checksum(
                artifact_path, ResourceChecksumAlgorithm.SHA256
            )
            if existing_sha256 != sha256_hex:
                raise UCSCResourceAcquisitionError(
                    f"concurrent UCSC cache publication produced an unexpected artifact "
                    f"for {url}"
                )

    retrieved_at = _format_timestamp(now())
    _write_url_index(
        index_path,
        source_url=url,
        retrieved_at=retrieved_at,
        sha256=sha256_hex,
        size_bytes=size_bytes,
        provider_checksum=provider_checksum,
        terms=terms,
    )
    return CachedResource(
        path=artifact_path,
        source_url=url,
        retrieved_at=retrieved_at,
        sha256=f"sha256:{sha256_hex}",
        size_bytes=size_bytes,
        provider_checksum=provider_checksum,
        terms=terms,
        cache_hit=False,
    )


def _read_ucsc_provider_md5(
    resource_url: str,
    open_url: URLopener,
) -> ProviderChecksum | None:
    filename = PurePosixPath(urlsplit(resource_url).path).name
    checksum_url = urljoin(resource_url, "md5sum.txt")
    try:
        with open_url(
            Request(checksum_url, headers={"User-Agent": _USER_AGENT})
        ) as response:
            payload = response.read().decode("utf-8")
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise UCSCResourceAcquisitionError(
            f"failed to read UCSC checksum metadata {checksum_url}: HTTP {exc.code}"
        ) from exc
    except (URLError, TimeoutError, OSError, UnicodeError) as exc:
        raise UCSCResourceAcquisitionError(
            f"failed to read UCSC checksum metadata {checksum_url}: {exc}"
        ) from exc

    for line in payload.splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) != 2:
            continue
        checksum, listed_name = fields
        listed_name = listed_name.strip().removeprefix("*")
        if listed_name != filename:
            continue
        if _MD5_RE.fullmatch(checksum) is None:
            raise UCSCResourceAcquisitionError(
                f"UCSC checksum metadata for {filename} contains a malformed MD5"
            )
        return ProviderChecksum(
            algorithm=ResourceChecksumAlgorithm.MD5,
            value=checksum.lower(),
            source_url=checksum_url,
        )

    return None


def _load_cached_resource_index_entry(
    root: Path,
    index_path: Path,
    *,
    source_url: str,
) -> _CachedResourceIndexEntry | None:
    """Validate cache metadata and exact artifact size without hashing its bytes."""

    if not index_path.exists():
        return None

    try:
        payload = json.loads(index_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None

    if not isinstance(payload, dict):
        return None
    if payload.get("schema_version") != _CACHE_SCHEMA_VERSION:
        return None
    if payload.get("source_url") != source_url:
        return None
    sha256_hex = payload.get("sha256")
    if (
        not isinstance(sha256_hex, str)
        or re.fullmatch(r"[0-9a-f]{64}", sha256_hex) is None
    ):
        return None

    stored_provider = payload.get("provider_checksum")
    provider_checksum = _provider_checksum_from_index(stored_provider)
    if stored_provider is not None and provider_checksum is None:
        return None

    size_bytes = payload.get("size_bytes")
    if not isinstance(size_bytes, int) or size_bytes < 0:
        return None

    artifact_path = _artifact_path(root, sha256_hex)
    if not artifact_path.is_file() or artifact_path.stat().st_size != size_bytes:
        return None

    retrieved_at = payload.get("retrieved_at")
    if not isinstance(retrieved_at, str) or not retrieved_at:
        return None

    return _CachedResourceIndexEntry(
        artifact_path=artifact_path,
        retrieved_at=retrieved_at,
        sha256_hex=sha256_hex,
        size_bytes=size_bytes,
        provider_checksum=provider_checksum,
    )


def _read_verified_cache_entry(
    root: Path,
    index_path: Path,
    *,
    source_url: str,
    terms: UCSCResourceTerms,
    verification: _CacheVerificationTracker | None = None,
    trusted_sha256_hexes: frozenset[str] = frozenset(),
) -> CachedResource | None:
    entry = _load_cached_resource_index_entry(
        root,
        index_path,
        source_url=source_url,
    )
    if entry is None:
        return None

    if entry.sha256_hex in trusted_sha256_hexes:
        return CachedResource(
            path=entry.artifact_path,
            source_url=source_url,
            retrieved_at=entry.retrieved_at,
            sha256=f"sha256:{entry.sha256_hex}",
            size_bytes=entry.size_bytes,
            provider_checksum=entry.provider_checksum,
            terms=terms,
            cache_hit=True,
        )

    checksum_progress = (
        verification.checksum_callback(index_path, size_bytes=entry.size_bytes)
        if verification is not None
        else None
    )
    actual_sha256 = compute_resource_checksum(
        entry.artifact_path,
        ResourceChecksumAlgorithm.SHA256,
        progress_callback=checksum_progress,
    )
    if verification is not None:
        verification.candidate_hashed(entry.size_bytes)
    if actual_sha256 != entry.sha256_hex:
        return None

    return CachedResource(
        path=entry.artifact_path,
        source_url=source_url,
        retrieved_at=entry.retrieved_at,
        sha256=f"sha256:{entry.sha256_hex}",
        size_bytes=entry.size_bytes,
        provider_checksum=entry.provider_checksum,
        terms=terms,
        cache_hit=True,
    )


def _response_content_length(response: _BinaryResponse) -> int | None:
    value = response.getheader("Content-Length")
    if value is None:
        return None
    try:
        size = int(value)
    except ValueError as exc:
        raise UCSCResourceAcquisitionError(
            f"invalid UCSC Content-Length header: {value!r}"
        ) from exc
    if size < 0:
        raise UCSCResourceAcquisitionError(
            f"invalid UCSC Content-Length header: {value!r}"
        )
    return size


def _optional_response_header(response: _BinaryResponse, name: str) -> str | None:
    value = response.getheader(name)
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None


def _provider_checksum_from_index(value: object) -> ProviderChecksum | None:
    if value is None:
        return None
    if not isinstance(value, dict):
        return None
    if value.get("algorithm") != ResourceChecksumAlgorithm.MD5.value:
        return None
    checksum = value.get("value")
    source_url = value.get("source_url")
    if not isinstance(checksum, str) or _MD5_RE.fullmatch(checksum) is None:
        return None
    if not isinstance(source_url, str) or not source_url:
        return None
    return ProviderChecksum(
        algorithm=ResourceChecksumAlgorithm.MD5,
        value=checksum.lower(),
        source_url=source_url,
    )


def _write_url_index(
    index_path: Path,
    *,
    source_url: str,
    retrieved_at: str,
    sha256: str,
    size_bytes: int,
    provider_checksum: ProviderChecksum | None,
    terms: UCSCResourceTerms,
) -> None:
    index_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": _CACHE_SCHEMA_VERSION,
        "source_url": source_url,
        "retrieved_at": retrieved_at,
        "sha256": sha256,
        "size_bytes": size_bytes,
        "provider_checksum": (
            None
            if provider_checksum is None
            else {
                "algorithm": provider_checksum.algorithm.value,
                "value": provider_checksum.value,
                "source_url": provider_checksum.source_url,
            }
        ),
        "terms": {
            "resource_class": terms.resource_class.value,
            "general_terms_url": terms.general_terms_url,
            "directory_terms_url": terms.directory_terms_url,
            "restricted_liftover_chain": terms.restricted_liftover_chain,
        },
    }
    fd, temp_name = tempfile.mkstemp(
        prefix=f"{index_path.stem}-", suffix=".tmp", dir=index_path.parent
    )
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, sort_keys=True, indent=2) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, index_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _artifact_path(root: Path, sha256_hex: str) -> Path:
    return root / "artifacts" / "sha256" / sha256_hex[:2] / sha256_hex


def _url_index_path(root: Path, url: str) -> Path:
    key = hashlib.sha256(url.encode("utf-8")).hexdigest()
    return root / "by-url" / f"{key}.json"


def _directory_url(
    parts: SplitResult, segments: list[str], directory_index: int
) -> str:
    path = "/" + "/".join(segments[: directory_index + 1]) + "/"
    return urlunsplit((parts.scheme, parts.netloc, path, "", ""))


def _validate_ucsc_resource_url(url: str) -> SplitResult:
    parts = urlsplit(url)
    if parts.scheme != "https" or parts.hostname not in _UCSC_HOSTS:
        raise ValueError(
            "UCSC resource URL must use https on an approved hgdownload host"
        )
    if parts.query or parts.fragment:
        raise ValueError("UCSC resource URL must not contain a query or fragment")
    if not parts.path.startswith("/goldenPath/") or parts.path.endswith("/"):
        raise ValueError("UCSC resource URL must identify a file under /goldenPath/")
    return parts


def _format_timestamp(value: datetime) -> str:
    if value.tzinfo is None:
        raise ValueError("cache retrieval timestamp must be timezone-aware")
    return value.astimezone(UTC).isoformat().replace("+00:00", "Z")


def _open_url(request: Request) -> _BinaryResponse:
    return cast(_BinaryResponse, urlopen(request, timeout=30))
