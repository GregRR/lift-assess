"""Acquire UCSC resources into an explicit content-addressed local cache.

Discovery answers which provider resources exist; this module handles the separate
act of planning and retrieving already-discovered UCSC URLs.  Retrieval is deliberately
explicit about provider terms, streams bytes without materializing large resources,
verifies UCSC-published MD5 metadata when available, and stores the resulting exact
artifact by liftAssess's canonical SHA-256 identity.

The cache root is always supplied by the caller.  liftAssess does not create a cache
inside the source tree and this module does not choose a user/OS default on behalf of
the future CLI.

Primary provider references checked 2026-08-13:

- UCSC data/software licensing: https://genome.ucsc.edu/license/
- restricted canFam3 liftOver-chain terms: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/
- canFam3/canFam4 comparative terms and files: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/
- comparative MD5 metadata: https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/md5sum.txt
- reciprocal-best MD5 metadata: https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/md5sum.txt

The canFam3/canFam4 comparison demonstrates why checksum lookup is exact-filename
based: its MD5 file covers some comparison artifacts but not ``canFam3.canFam4.net.gz``.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime
from enum import Enum
from pathlib import Path, PurePosixPath
from types import TracebackType
from typing import Protocol, Self, TypeAlias, cast
from urllib.error import HTTPError, URLError
from urllib.parse import SplitResult, urljoin, urlsplit, urlunsplit
from urllib.request import Request, urlopen

from .models import EvidenceAvailabilityTier
from .resource_identity import (
    ResourceChecksumAlgorithm,
    ResourceChecksumMismatchError,
    compute_resource_checksum,
)
from .resources import UCSCResourceBundle

ResourcePath: TypeAlias = str | os.PathLike[str]

_UCSC_HOSTS = frozenset({"hgdownload.soe.ucsc.edu", "hgdownload.gi.ucsc.edu"})
_UCSC_LICENSE_URL = "https://genome.ucsc.edu/license/"
_USER_AGENT = "liftAssess/0.0 resource-acquisition"
_CHUNK_SIZE = 1024 * 1024
_MD5_RE = re.compile(r"^[0-9a-fA-F]{32}$")
_CACHE_SCHEMA_VERSION = 1


class UCSCResourceAcquisitionError(RuntimeError):
    """A UCSC resource could not be acquired or cached reliably."""


class UCSCResourceTermsAcknowledgementRequired(UCSCResourceAcquisitionError):
    """Retrieval was requested without explicit acknowledgement of provider terms."""


class UCSCBundleAcquisitionPlanAcknowledgementRequired(UCSCResourceAcquisitionError):
    """Bundle retrieval was requested without acknowledging its explicit transfer plan."""


class UCSCResourceClass(str, Enum):
    """UCSC publication classes currently produced by the resource resolver."""

    COMPARATIVE = "COMPARATIVE"
    LIFTOVER_CHAIN = "LIFTOVER_CHAIN"


class UCSCBundleResourceRole(str, Enum):
    """One resource role in a discovered UCSC evidence bundle."""

    CHAIN = "CHAIN"
    NET = "NET"
    SYNTENIC_NET = "SYNTENIC_NET"
    RECIPROCAL_BEST_CHAIN = "RECIPROCAL_BEST_CHAIN"
    RECIPROCAL_BEST_NET = "RECIPROCAL_BEST_NET"


@dataclass(frozen=True)
class UCSCResourceTerms:
    """Provider terms references that must be surfaced before retrieval.

    ``restricted_liftover_chain`` distinguishes UCSC's dedicated
    ``liftOver/*.over.chain.gz`` files, for which UCSC currently states that
    downloading/using indicates EULA acceptance and that free use is limited to the
    described non-commercial/nonprofit cases.  Comparative resources are a distinct
    publication class and retain their own directory terms URL instead.
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
    def read(self, size: int = -1) -> bytes:
        ...

    def getheader(self, name: str, default: str | None = None) -> str | None:
        ...

    def __enter__(self) -> Self:
        ...

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        ...


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
    ) -> CachedResource:
        ...


def ucsc_resource_terms(url: str) -> UCSCResourceTerms:
    """Return the terms references for one resolver-produced UCSC resource URL."""

    parts = _validate_ucsc_resource_url(url)
    segments = parts.path.strip("/").split("/")
    filename = PurePosixPath(parts.path).name

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

    raise ValueError("unsupported UCSC resource URL outside comparative/liftOver paths")


def plan_ucsc_bundle_acquisition(
    bundle: UCSCResourceBundle,
) -> UCSCBundleAcquisitionPlan:
    """Build an inspectable no-network transfer plan from a discovered bundle.

    Planning does not create the cache or contact UCSC.  It enumerates the exact URLs
    required by the bundle's evidence-availability tier and records the applicable
    provider terms for each resource.  Size discovery is intentionally not guessed in
    this slice; callers must inspect this plan and explicitly acknowledge it before
    execution, and future CLI work can add measured remote-size metadata separately.
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


def acquire_ucsc_resource_bundle(
    plan: UCSCBundleAcquisitionPlan,
    cache_root: ResourcePath,
    *,
    transfer_plan_acknowledged: bool,
    terms_acknowledged: bool,
    refresh: bool = False,
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
    UCSC resources.  Remote size metadata and resumable large-file transfer remain a
    later milestone.
    """

    return _acquire_ucsc_resource_bundle(
        plan,
        cache_root,
        transfer_plan_acknowledged=transfer_plan_acknowledged,
        terms_acknowledged=terms_acknowledged,
        refresh=refresh,
        acquire_resource=acquire_ucsc_resource,
    )


def _acquire_ucsc_resource_bundle(
    plan: UCSCBundleAcquisitionPlan,
    cache_root: ResourcePath,
    *,
    transfer_plan_acknowledged: bool,
    terms_acknowledged: bool,
    refresh: bool,
    acquire_resource: _ResourceAcquirer,
) -> CachedUCSCResourceBundle:
    if not transfer_plan_acknowledged:
        raise UCSCBundleAcquisitionPlanAcknowledgementRequired(
            "UCSC bundle acquisition requires explicit acknowledgement of the "
            "pre-transfer plan before any resource acquisition is attempted"
        )

    acquired: dict[UCSCBundleResourceRole, CachedResource] = {}
    for item in plan.items:
        acquired[item.role] = acquire_resource(
            item.url,
            cache_root,
            terms_acknowledged=terms_acknowledged,
            refresh=refresh,
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


def _validate_resource_role_filename(
    role: UCSCBundleResourceRole, url: str
) -> None:
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


def acquire_ucsc_resource(
    url: str,
    cache_root: ResourcePath,
    *,
    terms_acknowledged: bool,
    refresh: bool = False,
) -> CachedResource:
    """Acquire one UCSC resource into a caller-selected content-addressed cache.

    The caller must explicitly acknowledge that they reviewed the applicable UCSC and
    directory-specific terms.  For restricted liftOver chains this acknowledgement is
    especially important because UCSC currently states that downloading or using the
    files indicates EULA acceptance.

    UCSC ``md5sum.txt`` metadata is verified when the resource's parent directory
    publishes an exact filename entry.  MD5 remains transfer-integrity metadata only;
    the cached artifact is named and reported by SHA-256. When no exact provider
    checksum was published at retrieval time, that metadata is retained with the URL
    index. A verified cache hit is intentionally usable offline and does not claim that
    the remote URL is unchanged; callers set ``refresh=True`` to contact UCSC and
    reacquire current bytes.
    """

    return _acquire_ucsc_resource(
        url,
        cache_root,
        terms_acknowledged=terms_acknowledged,
        refresh=refresh,
        open_url=_open_url,
        now=lambda: datetime.now(UTC),
    )


def _acquire_ucsc_resource(
    url: str,
    cache_root: ResourcePath,
    *,
    terms_acknowledged: bool,
    refresh: bool = False,
    open_url: URLopener,
    now: Clock,
) -> CachedResource:
    terms = ucsc_resource_terms(url)
    if not terms_acknowledged:
        restriction = (
            " UCSC identifies dedicated liftOver chain files as restricted and states "
            "that downloading/using them indicates EULA acceptance."
            if terms.restricted_liftover_chain
            else ""
        )
        raise UCSCResourceTermsAcknowledgementRequired(
            "UCSC resource retrieval requires explicit acknowledgement that the "
            "applicable provider/directory terms were reviewed and permit the intended "
            f"use.{restriction} Terms: {terms.general_terms_url} and "
            f"{terms.directory_terms_url}"
        )

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
    return _download_and_cache(
        url,
        root,
        index_path=index_path,
        provider_checksum=provider_checksum,
        terms=terms,
        open_url=open_url,
        now=now,
    )


def _download_and_cache(
    url: str,
    root: Path,
    *,
    index_path: Path,
    provider_checksum: ProviderChecksum | None,
    terms: UCSCResourceTerms,
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
                with open_url(Request(url, headers={"User-Agent": _USER_AGENT})) as response:
                    expected_size = _response_content_length(response)
                    size_bytes = 0
                    while chunk := response.read(_CHUNK_SIZE):
                        output.write(chunk)
                        size_bytes += len(chunk)
                        sha256.update(chunk)
                        if md5 is not None:
                            md5.update(chunk)
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

        artifact_path = _artifact_path(root, sha256_hex)
        artifact_path.parent.mkdir(parents=True, exist_ok=True)
        if artifact_path.exists():
            existing_sha256 = compute_resource_checksum(
                artifact_path, ResourceChecksumAlgorithm.SHA256
            )
            if existing_sha256 == sha256_hex:
                temp_path.unlink()
            else:
                os.replace(temp_path, artifact_path)
        else:
            os.replace(temp_path, artifact_path)

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
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _read_ucsc_provider_md5(
    resource_url: str,
    open_url: URLopener,
) -> ProviderChecksum | None:
    filename = PurePosixPath(urlsplit(resource_url).path).name
    checksum_url = urljoin(resource_url, "md5sum.txt")
    try:
        with open_url(Request(checksum_url, headers={"User-Agent": _USER_AGENT})) as response:
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


def _read_verified_cache_entry(
    root: Path,
    index_path: Path,
    *,
    source_url: str,
    terms: UCSCResourceTerms,
) -> CachedResource | None:
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
    if not isinstance(sha256_hex, str) or re.fullmatch(r"[0-9a-f]{64}", sha256_hex) is None:
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
    actual_sha256 = compute_resource_checksum(
        artifact_path, ResourceChecksumAlgorithm.SHA256
    )
    if actual_sha256 != sha256_hex:
        return None

    retrieved_at = payload.get("retrieved_at")
    if not isinstance(retrieved_at, str) or not retrieved_at:
        return None

    return CachedResource(
        path=artifact_path,
        source_url=source_url,
        retrieved_at=retrieved_at,
        sha256=f"sha256:{sha256_hex}",
        size_bytes=size_bytes,
        provider_checksum=provider_checksum,
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
        raise ValueError("UCSC resource URL must use https on an approved hgdownload host")
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
