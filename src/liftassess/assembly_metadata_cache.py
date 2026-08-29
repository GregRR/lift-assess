"""Acquire and load authoritative UCSC assembly-sequence metadata.

The metadata path is deliberately separate from alignment evidence acquisition.
``chromInfo`` and ``chromAlias`` are small UCSC database tables used for input
preflight; they are not mapping evidence and are not part of an evidence tier.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from io import BytesIO
from pathlib import Path
from re import fullmatch
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urljoin
from urllib.request import Request, urlopen

from .assembly_metadata import (
    AssemblySequenceCatalog,
    attach_ncbi_sequence_role_context,
    build_ucsc_assembly_sequence_catalog,
    parse_ucsc_assembly_description_accession,
)
from .models import (
    AssemblyIdentifier,
    ProvenanceIdentifier,
    ProvenanceIdentifierKind,
    ProvenanceSource,
)
from .resource_cache import (
    CachedResource,
    ResourcePath,
    acquire_ucsc_resource,
    load_cached_ucsc_resource,
)
from .resource_identity import sha256_hex_from_identifier
from .resource_stream import open_text_resource
from .resources import UCSCAssemblyMetadataResources

_USER_AGENT = "liftAssess assembly-role metadata"
_NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
_ROLE_CACHE_SCHEMA_VERSION = 1

_UCSC_GOLDEN_PATH = "https://hgdownload.soe.ucsc.edu/goldenPath/"
_UCSC_DB_PATTERN = r"[A-Za-z0-9_.-]+"


class AssemblyRoleMetadataAcquisitionError(RuntimeError):
    """Target-role metadata could not be acquired from UCSC/NCBI."""


class AssemblyRoleMetadataCacheIntegrityError(RuntimeError):
    """Cached target-role metadata does not match its recorded content identity."""


@dataclass(frozen=True)
class CachedAssemblyRoleArtifact:
    """One exact small metadata artifact retained by SHA-256."""

    path: Path
    source_url: str
    retrieved_at: str
    sha256: str
    size_bytes: int
    cache_hit: bool
    archive_member: str | None = None


@dataclass(frozen=True)
class CachedTargetAssemblyRoleMetadata:
    """Version binding plus exact NCBI sequence-report bytes for one UCSC database."""

    db: str
    assembly_accession: str
    assembly_description: CachedAssemblyRoleArtifact
    sequence_report: CachedAssemblyRoleArtifact

    def __post_init__(self) -> None:
        _validate_ucsc_db(self.db)
        expected_description = _ucsc_assembly_description_url(self.db)
        if self.assembly_description.source_url != expected_description:
            raise ValueError("cached assembly description does not match UCSC database")
        expected_report_url = _ncbi_sequence_report_download_url(
            self.assembly_accession
        )
        if self.sequence_report.source_url != expected_report_url:
            raise ValueError(
                "cached NCBI sequence report does not match assembly accession"
            )
        if self.sequence_report.archive_member != _ncbi_sequence_report_member(
            self.assembly_accession
        ):
            raise ValueError(
                "cached NCBI sequence report has unexpected archive member"
            )


@dataclass(frozen=True)
class CachedUCSCAssemblyMetadata:
    """Cached UCSC table artifacts used to define one assembly namespace."""

    db: str
    chrom_info: CachedResource
    chrom_alias: CachedResource | None = None

    def __post_init__(self) -> None:
        _validate_ucsc_db(self.db)
        expected = _metadata_urls(self.db)
        if self.chrom_info.source_url != expected.chrom_info_url:
            raise ValueError("cached chromInfo resource does not match UCSC database")
        if self.chrom_alias is not None and (
            self.chrom_alias.source_url != expected.chrom_alias_url
        ):
            raise ValueError("cached chromAlias resource does not match UCSC database")


def load_cached_ucsc_assembly_metadata(
    cache_root: ResourcePath,
    db: str,
) -> CachedUCSCAssemblyMetadata | None:
    """Load verified assembly metadata from cache without provider access.

    ``chromInfo`` is required. A cached ``chromAlias`` table is used when present;
    absence of alias bytes does not invalidate the canonical sequence namespace.
    """

    urls = _metadata_urls(db)
    chrom_info = load_cached_ucsc_resource(cache_root, urls.chrom_info_url)
    if chrom_info is None:
        return None
    chrom_alias = (
        load_cached_ucsc_resource(cache_root, urls.chrom_alias_url)
        if urls.chrom_alias_url is not None
        else None
    )
    return CachedUCSCAssemblyMetadata(
        db=db,
        chrom_info=chrom_info,
        chrom_alias=chrom_alias,
    )


def acquire_ucsc_assembly_metadata(
    resources: UCSCAssemblyMetadataResources,
    cache_root: ResourcePath,
    *,
    refresh: bool = False,
) -> CachedUCSCAssemblyMetadata:
    """Acquire discovered UCSC assembly tables into the shared content cache.

    UCSC's database download directory states that its files and tables are freely
    usable for any purpose, so these two metadata tables do not use the explicit
    evidence-resource terms acknowledgement gate.
    """

    chrom_info = acquire_ucsc_resource(
        resources.chrom_info_url,
        cache_root,
        terms_acknowledged=False,
        refresh=refresh,
    )
    chrom_alias = (
        acquire_ucsc_resource(
            resources.chrom_alias_url,
            cache_root,
            terms_acknowledged=False,
            refresh=refresh,
        )
        if resources.chrom_alias_url is not None
        else None
    )
    return CachedUCSCAssemblyMetadata(
        db=resources.db,
        chrom_info=chrom_info,
        chrom_alias=chrom_alias,
    )


def build_cached_ucsc_assembly_sequence_catalog(
    assembly: AssemblyIdentifier,
    metadata: CachedUCSCAssemblyMetadata,
) -> AssemblySequenceCatalog:
    """Parse a verified cached UCSC metadata pair into an immutable catalog."""

    if metadata.db != assembly.name and metadata.db not in assembly.aliases:
        raise ValueError(
            "assembly identifier does not represent cached UCSC metadata database"
        )

    sequence_provenance = _cached_metadata_provenance(
        metadata.chrom_info,
        label=f"UCSC {metadata.db} chromInfo assembly-sequence metadata",
    )
    expected_info_sha256 = sha256_hex_from_identifier(metadata.chrom_info.sha256)

    with open_text_resource(
        metadata.chrom_info.path,
        expected_sha256_hex=expected_info_sha256,
    ) as chrom_info_lines:
        if metadata.chrom_alias is None:
            return build_ucsc_assembly_sequence_catalog(
                assembly,
                chrom_info_lines,
                sequence_provenance=sequence_provenance,
            )

        alias_provenance = _cached_metadata_provenance(
            metadata.chrom_alias,
            label=f"UCSC {metadata.db} chromAlias sequence-alias metadata",
        )
        expected_alias_sha256 = sha256_hex_from_identifier(metadata.chrom_alias.sha256)
        with open_text_resource(
            metadata.chrom_alias.path,
            expected_sha256_hex=expected_alias_sha256,
        ) as chrom_alias_lines:
            return build_ucsc_assembly_sequence_catalog(
                assembly,
                chrom_info_lines,
                sequence_provenance=sequence_provenance,
                chrom_alias_lines=chrom_alias_lines,
                alias_provenance=alias_provenance,
            )


def _cached_metadata_provenance(
    resource: CachedResource,
    *,
    label: str,
) -> ProvenanceSource:
    identifier = ProvenanceIdentifier(
        kind=ProvenanceIdentifierKind.SHA256,
        value=resource.sha256,
    )
    return ProvenanceSource(
        source_id=f"file:{identifier.value}",
        label=label,
        identifiers=(identifier,),
    )


def load_cached_target_assembly_role_metadata(
    cache_root: ResourcePath,
    db: str,
) -> CachedTargetAssemblyRoleMetadata | None:
    """Load fully verified target-role metadata without provider access."""

    _validate_ucsc_db(db)
    description_url = _ucsc_assembly_description_url(db)
    description = _load_role_artifact(cache_root, description_url)
    if description is None:
        return None
    try:
        accession = parse_ucsc_assembly_description_accession(
            _read_verified_role_artifact_text(description)
        )
    except ValueError:
        return None
    report_url = _ncbi_sequence_report_download_url(accession)
    member = _ncbi_sequence_report_member(accession)
    report = _load_role_artifact(cache_root, report_url, archive_member=member)
    if report is None:
        return None
    return CachedTargetAssemblyRoleMetadata(
        db=db,
        assembly_accession=accession,
        assembly_description=description,
        sequence_report=report,
    )


def acquire_target_assembly_role_metadata(
    cache_root: ResourcePath,
    db: str,
    *,
    refresh: bool = False,
) -> CachedTargetAssemblyRoleMetadata:
    """Acquire the UCSC version binding and exact NCBI sequence-report artifact."""

    _validate_ucsc_db(db)
    description_url = _ucsc_assembly_description_url(db)
    description = _acquire_role_artifact(cache_root, description_url, refresh=refresh)
    try:
        accession = parse_ucsc_assembly_description_accession(
            _read_verified_role_artifact_text(description)
        )
    except ValueError as exc:
        raise AssemblyRoleMetadataAcquisitionError(
            f"UCSC assembly description for {db} does not state a versioned NCBI "
            "assembly accession"
        ) from exc
    report_url = _ncbi_sequence_report_download_url(accession)
    member = _ncbi_sequence_report_member(accession)
    report = _acquire_ncbi_sequence_report_artifact(
        cache_root,
        report_url,
        archive_member=member,
        refresh=refresh,
    )
    return CachedTargetAssemblyRoleMetadata(
        db=db,
        assembly_accession=accession,
        assembly_description=description,
        sequence_report=report,
    )


def attach_cached_target_role_context(
    catalog: AssemblySequenceCatalog,
    metadata: CachedTargetAssemblyRoleMetadata,
) -> AssemblySequenceCatalog:
    """Attach exact cached NCBI provider roles to an authoritative UCSC catalog."""

    if (
        metadata.db != catalog.assembly.name
        and metadata.db not in catalog.assembly.aliases
    ):
        raise ValueError(
            "target role metadata database does not match assembly catalog"
        )
    description_provenance = _cached_role_artifact_provenance(
        metadata.assembly_description,
        label=f"UCSC {metadata.db} assembly description",
    )
    role_provenance = _cached_role_artifact_provenance(
        metadata.sequence_report,
        label=(f"NCBI Datasets {metadata.assembly_accession} genome sequence report"),
        derived_from=(description_provenance,),
    )
    expected_report_sha256 = sha256_hex_from_identifier(metadata.sequence_report.sha256)
    with open_text_resource(
        metadata.sequence_report.path,
        expected_sha256_hex=expected_report_sha256,
    ) as lines:
        return attach_ncbi_sequence_role_context(
            catalog,
            lines,
            expected_assembly_accession=metadata.assembly_accession,
            role_provenance=role_provenance,
        )


def _read_verified_role_artifact_text(artifact: CachedAssemblyRoleArtifact) -> str:
    expected_sha256 = sha256_hex_from_identifier(artifact.sha256)
    with open_text_resource(
        artifact.path,
        expected_sha256_hex=expected_sha256,
    ) as lines:
        return "".join(lines)


def _cached_role_artifact_provenance(
    artifact: CachedAssemblyRoleArtifact,
    *,
    label: str,
    derived_from: tuple[ProvenanceSource, ...] = (),
) -> ProvenanceSource:
    identifier = ProvenanceIdentifier(
        kind=ProvenanceIdentifierKind.SHA256,
        value=artifact.sha256,
    )
    return ProvenanceSource(
        source_id=f"file:{identifier.value}",
        label=label,
        identifiers=(identifier,),
        derived_from=derived_from,
    )


def _ucsc_assembly_description_url(db: str) -> str:
    _validate_ucsc_db(db)
    return f"https://hgdownload.soe.ucsc.edu/gbdb/{db}/html/description.html"


def _ncbi_sequence_report_download_url(accession: str) -> str:
    if fullmatch(r"GC[AF]_[0-9]+\.[0-9]+", accession) is None:
        raise ValueError("NCBI assembly accession must be versioned GCA_ or GCF_")
    encoded = quote(accession, safe="")
    return (
        f"{_NCBI_DATASETS_BASE}/genome/accession/{encoded}/download?"
        "include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
    )


def _ncbi_sequence_report_member(accession: str) -> str:
    return f"ncbi_dataset/data/{accession}/sequence_report.jsonl"


def _role_cache_key(source_url: str, archive_member: str | None) -> str:
    identity = (
        source_url if archive_member is None else f"{source_url}#{archive_member}"
    )
    return hashlib.sha256(identity.encode("utf-8")).hexdigest()


def _role_index_path(
    cache_root: ResourcePath, source_url: str, archive_member: str | None
) -> Path:
    return (
        Path(cache_root)
        / "assembly-role"
        / "by-source"
        / (_role_cache_key(source_url, archive_member) + ".json")
    )


def _role_artifact_path(cache_root: ResourcePath, sha256_hex: str) -> Path:
    return (
        Path(cache_root)
        / "assembly-role"
        / "artifacts"
        / "sha256"
        / sha256_hex[:2]
        / sha256_hex
    )


def _load_role_artifact(
    cache_root: ResourcePath,
    source_url: str,
    *,
    archive_member: str | None = None,
) -> CachedAssemblyRoleArtifact | None:
    index_path = _role_index_path(cache_root, source_url, archive_member)
    if not index_path.is_file():
        return None
    try:
        payload = json.loads(index_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AssemblyRoleMetadataCacheIntegrityError(
            f"invalid target-role cache index for {source_url}"
        ) from exc
    if (
        not isinstance(payload, dict)
        or payload.get("schema_version") != _ROLE_CACHE_SCHEMA_VERSION
    ):
        raise AssemblyRoleMetadataCacheIntegrityError(
            f"unsupported target-role cache index for {source_url}"
        )
    expected_identity = (
        source_url if archive_member is None else f"{source_url}#{archive_member}"
    )
    recorded_identity = payload.get("source_identity")
    sha256_hex = payload.get("sha256")
    size_bytes = payload.get("size_bytes")
    retrieved_at = payload.get("retrieved_at")
    if (
        recorded_identity != expected_identity
        or not isinstance(sha256_hex, str)
        or fullmatch(r"[0-9a-f]{64}", sha256_hex) is None
        or not isinstance(size_bytes, int)
        or isinstance(size_bytes, bool)
        or size_bytes < 0
        or not isinstance(retrieved_at, str)
        or not retrieved_at
    ):
        raise AssemblyRoleMetadataCacheIntegrityError(
            f"invalid target-role cache metadata for {source_url}"
        )
    artifact_path = _role_artifact_path(cache_root, sha256_hex)
    if not artifact_path.is_file() or artifact_path.stat().st_size != size_bytes:
        raise AssemblyRoleMetadataCacheIntegrityError(
            f"target-role cache artifact is missing or has wrong size for {source_url}"
        )
    actual = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
    if actual != sha256_hex:
        raise AssemblyRoleMetadataCacheIntegrityError(
            f"target-role cache artifact SHA-256 mismatch for {source_url}"
        )
    return CachedAssemblyRoleArtifact(
        path=artifact_path,
        source_url=source_url,
        retrieved_at=retrieved_at,
        sha256=f"sha256:{sha256_hex}",
        size_bytes=size_bytes,
        cache_hit=True,
        archive_member=archive_member,
    )


def _acquire_role_artifact(
    cache_root: ResourcePath,
    source_url: str,
    *,
    refresh: bool,
) -> CachedAssemblyRoleArtifact:
    if not refresh:
        try:
            cached = _load_role_artifact(cache_root, source_url)
        except AssemblyRoleMetadataCacheIntegrityError:
            cached = None
        if cached is not None:
            return cached
    data = _fetch_url_bytes(source_url, accept="text/html,text/plain")
    return _publish_role_artifact(cache_root, source_url, data)


def _acquire_ncbi_sequence_report_artifact(
    cache_root: ResourcePath,
    source_url: str,
    *,
    archive_member: str,
    refresh: bool,
) -> CachedAssemblyRoleArtifact:
    if not refresh:
        try:
            cached = _load_role_artifact(
                cache_root, source_url, archive_member=archive_member
            )
        except AssemblyRoleMetadataCacheIntegrityError:
            cached = None
        if cached is not None:
            return cached
    package_bytes = _fetch_url_bytes(source_url, accept="application/zip")
    try:
        with zipfile.ZipFile(BytesIO(package_bytes)) as archive:
            report_bytes = archive.read(archive_member)
    except (KeyError, zipfile.BadZipFile) as exc:
        raise AssemblyRoleMetadataAcquisitionError(
            "NCBI Datasets package did not contain the expected sequence report"
        ) from exc
    return _publish_role_artifact(
        cache_root,
        source_url,
        report_bytes,
        archive_member=archive_member,
    )


def _fetch_url_bytes(source_url: str, *, accept: str) -> bytes:
    request = Request(
        source_url,
        headers={"User-Agent": _USER_AGENT, "Accept": accept},
        method="GET",
    )
    try:
        with urlopen(request, timeout=60) as response:
            return bytes(response.read())
    except HTTPError as exc:
        raise AssemblyRoleMetadataAcquisitionError(
            f"failed to acquire target-role metadata {source_url}: HTTP {exc.code}"
        ) from exc
    except (URLError, TimeoutError, OSError) as exc:
        raise AssemblyRoleMetadataAcquisitionError(
            f"failed to acquire target-role metadata {source_url}: {exc}"
        ) from exc


def _publish_role_artifact(
    cache_root: ResourcePath,
    source_url: str,
    data: bytes,
    *,
    archive_member: str | None = None,
) -> CachedAssemblyRoleArtifact:
    root = Path(cache_root)
    sha256_hex = hashlib.sha256(data).hexdigest()
    artifact_path = _role_artifact_path(root, sha256_hex)
    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    if (
        not artifact_path.is_file()
        or hashlib.sha256(artifact_path.read_bytes()).hexdigest() != sha256_hex
    ):
        fd, temp_name = tempfile.mkstemp(
            prefix="role-artifact-", dir=artifact_path.parent
        )
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(data)
            os.replace(temp_name, artifact_path)
        finally:
            if os.path.exists(temp_name):
                os.unlink(temp_name)
    retrieved_at = datetime.now(UTC).isoformat()
    source_identity = (
        source_url if archive_member is None else f"{source_url}#{archive_member}"
    )
    index_path = _role_index_path(root, source_url, archive_member)
    index_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": _ROLE_CACHE_SCHEMA_VERSION,
        "source_identity": source_identity,
        "sha256": sha256_hex,
        "size_bytes": len(data),
        "retrieved_at": retrieved_at,
    }
    fd, temp_name = tempfile.mkstemp(prefix="role-index-", dir=index_path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True)
            handle.write("\n")
        os.replace(temp_name, index_path)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
    return CachedAssemblyRoleArtifact(
        path=artifact_path,
        source_url=source_url,
        retrieved_at=retrieved_at,
        sha256=f"sha256:{sha256_hex}",
        size_bytes=len(data),
        cache_hit=False,
        archive_member=archive_member,
    )


def _metadata_urls(db: str) -> UCSCAssemblyMetadataResources:
    """Return canonical URLs only for cache lookup, never as existence evidence."""

    _validate_ucsc_db(db)
    base = urljoin(_UCSC_GOLDEN_PATH, f"{db}/database/")
    return UCSCAssemblyMetadataResources(
        db=db,
        chrom_info_url=urljoin(base, "chromInfo.txt.gz"),
        chrom_alias_url=urljoin(base, "chromAlias.txt.gz"),
    )


def _validate_ucsc_db(db: str) -> None:
    if not db or fullmatch(_UCSC_DB_PATTERN, db) is None:
        raise ValueError(
            "UCSC database identifier may contain only letters, digits, '.', '_', "
            "and '-'"
        )
