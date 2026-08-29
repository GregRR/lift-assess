"""Acquire and load authoritative UCSC assembly-sequence metadata.

The metadata path is deliberately separate from alignment evidence acquisition.
``chromInfo`` and ``chromAlias`` are small UCSC database tables used for input
preflight; they are not mapping evidence and are not part of an evidence tier.
"""

from __future__ import annotations

from dataclasses import dataclass
from re import fullmatch
from urllib.parse import urljoin

from .assembly_metadata import (
    AssemblySequenceCatalog,
    build_ucsc_assembly_sequence_catalog,
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

_UCSC_GOLDEN_PATH = "https://hgdownload.soe.ucsc.edu/goldenPath/"
_UCSC_DB_PATTERN = r"[A-Za-z0-9_.-]+"


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
