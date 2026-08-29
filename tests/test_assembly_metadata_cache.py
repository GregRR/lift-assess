from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCAssemblyMetadata,
    GenomicInterval,
    ProvenanceIdentifierKind,
    ResourceIdentityMismatchError,
    SourceIntervalPreflightState,
    build_cached_ucsc_assembly_sequence_catalog,
    preflight_source_interval,
    sha256_identifier_for_file,
    ucsc_resource_terms,
)

_DB = "canFam3"
_CHROM_INFO_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/chromInfo.txt.gz"
)
_CHROM_ALIAS_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/chromAlias.txt.gz"
)


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _cached_resource(path: Path, url: str) -> CachedResource:
    return CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-28T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=True,
    )


def _metadata(tmp_path: Path) -> CachedUCSCAssemblyMetadata:
    chrom_info = tmp_path / "chrom-info"
    chrom_alias = tmp_path / "chrom-alias"
    _write_gzip(
        chrom_info,
        "chr1\t1000\t/gbdb/canFam3/canFam3.2bit\n"
        "chrNoChain\t200\t/gbdb/canFam3/canFam3.2bit\n",
    )
    _write_gzip(chrom_alias, "1\tchr1\tassembly\n")
    return CachedUCSCAssemblyMetadata(
        db=_DB,
        chrom_info=_cached_resource(chrom_info, _CHROM_INFO_URL),
        chrom_alias=_cached_resource(chrom_alias, _CHROM_ALIAS_URL),
    )


def test_cached_metadata_builds_catalog_with_exact_file_provenance(
    tmp_path: Path,
) -> None:
    metadata = _metadata(tmp_path)
    assembly = AssemblyIdentifier(_DB, "UCSC")

    catalog = build_cached_ucsc_assembly_sequence_catalog(assembly, metadata)
    valid = preflight_source_interval(
        GenomicInterval(assembly, "chrNoChain", 10, 20), catalog
    )
    alias = preflight_source_interval(GenomicInterval(assembly, "1", 10, 20), catalog)

    assert valid.state is SourceIntervalPreflightState.VALID
    assert valid.sequence_length == 200
    assert valid.provenance_sources == (catalog.sequence_provenance,)
    sequence_identifier = catalog.sequence_provenance.identifiers[0]
    assert sequence_identifier.kind is ProvenanceIdentifierKind.SHA256
    assert sequence_identifier.value == metadata.chrom_info.sha256

    assert alias.state is SourceIntervalPreflightState.UNRECOGNIZED_SOURCE_SEQUENCE_NAME
    assert alias.suggested_sequence_name == "chr1"
    assert alias.provenance_sources == (
        catalog.sequence_provenance,
        catalog.alias_provenance,
    )
    assert catalog.alias_provenance is not None
    assert metadata.chrom_alias is not None
    assert catalog.alias_provenance.identifiers[0].value == metadata.chrom_alias.sha256


def test_cached_metadata_catalog_reverifies_raw_bytes_while_parsing(
    tmp_path: Path,
) -> None:
    metadata = _metadata(tmp_path)
    _write_gzip(
        metadata.chrom_info.path,
        "chr1\t1001\t/gbdb/canFam3/canFam3.2bit\n"
        "chrNoChain\t200\t/gbdb/canFam3/canFam3.2bit\n",
    )

    with pytest.raises(
        ResourceIdentityMismatchError, match="SHA256 provenance mismatch"
    ):
        build_cached_ucsc_assembly_sequence_catalog(
            AssemblyIdentifier(_DB, "UCSC"), metadata
        )


def test_cached_metadata_requires_matching_ucsc_database(tmp_path: Path) -> None:
    metadata = _metadata(tmp_path)

    with pytest.raises(ValueError, match="does not represent cached UCSC metadata"):
        build_cached_ucsc_assembly_sequence_catalog(
            AssemblyIdentifier("canFam4", "UCSC"), metadata
        )


def test_acquisition_uses_only_provider_discovered_metadata_urls(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    import liftassess.assembly_metadata_cache as metadata_cache
    from liftassess import UCSCAssemblyMetadataResources

    metadata = _metadata(tmp_path)
    assert metadata.chrom_alias is not None
    by_url = {
        metadata.chrom_info.source_url: metadata.chrom_info,
        metadata.chrom_alias.source_url: metadata.chrom_alias,
    }
    calls: list[tuple[str, bool]] = []

    def acquire(url: str, cache_root: object, **kwargs: object) -> CachedResource:
        del cache_root
        calls.append((url, bool(kwargs["terms_acknowledged"])))
        resource = by_url[url]
        assert resource is not None
        return resource

    monkeypatch.setattr(metadata_cache, "acquire_ucsc_resource", acquire)
    discovered = UCSCAssemblyMetadataResources(
        db=_DB,
        chrom_info_url=_CHROM_INFO_URL,
        chrom_alias_url=_CHROM_ALIAS_URL,
    )

    result = metadata_cache.acquire_ucsc_assembly_metadata(discovered, tmp_path)

    assert result == metadata
    assert calls == [(_CHROM_INFO_URL, False), (_CHROM_ALIAS_URL, False)]


def test_cached_metadata_loader_requires_chrom_info_but_alias_is_optional(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    import liftassess.assembly_metadata_cache as metadata_cache

    metadata = _metadata(tmp_path)
    seen: list[str] = []

    def load(cache_root: object, url: str) -> CachedResource | None:
        del cache_root
        seen.append(url)
        if url == _CHROM_INFO_URL:
            return metadata.chrom_info
        return None

    monkeypatch.setattr(metadata_cache, "load_cached_ucsc_resource", load)

    result = metadata_cache.load_cached_ucsc_assembly_metadata(tmp_path, _DB)

    assert result == CachedUCSCAssemblyMetadata(db=_DB, chrom_info=metadata.chrom_info)
    assert seen == [_CHROM_INFO_URL, _CHROM_ALIAS_URL]
