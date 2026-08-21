from __future__ import annotations

import gzip
import json
import sqlite3
from pathlib import Path

import pytest

from liftassess.chain_index import (
    ChainIndexCorruptionError,
    build_cached_chain_index,
    build_chain_index,
    chain_index_cache_path,
    load_cached_chain_index,
    load_chain_index,
)
from liftassess.models import AssemblyIdentifier, GenomicInterval
from liftassess.resource_cache import CachedResource, ucsc_resource_terms
from liftassess.resource_files import iter_chain_file
from liftassess.resource_identity import (
    ResourceIdentityMismatchError,
    sha256_hex_from_identifier,
    sha256_identifier_for_file,
)


def _assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier(name="source", provider="test")


def _chain_text() -> str:
    return """\
chain 10 chr1 200000 + 65530 65555 chrA 300000 - 1000 1028 1
10 5 8
10

chain 20 chr1 200000 + 70000 70010 chrB 300000 + 2000 2010 2
10

chain 30 chr2 100000 + 100 110 chrC 100000 + 300 310 3
10

"""


def _write_chain(path: Path, *, gzip_compressed: bool = False) -> None:
    if gzip_compressed:
        with gzip.open(path, "wt", encoding="ascii", newline="\n") as handle:
            handle.write(_chain_text())
    else:
        path.write_text(_chain_text(), encoding="ascii")


def _identifier_and_size(path: Path) -> tuple[str, int]:
    identifier = sha256_identifier_for_file(path)
    return identifier.value, path.stat().st_size


def test_sha256_identifier_helper_rejects_double_prefix() -> None:
    valid = "sha256:" + "a" * 64

    assert sha256_hex_from_identifier(valid) == "a" * 64
    with pytest.raises(ValueError, match="canonical"):
        sha256_hex_from_identifier("sha256:" + valid)


def test_chain_index_cache_path_uses_canonical_digest(tmp_path: Path) -> None:
    identifier = "sha256:" + "a" * 64

    path = chain_index_cache_path(tmp_path, identifier)

    assert path == (
        tmp_path / "derived" / "chain-index-v2" / "sha256" / "aa" / ("a" * 64)
    )


def test_build_and_query_chain_index_preserves_records_and_order(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)

    result = build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    )

    interval = GenomicInterval(
        assembly=_assembly(),
        sequence_name="chr1",
        start=65535,
        end=70005,
    )
    indexed = result.index.records_for_interval(interval)
    full = tuple(
        record
        for record in iter_chain_file(chain_path)
        if record.target_name == "chr1"
        and record.target_start < interval.end
        and record.target_end > interval.start
    )
    assert indexed == full
    assert [record.chain_id for record in indexed] == [1, 2]
    assert result.manifest.source_chain_sha256_identifier == identifier
    assert result.manifest.record_count == 3


def test_chain_crossing_bin_boundary_is_discoverable_from_both_bins(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "example.chain.gz"
    index_path = tmp_path / "index"
    _write_chain(chain_path, gzip_compressed=True)
    identifier, size_bytes = _identifier_and_size(chain_path)
    index = build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    ).index

    left = GenomicInterval(_assembly(), "chr1", 65535, 65536)
    right = GenomicInterval(_assembly(), "chr1", 65536, 65537)

    assert [record.chain_id for record in index.records_for_interval(left)] == [1]
    assert [record.chain_id for record in index.records_for_interval(right)] == [1]


def test_chain_index_preserves_full_traversal_source_bound_failure(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)
    index = build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    ).index
    interval = GenomicInterval(_assembly(), "chr1", 199_999, 200_001)

    with pytest.raises(ValueError, match="source interval exceeds"):
        index.records_for_interval(interval)


def test_build_rejects_source_identity_mismatch_without_publishing_index(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)

    with pytest.raises(
        ResourceIdentityMismatchError, match="SHA256 provenance mismatch"
    ):
        build_chain_index(
            chain_path,
            index_path,
            source_chain_sha256_identifier="sha256:" + "0" * 64,
            source_chain_size_bytes=chain_path.stat().st_size,
        )

    assert not index_path.exists()


def test_build_rejects_source_size_mismatch(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)

    with pytest.raises(
        ResourceIdentityMismatchError, match="source chain size mismatch"
    ):
        build_chain_index(
            chain_path,
            index_path,
            source_chain_sha256_identifier=identifier,
            source_chain_size_bytes=size_bytes + 1,
        )

    assert not index_path.exists()


def test_load_rejects_double_prefixed_source_identity_in_manifest(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)
    build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    )
    manifest_path = index_path / "manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["source_chain_sha256_identifier"] = "sha256:" + identifier
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ChainIndexCorruptionError, match="manifest is invalid"):
        load_chain_index(index_path)


def test_load_rejects_database_mutation(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)
    build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    )
    database_path = index_path / "index.sqlite3"
    with database_path.open("r+b") as handle:
        handle.seek(100)
        original = handle.read(1)
        handle.seek(100)
        handle.write(bytes([original[0] ^ 1]))

    with pytest.raises(ChainIndexCorruptionError, match="database SHA256 mismatch"):
        load_chain_index(index_path)


def test_fast_load_rejects_lookup_catalog_mutation(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)
    build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    )
    catalog_path = index_path / "lookup-catalog.json"
    with catalog_path.open("r+b") as handle:
        handle.seek(20)
        original = handle.read(1)
        handle.seek(20)
        handle.write(bytes([original[0] ^ 1]))

    with pytest.raises(ChainIndexCorruptionError, match="lookup catalog SHA256"):
        load_chain_index(index_path, verify_database=False)


def test_query_rejects_queried_bin_database_mutation_without_full_hash(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)
    build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    )
    database_path = index_path / "index.sqlite3"
    with sqlite3.connect(database_path) as connection:
        connection.execute(
            "DELETE FROM bin_memberships WHERE sequence_id = 0 AND record_id = 0"
        )
        connection.commit()

    index = load_chain_index(index_path, verify_database=False)
    interval = GenomicInterval(_assembly(), "chr1", 65535, 65536)

    with pytest.raises(ChainIndexCorruptionError, match="bin lookup integrity"):
        index.records_for_interval(interval)


def test_query_allows_unqueried_database_mutation_without_full_hash(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)
    build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    )
    database_path = index_path / "index.sqlite3"
    with sqlite3.connect(database_path) as connection:
        connection.execute(
            "DELETE FROM bin_memberships WHERE sequence_id = 1 AND record_id = 2"
        )
        connection.commit()

    index = load_chain_index(index_path, verify_database=False)
    interval = GenomicInterval(_assembly(), "chr1", 65535, 65536)

    assert [record.chain_id for record in index.records_for_interval(interval)] == [1]


def test_query_rejects_selected_block_mutation(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain"
    index_path = tmp_path / "index"
    _write_chain(chain_path)
    identifier, size_bytes = _identifier_and_size(chain_path)
    index = build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
    ).index
    block_path = index_path / "records.blocks"
    with block_path.open("r+b") as handle:
        original = handle.read(1)
        handle.seek(0)
        handle.write(bytes([original[0] ^ 1]))

    interval = GenomicInterval(_assembly(), "chr1", 65535, 65536)
    with pytest.raises(ChainIndexCorruptionError, match="block SHA256 mismatch"):
        index.records_for_interval(interval)


def test_cached_chain_index_helpers_use_resource_identity_path(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain.gz"
    cache_root = tmp_path / "cache"
    _write_chain(chain_path, gzip_compressed=True)
    identifier, size_bytes = _identifier_and_size(chain_path)
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
        "canFam3ToCanFam4.over.chain.gz"
    )
    resource = CachedResource(
        path=chain_path,
        source_url=url,
        retrieved_at="2026-08-20T00:00:00Z",
        sha256=identifier,
        size_bytes=size_bytes,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=True,
    )

    built = build_cached_chain_index(cache_root, resource)
    loaded = load_cached_chain_index(cache_root, resource)

    assert built.index.root == chain_index_cache_path(cache_root, identifier)
    assert loaded is not None
    assert loaded.manifest == built.manifest


def test_build_progress_reports_exact_raw_source_bytes(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain.gz"
    index_path = tmp_path / "index"
    _write_chain(chain_path, gzip_compressed=True)
    identifier, size_bytes = _identifier_and_size(chain_path)
    progress: list[int] = []

    build_chain_index(
        chain_path,
        index_path,
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=size_bytes,
        progress_callback=progress.append,
    )

    assert progress
    assert progress[-1] == size_bytes
    assert progress == sorted(progress)
