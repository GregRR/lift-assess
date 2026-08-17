from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest

from liftassess import (
    ProvenanceIdentifier,
    ProvenanceIdentifierKind,
    ProvenanceSource,
    ResourceChecksumAlgorithm,
    ResourceChecksumMismatchError,
    compute_resource_checksum,
    provenance_source_for_file,
    sha256_identifier_for_file,
    verify_resource_checksum,
)


def test_sha256_identifier_hashes_exact_compressed_file_bytes(tmp_path: Path) -> None:
    path = tmp_path / "example.chain.gz"
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write("chain bytes\n")

    expected = hashlib.sha256(path.read_bytes()).hexdigest()

    identifier = sha256_identifier_for_file(path)

    assert identifier.kind is ProvenanceIdentifierKind.SHA256
    assert identifier.value == f"sha256:{expected}"


def test_file_provenance_is_content_addressed_and_preserves_upstream_source(
    tmp_path: Path,
) -> None:
    first_path = tmp_path / "first.chain"
    second_path = tmp_path / "renamed.chain"
    first_path.write_bytes(b"identical bytes")
    second_path.write_bytes(b"identical bytes")
    alignment = ProvenanceSource("alignment", "shared alignment")

    first = provenance_source_for_file(
        first_path,
        label="first local name",
        derived_from=(alignment,),
    )
    second = provenance_source_for_file(
        second_path,
        label="second local name",
        derived_from=(alignment,),
    )

    assert first.source_id == second.source_id
    assert first.identifiers == second.identifiers
    assert first.derived_from == (alignment,)
    assert second.derived_from == (alignment,)


def test_different_file_bytes_produce_different_provenance_identity(
    tmp_path: Path,
) -> None:
    first_path = tmp_path / "first.chain"
    second_path = tmp_path / "second.chain"
    first_path.write_bytes(b"first bytes")
    second_path.write_bytes(b"second bytes")

    first = provenance_source_for_file(first_path, label="first", derived_from=())
    second = provenance_source_for_file(second_path, label="second", derived_from=())

    assert first.source_id != second.source_id
    assert first.identifiers != second.identifiers


def test_provider_md5_can_be_verified_without_becoming_provenance_identity(
    tmp_path: Path,
) -> None:
    path = tmp_path / "resource.chain.gz"
    path.write_bytes(b"provider resource bytes")
    expected = hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()

    actual = verify_resource_checksum(
        path,
        algorithm=ResourceChecksumAlgorithm.MD5,
        expected=expected.upper(),
    )
    provenance = provenance_source_for_file(path, label="resource", derived_from=())

    assert actual == expected
    assert provenance.identifiers[0].kind is ProvenanceIdentifierKind.SHA256
    assert provenance.identifiers[0].value.startswith("sha256:")


def test_provider_checksum_mismatch_fails_loudly(tmp_path: Path) -> None:
    path = tmp_path / "resource.chain"
    path.write_bytes(b"actual bytes")

    with pytest.raises(ResourceChecksumMismatchError, match="checksum mismatch"):
        verify_resource_checksum(
            path,
            algorithm=ResourceChecksumAlgorithm.MD5,
            expected="0" * 32,
        )


def test_expected_checksum_must_have_algorithm_specific_hex_shape(
    tmp_path: Path,
) -> None:
    path = tmp_path / "resource.chain"
    path.write_bytes(b"bytes")

    with pytest.raises(ValueError, match="32 hexadecimal characters"):
        verify_resource_checksum(
            path,
            algorithm=ResourceChecksumAlgorithm.MD5,
            expected="not-an-md5",
        )


def test_compute_sha256_returns_raw_lowercase_hex(tmp_path: Path) -> None:
    path = tmp_path / "resource.net"
    path.write_bytes(b"net bytes")

    checksum = compute_resource_checksum(path, ResourceChecksumAlgorithm.SHA256)

    assert checksum == hashlib.sha256(b"net bytes").hexdigest()
    assert not checksum.startswith("sha256:")


def test_compute_checksum_reports_exact_bytes_hashed(tmp_path: Path) -> None:
    path = tmp_path / "resource.net"
    data = b"abcdef" * 10
    path.write_bytes(data)
    progress: list[tuple[int, int]] = []

    checksum = compute_resource_checksum(
        path,
        ResourceChecksumAlgorithm.SHA256,
        progress_callback=lambda hashed, total: progress.append((hashed, total)),
    )

    assert checksum == hashlib.sha256(data).hexdigest()
    assert progress[0] == (0, len(data))
    assert progress[-1] == (len(data), len(data))
    assert all(total == len(data) for _, total in progress)
    assert [hashed for hashed, _ in progress] == sorted(
        hashed for hashed, _ in progress
    )


def test_sha256_provenance_identifier_requires_canonical_form() -> None:
    with pytest.raises(ValueError, match="canonical"):
        ProvenanceIdentifier(ProvenanceIdentifierKind.SHA256, "a" * 64)

    with pytest.raises(ValueError, match="canonical"):
        ProvenanceIdentifier(
            ProvenanceIdentifierKind.SHA256,
            "sha256:" + "A" * 64,
        )
