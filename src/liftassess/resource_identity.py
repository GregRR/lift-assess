"""Content identity and integrity checks for local evidence-resource files.

liftAssess uses SHA-256 over the exact file bytes as the canonical identity for
non-sequence resource artifacts such as chain and net files.  Provider-published
checksums may use other algorithms (for example UCSC currently publishes MD5 for
some download directories); those checksums are useful for transfer-integrity
verification but are not provenance identifiers in the v1 model.

These helpers hash the raw on-disk bytes, including compression bytes when the
resource is compressed.  A gzip file and an independently recompressed copy may
therefore have different file identities even when they decode to the same text.
That is intentional: this layer identifies the exact external artifact consumed,
not a normalized semantic representation of its contents.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Callable, Iterable
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Protocol, TypeAlias

from .models import (
    ProvenanceIdentifier,
    ProvenanceIdentifierKind,
    ProvenanceSource,
)

ResourcePath: TypeAlias = str | PathLike[str]
ResourceChecksumProgressCallback: TypeAlias = Callable[[int, int], None]

_HEX_RE = re.compile(r"^[0-9a-fA-F]+$")
_CHUNK_SIZE = 1024 * 1024


class ResourceChecksumAlgorithm(str, Enum):
    """Algorithms accepted for external file-integrity verification."""

    MD5 = "md5"
    SHA256 = "sha256"


class ResourceChecksumMismatchError(ValueError):
    """Raised when a local resource does not match an expected checksum."""


class ResourceIdentityMismatchError(ValueError):
    """Raised when parsed file bytes do not match their provenance identity."""


class _Digest(Protocol):
    def update(self, data: bytes) -> None: ...

    def hexdigest(self) -> str: ...


def compute_resource_checksum(
    path: ResourcePath,
    algorithm: ResourceChecksumAlgorithm,
    *,
    progress_callback: ResourceChecksumProgressCallback | None = None,
) -> str:
    """Return a lowercase hexadecimal checksum of the exact local file bytes.

    When supplied, ``progress_callback`` receives the cumulative raw bytes hashed and
    the exact on-disk file size.  The callback reports checksum work only; callers that
    compare the resulting digest with an expected identity remain responsible for
    deciding when integrity verification has actually succeeded.
    """

    resource_path = Path(path)
    total_bytes = resource_path.stat().st_size
    bytes_hashed = 0
    if progress_callback is not None:
        progress_callback(0, total_bytes)

    digest = _new_digest(algorithm)
    with resource_path.open("rb") as handle:
        while chunk := handle.read(_CHUNK_SIZE):
            digest.update(chunk)
            bytes_hashed += len(chunk)
            if progress_callback is not None:
                progress_callback(bytes_hashed, total_bytes)
    return digest.hexdigest()


def verify_resource_checksum(
    path: ResourcePath,
    *,
    algorithm: ResourceChecksumAlgorithm,
    expected: str,
) -> str:
    """Verify an external/provider checksum and return the computed checksum.

    This is an integrity check, not a provenance decision.  In particular, an MD5
    published by a provider can confirm that downloaded bytes match the provider's
    advertised file, but liftAssess still records SHA-256 as the file's canonical
    provenance identity.
    """

    expected_normalized = _normalize_expected_checksum(expected, algorithm)
    actual = compute_resource_checksum(path, algorithm)
    if actual != expected_normalized:
        raise ResourceChecksumMismatchError(
            f"{algorithm.value} checksum mismatch for {Path(path)}: "
            f"expected {expected_normalized}, got {actual}"
        )
    return actual


def sha256_identifier_for_file(path: ResourcePath) -> ProvenanceIdentifier:
    """Return the canonical v1 provenance identifier for one local file artifact."""

    checksum = compute_resource_checksum(path, ResourceChecksumAlgorithm.SHA256)
    return ProvenanceIdentifier(
        kind=ProvenanceIdentifierKind.SHA256,
        value=f"sha256:{checksum}",
    )


def provenance_source_for_file(
    path: ResourcePath,
    *,
    label: str,
    derived_from: Iterable[ProvenanceSource],
) -> ProvenanceSource:
    """Create a content-addressed provenance node for one exact local file.

    The structural ``source_id`` is derived from the same SHA-256 identifier stored
    on the node.  This makes identical bytes loaded through different paths converge
    on the same file-source identity instead of relying on a caller-chosen filename or
    label. ``derived_from`` is deliberately required because a file digest can identify
    the artifact but cannot infer the alignment, pipeline, or other upstream process
    that produced it; callers must make that relationship explicit, including choosing
    an empty tuple when no upstream source is claimed.
    """

    identifier = sha256_identifier_for_file(path)
    return ProvenanceSource(
        source_id=f"file:{identifier.value}",
        label=label,
        identifiers=(identifier,),
        derived_from=tuple(derived_from),
    )


def _sha256_checksum_from_file_provenance(source: ProvenanceSource) -> str:
    """Return the canonical file SHA-256 recorded by one provenance node.

    File-backed parsing requires exactly one SHA-256 identifier. The
    :func:`provenance_source_for_file` helper also content-addresses ``source_id``, but
    parsing verifies the typed digest itself so older/explicit logical source IDs do
    not silently substitute for exact artifact identity.
    """

    identifiers = tuple(
        identifier
        for identifier in source.identifiers
        if identifier.kind is ProvenanceIdentifierKind.SHA256
    )
    if len(identifiers) != 1:
        raise ValueError(
            "file provenance must contain exactly one canonical SHA256 identifier"
        )

    return identifiers[0].value.removeprefix("sha256:")


def _new_digest(algorithm: ResourceChecksumAlgorithm) -> _Digest:
    if algorithm is ResourceChecksumAlgorithm.MD5:
        # MD5 is supported only for compatibility with provider-published integrity
        # metadata, never as a liftAssess provenance identity or security primitive.
        return hashlib.md5(usedforsecurity=False)
    return hashlib.sha256()


def _normalize_expected_checksum(
    expected: str,
    algorithm: ResourceChecksumAlgorithm,
) -> str:
    value = expected.strip()
    expected_length = 32 if algorithm is ResourceChecksumAlgorithm.MD5 else 64
    if len(value) != expected_length or _HEX_RE.fullmatch(value) is None:
        raise ValueError(
            f"expected {algorithm.value} checksum must be exactly "
            f"{expected_length} hexadecimal characters"
        )
    return value.lower()
