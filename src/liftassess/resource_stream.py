"""Verified streaming of local plain-text and gzip resources.

This internal module owns the raw-byte hashing/decompression boundary shared by normal
resource parsing and derived local indexes.  Callers that supply an expected SHA-256
pass the raw 64-character hexadecimal digest explicitly; canonical provenance
identifiers remain ``sha256:<hex>`` at higher layers.
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

from .resource_identity import ResourceIdentityMismatchError

ResourcePath: TypeAlias = str | PathLike[str]
RawReadProgressCallback: TypeAlias = Callable[[int], None]

_CHUNK_SIZE = 1024 * 1024


class _Digest(Protocol):
    def update(self, data: bytes) -> None: ...

    def hexdigest(self) -> str: ...


class _HashingRawReader(io.RawIOBase):
    """Read one binary stream while hashing exactly the bytes returned upstream."""

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


@contextmanager
def open_text_resource(
    path: ResourcePath,
    *,
    expected_sha256_hex: str | None = None,
    progress_callback: RawReadProgressCallback | None = None,
    progress_interval_bytes: int = _CHUNK_SIZE,
) -> Iterator[TextIO]:
    """Open one local text resource and optionally verify its exact raw SHA-256.

    Compression is detected from the first two bytes of the same open file handle.
    When ``expected_sha256_hex`` is supplied, the exact raw bytes feeding the parser
    are hashed and verified after successful parser exhaustion.
    """

    resource_path = Path(path)
    raw = io.FileIO(resource_path, mode="r")
    try:
        prefix = raw.read(2)
        raw.seek(0)
        is_gzip = prefix == b"\x1f\x8b"

        digest: _Digest | None = (
            hashlib.sha256() if expected_sha256_hex is not None else None
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
            if digest is not None and expected_sha256_hex is not None:
                # Parser exhaustion should already have reached EOF. Draining is
                # defensive and includes any raw trailing bytes in the exact artifact
                # identity rather than only the decoded text payload.
                while text.read(_CHUNK_SIZE):
                    pass
                if is_gzip:
                    while buffered.read(_CHUNK_SIZE):
                        pass
                actual = digest.hexdigest()
                if actual != expected_sha256_hex:
                    raise ResourceIdentityMismatchError(
                        f"SHA256 provenance mismatch for {resource_path}: "
                        f"expected {expected_sha256_hex}, got {actual}"
                    )
                assert hashing_raw is not None
                hashing_raw.finish_progress()
        finally:
            text.close()
            if is_gzip:
                buffered.close()
    finally:
        raw.close()
