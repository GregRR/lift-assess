"""Reusable region-addressable indexes for UCSC chain resources.

Milestone 18 benchmarks selected a two-level local representation for large chain
resources: 64-KiB source-coordinate bin memberships plus encounter-order chain records
stored exactly once in independently compressed blocks.  This module implements that
representation without changing chain parsing or projection semantics.

An index is a derived acceleration artifact, never replacement scientific evidence.
Its manifest is bound to the canonical SHA-256 identifier and exact byte size of the
source chain.  Index construction hashes the exact raw source bytes while parsing them.
Normal index loading verifies a compact lookup-integrity catalog; each query validates
the exact genomic bins it reads against that catalog, and each selected compressed
payload block is verified before parsing. A full SQLite SHA-256 remains available for
explicit deep verification without becoming a per-query whole-file read.
"""

from __future__ import annotations

import hashlib
import io
import json
import shutil
import sqlite3
import struct
import tempfile
import zlib
from collections.abc import Callable, Mapping
from contextlib import closing
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Protocol

from .chain import ChainRecord, iter_chain_records
from .models import GenomicInterval
from .resource_cache import CachedResource
from .resource_identity import (
    ResourceChecksumAlgorithm,
    ResourceIdentityMismatchError,
    compute_resource_checksum,
    sha256_hex_from_identifier,
)
from .resource_stream import ResourcePath, open_text_resource

_INDEX_FORMAT = "liftassess-chain-index"
# v1 used whole-database verification; v2 adds authenticated query-local lookup.
_INDEX_SCHEMA_VERSION = 2
_LOOKUP_CATALOG_FORMAT = "liftassess-chain-index-lookup-catalog"
_LOOKUP_CATALOG_SCHEMA_VERSION = 1
DEFAULT_CHAIN_INDEX_BIN_WIDTH = 65_536
DEFAULT_CHAIN_INDEX_BLOCK_SIZE = 1024 * 1024
_CHAIN_BATCH_SIZE = 10_000
_MEMBERSHIP_BATCH_SIZE = 100_000


class ChainIndexError(RuntimeError):
    """Base class for reusable chain-index failures."""


class ChainIndexCorruptionError(ChainIndexError):
    """A derived chain index no longer matches its recorded integrity metadata."""


class _Digest(Protocol):
    def update(self, data: bytes, /) -> None: ...

    def hexdigest(self) -> str: ...


@dataclass
class _BinDigestAccumulator:
    count: int = 0
    digest: _Digest = field(default_factory=hashlib.sha256)

    def update(
        self,
        *,
        sequence_id: int,
        bin_id: int,
        record_id: int,
        target_start: int,
        target_end: int,
        block_id: int,
        block_offset: int,
        record_length: int,
    ) -> None:
        self.digest.update(
            _bin_integrity_payload(
                sequence_id=sequence_id,
                bin_id=bin_id,
                record_id=record_id,
                target_start=target_start,
                target_end=target_end,
                block_id=block_id,
                block_offset=block_offset,
                record_length=record_length,
            )
        )
        self.count += 1


@dataclass(frozen=True)
class _LookupSequence:
    sequence_id: int
    name: str
    minimum_target_size: int


@dataclass(frozen=True)
class _LookupBlock:
    block_id: int
    file_offset: int
    compressed_length: int
    uncompressed_length: int
    compressed_sha256_hex: str


@dataclass(frozen=True)
class _LookupBinIntegrity:
    count: int
    sha256_hex: str


@dataclass(frozen=True)
class _ChainIndexLookupCatalog:
    sequences_by_name: Mapping[str, _LookupSequence]
    blocks_by_id: Mapping[int, _LookupBlock]
    bins: Mapping[tuple[int, int], _LookupBinIntegrity]


@dataclass(frozen=True)
class ChainIndexManifest:
    """Validated metadata binding one derived index to one exact source chain."""

    source_chain_sha256_identifier: str
    source_chain_size_bytes: int
    bin_width: int
    block_size: int
    record_count: int
    membership_count: int
    block_count: int
    database_size_bytes: int
    block_store_size_bytes: int
    lookup_catalog_size_bytes: int
    database_sha256_identifier: str
    lookup_catalog_sha256_identifier: str

    def __post_init__(self) -> None:
        sha256_hex_from_identifier(self.source_chain_sha256_identifier)
        sha256_hex_from_identifier(self.database_sha256_identifier)
        sha256_hex_from_identifier(self.lookup_catalog_sha256_identifier)
        if self.source_chain_size_bytes < 0:
            raise ValueError("source chain size must not be negative")
        if self.bin_width <= 0:
            raise ValueError("chain index bin width must be positive")
        if self.block_size <= 0:
            raise ValueError("chain index block size must be positive")
        for name, value in (
            ("record count", self.record_count),
            ("membership count", self.membership_count),
            ("block count", self.block_count),
            ("database size", self.database_size_bytes),
            ("block-store size", self.block_store_size_bytes),
            ("lookup-catalog size", self.lookup_catalog_size_bytes),
        ):
            if value < 0:
                raise ValueError(f"chain index {name} must not be negative")

    def to_json_payload(self) -> dict[str, object]:
        return {
            "format": _INDEX_FORMAT,
            "schema_version": _INDEX_SCHEMA_VERSION,
            "source_chain_sha256_identifier": self.source_chain_sha256_identifier,
            "source_chain_size_bytes": self.source_chain_size_bytes,
            "bin_width": self.bin_width,
            "block_size": self.block_size,
            "record_count": self.record_count,
            "membership_count": self.membership_count,
            "block_count": self.block_count,
            "database_size_bytes": self.database_size_bytes,
            "block_store_size_bytes": self.block_store_size_bytes,
            "lookup_catalog_size_bytes": self.lookup_catalog_size_bytes,
            "database_sha256_identifier": self.database_sha256_identifier,
            "lookup_catalog_sha256_identifier": (self.lookup_catalog_sha256_identifier),
        }

    @classmethod
    def from_json_payload(cls, payload: object) -> ChainIndexManifest:
        if not isinstance(payload, dict):
            raise TypeError("chain index manifest must be a JSON object")
        if payload.get("format") != _INDEX_FORMAT:
            raise ValueError("unsupported chain index format")
        if payload.get("schema_version") != _INDEX_SCHEMA_VERSION:
            raise ValueError("unsupported chain index schema version")

        return cls(
            source_chain_sha256_identifier=_required_str(
                payload, "source_chain_sha256_identifier"
            ),
            source_chain_size_bytes=_required_int(payload, "source_chain_size_bytes"),
            bin_width=_required_int(payload, "bin_width"),
            block_size=_required_int(payload, "block_size"),
            record_count=_required_int(payload, "record_count"),
            membership_count=_required_int(payload, "membership_count"),
            block_count=_required_int(payload, "block_count"),
            database_size_bytes=_required_int(payload, "database_size_bytes"),
            block_store_size_bytes=_required_int(payload, "block_store_size_bytes"),
            lookup_catalog_size_bytes=_required_int(
                payload, "lookup_catalog_size_bytes"
            ),
            database_sha256_identifier=_required_str(
                payload, "database_sha256_identifier"
            ),
            lookup_catalog_sha256_identifier=_required_str(
                payload, "lookup_catalog_sha256_identifier"
            ),
        )


@dataclass(frozen=True)
class ChainIndexBuildResult:
    """One successfully built reusable index and its durable manifest."""

    index: ChainIndex
    manifest: ChainIndexManifest


@dataclass(frozen=True)
class ChainIndex:
    """One validated chain index ready for region-addressable queries."""

    root: Path
    manifest: ChainIndexManifest
    _lookup_catalog: _ChainIndexLookupCatalog

    @property
    def database_path(self) -> Path:
        return self.root / "index.sqlite3"

    @property
    def block_store_path(self) -> Path:
        return self.root / "records.blocks"

    @property
    def manifest_path(self) -> Path:
        return self.root / "manifest.json"

    @property
    def lookup_catalog_path(self) -> Path:
        return self.root / "lookup-catalog.json"

    def records_for_interval(
        self, interval: GenomicInterval
    ) -> tuple[ChainRecord, ...]:
        """Return exact overlapping chain records in original encounter order.

        Normal queries do not hash the complete SQLite database. Instead, they validate
        every queried bin's membership/record-locator rows against the authenticated
        lookup catalog loaded with the index, then verify each selected compressed
        record block before parsing it.
        """

        if interval.length == 0:
            raise ValueError(
                "zero-length source interval projection is not defined for "
                "liftAssess v1"
            )

        sequence = self._lookup_catalog.sequences_by_name.get(interval.sequence_name)
        if sequence is None:
            return ()
        if interval.end > sequence.minimum_target_size:
            raise ValueError("source interval exceeds chain target sequence bounds")

        first_bin = interval.start // self.manifest.bin_width
        last_bin = (interval.end - 1) // self.manifest.bin_width
        database_uri = f"file:{self.database_path}?mode=ro"
        try:
            with closing(sqlite3.connect(database_uri, uri=True)) as connection:
                rows = tuple(
                    connection.execute(
                        """
                        SELECT m.bin_id, c.record_id, c.target_start, c.target_end,
                               c.block_id, c.block_offset, c.record_length
                        FROM bin_memberships AS m
                        JOIN chains AS c ON c.record_id = m.record_id
                        WHERE m.sequence_id = ? AND m.bin_id BETWEEN ? AND ?
                        ORDER BY m.bin_id, c.record_id
                        """,
                        (sequence.sequence_id, first_bin, last_bin),
                    )
                )
        except sqlite3.DatabaseError as exc:
            raise ChainIndexCorruptionError(
                "chain index database lookup failed"
            ) from exc

        rows_by_bin: dict[int, list[tuple[int, int, int, int, int, int]]] = {}
        try:
            for (
                bin_id,
                record_id,
                target_start,
                target_end,
                block_id,
                block_offset,
                record_length,
            ) in rows:
                typed_row = (
                    _catalog_nonnegative_int(record_id),
                    _catalog_nonnegative_int(target_start),
                    _catalog_nonnegative_int(target_end),
                    _catalog_nonnegative_int(block_id),
                    _catalog_nonnegative_int(block_offset),
                    _catalog_positive_int(record_length),
                )
                rows_by_bin.setdefault(_catalog_nonnegative_int(bin_id), []).append(
                    typed_row
                )
        except (TypeError, ValueError) as exc:
            raise ChainIndexCorruptionError(
                "chain index database lookup returned malformed metadata"
            ) from exc

        authenticated_records: dict[int, tuple[int, int, int, int, int]] = {}
        for bin_id in range(first_bin, last_bin + 1):
            selected_rows = rows_by_bin.get(bin_id, [])
            expected = self._lookup_catalog.bins.get((sequence.sequence_id, bin_id))
            actual_digest = hashlib.sha256()
            try:
                for (
                    record_id,
                    target_start,
                    target_end,
                    block_id,
                    block_offset,
                    record_length,
                ) in selected_rows:
                    actual_digest.update(
                        _bin_integrity_payload(
                            sequence_id=sequence.sequence_id,
                            bin_id=bin_id,
                            record_id=record_id,
                            target_start=target_start,
                            target_end=target_end,
                            block_id=block_id,
                            block_offset=block_offset,
                            record_length=record_length,
                        )
                    )
            except (struct.error, OverflowError) as exc:
                raise ChainIndexCorruptionError(
                    "chain index database lookup returned invalid integer metadata"
                ) from exc

            if expected is None:
                if selected_rows:
                    raise ChainIndexCorruptionError(
                        "chain index bin lookup integrity mismatch"
                    )
            elif (
                len(selected_rows) != expected.count
                or actual_digest.hexdigest() != expected.sha256_hex
            ):
                raise ChainIndexCorruptionError(
                    "chain index bin lookup integrity mismatch"
                )

            for (
                record_id,
                target_start,
                target_end,
                block_id,
                block_offset,
                record_length,
            ) in selected_rows:
                previous = authenticated_records.get(record_id)
                current = (
                    target_start,
                    target_end,
                    block_id,
                    block_offset,
                    record_length,
                )
                if previous is not None and previous != current:
                    raise ChainIndexCorruptionError(
                        "chain index record metadata disagrees across bins"
                    )
                authenticated_records[record_id] = current

        selected: list[tuple[int, int, int, int, int, int]] = []
        for record_id, metadata in authenticated_records.items():
            target_start, target_end, block_id, block_offset, record_length = metadata
            if target_start < interval.end and target_end > interval.start:
                selected.append(
                    (
                        record_id,
                        target_start,
                        target_end,
                        block_id,
                        block_offset,
                        record_length,
                    )
                )
        selected.sort(key=lambda row: row[0])

        payload = bytearray()
        expected_geometry: list[tuple[int, int]] = []
        current_block_id: int | None = None
        current_block_payload = b""
        with self.block_store_path.open("rb") as block_store:
            for (
                _record_id,
                target_start,
                target_end,
                selected_block_id,
                block_offset,
                record_length,
            ) in selected:
                block = self._lookup_catalog.blocks_by_id.get(selected_block_id)
                if block is None:
                    raise ChainIndexCorruptionError(
                        "chain index record refers to an unknown compression block"
                    )
                if (
                    block_offset < 0
                    or record_length <= 0
                    or (block_offset + record_length > block.uncompressed_length)
                ):
                    raise ChainIndexCorruptionError(
                        "chain index record slice exceeds catalogued compression block"
                    )

                if selected_block_id != current_block_id:
                    block_store.seek(block.file_offset)
                    compressed = block_store.read(block.compressed_length)
                    if len(compressed) != block.compressed_length:
                        raise ChainIndexCorruptionError(
                            "chain index compressed block ended unexpectedly"
                        )
                    actual_sha256 = hashlib.sha256(compressed).hexdigest()
                    if actual_sha256 != block.compressed_sha256_hex:
                        raise ChainIndexCorruptionError(
                            "chain index compressed block SHA256 mismatch"
                        )
                    try:
                        current_block_payload = zlib.decompress(compressed)
                    except zlib.error as exc:
                        raise ChainIndexCorruptionError(
                            "chain index compressed block could not be decoded"
                        ) from exc
                    if len(current_block_payload) != block.uncompressed_length:
                        raise ChainIndexCorruptionError(
                            "chain index decoded block length does not match metadata"
                        )
                    current_block_id = selected_block_id

                record_bytes = current_block_payload[
                    block_offset : block_offset + record_length
                ]
                if len(record_bytes) != record_length:
                    raise ChainIndexCorruptionError(
                        "chain index record slice exceeds decoded compression block"
                    )
                payload.extend(record_bytes)
                expected_geometry.append((target_start, target_end))

        try:
            records = tuple(iter_chain_records(io.StringIO(payload.decode("ascii"))))
        except (UnicodeDecodeError, ValueError) as exc:
            raise ChainIndexCorruptionError(
                "chain index selected record payload is malformed"
            ) from exc
        if len(records) != len(selected):
            raise ChainIndexCorruptionError(
                "chain index selected record count does not match metadata"
            )
        for parsed_record, (target_start, target_end) in zip(
            records, expected_geometry
        ):
            if (
                parsed_record.target_name != interval.sequence_name
                or parsed_record.target_start != target_start
                or parsed_record.target_end != target_end
            ):
                raise ChainIndexCorruptionError(
                    "chain index selected record geometry does not match "
                    "lookup metadata"
                )
        return records


def chain_index_cache_path(
    cache_root: ResourcePath,
    source_chain_sha256_identifier: str,
) -> Path:
    """Return the deterministic cache directory for one exact source chain index."""

    sha256_hex = sha256_hex_from_identifier(source_chain_sha256_identifier)
    return (
        Path(cache_root)
        / "derived"
        / f"chain-index-v{_INDEX_SCHEMA_VERSION}"
        / "sha256"
        / sha256_hex[:2]
        / sha256_hex
    )


def build_cached_chain_index(
    cache_root: ResourcePath,
    resource: CachedResource,
    *,
    bin_width: int = DEFAULT_CHAIN_INDEX_BIN_WIDTH,
    block_size: int = DEFAULT_CHAIN_INDEX_BLOCK_SIZE,
    progress_callback: Callable[[int], None] | None = None,
) -> ChainIndexBuildResult:
    """Build the deterministic reusable index for one verified cached chain resource."""

    output_dir = chain_index_cache_path(cache_root, resource.sha256)
    return build_chain_index(
        resource.path,
        output_dir,
        source_chain_sha256_identifier=resource.sha256,
        source_chain_size_bytes=resource.size_bytes,
        bin_width=bin_width,
        block_size=block_size,
        progress_callback=progress_callback,
    )


def load_cached_chain_index(
    cache_root: ResourcePath,
    resource: CachedResource,
    *,
    verify_database: bool = False,
) -> ChainIndex | None:
    """Load the deterministic index for one cached chain resource when present."""

    index_dir = chain_index_cache_path(cache_root, resource.sha256)
    if not index_dir.is_dir():
        return None
    return load_chain_index(
        index_dir,
        expected_source_chain_sha256_identifier=resource.sha256,
        expected_source_chain_size_bytes=resource.size_bytes,
        verify_database=verify_database,
    )


def build_chain_index(
    chain_path: ResourcePath,
    output_dir: ResourcePath,
    *,
    source_chain_sha256_identifier: str,
    source_chain_size_bytes: int,
    bin_width: int = DEFAULT_CHAIN_INDEX_BIN_WIDTH,
    block_size: int = DEFAULT_CHAIN_INDEX_BLOCK_SIZE,
    progress_callback: Callable[[int], None] | None = None,
) -> ChainIndexBuildResult:
    """Build one atomic reusable index while verifying the exact source bytes.

    The destination must not already exist. The complete index is built in a private
    sibling directory and renamed into place only after source SHA-256 verification,
    database finalization, and manifest creation all succeed.
    """

    expected_sha256_hex = sha256_hex_from_identifier(source_chain_sha256_identifier)
    if source_chain_size_bytes < 0:
        raise ValueError("source chain size must not be negative")
    if bin_width <= 0:
        raise ValueError("chain index bin width must be positive")
    if block_size <= 0:
        raise ValueError("chain index block size must be positive")

    source_path = Path(chain_path)
    actual_size = source_path.stat().st_size
    if actual_size != source_chain_size_bytes:
        raise ResourceIdentityMismatchError(
            f"source chain size mismatch for {source_path}: "
            f"expected {source_chain_size_bytes}, got {actual_size}"
        )

    destination = Path(output_dir)
    if destination.exists():
        raise FileExistsError(f"chain index destination already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.tmp-", dir=destination.parent)
    )

    try:
        database = temporary / "index.sqlite3"
        blocks = temporary / "records.blocks"
        lookup_catalog = temporary / "lookup-catalog.json"
        manifest_path = temporary / "manifest.json"
        (
            record_count,
            membership_count,
            block_count,
            lookup_catalog_payload,
        ) = _build_index_files(
            source_path,
            database,
            blocks,
            expected_sha256_hex=expected_sha256_hex,
            bin_width=bin_width,
            block_size=block_size,
            progress_callback=progress_callback,
        )
        lookup_catalog.write_text(
            json.dumps(
                lookup_catalog_payload,
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        database_sha256_hex = compute_resource_checksum(
            database, ResourceChecksumAlgorithm.SHA256
        )
        lookup_catalog_sha256_hex = compute_resource_checksum(
            lookup_catalog, ResourceChecksumAlgorithm.SHA256
        )
        manifest = ChainIndexManifest(
            source_chain_sha256_identifier=source_chain_sha256_identifier,
            source_chain_size_bytes=source_chain_size_bytes,
            bin_width=bin_width,
            block_size=block_size,
            record_count=record_count,
            membership_count=membership_count,
            block_count=block_count,
            database_size_bytes=database.stat().st_size,
            block_store_size_bytes=blocks.stat().st_size,
            lookup_catalog_size_bytes=lookup_catalog.stat().st_size,
            database_sha256_identifier=f"sha256:{database_sha256_hex}",
            lookup_catalog_sha256_identifier=(f"sha256:{lookup_catalog_sha256_hex}"),
        )
        manifest_path.write_text(
            json.dumps(manifest.to_json_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.rename(destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    index = load_chain_index(
        destination,
        expected_source_chain_sha256_identifier=source_chain_sha256_identifier,
        verify_database=False,
    )
    return ChainIndexBuildResult(index=index, manifest=manifest)


def load_chain_index(
    index_dir: ResourcePath,
    *,
    expected_source_chain_sha256_identifier: str | None = None,
    expected_source_chain_size_bytes: int | None = None,
    verify_database: bool = True,
) -> ChainIndex:
    """Load and validate one existing reusable chain index.

    The compact lookup catalog is always checked for exact size and SHA-256 and is the
    per-query integrity authority for sequence metadata, bin membership/record
    locators, and compressed-block metadata. ``verify_database=True`` additionally
    performs the expensive whole-database SHA-256 check for deep maintenance
    verification.
    """

    root = Path(index_dir)
    manifest_path = root / "manifest.json"
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ChainIndexCorruptionError(
            "chain index manifest could not be read"
        ) from exc
    try:
        manifest = ChainIndexManifest.from_json_payload(payload)
    except (TypeError, ValueError) as exc:
        raise ChainIndexCorruptionError("chain index manifest is invalid") from exc

    if (
        expected_source_chain_sha256_identifier is not None
        and manifest.source_chain_sha256_identifier
        != expected_source_chain_sha256_identifier
    ):
        raise ChainIndexCorruptionError(
            "chain index source identity does not match the expected chain"
        )

    if (
        expected_source_chain_size_bytes is not None
        and manifest.source_chain_size_bytes != expected_source_chain_size_bytes
    ):
        raise ChainIndexCorruptionError(
            "chain index source size does not match the expected chain"
        )

    database = root / "index.sqlite3"
    blocks = root / "records.blocks"
    lookup_catalog_path = root / "lookup-catalog.json"
    if (
        not database.is_file()
        or database.stat().st_size != manifest.database_size_bytes
    ):
        raise ChainIndexCorruptionError(
            "chain index database is missing or has the wrong size"
        )
    if not blocks.is_file() or blocks.stat().st_size != manifest.block_store_size_bytes:
        raise ChainIndexCorruptionError(
            "chain index block store is missing or has the wrong size"
        )
    if (
        not lookup_catalog_path.is_file()
        or lookup_catalog_path.stat().st_size != manifest.lookup_catalog_size_bytes
    ):
        raise ChainIndexCorruptionError(
            "chain index lookup catalog is missing or has the wrong size"
        )

    expected_lookup_catalog_sha256_hex = sha256_hex_from_identifier(
        manifest.lookup_catalog_sha256_identifier
    )
    actual_lookup_catalog_sha256_hex = compute_resource_checksum(
        lookup_catalog_path, ResourceChecksumAlgorithm.SHA256
    )
    if actual_lookup_catalog_sha256_hex != expected_lookup_catalog_sha256_hex:
        raise ChainIndexCorruptionError("chain index lookup catalog SHA256 mismatch")
    lookup_catalog = _load_lookup_catalog(lookup_catalog_path, manifest)

    if verify_database:
        expected_database_sha256_hex = sha256_hex_from_identifier(
            manifest.database_sha256_identifier
        )
        actual_database_sha256_hex = compute_resource_checksum(
            database, ResourceChecksumAlgorithm.SHA256
        )
        if actual_database_sha256_hex != expected_database_sha256_hex:
            raise ChainIndexCorruptionError("chain index database SHA256 mismatch")

    return ChainIndex(
        root=root,
        manifest=manifest,
        _lookup_catalog=lookup_catalog,
    )


def _build_index_files(
    source_path: Path,
    database: Path,
    blocks: Path,
    *,
    expected_sha256_hex: str,
    bin_width: int,
    block_size: int,
    progress_callback: Callable[[int], None] | None,
) -> tuple[int, int, int, dict[str, object]]:
    sequence_ids: dict[str, int] = {}
    sequence_min_sizes: dict[str, int] = {}
    chain_rows: list[tuple[int, int, int, int, int, int, int]] = []
    membership_rows: list[tuple[int, int, int]] = []
    # One accumulator is retained per populated genomic bin, not per membership row.
    # Memory therefore scales with covered 65,536-bp bins rather than chain count.
    bin_digests: dict[tuple[int, int], _BinDigestAccumulator] = {}
    block_catalog: list[_LookupBlock] = []
    block_payload = bytearray()
    block_id = 0
    block_count = 0
    record_count = 0
    membership_count = 0

    with (
        closing(sqlite3.connect(database)) as connection,
        blocks.open("wb") as block_store,
        open_text_resource(
            source_path,
            expected_sha256_hex=expected_sha256_hex,
            progress_callback=progress_callback,
        ) as lines,
    ):
        _configure_database(connection)
        for record_id, chain in enumerate(iter_chain_records(lines)):
            sequence_id = sequence_ids.get(chain.target_name)
            if sequence_id is None:
                sequence_id = len(sequence_ids)
                sequence_ids[chain.target_name] = sequence_id
                sequence_min_sizes[chain.target_name] = chain.target_size
                connection.execute(
                    "INSERT INTO sequences VALUES (?, ?, ?)",
                    (sequence_id, chain.target_name, chain.target_size),
                )
            elif chain.target_size < sequence_min_sizes[chain.target_name]:
                # The full traversal checks bounds against every chain with this target
                # name before testing overlap. Using the minimum size preserves that
                # fail-fast behavior even if a malformed resource disagrees internally.
                sequence_min_sizes[chain.target_name] = chain.target_size
                connection.execute(
                    "UPDATE sequences SET minimum_target_size = ? "
                    "WHERE sequence_id = ?",
                    (chain.target_size, sequence_id),
                )

            payload = _format_chain_record(chain)
            if block_payload and len(block_payload) + len(payload) > block_size:
                block_catalog.append(
                    _write_compressed_block(
                        block_store,
                        connection,
                        block_id=block_id,
                        payload=bytes(block_payload),
                    )
                )
                block_count += 1
                block_id += 1
                block_payload.clear()

            block_offset = len(block_payload)
            block_payload.extend(payload)
            chain_rows.append(
                (
                    record_id,
                    sequence_id,
                    chain.target_start,
                    chain.target_end,
                    block_id,
                    block_offset,
                    len(payload),
                )
            )

            first_bin = chain.target_start // bin_width
            last_bin = (chain.target_end - 1) // bin_width
            for bin_id_value in range(first_bin, last_bin + 1):
                membership_rows.append((sequence_id, bin_id_value, record_id))
                bin_key = (sequence_id, bin_id_value)
                accumulator = bin_digests.get(bin_key)
                if accumulator is None:
                    accumulator = _BinDigestAccumulator()
                    bin_digests[bin_key] = accumulator
                accumulator.update(
                    sequence_id=sequence_id,
                    bin_id=bin_id_value,
                    record_id=record_id,
                    target_start=chain.target_start,
                    target_end=chain.target_end,
                    block_id=block_id,
                    block_offset=block_offset,
                    record_length=len(payload),
                )
                membership_count += 1

            record_count += 1
            if (
                len(chain_rows) >= _CHAIN_BATCH_SIZE
                or len(membership_rows) >= _MEMBERSHIP_BATCH_SIZE
            ):
                _flush_rows(connection, chain_rows, membership_rows)

        if block_payload:
            block_catalog.append(
                _write_compressed_block(
                    block_store,
                    connection,
                    block_id=block_id,
                    payload=bytes(block_payload),
                )
            )
            block_count += 1
        _flush_rows(connection, chain_rows, membership_rows)
        connection.commit()
        connection.execute("ANALYZE")
        connection.commit()

    sequence_catalog = [
        [sequence_id, name, sequence_min_sizes[name]]
        for name, sequence_id in sorted(sequence_ids.items(), key=lambda item: item[1])
    ]
    block_catalog_payload = [
        [
            block.block_id,
            block.file_offset,
            block.compressed_length,
            block.uncompressed_length,
            block.compressed_sha256_hex,
        ]
        for block in block_catalog
    ]
    bin_catalog = [
        [sequence_id, bin_id, accumulator.count, accumulator.digest.hexdigest()]
        for (sequence_id, bin_id), accumulator in sorted(bin_digests.items())
    ]
    lookup_catalog_payload: dict[str, object] = {
        "format": _LOOKUP_CATALOG_FORMAT,
        "schema_version": _LOOKUP_CATALOG_SCHEMA_VERSION,
        "sequences": sequence_catalog,
        "blocks": block_catalog_payload,
        "bins": bin_catalog,
    }
    return record_count, membership_count, block_count, lookup_catalog_payload


def _configure_database(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("PRAGMA temp_store = MEMORY")
    connection.executescript(
        """
        CREATE TABLE sequences (
            sequence_id INTEGER PRIMARY KEY,
            name TEXT NOT NULL UNIQUE,
            minimum_target_size INTEGER NOT NULL
        );
        CREATE TABLE chains (
            record_id INTEGER PRIMARY KEY,
            sequence_id INTEGER NOT NULL,
            target_start INTEGER NOT NULL,
            target_end INTEGER NOT NULL,
            block_id INTEGER NOT NULL,
            block_offset INTEGER NOT NULL,
            record_length INTEGER NOT NULL
        );
        CREATE TABLE blocks (
            block_id INTEGER PRIMARY KEY,
            file_offset INTEGER NOT NULL,
            compressed_length INTEGER NOT NULL,
            uncompressed_length INTEGER NOT NULL,
            compressed_sha256_hex TEXT NOT NULL
        );
        CREATE TABLE bin_memberships (
            sequence_id INTEGER NOT NULL,
            bin_id INTEGER NOT NULL,
            record_id INTEGER NOT NULL,
            PRIMARY KEY (sequence_id, bin_id, record_id)
        ) WITHOUT ROWID;
        """
    )


def _write_compressed_block(
    block_store: io.BufferedWriter,
    connection: sqlite3.Connection,
    *,
    block_id: int,
    payload: bytes,
) -> _LookupBlock:
    compressed = zlib.compress(payload)
    file_offset = block_store.tell()
    compressed_sha256_hex = hashlib.sha256(compressed).hexdigest()
    block_store.write(compressed)
    connection.execute(
        "INSERT INTO blocks VALUES (?, ?, ?, ?, ?)",
        (
            block_id,
            file_offset,
            len(compressed),
            len(payload),
            compressed_sha256_hex,
        ),
    )
    return _LookupBlock(
        block_id=block_id,
        file_offset=file_offset,
        compressed_length=len(compressed),
        uncompressed_length=len(payload),
        compressed_sha256_hex=compressed_sha256_hex,
    )


def _flush_rows(
    connection: sqlite3.Connection,
    chain_rows: list[tuple[int, int, int, int, int, int, int]],
    membership_rows: list[tuple[int, int, int]],
) -> None:
    if chain_rows:
        connection.executemany(
            "INSERT INTO chains VALUES (?, ?, ?, ?, ?, ?, ?)",
            chain_rows,
        )
        chain_rows.clear()
    if membership_rows:
        connection.executemany(
            "INSERT INTO bin_memberships VALUES (?, ?, ?)",
            membership_rows,
        )
        membership_rows.clear()


def _format_chain_record(record: ChainRecord) -> bytes:
    header = " ".join(
        (
            "chain",
            repr(record.score),
            record.target_name,
            str(record.target_size),
            record.target_strand.value,
            str(record.target_start),
            str(record.target_end),
            record.query_name,
            str(record.query_size),
            record.query_strand.value,
            str(record.query_start),
            str(record.query_end),
            str(record.chain_id),
        )
    )
    lines = [header]
    for block in record.blocks:
        if block.is_terminal:
            lines.append(str(block.size))
        else:
            target_gap, query_gap = block.gaps_after()
            lines.append(f"{block.size}\t{target_gap}\t{query_gap}")
    return ("\n".join(lines) + "\n\n").encode("ascii")


def _bin_integrity_payload(
    *,
    sequence_id: int,
    bin_id: int,
    record_id: int,
    target_start: int,
    target_end: int,
    block_id: int,
    block_offset: int,
    record_length: int,
) -> bytes:
    return struct.pack(
        ">8Q",
        sequence_id,
        bin_id,
        record_id,
        target_start,
        target_end,
        block_id,
        block_offset,
        record_length,
    )


def _load_lookup_catalog(
    path: Path,
    manifest: ChainIndexManifest,
) -> _ChainIndexLookupCatalog:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ChainIndexCorruptionError(
            "chain index lookup catalog could not be read"
        ) from exc
    if not isinstance(payload, dict):
        raise ChainIndexCorruptionError("chain index lookup catalog is invalid")
    if (
        payload.get("format") != _LOOKUP_CATALOG_FORMAT
        or payload.get("schema_version") != _LOOKUP_CATALOG_SCHEMA_VERSION
    ):
        raise ChainIndexCorruptionError("chain index lookup catalog is invalid")

    sequences_payload = payload.get("sequences")
    blocks_payload = payload.get("blocks")
    bins_payload = payload.get("bins")
    if not isinstance(sequences_payload, list) or not isinstance(blocks_payload, list):
        raise ChainIndexCorruptionError("chain index lookup catalog is invalid")
    if not isinstance(bins_payload, list):
        raise ChainIndexCorruptionError("chain index lookup catalog is invalid")
    if len(blocks_payload) != manifest.block_count:
        raise ChainIndexCorruptionError(
            "chain index lookup catalog block count does not match manifest"
        )

    sequences_by_name: dict[str, _LookupSequence] = {}
    sequence_ids: set[int] = set()
    try:
        for item in sequences_payload:
            if not isinstance(item, list):
                raise TypeError
            if len(item) != 3:
                raise ValueError
            sequence_id = _catalog_nonnegative_int(item[0])
            name = _catalog_nonempty_str(item[1])
            minimum_target_size = _catalog_nonnegative_int(item[2])
            if sequence_id in sequence_ids or name in sequences_by_name:
                raise ValueError
            sequence_ids.add(sequence_id)
            sequences_by_name[name] = _LookupSequence(
                sequence_id=sequence_id,
                name=name,
                minimum_target_size=minimum_target_size,
            )

        blocks_by_id: dict[int, _LookupBlock] = {}
        for item in blocks_payload:
            if not isinstance(item, list):
                raise TypeError
            if len(item) != 5:
                raise ValueError
            block_id = _catalog_nonnegative_int(item[0])
            file_offset = _catalog_nonnegative_int(item[1])
            compressed_length = _catalog_positive_int(item[2])
            uncompressed_length = _catalog_positive_int(item[3])
            compressed_sha256_hex = _catalog_sha256_hex(item[4])
            if block_id in blocks_by_id:
                raise ValueError
            if file_offset + compressed_length > manifest.block_store_size_bytes:
                raise ValueError
            blocks_by_id[block_id] = _LookupBlock(
                block_id=block_id,
                file_offset=file_offset,
                compressed_length=compressed_length,
                uncompressed_length=uncompressed_length,
                compressed_sha256_hex=compressed_sha256_hex,
            )

        bins: dict[tuple[int, int], _LookupBinIntegrity] = {}
        for item in bins_payload:
            if not isinstance(item, list):
                raise TypeError
            if len(item) != 4:
                raise ValueError
            sequence_id = _catalog_nonnegative_int(item[0])
            bin_id = _catalog_nonnegative_int(item[1])
            count = _catalog_positive_int(item[2])
            sha256_hex = _catalog_sha256_hex(item[3])
            key = (sequence_id, bin_id)
            if sequence_id not in sequence_ids or key in bins:
                raise ValueError
            bins[key] = _LookupBinIntegrity(count=count, sha256_hex=sha256_hex)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ChainIndexCorruptionError(
            "chain index lookup catalog is invalid"
        ) from exc

    if sum(entry.count for entry in bins.values()) != manifest.membership_count:
        raise ChainIndexCorruptionError(
            "chain index lookup catalog membership count does not match manifest"
        )

    return _ChainIndexLookupCatalog(
        sequences_by_name=MappingProxyType(sequences_by_name),
        blocks_by_id=MappingProxyType(blocks_by_id),
        bins=MappingProxyType(bins),
    )


def _catalog_nonnegative_int(value: object) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError
    if value < 0:
        raise ValueError
    return value


def _catalog_positive_int(value: object) -> int:
    result = _catalog_nonnegative_int(value)
    if result == 0:
        raise ValueError
    return result


def _catalog_nonempty_str(value: object) -> str:
    if not isinstance(value, str):
        raise TypeError
    if not value:
        raise ValueError
    return value


def _catalog_sha256_hex(value: object) -> str:
    result = _catalog_nonempty_str(value)
    if len(result) != 64 or any(
        character not in "0123456789abcdef" for character in result
    ):
        raise ValueError
    return result


def _required_str(payload: dict[Any, Any], key: str) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"chain index manifest field {key!r} must be a string")
    return value


def _required_int(payload: dict[Any, Any], key: str) -> int:
    value = payload.get(key)
    if not isinstance(value, int) or isinstance(value, bool):
        raise TypeError(f"chain index manifest field {key!r} must be an integer")
    return value
