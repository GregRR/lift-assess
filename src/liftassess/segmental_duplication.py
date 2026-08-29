"""Typed UCSC segmental-duplication context.

The UCSC ``genomicSuperDups`` table is treated as descriptive contextual evidence,
not as a mapping-quality score or a biological-correctness verdict.  Rows identify
putative paired genomic duplications published for one UCSC assembly.  liftAssess
keeps the exact assembly, paired intervals, orientation, provider UID, alignment
size, fraction matching, and content-addressed provenance needed to report overlaps
without collapsing them into a generic "difficult region" flag.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType

from .models import (
    AssemblyIdentifier,
    GenomicInterval,
    NormalizedCandidate,
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
from .resources import UCSCSegmentalDuplicationResource

_GENOMIC_SUPER_DUPS_FIELD_COUNT = 30


class UCSCSegmentalDuplicationCatalogError(RuntimeError):
    """Cached UCSC segmental-duplication bytes could not be verified or parsed."""


class SegmentalDuplicationCheckState(str, Enum):
    """Availability state for one assembly side of the typed context check."""

    ASSESSED = "ASSESSED"
    UNAVAILABLE = "UNAVAILABLE"
    NO_TARGET_PROJECTIONS = "NO_TARGET_PROJECTIONS"


@dataclass(frozen=True)
class UCSCSegmentalDuplicationRecord:
    """One exact row from UCSC's ``genomicSuperDups`` table."""

    interval: GenomicInterval
    paired_interval: GenomicInterval
    strand: str
    uid: int
    aligned_bases: int
    fraction_matching_bases: float

    def __post_init__(self) -> None:
        if self.interval.assembly != self.paired_interval.assembly:
            raise ValueError("segmental-duplication pair must use one assembly")
        if self.interval.length <= 0 or self.paired_interval.length <= 0:
            raise ValueError("segmental-duplication intervals must be non-empty")
        if self.strand not in {"+", "-"}:
            raise ValueError("segmental-duplication strand must be '+' or '-'")
        if self.uid < 0:
            raise ValueError("segmental-duplication UID cannot be negative")
        if self.aligned_bases <= 0:
            raise ValueError("segmental-duplication aligned bases must be positive")
        if not 0.0 <= self.fraction_matching_bases <= 1.0:
            raise ValueError(
                "segmental-duplication matching fraction must be between 0 and 1"
            )


@dataclass(frozen=True)
class UCSCSegmentalDuplicationCatalog:
    """One assembly's parsed UCSC segmental-duplication rows and provenance."""

    assembly: AssemblyIdentifier
    records: tuple[UCSCSegmentalDuplicationRecord, ...]
    provenance: ProvenanceSource
    _records_by_sequence: Mapping[str, tuple[UCSCSegmentalDuplicationRecord, ...]] = (
        field(init=False, repr=False, compare=False)
    )

    def __post_init__(self) -> None:
        by_sequence: dict[str, list[UCSCSegmentalDuplicationRecord]] = {}
        for record in self.records:
            if record.interval.assembly != self.assembly:
                raise ValueError(
                    "segmental-duplication catalog rows must match catalog assembly"
                )
            by_sequence.setdefault(record.interval.sequence_name, []).append(record)

        frozen: dict[str, tuple[UCSCSegmentalDuplicationRecord, ...]] = {}
        for sequence_name, records in by_sequence.items():
            records.sort(
                key=lambda item: (item.interval.start, item.interval.end, item.uid)
            )
            frozen[sequence_name] = tuple(records)
        object.__setattr__(self, "_records_by_sequence", MappingProxyType(frozen))

    def overlapping(
        self, interval: GenomicInterval
    ) -> tuple[UCSCSegmentalDuplicationOverlap, ...]:
        """Return exact table rows whose reference interval overlaps ``interval``."""

        if interval.assembly != self.assembly:
            raise ValueError(
                "segmental-duplication query interval must match catalog assembly"
            )
        records = self._records_by_sequence.get(interval.sequence_name, ())
        overlaps: list[UCSCSegmentalDuplicationOverlap] = []
        for record in records:
            if record.interval.start >= interval.end:
                break
            if record.interval.end <= interval.start:
                continue
            overlap_start = max(interval.start, record.interval.start)
            overlap_end = min(interval.end, record.interval.end)
            overlaps.append(
                UCSCSegmentalDuplicationOverlap(
                    record=record,
                    overlap_interval=GenomicInterval(
                        assembly=self.assembly,
                        sequence_name=interval.sequence_name,
                        start=overlap_start,
                        end=overlap_end,
                    ),
                )
            )
        return tuple(overlaps)


@dataclass(frozen=True)
class UCSCSegmentalDuplicationOverlap:
    """Exact overlap between a queried interval and one UCSC duplication row."""

    record: UCSCSegmentalDuplicationRecord
    overlap_interval: GenomicInterval

    def __post_init__(self) -> None:
        if self.overlap_interval.assembly != self.record.interval.assembly:
            raise ValueError("segmental-duplication overlap assembly must match row")
        if self.overlap_interval.sequence_name != self.record.interval.sequence_name:
            raise ValueError("segmental-duplication overlap sequence must match row")
        if (
            self.overlap_interval.start < self.record.interval.start
            or self.overlap_interval.end > self.record.interval.end
        ):
            raise ValueError(
                "segmental-duplication overlap must lie within row interval"
            )


@dataclass(frozen=True)
class CandidateSegmentalDuplicationOverlap:
    """One candidate's exact mapped target segments overlapping one duplication row."""

    candidate_id: str
    record: UCSCSegmentalDuplicationRecord
    overlap_intervals: tuple[GenomicInterval, ...]

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("candidate segmental-duplication overlap requires an ID")
        if not self.overlap_intervals:
            raise ValueError(
                "candidate segmental-duplication overlap requires at least one interval"
            )
        previous_end: int | None = None
        for interval in self.overlap_intervals:
            if interval.assembly != self.record.interval.assembly:
                raise ValueError(
                    "candidate overlap assembly must match duplication row"
                )
            if interval.sequence_name != self.record.interval.sequence_name:
                raise ValueError(
                    "candidate overlap sequence must match duplication row"
                )
            if (
                interval.start < self.record.interval.start
                or interval.end > self.record.interval.end
            ):
                raise ValueError("candidate overlap must lie within duplication row")
            if previous_end is not None and interval.start < previous_end:
                raise ValueError(
                    "candidate overlap intervals must be ordered/non-overlapping"
                )
            previous_end = interval.end


@dataclass(frozen=True)
class UCSCSegmentalDuplicationContextResult:
    """Typed source/target overlap observations for one mapping result."""

    source_state: SegmentalDuplicationCheckState
    target_state: SegmentalDuplicationCheckState
    source_overlaps: tuple[UCSCSegmentalDuplicationOverlap, ...] = ()
    target_overlaps: tuple[CandidateSegmentalDuplicationOverlap, ...] = ()
    source_provenance: ProvenanceSource | None = None
    target_provenance: ProvenanceSource | None = None

    def __post_init__(self) -> None:
        if self.source_state is SegmentalDuplicationCheckState.NO_TARGET_PROJECTIONS:
            raise ValueError("source segmental-duplication state cannot be no-target")
        if self.source_state is SegmentalDuplicationCheckState.ASSESSED:
            if self.source_provenance is None:
                raise ValueError("assessed source context requires provenance")
        elif self.source_overlaps or self.source_provenance is not None:
            raise ValueError("unavailable source context cannot carry observations")

        if self.target_state is SegmentalDuplicationCheckState.ASSESSED:
            if self.target_provenance is None:
                raise ValueError("assessed target context requires provenance")
        elif self.target_state is SegmentalDuplicationCheckState.NO_TARGET_PROJECTIONS:
            if self.target_overlaps or self.target_provenance is not None:
                raise ValueError("no-target context cannot carry target observations")
        elif self.target_overlaps or self.target_provenance is not None:
            raise ValueError("unavailable target context cannot carry observations")


def ucsc_segmental_duplication_table_url(db: str) -> str:
    """Return the canonical UCSC table-dump URL used as a cache identity.

    This function does not establish that the table exists. Online acquisition must use
    a separately discovered :class:`UCSCSegmentalDuplicationResource`.
    """

    allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-"
    if not db or any(character not in allowed for character in db):
        raise ValueError(f"invalid UCSC database name: {db!r}")
    return (
        "https://hgdownload.soe.ucsc.edu/goldenPath/"
        f"{db}/database/genomicSuperDups.txt.gz"
    )


def load_cached_ucsc_segmental_duplication_resource(
    cache_root: ResourcePath,
    db: str,
) -> CachedResource | None:
    """Load and verify an already-cached segmental-duplication table without network."""

    return load_cached_ucsc_resource(
        cache_root,
        ucsc_segmental_duplication_table_url(db),
    )


def acquire_ucsc_segmental_duplication_resource(
    resource: UCSCSegmentalDuplicationResource,
    cache_root: ResourcePath,
    *,
    refresh: bool = False,
) -> CachedResource:
    """Acquire a discovered freely-usable UCSC ``genomicSuperDups`` table."""

    return acquire_ucsc_resource(
        resource.url,
        cache_root,
        terms_acknowledged=False,
        refresh=refresh,
    )


def build_cached_ucsc_segmental_duplication_catalog(
    assembly: AssemblyIdentifier,
    resource: CachedResource,
) -> UCSCSegmentalDuplicationCatalog:
    """Parse one verified cached UCSC table into an in-memory interval catalog."""

    try:
        expected_sha256 = sha256_hex_from_identifier(resource.sha256)
        provenance = _segmental_duplication_provenance(assembly, resource)
        with open_text_resource(
            resource.path,
            expected_sha256_hex=expected_sha256,
        ) as lines:
            records = tuple(iter_ucsc_genomic_super_dups(lines, assembly=assembly))
        return UCSCSegmentalDuplicationCatalog(
            assembly=assembly,
            records=records,
            provenance=provenance,
        )
    except (OSError, ValueError) as exc:
        raise UCSCSegmentalDuplicationCatalogError(
            "cached UCSC segmental-duplication table could not be verified or parsed"
        ) from exc


def iter_ucsc_genomic_super_dups(
    lines: Iterable[str],
    *,
    assembly: AssemblyIdentifier,
) -> Iterable[UCSCSegmentalDuplicationRecord]:
    """Parse UCSC ``genomicSuperDups`` tab-separated rows.

    Coordinates are already 0-based, half-open in the UCSC table and therefore match
    liftAssess's canonical internal interval convention exactly.
    """

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != _GENOMIC_SUPER_DUPS_FIELD_COUNT:
            raise ValueError(
                "UCSC genomicSuperDups row must contain exactly "
                f"{_GENOMIC_SUPER_DUPS_FIELD_COUNT} tab-separated fields at line "
                f"{line_number}; found {len(fields)}"
            )

        chrom = fields[1]
        chrom_start = _parse_nonnegative_int(fields[2], line_number, "chromStart")
        chrom_end = _parse_nonnegative_int(fields[3], line_number, "chromEnd")
        strand = fields[6]
        other_chrom = fields[7]
        other_start = _parse_nonnegative_int(fields[8], line_number, "otherStart")
        other_end = _parse_nonnegative_int(fields[9], line_number, "otherEnd")
        uid = _parse_nonnegative_int(fields[11], line_number, "uid")
        aligned_bases = _parse_positive_int(fields[21], line_number, "alignB")
        fraction_matching = _parse_fraction(fields[26], line_number, "fracMatch")

        if chrom_end <= chrom_start:
            raise ValueError(
                "UCSC genomicSuperDups chromEnd must exceed chromStart at line "
                f"{line_number}"
            )
        if other_end <= other_start:
            raise ValueError(
                "UCSC genomicSuperDups otherEnd must exceed otherStart at line "
                f"{line_number}"
            )

        yield UCSCSegmentalDuplicationRecord(
            interval=GenomicInterval(
                assembly=assembly,
                sequence_name=chrom,
                start=chrom_start,
                end=chrom_end,
            ),
            paired_interval=GenomicInterval(
                assembly=assembly,
                sequence_name=other_chrom,
                start=other_start,
                end=other_end,
            ),
            strand=strand,
            uid=uid,
            aligned_bases=aligned_bases,
            fraction_matching_bases=fraction_matching,
        )


def build_ucsc_segmental_duplication_context(
    source_interval: GenomicInterval,
    candidates: tuple[NormalizedCandidate, ...],
    *,
    source_catalog: UCSCSegmentalDuplicationCatalog | None,
    target_catalog: UCSCSegmentalDuplicationCatalog | None,
    source_unavailable: bool,
    target_unavailable: bool,
) -> UCSCSegmentalDuplicationContextResult:
    """Build typed source/target overlap facts without changing mapping interpretation."""

    if source_catalog is not None and source_unavailable:
        raise ValueError("source context catalog and unavailable state are exclusive")
    if target_catalog is not None and target_unavailable:
        raise ValueError("target context catalog and unavailable state are exclusive")

    if source_catalog is None:
        if not source_unavailable:
            raise ValueError(
                "source segmental-duplication context state is unspecified"
            )
        source_state = SegmentalDuplicationCheckState.UNAVAILABLE
        source_overlaps: tuple[UCSCSegmentalDuplicationOverlap, ...] = ()
        source_provenance = None
    else:
        if source_catalog.assembly != source_interval.assembly:
            raise ValueError("source context catalog assembly must match source query")
        source_state = SegmentalDuplicationCheckState.ASSESSED
        source_overlaps = source_catalog.overlapping(source_interval)
        source_provenance = source_catalog.provenance

    if not candidates:
        if target_catalog is not None or target_unavailable:
            # Target context is not needed when no target projection exists. Do not
            # imply resource availability or failure for a side that was not queried.
            target_catalog = None
            target_unavailable = False
        target_state = SegmentalDuplicationCheckState.NO_TARGET_PROJECTIONS
        target_overlaps: tuple[CandidateSegmentalDuplicationOverlap, ...] = ()
        target_provenance = None
    elif target_catalog is None:
        if not target_unavailable:
            raise ValueError(
                "target segmental-duplication context state is unspecified"
            )
        target_state = SegmentalDuplicationCheckState.UNAVAILABLE
        target_overlaps = ()
        target_provenance = None
    else:
        target_assembly = candidates[0].target_interval.assembly
        if target_catalog.assembly != target_assembly:
            raise ValueError("target context catalog assembly must match candidates")
        if any(
            candidate.target_interval.assembly != target_assembly
            for candidate in candidates
        ):
            raise ValueError("all candidates must share one target assembly")
        target_state = SegmentalDuplicationCheckState.ASSESSED
        target_overlaps = _candidate_target_overlaps(candidates, target_catalog)
        target_provenance = target_catalog.provenance

    return UCSCSegmentalDuplicationContextResult(
        source_state=source_state,
        target_state=target_state,
        source_overlaps=source_overlaps,
        target_overlaps=target_overlaps,
        source_provenance=source_provenance,
        target_provenance=target_provenance,
    )


def _candidate_target_overlaps(
    candidates: tuple[NormalizedCandidate, ...],
    catalog: UCSCSegmentalDuplicationCatalog,
) -> tuple[CandidateSegmentalDuplicationOverlap, ...]:
    observations: list[CandidateSegmentalDuplicationOverlap] = []
    for candidate in candidates:
        by_record: dict[UCSCSegmentalDuplicationRecord, list[GenomicInterval]] = {}
        for segment in candidate.segments:
            for overlap in catalog.overlapping(segment.target_interval):
                by_record.setdefault(overlap.record, []).append(
                    overlap.overlap_interval
                )
        for record, intervals in by_record.items():
            intervals.sort(key=lambda item: (item.start, item.end))
            observations.append(
                CandidateSegmentalDuplicationOverlap(
                    candidate_id=candidate.candidate_id,
                    record=record,
                    overlap_intervals=tuple(intervals),
                )
            )
    observations.sort(
        key=lambda item: (
            item.candidate_id,
            item.record.interval.sequence_name,
            item.record.interval.start,
            item.record.interval.end,
            item.record.uid,
        )
    )
    return tuple(observations)


def _segmental_duplication_provenance(
    assembly: AssemblyIdentifier,
    resource: CachedResource,
) -> ProvenanceSource:
    identifier = ProvenanceIdentifier(
        kind=ProvenanceIdentifierKind.SHA256,
        value=resource.sha256,
    )
    return ProvenanceSource(
        source_id=f"file:{identifier.value}",
        label=f"UCSC {assembly.name} genomicSuperDups segmental-duplication table",
        identifiers=(identifier,),
    )


def _parse_nonnegative_int(value: str, line_number: int, field_name: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(
            f"invalid UCSC genomicSuperDups {field_name} at line {line_number}: "
            f"{value!r}"
        ) from exc
    if parsed < 0:
        raise ValueError(
            f"UCSC genomicSuperDups {field_name} cannot be negative at line "
            f"{line_number}"
        )
    return parsed


def _parse_positive_int(value: str, line_number: int, field_name: str) -> int:
    parsed = _parse_nonnegative_int(value, line_number, field_name)
    if parsed == 0:
        raise ValueError(
            f"UCSC genomicSuperDups {field_name} must be positive at line {line_number}"
        )
    return parsed


def _parse_fraction(value: str, line_number: int, field_name: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(
            f"invalid UCSC genomicSuperDups {field_name} at line {line_number}: "
            f"{value!r}"
        ) from exc
    if not 0.0 <= parsed <= 1.0:
        raise ValueError(
            f"UCSC genomicSuperDups {field_name} must be between 0 and 1 at line "
            f"{line_number}"
        )
    return parsed
