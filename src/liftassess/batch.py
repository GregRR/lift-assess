"""Batch input records and cross-record target-projection relationships.

Batch relationships are derived above candidate-level scientific evidence.  They
compare exact mapped target segments from distinct input records and deliberately
never use a candidate's target bounding span as if it were continuous alignment.
"""

from dataclasses import dataclass
from enum import Enum
from itertools import combinations

from .models import AssemblyIdentifier, GenomicInterval, NormalizedCandidate


class BatchTargetRelationshipKind(str, Enum):
    """Factual relationship between projections from two distinct input records."""

    EXACT_TARGET_COLLISION = "EXACT_TARGET_COLLISION"
    OVERLAPPING_TARGET_PROJECTIONS = "OVERLAPPING_TARGET_PROJECTIONS"


@dataclass(frozen=True)
class BatchInputRecord:
    """One non-empty source interval from a batch input stream.

    ``record_id`` is a deterministic per-input identifier. ``label`` preserves an
    optional user-facing name without requiring labels to be unique.
    """

    record_id: str
    source_interval: GenomicInterval
    source_line_number: int
    label: str | None = None

    def __post_init__(self) -> None:
        if not self.record_id:
            raise ValueError("batch record_id must not be empty")
        if self.source_interval.length <= 0:
            raise ValueError("batch source intervals must span at least one base")
        if self.source_line_number < 1:
            raise ValueError("batch source_line_number must be at least 1")
        if self.label == "":
            raise ValueError("batch label must not be empty when provided")


@dataclass(frozen=True)
class BatchRecordAssessment:
    """Candidate inventory associated with one batch input record."""

    record: BatchInputRecord
    candidates: tuple[NormalizedCandidate, ...]

    def __post_init__(self) -> None:
        candidate_ids = [candidate.candidate_id for candidate in self.candidates]
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError(
                "batch candidate IDs must be unique within one input record"
            )

        source = self.record.source_interval
        for candidate in self.candidates:
            for segment in candidate.segments:
                if (
                    segment.source_interval.assembly != source.assembly
                    or segment.source_interval.sequence_name != source.sequence_name
                    or segment.source_interval.start < source.start
                    or segment.source_interval.end > source.end
                ):
                    raise ValueError(
                        "batch candidates must map source segments contained within "
                        "their input record"
                    )


@dataclass(frozen=True)
class BatchTargetRelationship:
    """One exact cross-record relationship between two candidate projections."""

    kind: BatchTargetRelationshipKind
    left_record_id: str
    left_candidate_id: str
    right_record_id: str
    right_candidate_id: str
    target_assembly: AssemblyIdentifier
    target_sequence_name: str
    overlap_intervals: tuple[GenomicInterval, ...]

    def __post_init__(self) -> None:
        if not self.left_record_id or not self.right_record_id:
            raise ValueError("batch relationship record IDs must not be empty")
        if self.left_record_id == self.right_record_id:
            raise ValueError("batch relationships require two distinct input records")
        if not self.left_candidate_id or not self.right_candidate_id:
            raise ValueError("batch relationship candidate IDs must not be empty")
        if not self.target_sequence_name:
            raise ValueError("batch relationship target sequence must not be empty")
        if not self.overlap_intervals:
            raise ValueError("batch relationships require mapped target overlap")
        for interval in self.overlap_intervals:
            if (
                interval.assembly != self.target_assembly
                or interval.sequence_name != self.target_sequence_name
                or interval.length <= 0
            ):
                raise ValueError(
                    "batch relationship overlap intervals must share one non-empty "
                    "target assembly and sequence"
                )


@dataclass(frozen=True)
class BatchRelationshipResult:
    """Deterministic cross-record relationships for one assessed batch."""

    relationships: tuple[BatchTargetRelationship, ...]


def build_batch_target_relationships(
    assessments: tuple[BatchRecordAssessment, ...],
) -> BatchRelationshipResult:
    """Compare exact candidate target coverage across distinct input records.

    An exact collision requires the same covered target bases after adjacent target
    segments are canonicalized.  A non-identical positive intersection is reported
    separately as an overlapping target projection.  Bounding spans are never used
    for either relationship.
    """

    _validate_batch_record_ids(assessments)
    relationships: list[BatchTargetRelationship] = []

    for left, right in combinations(assessments, 2):
        for left_candidate in left.candidates:
            left_coverage = _canonical_target_coverage(left_candidate)
            for right_candidate in right.candidates:
                if (
                    left_candidate.target_interval.assembly
                    != right_candidate.target_interval.assembly
                    or left_candidate.target_interval.sequence_name
                    != right_candidate.target_interval.sequence_name
                ):
                    continue

                right_coverage = _canonical_target_coverage(right_candidate)
                overlaps = _target_intersections(left_coverage, right_coverage)
                if not overlaps:
                    continue

                kind = (
                    BatchTargetRelationshipKind.EXACT_TARGET_COLLISION
                    if left_coverage == right_coverage
                    else BatchTargetRelationshipKind.OVERLAPPING_TARGET_PROJECTIONS
                )
                relationships.append(
                    BatchTargetRelationship(
                        kind=kind,
                        left_record_id=left.record.record_id,
                        left_candidate_id=left_candidate.candidate_id,
                        right_record_id=right.record.record_id,
                        right_candidate_id=right_candidate.candidate_id,
                        target_assembly=left_candidate.target_interval.assembly,
                        target_sequence_name=(
                            left_candidate.target_interval.sequence_name
                        ),
                        overlap_intervals=overlaps,
                    )
                )

    return BatchRelationshipResult(relationships=tuple(relationships))


def _validate_batch_record_ids(
    assessments: tuple[BatchRecordAssessment, ...],
) -> None:
    record_ids = [assessment.record.record_id for assessment in assessments]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("batch record IDs must be unique within one assessment")


def _canonical_target_coverage(
    candidate: NormalizedCandidate,
) -> tuple[GenomicInterval, ...]:
    intervals = sorted(
        (segment.target_interval for segment in candidate.segments),
        key=lambda interval: (interval.start, interval.end),
    )
    merged: list[GenomicInterval] = []
    for interval in intervals:
        if merged and merged[-1].end == interval.start:
            previous = merged[-1]
            merged[-1] = GenomicInterval(
                assembly=previous.assembly,
                sequence_name=previous.sequence_name,
                start=previous.start,
                end=interval.end,
            )
        else:
            merged.append(interval)
    return tuple(merged)


def _target_intersections(
    left: tuple[GenomicInterval, ...],
    right: tuple[GenomicInterval, ...],
) -> tuple[GenomicInterval, ...]:
    intersections: list[GenomicInterval] = []
    left_index = 0
    right_index = 0

    while left_index < len(left) and right_index < len(right):
        left_interval = left[left_index]
        right_interval = right[right_index]
        overlap_start = max(left_interval.start, right_interval.start)
        overlap_end = min(left_interval.end, right_interval.end)
        if overlap_start < overlap_end:
            intersections.append(
                GenomicInterval(
                    assembly=left_interval.assembly,
                    sequence_name=left_interval.sequence_name,
                    start=overlap_start,
                    end=overlap_end,
                )
            )

        if left_interval.end <= right_interval.end:
            left_index += 1
        else:
            right_index += 1

    return tuple(intersections)
