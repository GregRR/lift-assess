"""Derived factual result profile for liftAssess.

This module is the boundary between scientific candidate/evidence data and result
rendering.  It derives deterministic factual states from already-computed mapping
geometry and evidence without assigning an aggregate verdict, confidence score, or
biological truth claim.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import assert_never

from .models import (
    ChainGapSummary,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    NormalizedCandidate,
    ReciprocalBestMembershipSummary,
)


class InputValidityState(str, Enum):
    """Preflight state available to the current result profile."""

    NOT_ASSESSED = "NOT_ASSESSED"


class ProjectionCountState(str, Enum):
    """How many chain-derived candidates were generated."""

    NONE = "NONE"
    ONE = "ONE"
    MULTIPLE = "MULTIPLE"


class SourceCoverageState(str, Enum):
    """Source coverage of one candidate, or maximum candidate-level coverage."""

    NONE = "NONE"
    COMPLETE = "COMPLETE"
    PARTIAL = "PARTIAL"


class OrientationState(str, Enum):
    """Orientation relationship across the current projection set."""

    NONE = "NONE"
    SAME = "SAME"
    REVERSE = "REVERSE"
    MIXED = "MIXED"


class TargetRoleState(str, Enum):
    """Target-sequence metadata state reserved for the metadata milestone."""

    NOT_ASSESSED = "NOT_ASSESSED"


class ReverseResultState(str, Enum):
    """Actual reverse-mapping state, distinct from reciprocal-best membership."""

    NOT_RUN = "NOT_RUN"


class QueryContextState(str, Enum):
    """Point/neighborhood comparison state reserved for a later milestone."""

    NOT_RUN = "NOT_RUN"


class ComparativeRelationshipState(str, Enum):
    """Cross-resource comparative synthesis state for the current slice."""

    NOT_ASSESSED = "NOT_ASSESSED"


class BatchRelationshipState(str, Enum):
    """Cross-record relationship state for a single-record result."""

    NOT_CHECKED = "NOT_CHECKED"


class ExternalContextState(str, Enum):
    """Typed external/context evidence state for the current slice."""

    NOT_ASSESSED = "NOT_ASSESSED"


class FactualHeadline(str, Enum):
    """Deterministic factual mapping headline tokens."""

    NO_CHAIN_PROJECTION = "NO_CHAIN_PROJECTION"
    ONE_COMPLETE_CHAIN_PROJECTION = "ONE_COMPLETE_CHAIN_PROJECTION"
    PARTIAL_SOURCE_COVERAGE = "PARTIAL_SOURCE_COVERAGE"
    PARTIAL_AND_FRAGMENTED_PROJECTION = "PARTIAL_AND_FRAGMENTED_PROJECTION"
    COMPLETE_BUT_DISCONTINUOUS_PROJECTION = "COMPLETE_BUT_DISCONTINUOUS_PROJECTION"
    MULTIPLE_CHAIN_PROJECTIONS = "MULTIPLE_CHAIN_PROJECTIONS"
    SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS = (
        "SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS"
    )


@dataclass(frozen=True)
class CandidateResultProfile:
    """Derived factual geometry for one normalized candidate."""

    candidate_id: str
    coverage_state: SourceCoverageState
    covered_source_bases: int
    source_bases: int
    uncovered_source_intervals: tuple[GenomicInterval, ...]
    exact_mapped_segment_count: int
    geometric_segment_count: int
    fragmented: bool
    target_discontinuous: bool
    target_bounding_span: GenomicInterval
    source_gap_intervals: tuple[GenomicInterval, ...]
    target_gap_intervals: tuple[GenomicInterval, ...]
    largest_uncovered_source_span_bases: int
    largest_source_gap_bases: int
    largest_target_gap_bases: int
    orientation: MappingOrientation


@dataclass(frozen=True)
class ResultScopeProfile:
    """Explicit boundaries for dimensions not yet assessed in this tranche."""

    target_role: TargetRoleState = TargetRoleState.NOT_ASSESSED
    reverse_result: ReverseResultState = ReverseResultState.NOT_RUN
    query_context: QueryContextState = QueryContextState.NOT_RUN
    comparative_relationship: ComparativeRelationshipState = (
        ComparativeRelationshipState.NOT_ASSESSED
    )
    batch_relationship: BatchRelationshipState = BatchRelationshipState.NOT_CHECKED
    external_context: ExternalContextState = ExternalContextState.NOT_ASSESSED
    named_variant_identity_assessed: bool = False
    gene_transcript_identity_assessed: bool = False
    downstream_workflow_assessed: bool = False


@dataclass(frozen=True)
class ResultProfile:
    """Orthogonal factual result profile derived from current scientific facts."""

    source_interval: GenomicInterval
    input_validity: InputValidityState
    projection_count: ProjectionCountState
    source_coverage: SourceCoverageState
    orientation: OrientationState
    maximum_candidate_covered_source_bases: int
    source_bases: int
    maximum_coverage_candidate_ids: tuple[str, ...]
    union_covered_source_bases: int
    candidate_profiles: tuple[CandidateResultProfile, ...]
    headline: FactualHeadline
    interpretation: str
    evidence_tier: EvidenceAvailabilityTier
    consumed_resource_roles: tuple[str, ...]
    scope: ResultScopeProfile = field(default_factory=ResultScopeProfile)


def build_result_profile(
    source_interval: GenomicInterval,
    candidates: tuple[NormalizedCandidate, ...],
    *,
    evidence_tier: EvidenceAvailabilityTier,
    consumed_resource_roles: tuple[str, ...] = (),
) -> ResultProfile:
    """Derive one deterministic factual profile from normalized scientific data."""

    if source_interval.length <= 0:
        raise ValueError("result profile requires a non-empty source interval")

    _validate_candidate_ids(candidates)
    _validate_distinct_candidate_geometries(candidates)
    candidate_profiles = tuple(
        _candidate_result_profile(
            source_interval,
            candidate,
            evidence_tier=evidence_tier,
        )
        for candidate in candidates
    )

    projection_count = _projection_count(len(candidates))
    if not candidate_profiles:
        maximum_covered = 0
        maximum_ids: tuple[str, ...] = ()
        coverage_state = SourceCoverageState.NONE
        union_covered = 0
    else:
        maximum_covered = max(
            candidate.covered_source_bases for candidate in candidate_profiles
        )
        maximum_ids = tuple(
            candidate.candidate_id
            for candidate in candidate_profiles
            if candidate.covered_source_bases == maximum_covered
        )
        coverage_state = (
            SourceCoverageState.COMPLETE
            if maximum_covered == source_interval.length
            else SourceCoverageState.PARTIAL
        )
        union_covered = _union_covered_source_bases(candidates)

    headline = _headline(
        projection_count,
        candidate_profiles,
        maximum_covered=maximum_covered,
        union_covered=union_covered,
        source_bases=source_interval.length,
    )
    return ResultProfile(
        source_interval=source_interval,
        input_validity=InputValidityState.NOT_ASSESSED,
        projection_count=projection_count,
        source_coverage=coverage_state,
        orientation=_orientation_state(candidates),
        maximum_candidate_covered_source_bases=maximum_covered,
        source_bases=source_interval.length,
        maximum_coverage_candidate_ids=maximum_ids,
        union_covered_source_bases=union_covered,
        candidate_profiles=candidate_profiles,
        headline=headline,
        interpretation=_interpretation(headline),
        evidence_tier=evidence_tier,
        consumed_resource_roles=consumed_resource_roles,
    )


def _candidate_result_profile(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
    *,
    evidence_tier: EvidenceAvailabilityTier,
) -> CandidateResultProfile:
    _validate_candidate_geometry(source_interval, candidate)

    coverage_observation = _single_observation(
        candidate,
        EvidenceKind.MAPPING_COVERAGE,
        required=True,
    )
    if coverage_observation is None or not isinstance(
        coverage_observation.value, MappingCoverageSummary
    ):
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage must use "
            "MappingCoverageSummary"
        )
    coverage = coverage_observation.value
    expected_covered_bases = sum(
        segment.source_interval.length for segment in candidate.segments
    )
    if coverage.source_bases != source_interval.length:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage source_bases "
            "does not match the assessed source interval"
        )
    if coverage.covered_source_bases != expected_covered_bases:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage does not match "
            "its normalized mapping segments"
        )
    expected_uncovered = _uncovered_source_intervals(source_interval, candidate)
    if coverage.uncovered_source_intervals != expected_uncovered:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage uncovered "
            "intervals do not match its normalized mapping segments"
        )

    gap_observation = _single_observation(
        candidate,
        EvidenceKind.CHAIN_GAPS,
        required=True,
    )
    if gap_observation is None or not isinstance(
        gap_observation.value, ChainGapSummary
    ):
        raise ValueError(
            f"candidate {candidate.candidate_id!r} chain gaps must use ChainGapSummary"
        )
    gaps = gap_observation.value
    source_gap_intervals = tuple(
        gap.source_gap_overlap
        for gap in gaps.gaps
        if gap.source_gap_overlap is not None
    )
    target_gap_intervals = tuple(
        gap.target_gap_interval
        for gap in gaps.gaps
        if gap.target_gap_interval is not None
    )

    reciprocal_observation = _single_observation(
        candidate,
        EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
        required=evidence_tier is EvidenceAvailabilityTier.COMPARATIVE,
    )
    if evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        if reciprocal_observation is not None:
            raise ValueError(
                "LIFTOVER-ONLY result profile cannot contain reciprocal-best evidence"
            )
    elif reciprocal_observation is not None:
        if not isinstance(
            reciprocal_observation.value, ReciprocalBestMembershipSummary
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} reciprocal-best evidence must "
                "use ReciprocalBestMembershipSummary"
            )
        reciprocal = reciprocal_observation.value
        if reciprocal.candidate_source_bases != expected_covered_bases:
            raise ValueError(
                f"candidate {candidate.candidate_id!r} reciprocal-best denominator "
                "does not match its normalized mapping segments"
            )
        if any(
            not _source_interval_is_covered_by_candidate(interval, candidate)
            for interval in reciprocal.covered_source_intervals
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} reciprocal-best covered "
                "intervals must lie within normalized mapping segments"
            )

    if coverage.status is MappingCoverageStatus.FULL:
        coverage_state = SourceCoverageState.COMPLETE
    elif coverage.status is MappingCoverageStatus.PARTIAL:
        coverage_state = SourceCoverageState.PARTIAL
    else:
        assert_never(coverage.status)
    canonical_segments = _canonical_mapping_segments(candidate)
    return CandidateResultProfile(
        candidate_id=candidate.candidate_id,
        coverage_state=coverage_state,
        covered_source_bases=coverage.covered_source_bases,
        source_bases=coverage.source_bases,
        uncovered_source_intervals=coverage.uncovered_source_intervals,
        exact_mapped_segment_count=len(candidate.segments),
        geometric_segment_count=len(canonical_segments),
        fragmented=len(canonical_segments) > 1,
        target_discontinuous=bool(target_gap_intervals),
        target_bounding_span=candidate.target_interval,
        source_gap_intervals=source_gap_intervals,
        target_gap_intervals=target_gap_intervals,
        largest_uncovered_source_span_bases=max(
            (interval.length for interval in coverage.uncovered_source_intervals),
            default=0,
        ),
        largest_source_gap_bases=max(
            (interval.length for interval in source_gap_intervals),
            default=0,
        ),
        largest_target_gap_bases=max(
            (interval.length for interval in target_gap_intervals),
            default=0,
        ),
        orientation=candidate.orientation,
    )


def _projection_count(candidate_count: int) -> ProjectionCountState:
    if candidate_count == 0:
        return ProjectionCountState.NONE
    if candidate_count == 1:
        return ProjectionCountState.ONE
    return ProjectionCountState.MULTIPLE


def _orientation_state(
    candidates: tuple[NormalizedCandidate, ...],
) -> OrientationState:
    orientations = {candidate.orientation for candidate in candidates}
    if not orientations:
        return OrientationState.NONE
    if orientations == {MappingOrientation.SAME}:
        return OrientationState.SAME
    if orientations == {MappingOrientation.REVERSE}:
        return OrientationState.REVERSE
    return OrientationState.MIXED


def _headline(
    projection_count: ProjectionCountState,
    candidate_profiles: tuple[CandidateResultProfile, ...],
    *,
    maximum_covered: int,
    union_covered: int,
    source_bases: int,
) -> FactualHeadline:
    if projection_count is ProjectionCountState.NONE:
        return FactualHeadline.NO_CHAIN_PROJECTION

    if projection_count is ProjectionCountState.ONE:
        candidate = candidate_profiles[0]
        if candidate.coverage_state is SourceCoverageState.PARTIAL:
            if candidate.fragmented:
                return FactualHeadline.PARTIAL_AND_FRAGMENTED_PROJECTION
            return FactualHeadline.PARTIAL_SOURCE_COVERAGE
        if candidate.target_discontinuous:
            return FactualHeadline.COMPLETE_BUT_DISCONTINUOUS_PROJECTION
        return FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION

    if maximum_covered < source_bases and union_covered > maximum_covered:
        return FactualHeadline.SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS
    return FactualHeadline.MULTIPLE_CHAIN_PROJECTIONS


def _interpretation(headline: FactualHeadline) -> str:
    if headline is FactualHeadline.NO_CHAIN_PROJECTION:
        return (
            "No consumed chain produced a projection for this source interval; "
            "this result does not establish why."
        )
    if headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION:
        return (
            "One chain projects every requested source base; this describes coordinate "
            "geometry, not biological identity or correctness."
        )
    if headline is FactualHeadline.PARTIAL_SOURCE_COVERAGE:
        return (
            "One chain projects only part of the requested source interval; uncovered "
            "source bases remain explicit."
        )
    if headline is FactualHeadline.PARTIAL_AND_FRAGMENTED_PROJECTION:
        return (
            "One chain projects only part of the source interval and the mapped "
            "portion is split across multiple exact segments."
        )
    if headline is FactualHeadline.COMPLETE_BUT_DISCONTINUOUS_PROJECTION:
        return (
            "Every requested source base projects, but target adjacency is not "
            "preserved across the mapped segments."
        )
    if headline is FactualHeadline.MULTIPLE_CHAIN_PROJECTIONS:
        return (
            "More than one chain projection exists; candidate encounter order is not a "
            "scientific rank and this result does not choose a biological locus."
        )
    if headline is FactualHeadline.SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS:
        return (
            "Different chain projections cover different portions of the source interval; "
            "they are not equivalent to multiple complete alternative mappings."
        )
    assert_never(headline)


def _single_observation(
    candidate: NormalizedCandidate,
    kind: EvidenceKind,
    *,
    required: bool,
) -> EvidenceObservation | None:
    matches = tuple(
        observation for observation in candidate.evidence if observation.kind is kind
    )
    if len(matches) > 1:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} has duplicate {kind.value} evidence"
        )
    if not matches:
        if required:
            raise ValueError(
                f"candidate {candidate.candidate_id!r} is missing {kind.value} evidence"
            )
        return None
    return matches[0]


def _validate_candidate_ids(candidates: tuple[NormalizedCandidate, ...]) -> None:
    candidate_ids = [candidate.candidate_id for candidate in candidates]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError("candidate IDs must be unique")


def _validate_distinct_candidate_geometries(
    candidates: tuple[NormalizedCandidate, ...],
) -> None:
    seen: dict[tuple[object, ...], str] = {}
    for candidate in candidates:
        key = _canonical_mapping_geometry(candidate)
        existing = seen.get(key)
        if existing is not None:
            raise ValueError(
                "distinct candidate IDs must not describe identical normalized "
                f"mapping geometry: {existing!r} and {candidate.candidate_id!r}"
            )
        seen[key] = candidate.candidate_id


def _canonical_mapping_segments(
    candidate: NormalizedCandidate,
) -> tuple[tuple[int, int, int, int], ...]:
    canonical_segments: list[list[int]] = []
    for segment in candidate.segments:
        current = [
            segment.source_interval.start,
            segment.source_interval.end,
            segment.target_interval.start,
            segment.target_interval.end,
        ]
        if canonical_segments and _segments_are_collinear_adjacent(
            canonical_segments[-1],
            current,
            orientation=candidate.orientation,
        ):
            previous = canonical_segments[-1]
            previous[1] = current[1]
            if candidate.orientation is MappingOrientation.SAME:
                previous[3] = current[3]
            else:
                previous[2] = current[2]
        else:
            canonical_segments.append(current)

    return tuple(
        (segment[0], segment[1], segment[2], segment[3])
        for segment in canonical_segments
    )


def _canonical_mapping_geometry(candidate: NormalizedCandidate) -> tuple[object, ...]:
    canonical_segments = _canonical_mapping_segments(candidate)

    first_source = candidate.segments[0].source_interval
    target = candidate.target_interval
    return (
        first_source.assembly,
        first_source.sequence_name,
        target.assembly,
        target.sequence_name,
        candidate.orientation,
        canonical_segments,
    )


def _segments_are_collinear_adjacent(
    previous: list[int],
    current: list[int],
    *,
    orientation: MappingOrientation,
) -> bool:
    if previous[1] != current[0]:
        return False
    if orientation is MappingOrientation.SAME:
        return previous[3] == current[2]
    return previous[2] == current[3]


def _validate_candidate_geometry(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
) -> None:
    for segment in candidate.segments:
        if (
            segment.source_interval.assembly != source_interval.assembly
            or segment.source_interval.sequence_name != source_interval.sequence_name
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} source geometry does not match "
                "the assessed source locus"
            )
        if (
            segment.source_interval.start < source_interval.start
            or segment.source_interval.end > source_interval.end
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} source geometry lies outside "
                "the assessed source locus"
            )


def _uncovered_source_intervals(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
) -> tuple[GenomicInterval, ...]:
    uncovered: list[GenomicInterval] = []
    cursor = source_interval.start
    for segment in candidate.segments:
        if cursor < segment.source_interval.start:
            uncovered.append(
                GenomicInterval(
                    assembly=source_interval.assembly,
                    sequence_name=source_interval.sequence_name,
                    start=cursor,
                    end=segment.source_interval.start,
                )
            )
        cursor = segment.source_interval.end
    if cursor < source_interval.end:
        uncovered.append(
            GenomicInterval(
                assembly=source_interval.assembly,
                sequence_name=source_interval.sequence_name,
                start=cursor,
                end=source_interval.end,
            )
        )
    return tuple(uncovered)


def _source_interval_is_covered_by_candidate(
    interval: GenomicInterval,
    candidate: NormalizedCandidate,
) -> bool:
    if interval.length <= 0:
        return False
    cursor = interval.start
    for segment in candidate.segments:
        source = segment.source_interval
        if (
            source.assembly != interval.assembly
            or source.sequence_name != interval.sequence_name
        ):
            continue
        if source.end <= cursor:
            continue
        if source.start > cursor:
            return False
        cursor = min(interval.end, source.end)
        if cursor == interval.end:
            return True
    return False


def _union_covered_source_bases(candidates: tuple[NormalizedCandidate, ...]) -> int:
    spans = sorted(
        (
            segment.source_interval.start,
            segment.source_interval.end,
        )
        for candidate in candidates
        for segment in candidate.segments
    )
    if not spans:
        return 0
    total = 0
    current_start, current_end = spans[0]
    for start, end in spans[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
            continue
        total += current_end - current_start
        current_start, current_end = start, end
    return total + current_end - current_start
