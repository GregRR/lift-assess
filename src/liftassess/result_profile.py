"""Derived factual result profile for liftAssess.

This module is the boundary between scientific candidate/evidence data and result
rendering.  It derives deterministic factual states from already-computed mapping
geometry and evidence without assigning an aggregate verdict, confidence score, or
biological truth claim.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import assert_never

from ._candidate_geometry import (
    canonical_mapping_segments,
    validate_distinct_candidate_geometries,
)
from .comparative_inventory import (
    FilteredAllChainComparisonResult,
    FilteredAllChainInventoryState,
)
from .comparative_relationship import (
    ComparativeEvidenceRelationshipResult,
    build_comparative_evidence_relationship,
)
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
from .query_context import (
    PointQueryContextResult,
    QueryContextNotRunReason,
    QueryContextState,
)
from .reverse_mapping import (
    CandidateReverseMappingResult,
    ReverseCheckState,
    ReverseOriginalSourceCoverageState,
    ReverseRelationshipState,
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


class ComparativeRelationshipState(str, Enum):
    """Cross-resource comparative synthesis state."""

    NOT_ASSESSED = "NOT_ASSESSED"
    NO_COMPETING_FULL_PLACEMENTS = "NO_COMPETING_FULL_PLACEMENTS"
    FAVORS_ONE_PLACEMENT = "FAVORS_ONE_PLACEMENT"
    DOES_NOT_SEPARATE_PLACEMENTS = "DOES_NOT_SEPARATE_PLACEMENTS"
    MIXED_CONFLICTING = "MIXED_CONFLICTING"


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


class QueryContextFinding(str, Enum):
    """Factual relationships between a point and its tested local context."""

    AGREES_WITH_POINT = "AGREES_WITH_POINT"
    NO_PROJECTION_AT_EITHER_SCALE = "NO_PROJECTION_AT_EITHER_SCALE"
    REVEALS_PARTIAL_COVERAGE = "REVEALS_PARTIAL_COVERAGE"
    REVEALS_FRAGMENTATION = "REVEALS_FRAGMENTATION"
    REVEALS_TARGET_DISCONTINUITY = "REVEALS_TARGET_DISCONTINUITY"
    CHANGES_WITH_QUERY_SCALE = "CHANGES_WITH_QUERY_SCALE"


_QUERY_CONTEXT_CHAIN_EVIDENCE_KINDS = frozenset(
    {
        EvidenceKind.CHAIN_GAPS,
        EvidenceKind.CHAIN_SCORE,
        EvidenceKind.MAPPING_COVERAGE,
    }
)


@dataclass(frozen=True)
class CandidateReverseMappingProfile:
    """Derived candidate-level facts from an actual reverse-mapping check."""

    check_state: ReverseCheckState
    relationship: ReverseRelationshipState | None
    original_source_bases: int
    original_source_covered_bases: int | None
    original_source_coverage: ReverseOriginalSourceCoverageState | None
    exact_original_geometry_return: bool | None
    reverse_projection_count: int | None
    segments_with_reverse_projection: int | None
    queried_target_segments: tuple[GenomicInterval, ...]


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
    reverse_mapping: CandidateReverseMappingProfile


@dataclass(frozen=True)
class ComparativePlacementProfile:
    """Renderer-facing categorical support facts for one all-chain placement."""

    candidate_id: str
    complete_source_coverage: bool
    retained_by_filtered_chain: bool
    depth1_top_net: bool
    full_reciprocal_best: bool


@dataclass(frozen=True)
class ComparativeRelationshipProfile:
    """Derived paired-inventory and categorical comparative facts."""

    state: ComparativeRelationshipState
    inventory_state: FilteredAllChainInventoryState | None
    favored_candidate_id: str | None
    additional_all_chain_candidate_ids: tuple[str, ...]
    placement_support: tuple[ComparativePlacementProfile, ...]


@dataclass(frozen=True)
class QueryContextProfile:
    """Derived facts for one automatic or explicitly sized point neighborhood."""

    check_state: QueryContextState
    findings: tuple[QueryContextFinding, ...]
    requested_window_bases: int | None
    tested_source_interval: GenomicInterval | None
    actual_window_bases: int | None
    not_run_reason: QueryContextNotRunReason | None
    projection_count: ProjectionCountState | None
    source_coverage: SourceCoverageState | None
    maximum_candidate_covered_source_bases: int | None
    union_covered_source_bases: int | None
    candidate_profiles: tuple[CandidateResultProfile, ...]
    headline: FactualHeadline | None
    point_and_local_context_map_together: bool


@dataclass(frozen=True)
class ResultScopeProfile:
    """Explicit scope states for orthogonal result dimensions."""

    target_role: TargetRoleState = TargetRoleState.NOT_ASSESSED
    reverse_result: ReverseCheckState = ReverseCheckState.NOT_RUN
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
    query_context: QueryContextProfile
    comparative_relationship: ComparativeRelationshipProfile
    scope: ResultScopeProfile = field(default_factory=ResultScopeProfile)


def build_result_profile(
    source_interval: GenomicInterval,
    candidates: tuple[NormalizedCandidate, ...],
    *,
    evidence_tier: EvidenceAvailabilityTier,
    consumed_resource_roles: tuple[str, ...] = (),
    reverse_mapping_results: tuple[CandidateReverseMappingResult, ...] | None = None,
    query_context_result: PointQueryContextResult | None = None,
    filtered_all_chain_comparison: FilteredAllChainComparisonResult | None = None,
    comparative_evidence_relationship: (
        ComparativeEvidenceRelationshipResult | None
    ) = None,
) -> ResultProfile:
    """Derive one deterministic factual profile from normalized scientific data."""

    if source_interval.length <= 0:
        raise ValueError("result profile requires a non-empty source interval")

    _validate_candidate_ids(candidates)
    validate_distinct_candidate_geometries(candidates)
    reverse_profiles = _reverse_mapping_profiles(candidates, reverse_mapping_results)
    candidate_profiles = tuple(
        _candidate_result_profile(
            source_interval,
            candidate,
            evidence_tier=evidence_tier,
            reverse_mapping=reverse_profile,
        )
        for candidate, reverse_profile in zip(candidates, reverse_profiles, strict=True)
    )

    projection_count = _projection_count(len(candidates))
    maximum_covered, maximum_ids, coverage_state, union_covered = (
        _aggregate_candidate_coverage(
            source_interval,
            candidates,
            candidate_profiles,
        )
    )

    headline = _headline(
        projection_count,
        candidate_profiles,
        maximum_covered=maximum_covered,
        union_covered=union_covered,
        source_bases=source_interval.length,
    )
    query_context = _query_context_profile(
        source_interval,
        candidates,
        candidate_profiles,
        evidence_tier=evidence_tier,
        query_context_result=query_context_result,
    )
    comparative_relationship = _comparative_relationship_profile(
        candidates,
        evidence_tier=evidence_tier,
        comparison=filtered_all_chain_comparison,
        relationship=comparative_evidence_relationship,
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
        query_context=query_context,
        comparative_relationship=comparative_relationship,
        scope=ResultScopeProfile(
            reverse_result=_reverse_scope_state(reverse_profiles),
            query_context=query_context.check_state,
            comparative_relationship=comparative_relationship.state,
        ),
    )


def _comparative_relationship_profile(
    candidates: tuple[NormalizedCandidate, ...],
    *,
    evidence_tier: EvidenceAvailabilityTier,
    comparison: FilteredAllChainComparisonResult | None,
    relationship: ComparativeEvidenceRelationshipResult | None,
) -> ComparativeRelationshipProfile:
    if comparison is None and relationship is None:
        return ComparativeRelationshipProfile(
            state=ComparativeRelationshipState.NOT_ASSESSED,
            inventory_state=None,
            favored_candidate_id=None,
            additional_all_chain_candidate_ids=(),
            placement_support=(),
        )
    if comparison is None or relationship is None:
        raise ValueError(
            "comparative result profile requires both paired inventory and "
            "categorical relationship"
        )
    if evidence_tier is not EvidenceAvailabilityTier.COMPARATIVE:
        raise ValueError(
            "comparative result profile requires COMPARATIVE evidence tier"
        )
    if comparison.all_chain_candidates != candidates:
        raise ValueError(
            "comparative result profile all-chain inventory must match candidates"
        )
    expected_relationship = build_comparative_evidence_relationship(comparison)
    if relationship != expected_relationship:
        raise ValueError(
            "comparative result profile relationship must match paired inventory "
            "and candidate evidence"
        )

    return ComparativeRelationshipProfile(
        state=ComparativeRelationshipState(relationship.relationship.value),
        inventory_state=comparison.relationship,
        favored_candidate_id=relationship.favored_candidate_id,
        additional_all_chain_candidate_ids=(
            comparison.additional_all_chain_candidate_ids
        ),
        placement_support=tuple(
            ComparativePlacementProfile(
                candidate_id=item.candidate_id,
                complete_source_coverage=item.complete_source_coverage,
                retained_by_filtered_chain=item.retained_by_filtered_chain,
                depth1_top_net=item.depth1_top_net,
                full_reciprocal_best=item.full_reciprocal_best,
            )
            for item in relationship.placement_support
        ),
    )


def _query_context_profile(
    source_interval: GenomicInterval,
    point_candidates: tuple[NormalizedCandidate, ...],
    point_candidate_profiles: tuple[CandidateResultProfile, ...],
    *,
    evidence_tier: EvidenceAvailabilityTier,
    query_context_result: PointQueryContextResult | None,
) -> QueryContextProfile:
    if query_context_result is None:
        return QueryContextProfile(
            check_state=QueryContextState.NOT_RUN,
            findings=(),
            requested_window_bases=None,
            tested_source_interval=None,
            actual_window_bases=None,
            not_run_reason=None,
            projection_count=None,
            source_coverage=None,
            maximum_candidate_covered_source_bases=None,
            union_covered_source_bases=None,
            candidate_profiles=(),
            headline=None,
            point_and_local_context_map_together=False,
        )

    if source_interval.length != 1:
        raise ValueError("point query context can only be attached to a one-base query")
    if query_context_result.check_state is QueryContextState.NOT_RUN:
        return QueryContextProfile(
            check_state=QueryContextState.NOT_RUN,
            findings=(),
            requested_window_bases=query_context_result.requested_window_bases,
            tested_source_interval=None,
            actual_window_bases=None,
            not_run_reason=query_context_result.not_run_reason,
            projection_count=None,
            source_coverage=None,
            maximum_candidate_covered_source_bases=None,
            union_covered_source_bases=None,
            candidate_profiles=(),
            headline=None,
            point_and_local_context_map_together=False,
        )

    tested_interval = query_context_result.tested_source_interval
    if tested_interval is None:
        raise ValueError("completed point context requires its tested source interval")
    if (
        tested_interval.assembly != source_interval.assembly
        or tested_interval.sequence_name != source_interval.sequence_name
        or tested_interval.start > source_interval.start
        or tested_interval.end < source_interval.end
    ):
        raise ValueError(
            "point context interval must contain the assessed source point"
        )

    context_candidates = query_context_result.candidates
    _validate_query_context_evidence_scope(context_candidates)
    _validate_candidate_ids(context_candidates)
    validate_distinct_candidate_geometries(context_candidates)
    context_reverse_profiles = _reverse_mapping_profiles(context_candidates, None)
    context_candidate_profiles = tuple(
        _candidate_result_profile(
            tested_interval,
            candidate,
            evidence_tier=evidence_tier,
            reverse_mapping=reverse_profile,
            require_comparative_evidence=False,
        )
        for candidate, reverse_profile in zip(
            context_candidates,
            context_reverse_profiles,
            strict=True,
        )
    )
    projection_count = _projection_count(len(context_candidates))
    maximum_covered, _, coverage_state, union_covered = _aggregate_candidate_coverage(
        tested_interval,
        context_candidates,
        context_candidate_profiles,
    )
    context_headline = _headline(
        projection_count,
        context_candidate_profiles,
        maximum_covered=maximum_covered,
        union_covered=union_covered,
        source_bases=tested_interval.length,
    )

    point_candidates_by_id = {
        candidate.candidate_id: candidate for candidate in point_candidates
    }
    context_candidates_by_id = {
        candidate.candidate_id: candidate for candidate in context_candidates
    }
    point_candidate_ids = frozenset(point_candidates_by_id)
    context_candidate_ids = frozenset(context_candidates_by_id)
    for candidate_id in point_candidate_ids & context_candidate_ids:
        _validate_point_context_candidate_geometry(
            point_candidates_by_id[candidate_id],
            context_candidates_by_id[candidate_id],
        )
    reveals_partial_coverage = any(
        candidate.coverage_state is SourceCoverageState.PARTIAL
        for candidate in context_candidate_profiles
    )
    reveals_fragmentation = any(
        candidate.fragmented for candidate in context_candidate_profiles
    )
    reveals_target_discontinuity = any(
        candidate.target_discontinuous for candidate in context_candidate_profiles
    )
    changes_with_scale = (
        point_candidate_ids != context_candidate_ids
        or reveals_partial_coverage
        or reveals_fragmentation
        or reveals_target_discontinuity
    )
    findings: list[QueryContextFinding] = []
    if not point_candidate_ids and not context_candidate_ids:
        findings.append(QueryContextFinding.NO_PROJECTION_AT_EITHER_SCALE)
    elif not changes_with_scale:
        findings.append(QueryContextFinding.AGREES_WITH_POINT)
    if reveals_partial_coverage:
        findings.append(QueryContextFinding.REVEALS_PARTIAL_COVERAGE)
    if reveals_fragmentation:
        findings.append(QueryContextFinding.REVEALS_FRAGMENTATION)
    if reveals_target_discontinuity:
        findings.append(QueryContextFinding.REVEALS_TARGET_DISCONTINUITY)
    if changes_with_scale:
        findings.append(QueryContextFinding.CHANGES_WITH_QUERY_SCALE)

    point_and_local_context_map_together = (
        len(point_candidate_profiles) == 1
        and len(context_candidate_profiles) == 1
        and point_candidate_profiles[0].candidate_id
        == context_candidate_profiles[0].candidate_id
        and context_candidate_profiles[0].coverage_state is SourceCoverageState.COMPLETE
        and not context_candidate_profiles[0].fragmented
        and not context_candidate_profiles[0].target_discontinuous
    )
    return QueryContextProfile(
        check_state=QueryContextState.RUN,
        findings=tuple(findings),
        requested_window_bases=query_context_result.requested_window_bases,
        tested_source_interval=tested_interval,
        actual_window_bases=tested_interval.length,
        not_run_reason=None,
        projection_count=projection_count,
        source_coverage=coverage_state,
        maximum_candidate_covered_source_bases=maximum_covered,
        union_covered_source_bases=union_covered,
        candidate_profiles=context_candidate_profiles,
        headline=context_headline,
        point_and_local_context_map_together=point_and_local_context_map_together,
    )


def _validate_query_context_evidence_scope(
    candidates: tuple[NormalizedCandidate, ...],
) -> None:
    invalid_kinds = sorted(
        {
            observation.kind.value
            for candidate in candidates
            for observation in candidate.evidence
            if observation.kind not in _QUERY_CONTEXT_CHAIN_EVIDENCE_KINDS
        }
    )
    if invalid_kinds:
        raise ValueError(
            "point query context can contain only forward-chain evidence; found "
            + ", ".join(invalid_kinds)
        )


def _validate_point_context_candidate_geometry(
    point_candidate: NormalizedCandidate,
    context_candidate: NormalizedCandidate,
) -> None:
    """Require a shared chain candidate to reproduce the point projection exactly."""

    mismatch = (
        "query-context candidate geometry must reproduce the point mapping "
        "for shared candidate IDs"
    )
    if point_candidate.orientation is not context_candidate.orientation:
        raise ValueError(mismatch)
    if point_candidate.mapping_provenance != context_candidate.mapping_provenance:
        raise ValueError(mismatch)
    if len(point_candidate.segments) != 1:
        raise ValueError(mismatch)

    point_segment = point_candidate.segments[0]
    point_source = point_segment.source_interval
    matching_context_segments = tuple(
        segment
        for segment in context_candidate.segments
        if segment.source_interval.assembly == point_source.assembly
        and segment.source_interval.sequence_name == point_source.sequence_name
        and segment.source_interval.start <= point_source.start
        and segment.source_interval.end >= point_source.end
    )
    if len(matching_context_segments) != 1:
        raise ValueError(mismatch)

    context_segment = matching_context_segments[0]
    context_source = context_segment.source_interval
    context_target = context_segment.target_interval
    offset_start = point_source.start - context_source.start
    offset_end = point_source.end - context_source.start
    if context_candidate.orientation is MappingOrientation.SAME:
        expected_target_start = context_target.start + offset_start
        expected_target_end = context_target.start + offset_end
    else:
        expected_target_start = context_target.end - offset_end
        expected_target_end = context_target.end - offset_start

    point_target = point_segment.target_interval
    if (
        point_target.assembly != context_target.assembly
        or point_target.sequence_name != context_target.sequence_name
        or point_target.start != expected_target_start
        or point_target.end != expected_target_end
    ):
        raise ValueError(mismatch)


def _reverse_mapping_profiles(
    candidates: tuple[NormalizedCandidate, ...],
    reverse_mapping_results: tuple[CandidateReverseMappingResult, ...] | None,
) -> tuple[CandidateReverseMappingProfile, ...]:
    if reverse_mapping_results is None:
        return tuple(
            CandidateReverseMappingProfile(
                check_state=ReverseCheckState.NOT_RUN,
                relationship=None,
                original_source_bases=sum(
                    segment.source_interval.length for segment in candidate.segments
                ),
                original_source_covered_bases=None,
                original_source_coverage=None,
                exact_original_geometry_return=None,
                reverse_projection_count=None,
                segments_with_reverse_projection=None,
                queried_target_segments=tuple(
                    segment.target_interval for segment in candidate.segments
                ),
            )
            for candidate in candidates
        )

    if len(reverse_mapping_results) != len(candidates):
        raise ValueError("reverse mapping results must match the candidate count")

    profiles: list[CandidateReverseMappingProfile] = []
    for candidate, result in zip(candidates, reverse_mapping_results, strict=True):
        if result.forward_candidate_id != candidate.candidate_id:
            raise ValueError(
                "reverse mapping results must preserve forward candidate order"
            )
        expected_original_source_segments = tuple(
            segment.source_interval for segment in candidate.segments
        )
        expected_queried_target_segments = tuple(
            segment.target_interval for segment in candidate.segments
        )
        if result.original_source_segments != expected_original_source_segments:
            raise ValueError(
                "reverse mapping original-source geometry must match forward candidate"
            )
        if result.queried_target_segments != expected_queried_target_segments:
            raise ValueError(
                "reverse mapping query geometry must match forward candidate"
            )
        if result.check_state is ReverseCheckState.RUN:
            profiles.append(
                CandidateReverseMappingProfile(
                    check_state=result.check_state,
                    relationship=result.relationship,
                    original_source_bases=result.original_source_bases,
                    original_source_covered_bases=result.original_source_covered_bases,
                    original_source_coverage=result.original_source_coverage,
                    exact_original_geometry_return=(
                        result.exact_original_geometry_return
                    ),
                    reverse_projection_count=result.reverse_projection_count,
                    segments_with_reverse_projection=(
                        result.segments_with_reverse_projection
                    ),
                    queried_target_segments=result.queried_target_segments,
                )
            )
        else:
            profiles.append(
                CandidateReverseMappingProfile(
                    check_state=result.check_state,
                    relationship=None,
                    original_source_bases=result.original_source_bases,
                    original_source_covered_bases=None,
                    original_source_coverage=None,
                    exact_original_geometry_return=None,
                    reverse_projection_count=None,
                    segments_with_reverse_projection=None,
                    queried_target_segments=result.queried_target_segments,
                )
            )
    return tuple(profiles)


def _reverse_scope_state(
    reverse_profiles: tuple[CandidateReverseMappingProfile, ...],
) -> ReverseCheckState:
    if not reverse_profiles:
        return ReverseCheckState.NOT_RUN
    states = {profile.check_state for profile in reverse_profiles}
    if len(states) != 1:
        raise ValueError(
            "reverse mapping check state must be uniform across candidates"
        )
    return next(iter(states))


def _candidate_result_profile(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
    *,
    evidence_tier: EvidenceAvailabilityTier,
    reverse_mapping: CandidateReverseMappingProfile,
    require_comparative_evidence: bool = True,
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
        required=(
            evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
            and require_comparative_evidence
        ),
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
    canonical_segments = canonical_mapping_segments(candidate)
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
        reverse_mapping=reverse_mapping,
    )


def _aggregate_candidate_coverage(
    source_interval: GenomicInterval,
    candidates: tuple[NormalizedCandidate, ...],
    candidate_profiles: tuple[CandidateResultProfile, ...],
) -> tuple[int, tuple[str, ...], SourceCoverageState, int]:
    if not candidate_profiles:
        return 0, (), SourceCoverageState.NONE, 0

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
    return (
        maximum_covered,
        maximum_ids,
        coverage_state,
        _union_covered_source_bases(candidates),
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
