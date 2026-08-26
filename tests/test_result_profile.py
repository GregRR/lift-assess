from __future__ import annotations

from dataclasses import replace

import pytest

from liftassess import (
    AssemblyIdentifier,
    ChainGap,
    ChainGapSummary,
    ComparativeEvidenceRelationshipResult,
    ComparativeRelationshipState,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    FactualHeadline,
    FilteredAllChainComparisonResult,
    FilteredAllChainInventoryState,
    GenomicInterval,
    InputValidityState,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    MappingSegment,
    NetHierarchySummary,
    NormalizedCandidate,
    OrientationState,
    PointQueryContextResult,
    ProjectionCountState,
    ProvenanceSource,
    QueryContextFinding,
    QueryContextState,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
    SourceCoverageState,
    build_comparative_evidence_relationship,
    build_filtered_all_chain_comparison,
    build_result_profile,
    reverse_mapping_unavailable,
)

SOURCE_ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")
TARGET_ASSEMBLY = AssemblyIdentifier("targetAsm", "test")
SOURCE = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 200)
ALIGNMENT = ProvenanceSource("alignment", "shared test alignment")
CHAIN = ProvenanceSource("chain", "chain resource", derived_from=(ALIGNMENT,))
FILTERED_CHAIN = ProvenanceSource(
    "filtered-chain", "filtered chain resource", derived_from=(ALIGNMENT,)
)
NET = ProvenanceSource("net", "net resource", derived_from=(ALIGNMENT,))
RBEST = ProvenanceSource("rbest", "rbest resource", derived_from=(ALIGNMENT,))
POINT = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 150, 151)
POINT_CONTEXT = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 201)


def _candidate(
    candidate_id: str,
    *,
    source_interval: GenomicInterval = SOURCE,
    source_spans: tuple[tuple[int, int], ...] = ((100, 200),),
    target_spans: tuple[tuple[int, int], ...] = ((1000, 1100),),
    orientation: MappingOrientation = MappingOrientation.SAME,
    target_gaps: tuple[tuple[int, int], ...] = (),
    reciprocal_best: ReciprocalBestMembershipStatus | None = None,
    covered_source_bases_override: int | None = None,
) -> NormalizedCandidate:
    segments = tuple(
        MappingSegment(
            GenomicInterval(
                source_interval.assembly,
                source_interval.sequence_name,
                source_start,
                source_end,
            ),
            GenomicInterval(TARGET_ASSEMBLY, "chrA", target_start, target_end),
        )
        for (source_start, source_end), (target_start, target_end) in zip(
            source_spans, target_spans, strict=True
        )
    )
    covered = sum(end - start for start, end in source_spans)
    uncovered: list[GenomicInterval] = []
    cursor = source_interval.start
    for start, end in source_spans:
        if cursor < start:
            uncovered.append(
                GenomicInterval(
                    source_interval.assembly,
                    source_interval.sequence_name,
                    cursor,
                    start,
                )
            )
        cursor = end
    if cursor < source_interval.end:
        uncovered.append(
            GenomicInterval(
                source_interval.assembly,
                source_interval.sequence_name,
                cursor,
                source_interval.end,
            )
        )

    coverage = EvidenceObservation(
        f"{candidate_id}:coverage",
        EvidenceKind.MAPPING_COVERAGE,
        MappingCoverageSummary(
            status=(
                MappingCoverageStatus.FULL
                if covered == source_interval.length
                else MappingCoverageStatus.PARTIAL
            ),
            covered_source_bases=(
                covered
                if covered_source_bases_override is None
                else covered_source_bases_override
            ),
            source_bases=source_interval.length,
            uncovered_source_intervals=tuple(uncovered),
        ),
        CHAIN,
    )
    gaps = EvidenceObservation(
        f"{candidate_id}:gaps",
        EvidenceKind.CHAIN_GAPS,
        ChainGapSummary(
            tuple(
                ChainGap(
                    source_boundary=source_spans[min(index + 1, len(source_spans) - 1)][
                        0
                    ],
                    target_gap_interval=GenomicInterval(
                        TARGET_ASSEMBLY, "chrA", target_start, target_end
                    ),
                )
                for index, (target_start, target_end) in enumerate(target_gaps)
            )
        ),
        CHAIN,
    )
    evidence: list[EvidenceObservation] = [coverage, gaps]
    if reciprocal_best is not None:
        if reciprocal_best is ReciprocalBestMembershipStatus.FULL:
            rbest_covered = covered
            rbest_intervals = tuple(segment.source_interval for segment in segments)
        elif reciprocal_best is ReciprocalBestMembershipStatus.PARTIAL:
            rbest_covered = covered - 1
            rbest_intervals = (
                GenomicInterval(
                    SOURCE_ASSEMBLY,
                    "chr1",
                    source_spans[0][0],
                    source_spans[0][1] - 1,
                ),
            )
        else:
            rbest_covered = 0
            rbest_intervals = ()
        evidence.append(
            EvidenceObservation(
                f"{candidate_id}:rbest",
                EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
                ReciprocalBestMembershipSummary(
                    status=reciprocal_best,
                    resource_completeness=(
                        ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
                    ),
                    chains_examined=1,
                    covered_source_bases=rbest_covered,
                    candidate_source_bases=covered,
                    covered_source_intervals=rbest_intervals,
                ),
                RBEST,
            )
        )

    target_start = min(start for start, _ in target_spans)
    target_end = max(end for _, end in target_spans)
    return NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=GenomicInterval(
            TARGET_ASSEMBLY, "chrA", target_start, target_end
        ),
        orientation=orientation,
        mapping_provenance=CHAIN,
        segments=segments,
        evidence=tuple(evidence),
    )


def _filtered_candidate(candidate: NormalizedCandidate) -> NormalizedCandidate:
    chain_evidence = tuple(
        replace(observation, provenance=FILTERED_CHAIN)
        for observation in candidate.evidence
        if observation.kind in {EvidenceKind.MAPPING_COVERAGE, EvidenceKind.CHAIN_GAPS}
    )
    return replace(
        candidate,
        candidate_id=f"filtered:{candidate.candidate_id}",
        mapping_provenance=FILTERED_CHAIN,
        evidence=chain_evidence,
    )


def _with_depth1_top_net(candidate: NormalizedCandidate) -> NormalizedCandidate:
    fill = ProvenanceSource(
        f"{candidate.candidate_id}:fill",
        "top net fill",
        derived_from=(NET,),
    )
    return replace(
        candidate,
        evidence=candidate.evidence
        + (
            EvidenceObservation(
                f"{candidate.candidate_id}:net:classification",
                EvidenceKind.NET_CLASSIFICATION,
                "top",
                fill,
            ),
            EvidenceObservation(
                f"{candidate.candidate_id}:net:hierarchy",
                EvidenceKind.NET_HIERARCHY,
                NetHierarchySummary(depth=1, source_fill_interval=SOURCE),
                fill,
            ),
        ),
    )


def _comparative_inputs(
    candidates: tuple[NormalizedCandidate, ...],
    *,
    filtered_candidate: NormalizedCandidate,
) -> tuple[
    FilteredAllChainComparisonResult,
    ComparativeEvidenceRelationshipResult,
]:
    comparison = build_filtered_all_chain_comparison(
        SOURCE,
        candidates,
        (filtered_candidate,),
        all_chain_provenance=CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )
    return comparison, build_comparative_evidence_relationship(comparison)


def test_no_projection_has_factual_no_projection_profile() -> None:
    profile = build_result_profile(
        SOURCE,
        (),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        consumed_resource_roles=("CHAIN",),
    )

    assert profile.input_validity is InputValidityState.NOT_ASSESSED
    assert profile.projection_count is ProjectionCountState.NONE
    assert profile.source_coverage is SourceCoverageState.NONE
    assert profile.orientation is OrientationState.NONE
    assert profile.maximum_candidate_covered_source_bases == 0
    assert profile.headline is FactualHeadline.NO_CHAIN_PROJECTION
    assert "does not establish why" in profile.interpretation


def test_one_complete_contiguous_projection_has_complete_profile() -> None:
    profile = build_result_profile(
        SOURCE,
        (_candidate("one"),),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        consumed_resource_roles=("CHAIN",),
    )

    candidate = profile.candidate_profiles[0]
    assert profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION
    assert profile.source_coverage is SourceCoverageState.COMPLETE
    assert profile.maximum_candidate_covered_source_bases == 100
    assert candidate.coverage_state is SourceCoverageState.COMPLETE
    assert candidate.exact_mapped_segment_count == 1
    assert candidate.geometric_segment_count == 1
    assert not candidate.fragmented
    assert not candidate.target_discontinuous


def test_partial_projection_reports_exact_uncovered_source_geometry() -> None:
    profile = build_result_profile(
        SOURCE,
        (
            _candidate(
                "partial",
                source_spans=((100, 190),),
                target_spans=((1000, 1090),),
            ),
        ),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    candidate = profile.candidate_profiles[0]
    assert profile.headline is FactualHeadline.PARTIAL_SOURCE_COVERAGE
    assert candidate.covered_source_bases == 90
    assert candidate.uncovered_source_intervals == (
        GenomicInterval(SOURCE_ASSEMBLY, "chr1", 190, 200),
    )
    assert candidate.largest_uncovered_source_span_bases == 10


def test_partial_fragmented_projection_reports_segments_and_target_gap() -> None:
    profile = build_result_profile(
        SOURCE,
        (
            _candidate(
                "split",
                source_spans=((100, 150), (160, 190)),
                target_spans=((1000, 1050), (1060, 1090)),
                target_gaps=((1050, 1060),),
            ),
        ),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    candidate = profile.candidate_profiles[0]
    assert profile.headline is FactualHeadline.PARTIAL_AND_FRAGMENTED_PROJECTION
    assert candidate.exact_mapped_segment_count == 2
    assert candidate.geometric_segment_count == 2
    assert candidate.fragmented
    assert candidate.target_discontinuous
    assert candidate.largest_target_gap_bases == 10


def test_complete_discontinuous_projection_keeps_bounding_span_distinct() -> None:
    profile = build_result_profile(
        SOURCE,
        (
            _candidate(
                "discontinuous",
                source_spans=((100, 150), (150, 200)),
                target_spans=((1000, 1050), (1100, 1150)),
                target_gaps=((1050, 1100),),
            ),
        ),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    candidate = profile.candidate_profiles[0]
    assert profile.headline is FactualHeadline.COMPLETE_BUT_DISCONTINUOUS_PROJECTION
    assert candidate.target_bounding_span == GenomicInterval(
        TARGET_ASSEMBLY, "chrA", 1000, 1150
    )
    assert candidate.largest_target_gap_bases == 50


def test_multiple_complete_projections_are_factual_multiplicity_without_rank() -> None:
    profile = build_result_profile(
        SOURCE,
        (
            _candidate("first"),
            _candidate("second", target_spans=((2000, 2100),)),
        ),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert profile.headline is FactualHeadline.MULTIPLE_CHAIN_PROJECTIONS
    assert profile.maximum_candidate_covered_source_bases == 100
    assert profile.maximum_coverage_candidate_ids == ("first", "second")


def test_source_split_across_projections_is_distinct_from_alternatives() -> None:
    profile = build_result_profile(
        SOURCE,
        (
            _candidate(
                "left",
                source_spans=((100, 150),),
                target_spans=((1000, 1050),),
            ),
            _candidate(
                "right",
                source_spans=((150, 200),),
                target_spans=((3000, 3050),),
            ),
        ),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert (
        profile.headline
        is FactualHeadline.SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS
    )
    assert profile.maximum_candidate_covered_source_bases == 50
    assert profile.union_covered_source_bases == 100


def test_orientation_state_reports_mixed_candidate_orientations() -> None:
    profile = build_result_profile(
        SOURCE,
        (
            _candidate("same"),
            _candidate(
                "reverse",
                target_spans=((2000, 2100),),
                orientation=MappingOrientation.REVERSE,
            ),
        ),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert profile.orientation is OrientationState.MIXED


def test_comparative_tier_preserves_rbest_without_aggregate_result() -> None:
    profile = build_result_profile(
        SOURCE,
        (_candidate("one", reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL),),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        consumed_resource_roles=("CHAIN", "NET", "RECIPROCAL_BEST_CHAIN"),
    )

    assert profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION
    assert profile.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
    assert profile.scope.comparative_relationship.value == "NOT_ASSESSED"


def test_comparative_relationship_profile_preserves_b14_style_support() -> None:
    favored = _with_depth1_top_net(
        _candidate(
            "favored",
            reciprocal_best=ReciprocalBestMembershipStatus.FULL,
        )
    )
    competitor = _candidate(
        "competitor",
        target_spans=((2000, 2100),),
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )
    comparison, relationship = _comparative_inputs(
        (favored, competitor),
        filtered_candidate=_filtered_candidate(favored),
    )

    profile = build_result_profile(
        SOURCE,
        (favored, competitor),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        filtered_all_chain_comparison=comparison,
        comparative_evidence_relationship=relationship,
    )

    comparative = profile.comparative_relationship
    assert comparative.state is ComparativeRelationshipState.FAVORS_ONE_PLACEMENT
    assert profile.scope.comparative_relationship is comparative.state
    assert (
        comparative.inventory_state
        is FilteredAllChainInventoryState.ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS
    )
    assert comparative.favored_candidate_id == "favored"
    assert comparative.additional_all_chain_candidate_ids == ("competitor",)
    assert tuple(item.candidate_id for item in comparative.placement_support) == (
        "favored",
        "competitor",
    )
    favored_support = comparative.placement_support[0]
    assert favored_support.complete_source_coverage
    assert favored_support.retained_by_filtered_chain
    assert favored_support.depth1_top_net
    assert favored_support.full_reciprocal_best


def test_zero_zero_comparative_inventory_is_structured_in_profile() -> None:
    comparison = build_filtered_all_chain_comparison(
        SOURCE,
        (),
        (),
        all_chain_provenance=CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )
    relationship = build_comparative_evidence_relationship(comparison)

    profile = build_result_profile(
        SOURCE,
        (),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        filtered_all_chain_comparison=comparison,
        comparative_evidence_relationship=relationship,
    )

    comparative = profile.comparative_relationship
    assert (
        comparative.state is ComparativeRelationshipState.NO_COMPETING_FULL_PLACEMENTS
    )
    assert (
        comparative.inventory_state
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    )
    assert comparative.placement_support == ()
    assert comparative.additional_all_chain_candidate_ids == ()


def test_comparative_profile_requires_inventory_and_relationship_together() -> None:
    candidate = _candidate(
        "one",
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )
    comparison, relationship = _comparative_inputs(
        (candidate,),
        filtered_candidate=_filtered_candidate(candidate),
    )

    with pytest.raises(ValueError, match="requires both paired inventory"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            filtered_all_chain_comparison=comparison,
        )
    with pytest.raises(ValueError, match="requires both paired inventory"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            comparative_evidence_relationship=relationship,
        )


def test_no_projection_comparative_profile_does_not_require_candidate_evidence() -> (
    None
):
    profile = build_result_profile(
        SOURCE,
        (),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        consumed_resource_roles=("CHAIN",),
    )

    assert profile.headline is FactualHeadline.NO_CHAIN_PROJECTION
    assert profile.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
    assert profile.consumed_resource_roles == ("CHAIN",)


def test_profile_rejects_duplicate_candidate_ids() -> None:
    first = _candidate("duplicate")
    second = _candidate(
        "duplicate",
        target_spans=((2000, 2100),),
    )

    with pytest.raises(ValueError, match="candidate IDs must be unique"):
        build_result_profile(
            SOURCE,
            (first, second),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_aggregate_profile_semantics_do_not_depend_on_candidate_order() -> None:
    first = _candidate(
        "first",
        source_spans=((100, 160),),
        target_spans=((1000, 1060),),
    )
    second = _candidate(
        "second",
        source_spans=((160, 200),),
        target_spans=((3000, 3040),),
        orientation=MappingOrientation.REVERSE,
    )

    forward = build_result_profile(
        SOURCE,
        (first, second),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )
    reversed_order = build_result_profile(
        SOURCE,
        (second, first),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert forward.headline is reversed_order.headline
    assert forward.source_coverage is reversed_order.source_coverage
    assert (
        forward.maximum_candidate_covered_source_bases
        == reversed_order.maximum_candidate_covered_source_bases
    )
    assert (
        forward.union_covered_source_bases == reversed_order.union_covered_source_bases
    )
    assert forward.orientation is reversed_order.orientation


def test_equivalent_geometry_is_rejected_after_partition_canonicalization() -> None:
    whole = _candidate("whole")
    split = _candidate(
        "split",
        source_spans=((100, 150), (150, 200)),
        target_spans=((1000, 1050), (1050, 1100)),
    )

    with pytest.raises(ValueError, match="identical normalized mapping geometry"):
        build_result_profile(
            SOURCE,
            (whole, split),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_adjacent_collinear_blocks_are_not_geometric_fragmentation() -> None:
    candidate = _candidate(
        "partitioned",
        source_spans=((100, 150), (150, 200)),
        target_spans=((1000, 1050), (1050, 1100)),
    )

    profile = build_result_profile(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    candidate_profile = profile.candidate_profiles[0]
    assert candidate_profile.exact_mapped_segment_count == 2
    assert candidate_profile.geometric_segment_count == 1
    assert not candidate_profile.fragmented
    assert not candidate_profile.target_discontinuous
    assert profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION


def test_large_target_gap_is_complete_but_discontinuous() -> None:
    candidate = _candidate(
        "large-gap",
        source_spans=((100, 150), (150, 200)),
        target_spans=((1000, 1050), (10_001_050, 10_001_100)),
        target_gaps=((1050, 10_001_050),),
    )

    profile = build_result_profile(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    candidate_profile = profile.candidate_profiles[0]
    assert candidate_profile.exact_mapped_segment_count == 2
    assert candidate_profile.geometric_segment_count == 2
    assert candidate_profile.fragmented
    assert candidate_profile.target_discontinuous
    assert candidate_profile.largest_target_gap_bases == 10_000_000
    assert profile.headline is FactualHeadline.COMPLETE_BUT_DISCONTINUOUS_PROJECTION


def test_union_coverage_deduplicates_overlapping_candidate_source_spans() -> None:
    first = _candidate(
        "first",
        source_spans=((100, 180),),
        target_spans=((1000, 1080),),
    )
    second = _candidate(
        "second",
        source_spans=((120, 200),),
        target_spans=((2000, 2080),),
    )

    profile = build_result_profile(
        SOURCE,
        (first, second),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert profile.maximum_candidate_covered_source_bases == 80
    assert profile.union_covered_source_bases == 100
    assert (
        profile.headline
        is FactualHeadline.SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS
    )


def test_overlapping_contained_projection_does_not_imply_source_split() -> None:
    first = _candidate(
        "first",
        source_spans=((100, 180),),
        target_spans=((1000, 1080),),
    )
    second = _candidate(
        "second",
        source_spans=((120, 160),),
        target_spans=((2000, 2040),),
    )

    profile = build_result_profile(
        SOURCE,
        (first, second),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert profile.maximum_candidate_covered_source_bases == 80
    assert profile.union_covered_source_bases == 80
    assert profile.headline is FactualHeadline.MULTIPLE_CHAIN_PROJECTIONS


def test_profile_rejects_coverage_that_disagrees_with_candidate_segments() -> None:
    candidate = NormalizedCandidate(
        candidate_id="bad",
        target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1080),
        orientation=MappingOrientation.SAME,
        mapping_provenance=CHAIN,
        segments=(
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 180),
                GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1080),
            ),
        ),
        evidence=(
            EvidenceObservation(
                "bad:coverage",
                EvidenceKind.MAPPING_COVERAGE,
                MappingCoverageSummary(
                    status=MappingCoverageStatus.PARTIAL,
                    covered_source_bases=90,
                    source_bases=100,
                    uncovered_source_intervals=(
                        GenomicInterval(SOURCE_ASSEMBLY, "chr1", 190, 200),
                    ),
                ),
                CHAIN,
            ),
            EvidenceObservation(
                "bad:gaps", EvidenceKind.CHAIN_GAPS, ChainGapSummary(), CHAIN
            ),
        ),
    )

    with pytest.raises(
        ValueError, match="does not match its normalized mapping segments"
    ):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_adjacent_reverse_orientation_partition_is_canonicalized() -> None:
    whole = _candidate(
        "whole",
        orientation=MappingOrientation.REVERSE,
    )
    split = _candidate(
        "split",
        source_spans=((100, 150), (150, 200)),
        target_spans=((1050, 1100), (1000, 1050)),
        orientation=MappingOrientation.REVERSE,
    )

    with pytest.raises(ValueError, match="identical normalized mapping geometry"):
        build_result_profile(
            SOURCE,
            (whole, split),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_same_target_bounds_with_different_internal_geometry_remain_distinct() -> None:
    first = _candidate(
        "first",
        source_spans=((100, 140), (160, 200)),
        target_spans=((1000, 1040), (1060, 1100)),
        target_gaps=((1040, 1060),),
    )
    second = _candidate(
        "second",
        source_spans=((100, 130), (150, 200)),
        target_spans=((1000, 1030), (1050, 1100)),
        target_gaps=((1030, 1050),),
    )

    profile = build_result_profile(
        SOURCE,
        (first, second),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert (
        profile.headline
        is FactualHeadline.SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS
    )
    assert tuple(
        candidate.target_bounding_span for candidate in profile.candidate_profiles
    ) == (
        GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1100),
        GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1100),
    )


def test_profile_rejects_duplicate_mapping_coverage_observations() -> None:
    candidate = _candidate("candidate")
    duplicate = replace(
        candidate.evidence[0],
        observation_id="candidate:coverage:duplicate",
    )
    candidate = replace(candidate, evidence=(*candidate.evidence, duplicate))

    with pytest.raises(ValueError, match="duplicate MAPPING_COVERAGE"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_comparative_profile_requires_reciprocal_best_for_each_candidate() -> None:
    with pytest.raises(ValueError, match="missing RECIPROCAL_BEST_MEMBERSHIP"):
        build_result_profile(
            SOURCE,
            (_candidate("candidate"),),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        )


def test_liftover_only_profile_rejects_reciprocal_best_evidence() -> None:
    candidate = _candidate(
        "candidate",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )

    with pytest.raises(ValueError, match="LIFTOVER-ONLY"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_profile_rejects_candidate_geometry_outside_source_interval() -> None:
    candidate = _candidate("candidate")
    candidate = replace(
        candidate,
        segments=(
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 99, 199),
                candidate.segments[0].target_interval,
            ),
        ),
    )

    with pytest.raises(ValueError, match="outside the assessed source locus"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_profile_rejects_wrong_uncovered_interval_geometry() -> None:
    candidate = _candidate(
        "partial",
        source_spans=((100, 199),),
        target_spans=((1000, 1099),),
    )
    coverage = candidate.evidence[0]
    assert isinstance(coverage.value, MappingCoverageSummary)
    inconsistent = replace(
        coverage,
        value=MappingCoverageSummary(
            status=MappingCoverageStatus.PARTIAL,
            covered_source_bases=99,
            source_bases=100,
            uncovered_source_intervals=(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 101),
            ),
        ),
    )
    candidate = replace(candidate, evidence=(inconsistent, *candidate.evidence[1:]))

    with pytest.raises(ValueError, match="uncovered intervals"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_profile_rejects_rbest_geometry_outside_candidate_segments() -> None:
    candidate = _candidate(
        "partial",
        source_spans=((100, 199),),
        target_spans=((1000, 1099),),
        reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL,
    )
    reciprocal = candidate.evidence[-1]
    assert isinstance(reciprocal.value, ReciprocalBestMembershipSummary)
    inconsistent = replace(
        reciprocal,
        value=ReciprocalBestMembershipSummary(
            status=ReciprocalBestMembershipStatus.PARTIAL,
            resource_completeness=(
                ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
            ),
            chains_examined=1,
            covered_source_bases=98,
            candidate_source_bases=99,
            covered_source_intervals=(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 102, 200),
            ),
        ),
    )
    candidate = replace(candidate, evidence=(*candidate.evidence[:-1], inconsistent))

    with pytest.raises(ValueError, match="within normalized mapping segments"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        )


def test_rbest_interval_may_span_adjacent_source_segments() -> None:
    candidate = _candidate(
        "split",
        source_spans=((100, 150), (150, 200)),
        target_spans=((1000, 1050), (1070, 1120)),
        target_gaps=((1050, 1070),),
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    reciprocal = candidate.evidence[-1]
    assert isinstance(reciprocal.value, ReciprocalBestMembershipSummary)
    spanning = replace(
        reciprocal,
        value=ReciprocalBestMembershipSummary(
            status=ReciprocalBestMembershipStatus.FULL,
            resource_completeness=(
                ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
            ),
            chains_examined=1,
            covered_source_bases=100,
            candidate_source_bases=100,
            covered_source_intervals=(SOURCE,),
        ),
    )
    candidate = replace(candidate, evidence=(*candidate.evidence[:-1], spanning))

    profile = build_result_profile(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert profile.headline is FactualHeadline.COMPLETE_BUT_DISCONTINUOUS_PROJECTION


def test_profile_exposes_exact_source_chain_gap_geometry() -> None:
    candidate = _candidate(
        "source-gap",
        source_spans=((100, 150), (160, 200)),
        target_spans=((1000, 1050), (1050, 1090)),
    )
    gaps = EvidenceObservation(
        "source-gap:gaps",
        EvidenceKind.CHAIN_GAPS,
        ChainGapSummary(
            (
                ChainGap(
                    source_boundary=150,
                    source_gap_overlap=GenomicInterval(
                        SOURCE_ASSEMBLY, "chr1", 150, 160
                    ),
                ),
            )
        ),
        CHAIN,
    )
    candidate = replace(candidate, evidence=(candidate.evidence[0], gaps))

    profile = build_result_profile(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    candidate_profile = profile.candidate_profiles[0]
    assert candidate_profile.source_gap_intervals == (
        GenomicInterval(SOURCE_ASSEMBLY, "chr1", 150, 160),
    )
    assert candidate_profile.largest_source_gap_bases == 10


def test_reverse_mapping_geometry_must_match_forward_candidate() -> None:
    candidate = _candidate("reverse-geometry")
    reverse = reverse_mapping_unavailable(candidate)
    mismatched = replace(
        reverse,
        original_source_segments=(GenomicInterval(SOURCE_ASSEMBLY, "chr1", 101, 201),),
    )

    with pytest.raises(ValueError, match="original-source geometry"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            reverse_mapping_results=(mismatched,),
        )


def test_reverse_mapping_candidate_id_must_match_forward_candidate() -> None:
    candidate = _candidate("reverse-id")
    reverse = reverse_mapping_unavailable(candidate)
    mismatched = replace(reverse, forward_candidate_id="different-candidate")

    with pytest.raises(ValueError, match="preserve forward candidate order"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            reverse_mapping_results=(mismatched,),
        )


def test_reverse_mapping_query_geometry_must_match_forward_candidate() -> None:
    candidate = _candidate("reverse-query-geometry")
    reverse = reverse_mapping_unavailable(candidate)
    mismatched = replace(
        reverse,
        queried_target_segments=(GenomicInterval(TARGET_ASSEMBLY, "chrA", 1001, 1101),),
    )

    with pytest.raises(ValueError, match="query geometry"):
        build_result_profile(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            reverse_mapping_results=(mismatched,),
        )


def test_point_context_zero_at_both_scales_is_not_called_agreement() -> None:
    profile = build_result_profile(
        POINT,
        (),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        query_context_result=PointQueryContextResult(
            check_state=QueryContextState.RUN,
            requested_window_bases=101,
            tested_source_interval=POINT_CONTEXT,
            candidates=(),
        ),
    )

    context = profile.query_context
    assert context.findings == (QueryContextFinding.NO_PROJECTION_AT_EITHER_SCALE,)
    assert context.projection_count is ProjectionCountState.NONE
    assert not context.point_and_local_context_map_together


def test_point_context_mapped_agreement_requires_real_candidate() -> None:
    point_candidate = _candidate(
        "shared",
        source_interval=POINT,
        source_spans=((150, 151),),
        target_spans=((1050, 1051),),
    )
    context_candidate = _candidate(
        "shared",
        source_interval=POINT_CONTEXT,
        source_spans=((100, 201),),
        target_spans=((1000, 1101),),
    )
    profile = build_result_profile(
        POINT,
        (point_candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        query_context_result=PointQueryContextResult(
            check_state=QueryContextState.RUN,
            requested_window_bases=101,
            tested_source_interval=POINT_CONTEXT,
            candidates=(context_candidate,),
        ),
    )

    context = profile.query_context
    assert context.findings == (QueryContextFinding.AGREES_WITH_POINT,)
    assert context.point_and_local_context_map_together


def test_point_context_profile_keeps_structural_findings_independent() -> None:
    point_candidate = _candidate(
        "shared",
        source_interval=POINT,
        source_spans=((150, 151),),
        target_spans=((1050, 1051),),
    )
    context_candidate = _candidate(
        "shared",
        source_interval=POINT_CONTEXT,
        source_spans=((100, 160), (161, 201)),
        target_spans=((1000, 1060), (1070, 1110)),
        target_gaps=((1060, 1070),),
    )
    profile = build_result_profile(
        POINT,
        (point_candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        query_context_result=PointQueryContextResult(
            check_state=QueryContextState.RUN,
            requested_window_bases=101,
            tested_source_interval=POINT_CONTEXT,
            candidates=(context_candidate,),
        ),
    )

    assert set(profile.query_context.findings) == {
        QueryContextFinding.REVEALS_PARTIAL_COVERAGE,
        QueryContextFinding.REVEALS_FRAGMENTATION,
        QueryContextFinding.REVEALS_TARGET_DISCONTINUITY,
        QueryContextFinding.CHANGES_WITH_QUERY_SCALE,
    }
