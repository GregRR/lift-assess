from __future__ import annotations

from dataclasses import replace

import pytest

from liftassess import (
    AssemblyIdentifier,
    AssessmentDecisionReason,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
    Verdict,
    assess_candidates,
)

SOURCE_ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")
TARGET_ASSEMBLY = AssemblyIdentifier("targetAsm", "test")
SOURCE = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 200)
ALIGNMENT = ProvenanceSource("alignment", "shared test alignment")
CHAIN = ProvenanceSource("chain", "chain file", derived_from=(ALIGNMENT,))
RBEST = ProvenanceSource("rbest", "reciprocal-best file", derived_from=(ALIGNMENT,))


def _candidate(
    candidate_id: str,
    *,
    covered_start: int = 100,
    covered_end: int = 200,
    target_start: int = 1000,
    reciprocal_best: ReciprocalBestMembershipStatus | None = None,
    reciprocal_best_completeness: ReciprocalBestResourceCompleteness = (
        ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
    ),
    extra_evidence: tuple[EvidenceObservation, ...] = (),
) -> NormalizedCandidate:
    covered_bases = covered_end - covered_start
    status = (
        MappingCoverageStatus.FULL
        if covered_start == SOURCE.start and covered_end == SOURCE.end
        else MappingCoverageStatus.PARTIAL
    )
    uncovered: list[GenomicInterval] = []
    if SOURCE.start < covered_start:
        uncovered.append(
            GenomicInterval(SOURCE_ASSEMBLY, "chr1", SOURCE.start, covered_start)
        )
    if covered_end < SOURCE.end:
        uncovered.append(
            GenomicInterval(SOURCE_ASSEMBLY, "chr1", covered_end, SOURCE.end)
        )

    coverage = EvidenceObservation(
        observation_id=f"{candidate_id}:coverage",
        kind=EvidenceKind.MAPPING_COVERAGE,
        value=MappingCoverageSummary(
            status=status,
            covered_source_bases=covered_bases,
            source_bases=SOURCE.length,
            uncovered_source_intervals=tuple(uncovered),
        ),
        provenance=CHAIN,
    )
    evidence: list[EvidenceObservation] = [coverage, *extra_evidence]

    if reciprocal_best is not None:
        reciprocal_intervals: tuple[GenomicInterval, ...]
        if reciprocal_best is ReciprocalBestMembershipStatus.FULL:
            reciprocal_covered = covered_bases
            reciprocal_intervals = (
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", covered_start, covered_end),
            )
        elif reciprocal_best is ReciprocalBestMembershipStatus.PARTIAL:
            reciprocal_covered = max(1, covered_bases - 1)
            reciprocal_intervals = (
                GenomicInterval(
                    SOURCE_ASSEMBLY,
                    "chr1",
                    covered_start,
                    covered_start + reciprocal_covered,
                ),
            )
        else:
            reciprocal_covered = 0
            reciprocal_intervals = ()

        evidence.append(
            EvidenceObservation(
                observation_id=f"{candidate_id}:rbest",
                kind=EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
                value=ReciprocalBestMembershipSummary(
                    status=reciprocal_best,
                    resource_completeness=reciprocal_best_completeness,
                    chains_examined=1,
                    covered_source_bases=reciprocal_covered,
                    candidate_source_bases=covered_bases,
                    covered_source_intervals=reciprocal_intervals,
                ),
                provenance=RBEST,
            )
        )

    target_end = target_start + covered_bases
    return NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=GenomicInterval(
            TARGET_ASSEMBLY, "chrA", target_start, target_end
        ),
        orientation=MappingOrientation.SAME,
        mapping_provenance=CHAIN,
        segments=(
            MappingSegment(
                source_interval=GenomicInterval(
                    SOURCE_ASSEMBLY, "chr1", covered_start, covered_end
                ),
                target_interval=GenomicInterval(
                    TARGET_ASSEMBLY, "chrA", target_start, target_end
                ),
            ),
        ),
        evidence=tuple(evidence),
    )


def test_no_candidates_is_indeterminate_for_either_evidence_tier() -> None:
    for tier in EvidenceAvailabilityTier:
        assessment = assess_candidates(SOURCE, (), evidence_tier=tier)
        assert assessment.verdict is Verdict.INDETERMINATE
        assert assessment.decision_reason is AssessmentDecisionReason.NO_CANDIDATES
        assert assessment.preferred_candidate_id is None


def test_assessor_terminal_cases_cover_complete_decision_reason_vocabulary() -> None:
    liftover_full = _candidate("lift-full")
    liftover_partial = _candidate("lift-partial", covered_end=199)
    comparative_full_full = _candidate(
        "comp-full-full",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    comparative_full_partial = _candidate(
        "comp-full-partial",
        reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL,
    )
    comparative_full_none = _candidate(
        "comp-full-none",
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )
    comparative_partial_full = _candidate(
        "comp-partial-full",
        covered_end=199,
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    comparative_partial_none = _candidate(
        "comp-partial-none",
        covered_end=199,
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )

    assessments = (
        assess_candidates(
            SOURCE,
            (),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        ),
        assess_candidates(
            SOURCE,
            (liftover_full, _candidate("lift-other", target_start=2000)),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        ),
        assess_candidates(
            SOURCE,
            (liftover_full,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        ),
        assess_candidates(
            SOURCE,
            (liftover_partial,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        ),
        assess_candidates(
            SOURCE,
            (
                comparative_full_full,
                _candidate(
                    "comp-other-material",
                    target_start=2000,
                    reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL,
                ),
            ),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        ),
        assess_candidates(
            SOURCE,
            (comparative_full_full,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        ),
        assess_candidates(
            SOURCE,
            (comparative_full_none,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        ),
        assess_candidates(
            SOURCE,
            (comparative_full_partial,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        ),
        assess_candidates(
            SOURCE,
            (comparative_partial_full,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        ),
        assess_candidates(
            SOURCE,
            (comparative_partial_none,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        ),
    )

    assert {assessment.decision_reason for assessment in assessments} == set(
        AssessmentDecisionReason
    )


def test_liftover_only_single_full_candidate_is_well_supported() -> None:
    candidate = _candidate("full")

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert assessment.verdict is Verdict.WELL_SUPPORTED
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.LIFTOVER_SINGLE_FULL_MAPPING
    )
    assert assessment.preferred_candidate_id == "full"
    assert assessment.supporting_evidence[0].observation_id == "full:coverage"


def test_liftover_only_single_partial_candidate_is_indeterminate() -> None:
    candidate = _candidate("partial", covered_end=199)

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert assessment.verdict is Verdict.INDETERMINATE
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.LIFTOVER_SINGLE_PARTIAL_MAPPING
    )
    assert assessment.preferred_candidate_id is None
    assert assessment.contradicting_evidence[0].observation_id == "partial:coverage"


def test_liftover_only_multiple_candidates_are_contested_without_ranking() -> None:
    full = _candidate("full")
    almost_full = _candidate("almost-full", covered_end=199, target_start=2000)

    assessment = assess_candidates(
        SOURCE,
        (full, almost_full),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert assessment.verdict is Verdict.CONTESTED
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.LIFTOVER_MULTIPLE_CANDIDATES
    )
    assert assessment.preferred_candidate_id is None
    assert {reference.candidate_id for reference in assessment.supporting_evidence} == {
        "full"
    }
    contradicting_candidate_ids = {
        reference.candidate_id for reference in assessment.contradicting_evidence
    }
    assert contradicting_candidate_ids == {"almost-full"}


def test_comparative_unique_full_and_full_rbest_candidate_is_well_supported() -> None:
    candidate = _candidate(
        "retained",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.WELL_SUPPORTED
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_FULL
    )
    assert assessment.preferred_candidate_id == "retained"
    supporting_ids = {
        reference.observation_id for reference in assessment.supporting_evidence
    }
    assert supporting_ids == {
        "retained:coverage",
        "retained:rbest",
    }


def test_comparative_single_full_mapping_with_partial_rbest_is_indeterminate() -> None:
    candidate = _candidate(
        "candidate", reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.INDETERMINATE
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_PARTIAL
    )
    assert assessment.preferred_candidate_id is None
    assert {
        reference.observation_id for reference in assessment.supporting_evidence
    } == {"candidate:coverage", "candidate:rbest"}
    assert {
        reference.observation_id for reference in assessment.contradicting_evidence
    } == {"candidate:rbest"}


def test_comparative_single_full_mapping_with_no_rbest_is_contested() -> None:
    candidate = _candidate(
        "candidate", reciprocal_best=ReciprocalBestMembershipStatus.NONE
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.CONTESTED
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_NONE
    )
    assert assessment.preferred_candidate_id is None
    assert {
        reference.observation_id for reference in assessment.supporting_evidence
    } == {"candidate:coverage"}
    assert {
        reference.observation_id for reference in assessment.contradicting_evidence
    } == {"candidate:rbest"}


def test_comparative_two_material_candidates_are_contested() -> None:
    retained = _candidate(
        "retained",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    competing_full = _candidate(
        "competing",
        target_start=2000,
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )

    assessment = assess_candidates(
        SOURCE,
        (retained, competing_full),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.CONTESTED
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.COMPARATIVE_MULTIPLE_MATERIAL_CANDIDATES
    )
    assert assessment.preferred_candidate_id is None


def test_comparative_two_fully_retained_candidates_are_contested() -> None:
    first = _candidate(
        "first",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    second = _candidate(
        "second",
        target_start=2000,
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )

    assessment = assess_candidates(
        SOURCE,
        (first, second),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.CONTESTED
    assert assessment.preferred_candidate_id is None


def test_comparative_partial_candidate_with_rbest_survival_is_material() -> None:
    retained = _candidate(
        "retained",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    partial = _candidate(
        "partial",
        covered_end=199,
        target_start=2000,
        reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL,
    )

    assessment = assess_candidates(
        SOURCE,
        (retained, partial),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.CONTESTED
    assert assessment.preferred_candidate_id is None


def test_comparative_partial_rbest_none_does_not_block_unique_retained() -> None:
    retained = _candidate(
        "retained",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    partial_rejected = _candidate(
        "partial-rejected",
        covered_end=199,
        target_start=2000,
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )

    assessment = assess_candidates(
        SOURCE,
        (retained, partial_rejected),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.WELL_SUPPORTED
    assert assessment.preferred_candidate_id == "retained"


def test_comparative_partial_only_is_indeterminate_even_with_full_rbest() -> None:
    candidate = _candidate(
        "partial",
        covered_end=199,
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.INDETERMINATE
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_PARTIAL
    )
    assert assessment.preferred_candidate_id is None


def test_comparative_partial_only_with_partial_rbest_is_indeterminate() -> None:
    candidate = _candidate(
        "partial",
        covered_end=199,
        reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL,
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.INDETERMINATE
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_PARTIAL
    )
    assert assessment.preferred_candidate_id is None


def test_comparative_partial_only_with_no_rbest_is_no_material_candidate() -> None:
    candidate = _candidate(
        "partial",
        covered_end=199,
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.INDETERMINATE
    assert (
        assessment.decision_reason
        is AssessmentDecisionReason.COMPARATIVE_NO_MATERIAL_CANDIDATE
    )


def test_three_material_candidates_are_contested() -> None:
    candidates = (
        _candidate(
            "first",
            target_start=1000,
            reciprocal_best=ReciprocalBestMembershipStatus.FULL,
        ),
        _candidate(
            "second",
            target_start=2000,
            reciprocal_best=ReciprocalBestMembershipStatus.PARTIAL,
        ),
        _candidate(
            "third",
            covered_end=199,
            target_start=3000,
            reciprocal_best=ReciprocalBestMembershipStatus.FULL,
        ),
    )

    assessment = assess_candidates(
        SOURCE,
        candidates,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.CONTESTED
    assert assessment.preferred_candidate_id is None


def test_distinct_candidate_ids_with_identical_geometry_are_rejected() -> None:
    first = _candidate("first")
    second = _candidate("second")

    with pytest.raises(ValueError, match="same canonical local mapping geometry"):
        assess_candidates(
            SOURCE,
            (first, second),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_adjacent_same_orientation_segment_partition_is_canonicalized() -> None:
    single = _candidate("single")
    split = replace(
        _candidate("split"),
        segments=(
            MappingSegment(
                source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 150),
                target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1050),
            ),
            MappingSegment(
                source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 150, 200),
                target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1050, 1100),
            ),
        ),
    )

    with pytest.raises(ValueError, match="same canonical local mapping geometry"):
        assess_candidates(
            SOURCE,
            (single, split),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_adjacent_reverse_orientation_segment_partition_is_canonicalized() -> None:
    single = replace(
        _candidate("single"),
        orientation=MappingOrientation.REVERSE,
    )
    split = replace(
        _candidate("split"),
        orientation=MappingOrientation.REVERSE,
        segments=(
            MappingSegment(
                source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 150),
                target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1050, 1100),
            ),
            MappingSegment(
                source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 150, 200),
                target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1050),
            ),
        ),
    )

    with pytest.raises(ValueError, match="same canonical local mapping geometry"):
        assess_candidates(
            SOURCE,
            (single, split),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_same_target_bounds_with_different_internal_geometry_remain_distinct() -> None:
    first_uncovered = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 140, 160)
    second_uncovered = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 130, 150)

    first = NormalizedCandidate(
        candidate_id="first",
        target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1100),
        orientation=MappingOrientation.SAME,
        mapping_provenance=CHAIN,
        segments=(
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 140),
                GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1040),
            ),
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 160, 200),
                GenomicInterval(TARGET_ASSEMBLY, "chrA", 1060, 1100),
            ),
        ),
        evidence=(
            EvidenceObservation(
                observation_id="first:coverage",
                kind=EvidenceKind.MAPPING_COVERAGE,
                value=MappingCoverageSummary(
                    status=MappingCoverageStatus.PARTIAL,
                    covered_source_bases=80,
                    source_bases=100,
                    uncovered_source_intervals=(first_uncovered,),
                ),
                provenance=CHAIN,
            ),
        ),
    )
    second = NormalizedCandidate(
        candidate_id="second",
        target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1100),
        orientation=MappingOrientation.SAME,
        mapping_provenance=CHAIN,
        segments=(
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 130),
                GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1030),
            ),
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 150, 200),
                GenomicInterval(TARGET_ASSEMBLY, "chrA", 1050, 1100),
            ),
        ),
        evidence=(
            EvidenceObservation(
                observation_id="second:coverage",
                kind=EvidenceKind.MAPPING_COVERAGE,
                value=MappingCoverageSummary(
                    status=MappingCoverageStatus.PARTIAL,
                    covered_source_bases=80,
                    source_bases=100,
                    uncovered_source_intervals=(second_uncovered,),
                ),
                provenance=CHAIN,
            ),
        ),
    )

    assessment = assess_candidates(
        SOURCE,
        (first, second),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert assessment.verdict is Verdict.CONTESTED
    assert assessment.supporting_evidence == ()
    assert {
        reference.candidate_id for reference in assessment.contradicting_evidence
    } == {"first", "second"}


@pytest.mark.parametrize(
    "completeness",
    [
        ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        ReciprocalBestResourceCompleteness.COMPLETE_CANDIDATE_SUBSET,
    ],
)
def test_exhaustive_rbest_completeness_bases_give_same_none_semantics(
    completeness: ReciprocalBestResourceCompleteness,
) -> None:
    candidate = _candidate(
        "candidate",
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
        reciprocal_best_completeness=completeness,
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.CONTESTED


def test_non_verdict_evidence_cannot_change_verdict() -> None:
    noisy_context = (
        EvidenceObservation(
            observation_id="context:score",
            kind=EvidenceKind.CHAIN_SCORE,
            value=999999999.0,
            provenance=CHAIN,
        ),
        EvidenceObservation(
            observation_id="context:qdup",
            kind=EvidenceKind.DUPLICATED_QUERY_BASES,
            value=999999999,
            provenance=CHAIN,
        ),
        EvidenceObservation(
            observation_id="context:classification",
            kind=EvidenceKind.NET_CLASSIFICATION,
            value="nonSyn",
            provenance=CHAIN,
        ),
    )
    candidate = _candidate("candidate", extra_evidence=noisy_context)

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    assert assessment.verdict is Verdict.WELL_SUPPORTED


def test_duplicate_verdict_driving_observation_is_rejected_instead_of_counted() -> None:
    candidate = _candidate("candidate")
    duplicate = replace(
        candidate.evidence[0],
        observation_id="candidate:coverage:duplicate",
    )
    candidate = replace(candidate, evidence=(*candidate.evidence, duplicate))

    with pytest.raises(ValueError, match="multiple MAPPING_COVERAGE"):
        assess_candidates(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_comparative_requires_reciprocal_best_for_every_candidate() -> None:
    candidate = _candidate("candidate")

    with pytest.raises(ValueError, match="missing required RECIPROCAL_BEST_MEMBERSHIP"):
        assess_candidates(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        )


def test_liftover_only_rejects_reciprocal_best_evidence() -> None:
    candidate = _candidate(
        "candidate",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )

    with pytest.raises(ValueError, match="LIFTOVER_ONLY"):
        assess_candidates(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_assessor_rejects_coverage_that_disagrees_with_candidate_segments() -> None:
    candidate = _candidate("candidate")
    coverage = candidate.evidence[0]
    assert isinstance(coverage.value, MappingCoverageSummary)
    inconsistent = replace(
        coverage,
        value=MappingCoverageSummary(
            status=MappingCoverageStatus.PARTIAL,
            covered_source_bases=99,
            source_bases=100,
            uncovered_source_intervals=(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 199, 200),
            ),
        ),
    )
    candidate = replace(candidate, evidence=(inconsistent, *candidate.evidence[1:]))

    with pytest.raises(ValueError, match="normalized mapping segments"):
        assess_candidates(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_assessor_rejects_candidate_geometry_outside_source_locus() -> None:
    candidate = _candidate("candidate")
    outside_source = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 99, 199)
    replacement_segment = MappingSegment(
        source_interval=outside_source,
        target_interval=candidate.segments[0].target_interval,
    )
    candidate = replace(candidate, segments=(replacement_segment,))

    with pytest.raises(ValueError, match="outside the assessed source interval"):
        assess_candidates(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_assessor_rejects_wrong_uncovered_interval_geometry() -> None:
    candidate = _candidate("partial", covered_end=199)
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
        assess_candidates(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )


def test_assessor_rejects_rbest_geometry_outside_candidate_segments() -> None:
    candidate = _candidate(
        "partial",
        covered_end=199,
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
        assess_candidates(
            SOURCE,
            (candidate,),
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        )


def test_candidate_order_does_not_change_comparative_verdict() -> None:
    retained = _candidate(
        "retained",
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    competing = _candidate(
        "competing",
        target_start=2000,
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )

    forward = assess_candidates(
        SOURCE,
        (retained, competing),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )
    reversed_order = assess_candidates(
        SOURCE,
        (competing, retained),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert forward.verdict is reversed_order.verdict is Verdict.CONTESTED
    assert forward.preferred_candidate_id is None
    assert reversed_order.preferred_candidate_id is None


def test_rbest_interval_may_span_adjacent_source_segments() -> None:
    coverage = EvidenceObservation(
        observation_id="split:coverage",
        kind=EvidenceKind.MAPPING_COVERAGE,
        value=MappingCoverageSummary(
            status=MappingCoverageStatus.FULL,
            covered_source_bases=100,
            source_bases=100,
        ),
        provenance=CHAIN,
    )
    reciprocal = EvidenceObservation(
        observation_id="split:rbest",
        kind=EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
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
        provenance=RBEST,
    )
    candidate = NormalizedCandidate(
        candidate_id="split",
        target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1120),
        orientation=MappingOrientation.SAME,
        mapping_provenance=CHAIN,
        segments=(
            MappingSegment(
                source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 150),
                target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1050),
            ),
            MappingSegment(
                source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 150, 200),
                target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1070, 1120),
            ),
        ),
        evidence=(coverage, reciprocal),
    )

    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    assert assessment.verdict is Verdict.WELL_SUPPORTED
    assert assessment.preferred_candidate_id == "split"
