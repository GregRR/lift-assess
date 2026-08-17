from __future__ import annotations

from liftassess import (
    AssemblyIdentifier,
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
    assess_candidates,
)
from liftassess.reporting import format_display_interval, render_assessment_summary

SOURCE_ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")
TARGET_ASSEMBLY = AssemblyIdentifier("targetAsm", "test")
SOURCE = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 200)
ALIGNMENT = ProvenanceSource("alignment", "shared test alignment")
CHAIN = ProvenanceSource("chain", "chain resource", derived_from=(ALIGNMENT,))
RBEST = ProvenanceSource("rbest", "reciprocal-best resource", derived_from=(ALIGNMENT,))


def _candidate(
    candidate_id: str,
    *,
    covered_end: int = 200,
    target_start: int = 1000,
    reciprocal_best: ReciprocalBestMembershipStatus | None = None,
) -> NormalizedCandidate:
    covered_bases = covered_end - SOURCE.start
    coverage_status = (
        MappingCoverageStatus.FULL
        if covered_end == SOURCE.end
        else MappingCoverageStatus.PARTIAL
    )
    uncovered = (
        ()
        if coverage_status is MappingCoverageStatus.FULL
        else (GenomicInterval(SOURCE_ASSEMBLY, "chr1", covered_end, SOURCE.end),)
    )
    evidence = [
        EvidenceObservation(
            observation_id=f"{candidate_id}:coverage",
            kind=EvidenceKind.MAPPING_COVERAGE,
            value=MappingCoverageSummary(
                status=coverage_status,
                covered_source_bases=covered_bases,
                source_bases=SOURCE.length,
                uncovered_source_intervals=uncovered,
            ),
            provenance=CHAIN,
        )
    ]

    if reciprocal_best is not None:
        reciprocal_intervals: tuple[GenomicInterval, ...]
        if reciprocal_best is ReciprocalBestMembershipStatus.FULL:
            reciprocal_covered = covered_bases
            reciprocal_intervals = (
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", SOURCE.start, covered_end),
            )
        elif reciprocal_best is ReciprocalBestMembershipStatus.PARTIAL:
            reciprocal_covered = covered_bases - 1
            reciprocal_intervals = (
                GenomicInterval(
                    SOURCE_ASSEMBLY,
                    "chr1",
                    SOURCE.start,
                    SOURCE.start + reciprocal_covered,
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
                    resource_completeness=(
                        ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
                    ),
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
                    SOURCE_ASSEMBLY, "chr1", SOURCE.start, covered_end
                ),
                target_interval=GenomicInterval(
                    TARGET_ASSEMBLY, "chrA", target_start, target_end
                ),
            ),
        ),
        evidence=tuple(evidence),
    )


def test_format_display_interval_makes_coordinate_convention_explicit() -> None:
    interval = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 12344, 12400)

    assert format_display_interval(interval) == "chr1:12345-12400 (1-based inclusive)"


def test_liftover_only_well_supported_summary_is_concise_and_explicit() -> None:
    assessment = assess_candidates(
        SOURCE,
        (_candidate("only"),),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    summary = render_assessment_summary(assessment)

    assert summary.splitlines()[:4] == [
        "Source locus: chr1:101-200 (1-based inclusive)",
        "Evidence availability: LIFTOVER-ONLY — chain mapping evidence only",
        "Assessment: WELL SUPPORTED",
        "Preferred candidate: chrA:1001-1100 (1-based inclusive; same orientation)",
    ]
    assert "Why: full source-locus mapping coverage" in summary
    assert summary.endswith("This does not establish biological correctness.")


def test_comparative_well_supported_summary_names_both_verdict_driving_states() -> None:
    assessment = assess_candidates(
        SOURCE,
        (_candidate("only", reciprocal_best=ReciprocalBestMembershipStatus.FULL),),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    summary = render_assessment_summary(assessment)

    assert (
        "Evidence availability: COMPARATIVE — mapping plus comparative evidence "
        "available" in summary
    )
    assert (
        "Why: full source-locus mapping coverage and full reciprocal-best membership"
        in summary
    )


def test_contested_multi_candidate_summary_does_not_imply_candidate_ranking() -> None:
    assessment = assess_candidates(
        SOURCE,
        (_candidate("first"), _candidate("second", target_start=2000)),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    summary = render_assessment_summary(assessment)

    assert "Assessment: CONTESTED" in summary
    assert "Candidates assessed: 2" in summary
    assert "Preferred candidate:" not in summary
    assert "chrA:1001-1100" not in summary
    assert "Why: multiple candidates retain material assessment evidence" in summary


def test_indeterminate_zero_candidate_summary_states_why() -> None:
    assessment = assess_candidates(
        SOURCE,
        (),
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    summary = render_assessment_summary(assessment)

    assert "Assessment: INDETERMINATE" in summary
    assert "Candidates assessed: 0" in summary
    assert "Why: no candidate mapping was generated for the requested locus" in summary


def test_indeterminate_partial_candidate_summary_reports_partial_mapping() -> None:
    assessment = assess_candidates(
        SOURCE,
        (_candidate("partial", covered_end=190),),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    summary = render_assessment_summary(assessment)

    assert "Assessment: INDETERMINATE" in summary
    assert "Why: candidate maps only part of the requested source locus" in summary


def test_split_candidate_summary_labels_target_interval_as_bounding_span() -> None:
    first_source = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 150)
    second_source = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 150, 200)
    first_target = GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1050)
    second_target = GenomicInterval(TARGET_ASSEMBLY, "chrA", 1100, 1150)
    candidate = NormalizedCandidate(
        candidate_id="split",
        target_interval=GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1150),
        orientation=MappingOrientation.SAME,
        mapping_provenance=CHAIN,
        segments=(
            MappingSegment(first_source, first_target),
            MappingSegment(second_source, second_target),
        ),
        evidence=(
            EvidenceObservation(
                observation_id="split:coverage",
                kind=EvidenceKind.MAPPING_COVERAGE,
                value=MappingCoverageSummary(
                    status=MappingCoverageStatus.FULL,
                    covered_source_bases=100,
                    source_bases=100,
                ),
                provenance=CHAIN,
            ),
        ),
    )
    assessment = assess_candidates(
        SOURCE,
        (candidate,),
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )

    summary = render_assessment_summary(assessment)

    assert (
        "Preferred candidate: chrA:1001-1150 (1-based inclusive; same orientation; "
        "bounding span of 2 mapped segments)" in summary
    )
