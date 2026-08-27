from __future__ import annotations

from dataclasses import replace

import pytest

from liftassess import (
    AssemblyIdentifier,
    ComparativeEvidenceRelationship,
    ComparativeEvidenceRelationshipResult,
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    MappingSegment,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
    build_comparative_evidence_relationship,
    build_filtered_all_chain_comparison,
)

SOURCE_ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")
TARGET_ASSEMBLY = AssemblyIdentifier("targetAsm", "test")
SOURCE_INTERVAL = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 110)
ALIGNMENT = ProvenanceSource("alignment", "shared alignment lineage")
ALL_CHAIN = ProvenanceSource("all-chain", "all-chain bytes", derived_from=(ALIGNMENT,))
FILTERED_CHAIN = ProvenanceSource(
    "filtered-chain",
    "filtered chain bytes",
    derived_from=(ALIGNMENT,),
)
NET = ProvenanceSource("net", "net bytes", derived_from=(ALIGNMENT,))
RBEST = ProvenanceSource("rbest", "reciprocal-best bytes", derived_from=(ALIGNMENT,))


def _candidate(
    candidate_id: str,
    *,
    target_sequence: str,
    target_start: int,
    full: bool = True,
    top_net: bool = False,
    full_rbest: bool = False,
    split_top_fill_evidence: bool = False,
) -> NormalizedCandidate:
    source_end = 110 if full else 109
    segment = MappingSegment(
        GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, source_end),
        GenomicInterval(
            TARGET_ASSEMBLY,
            target_sequence,
            target_start,
            target_start + (source_end - 100),
        ),
    )
    mapping = MappingCoverageSummary(
        status=MappingCoverageStatus.FULL if full else MappingCoverageStatus.PARTIAL,
        covered_source_bases=source_end - 100,
        source_bases=SOURCE_INTERVAL.length,
        uncovered_source_intervals=(
            () if full else (GenomicInterval(SOURCE_ASSEMBLY, "chr1", 109, 110),)
        ),
    )
    evidence: list[EvidenceObservation] = [
        EvidenceObservation(
            f"{candidate_id}:coverage",
            EvidenceKind.MAPPING_COVERAGE,
            mapping,
            ALL_CHAIN,
        )
    ]

    if top_net:
        classification_fill = ProvenanceSource(
            f"{candidate_id}:fill:classification",
            "top net fill",
            derived_from=(NET,),
        )
        hierarchy_fill = (
            ProvenanceSource(
                f"{candidate_id}:fill:hierarchy",
                "different net fill",
                derived_from=(NET,),
            )
            if split_top_fill_evidence
            else classification_fill
        )
        evidence.extend(
            (
                EvidenceObservation(
                    f"{candidate_id}:net:classification",
                    EvidenceKind.NET_CLASSIFICATION,
                    "top",
                    classification_fill,
                ),
                EvidenceObservation(
                    f"{candidate_id}:net:hierarchy",
                    EvidenceKind.NET_HIERARCHY,
                    NetHierarchySummary(
                        depth=1,
                        source_fill_interval=SOURCE_INTERVAL,
                    ),
                    hierarchy_fill,
                ),
            )
        )

    evidence.append(
        EvidenceObservation(
            f"{candidate_id}:rbest",
            EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
            ReciprocalBestMembershipSummary(
                status=(
                    ReciprocalBestMembershipStatus.FULL
                    if full_rbest
                    else ReciprocalBestMembershipStatus.NONE
                ),
                resource_completeness=(
                    ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
                ),
                chains_examined=1,
                covered_source_bases=source_end - 100 if full_rbest else 0,
                candidate_source_bases=source_end - 100,
                covered_source_intervals=(
                    (segment.source_interval,) if full_rbest else ()
                ),
            ),
            RBEST,
        )
    )

    return NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=segment.target_interval,
        orientation=MappingOrientation.SAME,
        mapping_provenance=ALL_CHAIN,
        segments=(segment,),
        evidence=tuple(evidence),
    )


def _filtered_for(
    candidate: NormalizedCandidate, candidate_id: str
) -> NormalizedCandidate:
    return replace(
        candidate,
        candidate_id=candidate_id,
        mapping_provenance=FILTERED_CHAIN,
        evidence=(),
    )


def _relationship(
    all_chain_candidates: tuple[NormalizedCandidate, ...],
    *,
    retained: NormalizedCandidate | None,
) -> ComparativeEvidenceRelationshipResult:
    filtered = () if retained is None else (_filtered_for(retained, "filtered"),)
    comparison = build_filtered_all_chain_comparison(
        SOURCE_INTERVAL,
        all_chain_candidates,
        filtered,
        all_chain_provenance=ALL_CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )
    return build_comparative_evidence_relationship(comparison)


def test_b14_style_pattern_favors_the_uniquely_supported_filtered_placement() -> None:
    retained = _candidate(
        "retained",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
        full_rbest=True,
    )
    extra_one = _candidate("extra-one", target_sequence="chr25", target_start=700)
    extra_two = _candidate("extra-two", target_sequence="chr13", target_start=900)

    result = _relationship((retained, extra_one, extra_two), retained=retained)

    assert result.relationship is ComparativeEvidenceRelationship.FAVORS_ONE_PLACEMENT
    assert result.favored_candidate_id == "retained"
    assert result.full_candidate_ids == ("retained", "extra-one", "extra-two")
    assert result.filtered_retained_full_candidate_ids == ("retained",)
    assert result.depth1_top_net_full_candidate_ids == ("retained",)
    assert result.full_rbest_full_candidate_ids == ("retained",)
    assert result.joint_top_net_full_rbest_candidate_ids == ("retained",)


def test_partial_support_does_not_become_hidden_weighting() -> None:
    retained = _candidate(
        "retained",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
    )
    competitor = _candidate("competitor", target_sequence="chr25", target_start=700)

    result = _relationship((retained, competitor), retained=retained)

    assert (
        result.relationship
        is ComparativeEvidenceRelationship.DOES_NOT_SEPARATE_PLACEMENTS
    )
    assert result.filtered_retained_full_candidate_ids == ("retained",)
    assert result.depth1_top_net_full_candidate_ids == ("retained",)
    assert result.full_rbest_full_candidate_ids == ()
    assert result.favored_candidate_id is None


def test_equivalent_joint_support_on_competitor_does_not_separate_placements() -> None:
    retained = _candidate(
        "retained",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
        full_rbest=True,
    )
    competitor = _candidate(
        "competitor",
        target_sequence="chr25",
        target_start=700,
        top_net=True,
        full_rbest=True,
    )

    result = _relationship((retained, competitor), retained=retained)

    assert (
        result.relationship
        is ComparativeEvidenceRelationship.DOES_NOT_SEPARATE_PLACEMENTS
    )
    assert result.favored_candidate_id is None
    assert result.joint_top_net_full_rbest_candidate_ids == (
        "retained",
        "competitor",
    )


def test_multiple_filtered_retained_placements_do_not_favor_one() -> None:
    retained_full = _candidate(
        "retained-full",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
        full_rbest=True,
    )
    retained_partial = _candidate(
        "retained-partial",
        target_sequence="chr7",
        target_start=600,
        full=False,
    )
    competitor = _candidate(
        "competitor",
        target_sequence="chr25",
        target_start=700,
    )
    comparison = build_filtered_all_chain_comparison(
        SOURCE_INTERVAL,
        (retained_full, retained_partial, competitor),
        (
            _filtered_for(retained_full, "filtered-full"),
            _filtered_for(retained_partial, "filtered-partial"),
        ),
        all_chain_provenance=ALL_CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )

    result = build_comparative_evidence_relationship(comparison)

    assert (
        result.relationship
        is ComparativeEvidenceRelationship.DOES_NOT_SEPARATE_PLACEMENTS
    )
    assert result.filtered_retained_full_candidate_ids == ("retained-full",)
    assert result.joint_top_net_full_rbest_candidate_ids == ("retained-full",)
    assert result.favored_candidate_id is None


def test_unique_joint_support_conflicting_with_filtered_retention_is_mixed() -> None:
    retained = _candidate("retained", target_sequence="chr5", target_start=500)
    competitor = _candidate(
        "competitor",
        target_sequence="chr25",
        target_start=700,
        top_net=True,
        full_rbest=True,
    )

    result = _relationship((retained, competitor), retained=retained)

    assert result.relationship is ComparativeEvidenceRelationship.MIXED_CONFLICTING
    assert result.favored_candidate_id is None
    assert result.filtered_retained_full_candidate_ids == ("retained",)
    assert result.joint_top_net_full_rbest_candidate_ids == ("competitor",)


def test_unique_top_net_and_full_rbest_on_different_placements_is_mixed() -> None:
    top = _candidate(
        "top",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
    )
    rbest = _candidate(
        "rbest",
        target_sequence="chr25",
        target_start=700,
        full_rbest=True,
    )

    result = _relationship((top, rbest), retained=None)

    assert result.relationship is ComparativeEvidenceRelationship.MIXED_CONFLICTING
    assert result.depth1_top_net_full_candidate_ids == ("top",)
    assert result.full_rbest_full_candidate_ids == ("rbest",)


def test_single_full_placement_is_not_forced_into_a_separation_category() -> None:
    retained = _candidate(
        "retained",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
        full_rbest=True,
    )

    result = _relationship((retained,), retained=retained)

    assert (
        result.relationship
        is ComparativeEvidenceRelationship.NO_COMPETING_FULL_PLACEMENTS
    )
    assert result.favored_candidate_id is None


def test_partial_extra_placement_is_not_a_competing_full_placement() -> None:
    retained = _candidate(
        "retained",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
        full_rbest=True,
    )
    partial = _candidate(
        "partial",
        target_sequence="chr25",
        target_start=700,
        full=False,
        top_net=True,
        full_rbest=True,
    )

    result = _relationship((retained, partial), retained=retained)

    assert (
        result.relationship
        is ComparativeEvidenceRelationship.NO_COMPETING_FULL_PLACEMENTS
    )
    assert result.full_candidate_ids == ("retained",)


def test_top_classification_and_depth_must_come_from_the_same_net_fill() -> None:
    retained = _candidate(
        "retained",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
        full_rbest=True,
        split_top_fill_evidence=True,
    )
    competitor = _candidate("competitor", target_sequence="chr25", target_start=700)

    result = _relationship((retained, competitor), retained=retained)

    assert (
        result.relationship
        is ComparativeEvidenceRelationship.DOES_NOT_SEPARATE_PLACEMENTS
    )
    assert result.depth1_top_net_full_candidate_ids == ()
    assert result.joint_top_net_full_rbest_candidate_ids == ()


def test_ali_qdup_and_chain_score_do_not_change_the_categorical_rule() -> None:
    retained = _candidate(
        "retained",
        target_sequence="chr5",
        target_start=500,
        top_net=True,
        full_rbest=True,
    )
    competitor = _candidate("competitor", target_sequence="chr25", target_start=700)
    noisy_competitor = replace(
        competitor,
        evidence=competitor.evidence
        + (
            EvidenceObservation("ali", EvidenceKind.ALIGNED_BASES, 10**12, NET),
            EvidenceObservation(
                "qdup",
                EvidenceKind.DUPLICATED_QUERY_BASES,
                0,
                NET,
            ),
            EvidenceObservation("score", EvidenceKind.CHAIN_SCORE, 10**15, ALL_CHAIN),
        ),
    )

    result = _relationship((retained, noisy_competitor), retained=retained)

    assert result.relationship is ComparativeEvidenceRelationship.FAVORS_ONE_PLACEMENT
    assert result.favored_candidate_id == "retained"


def test_missing_reciprocal_best_observation_is_an_invariant_failure() -> None:
    retained = _candidate("retained", target_sequence="chr5", target_start=500)
    invalid = replace(
        retained,
        evidence=tuple(
            observation
            for observation in retained.evidence
            if observation.kind is not EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
        ),
    )
    comparison = build_filtered_all_chain_comparison(
        SOURCE_INTERVAL,
        (invalid,),
        (_filtered_for(invalid, "filtered"),),
        all_chain_provenance=ALL_CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )

    with pytest.raises(ValueError, match="exactly one typed reciprocal-best"):
        build_comparative_evidence_relationship(comparison)
