import pytest

from liftassess import (
    AssemblyIdentifier,
    Assessment,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    EvidenceReference,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceIdentifier,
    ProvenanceIdentifierKind,
    ProvenanceSource,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
    Verdict,
)


@pytest.fixture
def source_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier(
        name="sourceAsm",
        provider="test-provider",
        accession="TEST_SOURCE_1",
        aliases=("source-assembly-alias",),
    )


@pytest.fixture
def target_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier(
        name="targetAsm",
        provider="test-provider",
        accession="TEST_TARGET_1",
        aliases=("target-assembly-alias",),
    )


def _segment(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    *,
    source_start: int = 10,
    target_start: int = 100,
    length: int = 50,
) -> MappingSegment:
    return MappingSegment(
        source_interval=GenomicInterval(
            source_assembly, "chr16", source_start, source_start + length
        ),
        target_interval=GenomicInterval(
            target_assembly, "chr16", target_start, target_start + length
        ),
    )


def test_interval_is_zero_based_and_half_open(
    source_assembly: AssemblyIdentifier,
) -> None:
    interval = GenomicInterval(
        assembly=source_assembly,
        sequence_name="chr16",
        start=0,
        end=10,
    )

    assert interval.length == 10
    assert interval.contains(0)
    assert interval.contains(9)
    assert not interval.contains(10)


def test_interval_rejects_negative_start_and_reversed_bounds(
    source_assembly: AssemblyIdentifier,
) -> None:
    with pytest.raises(ValueError, match="non-negative"):
        GenomicInterval(source_assembly, "chr16", -1, 10)

    with pytest.raises(ValueError, match="greater than or equal"):
        GenomicInterval(source_assembly, "chr16", 11, 10)


def test_candidate_can_carry_multiple_distinct_observations_from_one_source(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    alignment = ProvenanceSource(
        source_id="example-alignment",
        label="example source-to-target alignment",
    )
    chain_file = ProvenanceSource(
        source_id="example-chain-file",
        label="example chain file",
        identifiers=(
            ProvenanceIdentifier(
                ProvenanceIdentifierKind.SHA256,
                "sha256:" + "a" * 64,
            ),
        ),
        derived_from=(alignment,),
    )
    chain_score = EvidenceObservation(
        observation_id="chain-score",
        kind=EvidenceKind.CHAIN_SCORE,
        value=18432,
        provenance=chain_file,
    )
    duplicated_bases = EvidenceObservation(
        observation_id="qdup",
        kind=EvidenceKind.DUPLICATED_QUERY_BASES,
        value=42,
        provenance=chain_file,
    )

    candidate = NormalizedCandidate(
        candidate_id="candidate-a",
        target_interval=GenomicInterval(target_assembly, "chr16", 100, 150),
        orientation=MappingOrientation.SAME,
        mapping_provenance=chain_file,
        segments=(_segment(source_assembly, target_assembly),),
        evidence=(chain_score, duplicated_bases),
    )

    assert len(candidate.evidence) == 2
    assert candidate.evidence[0] is not candidate.evidence[1]
    assert candidate.evidence[0].provenance is candidate.evidence[1].provenance


def test_assessment_references_candidate_supporting_and_contradicting_evidence(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    source = ProvenanceSource("alignment", "alignment source")
    support = EvidenceObservation(
        "support",
        EvidenceKind.MAPPING_COVERAGE,
        "full",
        source,
    )
    contradiction = EvidenceObservation(
        "contradiction",
        EvidenceKind.CHAIN_GAPS,
        True,
        source,
    )
    candidate = NormalizedCandidate(
        candidate_id="candidate-a",
        target_interval=GenomicInterval(target_assembly, "chr16", 100, 150),
        orientation=MappingOrientation.SAME,
        mapping_provenance=source,
        segments=(_segment(source_assembly, target_assembly),),
        evidence=(support, contradiction),
    )

    assessment = Assessment(
        source_interval=GenomicInterval(source_assembly, "chr16", 10, 60),
        verdict=Verdict.CONTESTED,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        candidates=(candidate,),
        preferred_candidate_id="candidate-a",
        supporting_evidence=(EvidenceReference("candidate-a", "support"),),
        contradicting_evidence=(
            EvidenceReference("candidate-a", "contradiction"),
        ),
    )

    assert assessment.supporting_evidence[0].candidate_id == "candidate-a"
    assert assessment.supporting_evidence[0].observation_id == "support"
    assert assessment.contradicting_evidence[0].observation_id == "contradiction"


def test_verdict_enum_contains_exactly_the_three_v1_labels() -> None:
    assert set(Verdict) == {
        Verdict.WELL_SUPPORTED,
        Verdict.CONTESTED,
        Verdict.INDETERMINATE,
    }


def test_well_supported_assessment_is_representable(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    source = ProvenanceSource("alignment", "alignment source")
    observation = EvidenceObservation(
        "coverage",
        EvidenceKind.MAPPING_COVERAGE,
        "full",
        source,
    )
    candidate = NormalizedCandidate(
        candidate_id="candidate-a",
        target_interval=GenomicInterval(target_assembly, "chr16", 100, 150),
        orientation=MappingOrientation.SAME,
        mapping_provenance=source,
        segments=(_segment(source_assembly, target_assembly),),
        evidence=(observation,),
    )

    assessment = Assessment(
        source_interval=GenomicInterval(source_assembly, "chr16", 10, 60),
        verdict=Verdict.WELL_SUPPORTED,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        candidates=(candidate,),
        preferred_candidate_id="candidate-a",
        supporting_evidence=(EvidenceReference("candidate-a", "coverage"),),
    )

    assert assessment.verdict is Verdict.WELL_SUPPORTED


def test_distinct_sources_can_be_related_through_shared_upstream_provenance(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    alignment = ProvenanceSource(
        source_id="example-alignment",
        label="example source-to-target alignment",
    )
    chain_source = ProvenanceSource(
        source_id="chain-resource",
        label="chain resource",
        derived_from=(alignment,),
    )
    net_source = ProvenanceSource(
        source_id="net-resource",
        label="net resource",
        derived_from=(alignment,),
    )
    chain_observation = EvidenceObservation(
        "chain-score",
        EvidenceKind.CHAIN_SCORE,
        20000,
        chain_source,
    )
    net_observation = EvidenceObservation(
        "net-class",
        EvidenceKind.NET_CLASSIFICATION,
        "syn",
        net_source,
    )
    candidate = NormalizedCandidate(
        candidate_id="candidate-a",
        target_interval=GenomicInterval(target_assembly, "chr16", 100, 200),
        orientation=MappingOrientation.SAME,
        mapping_provenance=chain_source,
        segments=(
            _segment(
                source_assembly,
                target_assembly,
                source_start=10,
                target_start=100,
                length=100,
            ),
        ),
        evidence=(chain_observation, net_observation),
    )

    assert chain_observation.provenance is not net_observation.provenance
    assert chain_observation.provenance.derived_from == (alignment,)
    assert net_observation.provenance.derived_from == (alignment,)
    assert len(candidate.evidence) == 2


def test_dependent_evidence_can_remain_indeterminate(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    alignment = ProvenanceSource("alignment", "shared upstream alignment")
    chain_source = ProvenanceSource(
        "chain",
        "chain derived from alignment",
        derived_from=(alignment,),
    )
    net_source = ProvenanceSource(
        "net",
        "net derived from alignment",
        derived_from=(alignment,),
    )
    candidate = NormalizedCandidate(
        candidate_id="candidate-a",
        target_interval=GenomicInterval(target_assembly, "chr16", 100, 150),
        orientation=MappingOrientation.SAME,
        mapping_provenance=chain_source,
        segments=(_segment(source_assembly, target_assembly),),
        evidence=(
            EvidenceObservation(
                "chain-score",
                EvidenceKind.CHAIN_SCORE,
                19000,
                chain_source,
            ),
            EvidenceObservation(
                "net-class",
                EvidenceKind.NET_CLASSIFICATION,
                "syn",
                net_source,
            ),
        ),
    )

    assessment = Assessment(
        source_interval=GenomicInterval(source_assembly, "chr16", 10, 60),
        verdict=Verdict.INDETERMINATE,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        candidates=(candidate,),
        supporting_evidence=(
            EvidenceReference("candidate-a", "chain-score"),
            EvidenceReference("candidate-a", "net-class"),
        ),
    )

    assert assessment.verdict is Verdict.INDETERMINATE
    assert len(assessment.supporting_evidence) == 2
    upstream_sources = {
        observation.provenance.derived_from[0].source_id
        for observation in candidate.evidence
    }
    assert upstream_sources == {"alignment"}


def test_mapping_segment_requires_equal_nonzero_lengths(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    with pytest.raises(ValueError, match="at least one base"):
        MappingSegment(
            GenomicInterval(source_assembly, "chr16", 10, 10),
            GenomicInterval(target_assembly, "chr16", 100, 100),
        )

    with pytest.raises(ValueError, match="lengths must match"):
        MappingSegment(
            GenomicInterval(source_assembly, "chr16", 10, 20),
            GenomicInterval(target_assembly, "chr16", 100, 111),
        )


def test_candidate_target_interval_is_exact_segment_bounding_span(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    provenance = ProvenanceSource("chain", "chain resource")
    segments = (
        _segment(
            source_assembly,
            target_assembly,
            source_start=10,
            target_start=100,
            length=5,
        ),
        _segment(
            source_assembly,
            target_assembly,
            source_start=20,
            target_start=120,
            length=5,
        ),
    )

    candidate = NormalizedCandidate(
        candidate_id="candidate-a",
        target_interval=GenomicInterval(target_assembly, "chr16", 100, 125),
        orientation=MappingOrientation.SAME,
        mapping_provenance=provenance,
        segments=segments,
    )

    assert candidate.target_interval == GenomicInterval(
        target_assembly, "chr16", 100, 125
    )
    assert candidate.segments == segments

    with pytest.raises(ValueError, match="exactly bound"):
        NormalizedCandidate(
            candidate_id="candidate-b",
            target_interval=GenomicInterval(target_assembly, "chr16", 100, 130),
            orientation=MappingOrientation.SAME,
            mapping_provenance=provenance,
            segments=segments,
        )


def test_candidate_segments_must_be_source_ordered_and_nonoverlapping(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    provenance = ProvenanceSource("chain", "chain resource")
    segments = (
        _segment(
            source_assembly,
            target_assembly,
            source_start=20,
            target_start=120,
            length=5,
        ),
        _segment(
            source_assembly,
            target_assembly,
            source_start=10,
            target_start=100,
            length=5,
        ),
    )

    with pytest.raises(ValueError, match="ordered and non-overlapping"):
        NormalizedCandidate(
            candidate_id="candidate-a",
            target_interval=GenomicInterval(target_assembly, "chr16", 100, 125),
            orientation=MappingOrientation.SAME,
            mapping_provenance=provenance,
            segments=segments,
        )


def test_candidate_rejects_target_segments_inconsistent_with_orientation(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    provenance = ProvenanceSource("chain", "chain source")

    with pytest.raises(ValueError, match="same-orientation"):
        NormalizedCandidate(
            candidate_id="bad-same",
            target_interval=GenomicInterval(target_assembly, "chr16", 100, 145),
            orientation=MappingOrientation.SAME,
            mapping_provenance=provenance,
            segments=(
                _segment(
                    source_assembly,
                    target_assembly,
                    source_start=10,
                    target_start=120,
                    length=10,
                ),
                _segment(
                    source_assembly,
                    target_assembly,
                    source_start=30,
                    target_start=100,
                    length=10,
                ),
            ),
        )

    with pytest.raises(ValueError, match="reverse-orientation"):
        NormalizedCandidate(
            candidate_id="bad-reverse",
            target_interval=GenomicInterval(target_assembly, "chr16", 100, 140),
            orientation=MappingOrientation.REVERSE,
            mapping_provenance=provenance,
            segments=(
                _segment(
                    source_assembly,
                    target_assembly,
                    source_start=10,
                    target_start=100,
                    length=10,
                ),
                _segment(
                    source_assembly,
                    target_assembly,
                    source_start=30,
                    target_start=130,
                    length=10,
                ),
            ),
        )


def test_reciprocal_best_membership_summary_enforces_status_counts(
    source_assembly: AssemblyIdentifier,
) -> None:
    covered = (GenomicInterval(source_assembly, "chr16", 10, 20),)

    full = ReciprocalBestMembershipSummary(
        status=ReciprocalBestMembershipStatus.FULL,
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        chains_examined=1,
        covered_source_bases=10,
        candidate_source_bases=10,
        covered_source_intervals=covered,
    )
    assert full.status is ReciprocalBestMembershipStatus.FULL

    with pytest.raises(ValueError, match="full reciprocal-best"):
        ReciprocalBestMembershipSummary(
            status=ReciprocalBestMembershipStatus.FULL,
            resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
            chains_examined=1,
            covered_source_bases=10,
            candidate_source_bases=20,
            covered_source_intervals=covered,
        )

    with pytest.raises(ValueError, match="no reciprocal-best"):
        ReciprocalBestMembershipSummary(
            status=ReciprocalBestMembershipStatus.NONE,
            resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
            chains_examined=1,
            covered_source_bases=10,
            candidate_source_bases=20,
            covered_source_intervals=covered,
        )

    with pytest.raises(ValueError, match="partial reciprocal-best"):
        ReciprocalBestMembershipSummary(
            status=ReciprocalBestMembershipStatus.PARTIAL,
            resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
            chains_examined=1,
            covered_source_bases=0,
            candidate_source_bases=20,
        )


def test_reciprocal_best_membership_intervals_must_match_covered_bases(
    source_assembly: AssemblyIdentifier,
) -> None:
    with pytest.raises(ValueError, match="account for every covered"):
        ReciprocalBestMembershipSummary(
            status=ReciprocalBestMembershipStatus.PARTIAL,
            resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
            chains_examined=1,
            covered_source_bases=5,
            candidate_source_bases=20,
            covered_source_intervals=(
                GenomicInterval(source_assembly, "chr16", 10, 20),
            ),
        )


def test_reciprocal_best_membership_rejects_negative_scan_count(
    source_assembly: AssemblyIdentifier,
) -> None:
    with pytest.raises(ValueError, match="chains_examined"):
        ReciprocalBestMembershipSummary(
            status=ReciprocalBestMembershipStatus.NONE,
            resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
            chains_examined=-1,
            covered_source_bases=0,
            candidate_source_bases=20,
        )
