import pytest

from liftassess import (
    AssemblyIdentifier,
    Assessment,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    EvidenceReference,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceIdentifier,
    ProvenanceIdentifierKind,
    ProvenanceSource,
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
        mapping_provenance=chain_file,
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
        "candidate-a",
        GenomicInterval(target_assembly, "chr16", 100, 150),
        source,
        (support, contradiction),
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
        "candidate-a",
        GenomicInterval(target_assembly, "chr16", 100, 150),
        source,
        (observation,),
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
        "candidate-a",
        GenomicInterval(target_assembly, "chr16", 100, 200),
        chain_source,
        (chain_observation, net_observation),
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
        mapping_provenance=chain_source,
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
