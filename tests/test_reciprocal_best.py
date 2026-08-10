from collections.abc import Collection
from typing import cast

import pytest

from liftassess.chain import ChainBlock, ChainRecord, ChainStrand
from liftassess.models import (
    AssemblyIdentifier,
    EvidenceKind,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
)
from liftassess.projection import project_interval_through_chain
from liftassess.reciprocal_best import _annotate_candidate_with_reciprocal_best_chains


def _assemblies() -> tuple[AssemblyIdentifier, AssemblyIdentifier]:
    return (
        AssemblyIdentifier(name="source", provider="test"),
        AssemblyIdentifier(name="target", provider="test"),
    )


def _provenance() -> tuple[ProvenanceSource, ProvenanceSource]:
    alignment = ProvenanceSource(source_id="alignment", label="shared alignment")
    return (
        ProvenanceSource(
            source_id="chain-file",
            label="ordinary chain resource",
            derived_from=(alignment,),
        ),
        ProvenanceSource(
            source_id="rbest-chain-file",
            label="reciprocal-best chain resource",
            derived_from=(alignment,),
        ),
    )


def _chain(
    *,
    chain_id: int = 7,
    source_start: int = 100,
    query_start: int = 500,
    size: int = 20,
    query_strand: ChainStrand = ChainStrand.PLUS,
    query_size: int = 2000,
    source_name: str = "chr1",
    query_name: str = "chrA",
) -> ChainRecord:
    return ChainRecord(
        score=100,
        target_name=source_name,
        target_size=1000,
        target_strand=ChainStrand.PLUS,
        target_start=source_start,
        target_end=source_start + size,
        query_name=query_name,
        query_size=query_size,
        query_strand=query_strand,
        query_start=query_start,
        query_end=query_start + size,
        chain_id=chain_id,
        blocks=(ChainBlock(size=size),),
    )


def _candidate(
    chain: ChainRecord | None = None,
) -> tuple[
    NormalizedCandidate,
    ProvenanceSource,
    ProvenanceSource,
]:
    source_assembly, target_assembly = _assemblies()
    chain_provenance, rbest_provenance = _provenance()
    source_chain = chain or _chain()
    source_interval = GenomicInterval(
        assembly=source_assembly,
        sequence_name=source_chain.target_name,
        start=source_chain.target_start,
        end=source_chain.target_end,
    )
    candidate = project_interval_through_chain(
        source_interval,
        source_chain,
        target_assembly=target_assembly,
        mapping_provenance=chain_provenance,
    )
    assert candidate is not None
    return candidate, chain_provenance, rbest_provenance


def _membership(candidate: NormalizedCandidate) -> ReciprocalBestMembershipSummary:
    observations = [
        observation
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
    ]
    assert len(observations) == 1
    value = observations[0].value
    assert isinstance(value, ReciprocalBestMembershipSummary)
    return value


def test_full_membership_matches_exact_geometry_without_relying_on_chain_id() -> None:
    candidate, _, rbest_provenance = _candidate()
    reciprocal_chain = _chain(chain_id=999)

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(reciprocal_chain,),
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.FULL
    assert membership.covered_source_bases == 20
    assert membership.candidate_source_bases == 20
    assert [(interval.start, interval.end) for interval in membership.covered_source_intervals] == [
        (100, 120)
    ]
    assert annotated.evidence[-1].provenance is rbest_provenance


def test_partial_membership_preserves_exact_covered_source_interval() -> None:
    candidate, _, rbest_provenance = _candidate()
    reciprocal_chain = _chain(
        chain_id=88,
        source_start=105,
        query_start=505,
        size=10,
    )

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(reciprocal_chain,),
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.PARTIAL
    assert membership.covered_source_bases == 10
    assert [(interval.start, interval.end) for interval in membership.covered_source_intervals] == [
        (105, 115)
    ]


def test_none_membership_is_explicit_when_complete_resource_has_no_matching_geometry() -> None:
    candidate, _, rbest_provenance = _candidate()
    same_source_wrong_target = _chain(chain_id=88, query_start=700)

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(same_source_wrong_target,),
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.NONE
    assert (
        membership.resource_completeness
        is ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
    )
    assert membership.chains_examined == 1
    assert membership.covered_source_bases == 0
    assert membership.covered_source_intervals == ()


def test_two_reciprocal_best_fragments_can_combine_to_full_membership() -> None:
    candidate, _, rbest_provenance = _candidate()
    reciprocal_chains = (
        _chain(chain_id=80, source_start=100, query_start=500, size=10),
        _chain(chain_id=81, source_start=110, query_start=510, size=10),
    )

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=reciprocal_chains,
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.FULL
    assert membership.covered_source_bases == 20
    assert [(interval.start, interval.end) for interval in membership.covered_source_intervals] == [
        (100, 120)
    ]


def test_reverse_orientation_membership_uses_forward_reference_geometry() -> None:
    source_chain = _chain(
        query_strand=ChainStrand.MINUS,
        query_size=1000,
        query_start=100,
    )
    candidate, _, rbest_provenance = _candidate(source_chain)
    reciprocal_chain = _chain(
        chain_id=99,
        query_strand=ChainStrand.MINUS,
        query_size=1000,
        query_start=100,
    )

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(reciprocal_chain,),
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.FULL


def test_split_candidate_can_have_partial_membership_on_one_aligned_segment() -> None:
    source_chain = ChainRecord(
        score=100,
        target_name="chr1",
        target_size=1000,
        target_strand=ChainStrand.PLUS,
        target_start=100,
        target_end=130,
        query_name="chrA",
        query_size=2000,
        query_strand=ChainStrand.PLUS,
        query_start=500,
        query_end=520,
        chain_id=7,
        blocks=(
            ChainBlock(size=10, target_gap=10, query_gap=0),
            ChainBlock(size=10),
        ),
    )
    candidate, _, rbest_provenance = _candidate(source_chain)
    reciprocal_chain = _chain(chain_id=90, source_start=100, query_start=500, size=10)

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(reciprocal_chain,),
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.PARTIAL
    assert membership.covered_source_bases == 10
    assert membership.candidate_source_bases == 20


def test_reciprocal_best_internal_gap_is_not_bridged_in_membership_coverage() -> None:
    candidate, _, rbest_provenance = _candidate()
    reciprocal_chain = ChainRecord(
        score=100,
        target_name="chr1",
        target_size=1000,
        target_strand=ChainStrand.PLUS,
        target_start=100,
        target_end=120,
        query_name="chrA",
        query_size=2000,
        query_strand=ChainStrand.PLUS,
        query_start=500,
        query_end=520,
        chain_id=91,
        blocks=(
            ChainBlock(size=10, target_gap=3, query_gap=3),
            ChainBlock(size=7),
        ),
    )

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(reciprocal_chain,),
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.PARTIAL
    assert membership.covered_source_bases == 17
    assert [
        (interval.start, interval.end)
        for interval in membership.covered_source_intervals
    ] == [(100, 110), (113, 120)]


def test_wrong_sequence_or_orientation_does_not_count_as_membership() -> None:
    candidate, _, rbest_provenance = _candidate()
    reciprocal_chains = (
        _chain(chain_id=80, source_name="chr2"),
        _chain(chain_id=81, query_name="chrB"),
        _chain(
            chain_id=82,
            query_strand=ChainStrand.MINUS,
            query_size=2000,
            query_start=1480,
        ),
    )

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=reciprocal_chains,
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    assert _membership(annotated).status is ReciprocalBestMembershipStatus.NONE


def test_chains_examined_counts_only_candidate_sequence_pair_and_orientation() -> None:
    candidate, _, rbest_provenance = _candidate()
    reciprocal_chains = (
        _chain(chain_id=80, source_name="chr2"),
        _chain(chain_id=81, query_name="chrB"),
        _chain(
            chain_id=82,
            query_strand=ChainStrand.MINUS,
            query_size=2000,
            query_start=1480,
        ),
        _chain(chain_id=83),
    )

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=reciprocal_chains,
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.FULL
    assert membership.chains_examined == 1


def test_unrelated_reciprocal_best_provenance_is_rejected() -> None:
    candidate, _, _ = _candidate()
    unrelated = ProvenanceSource(source_id="other", label="unrelated alignment")

    with pytest.raises(ValueError, match="share an upstream source"):
        _annotate_candidate_with_reciprocal_best_chains(
            candidate,
            reciprocal_best_chains=(_chain(),),
            resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
            reciprocal_best_provenance=unrelated,
        )


def test_reverse_partial_membership_matches_subinterval_geometry() -> None:
    source_chain = _chain(
        query_strand=ChainStrand.MINUS,
        query_size=1000,
        query_start=100,
    )
    candidate, _, rbest_provenance = _candidate(source_chain)
    reciprocal_chain = _chain(
        chain_id=123,
        source_start=105,
        query_strand=ChainStrand.MINUS,
        query_size=1000,
        query_start=105,
        size=10,
    )

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(reciprocal_chain,),
        resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.PARTIAL
    assert membership.covered_source_bases == 10
    assert [(interval.start, interval.end) for interval in membership.covered_source_intervals] == [
        (105, 115)
    ]


def test_consumable_iterator_is_rejected_before_membership_is_computed() -> None:
    candidate, _, rbest_provenance = _candidate()
    consumable = iter((_chain(),))

    # cast simulates an untyped or incorrectly typed caller; runtime validation still
    # prevents a one-shot iterator from becoming a false exhaustive scan.
    with pytest.raises(TypeError, match="reusable materialized collection"):
        _annotate_candidate_with_reciprocal_best_chains(
            candidate,
            reciprocal_best_chains=cast(Collection[ChainRecord], consumable),
            resource_completeness=ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
            reciprocal_best_provenance=rbest_provenance,
        )


def test_complete_candidate_subset_basis_is_preserved_in_evidence() -> None:
    candidate, _, rbest_provenance = _candidate()

    annotated = _annotate_candidate_with_reciprocal_best_chains(
        candidate,
        reciprocal_best_chains=(),
        resource_completeness=(
            ReciprocalBestResourceCompleteness.COMPLETE_CANDIDATE_SUBSET
        ),
        reciprocal_best_provenance=rbest_provenance,
    )

    membership = _membership(annotated)
    assert membership.status is ReciprocalBestMembershipStatus.NONE
    assert (
        membership.resource_completeness
        is ReciprocalBestResourceCompleteness.COMPLETE_CANDIDATE_SUBSET
    )
    assert membership.chains_examined == 0
