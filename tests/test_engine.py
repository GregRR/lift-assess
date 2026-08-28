from collections.abc import Collection

import pytest

import liftassess.engine as engine_module
from liftassess.chain import ChainBlock, ChainRecord, ChainStrand
from liftassess.engine import build_ucsc_candidates
from liftassess.models import (
    AssemblyIdentifier,
    EvidenceKind,
    EvidenceValue,
    GenomicInterval,
    MappingOrientation,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
)
from liftassess.net import NetClassification, NetRecord, NetRecordKind
from liftassess.reciprocal_best import _annotate_candidate_with_reciprocal_best_chains


@pytest.fixture
def source_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier("sourceAsm", "test-provider")


@pytest.fixture
def target_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier("targetAsm", "test-provider")


@pytest.fixture
def provenance() -> tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource]:
    alignment = ProvenanceSource("alignment", "shared alignment")
    chain = ProvenanceSource("chain-file", "chain file", derived_from=(alignment,))
    net = ProvenanceSource("net-file", "net file", derived_from=(alignment,))
    rbest = ProvenanceSource("rbest-file", "rbest chain", derived_from=(alignment,))
    return chain, net, rbest


def _chain(
    *,
    chain_id: int,
    query_name: str,
    target_start: int = 100,
    query_start: int = 500,
    target_name: str = "chr1",
    query_strand: ChainStrand = ChainStrand.PLUS,
    blocks: tuple[ChainBlock, ...] = (ChainBlock(20),),
) -> ChainRecord:
    target_span = sum(block.size + (block.target_gap or 0) for block in blocks)
    query_span = sum(block.size + (block.query_gap or 0) for block in blocks)
    return ChainRecord(
        score=1000 + chain_id,
        target_name=target_name,
        target_size=5000,
        target_strand=ChainStrand.PLUS,
        target_start=target_start,
        target_end=target_start + target_span,
        query_name=query_name,
        query_size=5000,
        query_strand=query_strand,
        query_start=query_start,
        query_end=query_start + query_span,
        chain_id=chain_id,
        blocks=blocks,
    )


def _net_fill(*, chain_id: int, query_name: str, depth: int = 1) -> NetRecord:
    return NetRecord(
        target_name="chr1",
        target_sequence_size=5000,
        depth=depth,
        kind=NetRecordKind.FILL,
        target_start=100,
        target_span_size=20,
        query_name=query_name,
        orientation=MappingOrientation.SAME,
        query_start=500,
        query_span_size=20,
        chain_id=chain_id,
        aligned_bases=20,
        duplicated_query_bases=0,
        classification=NetClassification.SYNTENIC,
    )


def _observation_value(
    candidate: NormalizedCandidate, kind: EvidenceKind
) -> EvidenceValue:
    return next(item.value for item in candidate.evidence if item.kind is kind)


def test_chain_only_engine_builds_all_candidates_without_comparative_evidence(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, _ = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)

    candidates = build_ucsc_candidates(
        source,
        (_chain(chain_id=1, query_name="chrA"), _chain(chain_id=2, query_name="chrB")),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
    )

    assert [candidate.target_interval.sequence_name for candidate in candidates] == [
        "chrA",
        "chrB",
    ]
    assert all(
        EvidenceKind.NET_CLASSIFICATION
        not in {item.kind for item in candidate.evidence}
        for candidate in candidates
    )
    assert all(
        EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
        not in {item.kind for item in candidate.evidence}
        for candidate in candidates
    )


def test_one_shot_net_stream_annotates_multiple_candidates(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, net_provenance, _ = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    net_stream = (
        record
        for record in (
            _net_fill(chain_id=1, query_name="chrA"),
            _net_fill(chain_id=2, query_name="chrB"),
        )
    )

    candidates = build_ucsc_candidates(
        source,
        (_chain(chain_id=1, query_name="chrA"), _chain(chain_id=2, query_name="chrB")),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_records=net_stream,
        net_provenance=net_provenance,
    )

    assert len(candidates) == 2
    assert all(
        _observation_value(candidate, EvidenceKind.NET_CLASSIFICATION) == "syn"
        for candidate in candidates
    )


def test_one_shot_reciprocal_best_stream_annotates_multiple_candidates(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, rbest_provenance = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    reciprocal_stream = (
        chain
        for chain in (
            _chain(chain_id=101, query_name="chrA"),
            _chain(chain_id=202, query_name="chrB"),
            _chain(chain_id=303, query_name="chrUnused"),
        )
    )

    candidates = build_ucsc_candidates(
        source,
        (_chain(chain_id=1, query_name="chrA"), _chain(chain_id=2, query_name="chrB")),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        reciprocal_best_chains=reciprocal_stream,
        reciprocal_best_provenance=rbest_provenance,
        reciprocal_best_completeness=(
            ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
        ),
    )

    assert len(candidates) == 2
    for candidate in candidates:
        summary = _observation_value(candidate, EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP)
        assert isinstance(summary, ReciprocalBestMembershipSummary)
        assert summary.status is ReciprocalBestMembershipStatus.FULL
        assert (
            summary.resource_completeness
            is ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
        )
        assert summary.chains_examined == 1


def test_reciprocal_best_candidates_receive_only_matching_pair_chains(
    monkeypatch: pytest.MonkeyPatch,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, rbest_provenance = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    collection_sizes: list[int] = []
    original = _annotate_candidate_with_reciprocal_best_chains

    def counting_annotator(
        candidate: NormalizedCandidate,
        *,
        reciprocal_best_chains: Collection[ChainRecord],
        resource_completeness: ReciprocalBestResourceCompleteness,
        reciprocal_best_provenance: ProvenanceSource,
    ) -> NormalizedCandidate:
        assert isinstance(reciprocal_best_chains, list)
        collection_sizes.append(len(reciprocal_best_chains))
        return original(
            candidate,
            reciprocal_best_chains=reciprocal_best_chains,
            resource_completeness=resource_completeness,
            reciprocal_best_provenance=reciprocal_best_provenance,
        )

    monkeypatch.setattr(
        engine_module,
        "_annotate_candidate_with_reciprocal_best_chains",
        counting_annotator,
    )

    candidates = build_ucsc_candidates(
        source,
        (
            _chain(chain_id=1, query_name="chrA"),
            _chain(chain_id=2, query_name="chrB"),
        ),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        reciprocal_best_chains=(
            chain
            for chain in (
                _chain(chain_id=101, query_name="chrA"),
                _chain(chain_id=202, query_name="chrB"),
            )
        ),
        reciprocal_best_provenance=rbest_provenance,
        reciprocal_best_completeness=(
            ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
        ),
    )

    assert len(candidates) == 2
    assert collection_sizes == [1, 1]


def test_reverse_split_candidate_keeps_comparative_evidence_at_engine_boundary(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, net_provenance, rbest_provenance = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 125)
    split_blocks = (
        ChainBlock(size=10, target_gap=10, query_gap=15),
        ChainBlock(size=10),
    )
    candidate_chain = _chain(
        chain_id=1,
        query_name="chrA",
        query_strand=ChainStrand.MINUS,
        blocks=split_blocks,
    )
    net_fill = NetRecord(
        target_name="chr1",
        target_sequence_size=5000,
        depth=1,
        kind=NetRecordKind.FILL,
        target_start=100,
        target_span_size=30,
        query_name="chrA",
        orientation=MappingOrientation.REVERSE,
        query_start=500,
        query_span_size=35,
        chain_id=1,
        aligned_bases=20,
        duplicated_query_bases=0,
        classification=NetClassification.SYNTENIC,
    )
    reciprocal_chain = _chain(
        chain_id=101,
        query_name="chrA",
        query_strand=ChainStrand.MINUS,
        blocks=split_blocks,
    )

    (candidate,) = build_ucsc_candidates(
        source,
        (candidate_chain,),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_records=(record for record in (net_fill,)),
        net_provenance=net_provenance,
        reciprocal_best_chains=(chain for chain in (reciprocal_chain,)),
        reciprocal_best_provenance=rbest_provenance,
        reciprocal_best_completeness=(
            ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
        ),
    )

    assert candidate.orientation is MappingOrientation.REVERSE
    assert tuple(segment.source_interval.start for segment in candidate.segments) == (
        105,
        120,
    )
    assert tuple(
        (segment.target_interval.start, segment.target_interval.end)
        for segment in candidate.segments
    ) == ((4490, 4495), (4470, 4475))
    assert _observation_value(candidate, EvidenceKind.NET_CLASSIFICATION) == "syn"

    summary = _observation_value(candidate, EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP)
    assert isinstance(summary, ReciprocalBestMembershipSummary)
    assert summary.status is ReciprocalBestMembershipStatus.FULL
    assert (
        summary.resource_completeness
        is ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
    )
    assert summary.chains_examined == 1


def test_single_pass_rbest_preserves_internal_gap_as_partial_membership(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, rbest_provenance = provenance
    source = GenomicInterval(source_assembly, "chr1", 100, 120)
    candidate_chain = _chain(chain_id=1, query_name="chrA")
    rbest_chain = _chain(
        chain_id=101,
        query_name="chrA",
        blocks=(
            ChainBlock(size=10, target_gap=3, query_gap=3),
            ChainBlock(size=7),
        ),
    )

    (candidate,) = build_ucsc_candidates(
        source,
        (candidate_chain,),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        reciprocal_best_chains=(chain for chain in (rbest_chain,)),
        reciprocal_best_provenance=rbest_provenance,
        reciprocal_best_completeness=(
            ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
        ),
    )

    summary = _observation_value(candidate, EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP)
    assert isinstance(summary, ReciprocalBestMembershipSummary)
    assert summary.status is ReciprocalBestMembershipStatus.PARTIAL
    assert summary.covered_source_bases == 17
    assert summary.candidate_source_bases == 20
    assert [(item.start, item.end) for item in summary.covered_source_intervals] == [
        (100, 110),
        (113, 120),
    ]


def test_net_stream_retains_repeated_chain_id_fills_for_candidate(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, net_provenance, _ = provenance
    source = GenomicInterval(source_assembly, "chr1", 100, 120)
    records = (
        _net_fill(chain_id=1, query_name="chrA", depth=1),
        _net_fill(chain_id=1, query_name="chrA", depth=3),
    )

    (candidate,) = build_ucsc_candidates(
        source,
        (_chain(chain_id=1, query_name="chrA"),),
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_records=(record for record in records),
        net_provenance=net_provenance,
    )

    net_classes = [
        item
        for item in candidate.evidence
        if item.kind is EvidenceKind.NET_CLASSIFICATION
    ]
    hierarchy = [
        item for item in candidate.evidence if item.kind is EvidenceKind.NET_HIERARCHY
    ]
    assert len(net_classes) == 2
    assert len(hierarchy) == 2


def test_net_records_without_provenance_are_rejected(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, _ = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)

    with pytest.raises(ValueError, match="net records and net provenance"):
        build_ucsc_candidates(
            source,
            (_chain(chain_id=1, query_name="chrA"),),
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
            net_records=(),
        )


def test_net_provenance_without_records_is_rejected(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, net_provenance, _ = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)

    with pytest.raises(ValueError, match="net records and net provenance"):
        build_ucsc_candidates(
            source,
            (_chain(chain_id=1, query_name="chrA"),),
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
            net_provenance=net_provenance,
        )


def test_reciprocal_best_inputs_must_be_supplied_together(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, rbest_provenance = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)

    with pytest.raises(ValueError, match="reciprocal-best chains, provenance"):
        build_ucsc_candidates(
            source,
            (_chain(chain_id=1, query_name="chrA"),),
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
            reciprocal_best_chains=(),
            reciprocal_best_provenance=rbest_provenance,
        )


def test_duplicate_chain_candidate_identity_is_rejected(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, _ = provenance
    source = GenomicInterval(source_assembly, "chr1", 105, 115)

    with pytest.raises(ValueError, match="duplicate candidate identity"):
        build_ucsc_candidates(
            source,
            (
                _chain(chain_id=1, query_name="chrA"),
                _chain(chain_id=1, query_name="chrB"),
            ),
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
        )
