import pytest

from liftassess.chain import ChainBlock, ChainRecord, ChainStrand
from liftassess.models import (
    AssemblyIdentifier,
    EvidenceKind,
    GenomicInterval,
    MappingOrientation,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.projection import (
    iter_candidates_from_chains,
    project_interval_through_chain,
)


@pytest.fixture
def source_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier("sourceAsm", "test-provider")


@pytest.fixture
def target_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier("targetAsm", "test-provider")


@pytest.fixture
def chain_provenance() -> ProvenanceSource:
    return ProvenanceSource("chain-file", "test chain file")


def _chain(
    *,
    chain_id: int,
    blocks: tuple[ChainBlock, ...] = (ChainBlock(100),),
    target_name: str = "chr1",
    query_name: str = "chrA",
    query_strand: ChainStrand = ChainStrand.PLUS,
    target_start: int = 100,
    query_start: int = 500,
    score: float = 1000,
) -> ChainRecord:
    target_span = sum(block.size + (block.target_gap or 0) for block in blocks)
    query_span = sum(block.size + (block.query_gap or 0) for block in blocks)
    return ChainRecord(
        score=score,
        target_name=target_name,
        target_size=1000,
        target_strand=ChainStrand.PLUS,
        target_start=target_start,
        target_end=target_start + target_span,
        query_name=query_name,
        query_size=2000,
        query_strand=query_strand,
        query_start=query_start,
        query_end=query_start + query_span,
        chain_id=chain_id,
        blocks=blocks,
    )


def _project(
    source_interval: GenomicInterval,
    chain: ChainRecord,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> NormalizedCandidate | None:
    return project_interval_through_chain(
        source_interval,
        chain,
        target_assembly=target_assembly,
        mapping_provenance=chain_provenance,
    )


def test_projects_single_same_orientation_block(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    candidate = _project(
        GenomicInterval(source_assembly, "chr1", 120, 150),
        _chain(chain_id=7),
        target_assembly,
        chain_provenance,
    )

    assert candidate is not None
    assert candidate.candidate_id == "chain-file:chain:7"
    assert candidate.orientation is MappingOrientation.SAME
    assert candidate.target_interval == GenomicInterval(
        target_assembly, "chrA", 520, 550
    )
    assert len(candidate.segments) == 1
    assert candidate.segments[0].source_interval == GenomicInterval(
        source_assembly, "chr1", 120, 150
    )
    assert candidate.segments[0].target_interval == candidate.target_interval
    assert [observation.kind for observation in candidate.evidence] == [
        EvidenceKind.CHAIN_SCORE,
        EvidenceKind.MAPPING_COVERAGE,
        EvidenceKind.CHAIN_GAPS,
    ]
    assert candidate.evidence[0].value == 1000
    assert all(
        observation.provenance is chain_provenance for observation in candidate.evidence
    )


def test_distinct_chain_records_can_project_same_local_mapping(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    source = GenomicInterval(source_assembly, "chr1", 120, 150)
    first = _project(
        source,
        _chain(chain_id=31, target_start=100, query_start=500),
        target_assembly,
        chain_provenance,
    )
    second = _project(
        source,
        _chain(
            chain_id=32,
            blocks=(ChainBlock(120),),
            target_start=90,
            query_start=490,
        ),
        target_assembly,
        chain_provenance,
    )

    assert first is not None
    assert second is not None
    assert first.candidate_id != second.candidate_id
    assert first.orientation is second.orientation is MappingOrientation.SAME
    assert first.target_interval == second.target_interval
    assert first.segments == second.segments


def test_projects_reverse_query_coordinates_to_forward_reference(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    candidate = _project(
        GenomicInterval(source_assembly, "chr1", 120, 150),
        _chain(chain_id=8, query_strand=ChainStrand.MINUS),
        target_assembly,
        chain_provenance,
    )

    assert candidate is not None
    assert candidate.orientation is MappingOrientation.REVERSE
    assert candidate.target_interval == GenomicInterval(
        target_assembly, "chrA", 1450, 1480
    )


def test_split_mapping_remains_one_candidate_with_exact_segments(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        chain_id=9,
        blocks=(
            ChainBlock(size=10, target_gap=10, query_gap=15),
            ChainBlock(size=10),
        ),
    )
    candidate = _project(
        GenomicInterval(source_assembly, "chr1", 105, 125),
        chain,
        target_assembly,
        chain_provenance,
    )

    assert candidate is not None
    assert candidate.target_interval == GenomicInterval(
        target_assembly, "chrA", 505, 530
    )
    assert tuple(segment.source_interval for segment in candidate.segments) == (
        GenomicInterval(source_assembly, "chr1", 105, 110),
        GenomicInterval(source_assembly, "chr1", 120, 125),
    )
    assert tuple(segment.target_interval for segment in candidate.segments) == (
        GenomicInterval(target_assembly, "chrA", 505, 510),
        GenomicInterval(target_assembly, "chrA", 525, 530),
    )


def test_reverse_split_segments_stay_ordered_by_source_coordinates(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        chain_id=10,
        query_strand=ChainStrand.MINUS,
        blocks=(
            ChainBlock(size=10, target_gap=10, query_gap=15),
            ChainBlock(size=10),
        ),
    )
    candidate = _project(
        GenomicInterval(source_assembly, "chr1", 105, 125),
        chain,
        target_assembly,
        chain_provenance,
    )

    assert candidate is not None
    assert candidate.orientation is MappingOrientation.REVERSE
    assert tuple(segment.source_interval.start for segment in candidate.segments) == (
        105,
        120,
    )
    assert tuple(segment.target_interval for segment in candidate.segments) == (
        GenomicInterval(target_assembly, "chrA", 1490, 1495),
        GenomicInterval(target_assembly, "chrA", 1470, 1475),
    )
    assert candidate.target_interval == GenomicInterval(
        target_assembly, "chrA", 1470, 1495
    )


def test_returns_none_when_source_has_no_aligned_bases(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        chain_id=11,
        blocks=(
            ChainBlock(size=10, target_gap=10, query_gap=0),
            ChainBlock(size=10),
        ),
    )

    assert (
        _project(
            GenomicInterval(source_assembly, "chr1", 112, 118),
            chain,
            target_assembly,
            chain_provenance,
        )
        is None
    )


def test_returns_none_for_other_sequence_or_nonoverlapping_chain(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(chain_id=12)

    assert (
        _project(
            GenomicInterval(source_assembly, "chr2", 120, 150),
            chain,
            target_assembly,
            chain_provenance,
        )
        is None
    )
    assert (
        _project(
            GenomicInterval(source_assembly, "chr1", 10, 20),
            chain,
            target_assembly,
            chain_provenance,
        )
        is None
    )


def test_rejects_zero_length_source_interval_until_semantics_are_defined(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    with pytest.raises(ValueError, match="zero-length"):
        _project(
            GenomicInterval(source_assembly, "chr1", 120, 120),
            _chain(chain_id=13),
            target_assembly,
            chain_provenance,
        )


def test_rejects_source_interval_beyond_chain_target_sequence_size(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    with pytest.raises(ValueError, match="exceeds chain target sequence bounds"):
        _project(
            GenomicInterval(source_assembly, "chr1", 950, 1001),
            _chain(chain_id=14, target_start=900),
            target_assembly,
            chain_provenance,
        )


def test_iter_candidates_preserves_multiple_chain_mappings_without_ranking(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chains = (
        _chain(chain_id=21, query_name="chrA", score=200),
        _chain(chain_id=22, query_name="chrB", score=100),
        _chain(chain_id=23, target_name="chr2", query_name="chrC", score=999),
    )

    candidates = tuple(
        iter_candidates_from_chains(
            GenomicInterval(source_assembly, "chr1", 120, 150),
            chains,
            target_assembly=target_assembly,
            mapping_provenance=chain_provenance,
        )
    )

    assert [candidate.candidate_id for candidate in candidates] == [
        "chain-file:chain:21",
        "chain-file:chain:22",
    ]
    assert [candidate.target_interval.sequence_name for candidate in candidates] == [
        "chrA",
        "chrB",
    ]
