import pytest

from liftassess.chain import ChainBlock, ChainRecord, ChainStrand
from liftassess.evidence import _annotate_chain_mapping_structure
from liftassess.models import (
    AssemblyIdentifier,
    ChainGapSummary,
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.projection import project_interval_through_chain


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
    chain_id: int = 1,
    blocks: tuple[ChainBlock, ...] = (ChainBlock(20),),
    target_start: int = 100,
    query_start: int = 500,
    query_strand: ChainStrand = ChainStrand.PLUS,
) -> ChainRecord:
    target_span = sum(block.size + (block.target_gap or 0) for block in blocks)
    query_span = sum(block.size + (block.query_gap or 0) for block in blocks)
    return ChainRecord(
        score=1000,
        target_name="chr1",
        target_size=1000,
        target_strand=ChainStrand.PLUS,
        target_start=target_start,
        target_end=target_start + target_span,
        query_name="chrA",
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
) -> NormalizedCandidate:
    candidate = project_interval_through_chain(
        source_interval,
        chain,
        target_assembly=target_assembly,
        mapping_provenance=chain_provenance,
    )
    assert candidate is not None
    return candidate


def _observation(
    candidate: NormalizedCandidate, kind: EvidenceKind
) -> EvidenceObservation:
    return next(item for item in candidate.evidence if item.kind is kind)


def test_full_contiguous_mapping_records_full_coverage_and_no_gaps(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    candidate = _project(source, _chain(), target_assembly, chain_provenance)

    coverage = _observation(candidate, EvidenceKind.MAPPING_COVERAGE)
    gaps = _observation(candidate, EvidenceKind.CHAIN_GAPS)

    assert coverage.value == MappingCoverageSummary(
        status=MappingCoverageStatus.FULL,
        covered_source_bases=10,
        source_bases=10,
    )
    assert gaps.value == ChainGapSummary()
    assert coverage.provenance is chain_provenance
    assert gaps.provenance is chain_provenance
    assert candidate.evidence[0].kind is EvidenceKind.CHAIN_SCORE


def test_source_side_gap_makes_coverage_partial_and_records_exact_overlap(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        blocks=(
            ChainBlock(size=10, target_gap=10, query_gap=0),
            ChainBlock(size=10),
        )
    )
    source = GenomicInterval(source_assembly, "chr1", 105, 125)
    candidate = _project(source, chain, target_assembly, chain_provenance)

    coverage = _observation(candidate, EvidenceKind.MAPPING_COVERAGE).value
    gap_summary = _observation(candidate, EvidenceKind.CHAIN_GAPS).value

    assert isinstance(coverage, MappingCoverageSummary)
    assert coverage.status is MappingCoverageStatus.PARTIAL
    assert coverage.covered_source_bases == 10
    assert coverage.source_bases == 20
    assert coverage.uncovered_source_intervals == (
        GenomicInterval(source_assembly, "chr1", 110, 120),
    )

    assert isinstance(gap_summary, ChainGapSummary)
    assert len(gap_summary.gaps) == 1
    gap = gap_summary.gaps[0]
    assert gap.source_boundary == 110
    assert gap.source_gap_overlap == GenomicInterval(
        source_assembly, "chr1", 110, 120
    )
    assert gap.target_gap_interval is None


def test_destination_only_gap_does_not_make_source_coverage_partial(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        blocks=(
            ChainBlock(size=10, target_gap=0, query_gap=15),
            ChainBlock(size=10),
        )
    )
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    candidate = _project(source, chain, target_assembly, chain_provenance)

    coverage = _observation(candidate, EvidenceKind.MAPPING_COVERAGE).value
    gap_summary = _observation(candidate, EvidenceKind.CHAIN_GAPS).value

    assert isinstance(coverage, MappingCoverageSummary)
    assert coverage.status is MappingCoverageStatus.FULL
    assert coverage.covered_source_bases == coverage.source_bases == 10

    assert isinstance(gap_summary, ChainGapSummary)
    assert len(gap_summary.gaps) == 1
    gap = gap_summary.gaps[0]
    assert gap.source_boundary == 110
    assert gap.source_gap_overlap is None
    assert gap.target_gap_interval == GenomicInterval(
        target_assembly, "chrA", 510, 525
    )


def test_double_sided_gap_preserves_both_sides_without_interpretation(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        blocks=(
            ChainBlock(size=10, target_gap=10, query_gap=15),
            ChainBlock(size=10),
        )
    )
    source = GenomicInterval(source_assembly, "chr1", 105, 125)
    candidate = _project(source, chain, target_assembly, chain_provenance)

    gap_summary = _observation(candidate, EvidenceKind.CHAIN_GAPS).value

    assert isinstance(gap_summary, ChainGapSummary)
    assert gap_summary.gaps[0].source_gap_overlap == GenomicInterval(
        source_assembly, "chr1", 110, 120
    )
    assert gap_summary.gaps[0].target_gap_interval == GenomicInterval(
        target_assembly, "chrA", 510, 525
    )


def test_locus_starting_inside_chain_gap_is_not_misclassified_as_chain_edge(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        blocks=(
            ChainBlock(size=10, target_gap=10, query_gap=15),
            ChainBlock(size=10),
        )
    )
    source = GenomicInterval(source_assembly, "chr1", 115, 125)
    candidate = _project(source, chain, target_assembly, chain_provenance)

    coverage = _observation(candidate, EvidenceKind.MAPPING_COVERAGE).value
    gap_summary = _observation(candidate, EvidenceKind.CHAIN_GAPS).value

    assert isinstance(coverage, MappingCoverageSummary)
    assert coverage.status is MappingCoverageStatus.PARTIAL
    assert coverage.uncovered_source_intervals == (
        GenomicInterval(source_assembly, "chr1", 115, 120),
    )
    assert isinstance(gap_summary, ChainGapSummary)
    assert gap_summary.gaps[0].source_gap_overlap == GenomicInterval(
        source_assembly, "chr1", 115, 120
    )


def test_uncovered_locus_edge_outside_chain_is_not_reported_as_chain_gap(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    source = GenomicInterval(source_assembly, "chr1", 95, 105)
    candidate = _project(source, _chain(), target_assembly, chain_provenance)

    coverage = _observation(candidate, EvidenceKind.MAPPING_COVERAGE).value
    gap_summary = _observation(candidate, EvidenceKind.CHAIN_GAPS).value

    assert isinstance(coverage, MappingCoverageSummary)
    assert coverage.status is MappingCoverageStatus.PARTIAL
    assert coverage.uncovered_source_intervals == (
        GenomicInterval(source_assembly, "chr1", 95, 100),
    )
    assert gap_summary == ChainGapSummary()


def test_reverse_chain_gap_is_normalized_to_forward_target_coordinates(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    chain = _chain(
        query_strand=ChainStrand.MINUS,
        blocks=(
            ChainBlock(size=10, target_gap=0, query_gap=15),
            ChainBlock(size=10),
        ),
    )
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    candidate = _project(source, chain, target_assembly, chain_provenance)

    gap_summary = _observation(candidate, EvidenceKind.CHAIN_GAPS).value

    assert isinstance(gap_summary, ChainGapSummary)
    assert gap_summary.gaps[0].target_gap_interval == GenomicInterval(
        target_assembly, "chrA", 1475, 1490
    )


def test_annotation_rejects_candidate_with_wrong_chain_orientation(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    chain = _chain(query_strand=ChainStrand.MINUS)
    candidate = NormalizedCandidate(
        candidate_id="manual",
        target_interval=GenomicInterval(target_assembly, "chrA", 505, 515),
        orientation=MappingOrientation.SAME,
        mapping_provenance=chain_provenance,
        segments=(
            MappingSegment(
                source_interval=source,
                target_interval=GenomicInterval(target_assembly, "chrA", 505, 515),
            ),
        ),
    )

    with pytest.raises(ValueError, match="orientation"):
        _annotate_chain_mapping_structure(source, chain, candidate)


def test_annotation_rejects_candidate_from_different_chain(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    source = GenomicInterval(source_assembly, "chr1", 105, 115)
    candidate = _project(
        source,
        _chain(chain_id=1),
        target_assembly,
        chain_provenance,
    )

    with pytest.raises(ValueError, match="identity"):
        _annotate_chain_mapping_structure(
            source,
            _chain(chain_id=2),
            candidate,
        )
