"""Mechanical evidence extraction for chain-backed candidates.

Some chain evidence must be extracted while the raw ``ChainRecord`` is still
available. In particular, a requested locus may begin or end inside a chain gap,
and a mapped point may lie directly beside a source- or target-side gap. Those
relationships cannot always be reconstructed from normalized aligned segments
alone. This module records the source-specific observations on the normalized
candidate without assigning a verdict or interpreting biological meaning.
"""

from dataclasses import replace

from .chain import ChainRecord, chain_candidate_id
from .models import (
    ChainGap,
    ChainGapSummary,
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    NormalizedCandidate,
    PointGapBoundary,
    PointGapBoundaryPosition,
)


def _annotate_chain_mapping_structure(
    source_interval: GenomicInterval,
    chain: ChainRecord,
    candidate: NormalizedCandidate,
) -> NormalizedCandidate:
    """Attach source coverage and chain-gap observations to ``candidate``.

    Coverage answers only whether bases in the requested *source locus* are
    represented by exact mapping segments. Chain-gap evidence separately records
    source-side, target-side, or double-sided block gaps through that locus
    and exact gap adjacency for mapped 1-bp queries. Both observations retain the
    candidate's mapping provenance and therefore do not become independent
    evidence merely because they are distinct facts.
    """

    _validate_inputs(source_interval, chain, candidate)

    uncovered = _uncovered_source_intervals(source_interval, candidate)
    covered_source_bases = sum(
        segment.source_interval.length for segment in candidate.segments
    )
    coverage_status = (
        MappingCoverageStatus.FULL
        if covered_source_bases == source_interval.length
        else MappingCoverageStatus.PARTIAL
    )
    coverage = EvidenceObservation(
        observation_id=f"{candidate.candidate_id}:mapping-coverage",
        kind=EvidenceKind.MAPPING_COVERAGE,
        value=MappingCoverageSummary(
            status=coverage_status,
            covered_source_bases=covered_source_bases,
            source_bases=source_interval.length,
            uncovered_source_intervals=uncovered,
        ),
        provenance=candidate.mapping_provenance,
    )
    chain_gaps, point_gap_boundaries = _chain_gap_context(
        source_interval,
        chain,
        candidate,
    )
    gaps = EvidenceObservation(
        observation_id=f"{candidate.candidate_id}:chain-gaps",
        kind=EvidenceKind.CHAIN_GAPS,
        value=ChainGapSummary(
            gaps=chain_gaps,
            point_gap_boundaries=(
                point_gap_boundaries if source_interval.length == 1 else None
            ),
        ),
        provenance=candidate.mapping_provenance,
    )

    return replace(candidate, evidence=(*candidate.evidence, coverage, gaps))


def _validate_inputs(
    source_interval: GenomicInterval,
    chain: ChainRecord,
    candidate: NormalizedCandidate,
) -> None:
    if source_interval.length <= 0:
        raise ValueError("mapping-structure evidence requires a non-empty source locus")
    if source_interval.sequence_name != chain.target_name:
        raise ValueError("source locus sequence does not match chain target sequence")
    if candidate.target_interval.sequence_name != chain.query_name:
        raise ValueError(
            "candidate target sequence does not match chain query sequence"
        )
    if candidate.orientation is not chain.orientation:
        raise ValueError("candidate orientation does not match chain orientation")

    expected_candidate_id = chain_candidate_id(
        candidate.mapping_provenance.source_id, chain.chain_id
    )
    if candidate.candidate_id != expected_candidate_id:
        raise ValueError("candidate identity does not match chain provenance and ID")

    first_source = candidate.segments[0].source_interval
    if (
        first_source.assembly != source_interval.assembly
        or first_source.sequence_name != source_interval.sequence_name
    ):
        raise ValueError("candidate source segments do not match the source locus")

    for segment in candidate.segments:
        if (
            segment.source_interval.start < source_interval.start
            or segment.source_interval.end > source_interval.end
        ):
            raise ValueError("candidate source segment lies outside the source locus")


def _uncovered_source_intervals(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
) -> tuple[GenomicInterval, ...]:
    uncovered: list[GenomicInterval] = []
    cursor = source_interval.start

    for segment in candidate.segments:
        if cursor < segment.source_interval.start:
            uncovered.append(
                GenomicInterval(
                    assembly=source_interval.assembly,
                    sequence_name=source_interval.sequence_name,
                    start=cursor,
                    end=segment.source_interval.start,
                )
            )
        cursor = segment.source_interval.end

    if cursor < source_interval.end:
        uncovered.append(
            GenomicInterval(
                assembly=source_interval.assembly,
                sequence_name=source_interval.sequence_name,
                start=cursor,
                end=source_interval.end,
            )
        )

    return tuple(uncovered)


def _chain_gap_context(
    source_interval: GenomicInterval,
    chain: ChainRecord,
    candidate: NormalizedCandidate,
) -> tuple[tuple[ChainGap, ...], tuple[PointGapBoundary, ...]]:
    gaps: list[ChainGap] = []
    point_boundaries: list[PointGapBoundary] = []
    source_cursor = chain.target_start
    query_cursor = chain.query_start

    for block in chain.blocks:
        block_source_end = source_cursor + block.size
        block_query_end = query_cursor + block.size
        if block.is_terminal:
            break

        source_gap_bases, target_gap_bases = block.gaps_after()
        source_gap_start = block_source_end
        source_gap_end = source_gap_start + source_gap_bases

        source_overlap = _gap_source_overlap(
            source_interval,
            source_gap_start,
            source_gap_end,
        )
        query_only_gap_through_locus = (
            source_gap_bases == 0
            and target_gap_bases > 0
            and source_interval.start < source_gap_start < source_interval.end
        )
        point_at_internal_boundary = source_interval.length == 1 and (
            source_interval.end == source_gap_start
            or source_interval.start == source_gap_end
        )
        needs_gap_intervals = (
            source_overlap is not None
            or query_only_gap_through_locus
            or point_at_internal_boundary
        )
        source_gap = (
            _source_gap_interval(
                source_interval,
                source_gap_start,
                source_gap_end,
            )
            if point_at_internal_boundary
            else None
        )
        target_gap = (
            _target_gap_interval(
                chain,
                candidate,
                block_query_end,
                target_gap_bases,
            )
            if needs_gap_intervals
            else None
        )

        if source_overlap is not None or query_only_gap_through_locus:
            gaps.append(
                ChainGap(
                    source_boundary=source_gap_start,
                    source_gap_overlap=source_overlap,
                    target_gap_interval=target_gap,
                )
            )

        if point_at_internal_boundary:
            point_boundary = _point_gap_boundary(
                source_interval,
                candidate,
                source_gap,
                target_gap,
            )
            if point_boundary is not None:
                point_boundaries.append(point_boundary)

        source_cursor = source_gap_end
        query_cursor = block_query_end + target_gap_bases

    return tuple(gaps), tuple(point_boundaries)


def _source_gap_interval(
    source_interval: GenomicInterval,
    gap_start: int,
    gap_end: int,
) -> GenomicInterval | None:
    if gap_start == gap_end:
        return None
    return GenomicInterval(
        assembly=source_interval.assembly,
        sequence_name=source_interval.sequence_name,
        start=gap_start,
        end=gap_end,
    )


def _point_gap_boundary(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
    source_gap: GenomicInterval | None,
    target_gap: GenomicInterval | None,
) -> PointGapBoundary | None:
    if source_interval.length != 1:
        raise ValueError("point gap-boundary context requires a 1-bp source interval")
    if len(candidate.segments) != 1:
        raise ValueError("a mapped 1-bp query must contain exactly one mapping segment")

    source_position = _point_position_for_gap(source_interval, source_gap)
    target_position = _point_position_for_gap(
        candidate.segments[0].target_interval,
        target_gap,
    )
    if source_position is None and target_position is None:
        return None
    return PointGapBoundary(
        source_gap_interval=(source_gap if source_position is not None else None),
        source_position=source_position,
        target_gap_interval=(target_gap if target_position is not None else None),
        target_position=target_position,
    )


def _point_position_for_gap(
    point: GenomicInterval,
    gap: GenomicInterval | None,
) -> PointGapBoundaryPosition | None:
    if gap is None:
        return None
    if point.end == gap.start:
        return PointGapBoundaryPosition.BEFORE_GAP
    if point.start == gap.end:
        return PointGapBoundaryPosition.AFTER_GAP
    return None


def _gap_source_overlap(
    source_interval: GenomicInterval,
    gap_start: int,
    gap_end: int,
) -> GenomicInterval | None:
    if gap_start == gap_end:
        return None
    overlap_start = max(source_interval.start, gap_start)
    overlap_end = min(source_interval.end, gap_end)
    if overlap_start >= overlap_end:
        return None
    return GenomicInterval(
        assembly=source_interval.assembly,
        sequence_name=source_interval.sequence_name,
        start=overlap_start,
        end=overlap_end,
    )


def _target_gap_interval(
    chain: ChainRecord,
    candidate: NormalizedCandidate,
    query_gap_start: int,
    query_gap_bases: int,
) -> GenomicInterval | None:
    if query_gap_bases == 0:
        return None

    query_gap_end = query_gap_start + query_gap_bases
    target_start, target_end = chain.query_interval_to_forward(
        query_gap_start, query_gap_end
    )

    return GenomicInterval(
        assembly=candidate.target_interval.assembly,
        sequence_name=candidate.target_interval.sequence_name,
        start=target_start,
        end=target_end,
    )
