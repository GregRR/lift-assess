"""Mechanical evidence extraction for chain-backed candidates.

Some chain evidence must be extracted while the raw ``ChainRecord`` is still
available. In particular, a requested locus may begin or end inside a chain gap;
that distinction cannot always be reconstructed from normalized aligned segments
alone. This module records those source-specific observations on the normalized
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
)


def _annotate_chain_mapping_structure(
    source_interval: GenomicInterval,
    chain: ChainRecord,
    candidate: NormalizedCandidate,
) -> NormalizedCandidate:
    """Attach source coverage and chain-gap observations to ``candidate``.

    Coverage answers only whether bases in the requested *source locus* are
    represented by exact mapping segments. Chain-gap evidence separately records
    source-side, destination-side, or double-sided block gaps through that locus.
    Both observations retain the candidate's mapping provenance and therefore do
    not become independent evidence merely because they are distinct facts.
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
    gaps = EvidenceObservation(
        observation_id=f"{candidate.candidate_id}:chain-gaps",
        kind=EvidenceKind.CHAIN_GAPS,
        value=ChainGapSummary(
            gaps=_chain_gaps_through_locus(source_interval, chain, candidate)
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


def _chain_gaps_through_locus(
    source_interval: GenomicInterval,
    chain: ChainRecord,
    candidate: NormalizedCandidate,
) -> tuple[ChainGap, ...]:
    gaps: list[ChainGap] = []
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

        if source_overlap is not None or query_only_gap_through_locus:
            target_gap = _target_gap_interval(
                chain,
                candidate,
                block_query_end,
                target_gap_bases,
            )
            gaps.append(
                ChainGap(
                    source_boundary=source_gap_start,
                    source_gap_overlap=source_overlap,
                    target_gap_interval=target_gap,
                )
            )

        source_cursor = source_gap_end
        query_cursor = block_query_end + target_gap_bases

    return tuple(gaps)


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
