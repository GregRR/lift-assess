"""Mechanical projection of source intervals through parsed UCSC chains.

This module performs candidate generation only. It does not annotate candidates
with nets, rank candidates, interpret evidence, or assign assessment verdicts.
For UCSC liftOver map chains, the old/source assembly is the chain target side
(``t*`` fields) and the new/destination assembly is the query side (``q*``
fields).
"""

from collections.abc import Iterable, Iterator

from .chain import ChainRecord, chain_candidate_id
from .evidence import _annotate_chain_mapping_structure
from .models import (
    AssemblyIdentifier,
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
)


def project_interval_through_chain(
    source_interval: GenomicInterval,
    chain: ChainRecord,
    *,
    target_assembly: AssemblyIdentifier,
    mapping_provenance: ProvenanceSource,
) -> NormalizedCandidate | None:
    """Project one non-empty source interval through one parsed chain.

    The source interval is interpreted on the chain target/reference side,
    matching UCSC liftOver map-chain semantics. One chain produces at most one
    candidate. If the interval intersects several alignment blocks in that
    chain, those exact mappings remain separate ``MappingSegment`` objects on
    the same candidate.

    Return ``None`` when the source interval has no aligned bases in the chain.
    """

    if source_interval.length == 0:
        raise ValueError(
            "zero-length source interval projection is not defined for liftAssess v1"
        )
    if source_interval.sequence_name != chain.target_name:
        return None
    if source_interval.end > chain.target_size:
        raise ValueError("source interval exceeds chain target sequence bounds")
    if (
        source_interval.end <= chain.target_start
        or source_interval.start >= chain.target_end
    ):
        return None

    segments: list[MappingSegment] = []

    for block in chain.iter_aligned_blocks():
        overlap_start = max(source_interval.start, block.target_start)
        overlap_end = min(source_interval.end, block.target_end)
        if overlap_start >= overlap_end:
            continue

        target_start, target_end = block.query_interval_for_target_interval(
            overlap_start, overlap_end
        )
        segments.append(
            MappingSegment(
                source_interval=GenomicInterval(
                    assembly=source_interval.assembly,
                    sequence_name=source_interval.sequence_name,
                    start=overlap_start,
                    end=overlap_end,
                ),
                target_interval=GenomicInterval(
                    assembly=target_assembly,
                    sequence_name=chain.query_name,
                    start=target_start,
                    end=target_end,
                ),
            )
        )

    if not segments:
        return None

    candidate_id = chain_candidate_id(mapping_provenance.source_id, chain.chain_id)
    target_interval = GenomicInterval(
        assembly=target_assembly,
        sequence_name=chain.query_name,
        start=min(segment.target_interval.start for segment in segments),
        end=max(segment.target_interval.end for segment in segments),
    )
    chain_score = EvidenceObservation(
        observation_id=f"{candidate_id}:chain-score",
        kind=EvidenceKind.CHAIN_SCORE,
        value=chain.score,
        provenance=mapping_provenance,
    )

    candidate = NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=target_interval,
        orientation=chain.orientation,
        mapping_provenance=mapping_provenance,
        segments=tuple(segments),
        evidence=(chain_score,),
    )
    return _annotate_chain_mapping_structure(source_interval, chain, candidate)


def iter_candidates_from_chains(
    source_interval: GenomicInterval,
    chains: Iterable[ChainRecord],
    *,
    target_assembly: AssemblyIdentifier,
    mapping_provenance: ProvenanceSource,
) -> Iterator[NormalizedCandidate]:
    """Yield every chain-backed candidate with aligned source bases.

    Candidates are not ranked or filtered by score here. The function preserves
    the chain iterable's order and leaves comparative interpretation to later
    assessor logic.
    """

    for chain in chains:
        candidate = project_interval_through_chain(
            source_interval,
            chain,
            target_assembly=target_assembly,
            mapping_provenance=mapping_provenance,
        )
        if candidate is not None:
            yield candidate
