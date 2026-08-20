"""Mechanical UCSC net evidence annotation for chain-backed candidates.

This module lives on the candidate-generation side of the architecture boundary.
It matches parsed net fills to an already-normalized chain candidate, preserves
all relevant fill contexts, and attaches raw typed observations without ranking,
aggregating, or assigning an aggregate result verdict.
"""

from collections.abc import Iterable
from dataclasses import replace

from .chain import chain_candidate_id
from .models import (
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceSource,
)
from .net import NetRecord, NetRecordKind
from .provenance import shares_upstream_source


def _annotate_candidate_with_net_records(
    candidate: NormalizedCandidate,
    *,
    chain_id: int,
    net_records: Iterable[NetRecord],
    net_provenance: ProvenanceSource,
) -> NormalizedCandidate:
    """Attach raw net-fill observations relevant to ``candidate``.

    Matching is deliberately mechanical. A record must be a fill for the same
    chain, source/target sequence pair, and orientation, and its target-side span
    must overlap at least one exact aligned source segment on the candidate.
    Every matching fill is preserved separately; repeated chain IDs are never
    collapsed to a single net record.

    The net source must share provenance with the candidate's mapping source so
    net-derived observations cannot accidentally appear independent of the chain
    alignment from which the candidate was generated.
    """

    if chain_id < 0:
        raise ValueError("chain ID must be non-negative")

    expected_candidate_id = chain_candidate_id(
        candidate.mapping_provenance.source_id, chain_id
    )
    if candidate.candidate_id != expected_candidate_id:
        raise ValueError(
            "candidate identity does not match mapping provenance and chain ID"
        )

    if not shares_upstream_source(candidate.mapping_provenance, net_provenance):
        raise ValueError(
            "net provenance must share an upstream source with candidate mapping provenance"
        )

    observations = list(candidate.evidence)
    # File-encounter counter used only to make observation IDs unique. It is not rank.
    matched_fill_count = 0

    for record in net_records:
        if not _fill_matches_candidate(record, candidate, chain_id):
            continue

        matched_fill_count += 1
        fill_provenance = _fill_provenance(net_provenance, record)
        fill_key = f"{candidate.candidate_id}:net-fill:{matched_fill_count}"

        if record.aligned_bases is not None:
            observations.append(
                EvidenceObservation(
                    observation_id=f"{fill_key}:ali",
                    kind=EvidenceKind.ALIGNED_BASES,
                    value=record.aligned_bases,
                    provenance=fill_provenance,
                )
            )
        if record.duplicated_query_bases is not None:
            observations.append(
                EvidenceObservation(
                    observation_id=f"{fill_key}:qdup",
                    kind=EvidenceKind.DUPLICATED_QUERY_BASES,
                    value=record.duplicated_query_bases,
                    provenance=fill_provenance,
                )
            )
        if record.classification is not None:
            observations.append(
                EvidenceObservation(
                    observation_id=f"{fill_key}:classification",
                    kind=EvidenceKind.NET_CLASSIFICATION,
                    value=record.classification.value,
                    provenance=fill_provenance,
                )
            )

        source_segment = candidate.segments[0].source_interval
        observations.append(
            EvidenceObservation(
                observation_id=f"{fill_key}:hierarchy",
                kind=EvidenceKind.NET_HIERARCHY,
                value=NetHierarchySummary(
                    depth=record.depth,
                    source_fill_interval=GenomicInterval(
                        assembly=source_segment.assembly,
                        sequence_name=record.target_name,
                        start=record.target_start,
                        end=record.target_end,
                    ),
                ),
                provenance=fill_provenance,
            )
        )

    if matched_fill_count == 0:
        return candidate
    return replace(candidate, evidence=tuple(observations))


def _fill_matches_candidate(
    record: NetRecord,
    candidate: NormalizedCandidate,
    chain_id: int,
) -> bool:
    if record.kind is not NetRecordKind.FILL:
        return False
    if record.chain_id != chain_id:
        return False

    source_segment = candidate.segments[0].source_interval
    if record.target_name != source_segment.sequence_name:
        return False
    if record.query_name != candidate.target_interval.sequence_name:
        return False
    if record.orientation is not candidate.orientation:
        return False

    return any(
        _overlaps(
            segment.source_interval.start,
            segment.source_interval.end,
            record.target_start,
            record.target_end,
        )
        for segment in candidate.segments
    )


def _overlaps(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


def _fill_provenance(
    net_provenance: ProvenanceSource,
    record: NetRecord,
) -> ProvenanceSource:
    chain_text = "none" if record.chain_id is None else str(record.chain_id)
    record_identity = (
        f"chain-{chain_text}:t-{record.target_name}-{record.target_start}-"
        f"{record.target_end}:q-{record.query_name}-{record.orientation.value}-"
        f"{record.query_start}-{record.query_end}:depth-{record.depth}"
    )
    return ProvenanceSource(
        source_id=f"{net_provenance.source_id}:fill:{record_identity}",
        label=(
            f"net fill chain {chain_text} {record.target_name}:"
            f"{record.target_start}-{record.target_end} depth {record.depth}"
        ),
        derived_from=(net_provenance,),
    )
