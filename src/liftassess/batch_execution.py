"""Index-only batch execution over one exact cached UCSC chain resource.

This layer deliberately executes only chain-backed candidate generation.  It exists so
batch orchestration can reuse the prepared region-addressable chain index without
silently falling back to a whole-resource traversal.  Comparative net/reciprocal-best
batch evidence remains a later M22 slice rather than being approximated here.
"""

from dataclasses import dataclass

from .batch import (
    BatchInputRecord,
    BatchRecordAssessment,
    BatchRelationshipResult,
    build_batch_target_relationships,
)
from .chain_index import ChainIndex
from .models import AssemblyIdentifier, EvidenceAvailabilityTier, ProvenanceSource
from .resource_cache import CachedResource, CachedUCSCChainResource
from .resource_files import build_ucsc_chain_candidates_for_intervals_from_cached_chain
from .resource_identity import sha256_hex_from_identifier


@dataclass(frozen=True)
class IndexedChainBatchResult:
    """Chain-only indexed execution result for one source/target resource pair."""

    source_db: str
    target_db: str
    evidence_tier: EvidenceAvailabilityTier
    chain_sha256_identifier: str
    chain_resource: CachedResource
    alignment_provenance: ProvenanceSource
    record_assessments: tuple[BatchRecordAssessment, ...]
    relationships: BatchRelationshipResult

    def __post_init__(self) -> None:
        if not self.source_db or not self.target_db:
            raise ValueError("batch UCSC database identifiers must not be empty")
        sha256_hex_from_identifier(self.chain_sha256_identifier)
        if self.chain_resource.sha256 != self.chain_sha256_identifier:
            raise ValueError("batch chain resource identity must match result identity")
        record_ids = [item.record.record_id for item in self.record_assessments]
        if len(record_ids) != len(set(record_ids)):
            raise ValueError("batch result record IDs must be unique")


def run_indexed_chain_batch(
    records: tuple[BatchInputRecord, ...],
    chain_context: CachedUCSCChainResource,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    chain_index: ChainIndex | None,
) -> IndexedChainBatchResult:
    """Project a batch through one prepared chain index and derive relationships.

    A prepared index is mandatory.  Batch execution must never respond to a missing or
    unusable index by reparsing the complete chain once per row (or even once for the
    batch) behind the caller's back.  Index preparation remains an explicit cache-only
    operation at the existing preparation boundary.

    The underlying multi-interval adapter preserves ordinary chain candidate semantics
    and exact chain provenance.  This slice intentionally does not attach net or
    reciprocal-best observations, even when ``chain_context`` names the COMPARATIVE
    publication class.
    """

    if not records:
        raise ValueError("indexed batch execution requires at least one input record")
    if chain_index is None:
        raise ValueError("indexed batch execution requires a prepared chain index")

    candidate_sets = build_ucsc_chain_candidates_for_intervals_from_cached_chain(
        (record.source_interval for record in records),
        chain_context,
        target_assembly=target_assembly,
        alignment_provenance=alignment_provenance,
        chain_index=chain_index,
    )
    if len(candidate_sets) != len(records):
        raise ValueError(
            "batch candidate result count must match the input record count"
        )

    record_assessments = tuple(
        BatchRecordAssessment(record=record, candidates=candidates)
        for record, candidates in zip(records, candidate_sets)
    )
    relationships = build_batch_target_relationships(record_assessments)
    return IndexedChainBatchResult(
        source_db=chain_context.source_db,
        target_db=chain_context.target_db,
        evidence_tier=chain_context.evidence_tier,
        chain_sha256_identifier=chain_context.chain.sha256,
        chain_resource=chain_context.chain,
        alignment_provenance=alignment_provenance,
        record_assessments=record_assessments,
        relationships=relationships,
    )
