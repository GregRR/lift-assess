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
    BatchRecordPointContext,
    BatchRelationshipResult,
    build_batch_target_relationships,
)
from .chain_index import ChainIndex
from .models import (
    AssemblyIdentifier,
    EvidenceAvailabilityTier,
    GenomicInterval,
    ProvenanceSource,
)
from .query_context import (
    DEFAULT_POINT_CONTEXT_BASES,
    PointQueryContextResult,
    QueryContextNotRunReason,
    QueryContextState,
    build_centered_point_context_interval,
    point_context_not_run,
)
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
    point_context_window_bases: int
    point_context_records: tuple[BatchRecordPointContext, ...]
    point_context_relationships: BatchRelationshipResult

    def __post_init__(self) -> None:
        if not self.source_db or not self.target_db:
            raise ValueError("batch UCSC database identifiers must not be empty")
        sha256_hex_from_identifier(self.chain_sha256_identifier)
        if self.chain_resource.sha256 != self.chain_sha256_identifier:
            raise ValueError("batch chain resource identity must match result identity")
        record_ids = [item.record.record_id for item in self.record_assessments]
        if len(record_ids) != len(set(record_ids)):
            raise ValueError("batch result record IDs must be unique")
        if (
            self.point_context_window_bases < 3
            or self.point_context_window_bases % 2 == 0
        ):
            raise ValueError(
                "batch point-context window must be an odd value of at least 3"
            )
        context_record_ids = [
            item.record.record_id for item in self.point_context_records
        ]
        if context_record_ids != record_ids:
            raise ValueError(
                "batch point-context records must preserve the input record order"
            )


def run_indexed_chain_batch(
    records: tuple[BatchInputRecord, ...],
    chain_context: CachedUCSCChainResource,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    chain_index: ChainIndex | None,
    point_context_window_bases: int = DEFAULT_POINT_CONTEXT_BASES,
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
    point_context_records, point_context_relationships = _run_batch_point_context(
        records,
        chain_context,
        target_assembly=target_assembly,
        alignment_provenance=alignment_provenance,
        chain_index=chain_index,
        requested_window_bases=point_context_window_bases,
    )
    return IndexedChainBatchResult(
        source_db=chain_context.source_db,
        target_db=chain_context.target_db,
        evidence_tier=chain_context.evidence_tier,
        chain_sha256_identifier=chain_context.chain.sha256,
        chain_resource=chain_context.chain,
        alignment_provenance=alignment_provenance,
        record_assessments=record_assessments,
        relationships=relationships,
        point_context_window_bases=point_context_window_bases,
        point_context_records=point_context_records,
        point_context_relationships=point_context_relationships,
    )


def _run_batch_point_context(
    records: tuple[BatchInputRecord, ...],
    chain_context: CachedUCSCChainResource,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    chain_index: ChainIndex,
    requested_window_bases: int,
) -> tuple[tuple[BatchRecordPointContext, ...], BatchRelationshipResult]:
    """Assess point neighborhoods from the same prepared index in one batch slice."""

    runnable_positions: list[int] = []
    runnable_intervals: list[GenomicInterval] = []
    contexts: list[BatchRecordPointContext | None] = [None] * len(records)

    for position, record in enumerate(records):
        if record.source_interval.length != 1:
            contexts[position] = BatchRecordPointContext(
                record=record, context_result=None
            )
            continue

        source_bound = chain_index.source_sequence_query_bound(
            record.source_interval.sequence_name
        )
        if source_bound is None:
            contexts[position] = BatchRecordPointContext(
                record=record,
                context_result=point_context_not_run(
                    requested_window_bases=requested_window_bases,
                    reason=QueryContextNotRunReason.SOURCE_BOUNDS_UNAVAILABLE,
                ),
            )
            continue

        context_interval = build_centered_point_context_interval(
            record.source_interval,
            requested_window_bases=requested_window_bases,
            source_sequence_query_bound=source_bound,
        )
        runnable_positions.append(position)
        runnable_intervals.append(context_interval)

    context_candidate_sets = (
        build_ucsc_chain_candidates_for_intervals_from_cached_chain(
            runnable_intervals,
            chain_context,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
            chain_index=chain_index,
        )
        if runnable_intervals
        else ()
    )
    if len(context_candidate_sets) != len(runnable_positions):
        raise ValueError(
            "batch point-context candidate result count must match runnable points"
        )

    for position, candidates in zip(runnable_positions, context_candidate_sets):
        record = records[position]
        source_bound = chain_index.source_sequence_query_bound(
            record.source_interval.sequence_name
        )
        if source_bound is None:
            raise ValueError(
                "batch point-context source bound changed during execution"
            )
        tested_interval = build_centered_point_context_interval(
            record.source_interval,
            requested_window_bases=requested_window_bases,
            source_sequence_query_bound=source_bound,
        )
        contexts[position] = BatchRecordPointContext(
            record=record,
            context_result=PointQueryContextResult(
                check_state=QueryContextState.RUN,
                requested_window_bases=requested_window_bases,
                tested_source_interval=tested_interval,
                candidates=candidates,
            ),
        )

    finalized = tuple(context for context in contexts if context is not None)
    if len(finalized) != len(records):
        raise ValueError("batch point-context result count must match input records")

    context_assessments: list[BatchRecordAssessment] = []
    for item in finalized:
        result = item.context_result
        if result is None or result.check_state is not QueryContextState.RUN:
            continue
        completed_interval = result.tested_source_interval
        if completed_interval is None:
            raise ValueError("completed batch point context requires a tested interval")
        context_assessments.append(
            BatchRecordAssessment(
                record=BatchInputRecord(
                    record_id=item.record.record_id,
                    source_interval=completed_interval,
                    source_line_number=item.record.source_line_number,
                    label=item.record.label,
                ),
                candidates=result.candidates,
            )
        )

    return finalized, build_batch_target_relationships(tuple(context_assessments))
