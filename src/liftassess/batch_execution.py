"""Index-required batch execution over exact cached UCSC resources.

Candidate generation always uses one prepared region-addressable chain index and never
falls back to a whole-chain traversal. For a complete COMPARATIVE bundle, submitted-row
net and reciprocal-best evidence are attached with one shared pass over each resource.
Automatic point-context candidates intentionally remain forward-chain-only.
"""

from dataclasses import dataclass

from .assembly_metadata import (
    SourceIntervalPreflightResult,
    SourceIntervalPreflightState,
)
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
    ProvenanceIdentifierKind,
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
from .resource_cache import (
    CachedResource,
    CachedUCSCChainResource,
    CachedUCSCResourceBundle,
)
from .resource_files import (
    ResourceReadProgressCallback,
    attach_ucsc_comparative_evidence_to_candidate_sets_from_cached_bundle,
    build_ucsc_chain_candidates_for_intervals_from_cached_chain,
)
from .resource_identity import sha256_hex_from_identifier


@dataclass(frozen=True)
class IndexedChainBatchResult:
    """Indexed batch result for one exact UCSC publication class."""

    source_db: str
    target_db: str
    evidence_tier: EvidenceAvailabilityTier
    chain_sha256_identifier: str
    chain_resource: CachedResource
    net_resource: CachedResource | None
    reciprocal_best_chain_resource: CachedResource | None
    comparative_evidence_consumed: bool
    alignment_provenance: ProvenanceSource
    record_assessments: tuple[BatchRecordAssessment, ...]
    relationships: BatchRelationshipResult
    point_context_window_bases: int
    point_context_records: tuple[BatchRecordPointContext, ...]
    point_context_relationships: BatchRelationshipResult
    source_preflights: tuple[SourceIntervalPreflightResult, ...] | None = None
    source_preflight_resources: tuple[CachedResource, ...] = ()

    def __post_init__(self) -> None:
        if not self.source_db or not self.target_db:
            raise ValueError("batch UCSC database identifiers must not be empty")
        sha256_hex_from_identifier(self.chain_sha256_identifier)
        if self.chain_resource.sha256 != self.chain_sha256_identifier:
            raise ValueError("batch chain resource identity must match result identity")
        comparative_resources = (
            self.net_resource,
            self.reciprocal_best_chain_resource,
        )
        if self.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
            if any(resource is None for resource in comparative_resources):
                raise ValueError(
                    "COMPARATIVE batch results require net and reciprocal-best "
                    "resources"
                )
        elif any(resource is not None for resource in comparative_resources):
            raise ValueError(
                "LIFTOVER_ONLY batch results cannot carry comparative resources"
            )
        if (
            self.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
            and self.comparative_evidence_consumed
        ):
            raise ValueError(
                "LIFTOVER_ONLY batch results cannot consume comparative evidence"
            )
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
        _validate_batch_source_preflights(
            tuple(item.record for item in self.record_assessments),
            self.source_preflights,
            self.source_preflight_resources,
        )


def run_indexed_chain_batch(
    records: tuple[BatchInputRecord, ...],
    chain_context: CachedUCSCChainResource,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    chain_index: ChainIndex | None,
    comparative_bundle: CachedUCSCResourceBundle | None = None,
    progress_callback: ResourceReadProgressCallback | None = None,
    point_context_window_bases: int = DEFAULT_POINT_CONTEXT_BASES,
    source_preflights: tuple[SourceIntervalPreflightResult, ...] | None = None,
    source_preflight_resources: tuple[CachedResource, ...] = (),
) -> IndexedChainBatchResult:
    """Project a batch through one prepared chain index and derive relationships.

    A prepared index is mandatory.  Batch execution must never respond to a missing or
    unusable index by reparsing the complete chain once per row (or even once for the
    batch) behind the caller's back.  Index preparation remains an explicit cache-only
    operation at the existing preparation boundary.

    The underlying multi-interval adapter preserves ordinary chain candidate semantics
    and exact chain provenance. When ``comparative_bundle`` is supplied, its ordinary
    net and reciprocal-best chain are each consumed once across the complete submitted
    candidate collection. Point-context candidates remain forward-chain-only in this
    slice and do not inherit those observations.
    """

    if not records:
        raise ValueError("indexed batch execution requires at least one input record")
    _validate_batch_source_preflights(
        records,
        source_preflights,
        source_preflight_resources,
    )
    if chain_index is None:
        raise ValueError("indexed batch execution requires a prepared chain index")
    if comparative_bundle is not None:
        if chain_context.evidence_tier is not EvidenceAvailabilityTier.COMPARATIVE:
            raise ValueError(
                "comparative batch evidence requires a COMPARATIVE chain context"
            )
        if (
            comparative_bundle.source_db != chain_context.source_db
            or comparative_bundle.target_db != chain_context.target_db
            or comparative_bundle.evidence_tier
            is not EvidenceAvailabilityTier.COMPARATIVE
            or comparative_bundle.chain.sha256 != chain_context.chain.sha256
        ):
            raise ValueError(
                "comparative batch bundle must match the selected exact chain resource"
            )
    elif chain_context.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        raise ValueError(
            "COMPARATIVE batch execution requires the complete cached "
            "comparative bundle"
        )

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
    comparative_evidence_consumed = bool(
        comparative_bundle is not None and any(candidate_sets)
    )
    if comparative_evidence_consumed:
        assert comparative_bundle is not None
        candidate_sets = (
            attach_ucsc_comparative_evidence_to_candidate_sets_from_cached_bundle(
                candidate_sets,
                comparative_bundle,
                alignment_provenance=alignment_provenance,
                progress_callback=progress_callback,
            )
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
        source_preflights=source_preflights,
    )
    return IndexedChainBatchResult(
        source_db=chain_context.source_db,
        target_db=chain_context.target_db,
        evidence_tier=chain_context.evidence_tier,
        chain_sha256_identifier=chain_context.chain.sha256,
        chain_resource=chain_context.chain,
        net_resource=(
            comparative_bundle.net if comparative_bundle is not None else None
        ),
        reciprocal_best_chain_resource=(
            comparative_bundle.reciprocal_best_chain
            if comparative_bundle is not None
            else None
        ),
        comparative_evidence_consumed=comparative_evidence_consumed,
        alignment_provenance=alignment_provenance,
        record_assessments=record_assessments,
        relationships=relationships,
        point_context_window_bases=point_context_window_bases,
        point_context_records=point_context_records,
        point_context_relationships=point_context_relationships,
        source_preflights=source_preflights,
        source_preflight_resources=source_preflight_resources,
    )


def _validate_batch_source_preflights(
    records: tuple[BatchInputRecord, ...],
    source_preflights: tuple[SourceIntervalPreflightResult, ...] | None,
    source_preflight_resources: tuple[CachedResource, ...],
) -> None:
    """Require valid, record-aligned preflight facts when supplied."""

    if source_preflights is None:
        if source_preflight_resources:
            raise ValueError(
                "batch source preflight resources require source preflight facts"
            )
        return
    if len(source_preflights) != len(records):
        raise ValueError(
            "batch source preflight result count must match the input record count"
        )

    expected_sha256: set[str] = set()
    for record, preflight in zip(records, source_preflights, strict=True):
        if preflight.source_interval != record.source_interval:
            raise ValueError(
                "batch source preflight intervals must preserve input record order"
            )
        if preflight.state is not SourceIntervalPreflightState.VALID:
            raise ValueError(
                "indexed batch scientific assessment requires valid source preflight"
            )
        expected_sha256.update(
            identifier.value
            for source in preflight.provenance_sources
            for identifier in source.identifiers
            if identifier.kind is ProvenanceIdentifierKind.SHA256
        )

    actual_sha256 = {resource.sha256 for resource in source_preflight_resources}
    if actual_sha256 != expected_sha256:
        raise ValueError(
            "batch source preflight resources must match authoritative provenance"
        )


def _batch_source_sequence_bound(
    record: BatchInputRecord,
    *,
    position: int,
    chain_index: ChainIndex,
    source_preflights: tuple[SourceIntervalPreflightResult, ...] | None,
) -> int | None:
    """Return authoritative source bounds when preflight facts are available."""

    if source_preflights is not None:
        preflight = source_preflights[position]
        if preflight.source_interval != record.source_interval:
            raise ValueError(
                "batch source preflight intervals must preserve input record order"
            )
        return preflight.sequence_length
    return chain_index.source_sequence_query_bound(record.source_interval.sequence_name)


def _run_batch_point_context(
    records: tuple[BatchInputRecord, ...],
    chain_context: CachedUCSCChainResource,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    chain_index: ChainIndex,
    requested_window_bases: int,
    source_preflights: tuple[SourceIntervalPreflightResult, ...] | None,
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

        source_bound = _batch_source_sequence_bound(
            record,
            position=position,
            chain_index=chain_index,
            source_preflights=source_preflights,
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
        source_bound = _batch_source_sequence_bound(
            record,
            position=position,
            chain_index=chain_index,
            source_preflights=source_preflights,
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
