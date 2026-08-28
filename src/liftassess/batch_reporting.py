"""Human and machine reporting for indexed batch assessment."""

import json

from .batch import (
    BatchRecordPointContext,
    BatchTargetRelationship,
    BatchTargetRelationshipKind,
)
from .batch_execution import IndexedChainBatchResult
from .models import (
    EvidenceAvailabilityTier,
    EvidenceKind,
    GenomicInterval,
    NormalizedCandidate,
    ReciprocalBestMembershipSummary,
)
from .query_context import PointQueryContextResult, QueryContextState
from .reporting import (
    assembly_json_payload,
    candidate_json_payload,
    interval_json_payload,
    provenance_source_json_payload,
)

_BIOLOGICAL_CORRECTNESS_CAVEAT = "This does not establish biological correctness."
_JSON_SCHEMA_VERSION = 2
_JSON_INTERVAL_COORDINATE_SYSTEM = "0-based-half-open"
_DEFAULT_BATCH_RECORD_PREVIEW_LIMIT = 20
_DEFAULT_INLINE_PROJECTION_LIMIT = 4


def render_indexed_chain_batch_summary(result: IndexedChainBatchResult) -> str:
    """Render a compact factual summary for indexed batch assessment."""

    record_count = len(result.record_assessments)
    records_with_projection = sum(
        bool(item.candidates) for item in result.record_assessments
    )
    candidate_count = sum(len(item.candidates) for item in result.record_assessments)
    relationship_counts = _batch_relationship_counts(result)
    point_count, context_run_count, context_not_run_count = _point_context_counts(
        result
    )
    lines = [
        "* BATCH CHAIN PROJECTIONS *",
        "Input records:",
        f"    {record_count}",
        "Records with indexed candidate(s):",
        f"    {records_with_projection}",
        "Records with zero indexed candidates:",
        f"    {record_count - records_with_projection}",
        "Candidate projections:",
        f"    {candidate_count}",
        "Input target relationships:",
        (
            "    none"
            if not result.relationships.relationships
            else "    "
            + ", ".join(
                f"{kind.value}={relationship_counts[kind]}"
                for kind in BatchTargetRelationshipKind
                if relationship_counts[kind]
            )
        ),
    ]
    if point_count:
        lines.extend(
            [
                "Point context:",
                (
                    f"    requested={result.point_context_window_bases} bp; "
                    f"point records={point_count}; run={context_run_count}; "
                    f"not run={context_not_run_count}"
                ),
                "Point-context target relationships:",
                _render_context_relationship_counts(result),
            ]
        )
    lines.extend(
        [
            "Evidence:",
            _batch_evidence_summary(result),
            "Chain resource:",
            f"    {result.chain_sha256_identifier}",
            "Scope:",
            _batch_scope_summary(result),
        ]
    )

    preview = result.record_assessments[:_DEFAULT_BATCH_RECORD_PREVIEW_LIMIT]
    if preview:
        lines.append("Records:")
        for index, assessment in enumerate(preview):
            record = assessment.record
            context_result = result.point_context_records[index].context_result
            label = f" [{record.label}]" if record.label is not None else ""
            lines.append(f"    {record.record_id}{label}")
            lines.append(
                "        Source: " + _format_half_open_interval(record.source_interval)
            )
            if not assessment.candidates:
                lines.append("        Projections: 0")
            else:
                lines.append(f"        Projections: {len(assessment.candidates)}")
                for candidate in assessment.candidates[
                    :_DEFAULT_INLINE_PROJECTION_LIMIT
                ]:
                    lines.append(
                        f"            {candidate.candidate_id}: target bounding span "
                        + _format_half_open_interval(candidate.target_interval)
                        + f"; orientation={candidate.orientation.value}"
                        + _candidate_comparative_suffix(candidate)
                    )
                omitted = len(assessment.candidates) - _DEFAULT_INLINE_PROJECTION_LIMIT
                if omitted > 0:
                    lines.append(f"            ... {omitted} more projection(s)")
            _append_point_context_preview(lines, context_result)

    omitted_records = record_count - len(preview)
    if omitted_records > 0:
        lines.append(f"    ... {omitted_records} more record(s); use --json for all")
    return "\n".join(lines)


def render_indexed_chain_batch_json(result: IndexedChainBatchResult) -> str:
    """Render schema-v2 machine-readable indexed batch output."""

    chain_provenance_source_id = f"file:{result.chain_sha256_identifier}"
    payload: dict[str, object] = {
        "schema_version": _JSON_SCHEMA_VERSION,
        "report_type": "liftassess.ucsc_batch_result",
        "semantics": {
            "interval_coordinates": _JSON_INTERVAL_COORDINATE_SYSTEM,
            "input_record_order": "input_order",
            "candidate_order": "reproducibility_only_not_rank",
            "result_dimensions": "orthogonal_not_votes",
            "batch_relationships": "cross_record_exact_mapped_target_segments",
            "point_context_relationships": (
                "same_relationship_geometry_at_explicit_point_context_scale"
            ),
            "evidence_scope": (
                "indexed_chain_plus_shared_comparative"
                if result.comparative_evidence_consumed
                else "indexed_chain_only"
            ),
            "provenance_edges": "dependence_not_independent_confirmation",
        },
        "ucsc_database_pair": {
            "source_db": result.source_db,
            "target_db": result.target_db,
        },
        "evidence": {
            "resource_publication_class": result.evidence_tier.value,
            "assessment_scope": (
                "CHAIN_NET_RECIPROCAL_BEST"
                if result.comparative_evidence_consumed
                else "CHAIN_ONLY"
            ),
            "comparative_net_reciprocal_best": (
                "ASSESSED_FOR_SUBMITTED_RECORDS"
                if result.comparative_evidence_consumed
                else "NOT_USED_NO_SUBMITTED_CANDIDATES"
                if result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
                else "NOT_ASSESSED"
            ),
        },
        "records": [
            {
                "record_id": assessment.record.record_id,
                "source_line_number": assessment.record.source_line_number,
                "label": assessment.record.label,
                "source_interval": interval_json_payload(
                    assessment.record.source_interval
                ),
                "candidates": [
                    candidate_json_payload(candidate)
                    for candidate in assessment.candidates
                ],
            }
            for assessment in result.record_assessments
        ],
        "relationships": [
            {
                "kind": relationship.kind.value,
                "left_record_id": relationship.left_record_id,
                "left_candidate_id": relationship.left_candidate_id,
                "right_record_id": relationship.right_record_id,
                "right_candidate_id": relationship.right_candidate_id,
                "target_assembly": assembly_json_payload(relationship.target_assembly),
                "target_sequence_name": relationship.target_sequence_name,
                "overlap_intervals": [
                    interval_json_payload(interval)
                    for interval in relationship.overlap_intervals
                ],
            }
            for relationship in result.relationships.relationships
        ],
        "point_context": {
            "requested_window_bases": result.point_context_window_bases,
            "records": [
                _point_context_record_json(
                    item, requested_window_bases=result.point_context_window_bases
                )
                for item in result.point_context_records
            ],
            "relationships": [
                _context_relationship_json_payload(relationship)
                for relationship in result.point_context_relationships.relationships
            ],
        },
        "resource": _batch_chain_resource_json(result, chain_provenance_source_id),
        "comparative_resources": _batch_comparative_resources_json(result),
        "provenance": {
            "alignment_source_id": result.alignment_provenance.source_id,
            "sources": _batch_provenance_sources(result, chain_provenance_source_id),
        },
        "scope": {
            "authoritative_source_sequence_preflight": "NOT_ASSESSED",
            "cross_record_target_relationships": "ASSESSED",
            "reverse_mapping": "NOT_ASSESSED",
            "point_context": _point_context_scope(result),
            "point_context_comparative_evidence": "NOT_ASSESSED",
            "filtered_all_chain_comparison": "NOT_ASSESSED",
            "comparative_relationship_interpretation": "NOT_ASSESSED",
            "target_role": "NOT_ASSESSED",
            "named_variant_identity": "NOT_ASSESSED",
            "gene_transcript_identity": "NOT_ASSESSED",
        },
        "caveat": _BIOLOGICAL_CORRECTNESS_CAVEAT,
    }
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False)


def _batch_relationship_counts(
    result: IndexedChainBatchResult,
) -> dict[BatchTargetRelationshipKind, int]:
    counts = {kind: 0 for kind in BatchTargetRelationshipKind}
    for relationship in result.relationships.relationships:
        counts[relationship.kind] += 1
    return counts


def _point_context_counts(result: IndexedChainBatchResult) -> tuple[int, int, int]:
    point_count = 0
    run_count = 0
    not_run_count = 0
    for item in result.point_context_records:
        context = item.context_result
        if context is None:
            continue
        point_count += 1
        if context.check_state is QueryContextState.RUN:
            run_count += 1
        else:
            not_run_count += 1
    return point_count, run_count, not_run_count


def _render_context_relationship_counts(result: IndexedChainBatchResult) -> str:
    counts = {kind: 0 for kind in BatchTargetRelationshipKind}
    for relationship in result.point_context_relationships.relationships:
        counts[relationship.kind] += 1
    if not result.point_context_relationships.relationships:
        return "    none"
    parts: list[str] = []
    exact_count = counts[BatchTargetRelationshipKind.EXACT_TARGET_COLLISION]
    if exact_count:
        parts.append(f"NEIGHBORHOOD_LEVEL_TARGET_COLLISION={exact_count}")
    overlap_count = counts[BatchTargetRelationshipKind.OVERLAPPING_TARGET_PROJECTIONS]
    if overlap_count:
        parts.append(f"OVERLAPPING_TARGET_PROJECTIONS={overlap_count}")
    return "    " + ", ".join(parts)


def _append_point_context_preview(
    lines: list[str],
    context: PointQueryContextResult | None,
) -> None:
    if context is None:
        return
    if context.check_state is not QueryContextState.RUN:
        reason = context.not_run_reason
        reason_text = reason.value if reason is not None else "UNKNOWN"
        lines.append(f"        Point context: NOT RUN ({reason_text})")
        return
    tested = context.tested_source_interval
    if tested is None:
        raise ValueError("completed batch point context requires a tested interval")
    lines.append(
        "        Point context: "
        + _format_half_open_interval(tested)
        + f"; projections={len(context.candidates)}"
    )
    for candidate in context.candidates[:_DEFAULT_INLINE_PROJECTION_LIMIT]:
        covered = sum(segment.source_interval.length for segment in candidate.segments)
        lines.append(
            f"            {candidate.candidate_id}: target bounding span "
            + _format_half_open_interval(candidate.target_interval)
            + f"; source coverage={covered}/{tested.length}"
        )
    omitted = len(context.candidates) - _DEFAULT_INLINE_PROJECTION_LIMIT
    if omitted > 0:
        lines.append(f"            ... {omitted} more context projection(s)")


def _point_context_record_json(
    item: BatchRecordPointContext,
    *,
    requested_window_bases: int,
) -> dict[str, object]:
    context = item.context_result
    if context is None:
        return {
            "record_id": item.record.record_id,
            "state": "NOT_APPLICABLE",
            "requested_window_bases": requested_window_bases,
            "not_run_reason": "SOURCE_INTERVAL_IS_NOT_ONE_BASE",
        }
    payload: dict[str, object] = {
        "record_id": item.record.record_id,
        "state": context.check_state.value,
        "requested_window_bases": context.requested_window_bases,
        "not_run_reason": (
            context.not_run_reason.value if context.not_run_reason is not None else None
        ),
    }
    if context.check_state is QueryContextState.RUN:
        tested = context.tested_source_interval
        if tested is None:
            raise ValueError("completed batch point context requires a tested interval")
        payload["tested_source_interval"] = interval_json_payload(tested)
        payload["candidates"] = [
            candidate_json_payload(candidate) for candidate in context.candidates
        ]
    return payload


def _context_relationship_json_payload(
    relationship: BatchTargetRelationship,
) -> dict[str, object]:
    kind = relationship.kind.value
    if relationship.kind is BatchTargetRelationshipKind.EXACT_TARGET_COLLISION:
        kind = "NEIGHBORHOOD_LEVEL_TARGET_COLLISION"
    return {
        "kind": kind,
        "left_record_id": relationship.left_record_id,
        "left_candidate_id": relationship.left_candidate_id,
        "right_record_id": relationship.right_record_id,
        "right_candidate_id": relationship.right_candidate_id,
        "target_assembly": assembly_json_payload(relationship.target_assembly),
        "target_sequence_name": relationship.target_sequence_name,
        "overlap_intervals": [
            interval_json_payload(interval)
            for interval in relationship.overlap_intervals
        ],
    }


def _point_context_scope(result: IndexedChainBatchResult) -> str:
    point_count, run_count, not_run_count = _point_context_counts(result)
    if point_count == 0:
        return "NOT_APPLICABLE"
    if not_run_count == 0:
        return "ASSESSED_FOR_ALL_POINT_RECORDS"
    if run_count == 0:
        return "NOT_RUN_FOR_POINT_RECORDS"
    return "PARTIALLY_ASSESSED_FOR_POINT_RECORDS"


def _batch_scope_summary(result: IndexedChainBatchResult) -> str:
    point_count, run_count, not_run_count = _point_context_counts(result)
    if not point_count:
        context_text = (
            "Point context is not applicable because the batch has no 1-bp rows."
        )
    elif not not_run_count:
        context_text = (
            f"Automatic {result.point_context_window_bases}-bp point context is "
            "assessed from the same prepared chain index for every 1-bp row."
        )
    else:
        context_text = (
            f"Automatic {result.point_context_window_bases}-bp point context ran for "
            f"{run_count}/{point_count} point rows; unavailable indexed source bounds "
            "are reported per record."
        )
    return (
        "    Cross-record target relationships are derived from exact mapped target "
        "segments. "
        + context_text
        + (
            " Submitted-row net and reciprocal-best evidence are assessed once across "
            "the batch; point-context candidates remain chain-only."
            if result.comparative_evidence_consumed
            else (
                " Net/reciprocal-best resources were available but not consumed "
                "because no submitted candidate was generated; point-context remains "
                "chain-only."
            )
            if result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
            else (
                " Net/reciprocal-best evidence is not assessed for "
                "LIFTOVER-ONLY batches."
            )
        )
        + (
            " Filtered-vs-all-chain comparison and categorical comparative "
            "relationship interpretation are not assessed in this batch slice."
            if result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
            else ""
        )
        + " Authoritative assembly-sequence name/alias preflight, reverse, "
        "target-role, named-variant, and gene/transcript evidence are not assessed "
        "in this batch slice."
    )


def _batch_evidence_summary(result: IndexedChainBatchResult) -> str:
    tier = result.evidence_tier.value.replace("_", "-")
    if result.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        if result.comparative_evidence_consumed:
            return (
                f"    {tier} publication class; indexed all-chain plus one shared net "
                "scan and one shared reciprocal-best-chain scan for submitted rows"
            )
        return (
            f"    {tier} publication class; indexed all-chain; net/reciprocal-best "
            "not used because no submitted candidate was generated"
        )
    return f"    {tier} publication class; indexed chain only"


def _candidate_comparative_suffix(candidate: NormalizedCandidate) -> str:
    classifications = [
        str(observation.value)
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.NET_CLASSIFICATION
    ]
    reciprocal = [
        observation.value
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
        and isinstance(observation.value, ReciprocalBestMembershipSummary)
    ]
    if not classifications and not reciprocal:
        return ""
    parts: list[str] = []
    if classifications:
        parts.append("net=" + ",".join(classifications))
    if reciprocal:
        parts.append("reciprocal-best=" + reciprocal[-1].status.value)
    return "; " + "; ".join(parts)


def _batch_comparative_resources_json(
    result: IndexedChainBatchResult,
) -> list[dict[str, object]]:
    if result.evidence_tier is not EvidenceAvailabilityTier.COMPARATIVE:
        return []
    resources = (
        ("NET", result.net_resource),
        ("RECIPROCAL_BEST_CHAIN", result.reciprocal_best_chain_resource),
    )
    payload: list[dict[str, object]] = []
    for role, resource in resources:
        if resource is None:
            raise ValueError("COMPARATIVE batch result is missing a consumed resource")
        payload.append(
            {
                "role": role,
                "consumed_by_engine": result.comparative_evidence_consumed,
                "source_url": resource.source_url,
                "cache_path": str(resource.path),
                "retrieved_at": resource.retrieved_at,
                "size_bytes": resource.size_bytes,
                "sha256": resource.sha256,
            }
        )
    return payload


def _batch_provenance_sources(
    result: IndexedChainBatchResult,
    chain_provenance_source_id: str,
) -> list[dict[str, object]]:
    sources: list[dict[str, object]] = [
        provenance_source_json_payload(result.alignment_provenance),
        {
            "source_id": chain_provenance_source_id,
            "label": f"UCSC {result.source_db}→{result.target_db} chain resource",
            "identifiers": [
                {"kind": "SHA256", "value": result.chain_sha256_identifier}
            ],
            "derived_from_source_ids": [result.alignment_provenance.source_id],
        },
    ]
    if result.comparative_evidence_consumed:
        for label, resource in (
            ("net resource", result.net_resource),
            ("reciprocal-best chain resource", result.reciprocal_best_chain_resource),
        ):
            if resource is None:
                raise ValueError(
                    "COMPARATIVE batch result is missing comparative provenance"
                )
            sources.append(
                {
                    "source_id": f"file:{resource.sha256}",
                    "label": f"UCSC {result.source_db}→{result.target_db} {label}",
                    "identifiers": [{"kind": "SHA256", "value": resource.sha256}],
                    "derived_from_source_ids": [result.alignment_provenance.source_id],
                }
            )
    return sources


def _format_half_open_interval(interval: GenomicInterval) -> str:
    return (
        f"{interval.sequence_name}:{interval.start}-{interval.end} (0-based half-open)"
    )


def _batch_chain_resource_json(
    result: IndexedChainBatchResult,
    file_provenance_source_id: str,
) -> dict[str, object]:
    resource = result.chain_resource
    checksum = resource.provider_checksum
    return {
        "role": "CHAIN",
        "consumed_by_engine": True,
        "file_provenance_source_id": file_provenance_source_id,
        "source_url": resource.source_url,
        "cache_path": str(resource.path),
        "retrieved_at": resource.retrieved_at,
        "size_bytes": resource.size_bytes,
        "sha256": resource.sha256,
        "cache_hit_at_acquisition": resource.cache_hit,
        "provider_checksum": (
            {
                "algorithm": checksum.algorithm.value,
                "value": checksum.value,
                "source_url": checksum.source_url,
            }
            if checksum is not None
            else None
        ),
        "terms": {
            "resource_class": resource.terms.resource_class.value,
            "general_terms_url": resource.terms.general_terms_url,
            "directory_terms_url": resource.terms.directory_terms_url,
            "restricted_liftover_chain": resource.terms.restricted_liftover_chain,
        },
    }
