"""Human and machine reporting for indexed chain-only batch assessment."""

import json

from .batch import BatchTargetRelationshipKind
from .batch_execution import IndexedChainBatchResult
from .models import GenomicInterval
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
    """Render a compact factual summary for indexed chain-only batch assessment."""

    record_count = len(result.record_assessments)
    records_with_projection = sum(
        bool(item.candidates) for item in result.record_assessments
    )
    candidate_count = sum(len(item.candidates) for item in result.record_assessments)
    relationship_counts = _batch_relationship_counts(result)
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
        "Target relationships:",
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
        "Evidence:",
        (
            f"    {result.evidence_tier.value.replace('_', '-')} publication class; "
            "indexed chain only"
        ),
        "Chain resource:",
        f"    {result.chain_sha256_identifier}",
        "Scope:",
        (
            "    Cross-record target relationships are derived from exact mapped "
            "target segments. Authoritative assembly-sequence name/alias preflight, "
            "net/reciprocal-best, reverse, point-context, target-role, named-variant, "
            "and gene/transcript evidence are not assessed in this batch slice."
        ),
    ]

    preview = result.record_assessments[:_DEFAULT_BATCH_RECORD_PREVIEW_LIMIT]
    if preview:
        lines.append("Records:")
        for assessment in preview:
            record = assessment.record
            label = f" [{record.label}]" if record.label is not None else ""
            lines.append(f"    {record.record_id}{label}")
            lines.append(
                "        Source: " + _format_half_open_interval(record.source_interval)
            )
            if not assessment.candidates:
                lines.append("        Projections: 0")
                continue
            lines.append(f"        Projections: {len(assessment.candidates)}")
            for candidate in assessment.candidates[:_DEFAULT_INLINE_PROJECTION_LIMIT]:
                lines.append(
                    f"            {candidate.candidate_id}: target bounding span "
                    + _format_half_open_interval(candidate.target_interval)
                    + f"; orientation={candidate.orientation.value}"
                )
            omitted = len(assessment.candidates) - _DEFAULT_INLINE_PROJECTION_LIMIT
            if omitted > 0:
                lines.append(f"            ... {omitted} more projection(s)")

    omitted_records = record_count - len(preview)
    if omitted_records > 0:
        lines.append(f"    ... {omitted_records} more record(s); use --json for all")
    return "\n".join(lines)


def render_indexed_chain_batch_json(result: IndexedChainBatchResult) -> str:
    """Render schema-v2 machine-readable indexed chain-only batch output."""

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
            "evidence_scope": "indexed_chain_only",
            "provenance_edges": "dependence_not_independent_confirmation",
        },
        "ucsc_database_pair": {
            "source_db": result.source_db,
            "target_db": result.target_db,
        },
        "evidence": {
            "resource_publication_class": result.evidence_tier.value,
            "assessment_scope": "CHAIN_ONLY",
            "comparative_net_reciprocal_best": "NOT_ASSESSED",
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
        "resource": _batch_chain_resource_json(result, chain_provenance_source_id),
        "provenance": {
            "alignment_source_id": result.alignment_provenance.source_id,
            "sources": [
                provenance_source_json_payload(result.alignment_provenance),
                {
                    "source_id": chain_provenance_source_id,
                    "label": (
                        f"UCSC {result.source_db}→{result.target_db} chain resource"
                    ),
                    "identifiers": [
                        {"kind": "SHA256", "value": result.chain_sha256_identifier}
                    ],
                    "derived_from_source_ids": [result.alignment_provenance.source_id],
                },
            ],
        },
        "scope": {
            "authoritative_source_sequence_preflight": "NOT_ASSESSED",
            "cross_record_target_relationships": "ASSESSED",
            "reverse_mapping": "NOT_ASSESSED",
            "point_context": "NOT_ASSESSED",
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
