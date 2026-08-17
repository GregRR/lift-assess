"""Human-readable and machine-readable assessment reporting.

The renderer consumes an already-computed :class:`~liftassess.models.Assessment`.
It never assigns or modifies a verdict, counts evidence as a score, or infers
biological correctness.  Its job is presentation: expose evidence availability
before interpretation and preserve coordinate semantics explicitly. Human-readable
renderers convert canonical internal coordinates to a labeled 1-based inclusive display
convention; JSON retains canonical 0-based half-open intervals for exact machine use.

A candidate's ``target_interval`` is a bounding span.  For split mappings, the
summary therefore labels that span explicitly and reports the number of mapped
segments rather than implying that every base inside the span aligned.
"""

import json
from typing import assert_never

from .chain import chain_id_from_candidate_id
from .models import (
    Assessment,
    AssessmentDecisionReason,
    ChainGapSummary,
    EvidenceAvailabilityTier,
    EvidenceObservation,
    EvidenceReference,
    EvidenceValue,
    GenomicInterval,
    MappingCoverageSummary,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipSummary,
    Verdict,
)
from .orchestration import UCSCAssessmentReport, UCSCAssessmentResource
from .resource_cache import CachedResource

_BIOLOGICAL_CORRECTNESS_CAVEAT = "This does not establish biological correctness."
_JSON_SCHEMA_VERSION = 1
_JSON_REPORT_TYPE = "liftassess.ucsc_assessment"
_JSON_INTERVAL_COORDINATE_SYSTEM = "0-based-half-open"


def format_display_interval(interval: GenomicInterval) -> str:
    """Format a canonical interval as explicit 1-based inclusive coordinates."""

    if interval.length <= 0:
        raise ValueError("displayed genomic intervals must span at least one base")
    return (
        f"{interval.sequence_name}:{interval.start + 1}-{interval.end} "
        "(1-based inclusive)"
    )


def render_assessment_summary(assessment: Assessment) -> str:
    """Render the concise default v1 assessment summary.

    Evidence availability is intentionally emitted before the verdict so readers do
    not mistake a richer evidence tier for stronger confidence.  The final caveat is
    unconditional because every v1 verdict describes evidentiary support rather than
    biological truth.
    """

    lines = [
        f"Source locus: {format_display_interval(assessment.source_interval)}",
        f"Evidence availability: {_evidence_tier_text(assessment.evidence_tier)}",
        f"Assessment: {_verdict_text(assessment.verdict)}",
    ]

    preferred = _preferred_candidate(assessment)
    if preferred is not None:
        lines.append(f"Preferred candidate: {_candidate_text(preferred)}")
        other_count = len(assessment.candidates) - 1
        if other_count:
            lines.append(f"Other candidates assessed: {other_count}")
    elif len(assessment.candidates) == 1:
        lines.append(f"Candidate: {_candidate_text(assessment.candidates[0])}")
    else:
        lines.append(f"Candidates assessed: {len(assessment.candidates)}")

    lines.append(f"Why: {_decision_reason_text(assessment.decision_reason)}")
    if assessment.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        lines.append(
            "Comparative observations are not assumed to be independent; dependency "
            "provenance is available with --details or --json."
        )
    lines.extend(("", _BIOLOGICAL_CORRECTNESS_CAVEAT))
    return "\n".join(lines)


def _evidence_tier_text(tier: EvidenceAvailabilityTier) -> str:
    if tier is EvidenceAvailabilityTier.COMPARATIVE:
        return "COMPARATIVE — mapping plus comparative evidence available"
    if tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return "LIFTOVER-ONLY — chain mapping evidence only"
    raise ValueError(f"unsupported evidence availability tier: {tier!r}")


def _verdict_text(verdict: Verdict) -> str:
    return verdict.value.replace("_", " ")


def _preferred_candidate(assessment: Assessment) -> NormalizedCandidate | None:
    if assessment.preferred_candidate_id is None:
        return None
    for candidate in assessment.candidates:
        if candidate.candidate_id == assessment.preferred_candidate_id:
            return candidate
    # Assessment itself enforces referential integrity, so reaching this path would
    # indicate that an invalid object bypassed normal dataclass construction.
    raise ValueError("assessment preferred candidate is not present in candidates")


def _candidate_text(candidate: NormalizedCandidate) -> str:
    interval = candidate.target_interval
    coordinate_text = f"{interval.sequence_name}:{interval.start + 1}-{interval.end}"
    details = [
        "1-based inclusive",
        f"{candidate.orientation.value.lower()} orientation",
    ]
    if len(candidate.segments) > 1:
        details.append(f"bounding span of {len(candidate.segments)} mapped segments")
    return f"{coordinate_text} ({'; '.join(details)})"


def _decision_reason_text(reason: AssessmentDecisionReason) -> str:
    """Render the assessor-owned terminal decision reason without re-assessing."""

    match reason:
        case AssessmentDecisionReason.NO_CANDIDATES:
            return "no candidate mapping was generated for the requested locus"
        case AssessmentDecisionReason.LIFTOVER_MULTIPLE_CANDIDATES:
            return "multiple chain-derived candidate mappings remain"
        case AssessmentDecisionReason.LIFTOVER_SINGLE_FULL_MAPPING:
            return "full source-locus mapping coverage"
        case AssessmentDecisionReason.LIFTOVER_SINGLE_PARTIAL_MAPPING:
            return "candidate maps only part of the requested source locus"
        case AssessmentDecisionReason.COMPARATIVE_MULTIPLE_MATERIAL_CANDIDATES:
            return "multiple candidates retain material assessment evidence"
        case AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_FULL:
            return (
                "full source-locus mapping coverage and full reciprocal-best membership"
            )
        case AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_NONE:
            return (
                "full source-locus mapping coverage with reciprocal-best membership "
                "NONE"
            )
        case AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_PARTIAL:
            return (
                "full source-locus mapping coverage with reciprocal-best membership "
                "PARTIAL"
            )
        case AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_PARTIAL:
            return (
                "one material candidate maps only part of the requested source locus; "
                "other raw candidates are not material under the v1 comparative rule"
            )
        case AssessmentDecisionReason.COMPARATIVE_NO_MATERIAL_CANDIDATE:
            return "no candidate retains material comparative evidence"
    assert_never(reason)


def render_assessment_details(report: UCSCAssessmentReport) -> str:
    """Render the auditable human-readable v1 evidence dossier.

    Unlike the concise summary, this renderer exposes every normalized mapping
    segment, evidence observation, retrieval artifact, and provenance dependency in
    the report. Candidate order is preserved for reproducibility but deliberately not
    numbered or described as rank: v1 does not infer preference from encounter order
    or chain score. Supporting/contradicting labels mirror the assessor's categorical
    evidence references and must not be counted as votes or converted to a score.
    """

    assessment = report.assessment
    supporting = set(assessment.supporting_evidence)
    contradicting = set(assessment.contradicting_evidence)
    lines = [
        "Detailed evidence dossier",
        f"UCSC database pair: {report.source_db} -> {report.target_db}",
        f"Source locus: {format_display_interval(assessment.source_interval)}",
        f"Evidence availability: {_evidence_tier_text(assessment.evidence_tier)}",
        f"Assessment: {_verdict_text(assessment.verdict)}",
        f"Decision reason: {assessment.decision_reason.value}",
        f"Candidates assessed: {len(assessment.candidates)}",
    ]

    preferred = _preferred_candidate(assessment)
    if preferred is None:
        lines.append("Preferred candidate: none")
    else:
        lines.append(f"Preferred candidate: {_candidate_heading(preferred)}")

    lines.extend(
        (
            (
                "Role note: supporting/contradicting are categorical roles, not additive "
                "scores."
            ),
            "",
            "Candidates",
            (
                "Candidate order is preserved for reproducibility and does not indicate rank or preference."
            ),
        )
    )

    if not assessment.candidates:
        lines.append("  none")
    for candidate in assessment.candidates:
        lines.extend(
            _candidate_detail_lines(
                candidate,
                supporting=supporting,
                contradicting=contradicting,
            )
        )

    lines.extend(("", "Resources"))
    for assessment_resource in report.resources:
        resource = assessment_resource.resource
        consumption = (
            "consumed" if assessment_resource.consumed_by_engine else "not consumed"
        )
        lines.extend(
            (
                f"{assessment_resource.role.value} [{consumption}]",
                f"  Source URL: {resource.source_url}",
                f"  Cache path: {resource.path}",
                f"  Retrieved at: {resource.retrieved_at}",
                f"  Size: {resource.size_bytes} bytes",
                f"  SHA-256: {resource.sha256}",
                f"  Cache hit at acquisition: {'yes' if resource.cache_hit else 'no'}",
                f"  Provider checksum: {_provider_checksum_text(resource)}",
                f"  General terms: {resource.terms.general_terms_url}",
                f"  Directory terms: {resource.terms.directory_terms_url}",
                "  Evidence provenance: "
                + (
                    assessment_resource.file_provenance.source_id
                    if assessment_resource.file_provenance is not None
                    else "none (artifact was retrieval context only)"
                ),
            )
        )

    provenance_sources = _report_provenance_sources(report)
    lines.extend(("", "Provenance dependencies"))
    for source in provenance_sources:
        lines.append(source.source_id)
        lines.append(f"  Label: {source.label}")
        identifiers = ", ".join(
            f"{identifier.kind.value}={identifier.value}"
            for identifier in source.identifiers
        )
        lines.append(f"  Identifiers: {identifiers or 'none'}")
        parents = ", ".join(parent.source_id for parent in source.derived_from)
        lines.append(f"  Derived from: {parents or 'none'}")

    lines.extend(
        (
            "",
            (
                "Dependency note: provenance edges record shared upstream dependence; "
                "they do not establish independent confirmation."
            ),
            "",
            _BIOLOGICAL_CORRECTNESS_CAVEAT,
        )
    )
    return "\n".join(lines)


def render_assessment_json(report: UCSCAssessmentReport) -> str:
    """Render one deterministic machine-readable v1 UCSC assessment report.

    JSON mirrors the already-computed assessment/report model rather than introducing
    a second interpretation layer. All genomic intervals use liftAssess's canonical
    0-based, half-open coordinates; the schema labels that convention explicitly.
    Candidate array order is retained only for reproducibility, evidence roles remain
    categorical, and provenance edges record dependence rather than independence.
    """

    assessment = report.assessment
    supporting = set(assessment.supporting_evidence)
    contradicting = set(assessment.contradicting_evidence)
    payload: dict[str, object] = {
        "schema_version": _JSON_SCHEMA_VERSION,
        "report_type": _JSON_REPORT_TYPE,
        "semantics": {
            "interval_coordinates": _JSON_INTERVAL_COORDINATE_SYSTEM,
            "candidate_order": "reproducibility_only_not_rank",
            "evidence_roles": "categorical_not_additive",
            "provenance_edges": "dependence_not_independent_confirmation",
        },
        "ucsc_database_pair": {
            "source_db": report.source_db,
            "target_db": report.target_db,
        },
        "assessment": {
            "source_interval": _interval_json(assessment.source_interval),
            "evidence_tier": assessment.evidence_tier.value,
            "verdict": assessment.verdict.value,
            "decision_reason": assessment.decision_reason.value,
            "preferred_candidate_id": assessment.preferred_candidate_id,
            "supporting_evidence": [
                _evidence_reference_json(reference)
                for reference in assessment.supporting_evidence
            ],
            "contradicting_evidence": [
                _evidence_reference_json(reference)
                for reference in assessment.contradicting_evidence
            ],
            "candidates": [
                _candidate_json(
                    candidate,
                    supporting=supporting,
                    contradicting=contradicting,
                )
                for candidate in assessment.candidates
            ],
        },
        "resources": [
            _assessment_resource_json(assessment_resource)
            for assessment_resource in report.resources
        ],
        "provenance": {
            "alignment_source_id": report.alignment_provenance.source_id,
            "sources": [
                _provenance_source_json(source)
                for source in _report_provenance_sources(report)
            ],
        },
        "caveat": _BIOLOGICAL_CORRECTNESS_CAVEAT,
    }
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False)


def _interval_json(interval: GenomicInterval) -> dict[str, object]:
    return {
        "assembly": {
            "name": interval.assembly.name,
            "provider": interval.assembly.provider,
            "accession": interval.assembly.accession,
            "aliases": list(interval.assembly.aliases),
        },
        "sequence_name": interval.sequence_name,
        "start": interval.start,
        "end": interval.end,
        "coordinate_system": _JSON_INTERVAL_COORDINATE_SYSTEM,
    }


def _evidence_reference_json(reference: EvidenceReference) -> dict[str, str]:
    return {
        "candidate_id": reference.candidate_id,
        "observation_id": reference.observation_id,
    }


def _candidate_json(
    candidate: NormalizedCandidate,
    *,
    supporting: set[EvidenceReference],
    contradicting: set[EvidenceReference],
) -> dict[str, object]:
    return {
        "candidate_id": candidate.candidate_id,
        "ucsc_chain_id": _chain_id_from_candidate(candidate),
        "orientation": candidate.orientation.value,
        "target_bounding_interval": _interval_json(candidate.target_interval),
        "mapping_provenance_source_id": candidate.mapping_provenance.source_id,
        "segments": [
            {
                "source_interval": _interval_json(segment.source_interval),
                "target_interval": _interval_json(segment.target_interval),
            }
            for segment in candidate.segments
        ],
        "evidence": [
            _evidence_observation_json(
                candidate.candidate_id,
                observation,
                supporting=supporting,
                contradicting=contradicting,
            )
            for observation in candidate.evidence
        ],
    }


def _evidence_observation_json(
    candidate_id: str,
    observation: EvidenceObservation,
    *,
    supporting: set[EvidenceReference],
    contradicting: set[EvidenceReference],
) -> dict[str, object]:
    reference = EvidenceReference(candidate_id, observation.observation_id)
    return {
        "observation_id": observation.observation_id,
        "kind": observation.kind.value,
        "assessment_role": _evidence_role(
            reference,
            supporting=supporting,
            contradicting=contradicting,
        ),
        "value": _evidence_value_json(observation.value),
        "provenance_source_id": observation.provenance.source_id,
    }


def _evidence_value_json(value: EvidenceValue) -> dict[str, object]:
    if isinstance(value, MappingCoverageSummary):
        return {
            "type": "MAPPING_COVERAGE_SUMMARY",
            "status": value.status.value,
            "covered_source_bases": value.covered_source_bases,
            "source_bases": value.source_bases,
            "uncovered_source_intervals": [
                _interval_json(interval)
                for interval in value.uncovered_source_intervals
            ],
        }

    if isinstance(value, ChainGapSummary):
        return {
            "type": "CHAIN_GAP_SUMMARY",
            "gaps": [
                {
                    "source_boundary_0_based": gap.source_boundary,
                    "source_gap_overlap": (
                        _interval_json(gap.source_gap_overlap)
                        if gap.source_gap_overlap is not None
                        else None
                    ),
                    "target_gap_interval": (
                        _interval_json(gap.target_gap_interval)
                        if gap.target_gap_interval is not None
                        else None
                    ),
                }
                for gap in value.gaps
            ],
        }

    if isinstance(value, NetHierarchySummary):
        return {
            "type": "NET_HIERARCHY_SUMMARY",
            "depth": value.depth,
            "source_fill_interval": _interval_json(value.source_fill_interval),
        }

    if isinstance(value, ReciprocalBestMembershipSummary):
        return {
            "type": "RECIPROCAL_BEST_MEMBERSHIP_SUMMARY",
            "status": value.status.value,
            "resource_completeness": value.resource_completeness.value,
            "chains_examined": value.chains_examined,
            "covered_source_bases": value.covered_source_bases,
            "candidate_source_bases": value.candidate_source_bases,
            "covered_source_intervals": [
                _interval_json(interval) for interval in value.covered_source_intervals
            ],
        }

    if isinstance(value, (str, int, float, bool)):
        return {"type": "SCALAR", "value": value}

    raise TypeError(f"unsupported evidence value for JSON reporting: {type(value)!r}")


def _assessment_resource_json(
    assessment_resource: UCSCAssessmentResource,
) -> dict[str, object]:
    role = assessment_resource.role
    resource = assessment_resource.resource
    file_provenance = assessment_resource.file_provenance
    checksum = resource.provider_checksum
    return {
        "role": role.value,
        "consumed_by_engine": assessment_resource.consumed_by_engine,
        "file_provenance_source_id": (
            file_provenance.source_id if file_provenance is not None else None
        ),
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


def _provenance_source_json(source: ProvenanceSource) -> dict[str, object]:
    return {
        "source_id": source.source_id,
        "label": source.label,
        "identifiers": [
            {"kind": identifier.kind.value, "value": identifier.value}
            for identifier in source.identifiers
        ],
        "derived_from_source_ids": [parent.source_id for parent in source.derived_from],
    }


def _candidate_detail_lines(
    candidate: NormalizedCandidate,
    *,
    supporting: set[EvidenceReference],
    contradicting: set[EvidenceReference],
) -> list[str]:
    lines = [
        _candidate_heading(candidate),
        f"  Candidate ID: {candidate.candidate_id}",
        f"  Target: {_candidate_text(candidate)}",
        f"  Mapping provenance: {candidate.mapping_provenance.source_id}",
        f"  Exact mapped segments ({len(candidate.segments)}):",
    ]
    for segment in candidate.segments:
        lines.append(
            "    "
            f"{format_display_interval(segment.source_interval)} -> "
            f"{format_display_interval(segment.target_interval)}"
        )

    lines.append(f"  Evidence observations ({len(candidate.evidence)}):")
    for observation in candidate.evidence:
        reference = EvidenceReference(
            candidate.candidate_id, observation.observation_id
        )
        role = _evidence_role_text(
            reference,
            supporting=supporting,
            contradicting=contradicting,
        )
        value_lines = _evidence_value_lines(observation)
        lines.append(f"    {observation.kind.value} [{role}]: {value_lines[0]}")
        lines.extend(f"      {line}" for line in value_lines[1:])
        lines.append(f"      provenance: {observation.provenance.source_id}")
    return lines


def _candidate_heading(candidate: NormalizedCandidate) -> str:
    chain_id = _chain_id_from_candidate(candidate)
    if chain_id is None:
        return f"Candidate {candidate.candidate_id}"
    return f"Chain {chain_id}"


def _chain_id_from_candidate(candidate: NormalizedCandidate) -> int | None:
    """Extract the UCSC chain ID embedded by the one v1 candidate engine."""

    return chain_id_from_candidate_id(candidate.candidate_id)


def _evidence_role(
    reference: EvidenceReference,
    *,
    supporting: set[EvidenceReference],
    contradicting: set[EvidenceReference],
) -> str:
    """Return the one categorical assessment role shared by text and JSON output."""

    is_supporting = reference in supporting
    is_contradicting = reference in contradicting
    if is_supporting and is_contradicting:
        return "SUPPORTING_AND_CONTRADICTING"
    if is_supporting:
        return "SUPPORTING"
    if is_contradicting:
        return "CONTRADICTING"
    return "CONTEXT"


def _evidence_role_text(
    reference: EvidenceReference,
    *,
    supporting: set[EvidenceReference],
    contradicting: set[EvidenceReference],
) -> str:
    role = _evidence_role(
        reference,
        supporting=supporting,
        contradicting=contradicting,
    )
    return {
        "SUPPORTING": "supporting",
        "CONTRADICTING": "contradicting",
        "SUPPORTING_AND_CONTRADICTING": "supporting + contradicting",
        "CONTEXT": "context",
    }[role]


def _evidence_value_lines(observation: EvidenceObservation) -> list[str]:
    value = observation.value
    if isinstance(value, MappingCoverageSummary):
        lines = [
            (
                f"{value.status.value}; {value.covered_source_bases}/{value.source_bases} "
                "source bases covered"
            )
        ]
        if value.uncovered_source_intervals:
            lines.append(
                "uncovered source intervals: "
                + ", ".join(
                    format_display_interval(interval)
                    for interval in value.uncovered_source_intervals
                )
            )
        else:
            lines.append("uncovered source intervals: none")
        return lines

    if isinstance(value, ChainGapSummary):
        lines = [f"{len(value.gaps)} chain gap(s) through the requested locus"]
        for gap in value.gaps:
            source_gap = (
                format_display_interval(gap.source_gap_overlap)
                if gap.source_gap_overlap is not None
                else "none"
            )
            target_gap = (
                format_display_interval(gap.target_gap_interval)
                if gap.target_gap_interval is not None
                else "none"
            )
            lines.append(
                f"source boundary={gap.source_boundary} (0-based boundary); "
                f"source gap={source_gap}; target gap={target_gap}"
            )
        return lines

    if isinstance(value, NetHierarchySummary):
        return [
            (
                f"depth={value.depth}; fill span="
                f"{format_display_interval(value.source_fill_interval)}"
            )
        ]

    if isinstance(value, ReciprocalBestMembershipSummary):
        lines = [
            (
                f"{value.status.value}; {value.covered_source_bases}/"
                f"{value.candidate_source_bases} candidate mapped source bases covered; "
                f"completeness={value.resource_completeness.value}; "
                f"chains examined={value.chains_examined}"
            )
        ]
        if value.covered_source_intervals:
            lines.append(
                "covered source intervals: "
                + ", ".join(
                    format_display_interval(interval)
                    for interval in value.covered_source_intervals
                )
            )
        else:
            lines.append("covered source intervals: none")
        return lines

    return [str(value)]


def _provider_checksum_text(resource: CachedResource) -> str:
    # Kept structural rather than provider-specific so future resource providers can
    # preserve their own checksum metadata without the renderer interpreting it.
    checksum = resource.provider_checksum
    if checksum is None:
        return "none"
    return f"{checksum.algorithm.value}:{checksum.value} (from {checksum.source_url})"


def _report_provenance_sources(
    report: UCSCAssessmentReport,
) -> tuple[ProvenanceSource, ...]:
    roots: list[ProvenanceSource] = [report.alignment_provenance]
    for assessment_resource in report.resources:
        if assessment_resource.file_provenance is not None:
            roots.append(assessment_resource.file_provenance)
    for candidate in report.assessment.candidates:
        roots.append(candidate.mapping_provenance)
        roots.extend(observation.provenance for observation in candidate.evidence)

    by_id: dict[str, ProvenanceSource] = {}
    pending = list(roots)
    while pending:
        source = pending.pop()
        existing = by_id.get(source.source_id)
        if existing is not None:
            if _provenance_definition(existing) != _provenance_definition(source):
                raise ValueError(
                    "provenance source ID refers to conflicting source definitions"
                )
            continue
        by_id[source.source_id] = source
        pending.extend(source.derived_from)
    return tuple(by_id[source_id] for source_id in sorted(by_id))


def _provenance_definition(source: ProvenanceSource) -> tuple[object, ...]:
    """Return one cycle-safe structural definition for a provenance source ID."""

    return (
        source.label,
        source.identifiers,
        tuple(parent.source_id for parent in source.derived_from),
    )
