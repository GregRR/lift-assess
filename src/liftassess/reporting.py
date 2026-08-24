"""Human-readable and machine-readable factual result reporting.

Both renderers consume the same derived result profile. Exact candidate geometry,
evidence, resource identity, and provenance remain available alongside the profile;
no renderer assigns an aggregate verdict, confidence score, candidate rank, or
biological correctness claim.
"""

import json

from .chain import chain_id_from_candidate_id
from .models import (
    AssemblyIdentifier,
    ChainGapSummary,
    EvidenceAvailabilityTier,
    EvidenceObservation,
    EvidenceValue,
    GenomicInterval,
    MappingCoverageSummary,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipSummary,
)
from .orchestration import UCSCAssessmentReport, UCSCAssessmentResource
from .resource_cache import CachedResource
from .result_profile import (
    CandidateResultProfile,
    CandidateReverseMappingProfile,
    FactualHeadline,
    ResultProfile,
    SourceCoverageState,
)
from .reverse_mapping import (
    CandidateReverseMappingResult,
    ReverseCheckState,
    ReverseRelationshipState,
)

_BIOLOGICAL_CORRECTNESS_CAVEAT = "This does not establish biological correctness."
_JSON_SCHEMA_VERSION = 2
_JSON_REPORT_TYPE = "liftassess.ucsc_result"
_JSON_INTERVAL_COORDINATE_SYSTEM = "0-based-half-open"
_DEFAULT_INLINE_PROJECTION_LIMIT = 4


def format_display_interval(interval: GenomicInterval) -> str:
    """Format a canonical interval as explicit 1-based inclusive coordinates."""

    if interval.length <= 0:
        raise ValueError("displayed genomic intervals must span at least one base")
    return (
        f"{interval.sequence_name}:{interval.start + 1}-{interval.end} "
        "(1-based inclusive)"
    )


def render_assessment_summary(report: UCSCAssessmentReport) -> str:
    """Render the progressive-disclosure default factual result summary."""

    profile = report.result_profile
    lines = [
        _headline_text(profile.headline),
        f"Source: {format_display_interval(profile.source_interval)}",
    ]

    if not profile.candidate_profiles:
        lines.append("Chain projections: 0")
    elif len(profile.candidate_profiles) == 1:
        candidate_profile = profile.candidate_profiles[0]
        candidate = report.candidates[0]
        lines.extend(_single_candidate_summary_lines(candidate, candidate_profile))
    else:
        lines.extend(_multiple_candidate_summary_lines(report, profile))

    lines.append(f"Evidence: {_evidence_summary(profile)}")
    lines.append(f"Interpretation: {profile.interpretation}")
    lines.append(
        "Scope: coordinate projection/structure assessed; named-variant and "
        "gene/transcript identity not assessed."
    )
    lines.append(
        "Details: use --details for the full profile/evidence or --json for schema v2."
    )
    lines.append(_BIOLOGICAL_CORRECTNESS_CAVEAT)
    return "\n".join(lines)


def _single_candidate_summary_lines(
    candidate: NormalizedCandidate,
    profile: CandidateResultProfile,
) -> list[str]:
    lines = [
        f"Source coverage: {profile.covered_source_bases}/{profile.source_bases} source bases",
        f"Target: {_candidate_text(candidate, profile)}",
    ]
    if profile.fragmented or profile.coverage_state is SourceCoverageState.PARTIAL:
        lines.append(f"Geometric mapped segments: {profile.geometric_segment_count}")
    if profile.uncovered_source_intervals:
        lines.append(
            "Uncovered source: "
            + ", ".join(
                format_display_interval(interval)
                for interval in profile.uncovered_source_intervals
            )
        )
    if profile.target_gap_intervals:
        lines.append(
            "Target gaps: "
            + ", ".join(
                format_display_interval(interval)
                for interval in profile.target_gap_intervals
            )
        )
    lines.append(f"Reverse mapping: {_reverse_summary_text(profile.reverse_mapping)}")
    return lines


def _multiple_candidate_summary_lines(
    report: UCSCAssessmentReport,
    profile: ResultProfile,
) -> list[str]:
    lines = [
        (
            "Maximum candidate source coverage: "
            f"{profile.maximum_candidate_covered_source_bases}/"
            f"{profile.source_bases} bases"
        ),
        f"Chain projections: {len(profile.candidate_profiles)}",
        "Projection order: reproducibility only; not rank.",
    ]
    if (
        profile.union_covered_source_bases
        != profile.maximum_candidate_covered_source_bases
    ):
        lines.append(
            "Source bases represented across all projections: "
            f"{profile.union_covered_source_bases}/{profile.source_bases}"
        )
    if len(profile.candidate_profiles) <= _DEFAULT_INLINE_PROJECTION_LIMIT:
        lines.append("Projection details:")
        for candidate, candidate_profile in zip(
            report.candidates,
            profile.candidate_profiles,
            strict=True,
        ):
            lines.append(
                "  - "
                f"{_candidate_text(candidate, candidate_profile)}; "
                f"coverage {candidate_profile.covered_source_bases}/"
                f"{candidate_profile.source_bases}; "
                f"geometric segments {candidate_profile.geometric_segment_count}; "
                f"reverse {_reverse_summary_text(candidate_profile.reverse_mapping)}"
            )
        return lines

    segment_counts = [
        candidate.geometric_segment_count for candidate in profile.candidate_profiles
    ]
    target_sequence_count = len(
        {candidate.target_interval.sequence_name for candidate in report.candidates}
    )
    lines.extend(
        (
            f"Target sequences represented: {target_sequence_count}",
            f"Projection orientations: {profile.orientation.value}",
            _segment_count_summary(segment_counts),
            (
                "Projections at maximum source coverage: "
                f"{len(profile.maximum_coverage_candidate_ids)}"
            ),
            _reverse_set_summary(profile),
            (
                "Projection details: omitted from default output for this candidate "
                "set; use --details or --json for every projection."
            ),
        )
    )
    return lines


def _segment_count_summary(segment_counts: list[int]) -> str:
    minimum = min(segment_counts)
    maximum = max(segment_counts)
    if minimum == maximum:
        return f"Geometric mapped segments per projection: {minimum}"
    return f"Geometric mapped segments per projection: {minimum}-{maximum}"


def _reverse_summary_text(profile: CandidateReverseMappingProfile) -> str:
    if profile.check_state is ReverseCheckState.NOT_RUN:
        return "not run"
    if profile.check_state is ReverseCheckState.UNAVAILABLE:
        return "unavailable from the current prepared reverse resources"

    relationship = profile.relationship
    if relationship is ReverseRelationshipState.NO_PROJECTION:
        return "completed; no reverse chain projection"
    if relationship is ReverseRelationshipState.ELSEWHERE_ONLY:
        return "returns only to a different source locus"
    if relationship is ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE:
        return "returns to the original source locus and elsewhere"
    if relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY:
        if profile.exact_original_geometry_return:
            return "exactly reconstructs the original aligned source geometry"
        assert profile.original_source_covered_bases is not None
        return (
            "returns only to the original source locus; recovered "
            f"{profile.original_source_covered_bases}/{profile.original_source_bases} "
            "aligned source bases"
        )
    raise ValueError("completed reverse mapping requires a relationship state")


def _reverse_set_summary(profile: ResultProfile) -> str:
    states = [candidate.reverse_mapping for candidate in profile.candidate_profiles]
    if not states:
        return "Reverse mapping: not run"
    check_state = states[0].check_state
    if check_state is ReverseCheckState.NOT_RUN:
        return "Reverse mapping: not run"
    if check_state is ReverseCheckState.UNAVAILABLE:
        return (
            "Reverse mapping: unavailable from the current prepared reverse resources"
        )

    exact = sum(bool(item.exact_original_geometry_return) for item in states)
    original_only_nonexact = sum(
        item.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY
        and not item.exact_original_geometry_return
        for item in states
    )
    elsewhere_only = sum(
        item.relationship is ReverseRelationshipState.ELSEWHERE_ONLY for item in states
    )
    mixed = sum(
        item.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE
        for item in states
    )
    no_projection = sum(
        item.relationship is ReverseRelationshipState.NO_PROJECTION for item in states
    )
    return (
        "Reverse mapping: "
        f"{exact} exact original-geometry return(s); "
        f"{original_only_nonexact} other original-only return(s); "
        f"{elsewhere_only} elsewhere-only; {mixed} original+elsewhere; "
        f"{no_projection} no-projection"
    )


def _headline_text(headline: FactualHeadline) -> str:
    return headline.value.replace("_", " ")


def _evidence_summary(profile: ResultProfile) -> str:
    roles = ", ".join(profile.consumed_resource_roles) or "none"
    if profile.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return f"LIFTOVER-ONLY — consumed {roles}; chain mapping evidence only"
    if profile.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        return (
            f"COMPARATIVE — consumed {roles}; UCSC-derived observations may share "
            "upstream alignment lineage and are not independent votes"
        )
    raise ValueError(
        f"unsupported evidence availability tier: {profile.evidence_tier!r}"
    )


def _candidate_text(
    candidate: NormalizedCandidate,
    profile: CandidateResultProfile,
) -> str:
    interval = candidate.target_interval
    coordinate_text = f"{interval.sequence_name}:{interval.start + 1}-{interval.end}"
    details = [
        "1-based inclusive",
        f"{candidate.orientation.value.lower()} orientation",
    ]
    if profile.geometric_segment_count > 1:
        details.append(
            "bounding span of "
            f"{profile.geometric_segment_count} geometric mapped segments"
        )
    return f"{coordinate_text} ({'; '.join(details)})"


def render_assessment_details(report: UCSCAssessmentReport) -> str:
    """Render the complete factual profile, evidence, resources, and provenance."""

    profile = report.result_profile
    lines = [
        "Detailed factual result dossier",
        f"UCSC database pair: {report.source_db} -> {report.target_db}",
        f"Source locus: {format_display_interval(report.source_interval)}",
        f"Headline: {_headline_text(profile.headline)}",
        f"Interpretation: {profile.interpretation}",
        f"Input validity preflight: {profile.input_validity.value}",
        f"Projection count: {profile.projection_count.value}",
        f"Projection orientation: {profile.orientation.value}",
        (
            "Maximum candidate source coverage: "
            f"{profile.maximum_candidate_covered_source_bases}/{profile.source_bases}"
        ),
        (
            "Union source coverage across candidates: "
            f"{profile.union_covered_source_bases}/{profile.source_bases}"
        ),
        f"Evidence availability: {profile.evidence_tier.value.replace('_', '-')}",
        "Consumed resource roles: "
        + (", ".join(profile.consumed_resource_roles) or "none"),
        "",
        "Current scope boundaries",
        f"  Target role: {profile.scope.target_role.value}",
        f"  Actual reverse mapping: {profile.scope.reverse_result.value}",
        f"  Point/neighborhood context: {profile.scope.query_context.value}",
        (
            "  Comparative relationship synthesis: "
            f"{profile.scope.comparative_relationship.value}"
        ),
        f"  Batch relationships: {profile.scope.batch_relationship.value}",
        f"  Typed external context: {profile.scope.external_context.value}",
        "  Named variant / rsID identity: NOT ASSESSED",
        "  Gene / transcript identity: NOT ASSESSED",
        "  File / downstream workflow: NOT ASSESSED",
        "",
        "Candidates",
        (
            "Candidate order is preserved for reproducibility and does not indicate "
            "rank or preference."
        ),
    ]

    if not report.candidates:
        lines.append("  none")
    for candidate, candidate_profile in zip(
        report.candidates,
        profile.candidate_profiles,
        strict=True,
    ):
        lines.extend(_candidate_detail_lines(candidate, candidate_profile))

    if report.reverse_mapping_results is not None:
        lines.extend(("", "Reverse mapping results"))
        if not report.reverse_mapping_results:
            lines.append("  none")
        for result in report.reverse_mapping_results:
            lines.extend(_reverse_mapping_detail_lines(result))

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
                (
                    "  File provenance: "
                    + (
                        assessment_resource.file_provenance.source_id
                        if assessment_resource.file_provenance is not None
                        else "none"
                    )
                ),
            )
        )

    if report.reverse_mapping_resource is not None:
        lines.extend(("", "Reverse mapping resource"))
        reverse_resource = report.reverse_mapping_resource
        resource = reverse_resource.resource
        lines.extend(
            (
                (
                    f"{report.target_db}->{report.source_db} "
                    f"{reverse_resource.role.value} [consumed]"
                ),
                f"  Source URL: {resource.source_url}",
                f"  Cache path: {resource.path}",
                f"  Retrieved at: {resource.retrieved_at}",
                f"  Size: {resource.size_bytes} bytes",
                f"  SHA-256: {resource.sha256}",
                f"  Provider checksum: {_provider_checksum_text(resource)}",
                (
                    "  File provenance: "
                    + (
                        reverse_resource.file_provenance.source_id
                        if reverse_resource.file_provenance is not None
                        else "none"
                    )
                ),
            )
        )

    lines.extend(("", "Provenance dependency graph"))
    for source in _report_provenance_sources(report):
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
    """Render the schema-v2 factual result report."""

    payload: dict[str, object] = {
        "schema_version": _JSON_SCHEMA_VERSION,
        "report_type": _JSON_REPORT_TYPE,
        "semantics": {
            "interval_coordinates": _JSON_INTERVAL_COORDINATE_SYSTEM,
            "candidate_order": "reproducibility_only_not_rank",
            "result_dimensions": "orthogonal_not_votes",
            "provenance_edges": "dependence_not_independent_confirmation",
        },
        "ucsc_database_pair": {
            "source_db": report.source_db,
            "target_db": report.target_db,
        },
        "source_assembly": _assembly_json(report.source_interval.assembly),
        "target_assembly": _assembly_json(report.target_assembly),
        "source_interval": _interval_json(report.source_interval),
        "result_profile": _result_profile_json(report.result_profile),
        "candidates": [_candidate_json(candidate) for candidate in report.candidates],
        "reverse_mapping": _reverse_mapping_json(report),
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


def _result_profile_json(profile: ResultProfile) -> dict[str, object]:
    return {
        "input_validity": profile.input_validity.value,
        "headline": profile.headline.value,
        "interpretation": profile.interpretation,
        "projection_count": profile.projection_count.value,
        "orientation": profile.orientation.value,
        "source_coverage": {
            "state": profile.source_coverage.value,
            "maximum_candidate_covered_source_bases": (
                profile.maximum_candidate_covered_source_bases
            ),
            "source_bases": profile.source_bases,
            "maximum_coverage_candidate_ids": list(
                profile.maximum_coverage_candidate_ids
            ),
            "union_covered_source_bases": profile.union_covered_source_bases,
        },
        "evidence": {
            "tier": profile.evidence_tier.value,
            "consumed_resource_roles": list(profile.consumed_resource_roles),
        },
        "candidate_profiles": [
            _candidate_profile_json(candidate)
            for candidate in profile.candidate_profiles
        ],
        "scope": {
            "target_role": profile.scope.target_role.value,
            "actual_reverse_mapping": profile.scope.reverse_result.value,
            "query_context": profile.scope.query_context.value,
            "comparative_relationship": profile.scope.comparative_relationship.value,
            "batch_relationship": profile.scope.batch_relationship.value,
            "external_context": profile.scope.external_context.value,
            "named_variant_identity_assessed": (
                profile.scope.named_variant_identity_assessed
            ),
            "gene_transcript_identity_assessed": (
                profile.scope.gene_transcript_identity_assessed
            ),
            "downstream_workflow_assessed": profile.scope.downstream_workflow_assessed,
        },
    }


def _candidate_profile_json(profile: CandidateResultProfile) -> dict[str, object]:
    return {
        "candidate_id": profile.candidate_id,
        "source_coverage": {
            "state": profile.coverage_state.value,
            "covered_source_bases": profile.covered_source_bases,
            "source_bases": profile.source_bases,
            "uncovered_source_intervals": [
                _interval_json(interval)
                for interval in profile.uncovered_source_intervals
            ],
            "largest_uncovered_source_span_bases": (
                profile.largest_uncovered_source_span_bases
            ),
        },
        "geometry": {
            "exact_mapped_segment_count": profile.exact_mapped_segment_count,
            "geometric_segment_count": profile.geometric_segment_count,
            "fragmented": profile.fragmented,
            "target_discontinuous": profile.target_discontinuous,
            "target_bounding_span": _interval_json(profile.target_bounding_span),
            "source_gap_intervals": [
                _interval_json(interval) for interval in profile.source_gap_intervals
            ],
            "target_gap_intervals": [
                _interval_json(interval) for interval in profile.target_gap_intervals
            ],
            "largest_source_gap_bases": profile.largest_source_gap_bases,
            "largest_target_gap_bases": profile.largest_target_gap_bases,
        },
        "orientation": profile.orientation.value,
        "reverse_mapping": _reverse_profile_json(profile.reverse_mapping),
    }


def _reverse_profile_json(
    profile: CandidateReverseMappingProfile,
) -> dict[str, object]:
    return {
        "check_state": profile.check_state.value,
        "relationship": (
            profile.relationship.value if profile.relationship is not None else None
        ),
        "original_source_bases": profile.original_source_bases,
        "original_source_covered_bases": profile.original_source_covered_bases,
        "original_source_coverage": (
            profile.original_source_coverage.value
            if profile.original_source_coverage is not None
            else None
        ),
        "exact_original_geometry_return": profile.exact_original_geometry_return,
        "reverse_projection_count": profile.reverse_projection_count,
        "segments_with_reverse_projection": profile.segments_with_reverse_projection,
        "queried_target_segments": [
            _interval_json(interval) for interval in profile.queried_target_segments
        ],
    }


def _reverse_mapping_json(report: UCSCAssessmentReport) -> dict[str, object]:
    results = report.reverse_mapping_results
    return {
        "check_state": report.result_profile.scope.reverse_result.value,
        "reverse_database_pair": (
            {
                "source_db": report.target_db,
                "target_db": report.source_db,
            }
            if results is not None
            and any(result.check_state is ReverseCheckState.RUN for result in results)
            else None
        ),
        "resource": (
            _assessment_resource_json(report.reverse_mapping_resource)
            if report.reverse_mapping_resource is not None
            else None
        ),
        "candidate_results": (
            [_candidate_reverse_mapping_json(result) for result in results]
            if results is not None
            else []
        ),
    }


def _candidate_reverse_mapping_json(
    result: CandidateReverseMappingResult,
) -> dict[str, object]:
    return {
        "forward_candidate_id": result.forward_candidate_id,
        "check_state": result.check_state.value,
        "relationship": (
            result.relationship.value if result.relationship is not None else None
        ),
        "original_source_bases": result.original_source_bases,
        "original_source_covered_bases": (
            result.original_source_covered_bases
            if result.check_state is ReverseCheckState.RUN
            else None
        ),
        "original_source_coverage": (
            result.original_source_coverage.value
            if result.check_state is ReverseCheckState.RUN
            else None
        ),
        "exact_original_geometry_return": (
            result.exact_original_geometry_return
            if result.check_state is ReverseCheckState.RUN
            else None
        ),
        "queried_target_segments": [
            _interval_json(interval) for interval in result.queried_target_segments
        ],
        "segment_results": [
            {
                "queried_target_segment": _interval_json(
                    segment_result.queried_target_segment
                ),
                "expected_original_source_segment": _interval_json(
                    segment_result.expected_original_source_segment
                ),
                "reverse_candidates": [
                    _candidate_json(candidate)
                    for candidate in segment_result.candidates
                ],
            }
            for segment_result in result.segment_results
        ],
    }


def _assembly_json(assembly: AssemblyIdentifier) -> dict[str, object]:
    # Kept as a helper so assembly serialization is identical inside/outside intervals.
    return {
        "name": assembly.name,
        "provider": assembly.provider,
        "accession": assembly.accession,
        "aliases": list(assembly.aliases),
    }


def _interval_json(interval: GenomicInterval) -> dict[str, object]:
    return {
        "assembly": _assembly_json(interval.assembly),
        "sequence_name": interval.sequence_name,
        "start": interval.start,
        "end": interval.end,
        "coordinate_system": _JSON_INTERVAL_COORDINATE_SYSTEM,
    }


def _candidate_json(candidate: NormalizedCandidate) -> dict[str, object]:
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
            _evidence_observation_json(observation)
            for observation in candidate.evidence
        ],
    }


def _evidence_observation_json(observation: EvidenceObservation) -> dict[str, object]:
    return {
        "observation_id": observation.observation_id,
        "kind": observation.kind.value,
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
    profile: CandidateResultProfile,
) -> list[str]:
    lines = [
        _candidate_heading(candidate),
        f"  Candidate ID: {candidate.candidate_id}",
        f"  Target: {_candidate_text(candidate, profile)}",
        (
            f"  Source coverage: {profile.covered_source_bases}/{profile.source_bases} "
            f"({profile.coverage_state.value})"
        ),
        f"  Geometric mapped segments: {profile.geometric_segment_count}",
        f"  Fragmented: {'yes' if profile.fragmented else 'no'}",
        f"  Target discontinuous: {'yes' if profile.target_discontinuous else 'no'}",
        (
            "  Largest uncovered source span: "
            f"{profile.largest_uncovered_source_span_bases} bases"
        ),
        f"  Largest source chain gap: {profile.largest_source_gap_bases} bases",
        f"  Largest target gap: {profile.largest_target_gap_bases} bases",
        f"  Mapping provenance: {candidate.mapping_provenance.source_id}",
        (
            "  Exact chain-derived mapped segments "
            f"({profile.exact_mapped_segment_count}):"
        ),
    ]
    for segment in candidate.segments:
        lines.append(
            "    "
            f"{format_display_interval(segment.source_interval)} -> "
            f"{format_display_interval(segment.target_interval)}"
        )

    if profile.uncovered_source_intervals:
        lines.append("  Uncovered source intervals:")
        lines.extend(
            f"    {format_display_interval(interval)}"
            for interval in profile.uncovered_source_intervals
        )
    else:
        lines.append("  Uncovered source intervals: none")

    if profile.source_gap_intervals:
        lines.append("  Source chain-gap intervals:")
        lines.extend(
            f"    {format_display_interval(interval)}"
            for interval in profile.source_gap_intervals
        )
    else:
        lines.append("  Source chain-gap intervals: none")

    if profile.target_gap_intervals:
        lines.append("  Target gap intervals:")
        lines.extend(
            f"    {format_display_interval(interval)}"
            for interval in profile.target_gap_intervals
        )
    else:
        lines.append("  Target gap intervals: none")

    reverse = profile.reverse_mapping
    lines.append(f"  Reverse mapping check: {reverse.check_state.value}")
    if reverse.check_state is ReverseCheckState.RUN:
        assert reverse.relationship is not None
        assert reverse.original_source_covered_bases is not None
        assert reverse.original_source_coverage is not None
        assert reverse.exact_original_geometry_return is not None
        assert reverse.reverse_projection_count is not None
        assert reverse.segments_with_reverse_projection is not None
        lines.extend(
            (
                f"  Reverse relationship: {reverse.relationship.value}",
                (
                    "  Reverse original-source coverage: "
                    f"{reverse.original_source_covered_bases}/"
                    f"{reverse.original_source_bases} "
                    f"({reverse.original_source_coverage.value})"
                ),
                (
                    "  Exact original aligned geometry reconstructed: "
                    f"{'yes' if reverse.exact_original_geometry_return else 'no'}"
                ),
                f"  Reverse projections: {reverse.reverse_projection_count}",
                (
                    "  Forward target segments with reverse projection: "
                    f"{reverse.segments_with_reverse_projection}/"
                    f"{len(reverse.queried_target_segments)}"
                ),
            )
        )

    lines.append(f"  Evidence observations ({len(candidate.evidence)}):")
    for observation in candidate.evidence:
        value_lines = _evidence_value_lines(observation)
        lines.append(f"    {observation.kind.value}: {value_lines[0]}")
        lines.extend(f"      {line}" for line in value_lines[1:])
        lines.append(f"      provenance: {observation.provenance.source_id}")
    return lines


def _reverse_mapping_detail_lines(
    result: CandidateReverseMappingResult,
) -> list[str]:
    lines = [
        f"Candidate {result.forward_candidate_id}",
        f"  Check state: {result.check_state.value}",
    ]
    if result.check_state is not ReverseCheckState.RUN:
        return lines

    assert result.relationship is not None
    lines.extend(
        (
            f"  Relationship: {result.relationship.value}",
            (
                "  Original aligned source coverage: "
                f"{result.original_source_covered_bases}/"
                f"{result.original_source_bases} "
                f"({result.original_source_coverage.value})"
            ),
            (
                "  Exact original aligned geometry reconstructed: "
                f"{'yes' if result.exact_original_geometry_return else 'no'}"
            ),
        )
    )
    for index, segment_result in enumerate(result.segment_results, start=1):
        lines.append(
            f"  Segment {index} reverse query: "
            f"{format_display_interval(segment_result.queried_target_segment)}"
        )
        expected_source = format_display_interval(
            segment_result.expected_original_source_segment
        )
        lines.append(f"    Expected original source: {expected_source}")
        if not segment_result.candidates:
            lines.append("    Reverse projections: none")
            continue
        lines.append(f"    Reverse projections: {len(segment_result.candidates)}")
        for candidate in segment_result.candidates:
            lines.append(
                "      "
                f"{candidate.candidate_id}: "
                f"{format_display_interval(candidate.target_interval)}"
            )
    return lines


def _candidate_heading(candidate: NormalizedCandidate) -> str:
    chain_id = _chain_id_from_candidate(candidate)
    if chain_id is None:
        return f"Candidate {candidate.candidate_id}"
    return f"Chain {chain_id}"


def _chain_id_from_candidate(candidate: NormalizedCandidate) -> int | None:
    return chain_id_from_candidate_id(candidate.candidate_id)


def _evidence_value_lines(observation: EvidenceObservation) -> list[str]:
    value = observation.value
    if isinstance(value, MappingCoverageSummary):
        lines = [
            (
                f"{value.status.value}; {value.covered_source_bases}/"
                f"{value.source_bases} source bases covered"
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
                f"depth={value.depth}; "
                f"fill span={format_display_interval(value.source_fill_interval)}"
            )
        ]

    if isinstance(value, ReciprocalBestMembershipSummary):
        lines = [
            (
                f"{value.status.value}; {value.covered_source_bases}/"
                f"{value.candidate_source_bases} "
                "candidate mapped source bases covered; "
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
    for candidate in report.candidates:
        roots.append(candidate.mapping_provenance)
        roots.extend(observation.provenance for observation in candidate.evidence)
    if report.reverse_alignment_provenance is not None:
        roots.append(report.reverse_alignment_provenance)
    if (
        report.reverse_mapping_resource is not None
        and report.reverse_mapping_resource.file_provenance is not None
    ):
        roots.append(report.reverse_mapping_resource.file_provenance)
    if report.reverse_mapping_results is not None:
        for result in report.reverse_mapping_results:
            for segment_result in result.segment_results:
                for candidate in segment_result.candidates:
                    roots.append(candidate.mapping_provenance)
                    roots.extend(
                        observation.provenance for observation in candidate.evidence
                    )

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
    return (
        source.label,
        source.identifiers,
        tuple(parent.source_id for parent in source.derived_from),
    )
