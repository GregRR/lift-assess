"""Concise human-readable assessment reporting.

The renderer consumes an already-computed :class:`~liftassess.models.Assessment`.
It never assigns or modifies a verdict, counts evidence as a score, or infers
biological correctness.  Its job is presentation: expose evidence availability
before interpretation, convert canonical internal coordinates back to an explicit
1-based inclusive display convention, and summarize only the narrow evidence states
that drive the v1 assessor.

A candidate's ``target_interval`` is a bounding span.  For split mappings, the
summary therefore labels that span explicitly and reports the number of mapped
segments rather than implying that every base inside the span aligned.
"""

from .models import (
    Assessment,
    ChainGapSummary,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    EvidenceReference,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipSummary,
    Verdict,
)
from .orchestration import UCSCAssessmentReport
from .resource_cache import CachedResource

_BIOLOGICAL_CORRECTNESS_CAVEAT = "This does not establish biological correctness."


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

    lines.append(f"Why: {_verdict_basis_text(assessment)}")
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


def _verdict_basis_text(assessment: Assessment) -> str:
    if not assessment.candidates:
        return "no candidate mapping was generated for the requested locus"

    if len(assessment.candidates) > 1:
        if assessment.verdict is Verdict.CONTESTED:
            return "multiple candidates retain material assessment evidence"
        if assessment.verdict is Verdict.INDETERMINATE:
            return (
                "available evidence does not materially distinguish the candidate "
                "mappings"
            )

    candidate = _preferred_candidate(assessment)
    if candidate is None and len(assessment.candidates) == 1:
        candidate = assessment.candidates[0]

    if candidate is not None:
        coverage = _mapping_coverage(candidate)
        reciprocal = _reciprocal_best(candidate)

        if assessment.verdict is Verdict.WELL_SUPPORTED:
            if reciprocal is None:
                return "full source-locus mapping coverage"
            return (
                "full source-locus mapping coverage and full reciprocal-best membership"
            )

        if coverage is not None and coverage.status is MappingCoverageStatus.PARTIAL:
            return "candidate maps only part of the requested source locus"

        if reciprocal is not None:
            return (
                "full source-locus mapping coverage with reciprocal-best membership "
                f"{reciprocal.status.value}"
            )

    if assessment.verdict is Verdict.CONTESTED:
        return "available verdict-driving evidence materially disagrees"
    return "available verdict-driving evidence is insufficient or non-discriminating"


def _mapping_coverage(candidate: NormalizedCandidate) -> MappingCoverageSummary | None:
    for observation in candidate.evidence:
        if observation.kind is EvidenceKind.MAPPING_COVERAGE and isinstance(
            observation.value, MappingCoverageSummary
        ):
            return observation.value
    return None


def _reciprocal_best(
    candidate: NormalizedCandidate,
) -> ReciprocalBestMembershipSummary | None:
    for observation in candidate.evidence:
        if observation.kind is EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP and isinstance(
            observation.value, ReciprocalBestMembershipSummary
        ):
            return observation.value
    return None


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

    _, separator, chain_text = candidate.candidate_id.rpartition(":chain:")
    if not separator:
        return None
    try:
        chain_id = int(chain_text)
    except ValueError:
        return None
    return chain_id if chain_id >= 0 else None


def _evidence_role_text(
    reference: EvidenceReference,
    *,
    supporting: set[EvidenceReference],
    contradicting: set[EvidenceReference],
) -> str:
    is_supporting = reference in supporting
    is_contradicting = reference in contradicting
    if is_supporting and is_contradicting:
        return "supporting + contradicting"
    if is_supporting:
        return "supporting"
    if is_contradicting:
        return "contradicting"
    return "context"


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
