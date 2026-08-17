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
    EvidenceAvailabilityTier,
    EvidenceKind,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    NormalizedCandidate,
    ReciprocalBestMembershipSummary,
    Verdict,
)

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
