"""Deterministic v1 assessment rules for normalized liftOver candidates.

The assessor consumes normalized candidates and their provenance-aware evidence. It
never inspects how the candidates were generated and deliberately does not calculate
a numeric support score.

v1 verdict-driving evidence is intentionally narrow:

* source-locus mapping coverage (``FULL`` vs ``PARTIAL``); and
* reciprocal-best membership (``FULL``/``PARTIAL``/``NONE``) when the evidence tier
  is comparative.

Raw chain score, net ``ali``/``qDup``, net classification, and hierarchy remain
important report context, but v1 does not assign them a monotonic support direction or
an arbitrary threshold. Keeping them out of verdict arithmetic avoids turning
length-dependent or fill-scoped measurements into an undocumented confidence score.

Reciprocal-best membership may share the same upstream alignment as the chain mapping.
The rules therefore use its categorical self-consistency state directly; they never
count it as an independent vote or increase support because several observations share
the same provenance. Because v1 has no additive evidence voting, the verdict logic does
not branch on provenance independence itself; provenance remains attached for audit and
reporting. Any future rule that aggregates evidence sources must explicitly account for
shared upstream provenance before treating observations as independent support.
"""

from dataclasses import dataclass

from .models import (
    Assessment,
    AssessmentDecisionReason,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    EvidenceReference,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    NormalizedCandidate,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    Verdict,
)


@dataclass(frozen=True)
class _VerdictEvidence:
    """Typed verdict-driving observations for one candidate."""

    candidate: NormalizedCandidate
    coverage_observation: EvidenceObservation
    coverage: MappingCoverageSummary
    reciprocal_best_observation: EvidenceObservation | None
    reciprocal_best: ReciprocalBestMembershipSummary | None

    @property
    def is_full_mapping(self) -> bool:
        return self.coverage.status is MappingCoverageStatus.FULL

    @property
    def has_reciprocal_best_support(self) -> bool:
        return self.reciprocal_best is not None and self.reciprocal_best.status in {
            ReciprocalBestMembershipStatus.FULL,
            ReciprocalBestMembershipStatus.PARTIAL,
        }

    @property
    def is_fully_retained(self) -> bool:
        return (
            self.is_full_mapping
            and self.reciprocal_best is not None
            and self.reciprocal_best.status is ReciprocalBestMembershipStatus.FULL
        )

    @property
    def is_material_comparative_candidate(self) -> bool:
        # A full source-locus mapping is material even if reciprocal-best netting
        # rejects it: that disagreement is itself assessment-relevant. A partial
        # mapping becomes material only when at least some of its mapped geometry
        # survives the complete reciprocal-best resource.
        return self.is_full_mapping or self.has_reciprocal_best_support


def assess_candidates(
    source_interval: GenomicInterval,
    candidates: tuple[NormalizedCandidate, ...],
    *,
    evidence_tier: EvidenceAvailabilityTier,
) -> Assessment:
    """Return the deterministic v1 assessment for ``candidates``.

    Decision policy
    ---------------
    ``LIFTOVER_ONLY``
        One and only one full source-locus mapping is ``WELL_SUPPORTED``. Multiple
        chain-derived candidates are ``CONTESTED`` because the sparse tier has no
        comparative evidence with which to dismiss an alternative. No candidate, or
        one partial candidate, is ``INDETERMINATE``.

    ``COMPARATIVE``
        A candidate is *fully retained* only when source-locus coverage and
        reciprocal-best membership are both ``FULL``. A candidate remains materially
        in play when it either maps the full source locus or has non-``NONE``
        reciprocal-best membership. Exactly one fully retained candidate with no
        other material candidate is ``WELL_SUPPORTED``. Two or more material
        candidates are ``CONTESTED``. A sole full mapping with reciprocal-best
        membership ``NONE`` is also ``CONTESTED`` because the available mapping and
        self-consistency evidence materially disagree. A sole full mapping with
        ``PARTIAL`` reciprocal-best membership is ``INDETERMINATE``: the categorical
        state establishes mixed evidence but, without a magnitude threshold, does not
        establish that the disagreement is material. Remaining cases are likewise
        ``INDETERMINATE``.

    These are categorical rules, not calibrated biological truth criteria. A
    ``WELL_SUPPORTED`` verdict still does not mean the preferred candidate is correct.
    """

    if source_interval.length <= 0:
        raise ValueError("assessment requires a non-empty source interval")

    _validate_candidate_ids(candidates)
    _validate_distinct_candidate_geometries(candidates)
    if not candidates:
        return Assessment(
            source_interval=source_interval,
            verdict=Verdict.INDETERMINATE,
            evidence_tier=evidence_tier,
            decision_reason=AssessmentDecisionReason.NO_CANDIDATES,
            candidates=(),
        )

    profiles = tuple(
        _verdict_evidence_for_candidate(
            source_interval,
            candidate,
            evidence_tier=evidence_tier,
        )
        for candidate in candidates
    )

    if evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return _assess_liftover_only(source_interval, candidates, profiles)
    return _assess_comparative(source_interval, candidates, profiles)


def _assess_liftover_only(
    source_interval: GenomicInterval,
    candidates: tuple[NormalizedCandidate, ...],
    profiles: tuple[_VerdictEvidence, ...],
) -> Assessment:
    if len(profiles) > 1:
        supporting, contradicting = _mapping_coverage_references(profiles)
        return Assessment(
            source_interval=source_interval,
            verdict=Verdict.CONTESTED,
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            decision_reason=AssessmentDecisionReason.LIFTOVER_MULTIPLE_CANDIDATES,
            candidates=candidates,
            supporting_evidence=supporting,
            contradicting_evidence=contradicting,
        )

    profile = profiles[0]
    coverage_reference = _reference(profile.candidate, profile.coverage_observation)
    if profile.is_full_mapping:
        return Assessment(
            source_interval=source_interval,
            verdict=Verdict.WELL_SUPPORTED,
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            decision_reason=AssessmentDecisionReason.LIFTOVER_SINGLE_FULL_MAPPING,
            candidates=candidates,
            preferred_candidate_id=profile.candidate.candidate_id,
            supporting_evidence=(coverage_reference,),
        )

    return Assessment(
        source_interval=source_interval,
        verdict=Verdict.INDETERMINATE,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        decision_reason=AssessmentDecisionReason.LIFTOVER_SINGLE_PARTIAL_MAPPING,
        candidates=candidates,
        contradicting_evidence=(coverage_reference,),
    )


def _assess_comparative(
    source_interval: GenomicInterval,
    candidates: tuple[NormalizedCandidate, ...],
    profiles: tuple[_VerdictEvidence, ...],
) -> Assessment:
    material_profiles = tuple(
        profile for profile in profiles if profile.is_material_comparative_candidate
    )
    fully_retained = tuple(profile for profile in profiles if profile.is_fully_retained)

    if len(material_profiles) >= 2:
        supporting, contradicting = _comparative_references(material_profiles)
        return Assessment(
            source_interval=source_interval,
            verdict=Verdict.CONTESTED,
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            decision_reason=(
                AssessmentDecisionReason.COMPARATIVE_MULTIPLE_MATERIAL_CANDIDATES
            ),
            candidates=candidates,
            supporting_evidence=supporting,
            contradicting_evidence=contradicting,
        )

    if len(fully_retained) == 1:
        # The >=2-material branch above guarantees that this single fully retained
        # candidate is the only material candidate.
        profile = fully_retained[0]
        supporting, _ = _comparative_references((profile,))
        return Assessment(
            source_interval=source_interval,
            verdict=Verdict.WELL_SUPPORTED,
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            decision_reason=(
                AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_FULL
            ),
            candidates=candidates,
            preferred_candidate_id=profile.candidate.candidate_id,
            supporting_evidence=supporting,
        )

    full_profiles = tuple(profile for profile in profiles if profile.is_full_mapping)
    if full_profiles:
        # At most one material candidate can remain here. NONE is an exhaustive
        # categorical rejection of the candidate's mapped geometry under the stated
        # reciprocal-best completeness basis, so it is a material contradiction.
        # PARTIAL is deliberately weaker: without a magnitude threshold, v1 cannot
        # call the mixed state a material disagreement rather than insufficient
        # discrimination.
        profile = full_profiles[0]
        supporting, contradicting = _comparative_references((profile,))
        if profile.reciprocal_best is None:
            raise AssertionError(
                "comparative assessment requires reciprocal-best evidence"
            )
        if profile.reciprocal_best.status is ReciprocalBestMembershipStatus.NONE:
            verdict = Verdict.CONTESTED
            decision_reason = (
                AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_NONE
            )
        elif profile.reciprocal_best.status is ReciprocalBestMembershipStatus.PARTIAL:
            verdict = Verdict.INDETERMINATE
            decision_reason = (
                AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_FULL_RBEST_PARTIAL
            )
        else:
            raise AssertionError(
                "fully retained comparative candidate must have returned earlier"
            )
        return Assessment(
            source_interval=source_interval,
            verdict=verdict,
            evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
            decision_reason=decision_reason,
            candidates=candidates,
            supporting_evidence=supporting,
            contradicting_evidence=contradicting,
        )

    # The >=2-material branch already returned, and all remaining candidates are
    # partial because the full-candidate branch also returned. The existing
    # materiality predicate therefore leaves exactly two possible terminal states:
    # one material partial candidate, or no material candidate at all.
    if len(material_profiles) == 1:
        decision_reason = AssessmentDecisionReason.COMPARATIVE_SOLE_MATERIAL_PARTIAL
    else:
        if material_profiles:
            raise AssertionError(
                "comparative fallback cannot contain multiple material candidates"
            )
        decision_reason = AssessmentDecisionReason.COMPARATIVE_NO_MATERIAL_CANDIDATE

    supporting, contradicting = _comparative_references(material_profiles)
    return Assessment(
        source_interval=source_interval,
        verdict=Verdict.INDETERMINATE,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        decision_reason=decision_reason,
        candidates=candidates,
        supporting_evidence=supporting,
        contradicting_evidence=contradicting,
    )


def _verdict_evidence_for_candidate(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
    *,
    evidence_tier: EvidenceAvailabilityTier,
) -> _VerdictEvidence:
    _validate_candidate_geometry(source_interval, candidate)

    coverage_observation = _single_observation(
        candidate,
        EvidenceKind.MAPPING_COVERAGE,
        required=True,
    )
    if coverage_observation is None or not isinstance(
        coverage_observation.value, MappingCoverageSummary
    ):
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage must use "
            "MappingCoverageSummary"
        )
    coverage = coverage_observation.value
    expected_covered_bases = sum(
        segment.source_interval.length for segment in candidate.segments
    )
    if coverage.source_bases != source_interval.length:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage source_bases "
            "does not match the assessed source interval"
        )
    if coverage.covered_source_bases != expected_covered_bases:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage does not match "
            "its normalized mapping segments"
        )
    expected_uncovered = _uncovered_source_intervals(source_interval, candidate)
    if coverage.uncovered_source_intervals != expected_uncovered:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} mapping coverage uncovered "
            "intervals "
            "do not match its normalized mapping segments"
        )

    reciprocal_observation = _single_observation(
        candidate,
        EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
        required=evidence_tier is EvidenceAvailabilityTier.COMPARATIVE,
    )
    reciprocal: ReciprocalBestMembershipSummary | None = None
    if reciprocal_observation is not None:
        if not isinstance(
            reciprocal_observation.value, ReciprocalBestMembershipSummary
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} reciprocal-best evidence must "
                "use ReciprocalBestMembershipSummary"
            )
        reciprocal = reciprocal_observation.value
        if reciprocal.candidate_source_bases != expected_covered_bases:
            raise ValueError(
                f"candidate {candidate.candidate_id!r} reciprocal-best denominator "
                "does not match its normalized mapping segments"
            )
        if any(
            not _source_interval_is_covered_by_candidate(interval, candidate)
            for interval in reciprocal.covered_source_intervals
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} reciprocal-best covered "
                "intervals "
                "must lie within normalized mapping segments"
            )

    if (
        evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
        and reciprocal_observation is not None
    ):
        raise ValueError(
            "LIFTOVER_ONLY assessment cannot contain reciprocal-best evidence"
        )

    return _VerdictEvidence(
        candidate=candidate,
        coverage_observation=coverage_observation,
        coverage=coverage,
        reciprocal_best_observation=reciprocal_observation,
        reciprocal_best=reciprocal,
    )


def _single_observation(
    candidate: NormalizedCandidate,
    kind: EvidenceKind,
    *,
    required: bool,
) -> EvidenceObservation | None:
    matches = tuple(
        observation for observation in candidate.evidence if observation.kind is kind
    )
    if len(matches) > 1:
        raise ValueError(
            f"candidate {candidate.candidate_id!r} contains multiple {kind.value} "
            "observations; verdict-driving evidence must be unambiguous"
        )
    if not matches:
        if required:
            raise ValueError(
                f"candidate {candidate.candidate_id!r} is missing required "
                f"{kind.value} "
                "evidence"
            )
        return None
    return matches[0]


def _validate_candidate_ids(candidates: tuple[NormalizedCandidate, ...]) -> None:
    candidate_ids = tuple(candidate.candidate_id for candidate in candidates)
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError("candidate IDs must be unique before assessment")


def _validate_distinct_candidate_geometries(
    candidates: tuple[NormalizedCandidate, ...],
) -> None:
    """Reject multiple record IDs for one canonical local mapping hypothesis.

    Candidate IDs preserve source-record identity, but source-record identity is
    not itself evidence of a distinct placement. Different records can project the
    assessed locus to exactly the same source-to-target mapping. Counting those
    records separately would manufacture candidate multiplicity and can force a
    false ``CONTESTED`` verdict.

    v1 refuses to choose or silently merge such records because doing so would
    require explicit provenance/evidence aggregation semantics. Callers must
    consolidate equivalent mapping records before assessment.
    """

    candidate_id_by_geometry: dict[tuple[object, ...], str] = {}
    for candidate in candidates:
        geometry = _canonical_mapping_geometry(candidate)
        existing_id = candidate_id_by_geometry.get(geometry)
        if existing_id is not None:
            raise ValueError(
                f"candidate IDs {existing_id!r} and {candidate.candidate_id!r} "
                "describe the same canonical local mapping geometry; consolidate "
                "equivalent mapping records before assessment"
            )
        candidate_id_by_geometry[geometry] = candidate.candidate_id


def _canonical_mapping_geometry(candidate: NormalizedCandidate) -> tuple[object, ...]:
    """Return hypothesis-level geometry with adjacent collinear segments merged."""

    first_source = candidate.segments[0].source_interval
    target = candidate.target_interval
    merged: list[list[int]] = []

    for segment in candidate.segments:
        source = segment.source_interval
        segment_target = segment.target_interval
        if merged and _segments_are_collinear_adjacent(
            merged[-1],
            source.start,
            segment_target.start,
            segment_target.end,
            candidate.orientation,
        ):
            previous = merged[-1]
            previous[1] = source.end
            if candidate.orientation is MappingOrientation.SAME:
                previous[3] = segment_target.end
            else:
                previous[2] = segment_target.start
            continue

        merged.append(
            [
                source.start,
                source.end,
                segment_target.start,
                segment_target.end,
            ]
        )

    return (
        candidate.orientation,
        first_source.assembly,
        first_source.sequence_name,
        target.assembly,
        target.sequence_name,
        tuple(tuple(segment) for segment in merged),
    )


def _segments_are_collinear_adjacent(
    previous: list[int],
    source_start: int,
    target_start: int,
    target_end: int,
    orientation: MappingOrientation,
) -> bool:
    """Return whether two source-ordered segments form one continuous map."""

    if previous[1] != source_start:
        return False
    if orientation is MappingOrientation.SAME:
        return previous[3] == target_start
    return target_end == previous[2]


def _validate_candidate_geometry(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
) -> None:
    for segment in candidate.segments:
        if (
            segment.source_interval.assembly != source_interval.assembly
            or segment.source_interval.sequence_name != source_interval.sequence_name
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} source geometry does not match "
                "the assessed source sequence"
            )
        if (
            segment.source_interval.start < source_interval.start
            or segment.source_interval.end > source_interval.end
        ):
            raise ValueError(
                f"candidate {candidate.candidate_id!r} source geometry lies outside "
                "the assessed source interval"
            )


def _uncovered_source_intervals(
    source_interval: GenomicInterval,
    candidate: NormalizedCandidate,
) -> tuple[GenomicInterval, ...]:
    uncovered: list[GenomicInterval] = []
    cursor = source_interval.start

    for segment in candidate.segments:
        if cursor < segment.source_interval.start:
            uncovered.append(
                GenomicInterval(
                    assembly=source_interval.assembly,
                    sequence_name=source_interval.sequence_name,
                    start=cursor,
                    end=segment.source_interval.start,
                )
            )
        cursor = segment.source_interval.end

    if cursor < source_interval.end:
        uncovered.append(
            GenomicInterval(
                assembly=source_interval.assembly,
                sequence_name=source_interval.sequence_name,
                start=cursor,
                end=source_interval.end,
            )
        )
    return tuple(uncovered)


def _source_interval_is_covered_by_candidate(
    interval: GenomicInterval,
    candidate: NormalizedCandidate,
) -> bool:
    if interval.length <= 0:
        return False

    first_source = candidate.segments[0].source_interval
    if (
        interval.assembly != first_source.assembly
        or interval.sequence_name != first_source.sequence_name
    ):
        return False

    cursor = interval.start
    for segment in candidate.segments:
        source = segment.source_interval
        if source.end <= cursor:
            continue
        if source.start > cursor:
            return False
        cursor = max(cursor, source.end)
        if cursor >= interval.end:
            return True
    return False


def _mapping_coverage_references(
    profiles: tuple[_VerdictEvidence, ...],
) -> tuple[tuple[EvidenceReference, ...], tuple[EvidenceReference, ...]]:
    """Classify mapping-coverage observations by their categorical direction."""

    supporting: list[EvidenceReference] = []
    contradicting: list[EvidenceReference] = []
    for profile in profiles:
        reference = _reference(profile.candidate, profile.coverage_observation)
        if profile.is_full_mapping:
            supporting.append(reference)
        else:
            contradicting.append(reference)
    return tuple(supporting), tuple(contradicting)


def _comparative_references(
    profiles: tuple[_VerdictEvidence, ...],
) -> tuple[tuple[EvidenceReference, ...], tuple[EvidenceReference, ...]]:
    supporting: list[EvidenceReference] = []
    contradicting: list[EvidenceReference] = []

    for profile in profiles:
        coverage_reference = _reference(profile.candidate, profile.coverage_observation)
        if profile.is_full_mapping:
            supporting.append(coverage_reference)
        else:
            contradicting.append(coverage_reference)

        if (
            profile.reciprocal_best_observation is None
            or profile.reciprocal_best is None
        ):
            continue
        reciprocal_reference = _reference(
            profile.candidate, profile.reciprocal_best_observation
        )
        if profile.reciprocal_best.status is ReciprocalBestMembershipStatus.FULL:
            supporting.append(reciprocal_reference)
        elif profile.reciprocal_best.status is ReciprocalBestMembershipStatus.PARTIAL:
            # PARTIAL is mixed evidence: some geometry survives, but not all of it.
            # Preserve both roles rather than pretending it is wholly favorable or
            # wholly unfavorable.
            supporting.append(reciprocal_reference)
            contradicting.append(reciprocal_reference)
        else:
            contradicting.append(reciprocal_reference)

    return tuple(supporting), tuple(contradicting)


def _reference(
    candidate: NormalizedCandidate,
    observation: EvidenceObservation,
) -> EvidenceReference:
    return EvidenceReference(
        candidate_id=candidate.candidate_id,
        observation_id=observation.observation_id,
    )
