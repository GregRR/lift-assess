"""UCSC reciprocal-best membership evidence for normalized candidates.

UCSC's standard reciprocal-best resources are derived from the original chained
alignment by swapping direction, netting, subsetting the reciprocal-best net,
and swapping the retained chains back. They are therefore a self-consistency
check on the same upstream alignment, not an independent alignment experiment.

Implementation basis: UCSC kent ``src/hg/utils/automation/doRecipBest.pl``
(verified 2026-08-09). The reciprocal-best chain file can retain only part of an
original chain. This module therefore measures exact candidate-aligned source
coverage as FULL, PARTIAL, or NONE instead of reducing membership to a boolean.

Matching only accumulates confirming geometry, so an incomplete resource can
under-report membership (FULL -> PARTIAL/NONE or PARTIAL -> NONE) but cannot
create a false FULL. v1 nevertheless requires a complete resource or complete
candidate-relevant subset for every emitted membership state. Keeping one
conservative contract prevents callers from treating results from arbitrary
partial scans as comparable to exhaustive evidence.
"""

from collections.abc import Collection, Iterable
from dataclasses import replace

from .chain import ChainRecord
from .models import (
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
)
from .provenance import shares_upstream_source


def _annotate_candidate_with_reciprocal_best_chains(
    candidate: NormalizedCandidate,
    *,
    reciprocal_best_chains: Collection[ChainRecord],
    resource_completeness: ReciprocalBestResourceCompleteness,
    reciprocal_best_provenance: ProvenanceSource,
) -> NormalizedCandidate:
    """Attach locus-specific reciprocal-best membership to ``candidate``.

    ``reciprocal_best_chains`` must be a reusable, materialized collection. A bare
    generator is rejected because reusing a partially consumed iterator could turn
    real reciprocal-best matches into false ``PARTIAL`` or ``NONE`` observations.

    ``resource_completeness`` is an explicit caller claim that the collection is
    either the complete reciprocal-best resource for this assembly direction or a
    complete candidate-relevant subset derived from it. The claim is preserved in
    the evidence output; this function cannot independently prove that an external
    resource was completely loaded.

    Matching is based on exact source-to-target geometry, not chain ID. UCSC's
    reciprocal-best pipeline subsets and stitches chains, so chain identity is
    less informative here than whether the candidate's actual mapped bases are
    present in the reciprocal-best resource.
    """

    # Absence is evidence here, so a one-shot iterable is unsafe: exhaustion can
    # otherwise manufacture false PARTIAL/NONE results on a later candidate.
    if not isinstance(reciprocal_best_chains, Collection):
        raise TypeError(
            "reciprocal_best_chains must be a reusable materialized collection"
        )

    if not shares_upstream_source(
        candidate.mapping_provenance, reciprocal_best_provenance
    ):
        raise ValueError(
            "reciprocal-best provenance must share an upstream source with "
            "candidate mapping provenance"
        )

    covered_ranges: list[tuple[int, int]] = []
    source_sequence = candidate.segments[0].source_interval.sequence_name
    target_sequence = candidate.target_interval.sequence_name
    candidate_pair_chains_examined = 0

    for chain in reciprocal_best_chains:
        if chain.target_name != source_sequence:
            continue
        if chain.query_name != target_sequence:
            continue
        if chain.orientation is not candidate.orientation:
            continue

        # This is the audit count stored as ``chains_examined`` below. Count only
        # chains for the candidate's source/target sequence pair and orientation;
        # whole-genome chains for unrelated chromosomes would otherwise inflate the
        # number without saying anything about this candidate. It still does not
        # prove resource completeness or imply evidence strength.
        candidate_pair_chains_examined += 1

        for block in chain.iter_aligned_blocks():
            for segment in candidate.segments:
                overlap_start = max(segment.source_interval.start, block.target_start)
                overlap_end = min(segment.source_interval.end, block.target_end)
                if overlap_start >= overlap_end:
                    continue

                reciprocal_target = block.query_interval_for_target_interval(
                    overlap_start, overlap_end
                )
                candidate_target = _candidate_target_interval_for_source_overlap(
                    segment,
                    candidate.orientation,
                    overlap_start,
                    overlap_end,
                )
                if reciprocal_target == candidate_target:
                    covered_ranges.append((overlap_start, overlap_end))

    merged_ranges = _merge_ranges(covered_ranges)
    source_assembly = candidate.segments[0].source_interval.assembly
    covered_intervals = tuple(
        GenomicInterval(
            assembly=source_assembly,
            sequence_name=source_sequence,
            start=start,
            end=end,
        )
        for start, end in merged_ranges
    )
    covered_source_bases = sum(interval.length for interval in covered_intervals)
    candidate_source_bases = sum(
        segment.source_interval.length for segment in candidate.segments
    )

    if covered_source_bases == 0:
        status = ReciprocalBestMembershipStatus.NONE
    elif covered_source_bases == candidate_source_bases:
        status = ReciprocalBestMembershipStatus.FULL
    else:
        status = ReciprocalBestMembershipStatus.PARTIAL

    observation = EvidenceObservation(
        observation_id=f"{candidate.candidate_id}:reciprocal-best-membership",
        kind=EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
        value=ReciprocalBestMembershipSummary(
            status=status,
            resource_completeness=resource_completeness,
            chains_examined=candidate_pair_chains_examined,
            covered_source_bases=covered_source_bases,
            candidate_source_bases=candidate_source_bases,
            covered_source_intervals=covered_intervals,
        ),
        provenance=reciprocal_best_provenance,
    )
    return replace(candidate, evidence=(*candidate.evidence, observation))


def _candidate_target_interval_for_source_overlap(
    segment: MappingSegment,
    orientation: MappingOrientation,
    source_start: int,
    source_end: int,
) -> tuple[int, int]:
    """Map a contained source subinterval through one normalized mapping segment."""

    if (
        source_start < segment.source_interval.start
        or source_end > segment.source_interval.end
        or source_end <= source_start
    ):
        raise ValueError("source overlap lies outside mapping segment")

    offset_start = source_start - segment.source_interval.start
    offset_end = source_end - segment.source_interval.start
    if orientation is MappingOrientation.SAME:
        return (
            segment.target_interval.start + offset_start,
            segment.target_interval.start + offset_end,
        )
    return (
        segment.target_interval.end - offset_end,
        segment.target_interval.end - offset_start,
    )


def _merge_ranges(ranges: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Return the union of half-open source ranges in coordinate order."""

    ordered = sorted(ranges)
    if not ordered:
        return ()

    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        previous = merged[-1]
        if start <= previous[1]:
            previous[1] = max(previous[1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)
