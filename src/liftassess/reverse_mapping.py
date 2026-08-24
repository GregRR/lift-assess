"""Structured actual reverse-mapping context.

Reverse mapping is a separate analysis from UCSC reciprocal-best membership.  This
module compares explicit reverse candidate sets with the exact target segments of one
forward candidate.  It records factual geometry only; it does not assign an aggregate
verdict, confidence score, or biological interpretation.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .chain_index import ChainIndex
from .models import GenomicInterval, NormalizedCandidate, ProvenanceSource
from .resource_cache import CachedUCSCChainResource, CachedUCSCResourceBundle
from .resource_files import (
    ResourceReadProgressCallback,
    build_ucsc_chain_candidates_for_intervals_from_cached_chain,
)


class ReverseCheckState(str, Enum):
    """Whether actual reverse assessment was performed."""

    NOT_RUN = "NOT_RUN"
    UNAVAILABLE = "UNAVAILABLE"
    RUN = "RUN"


class ReverseRelationshipState(str, Enum):
    """Where observed reverse projections land relative to the original source."""

    NO_PROJECTION = "NO_PROJECTION"
    ORIGINAL_SOURCE_ONLY = "ORIGINAL_SOURCE_ONLY"
    ELSEWHERE_ONLY = "ELSEWHERE_ONLY"
    ORIGINAL_SOURCE_AND_ELSEWHERE = "ORIGINAL_SOURCE_AND_ELSEWHERE"


class ReverseOriginalSourceCoverageState(str, Enum):
    """How much of the forward candidate's aligned source geometry is returned."""

    NONE = "NONE"
    PARTIAL = "PARTIAL"
    COMPLETE = "COMPLETE"


@dataclass(frozen=True)
class ReverseSegmentResult:
    """Actual reverse candidates for one exact forward target segment."""

    queried_target_segment: GenomicInterval
    expected_original_source_segment: GenomicInterval
    candidates: tuple[NormalizedCandidate, ...]

    def __post_init__(self) -> None:
        if self.queried_target_segment.length <= 0:
            raise ValueError("reverse query segment must span at least one base")
        if (
            self.expected_original_source_segment.length
            != self.queried_target_segment.length
        ):
            raise ValueError(
                "reverse query and expected original source segments must have "
                "equal length"
            )
        _validate_unique_candidate_ids(self.candidates)
        for candidate in self.candidates:
            _validate_reverse_candidate_source_geometry(
                self.queried_target_segment,
                candidate,
            )


@dataclass(frozen=True)
class CandidateReverseMappingResult:
    """Orthogonal reverse-mapping facts for one forward candidate.

    ``original_source_segments`` and ``queried_target_segments`` are retained even when
    the check is not run or unavailable so the exact intended reverse geometry remains
    explicit.  A run contains one ``ReverseSegmentResult`` per forward mapping segment.
    """

    forward_candidate_id: str
    check_state: ReverseCheckState
    original_source_segments: tuple[GenomicInterval, ...]
    queried_target_segments: tuple[GenomicInterval, ...]
    segment_results: tuple[ReverseSegmentResult, ...] = ()

    def __post_init__(self) -> None:
        if not self.forward_candidate_id:
            raise ValueError("forward_candidate_id must not be empty")
        if not self.original_source_segments:
            raise ValueError(
                "reverse mapping requires at least one original source segment"
            )
        if len(self.original_source_segments) != len(self.queried_target_segments):
            raise ValueError(
                "original source and queried target segment counts must match"
            )
        for source_segment, target_segment in zip(
            self.original_source_segments,
            self.queried_target_segments,
            strict=True,
        ):
            if source_segment.length != target_segment.length:
                raise ValueError(
                    "paired original source and queried target segments must have "
                    "equal length"
                )

        if self.check_state is ReverseCheckState.RUN:
            if len(self.segment_results) != len(self.queried_target_segments):
                raise ValueError(
                    "a completed reverse run requires one result per queried target "
                    "segment"
                )
            for expected_source, queried_target, result in zip(
                self.original_source_segments,
                self.queried_target_segments,
                self.segment_results,
                strict=True,
            ):
                if result.expected_original_source_segment != expected_source:
                    raise ValueError(
                        "reverse segment result does not match original source geometry"
                    )
                if result.queried_target_segment != queried_target:
                    raise ValueError(
                        "reverse segment result does not match queried target geometry"
                    )
        elif self.segment_results:
            raise ValueError(
                "reverse segment results are only valid when reverse assessment was run"
            )

    @property
    def original_source_bases(self) -> int:
        """Number of forward-aligned source bases eligible for reverse comparison."""

        return sum(segment.length for segment in self.original_source_segments)

    @property
    def reverse_projection_count(self) -> int:
        """Total reverse candidates across all exact target-segment queries."""

        return sum(len(result.candidates) for result in self.segment_results)

    @property
    def segments_with_reverse_projection(self) -> int:
        """Number of exact forward target segments with at least one reverse result."""

        return sum(bool(result.candidates) for result in self.segment_results)

    @property
    def original_source_covered_bases(self) -> int:
        """Union of original aligned source bases touched by reverse projections."""

        if self.check_state is not ReverseCheckState.RUN:
            return 0
        return sum(_original_overlap_bases(result) for result in self.segment_results)

    @property
    def original_source_coverage(self) -> ReverseOriginalSourceCoverageState:
        """Coverage of original aligned source geometry by observed reverse returns."""

        covered = self.original_source_covered_bases
        if covered == 0:
            return ReverseOriginalSourceCoverageState.NONE
        if covered == self.original_source_bases:
            return ReverseOriginalSourceCoverageState.COMPLETE
        return ReverseOriginalSourceCoverageState.PARTIAL

    @property
    def exact_original_geometry_return(self) -> bool:
        """Whether every forward target segment reconstructs its source segment."""

        if self.check_state is not ReverseCheckState.RUN:
            return False
        return all(
            _segment_has_exact_original_return(result)
            for result in self.segment_results
        )

    @property
    def relationship(self) -> ReverseRelationshipState | None:
        """Categorical locus relationship for a completed reverse check.

        Exactness and amount of returned original-source geometry remain separate facts;
        this relationship only records whether returned bases touch the corresponding
        original source geometry and whether any returned bases land elsewhere.
        """

        if self.check_state is not ReverseCheckState.RUN:
            return None
        if self.reverse_projection_count == 0:
            return ReverseRelationshipState.NO_PROJECTION

        original_present = self.original_source_covered_bases > 0
        elsewhere_present = any(
            _candidate_has_elsewhere_geometry(result, candidate)
            for result in self.segment_results
            for candidate in result.candidates
        )
        if original_present and elsewhere_present:
            return ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE
        if original_present:
            return ReverseRelationshipState.ORIGINAL_SOURCE_ONLY
        return ReverseRelationshipState.ELSEWHERE_ONLY


def build_reverse_mapping_results_from_cached_chain(
    forward_candidates: tuple[NormalizedCandidate, ...],
    reverse_chain: CachedUCSCChainResource,
    *,
    reverse_alignment_provenance: ProvenanceSource,
    progress_callback: ResourceReadProgressCallback | None = None,
    chain_index: ChainIndex | None = None,
) -> tuple[CandidateReverseMappingResult, ...]:
    """Run actual chain-only reverse mapping for exact forward target segments."""

    if not forward_candidates:
        return ()

    original_source_assembly = (
        forward_candidates[0].segments[0].source_interval.assembly
    )
    queried_segments: list[GenomicInterval] = []
    segment_counts: list[int] = []
    for candidate in forward_candidates:
        if any(
            segment.source_interval.assembly != original_source_assembly
            for segment in candidate.segments
        ):
            raise ValueError(
                "forward candidates for one reverse run must share the original "
                "source assembly"
            )
        segment_counts.append(len(candidate.segments))
        queried_segments.extend(
            segment.target_interval for segment in candidate.segments
        )

    reverse_candidates = build_ucsc_chain_candidates_for_intervals_from_cached_chain(
        queried_segments,
        reverse_chain,
        target_assembly=original_source_assembly,
        alignment_provenance=reverse_alignment_provenance,
        progress_callback=progress_callback,
        chain_index=chain_index,
    )

    results: list[CandidateReverseMappingResult] = []
    offset = 0
    for candidate, segment_count in zip(
        forward_candidates, segment_counts, strict=True
    ):
        candidate_reverse_sets = reverse_candidates[offset : offset + segment_count]
        results.append(
            build_candidate_reverse_mapping_result(
                candidate,
                candidate_reverse_sets,
            )
        )
        offset += segment_count
    return tuple(results)


def build_reverse_mapping_results_from_cached_bundle(
    forward_candidates: tuple[NormalizedCandidate, ...],
    reverse_bundle: CachedUCSCResourceBundle,
    *,
    reverse_alignment_provenance: ProvenanceSource,
    progress_callback: ResourceReadProgressCallback | None = None,
    chain_index: ChainIndex | None = None,
) -> tuple[CandidateReverseMappingResult, ...]:
    """Compatibility wrapper that consumes only the reverse bundle's chain."""

    return build_reverse_mapping_results_from_cached_chain(
        forward_candidates,
        CachedUCSCChainResource(
            source_db=reverse_bundle.source_db,
            target_db=reverse_bundle.target_db,
            evidence_tier=reverse_bundle.evidence_tier,
            chain=reverse_bundle.chain,
        ),
        reverse_alignment_provenance=reverse_alignment_provenance,
        progress_callback=progress_callback,
        chain_index=chain_index,
    )


def reverse_mapping_not_run(
    forward_candidate: NormalizedCandidate,
) -> CandidateReverseMappingResult:
    """Represent an explicitly unperformed reverse check for one forward candidate."""

    return _empty_reverse_mapping(forward_candidate, ReverseCheckState.NOT_RUN)


def reverse_mapping_unavailable(
    forward_candidate: NormalizedCandidate,
) -> CandidateReverseMappingResult:
    """Represent unavailable reverse resources for one forward candidate."""

    return _empty_reverse_mapping(forward_candidate, ReverseCheckState.UNAVAILABLE)


def build_candidate_reverse_mapping_result(
    forward_candidate: NormalizedCandidate,
    reverse_candidates_by_segment: tuple[tuple[NormalizedCandidate, ...], ...],
) -> CandidateReverseMappingResult:
    """Compare actual reverse candidates with one forward candidate's exact geometry.

    ``reverse_candidates_by_segment`` must contain one candidate tuple for each exact
    forward mapping segment, in forward-segment order.  The reverse query for a
    fragmented candidate is therefore the set of exact target segments, never the
    candidate's target bounding span across unaligned gaps.
    """

    if len(reverse_candidates_by_segment) != len(forward_candidate.segments):
        raise ValueError(
            "reverse candidate sets must match the forward candidate segment count"
        )

    segment_results = tuple(
        ReverseSegmentResult(
            queried_target_segment=forward_segment.target_interval,
            expected_original_source_segment=forward_segment.source_interval,
            candidates=reverse_candidates,
        )
        for forward_segment, reverse_candidates in zip(
            forward_candidate.segments,
            reverse_candidates_by_segment,
            strict=True,
        )
    )
    return CandidateReverseMappingResult(
        forward_candidate_id=forward_candidate.candidate_id,
        check_state=ReverseCheckState.RUN,
        original_source_segments=tuple(
            segment.source_interval for segment in forward_candidate.segments
        ),
        queried_target_segments=tuple(
            segment.target_interval for segment in forward_candidate.segments
        ),
        segment_results=segment_results,
    )


def _empty_reverse_mapping(
    forward_candidate: NormalizedCandidate,
    check_state: ReverseCheckState,
) -> CandidateReverseMappingResult:
    if check_state is ReverseCheckState.RUN:
        raise ValueError("completed reverse runs require explicit segment results")
    return CandidateReverseMappingResult(
        forward_candidate_id=forward_candidate.candidate_id,
        check_state=check_state,
        original_source_segments=tuple(
            segment.source_interval for segment in forward_candidate.segments
        ),
        queried_target_segments=tuple(
            segment.target_interval for segment in forward_candidate.segments
        ),
    )


def _validate_unique_candidate_ids(candidates: tuple[NormalizedCandidate, ...]) -> None:
    candidate_ids = [candidate.candidate_id for candidate in candidates]
    if len(candidate_ids) != len(set(candidate_ids)):
        raise ValueError(
            "reverse candidate IDs must be unique within one segment query"
        )


def _validate_reverse_candidate_source_geometry(
    queried_target_segment: GenomicInterval,
    candidate: NormalizedCandidate,
) -> None:
    for segment in candidate.segments:
        source = segment.source_interval
        if (
            source.assembly != queried_target_segment.assembly
            or source.sequence_name != queried_target_segment.sequence_name
            or source.start < queried_target_segment.start
            or source.end > queried_target_segment.end
        ):
            raise ValueError(
                "reverse candidate source segments must lie within the exact queried "
                "forward target segment"
            )


def _original_overlap_bases(result: ReverseSegmentResult) -> int:
    expected = result.expected_original_source_segment
    overlaps: list[tuple[int, int]] = []
    for candidate in result.candidates:
        for segment in candidate.segments:
            target = segment.target_interval
            if (
                target.assembly != expected.assembly
                or target.sequence_name != expected.sequence_name
            ):
                continue
            start = max(target.start, expected.start)
            end = min(target.end, expected.end)
            if start < end:
                overlaps.append((start, end))
    return _union_span_length(overlaps)


def _candidate_has_elsewhere_geometry(
    result: ReverseSegmentResult,
    candidate: NormalizedCandidate,
) -> bool:
    expected = result.expected_original_source_segment
    return any(
        segment.target_interval.assembly != expected.assembly
        or segment.target_interval.sequence_name != expected.sequence_name
        or segment.target_interval.start < expected.start
        or segment.target_interval.end > expected.end
        for segment in candidate.segments
    )


def _segment_has_exact_original_return(result: ReverseSegmentResult) -> bool:
    return any(
        _candidate_fully_covers_query(result.queried_target_segment, candidate)
        and _candidate_exactly_covers_expected_source(
            result.expected_original_source_segment,
            candidate,
        )
        for candidate in result.candidates
    )


def _candidate_fully_covers_query(
    query: GenomicInterval,
    candidate: NormalizedCandidate,
) -> bool:
    spans = [
        (segment.source_interval.start, segment.source_interval.end)
        for segment in candidate.segments
    ]
    return _union_span_length(spans) == query.length


def _candidate_exactly_covers_expected_source(
    expected: GenomicInterval,
    candidate: NormalizedCandidate,
) -> bool:
    spans: list[tuple[int, int]] = []
    for segment in candidate.segments:
        target = segment.target_interval
        if (
            target.assembly != expected.assembly
            or target.sequence_name != expected.sequence_name
        ):
            return False
        if target.start < expected.start or target.end > expected.end:
            return False
        spans.append((target.start, target.end))
    return _union_span_length(spans) == expected.length


def _union_span_length(spans: list[tuple[int, int]]) -> int:
    if not spans:
        return 0
    ordered = sorted(spans)
    total = 0
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start > current_end:
            total += current_end - current_start
            current_start, current_end = start, end
        else:
            current_end = max(current_end, end)
    return total + current_end - current_start
