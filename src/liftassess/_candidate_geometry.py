"""Shared canonical geometry helpers for normalized mapping candidates."""

from __future__ import annotations

from .models import MappingOrientation, NormalizedCandidate


def canonical_mapping_segments(
    candidate: NormalizedCandidate,
) -> tuple[tuple[int, int, int, int], ...]:
    """Canonicalize collinear adjacent segment partitions for geometry comparison."""

    canonical_segments: list[list[int]] = []
    for segment in candidate.segments:
        current = [
            segment.source_interval.start,
            segment.source_interval.end,
            segment.target_interval.start,
            segment.target_interval.end,
        ]
        if canonical_segments and _segments_are_collinear_adjacent(
            canonical_segments[-1],
            current,
            orientation=candidate.orientation,
        ):
            previous = canonical_segments[-1]
            previous[1] = current[1]
            if candidate.orientation is MappingOrientation.SAME:
                previous[3] = current[3]
            else:
                previous[2] = current[2]
        else:
            canonical_segments.append(current)

    return tuple(
        (segment[0], segment[1], segment[2], segment[3])
        for segment in canonical_segments
    )


def canonical_mapping_geometry(candidate: NormalizedCandidate) -> tuple[object, ...]:
    """Return hypothesis-level local geometry independent of chain record identity."""

    first_source = candidate.segments[0].source_interval
    target = candidate.target_interval
    return (
        first_source.assembly,
        first_source.sequence_name,
        target.assembly,
        target.sequence_name,
        candidate.orientation,
        canonical_mapping_segments(candidate),
    )


def validate_distinct_candidate_geometries(
    candidates: tuple[NormalizedCandidate, ...],
) -> None:
    """Reject distinct IDs that represent the same canonical mapping hypothesis."""

    seen: dict[tuple[object, ...], str] = {}
    for candidate in candidates:
        key = canonical_mapping_geometry(candidate)
        existing = seen.get(key)
        if existing is not None:
            raise ValueError(
                "distinct candidate IDs must not describe identical normalized "
                f"mapping geometry: {existing!r} and {candidate.candidate_id!r}"
            )
        seen[key] = candidate.candidate_id


def _segments_are_collinear_adjacent(
    previous: list[int],
    current: list[int],
    *,
    orientation: MappingOrientation,
) -> bool:
    if previous[1] != current[0]:
        return False
    if orientation is MappingOrientation.SAME:
        return previous[3] == current[2]
    return previous[2] == current[3]
