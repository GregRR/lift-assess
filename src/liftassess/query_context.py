"""Point-query local-context execution model.

Milestone 20 keeps neighborhood execution separate from interpretation.  This module
represents whether a point-context check ran, the exact source window tested, and the
chain-derived candidates produced for that window.  The factual comparison with the
original point belongs to :mod:`liftassess.result_profile`.
"""

from dataclasses import dataclass
from enum import Enum

from .models import GenomicInterval, NormalizedCandidate

DEFAULT_POINT_CONTEXT_BASES = 101


class QueryContextState(str, Enum):
    """Whether point-neighborhood chain geometry was actually assessed."""

    NOT_RUN = "NOT_RUN"
    RUN = "RUN"


class QueryContextNotRunReason(str, Enum):
    """Why automatic point-context execution did not run."""

    INDEX_UNAVAILABLE = "INDEX_UNAVAILABLE"
    INDEX_UNUSABLE = "INDEX_UNUSABLE"
    SOURCE_BOUNDS_UNAVAILABLE = "SOURCE_BOUNDS_UNAVAILABLE"


@dataclass(frozen=True)
class PointQueryContextResult:
    """Raw chain-geometry result for one requested point-context window."""

    check_state: QueryContextState
    requested_window_bases: int
    tested_source_interval: GenomicInterval | None = None
    candidates: tuple[NormalizedCandidate, ...] = ()
    not_run_reason: QueryContextNotRunReason | None = None

    def __post_init__(self) -> None:
        _validate_context_window_bases(self.requested_window_bases)
        if self.check_state is QueryContextState.RUN:
            if self.tested_source_interval is None:
                raise ValueError("completed point context requires a tested interval")
            if self.tested_source_interval.length <= 0:
                raise ValueError("point context tested interval must not be empty")
            if self.tested_source_interval.length > self.requested_window_bases:
                raise ValueError(
                    "point context tested interval cannot exceed the requested window"
                )
            if self.not_run_reason is not None:
                raise ValueError(
                    "completed point context cannot carry a not-run reason"
                )
            return

        if self.tested_source_interval is not None or self.candidates:
            raise ValueError(
                "unperformed point context cannot carry assessment results"
            )
        if self.not_run_reason is None:
            raise ValueError("unperformed point context requires a not-run reason")


def build_centered_point_context_interval(
    point_interval: GenomicInterval,
    *,
    requested_window_bases: int = DEFAULT_POINT_CONTEXT_BASES,
    source_sequence_query_bound: int,
) -> GenomicInterval:
    """Return a boundary-clipped odd-width context window around one source base."""

    if point_interval.length != 1:
        raise ValueError("point context requires a one-base source interval")
    _validate_context_window_bases(requested_window_bases)
    if source_sequence_query_bound <= 0:
        raise ValueError("source sequence query bound must be positive")
    if point_interval.end > source_sequence_query_bound:
        raise ValueError("point interval exceeds the source sequence query bound")

    flank = requested_window_bases // 2
    return GenomicInterval(
        assembly=point_interval.assembly,
        sequence_name=point_interval.sequence_name,
        start=max(0, point_interval.start - flank),
        end=min(source_sequence_query_bound, point_interval.end + flank),
    )


def point_context_not_run(
    *,
    requested_window_bases: int,
    reason: QueryContextNotRunReason,
) -> PointQueryContextResult:
    """Build a truthful unperformed point-context result."""

    return PointQueryContextResult(
        check_state=QueryContextState.NOT_RUN,
        requested_window_bases=requested_window_bases,
        not_run_reason=reason,
    )


def _validate_context_window_bases(window_bases: int) -> None:
    if window_bases < 3:
        raise ValueError("point context window must be at least 3 bases")
    if window_bases % 2 == 0:
        raise ValueError("point context window must contain an odd number of bases")
