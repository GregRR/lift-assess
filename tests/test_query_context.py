from __future__ import annotations

import pytest

from liftassess import (
    AssemblyIdentifier,
    GenomicInterval,
    PointQueryContextResult,
    QueryContextNotRunReason,
    QueryContextState,
    build_centered_point_context_interval,
    point_context_not_run,
)

ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")


def test_centered_point_context_uses_requested_odd_window() -> None:
    point = GenomicInterval(ASSEMBLY, "chr1", 100, 101)

    context = build_centered_point_context_interval(
        point,
        requested_window_bases=101,
        source_sequence_query_bound=1000,
    )

    assert context == GenomicInterval(ASSEMBLY, "chr1", 50, 151)


def test_centered_point_context_clips_at_source_boundary() -> None:
    point = GenomicInterval(ASSEMBLY, "chr1", 0, 1)

    context = build_centered_point_context_interval(
        point,
        requested_window_bases=101,
        source_sequence_query_bound=1000,
    )

    assert context == GenomicInterval(ASSEMBLY, "chr1", 0, 51)


def test_centered_point_context_clips_at_upper_source_boundary() -> None:
    point = GenomicInterval(ASSEMBLY, "chr1", 999, 1000)

    context = build_centered_point_context_interval(
        point,
        requested_window_bases=101,
        source_sequence_query_bound=1000,
    )

    assert context == GenomicInterval(ASSEMBLY, "chr1", 949, 1000)


def test_context_window_must_be_odd_and_nontrivial() -> None:
    point = GenomicInterval(ASSEMBLY, "chr1", 100, 101)

    with pytest.raises(ValueError, match="odd"):
        build_centered_point_context_interval(
            point,
            requested_window_bases=100,
            source_sequence_query_bound=1000,
        )
    with pytest.raises(ValueError, match="at least 3"):
        point_context_not_run(
            requested_window_bases=1,
            reason=QueryContextNotRunReason.INDEX_UNAVAILABLE,
        )


def test_unperformed_context_cannot_claim_results() -> None:
    with pytest.raises(ValueError, match="cannot carry assessment results"):
        PointQueryContextResult(
            check_state=QueryContextState.NOT_RUN,
            requested_window_bases=101,
            tested_source_interval=GenomicInterval(ASSEMBLY, "chr1", 50, 151),
            not_run_reason=QueryContextNotRunReason.INDEX_UNAVAILABLE,
        )
