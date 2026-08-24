from __future__ import annotations

import pytest

from liftassess import (
    AssemblyIdentifier,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
    ReverseCheckState,
    ReverseOriginalSourceCoverageState,
    ReverseRelationshipState,
    build_candidate_reverse_mapping_result,
    reverse_mapping_not_run,
    reverse_mapping_unavailable,
)

SOURCE_ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")
TARGET_ASSEMBLY = AssemblyIdentifier("targetAsm", "test")
CHAIN = ProvenanceSource("chain", "test chain")
REVERSE_CHAIN = ProvenanceSource("reverse-chain", "test reverse chain")


def _candidate(
    candidate_id: str,
    *,
    source_assembly: AssemblyIdentifier,
    source_sequence: str,
    source_spans: tuple[tuple[int, int], ...],
    target_assembly: AssemblyIdentifier,
    target_sequence: str,
    target_spans: tuple[tuple[int, int], ...],
    provenance: ProvenanceSource,
) -> NormalizedCandidate:
    segments = tuple(
        MappingSegment(
            GenomicInterval(source_assembly, source_sequence, source_start, source_end),
            GenomicInterval(target_assembly, target_sequence, target_start, target_end),
        )
        for (source_start, source_end), (target_start, target_end) in zip(
            source_spans,
            target_spans,
            strict=True,
        )
    )
    return NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=GenomicInterval(
            target_assembly,
            target_sequence,
            min(start for start, _ in target_spans),
            max(end for _, end in target_spans),
        ),
        orientation=MappingOrientation.SAME,
        mapping_provenance=provenance,
        segments=segments,
    )


def _forward_single() -> NormalizedCandidate:
    return _candidate(
        "forward",
        source_assembly=SOURCE_ASSEMBLY,
        source_sequence="chr1",
        source_spans=((100, 200),),
        target_assembly=TARGET_ASSEMBLY,
        target_sequence="chrA",
        target_spans=((1000, 1100),),
        provenance=CHAIN,
    )


def _reverse(
    candidate_id: str,
    *,
    source_spans: tuple[tuple[int, int], ...] = ((1000, 1100),),
    target_sequence: str = "chr1",
    target_spans: tuple[tuple[int, int], ...] = ((100, 200),),
) -> NormalizedCandidate:
    return _candidate(
        candidate_id,
        source_assembly=TARGET_ASSEMBLY,
        source_sequence="chrA",
        source_spans=source_spans,
        target_assembly=SOURCE_ASSEMBLY,
        target_sequence=target_sequence,
        target_spans=target_spans,
        provenance=REVERSE_CHAIN,
    )


def test_not_run_and_unavailable_keep_intended_exact_reverse_geometry() -> None:
    forward = _forward_single()

    not_run = reverse_mapping_not_run(forward)
    unavailable = reverse_mapping_unavailable(forward)

    assert not_run.check_state is ReverseCheckState.NOT_RUN
    assert unavailable.check_state is ReverseCheckState.UNAVAILABLE
    assert not_run.relationship is None
    assert unavailable.relationship is None
    assert not_run.queried_target_segments == (forward.segments[0].target_interval,)
    assert not_run.original_source_segments == (forward.segments[0].source_interval,)


def test_exact_reverse_return_is_complete_original_only_geometry() -> None:
    result = build_candidate_reverse_mapping_result(
        _forward_single(),
        ((_reverse("reverse-exact"),),),
    )

    assert result.check_state is ReverseCheckState.RUN
    assert result.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY
    assert result.reverse_projection_count == 1
    assert result.segments_with_reverse_projection == 1
    assert result.original_source_covered_bases == 100
    assert (
        result.original_source_coverage is ReverseOriginalSourceCoverageState.COMPLETE
    )
    assert result.exact_original_geometry_return


def test_reverse_return_elsewhere_is_not_promoted_to_original_source() -> None:
    result = build_candidate_reverse_mapping_result(
        _forward_single(),
        (
            (
                _reverse(
                    "reverse-other",
                    target_sequence="chr9",
                    target_spans=((500, 600),),
                ),
            ),
        ),
    )

    assert result.relationship is ReverseRelationshipState.ELSEWHERE_ONLY
    assert result.original_source_covered_bases == 0
    assert result.original_source_coverage is ReverseOriginalSourceCoverageState.NONE
    assert not result.exact_original_geometry_return


def test_completed_reverse_run_can_have_no_projection() -> None:
    result = build_candidate_reverse_mapping_result(_forward_single(), ((),))

    assert result.relationship is ReverseRelationshipState.NO_PROJECTION
    assert result.reverse_projection_count == 0
    assert result.segments_with_reverse_projection == 0
    assert result.original_source_coverage is ReverseOriginalSourceCoverageState.NONE


def test_reverse_results_can_include_original_and_elsewhere_returns() -> None:
    result = build_candidate_reverse_mapping_result(
        _forward_single(),
        (
            (
                _reverse("reverse-exact"),
                _reverse(
                    "reverse-other",
                    target_sequence="chr9",
                    target_spans=((500, 600),),
                ),
            ),
        ),
    )

    assert result.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE
    assert (
        result.original_source_coverage is ReverseOriginalSourceCoverageState.COMPLETE
    )
    assert result.exact_original_geometry_return
    assert result.reverse_projection_count == 2


def test_partial_original_return_stays_partial_without_elsewhere_return() -> None:
    result = build_candidate_reverse_mapping_result(
        _forward_single(),
        (
            (
                _reverse(
                    "reverse-partial",
                    source_spans=((1000, 1090),),
                    target_spans=((100, 190),),
                ),
            ),
        ),
    )

    assert result.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY
    assert result.original_source_covered_bases == 90
    assert result.original_source_coverage is ReverseOriginalSourceCoverageState.PARTIAL
    assert not result.exact_original_geometry_return


def test_fragmented_forward_candidate_reverses_segments_not_bounding_span() -> None:
    forward = _candidate(
        "forward-split",
        source_assembly=SOURCE_ASSEMBLY,
        source_sequence="chr1",
        source_spans=((100, 140), (160, 200)),
        target_assembly=TARGET_ASSEMBLY,
        target_sequence="chrA",
        target_spans=((1000, 1040), (1060, 1100)),
        provenance=CHAIN,
    )
    reverse_first = _candidate(
        "reverse-first",
        source_assembly=TARGET_ASSEMBLY,
        source_sequence="chrA",
        source_spans=((1000, 1040),),
        target_assembly=SOURCE_ASSEMBLY,
        target_sequence="chr1",
        target_spans=((100, 140),),
        provenance=REVERSE_CHAIN,
    )
    reverse_second = _candidate(
        "reverse-second",
        source_assembly=TARGET_ASSEMBLY,
        source_sequence="chrA",
        source_spans=((1060, 1100),),
        target_assembly=SOURCE_ASSEMBLY,
        target_sequence="chr1",
        target_spans=((160, 200),),
        provenance=REVERSE_CHAIN,
    )

    result = build_candidate_reverse_mapping_result(
        forward,
        ((reverse_first,), (reverse_second,)),
    )

    assert result.queried_target_segments == (
        GenomicInterval(TARGET_ASSEMBLY, "chrA", 1000, 1040),
        GenomicInterval(TARGET_ASSEMBLY, "chrA", 1060, 1100),
    )
    assert forward.target_interval == GenomicInterval(
        TARGET_ASSEMBLY,
        "chrA",
        1000,
        1100,
    )
    assert result.original_source_bases == 80
    assert result.original_source_covered_bases == 80
    assert (
        result.original_source_coverage is ReverseOriginalSourceCoverageState.COMPLETE
    )
    assert result.exact_original_geometry_return


def test_reverse_candidate_source_must_stay_inside_exact_segment_query() -> None:
    with pytest.raises(
        ValueError,
        match="reverse candidate source segments must lie within the exact queried",
    ):
        build_candidate_reverse_mapping_result(
            _forward_single(),
            (
                (
                    _reverse(
                        "reverse-outside",
                        source_spans=((999, 1100),),
                        target_spans=((100, 201),),
                    ),
                ),
            ),
        )


def test_reverse_candidates_must_represent_distinct_canonical_geometry() -> None:
    duplicate_geometry = _reverse(
        "reverse-partitioned",
        source_spans=((1000, 1050), (1050, 1100)),
        target_spans=((100, 150), (150, 200)),
    )

    with pytest.raises(ValueError, match="identical normalized mapping geometry"):
        build_candidate_reverse_mapping_result(
            _forward_single(),
            ((_reverse("reverse-single"), duplicate_geometry),),
        )
