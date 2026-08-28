from liftassess.batch import (
    BatchInputRecord,
    BatchRecordAssessment,
    BatchTargetRelationshipKind,
    build_batch_target_relationships,
)
from liftassess.models import (
    AssemblyIdentifier,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
)

SOURCE = AssemblyIdentifier(name="hg38", provider="UCSC")
TARGET = AssemblyIdentifier(name="hg19", provider="UCSC")
PROVENANCE = ProvenanceSource(source_id="chain", label="chain")


def _record(record_id: str, start: int, end: int) -> BatchInputRecord:
    return BatchInputRecord(
        record_id=record_id,
        source_interval=GenomicInterval(SOURCE, "chr1", start, end),
        source_line_number=int(record_id.removeprefix("row-")),
    )


def _candidate(
    candidate_id: str,
    *,
    source_start: int,
    target_segments: tuple[tuple[int, int], ...],
    orientation: MappingOrientation = MappingOrientation.SAME,
) -> NormalizedCandidate:
    source_cursor = source_start
    segments: list[MappingSegment] = []
    ordered_target_segments = (
        target_segments
        if orientation is MappingOrientation.SAME
        else tuple(reversed(target_segments))
    )
    for target_start, target_end in ordered_target_segments:
        length = target_end - target_start
        segments.append(
            MappingSegment(
                source_interval=GenomicInterval(
                    SOURCE,
                    "chr1",
                    source_cursor,
                    source_cursor + length,
                ),
                target_interval=GenomicInterval(
                    TARGET,
                    "chr2",
                    target_start,
                    target_end,
                ),
            )
        )
        source_cursor += length

    return NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=GenomicInterval(
            TARGET,
            "chr2",
            min(start for start, _ in target_segments),
            max(end for _, end in target_segments),
        ),
        orientation=orientation,
        mapping_provenance=PROVENANCE,
        segments=tuple(segments),
    )


def test_batch_relationships_detect_exact_target_collision_across_partitioning() -> (
    None
):
    left = BatchRecordAssessment(
        record=_record("row-1", 0, 10),
        candidates=(
            _candidate(
                "left",
                source_start=0,
                target_segments=((100, 105), (105, 110)),
            ),
        ),
    )
    right = BatchRecordAssessment(
        record=_record("row-2", 20, 30),
        candidates=(
            _candidate("right", source_start=20, target_segments=((100, 110),)),
        ),
    )

    result = build_batch_target_relationships((left, right))

    assert len(result.relationships) == 1
    relationship = result.relationships[0]
    assert relationship.kind is BatchTargetRelationshipKind.EXACT_TARGET_COLLISION
    assert relationship.overlap_intervals == (
        GenomicInterval(TARGET, "chr2", 100, 110),
    )


def test_batch_relationships_detect_overlapping_but_offset_projection() -> None:
    left = BatchRecordAssessment(
        record=_record("row-1", 0, 10),
        candidates=(_candidate("left", source_start=0, target_segments=((100, 110),)),),
    )
    right = BatchRecordAssessment(
        record=_record("row-2", 20, 30),
        candidates=(
            _candidate("right", source_start=20, target_segments=((105, 115),)),
        ),
    )

    result = build_batch_target_relationships((left, right))

    assert result.relationships[0].kind is (
        BatchTargetRelationshipKind.OVERLAPPING_TARGET_PROJECTIONS
    )
    assert result.relationships[0].overlap_intervals == (
        GenomicInterval(TARGET, "chr2", 105, 110),
    )


def test_batch_relationships_do_not_use_fragmented_bounding_span_as_coverage() -> None:
    left = BatchRecordAssessment(
        record=_record("row-1", 0, 10),
        candidates=(
            _candidate(
                "fragmented",
                source_start=0,
                target_segments=((100, 105), (200, 205)),
            ),
        ),
    )
    right = BatchRecordAssessment(
        record=_record("row-2", 20, 30),
        candidates=(
            _candidate("inside-gap", source_start=20, target_segments=((150, 160),)),
        ),
    )

    result = build_batch_target_relationships((left, right))

    assert result.relationships == ()


def test_batch_relationships_can_overlap_one_fragment_without_becoming_exact() -> None:
    left = BatchRecordAssessment(
        record=_record("row-1", 0, 10),
        candidates=(
            _candidate(
                "fragmented",
                source_start=0,
                target_segments=((100, 105), (200, 205)),
            ),
        ),
    )
    right = BatchRecordAssessment(
        record=_record("row-2", 20, 25),
        candidates=(
            _candidate(
                "partial-overlap",
                source_start=20,
                target_segments=((102, 107),),
            ),
        ),
    )

    result = build_batch_target_relationships((left, right))

    assert result.relationships[0].kind is (
        BatchTargetRelationshipKind.OVERLAPPING_TARGET_PROJECTIONS
    )
    assert result.relationships[0].overlap_intervals == (
        GenomicInterval(TARGET, "chr2", 102, 105),
    )


def test_batch_relationships_ignore_different_target_sequences() -> None:
    left_candidate = _candidate("left", source_start=0, target_segments=((100, 110),))
    right_candidate = NormalizedCandidate(
        candidate_id="right",
        target_interval=GenomicInterval(TARGET, "chr3", 100, 110),
        orientation=MappingOrientation.SAME,
        mapping_provenance=PROVENANCE,
        segments=(
            MappingSegment(
                source_interval=GenomicInterval(SOURCE, "chr1", 20, 30),
                target_interval=GenomicInterval(TARGET, "chr3", 100, 110),
            ),
        ),
    )

    result = build_batch_target_relationships(
        (
            BatchRecordAssessment(_record("row-1", 0, 10), (left_candidate,)),
            BatchRecordAssessment(_record("row-2", 20, 30), (right_candidate,)),
        )
    )

    assert result.relationships == ()


def test_batch_relationships_ignore_relationships_within_one_record() -> None:
    assessment = BatchRecordAssessment(
        record=_record("row-1", 0, 10),
        candidates=(
            _candidate("a", source_start=0, target_segments=((100, 110),)),
            _candidate("b", source_start=0, target_segments=((105, 115),)),
        ),
    )

    result = build_batch_target_relationships((assessment,))

    assert result.relationships == ()


def test_batch_record_assessment_rejects_candidate_from_outside_record() -> None:
    record = _record("row-1", 0, 5)
    candidate = _candidate("outside", source_start=0, target_segments=((100, 110),))

    try:
        BatchRecordAssessment(record=record, candidates=(candidate,))
    except ValueError as exc:
        assert "contained within" in str(exc)
    else:
        raise AssertionError("expected out-of-record candidate to be rejected")
