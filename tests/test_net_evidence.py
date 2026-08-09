import pytest

from liftassess.chain import chain_candidate_id
from liftassess.models import (
    AssemblyIdentifier,
    EvidenceKind,
    EvidenceObservation,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.net import NetClassification, NetRecord, NetRecordKind
from liftassess.net_evidence import _annotate_candidate_with_net_records


@pytest.fixture
def source_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier("sourceAsm", "test-provider")


@pytest.fixture
def target_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier("targetAsm", "test-provider")


@pytest.fixture
def alignment_provenance() -> ProvenanceSource:
    return ProvenanceSource("alignment", "upstream alignment")


@pytest.fixture
def chain_provenance(alignment_provenance: ProvenanceSource) -> ProvenanceSource:
    return ProvenanceSource(
        "chain-file", "chain file", derived_from=(alignment_provenance,)
    )


@pytest.fixture
def net_provenance(alignment_provenance: ProvenanceSource) -> ProvenanceSource:
    return ProvenanceSource(
        "net-file", "net file", derived_from=(alignment_provenance,)
    )


def _candidate(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    *,
    chain_id: int = 7,
    orientation: MappingOrientation = MappingOrientation.SAME,
    source_segments: tuple[tuple[int, int], ...] = ((100, 120),),
) -> NormalizedCandidate:
    target_segments: list[MappingSegment] = []
    for index, (source_start, source_end) in enumerate(source_segments):
        length = source_end - source_start
        if orientation is MappingOrientation.SAME:
            target_start = 500 + index * 100
        else:
            target_start = 700 - index * 100 - length
        target_segments.append(
            MappingSegment(
                source_interval=GenomicInterval(
                    source_assembly, "chr1", source_start, source_end
                ),
                target_interval=GenomicInterval(
                    target_assembly,
                    "chrA",
                    target_start,
                    target_start + length,
                ),
            )
        )

    return NormalizedCandidate(
        candidate_id=chain_candidate_id(chain_provenance.source_id, chain_id),
        target_interval=GenomicInterval(
            target_assembly,
            "chrA",
            min(segment.target_interval.start for segment in target_segments),
            max(segment.target_interval.end for segment in target_segments),
        ),
        orientation=orientation,
        mapping_provenance=chain_provenance,
        segments=tuple(target_segments),
        evidence=(
            EvidenceObservation(
                observation_id=f"{chain_provenance.source_id}:chain:{chain_id}:score",
                kind=EvidenceKind.CHAIN_SCORE,
                value=1000,
                provenance=chain_provenance,
            ),
        ),
    )


def _fill(
    *,
    chain_id: int = 7,
    target_name: str = "chr1",
    target_start: int = 95,
    target_span_size: int = 40,
    query_name: str = "chrA",
    orientation: MappingOrientation = MappingOrientation.SAME,
    depth: int = 1,
    aligned_bases: int | None = 31,
    duplicated_query_bases: int | None = 4,
    classification: NetClassification | None = NetClassification.TOP,
) -> NetRecord:
    return NetRecord(
        target_name=target_name,
        target_sequence_size=1000,
        depth=depth,
        kind=NetRecordKind.FILL,
        target_start=target_start,
        target_span_size=target_span_size,
        query_name=query_name,
        orientation=orientation,
        query_start=400,
        query_span_size=target_span_size,
        chain_id=chain_id,
        score=1000,
        aligned_bases=aligned_bases,
        duplicated_query_bases=duplicated_query_bases,
        classification=classification,
    )


def test_matching_fill_attaches_raw_metrics_and_hierarchy_with_shared_fill_provenance(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(source_assembly, target_assembly, chain_provenance)

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=(_fill(depth=3),),
        net_provenance=net_provenance,
    )

    assert [observation.kind for observation in annotated.evidence] == [
        EvidenceKind.CHAIN_SCORE,
        EvidenceKind.ALIGNED_BASES,
        EvidenceKind.DUPLICATED_QUERY_BASES,
        EvidenceKind.NET_CLASSIFICATION,
        EvidenceKind.NET_HIERARCHY,
    ]
    assert [observation.value for observation in annotated.evidence[1:4]] == [
        31,
        4,
        "top",
    ]

    hierarchy = annotated.evidence[4].value
    assert isinstance(hierarchy, NetHierarchySummary)
    assert hierarchy.depth == 3
    assert hierarchy.source_fill_interval == GenomicInterval(
        source_assembly, "chr1", 95, 135
    )

    assert annotated.evidence[0].provenance is chain_provenance
    fill_source = annotated.evidence[1].provenance
    assert all(
        observation.provenance is fill_source for observation in annotated.evidence[1:]
    )
    assert fill_source.derived_from == (net_provenance,)
    assert net_provenance.derived_from == chain_provenance.derived_from


def test_repeated_chain_id_matches_by_actual_source_segment_overlap_not_id_alone(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(
        source_assembly,
        target_assembly,
        chain_provenance,
        source_segments=((300, 320),),
    )
    records = (
        _fill(
            chain_id=7,
            target_start=90,
            target_span_size=50,
            depth=3,
            aligned_bases=11,
            classification=NetClassification.NON_SYNTENIC,
        ),
        _fill(
            chain_id=7,
            target_start=290,
            target_span_size=50,
            depth=1,
            aligned_bases=37,
            classification=NetClassification.TOP,
        ),
    )

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=records,
        net_provenance=net_provenance,
    )

    aligned = [
        observation.value
        for observation in annotated.evidence
        if observation.kind is EvidenceKind.ALIGNED_BASES
    ]
    classifications = [
        observation.value
        for observation in annotated.evidence
        if observation.kind is EvidenceKind.NET_CLASSIFICATION
    ]
    hierarchy = [
        observation.value
        for observation in annotated.evidence
        if observation.kind is EvidenceKind.NET_HIERARCHY
    ]

    assert aligned == [37]
    assert classifications == ["top"]
    assert len(hierarchy) == 1
    assert isinstance(hierarchy[0], NetHierarchySummary)
    assert hierarchy[0].depth == 1


def test_one_candidate_can_preserve_multiple_relevant_fills_for_same_chain(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(
        source_assembly,
        target_assembly,
        chain_provenance,
        source_segments=((100, 110), (300, 310)),
    )
    records = (
        _fill(
            target_start=95,
            target_span_size=20,
            depth=3,
            aligned_bases=9,
            classification=NetClassification.NON_SYNTENIC,
        ),
        _fill(
            target_start=295,
            target_span_size=20,
            depth=1,
            aligned_bases=17,
            classification=NetClassification.TOP,
        ),
    )

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=records,
        net_provenance=net_provenance,
    )

    assert [
        observation.value
        for observation in annotated.evidence
        if observation.kind is EvidenceKind.ALIGNED_BASES
    ] == [9, 17]
    assert [
        observation.value
        for observation in annotated.evidence
        if observation.kind is EvidenceKind.NET_CLASSIFICATION
    ] == ["nonSyn", "top"]

    hierarchy_observations = [
        observation
        for observation in annotated.evidence
        if observation.kind is EvidenceKind.NET_HIERARCHY
    ]
    assert [
        observation.value.depth
        for observation in hierarchy_observations
        if isinstance(observation.value, NetHierarchySummary)
    ] == [3, 1]
    assert hierarchy_observations[0].provenance is not hierarchy_observations[1].provenance
    assert (
        hierarchy_observations[0].provenance.derived_from
        == hierarchy_observations[1].provenance.derived_from
        == (net_provenance,)
    )


def test_fill_overlapping_only_unaligned_space_between_segments_is_not_matched(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(
        source_assembly,
        target_assembly,
        chain_provenance,
        source_segments=((100, 110), (120, 130)),
    )

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=(_fill(target_start=110, target_span_size=10),),
        net_provenance=net_provenance,
    )

    assert annotated is candidate


def test_nonmatching_kind_sequence_orientation_or_chain_are_ignored(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(source_assembly, target_assembly, chain_provenance)
    records = (
        NetRecord(
            target_name="chr1",
            target_sequence_size=1000,
            depth=1,
            kind=NetRecordKind.GAP,
            target_start=95,
            target_span_size=40,
            query_name="chrA",
            orientation=MappingOrientation.SAME,
            query_start=400,
            query_span_size=40,
            chain_id=7,
        ),
        _fill(target_name="chr2"),
        _fill(chain_id=8),
        _fill(query_name="chrB"),
        _fill(orientation=MappingOrientation.REVERSE),
    )

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=records,
        net_provenance=net_provenance,
    )

    assert annotated is candidate


def test_fill_without_optional_metrics_still_records_hierarchy(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(source_assembly, target_assembly, chain_provenance)

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=(
            _fill(
                aligned_bases=None,
                duplicated_query_bases=None,
                classification=None,
            ),
        ),
        net_provenance=net_provenance,
    )

    assert [observation.kind for observation in annotated.evidence] == [
        EvidenceKind.CHAIN_SCORE,
        EvidenceKind.NET_HIERARCHY,
    ]


def test_rejects_chain_id_that_does_not_match_candidate_identity(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(source_assembly, target_assembly, chain_provenance)

    with pytest.raises(ValueError, match="candidate identity"):
        _annotate_candidate_with_net_records(
            candidate,
            chain_id=8,
            net_records=(_fill(),),
            net_provenance=net_provenance,
        )


def test_rejects_net_provenance_with_no_shared_upstream_source(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(source_assembly, target_assembly, chain_provenance)
    unrelated_net = ProvenanceSource("net-file", "unrelated net file")

    with pytest.raises(ValueError, match="share an upstream source"):
        _annotate_candidate_with_net_records(
            candidate,
            chain_id=7,
            net_records=(_fill(),),
            net_provenance=unrelated_net,
        )


def test_shared_upstream_provenance_handles_diamond_graph(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    root = ProvenanceSource("alignment", "upstream alignment")
    left = ProvenanceSource("left-branch", "left branch", derived_from=(root,))
    right = ProvenanceSource("right-branch", "right branch", derived_from=(root,))
    chain_provenance = ProvenanceSource(
        "chain-file",
        "chain file",
        derived_from=(left, right),
    )
    net_provenance = ProvenanceSource(
        "net-file",
        "net file",
        derived_from=(right,),
    )
    candidate = _candidate(source_assembly, target_assembly, chain_provenance)

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=(_fill(),),
        net_provenance=net_provenance,
    )

    assert len(annotated.evidence) > len(candidate.evidence)


def test_reverse_orientation_fill_matches_reverse_candidate(
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_provenance: ProvenanceSource,
) -> None:
    candidate = _candidate(
        source_assembly,
        target_assembly,
        chain_provenance,
        orientation=MappingOrientation.REVERSE,
    )

    annotated = _annotate_candidate_with_net_records(
        candidate,
        chain_id=7,
        net_records=(
            _fill(
                orientation=MappingOrientation.REVERSE,
                classification=NetClassification.INVERSION,
            ),
        ),
        net_provenance=net_provenance,
    )

    assert [
        observation.value
        for observation in annotated.evidence
        if observation.kind is EvidenceKind.NET_CLASSIFICATION
    ] == ["inv"]
