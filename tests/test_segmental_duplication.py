import gzip
from pathlib import Path

import pytest

from liftassess import ExternalContextState, build_external_context_profile
from liftassess.models import (
    AssemblyIdentifier,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
)
from liftassess.resource_cache import CachedResource, ucsc_resource_terms
from liftassess.resource_identity import (
    ResourceIdentityMismatchError,
    sha256_identifier_for_file,
)
from liftassess.resources import UCSCSegmentalDuplicationResource
from liftassess.segmental_duplication import (
    SegmentalDuplicationCheckState,
    UCSCSegmentalDuplicationCatalog,
    acquire_ucsc_segmental_duplication_resource,
    build_cached_ucsc_segmental_duplication_catalog,
    build_ucsc_segmental_duplication_context,
    iter_ucsc_genomic_super_dups,
    ucsc_segmental_duplication_table_url,
)

ASSEMBLY = AssemblyIdentifier(name="hg38", provider="UCSC")
SOURCE_ASSEMBLY = AssemblyIdentifier(name="hg19", provider="UCSC")

_SAMPLE_ROW = (
    "585\tchr1\t10000\t87112\tchr15:101906152\t0\t-\tchr15\t101906152\t"
    "101981189\t75037\t11764\t1000\tN/A\tN/A\tN/A\tN/A\t"
    "align_both/0009/both0046049\t77880\t71\t3611\t74269\t73743\t526\t"
    "331\t195\t0.992918\t0.991969\t0.00711601\t0.00711937\n"
)


def _record_row(
    *,
    chrom: str,
    start: int,
    end: int,
    other_chrom: str,
    other_start: int,
    other_end: int,
    uid: int,
    strand: str = "+",
    frac_match: float = 0.99,
) -> str:
    fields = [
        "0",
        chrom,
        str(start),
        str(end),
        f"{other_chrom}:{other_start}",
        "0",
        strand,
        other_chrom,
        str(other_start),
        str(other_end),
        str(other_end),
        str(uid),
        "1000",
        "N/A",
        "N/A",
        "N/A",
        "N/A",
        "alignment",
        str(end - start),
        "0",
        "0",
        str(end - start),
        str(end - start),
        "0",
        "0",
        "0",
        str(frac_match),
        str(frac_match),
        "0",
        "0",
    ]
    assert len(fields) == 30
    return "\t".join(fields) + "\n"


def _catalog(
    assembly: AssemblyIdentifier,
    rows: tuple[str, ...],
    *,
    source_id: str,
) -> UCSCSegmentalDuplicationCatalog:
    return UCSCSegmentalDuplicationCatalog(
        assembly=assembly,
        records=tuple(iter_ucsc_genomic_super_dups(rows, assembly=assembly)),
        provenance=ProvenanceSource(source_id=source_id, label=source_id),
    )


def _candidate(*segments: MappingSegment) -> NormalizedCandidate:
    target_start = min(segment.target_interval.start for segment in segments)
    target_end = max(segment.target_interval.end for segment in segments)
    return NormalizedCandidate(
        candidate_id="chain:1",
        target_interval=GenomicInterval(
            assembly=ASSEMBLY,
            sequence_name="chr1",
            start=target_start,
            end=target_end,
        ),
        orientation=MappingOrientation.SAME,
        mapping_provenance=ProvenanceSource(source_id="chain", label="chain"),
        segments=segments,
    )


def test_parses_ucsc_genomic_super_dups_schema_fields() -> None:
    (record,) = tuple(iter_ucsc_genomic_super_dups((_SAMPLE_ROW,), assembly=ASSEMBLY))

    assert record.interval == GenomicInterval(ASSEMBLY, "chr1", 10000, 87112)
    assert record.paired_interval == GenomicInterval(
        ASSEMBLY, "chr15", 101906152, 101981189
    )
    assert record.strand == "-"
    assert record.uid == 11764
    assert record.aligned_bases == 74269
    assert record.fraction_matching_bases == pytest.approx(0.992918)


def test_rejects_malformed_ucsc_genomic_super_dups_row() -> None:
    with pytest.raises(ValueError, match="exactly 30"):
        tuple(iter_ucsc_genomic_super_dups(("chr1\t1\t2\n",), assembly=ASSEMBLY))


def test_catalog_returns_exact_half_open_overlaps() -> None:
    catalog = _catalog(
        ASSEMBLY,
        (
            _record_row(
                chrom="chr1",
                start=100,
                end=200,
                other_chrom="chr2",
                other_start=500,
                other_end=600,
                uid=1,
            ),
            _record_row(
                chrom="chr1",
                start=300,
                end=400,
                other_chrom="chr3",
                other_start=700,
                other_end=800,
                uid=2,
            ),
        ),
        source_id="segdup:hg38",
    )

    overlaps = catalog.overlapping(GenomicInterval(ASSEMBLY, "chr1", 150, 320))

    assert [item.record.uid for item in overlaps] == [1, 2]
    assert [item.overlap_interval for item in overlaps] == [
        GenomicInterval(ASSEMBLY, "chr1", 150, 200),
        GenomicInterval(ASSEMBLY, "chr1", 300, 320),
    ]


def test_target_context_queries_exact_mapped_segments_not_bounding_gap() -> None:
    source_catalog = _catalog(SOURCE_ASSEMBLY, (), source_id="segdup:hg19")
    target_catalog = _catalog(
        ASSEMBLY,
        (
            _record_row(
                chrom="chr1",
                start=150,
                end=160,
                other_chrom="chr2",
                other_start=500,
                other_end=510,
                uid=1,
            ),
            _record_row(
                chrom="chr1",
                start=450,
                end=460,
                other_chrom="chr3",
                other_start=700,
                other_end=710,
                uid=2,
            ),
        ),
        source_id="segdup:hg38",
    )
    source = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 1000, 1020)
    candidate = _candidate(
        MappingSegment(
            source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 1000, 1010),
            target_interval=GenomicInterval(ASSEMBLY, "chr1", 100, 110),
        ),
        MappingSegment(
            source_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 1010, 1020),
            target_interval=GenomicInterval(ASSEMBLY, "chr1", 400, 410),
        ),
    )

    result = build_ucsc_segmental_duplication_context(
        source,
        (candidate,),
        source_catalog=source_catalog,
        target_catalog=target_catalog,
        source_unavailable=False,
        target_unavailable=False,
    )

    assert result.source_state is SegmentalDuplicationCheckState.ASSESSED
    assert result.target_state is SegmentalDuplicationCheckState.ASSESSED
    assert result.source_overlaps == ()
    assert result.target_overlaps == ()


def test_context_records_source_and_target_overlaps_without_scoring_mapping() -> None:
    source_catalog = _catalog(
        SOURCE_ASSEMBLY,
        (
            _record_row(
                chrom="chr1",
                start=995,
                end=1010,
                other_chrom="chr5",
                other_start=2000,
                other_end=2015,
                uid=10,
                frac_match=0.995,
            ),
        ),
        source_id="segdup:hg19",
    )
    target_catalog = _catalog(
        ASSEMBLY,
        (
            _record_row(
                chrom="chr1",
                start=105,
                end=115,
                other_chrom="chr7",
                other_start=3000,
                other_end=3010,
                uid=20,
                frac_match=0.999,
            ),
        ),
        source_id="segdup:hg38",
    )
    source = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 1000, 1020)
    candidate = _candidate(
        MappingSegment(
            source_interval=source,
            target_interval=GenomicInterval(ASSEMBLY, "chr1", 100, 120),
        )
    )

    result = build_ucsc_segmental_duplication_context(
        source,
        (candidate,),
        source_catalog=source_catalog,
        target_catalog=target_catalog,
        source_unavailable=False,
        target_unavailable=False,
    )

    assert len(result.source_overlaps) == 1
    assert result.source_overlaps[0].overlap_interval == GenomicInterval(
        SOURCE_ASSEMBLY, "chr1", 1000, 1010
    )
    assert len(result.target_overlaps) == 1
    assert result.target_overlaps[0].candidate_id == candidate.candidate_id
    assert result.target_overlaps[0].overlap_intervals == (
        GenomicInterval(ASSEMBLY, "chr1", 105, 115),
    )


def test_no_projection_target_context_is_not_reported_unavailable() -> None:
    source_catalog = _catalog(SOURCE_ASSEMBLY, (), source_id="segdup:hg19")
    source = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 1000, 1020)

    result = build_ucsc_segmental_duplication_context(
        source,
        (),
        source_catalog=source_catalog,
        target_catalog=None,
        source_unavailable=False,
        target_unavailable=True,
    )

    assert result.source_state is SegmentalDuplicationCheckState.ASSESSED
    assert result.target_state is SegmentalDuplicationCheckState.NO_TARGET_PROJECTIONS


def test_cached_catalog_verifies_exact_resource_bytes(tmp_path: Path) -> None:
    path = tmp_path / "genomicSuperDups.txt.gz"
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(_SAMPLE_ROW)
    url = ucsc_segmental_duplication_table_url("hg38")
    resource = CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-28T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=True,
    )

    catalog = build_cached_ucsc_segmental_duplication_catalog(ASSEMBLY, resource)
    assert len(catalog.records) == 1
    assert catalog.provenance.source_id == f"file:{resource.sha256}"

    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(_SAMPLE_ROW.replace("\t11764\t", "\t11765\t"))
    with pytest.raises(ResourceIdentityMismatchError):
        build_cached_ucsc_segmental_duplication_catalog(ASSEMBLY, resource)


def test_acquisition_uses_discovered_exact_url_without_terms_ack_gate(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    discovered = UCSCSegmentalDuplicationResource(
        db="hg38",
        url=ucsc_segmental_duplication_table_url("hg38"),
    )
    expected = object()
    calls: list[tuple[object, ...]] = []

    def fake_acquire(url: str, cache_root: object, **kwargs: object) -> object:
        calls.append((url, cache_root, kwargs))
        return expected

    monkeypatch.setattr(
        "liftassess.segmental_duplication.acquire_ucsc_resource",
        fake_acquire,
    )

    result = acquire_ucsc_segmental_duplication_resource(discovered, tmp_path)

    assert result is expected
    assert calls == [
        (
            discovered.url,
            tmp_path,
            {"terms_acknowledged": False, "refresh": False},
        )
    ]


def test_external_context_scope_distinguishes_partial_and_unavailable() -> None:
    source = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 1000, 1020)
    target = _candidate(
        MappingSegment(
            source_interval=source,
            target_interval=GenomicInterval(ASSEMBLY, "chr1", 100, 120),
        )
    )
    source_catalog = _catalog(SOURCE_ASSEMBLY, (), source_id="segdup:hg19")

    partial_result = build_ucsc_segmental_duplication_context(
        source,
        (target,),
        source_catalog=source_catalog,
        target_catalog=None,
        source_unavailable=False,
        target_unavailable=True,
    )
    assert (
        build_external_context_profile(partial_result).state
        is ExternalContextState.PARTIALLY_ASSESSED
    )

    unavailable_result = build_ucsc_segmental_duplication_context(
        source,
        (target,),
        source_catalog=None,
        target_catalog=None,
        source_unavailable=True,
        target_unavailable=True,
    )
    assert (
        build_external_context_profile(unavailable_result).state
        is ExternalContextState.UNAVAILABLE
    )
