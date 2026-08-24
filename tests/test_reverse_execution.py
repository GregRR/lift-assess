from __future__ import annotations

import gzip
from collections.abc import Iterator
from io import StringIO
from pathlib import Path

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCResourceBundle,
    EvidenceAvailabilityTier,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
    ReverseOriginalSourceCoverageState,
    ReverseRelationshipState,
    build_cached_chain_index,
    build_reverse_mapping_results_from_cached_bundle,
    sha256_identifier_for_file,
    ucsc_resource_terms,
)
from liftassess.chain import ChainRecord, iter_chain_records
from liftassess.engine import build_ucsc_chain_candidates_for_intervals
from liftassess.resource_files import (
    build_ucsc_chain_candidates_for_intervals_from_cached_bundle,
)

SOURCE_ASSEMBLY = AssemblyIdentifier("canFam3", "UCSC")
TARGET_ASSEMBLY = AssemblyIdentifier("canFam4", "UCSC")
FORWARD_PROVENANCE = ProvenanceSource("forward-chain", "forward chain")
REVERSE_ALIGNMENT = ProvenanceSource("reverse-alignment", "reverse UCSC alignment")


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="ascii", newline="\n") as handle:
        handle.write(text)


def _cached_resource(path: Path, url: str) -> CachedResource:
    return CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-23T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=False,
    )


def _reverse_liftover_bundle(
    tmp_path: Path, chain_text: str
) -> CachedUCSCResourceBundle:
    chain_path = tmp_path / "reverse-chain"
    _write_gzip(chain_path, chain_text)
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/"
        "canFam4ToCanFam3.over.chain.gz"
    )
    return CachedUCSCResourceBundle(
        source_db="canFam4",
        target_db="canFam3",
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=_cached_resource(chain_path, url),
    )


def _forward_candidate(
    *,
    source_spans: tuple[tuple[int, int], ...] = ((100, 110),),
    target_spans: tuple[tuple[int, int], ...] = ((500, 510),),
) -> NormalizedCandidate:
    segments = tuple(
        MappingSegment(
            GenomicInterval(SOURCE_ASSEMBLY, "chr1", source_start, source_end),
            GenomicInterval(TARGET_ASSEMBLY, "chrA", target_start, target_end),
        )
        for (source_start, source_end), (target_start, target_end) in zip(
            source_spans, target_spans, strict=True
        )
    )
    return NormalizedCandidate(
        candidate_id="forward",
        target_interval=GenomicInterval(
            TARGET_ASSEMBLY,
            "chrA",
            min(start for start, _ in target_spans),
            max(end for _, end in target_spans),
        ),
        orientation=MappingOrientation.SAME,
        mapping_provenance=FORWARD_PROVENANCE,
        segments=segments,
    )


class _OneShotChains:
    def __init__(self, records: tuple[ChainRecord, ...]) -> None:
        self.records = records
        self.iteration_count = 0

    def __iter__(self) -> Iterator[ChainRecord]:
        self.iteration_count += 1
        if self.iteration_count > 1:
            raise AssertionError("chain records were traversed more than once")
        yield from self.records


def test_chain_multi_interval_engine_consumes_chain_stream_once() -> None:
    records = tuple(
        iter_chain_records(
            StringIO(
                "chain 10 chrA 1000 + 100 110 chr1 1000 + 500 510 1\n"
                "10\n\n"
                "chain 20 chrA 1000 + 200 210 chr2 1000 + 600 610 2\n"
                "10\n\n"
            )
        )
    )
    chains = _OneShotChains(records)
    provenance = ProvenanceSource("reverse-chain", "reverse chain")

    candidates = build_ucsc_chain_candidates_for_intervals(
        (
            GenomicInterval(TARGET_ASSEMBLY, "chrA", 102, 108),
            GenomicInterval(TARGET_ASSEMBLY, "chrA", 202, 208),
        ),
        chains,
        target_assembly=SOURCE_ASSEMBLY,
        chain_provenance=provenance,
    )

    assert chains.iteration_count == 1
    assert tuple(candidate.candidate_id for candidate in candidates[0]) == (
        "reverse-chain:chain:1",
    )
    assert tuple(candidate.candidate_id for candidate in candidates[1]) == (
        "reverse-chain:chain:2",
    )


def test_actual_reverse_execution_matches_indexed_and_shared_traversal(
    tmp_path: Path,
) -> None:
    bundle = _reverse_liftover_bundle(
        tmp_path,
        "chain 100 chrA 2000 + 500 510 chr1 1000 + 100 110 1\n"
        "10\n\n"
        "chain 90 chrA 2000 + 500 510 chr9 1000 + 700 710 2\n"
        "10\n\n",
    )
    forward = _forward_candidate()

    shared = build_reverse_mapping_results_from_cached_bundle(
        (forward,),
        bundle,
        reverse_alignment_provenance=REVERSE_ALIGNMENT,
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index
    indexed = build_reverse_mapping_results_from_cached_bundle(
        (forward,),
        bundle,
        reverse_alignment_provenance=REVERSE_ALIGNMENT,
        chain_index=index,
    )

    assert indexed == shared
    (result,) = indexed
    assert result.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE
    assert result.original_source_covered_bases == 10
    assert (
        result.original_source_coverage is ReverseOriginalSourceCoverageState.COMPLETE
    )
    assert result.exact_original_geometry_return
    assert result.reverse_projection_count == 2
    for candidate in result.segment_results[0].candidates:
        assert candidate.mapping_provenance.identifiers[0].value == bundle.chain.sha256
        assert candidate.mapping_provenance.derived_from == (REVERSE_ALIGNMENT,)


def test_actual_reverse_execution_queries_fragmented_segments_not_bounding_span(
    tmp_path: Path,
) -> None:
    bundle = _reverse_liftover_bundle(
        tmp_path,
        "chain 100 chrA 2000 + 1000 1040 chr1 1000 + 100 140 1\n"
        "40\n\n"
        "chain 100 chrA 2000 + 1060 1100 chr1 1000 + 160 200 2\n"
        "40\n\n"
        "chain 100 chrA 2000 + 1045 1055 chr9 1000 + 700 710 3\n"
        "10\n\n",
    )
    forward = _forward_candidate(
        source_spans=((100, 140), (160, 200)),
        target_spans=((1000, 1040), (1060, 1100)),
    )

    (result,) = build_reverse_mapping_results_from_cached_bundle(
        (forward,),
        bundle,
        reverse_alignment_provenance=REVERSE_ALIGNMENT,
    )

    assert result.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY
    assert result.reverse_projection_count == 2
    assert result.exact_original_geometry_return
    assert tuple(
        segment.queried_target_segment for segment in result.segment_results
    ) == tuple(segment.target_interval for segment in forward.segments)


def test_chain_only_cached_multi_interval_path_ignores_comparative_resources(
    tmp_path: Path,
) -> None:
    chain_path = tmp_path / "chain"
    _write_gzip(
        chain_path,
        "chain 100 chrA 2000 + 500 510 chr1 1000 + 100 110 1\n10\n\n",
    )
    invalid_paths = {
        name: tmp_path / name for name in ("net", "syn-net", "rbest-chain", "rbest-net")
    }
    for path in invalid_paths.values():
        path.write_bytes(b"deliberately not parseable")

    forward_base = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
    reciprocal_base = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/reciprocalBest/"
    )
    bundle = CachedUCSCResourceBundle(
        source_db="canFam4",
        target_db="canFam3",
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain=_cached_resource(
            chain_path, f"{forward_base}canFam4.canFam3.all.chain.gz"
        ),
        net=_cached_resource(
            invalid_paths["net"], f"{forward_base}canFam4.canFam3.net.gz"
        ),
        syntenic_net=_cached_resource(
            invalid_paths["syn-net"], f"{forward_base}canFam4.canFam3.syn.net.gz"
        ),
        reciprocal_best_chain=_cached_resource(
            invalid_paths["rbest-chain"],
            f"{reciprocal_base}canFam4.canFam3.rbest.chain.gz",
        ),
        reciprocal_best_net=_cached_resource(
            invalid_paths["rbest-net"],
            f"{reciprocal_base}canFam4.canFam3.rbest.net.gz",
        ),
    )

    candidates = build_ucsc_chain_candidates_for_intervals_from_cached_bundle(
        (GenomicInterval(TARGET_ASSEMBLY, "chrA", 502, 508),),
        bundle,
        target_assembly=SOURCE_ASSEMBLY,
        alignment_provenance=REVERSE_ALIGNMENT,
    )

    assert len(candidates) == 1
    assert len(candidates[0]) == 1
