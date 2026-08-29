from __future__ import annotations

import gzip
from dataclasses import replace
from pathlib import Path

import pytest

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCChainResource,
    CachedUCSCResourceBundle,
    ComparativeEvidenceRelationship,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    FactualHeadline,
    FilteredAllChainCorrespondenceError,
    FilteredAllChainInventoryState,
    GenomicInterval,
    InputValidityState,
    MappingSegment,
    ProvenanceSource,
    QueryContextFinding,
    QueryContextNotRunReason,
    QueryContextState,
    SourceIntervalPreflightResult,
    UCSCAssessmentReport,
    UCSCBundleResourceRole,
    assess_ucsc_cached_bundle,
    attach_filtered_all_chain_comparison,
    attach_point_query_context,
    attach_query_context_result,
    attach_reverse_mapping_results,
    build_cached_chain_index,
    build_ucsc_assembly_sequence_catalog,
    preflight_source_interval,
    provenance_source_for_file,
    reverse_mapping_unavailable,
    sha256_identifier_for_file,
    ucsc_resource_terms,
)
from liftassess.reporting import render_assessment_summary


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _chain_text(
    *,
    chain_id: int = 1,
    start: int = 100,
    length: int = 20,
    target_start: int = 500,
) -> str:
    end = start + length
    target_end = target_start + length
    return (
        f"chain 100 chr1 1000 + {start} {end} chrA 2000 + "
        f"{target_start} {target_end} {chain_id}\n{length}\n\n"
    )


def _net_text() -> str:
    return (
        "net chr1 1000\n"
        " fill 100 20 chrA + 500 20 id 1 score 100 ali 20 qDup 0 type syn\n"
    )


def _cached_resource(path: Path, url: str) -> CachedResource:
    return CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-16T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=False,
    )


def _comparative_bundle(
    tmp_path: Path,
    *,
    chain_start: int = 100,
    chain_text: str | None = None,
    reciprocal_chain_text: str | None = None,
    valid_net: bool = True,
) -> CachedUCSCResourceBundle:
    chain_path = tmp_path / "chain-artifact"
    net_path = tmp_path / "net-artifact"
    syn_net_path = tmp_path / "syn-net-artifact"
    reciprocal_chain_path = tmp_path / "rbest-chain-artifact"
    reciprocal_net_path = tmp_path / "rbest-net-artifact"

    _write_gzip(
        chain_path,
        chain_text or _chain_text(chain_id=1, start=chain_start),
    )
    if valid_net:
        _write_gzip(net_path, _net_text())
    else:
        net_path.write_bytes(b"not a net")
    _write_gzip(
        reciprocal_chain_path,
        reciprocal_chain_text or _chain_text(chain_id=101),
    )
    # These resources are present in the acquired bundle but are deliberately not
    # parseable by the current engine. A successful report must not imply otherwise.
    syn_net_path.write_bytes(b"retrieval context only")
    reciprocal_net_path.write_bytes(b"retrieval context only")

    forward = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reciprocal = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/"
    )
    return CachedUCSCResourceBundle(
        source_db="canFam3",
        target_db="canFam4",
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain=_cached_resource(
            chain_path,
            f"{forward}canFam3.canFam4.all.chain.gz",
        ),
        net=_cached_resource(net_path, f"{forward}canFam3.canFam4.net.gz"),
        syntenic_net=_cached_resource(
            syn_net_path,
            f"{forward}canFam3.canFam4.syn.net.gz",
        ),
        reciprocal_best_chain=_cached_resource(
            reciprocal_chain_path,
            f"{reciprocal}canFam3.canFam4.rbest.chain.gz",
        ),
        reciprocal_best_net=_cached_resource(
            reciprocal_net_path,
            f"{reciprocal}canFam3.canFam4.rbest.net.gz",
        ),
    )


def _liftover_bundle(tmp_path: Path, chain_text: str) -> CachedUCSCResourceBundle:
    chain_path = tmp_path / "lift-chain"
    _write_gzip(chain_path, chain_text)
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
        "canFam3ToCanFam4.over.chain.gz"
    )
    return CachedUCSCResourceBundle(
        source_db="canFam3",
        target_db="canFam4",
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=_cached_resource(chain_path, url),
    )


def _source_preflight(
    tmp_path: Path,
    source: AssemblyIdentifier,
    *,
    sequence_name: str,
    sequence_length: int,
    start: int,
    end: int,
) -> tuple[GenomicInterval, SourceIntervalPreflightResult, CachedResource]:
    path = tmp_path / f"{sequence_name}-chromInfo"
    _write_gzip(
        path,
        f"{sequence_name}\t{sequence_length}\t/gbdb/canFam3/canFam3.2bit\n",
    )
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/chromInfo.txt.gz"
    resource = _cached_resource(path, url)
    provenance = provenance_source_for_file(
        path,
        label="UCSC canFam3 chromInfo assembly-sequence metadata",
        derived_from=(),
    )
    catalog = build_ucsc_assembly_sequence_catalog(
        source,
        (f"{sequence_name}\t{sequence_length}\t/gbdb/canFam3/canFam3.2bit\n",),
        sequence_provenance=provenance,
    )
    interval = GenomicInterval(source, sequence_name, start, end)
    return interval, preflight_source_interval(interval, catalog), resource


def _consumed_roles(report: UCSCAssessmentReport) -> set[UCSCBundleResourceRole]:
    return {
        resource.role for resource in report.resources if resource.consumed_by_engine
    }


def _assemblies() -> tuple[AssemblyIdentifier, AssemblyIdentifier]:
    return (
        AssemblyIdentifier(name="canFam3", provider="UCSC"),
        AssemblyIdentifier(name="canFam4", provider="UCSC"),
    )


def test_filtered_all_chain_comparison_uses_prepared_filtered_index(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    source_interval = GenomicInterval(source, "chr1", 105, 115)
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    report = assess_ucsc_cached_bundle(
        source_interval,
        _comparative_bundle(tmp_path),
        target_assembly=target,
        alignment_provenance=alignment,
    )
    filtered_bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=91),
    )
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    filtered_index = build_cached_chain_index(
        tmp_path / "filtered-index-cache",
        filtered_bundle.chain,
    ).index

    def fail_full_scan(*args: object, **kwargs: object) -> tuple[object, ...]:
        del args, kwargs
        raise AssertionError("filtered comparison must not start a whole-chain scan")

    monkeypatch.setattr(
        "liftassess.resource_files._iter_chain_file_with_provenance",
        fail_full_scan,
    )

    enriched = attach_filtered_all_chain_comparison(
        report,
        filtered_chain=filtered_chain,
        filtered_chain_index=filtered_index,
    )

    comparison = enriched.filtered_all_chain_comparison
    assert comparison is not None
    assert (
        comparison.relationship
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    )
    assert len(comparison.filtered_candidates) == 1
    assert comparison.candidate_matches[0].all_chain_candidate_id == (
        report.candidates[0].candidate_id
    )
    assert enriched.filtered_chain_comparison_resource is not None
    assert enriched.comparative_evidence_relationship is not None
    assert (
        enriched.comparative_evidence_relationship.relationship
        is ComparativeEvidenceRelationship.NO_COMPETING_FULL_PLACEMENTS
    )
    assert (
        comparison.filtered_chain_provenance.derived_from
        == comparison.all_chain_provenance.derived_from
        == (alignment,)
    )
    assert (
        enriched.result_profile.scope.comparative_relationship.value
        == "NO_COMPETING_FULL_PLACEMENTS"
    )
    assert (
        enriched.result_profile.comparative_relationship.inventory_state
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    )
    assert enriched.result_profile.comparative_relationship.favored_candidate_id is None
    first_support = enriched.result_profile.comparative_relationship.placement_support[
        0
    ]
    assert first_support.candidate_id == report.candidates[0].candidate_id


def test_filtered_all_chain_comparison_reports_net_clipped_geometry_as_unpairable(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    source_interval = GenomicInterval(source, "chr1", 100, 120)
    alignment = ProvenanceSource("alignment", "shared UCSC dependency group")
    report = assess_ucsc_cached_bundle(
        source_interval,
        _comparative_bundle(tmp_path),
        target_assembly=target,
        alignment_provenance=alignment,
    )
    filtered_bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=1, start=105, length=10, target_start=505),
    )
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    filtered_index = build_cached_chain_index(
        tmp_path / "filtered-clipped-index",
        filtered_bundle.chain,
    ).index

    with pytest.raises(
        FilteredAllChainCorrespondenceError,
        match="cannot be paired to identical all-chain geometry",
    ):
        attach_filtered_all_chain_comparison(
            report,
            filtered_chain=filtered_chain,
            filtered_chain_index=filtered_index,
        )


def test_comparative_profile_survives_reverse_context_rebuild(tmp_path: Path) -> None:
    source, target = _assemblies()
    source_interval = GenomicInterval(source, "chr1", 105, 115)
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    report = assess_ucsc_cached_bundle(
        source_interval,
        _comparative_bundle(tmp_path),
        target_assembly=target,
        alignment_provenance=alignment,
    )
    filtered_bundle = _liftover_bundle(tmp_path, _chain_text(chain_id=91))
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    filtered_index = build_cached_chain_index(
        tmp_path / "filtered-reverse-preservation-index",
        filtered_bundle.chain,
    ).index
    compared = attach_filtered_all_chain_comparison(
        report,
        filtered_chain=filtered_chain,
        filtered_chain_index=filtered_index,
    )

    enriched = attach_reverse_mapping_results(
        compared,
        tuple(
            reverse_mapping_unavailable(candidate) for candidate in compared.candidates
        ),
    )

    assert (
        enriched.result_profile.scope.comparative_relationship.value
        == "NO_COMPETING_FULL_PLACEMENTS"
    )
    assert (
        enriched.result_profile.comparative_relationship.inventory_state
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    )


def test_comparative_profile_survives_point_context_rebuild(tmp_path: Path) -> None:
    source, target = _assemblies()
    source_interval = GenomicInterval(source, "chr1", 105, 106)
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    bundle = _comparative_bundle(tmp_path)
    report = assess_ucsc_cached_bundle(
        source_interval,
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )
    filtered_bundle = _liftover_bundle(tmp_path, _chain_text(chain_id=91))
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    filtered_index = build_cached_chain_index(
        tmp_path / "filtered-context-preservation-index",
        filtered_bundle.chain,
    ).index
    compared = attach_filtered_all_chain_comparison(
        report,
        filtered_chain=filtered_chain,
        filtered_chain_index=filtered_index,
    )
    all_chain = CachedUCSCChainResource(
        source_db=bundle.source_db,
        target_db=bundle.target_db,
        evidence_tier=bundle.evidence_tier,
        chain=bundle.chain,
    )
    all_chain_index = build_cached_chain_index(
        tmp_path / "all-chain-context-preservation-index",
        bundle.chain,
    ).index

    enriched = attach_point_query_context(
        compared,
        chain_context=all_chain,
        chain_index=all_chain_index,
    )

    assert (
        enriched.result_profile.scope.comparative_relationship.value
        == "NO_COMPETING_FULL_PLACEMENTS"
    )
    assert (
        enriched.result_profile.comparative_relationship.inventory_state
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    )


def test_filtered_all_chain_comparison_requires_prepared_filtered_index(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        _comparative_bundle(tmp_path),
        target_assembly=target,
        alignment_provenance=ProvenanceSource(
            "alignment",
            "shared UCSC alignment",
        ),
    )
    filtered_bundle = _liftover_bundle(tmp_path, _chain_text(chain_id=91))
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    executed = False

    def fail_if_executed(*args: object, **kwargs: object) -> tuple[object, ...]:
        nonlocal executed
        del args, kwargs
        executed = True
        raise AssertionError("comparison execution must not start without an index")

    monkeypatch.setattr(
        (
            "liftassess.orchestration."
            "build_ucsc_chain_candidates_for_intervals_from_cached_chain"
        ),
        fail_if_executed,
    )

    with pytest.raises(ValueError, match="requires a prepared filtered-chain index"):
        attach_filtered_all_chain_comparison(
            report,
            filtered_chain=filtered_chain,
            filtered_chain_index=None,
        )

    assert executed is False


def test_cached_comparative_bundle_runs_engine_and_profile_end_to_end(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    source_interval = GenomicInterval(source, "chr1", 105, 115)
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    bundle = _comparative_bundle(tmp_path)

    report = assess_ucsc_cached_bundle(
        source_interval,
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )

    assert report.source_db == "canFam3"
    assert report.target_db == "canFam4"
    assert report.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE
    assert (
        report.result_profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION
    )
    assert report.result_profile.maximum_coverage_candidate_ids
    assert len(report.candidates) == 1
    candidate = report.candidates[0]
    assert any(
        observation.kind is EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
        for observation in candidate.evidence
    )


def test_assessment_with_chain_index_matches_full_report(tmp_path: Path) -> None:
    source, target = _assemblies()
    bundle = _comparative_bundle(tmp_path)
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    interval = GenomicInterval(source, "chr1", 105, 115)
    full = assess_ucsc_cached_bundle(
        interval,
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    indexed = assess_ucsc_cached_bundle(
        interval,
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
        chain_index=index,
    )

    assert indexed == full


def test_report_distinguishes_consumed_evidence_inputs_from_bundle_context(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    bundle = _comparative_bundle(tmp_path)

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )

    consumed = {
        resource.role for resource in report.resources if resource.consumed_by_engine
    }
    assert consumed == {
        UCSCBundleResourceRole.CHAIN,
        UCSCBundleResourceRole.NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
    }
    context_only = {
        resource.role
        for resource in report.resources
        if not resource.consumed_by_engine
    }
    assert context_only == {
        UCSCBundleResourceRole.SYNTENIC_NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_NET,
    }
    assert all(
        resource.file_provenance is None
        for resource in report.resources
        if not resource.consumed_by_engine
    )
    expected_resources = {
        UCSCBundleResourceRole.CHAIN: bundle.chain,
        UCSCBundleResourceRole.NET: bundle.net,
        UCSCBundleResourceRole.SYNTENIC_NET: bundle.syntenic_net,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN: bundle.reciprocal_best_chain,
        UCSCBundleResourceRole.RECIPROCAL_BEST_NET: bundle.reciprocal_best_net,
    }
    assert all(
        resource.resource is expected_resources[resource.role]
        for resource in report.resources
    )

    candidate = report.candidates[0]
    chain_resource = report.resources[0]
    assert chain_resource.file_provenance == candidate.mapping_provenance

    net_observation = next(
        observation
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.NET_CLASSIFICATION
    )
    net_resource = next(
        resource
        for resource in report.resources
        if resource.role is UCSCBundleResourceRole.NET
    )
    assert net_resource.file_provenance == net_observation.provenance.derived_from[0]

    reciprocal_observation = next(
        observation
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
    )
    reciprocal_resource = next(
        resource
        for resource in report.resources
        if resource.role is UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN
    )
    assert reciprocal_resource.file_provenance == reciprocal_observation.provenance


def test_no_candidate_report_does_not_claim_comparative_files_were_consumed(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    # The invalid net bytes are content-addressed honestly in the cached resource.
    # If orchestration tried to parse them despite having no candidates, profile
    # construction would fail instead of returning a report.
    bundle = _comparative_bundle(tmp_path, chain_start=200, valid_net=False)

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )

    assert report.result_profile.headline is FactualHeadline.NO_CHAIN_PROJECTION
    assert report.candidates == ()
    assert {
        resource.role for resource in report.resources if resource.consumed_by_engine
    } == {UCSCBundleResourceRole.CHAIN}


def test_liftover_only_report_keeps_tier_separate_from_result_state(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(tmp_path, _chain_text(chain_id=17))
    url = bundle.chain.source_url

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource(
            "alignment",
            "upstream UCSC alignment",
        ),
    )

    assert report.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
    assert (
        report.result_profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION
    )
    assert len(report.resources) == 1
    assert report.resources[0].consumed_by_engine
    assert report.resources[0].resource is bundle.chain
    assert report.resources[0].resource.source_url == url
    assert report.resources[0].resource.retrieved_at == "2026-08-16T00:00:00Z"
    assert report.resources[0].resource.terms.restricted_liftover_chain


def test_comparative_consumption_is_independent_of_multiple_projection_state(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _comparative_bundle(
        tmp_path,
        chain_text=(
            _chain_text(chain_id=1, target_start=500)
            + _chain_text(chain_id=2, target_start=700)
        ),
        reciprocal_chain_text=(
            _chain_text(chain_id=101, target_start=500)
            + _chain_text(chain_id=102, target_start=700)
        ),
    )

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )

    assert report.result_profile.headline is FactualHeadline.MULTIPLE_CHAIN_PROJECTIONS
    assert _consumed_roles(report) == {
        UCSCBundleResourceRole.CHAIN,
        UCSCBundleResourceRole.NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
    }


def test_comparative_consumption_is_independent_of_rbest_membership_state(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    # The reciprocal-best resource confirms only half of the assessed source bases.
    bundle = _comparative_bundle(
        tmp_path,
        reciprocal_chain_text=_chain_text(
            chain_id=101,
            start=105,
            length=5,
            target_start=505,
        ),
    )

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )

    assert (
        report.result_profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION
    )
    assert _consumed_roles(report) == {
        UCSCBundleResourceRole.CHAIN,
        UCSCBundleResourceRole.NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
    }


def test_liftover_only_multi_candidate_report_marks_only_chain_consumed(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=1, target_start=500)
        + _chain_text(chain_id=2, target_start=700),
    )

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "upstream UCSC alignment"),
    )

    assert report.result_profile.headline is FactualHeadline.MULTIPLE_CHAIN_PROJECTIONS
    assert _consumed_roles(report) == {UCSCBundleResourceRole.CHAIN}


def test_liftover_only_zero_candidate_report_marks_only_chain_consumed(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=1, start=200),
    )

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "upstream UCSC alignment"),
    )

    assert report.result_profile.headline is FactualHeadline.NO_CHAIN_PROJECTION
    assert report.candidates == ()
    assert _consumed_roles(report) == {UCSCBundleResourceRole.CHAIN}


def test_zero_candidate_progress_never_consumes_comparative_evidence_resources(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _comparative_bundle(tmp_path, chain_start=200)
    events: list[tuple[UCSCBundleResourceRole, int, int]] = []

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
        progress_callback=lambda role, read, total: events.append((role, read, total)),
    )

    assert report.candidates == ()
    assert {role for role, _, _ in events} == {UCSCBundleResourceRole.CHAIN}
    assert events[-1] == (
        UCSCBundleResourceRole.CHAIN,
        bundle.chain.size_bytes,
        bundle.chain.size_bytes,
    )


def test_cached_bundle_progress_reports_exact_consumed_raw_byte_totals(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _comparative_bundle(tmp_path)
    events: list[tuple[UCSCBundleResourceRole, int, int]] = []

    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 115),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
        progress_callback=lambda role, read, total: events.append((role, read, total)),
    )

    assert report.candidates
    assert bundle.net is not None
    assert bundle.reciprocal_best_chain is not None
    final_by_role = {role: (read, total) for role, read, total in events}
    assert final_by_role == {
        UCSCBundleResourceRole.CHAIN: (
            bundle.chain.size_bytes,
            bundle.chain.size_bytes,
        ),
        UCSCBundleResourceRole.NET: (bundle.net.size_bytes, bundle.net.size_bytes),
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN: (
            bundle.reciprocal_best_chain.size_bytes,
            bundle.reciprocal_best_chain.size_bytes,
        ),
    }
    assert all(0 < read <= total for _, read, total in events)


def test_indexed_point_context_maps_clean_window_without_reusing_comparative_evidence(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=0, length=300, target_start=500),
    )
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    context = enriched.result_profile.query_context
    assert context.check_state is QueryContextState.RUN
    assert context.tested_source_interval == GenomicInterval(source, "chr1", 50, 151)
    assert context.actual_window_bases == 101
    assert context.findings == (QueryContextFinding.AGREES_WITH_POINT,)
    assert context.point_and_local_context_map_together
    assert len(context.candidate_profiles) == 1
    assert enriched.query_context_result is not None
    assert {
        observation.kind
        for observation in enriched.query_context_result.candidates[0].evidence
    } == {
        EvidenceKind.CHAIN_SCORE,
        EvidenceKind.MAPPING_COVERAGE,
        EvidenceKind.CHAIN_GAPS,
    }


def test_indexed_point_context_with_no_projection_at_either_scale_is_explicit(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=0, length=300, target_start=500),
    )
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 500, 501),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    context = enriched.result_profile.query_context
    assert report.candidates == ()
    assert context.check_state is QueryContextState.RUN
    assert context.candidate_profiles == ()
    assert context.findings == (QueryContextFinding.NO_PROJECTION_AT_EITHER_SCALE,)
    assert not context.point_and_local_context_map_together
    summary = render_assessment_summary(enriched)
    assert "no chain projection was found for the point" in summary
    assert "agrees with the point-level chain result" not in summary


def test_point_context_without_index_is_explicitly_not_run(tmp_path: Path) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=0, length=300, target_start=500),
    )
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=None,
    )

    context = enriched.result_profile.query_context
    assert context.check_state is QueryContextState.NOT_RUN
    assert context.not_run_reason is QueryContextNotRunReason.INDEX_UNAVAILABLE
    assert enriched.query_context_result is not None
    assert enriched.query_context_result.candidates == ()


def test_valid_source_preflight_is_preserved_in_completed_report(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    interval, preflight, metadata_resource = _source_preflight(
        tmp_path,
        source,
        sequence_name="chr1",
        sequence_length=1000,
        start=100,
        end=120,
    )
    bundle = _liftover_bundle(tmp_path, _chain_text(chain_id=17))

    report = assess_ucsc_cached_bundle(
        interval,
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
        source_preflight=preflight,
        source_preflight_resources=(metadata_resource,),
    )

    assert report.source_preflight == preflight
    assert report.source_preflight_resources == (metadata_resource,)
    assert report.result_profile.input_validity is InputValidityState.VALID


def test_authoritative_bounds_allow_context_for_valid_sequence_absent_from_chain(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    interval, preflight, metadata_resource = _source_preflight(
        tmp_path,
        source,
        sequence_name="chrNoChain",
        sequence_length=200,
        start=100,
        end=101,
    )
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=0, length=200, target_start=500),
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index
    report = assess_ucsc_cached_bundle(
        interval,
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
        source_preflight=preflight,
        source_preflight_resources=(metadata_resource,),
        chain_index=index,
    )

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    context = enriched.result_profile.query_context
    assert report.candidates == ()
    assert context.check_state is QueryContextState.RUN
    assert context.not_run_reason is None
    assert context.findings == (QueryContextFinding.NO_PROJECTION_AT_EITHER_SCALE,)
    assert enriched.query_context_result is not None
    assert enriched.query_context_result.tested_source_interval == GenomicInterval(
        source, "chrNoChain", 50, 151
    )


def test_point_context_without_indexed_source_bound_is_explicitly_not_run(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=0, length=300, target_start=500),
    )
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chrMissing", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    context = enriched.result_profile.query_context
    assert context.check_state is QueryContextState.NOT_RUN
    assert context.not_run_reason is QueryContextNotRunReason.SOURCE_BOUNDS_UNAVAILABLE
    assert enriched.query_context_result is not None
    assert enriched.query_context_result.candidates == ()


def test_point_context_reports_partial_local_geometry_as_scale_change(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=90, length=21, target_start=500),
    )
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    context = enriched.result_profile.query_context
    assert context.check_state is QueryContextState.RUN
    assert set(context.findings) == {
        QueryContextFinding.REVEALS_PARTIAL_COVERAGE,
        QueryContextFinding.CHANGES_WITH_QUERY_SCALE,
    }
    assert not context.point_and_local_context_map_together
    assert context.maximum_candidate_covered_source_bases == 21
    assert context.actual_window_bases == 101


def test_point_context_keeps_structural_findings_separate(tmp_path: Path) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        ("chain 100 chr1 1000 + 50 151 chrA 2000 + 500 610 17\n70 1 10\n30\n\n"),
    )
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    context = enriched.result_profile.query_context
    assert set(context.findings) == {
        QueryContextFinding.REVEALS_PARTIAL_COVERAGE,
        QueryContextFinding.REVEALS_FRAGMENTATION,
        QueryContextFinding.REVEALS_TARGET_DISCONTINUITY,
        QueryContextFinding.CHANGES_WITH_QUERY_SCALE,
    }


def test_comparative_point_context_uses_all_chain_geometry_without_net_rbest_repass(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _comparative_bundle(
        tmp_path,
        chain_text=_chain_text(chain_id=1, start=0, length=300, target_start=500),
        reciprocal_chain_text=_chain_text(
            chain_id=101,
            start=0,
            length=300,
            target_start=500,
        ),
    )
    alignment = ProvenanceSource("alignment", "shared UCSC alignment")
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 105, 106),
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )
    assert any(
        observation.kind is EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
        for observation in report.candidates[0].evidence
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    assert enriched.result_profile.query_context.check_state is QueryContextState.RUN
    assert enriched.query_context_result is not None
    context_evidence_kinds = {
        observation.kind
        for observation in enriched.query_context_result.candidates[0].evidence
    }
    assert EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP not in context_evidence_kinds
    assert EvidenceKind.NET_CLASSIFICATION not in context_evidence_kinds


def test_query_context_rejects_non_chain_evidence_kind(tmp_path: Path) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=0, length=300, target_start=500),
    )
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index
    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )
    assert enriched.query_context_result is not None
    context_candidate = enriched.query_context_result.candidates[0]
    invalid_candidate = replace(
        context_candidate,
        evidence=context_candidate.evidence
        + (
            EvidenceObservation(
                observation_id="unexpected-net-classification",
                kind=EvidenceKind.NET_CLASSIFICATION,
                value="top",
                provenance=context_candidate.mapping_provenance,
            ),
        ),
    )
    invalid_context = replace(
        enriched.query_context_result,
        candidates=(invalid_candidate,),
    )

    with pytest.raises(ValueError, match="only forward-chain evidence"):
        attach_query_context_result(report, invalid_context)


def test_query_context_rejects_shared_candidate_id_with_shifted_point_geometry(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        _chain_text(chain_id=17, start=0, length=300, target_start=500),
    )
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index
    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )
    assert enriched.query_context_result is not None
    context_candidate = enriched.query_context_result.candidates[0]

    shifted_segments = tuple(
        MappingSegment(
            source_interval=segment.source_interval,
            target_interval=GenomicInterval(
                assembly=segment.target_interval.assembly,
                sequence_name=segment.target_interval.sequence_name,
                start=segment.target_interval.start + 1,
                end=segment.target_interval.end + 1,
            ),
        )
        for segment in context_candidate.segments
    )
    shifted_candidate = replace(
        context_candidate,
        target_interval=GenomicInterval(
            assembly=context_candidate.target_interval.assembly,
            sequence_name=context_candidate.target_interval.sequence_name,
            start=context_candidate.target_interval.start + 1,
            end=context_candidate.target_interval.end + 1,
        ),
        segments=shifted_segments,
    )
    mismatched_context = replace(
        enriched.query_context_result,
        candidates=(shifted_candidate,),
    )

    with pytest.raises(ValueError, match="reproduce the point mapping"):
        attach_query_context_result(report, mismatched_context)


def test_reverse_orientation_point_context_reproduces_point_geometry(
    tmp_path: Path,
) -> None:
    source, target = _assemblies()
    bundle = _liftover_bundle(
        tmp_path,
        ("chain 100 chr1 1000 + 0 300 chrA 2000 - 500 800 17\n300\n\n"),
    )
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 101),
        bundle,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("alignment", "shared UCSC alignment"),
    )
    index = build_cached_chain_index(tmp_path / "cache", bundle.chain).index

    enriched = attach_point_query_context(
        report,
        chain_context=CachedUCSCChainResource(
            source_db=bundle.source_db,
            target_db=bundle.target_db,
            evidence_tier=bundle.evidence_tier,
            chain=bundle.chain,
        ),
        chain_index=index,
    )

    assert report.candidates[0].target_interval == GenomicInterval(
        target, "chrA", 1399, 1400
    )
    assert enriched.query_context_result is not None
    assert enriched.query_context_result.candidates[
        0
    ].target_interval == GenomicInterval(target, "chrA", 1349, 1450)
    context = enriched.result_profile.query_context
    assert context.findings == (QueryContextFinding.AGREES_WITH_POINT,)
    assert context.point_and_local_context_map_together
