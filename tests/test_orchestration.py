from __future__ import annotations

import gzip
from pathlib import Path

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCResourceBundle,
    EvidenceAvailabilityTier,
    EvidenceKind,
    FactualHeadline,
    GenomicInterval,
    ProvenanceSource,
    UCSCAssessmentReport,
    UCSCBundleResourceRole,
    assess_ucsc_cached_bundle,
    build_cached_chain_index,
    sha256_identifier_for_file,
    ucsc_resource_terms,
)


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


def _consumed_roles(report: UCSCAssessmentReport) -> set[UCSCBundleResourceRole]:
    return {
        resource.role for resource in report.resources if resource.consumed_by_engine
    }


def _assemblies() -> tuple[AssemblyIdentifier, AssemblyIdentifier]:
    return (
        AssemblyIdentifier(name="canFam3", provider="UCSC"),
        AssemblyIdentifier(name="canFam4", provider="UCSC"),
    )


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
