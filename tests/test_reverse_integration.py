from __future__ import annotations

import gzip
import json
from dataclasses import replace
from pathlib import Path

import pytest

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCResourceBundle,
    EvidenceAvailabilityTier,
    GenomicInterval,
    ProvenanceSource,
    ReverseCheckState,
    ReverseRelationshipState,
    UCSCAssessmentReport,
    assess_ucsc_cached_bundle,
    attach_reverse_mapping_context,
    sha256_identifier_for_file,
    ucsc_resource_terms,
)
from liftassess.reporting import render_assessment_details, render_assessment_json


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _resource(path: Path, url: str) -> CachedResource:
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


def _bundles(
    tmp_path: Path,
) -> tuple[CachedUCSCResourceBundle, CachedUCSCResourceBundle]:
    forward_path = tmp_path / "forward-chain"
    reverse_path = tmp_path / "reverse-chain"
    _write_gzip(
        forward_path,
        "chain 100 chr1 1000 + 100 110 chrA 2000 + 500 510 1\n10\n\n",
    )
    _write_gzip(
        reverse_path,
        "chain 100 chrA 2000 + 500 510 chr1 1000 + 100 110 2\n10\n\n",
    )
    forward_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
        "canFam3ToCanFam4.over.chain.gz"
    )
    reverse_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/"
        "canFam4ToCanFam3.over.chain.gz"
    )
    return (
        CachedUCSCResourceBundle(
            source_db="canFam3",
            target_db="canFam4",
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            chain=_resource(forward_path, forward_url),
        ),
        CachedUCSCResourceBundle(
            source_db="canFam4",
            target_db="canFam3",
            evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
            chain=_resource(reverse_path, reverse_url),
        ),
    )


def _forward_report(
    tmp_path: Path,
) -> tuple[UCSCAssessmentReport, CachedUCSCResourceBundle]:
    source = AssemblyIdentifier(name="canFam3", provider="UCSC")
    target = AssemblyIdentifier(name="canFam4", provider="UCSC")
    forward, reverse = _bundles(tmp_path)
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 100, 110),
        forward,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("forward", "forward lineage"),
    )
    return report, reverse


def test_reverse_context_populates_candidate_profile_and_resource(
    tmp_path: Path,
) -> None:
    report, reverse = _forward_report(tmp_path)
    reverse_lineage = ProvenanceSource("reverse", "reverse lineage")

    enriched = attach_reverse_mapping_context(
        report,
        reverse_bundle=reverse,
        reverse_alignment_provenance=reverse_lineage,
    )

    assert enriched.result_profile.scope.reverse_result is ReverseCheckState.RUN
    reverse_profile = enriched.result_profile.candidate_profiles[0].reverse_mapping
    assert reverse_profile.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY
    assert reverse_profile.original_source_covered_bases == 10
    assert reverse_profile.exact_original_geometry_return is True
    assert enriched.reverse_mapping_resource is not None
    assert enriched.reverse_mapping_resource.resource is reverse.chain
    assert enriched.reverse_mapping_resource.file_provenance is not None
    assert enriched.reverse_mapping_resource.file_provenance.derived_from == (
        reverse_lineage,
    )


def test_reverse_context_rejects_candidate_provenance_outside_consumed_chain(
    tmp_path: Path,
) -> None:
    report, reverse = _forward_report(tmp_path)
    enriched = attach_reverse_mapping_context(
        report,
        reverse_bundle=reverse,
        reverse_alignment_provenance=ProvenanceSource("reverse", "reverse lineage"),
    )
    assert enriched.reverse_mapping_results is not None
    result = enriched.reverse_mapping_results[0]
    segment_result = result.segment_results[0]
    reverse_candidate = segment_result.candidates[0]
    mismatched_candidate = replace(
        reverse_candidate,
        mapping_provenance=ProvenanceSource("wrong", "wrong reverse chain"),
    )
    mismatched_result = replace(
        result,
        segment_results=(replace(segment_result, candidates=(mismatched_candidate,)),),
    )

    with pytest.raises(
        ValueError,
        match="reverse candidate mapping provenance must identify",
    ):
        replace(enriched, reverse_mapping_results=(mismatched_result,))


def test_reverse_context_marks_cached_reverse_resources_unavailable(
    tmp_path: Path,
) -> None:
    report, _ = _forward_report(tmp_path)

    enriched = attach_reverse_mapping_context(report, reverse_bundle=None)

    assert enriched.result_profile.scope.reverse_result is ReverseCheckState.UNAVAILABLE
    assert enriched.reverse_mapping_resource is None
    assert enriched.reverse_alignment_provenance is None
    assert enriched.reverse_mapping_results is not None
    assert (
        enriched.reverse_mapping_results[0].check_state is ReverseCheckState.UNAVAILABLE
    )
    reverse_profile = enriched.result_profile.candidate_profiles[0].reverse_mapping
    assert reverse_profile.original_source_covered_bases is None
    assert reverse_profile.exact_original_geometry_return is None


def test_reverse_context_is_rendered_in_details_and_schema_v2_json(
    tmp_path: Path,
) -> None:
    report, reverse = _forward_report(tmp_path)
    enriched = attach_reverse_mapping_context(
        report,
        reverse_bundle=reverse,
        reverse_alignment_provenance=ProvenanceSource("reverse", "reverse lineage"),
    )

    details = render_assessment_details(enriched)
    assert "Actual reverse mapping: RUN" in details
    assert "Reverse relationship: ORIGINAL_SOURCE_ONLY" in details
    assert "Exact original aligned geometry reconstructed: yes" in details
    assert "Reverse mapping resource" in details

    payload = json.loads(render_assessment_json(enriched))
    assert payload["schema_version"] == 2
    assert payload["result_profile"]["scope"]["actual_reverse_mapping"] == "RUN"
    candidate_profile = payload["result_profile"]["candidate_profiles"][0]
    assert (
        candidate_profile["reverse_mapping"]["relationship"] == "ORIGINAL_SOURCE_ONLY"
    )
    assert (
        candidate_profile["reverse_mapping"]["exact_original_geometry_return"] is True
    )
    reverse_payload = payload["reverse_mapping"]
    assert reverse_payload["reverse_database_pair"] == {
        "source_db": "canFam4",
        "target_db": "canFam3",
    }
    assert reverse_payload["resource"]["sha256"] == reverse.chain.sha256
    assert (
        reverse_payload["candidate_results"][0]["relationship"]
        == "ORIGINAL_SOURCE_ONLY"
    )
    provenance_ids = {item["source_id"] for item in payload["provenance"]["sources"]}
    assert "reverse" in provenance_ids


def test_reverse_context_rejects_invalid_context_before_execution(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    report, reverse = _forward_report(tmp_path)
    executed = False

    def fail_if_executed(*args: object, **kwargs: object) -> tuple[object, ...]:
        nonlocal executed
        del args, kwargs
        executed = True
        raise AssertionError("reverse execution should not start before validation")

    monkeypatch.setattr(
        "liftassess.orchestration.build_reverse_mapping_results_from_cached_bundle",
        fail_if_executed,
    )
    base = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/"
    reciprocal_base = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/reciprocalBest/"
    )
    wrong_tier = CachedUCSCResourceBundle(
        source_db="canFam4",
        target_db="canFam3",
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain=_resource(reverse.chain.path, f"{base}canFam4.canFam3.all.chain.gz"),
        net=_resource(reverse.chain.path, f"{base}canFam4.canFam3.net.gz"),
        syntenic_net=_resource(reverse.chain.path, f"{base}canFam4.canFam3.syn.net.gz"),
        reciprocal_best_chain=_resource(
            reverse.chain.path,
            f"{reciprocal_base}canFam4.canFam3.rbest.chain.gz",
        ),
        reciprocal_best_net=_resource(
            reverse.chain.path,
            f"{reciprocal_base}canFam4.canFam3.rbest.net.gz",
        ),
    )

    with pytest.raises(
        ValueError,
        match="reverse bundle publication class must match",
    ):
        attach_reverse_mapping_context(
            report,
            reverse_bundle=wrong_tier,
            reverse_alignment_provenance=ProvenanceSource("reverse", "reverse lineage"),
        )

    assert executed is False


def test_reverse_context_rejects_duplicate_attach_before_execution(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    report, reverse = _forward_report(tmp_path)
    enriched = attach_reverse_mapping_context(
        report,
        reverse_bundle=reverse,
        reverse_alignment_provenance=ProvenanceSource("reverse", "reverse lineage"),
    )

    monkeypatch.setattr(
        "liftassess.orchestration.build_reverse_mapping_results_from_cached_bundle",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("duplicate attach must fail before reverse execution")
        ),
    )

    with pytest.raises(ValueError, match="reverse mapping context is already attached"):
        attach_reverse_mapping_context(
            enriched,
            reverse_bundle=reverse,
            reverse_alignment_provenance=ProvenanceSource("reverse", "reverse lineage"),
        )


def test_reverse_context_leaves_zero_candidate_report_not_run(tmp_path: Path) -> None:
    source = AssemblyIdentifier(name="canFam3", provider="UCSC")
    target = AssemblyIdentifier(name="canFam4", provider="UCSC")
    forward, reverse = _bundles(tmp_path)
    report = assess_ucsc_cached_bundle(
        GenomicInterval(source, "chr1", 200, 210),
        forward,
        target_assembly=target,
        alignment_provenance=ProvenanceSource("forward", "forward lineage"),
    )
    assert not report.candidates

    enriched = attach_reverse_mapping_context(
        report,
        reverse_bundle=reverse,
        reverse_alignment_provenance=ProvenanceSource("reverse", "reverse lineage"),
    )

    assert enriched is report
    assert enriched.reverse_mapping_results is None
    assert enriched.reverse_mapping_resource is None
    assert enriched.result_profile.scope.reverse_result is ReverseCheckState.NOT_RUN
