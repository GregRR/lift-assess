from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from liftassess import (
    AssemblyIdentifier,
    AssemblySequenceCatalog,
    AssemblySequenceMetadata,
    AssemblySequenceRoleContext,
    CachedAssemblyRoleArtifact,
    CachedResource,
    CachedTargetAssemblyRoleMetadata,
    ChainGap,
    ChainGapSummary,
    ComparativeEvidenceRelationship,
    EvidenceAvailabilityTier,
    EvidenceKind,
    EvidenceObservation,
    FactualHeadline,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    MappingSegment,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceIdentifier,
    ProvenanceIdentifierKind,
    ProvenanceSource,
    ProviderChecksum,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
    ResourceChecksumAlgorithm,
    TargetRoleState,
    UCSCAssessmentReport,
    UCSCAssessmentResource,
    UCSCBundleResourceRole,
    build_comparative_evidence_relationship,
    build_filtered_all_chain_comparison,
    build_result_profile,
    reporting,
    ucsc_resource_terms,
)
from liftassess.chain import chain_candidate_id
from liftassess.reporting import (
    format_display_interval,
    render_assessment_details,
    render_assessment_summary,
)

SOURCE_ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")
TARGET_ASSEMBLY = AssemblyIdentifier("targetAsm", "test")
SOURCE = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 200)
ALIGNMENT = ProvenanceSource("alignment", "shared test alignment")
CHAIN = ProvenanceSource(
    "chain-file",
    "chain resource",
    identifiers=(
        ProvenanceIdentifier(ProvenanceIdentifierKind.SHA256, "sha256:" + "a" * 64),
    ),
    derived_from=(ALIGNMENT,),
)
NET = ProvenanceSource(
    "net-file",
    "net resource",
    identifiers=(
        ProvenanceIdentifier(ProvenanceIdentifierKind.SHA256, "sha256:" + "b" * 64),
    ),
    derived_from=(ALIGNMENT,),
)
RBEST = ProvenanceSource(
    "rbest-file",
    "reciprocal-best resource",
    identifiers=(
        ProvenanceIdentifier(ProvenanceIdentifierKind.SHA256, "sha256:" + "c" * 64),
    ),
    derived_from=(ALIGNMENT,),
)
FILTERED_CHAIN = ProvenanceSource(
    "filtered-chain-file",
    "ordinary filtered liftOver chain resource",
    identifiers=(
        ProvenanceIdentifier(ProvenanceIdentifierKind.SHA256, "sha256:" + "f" * 64),
    ),
    derived_from=(ALIGNMENT,),
)


def _candidate(
    chain_id: int,
    *,
    source_spans: tuple[tuple[int, int], ...] = ((100, 200),),
    target_spans: tuple[tuple[int, int], ...] = ((1000, 1100),),
    orientation: MappingOrientation = MappingOrientation.SAME,
    target_gaps: tuple[tuple[int, int], ...] = (),
    reciprocal_best: ReciprocalBestMembershipStatus | None = None,
    extra_evidence: tuple[EvidenceObservation, ...] = (),
) -> NormalizedCandidate:
    candidate_id = chain_candidate_id(CHAIN.source_id, chain_id)
    segments = tuple(
        MappingSegment(
            GenomicInterval(SOURCE_ASSEMBLY, "chr1", source_start, source_end),
            GenomicInterval(TARGET_ASSEMBLY, "chrA", target_start, target_end),
        )
        for (source_start, source_end), (target_start, target_end) in zip(
            source_spans, target_spans, strict=True
        )
    )
    covered = sum(end - start for start, end in source_spans)
    uncovered: list[GenomicInterval] = []
    cursor = SOURCE.start
    for start, end in source_spans:
        if cursor < start:
            uncovered.append(GenomicInterval(SOURCE_ASSEMBLY, "chr1", cursor, start))
        cursor = end
    if cursor < SOURCE.end:
        uncovered.append(GenomicInterval(SOURCE_ASSEMBLY, "chr1", cursor, SOURCE.end))

    evidence: list[EvidenceObservation] = [
        EvidenceObservation(
            f"{candidate_id}:coverage",
            EvidenceKind.MAPPING_COVERAGE,
            MappingCoverageSummary(
                status=(
                    MappingCoverageStatus.FULL
                    if covered == SOURCE.length
                    else MappingCoverageStatus.PARTIAL
                ),
                covered_source_bases=covered,
                source_bases=SOURCE.length,
                uncovered_source_intervals=tuple(uncovered),
            ),
            CHAIN,
        ),
        EvidenceObservation(
            f"{candidate_id}:gaps",
            EvidenceKind.CHAIN_GAPS,
            ChainGapSummary(
                tuple(
                    ChainGap(
                        source_boundary=source_spans[
                            min(index + 1, len(source_spans) - 1)
                        ][0],
                        target_gap_interval=GenomicInterval(
                            TARGET_ASSEMBLY, "chrA", gap_start, gap_end
                        ),
                    )
                    for index, (gap_start, gap_end) in enumerate(target_gaps)
                )
            ),
            CHAIN,
        ),
    ]
    evidence.extend(extra_evidence)
    if reciprocal_best is not None:
        if reciprocal_best is ReciprocalBestMembershipStatus.FULL:
            rbest_covered = covered
            intervals = tuple(segment.source_interval for segment in segments)
        elif reciprocal_best is ReciprocalBestMembershipStatus.PARTIAL:
            rbest_covered = covered - 1
            intervals = (
                GenomicInterval(
                    SOURCE_ASSEMBLY,
                    "chr1",
                    source_spans[0][0],
                    source_spans[0][1] - 1,
                ),
            )
        else:
            rbest_covered = 0
            intervals = ()
        evidence.append(
            EvidenceObservation(
                f"{candidate_id}:rbest",
                EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP,
                ReciprocalBestMembershipSummary(
                    status=reciprocal_best,
                    resource_completeness=(
                        ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
                    ),
                    chains_examined=2,
                    covered_source_bases=rbest_covered,
                    candidate_source_bases=covered,
                    covered_source_intervals=intervals,
                ),
                RBEST,
            )
        )

    return NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=GenomicInterval(
            TARGET_ASSEMBLY,
            "chrA",
            min(start for start, _ in target_spans),
            max(end for _, end in target_spans),
        ),
        orientation=orientation,
        mapping_provenance=CHAIN,
        segments=segments,
        evidence=tuple(evidence),
    )


def _cached_resource(
    role: UCSCBundleResourceRole,
    *,
    digest_char: str,
) -> CachedResource:
    return CachedResource(
        path=Path(f"/cache/{role.value.lower()}.gz"),
        source_url=f"https://example.test/{role.value.lower()}.gz",
        retrieved_at="2026-08-19T00:00:00Z",
        sha256="sha256:" + digest_char * 64,
        size_bytes=123,
        provider_checksum=None,
        terms=ucsc_resource_terms(
            "https://hgdownload.soe.ucsc.edu/goldenPath/sourceAsm/liftOver/"
            "sourceAsmToTargetAsm.over.chain.gz"
        ),
        cache_hit=True,
    )


def _resource(
    role: UCSCBundleResourceRole,
    *,
    consumed: bool,
    provenance: ProvenanceSource | None,
    digest_char: str,
) -> UCSCAssessmentResource:
    return UCSCAssessmentResource(
        role=role,
        resource=_cached_resource(role, digest_char=digest_char),
        consumed_by_engine=consumed,
        file_provenance=provenance,
    )


def _report(
    candidates: tuple[NormalizedCandidate, ...],
    *,
    tier: EvidenceAvailabilityTier = EvidenceAvailabilityTier.LIFTOVER_ONLY,
) -> UCSCAssessmentReport:
    resources: tuple[UCSCAssessmentResource, ...]
    if tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        resources = (
            _resource(
                UCSCBundleResourceRole.CHAIN,
                consumed=True,
                provenance=CHAIN,
                digest_char="a",
            ),
        )
    else:
        resources = (
            _resource(
                UCSCBundleResourceRole.CHAIN,
                consumed=True,
                provenance=CHAIN,
                digest_char="a",
            ),
            _resource(
                UCSCBundleResourceRole.NET,
                consumed=bool(candidates),
                provenance=NET if candidates else None,
                digest_char="b",
            ),
            _resource(
                UCSCBundleResourceRole.SYNTENIC_NET,
                consumed=False,
                provenance=None,
                digest_char="d",
            ),
            _resource(
                UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
                consumed=bool(candidates),
                provenance=RBEST if candidates else None,
                digest_char="c",
            ),
            _resource(
                UCSCBundleResourceRole.RECIPROCAL_BEST_NET,
                consumed=False,
                provenance=None,
                digest_char="e",
            ),
        )
    consumed_roles = tuple(
        resource.role.value for resource in resources if resource.consumed_by_engine
    )
    profile = build_result_profile(
        SOURCE,
        candidates,
        evidence_tier=tier,
        consumed_resource_roles=consumed_roles,
    )
    return UCSCAssessmentReport(
        source_interval=SOURCE,
        target_assembly=TARGET_ASSEMBLY,
        candidates=candidates,
        evidence_tier=tier,
        result_profile=profile,
        source_db="sourceAsm",
        target_db="targetAsm",
        alignment_provenance=ALIGNMENT,
        resources=resources,
    )


def _with_depth1_top_net(candidate: NormalizedCandidate) -> NormalizedCandidate:
    fill = ProvenanceSource(
        f"{candidate.candidate_id}:top-fill",
        "depth-1 top net fill",
        derived_from=(NET,),
    )
    return replace(
        candidate,
        evidence=candidate.evidence
        + (
            EvidenceObservation(
                f"{candidate.candidate_id}:net-classification",
                EvidenceKind.NET_CLASSIFICATION,
                "top",
                fill,
            ),
            EvidenceObservation(
                f"{candidate.candidate_id}:net-hierarchy",
                EvidenceKind.NET_HIERARCHY,
                NetHierarchySummary(depth=1, source_fill_interval=SOURCE),
                fill,
            ),
        ),
    )


def _filtered_candidate(candidate: NormalizedCandidate) -> NormalizedCandidate:
    evidence = tuple(
        replace(observation, provenance=FILTERED_CHAIN)
        for observation in candidate.evidence
        if observation.kind in {EvidenceKind.MAPPING_COVERAGE, EvidenceKind.CHAIN_GAPS}
    )
    return replace(
        candidate,
        candidate_id=f"filtered:{candidate.candidate_id}",
        mapping_provenance=FILTERED_CHAIN,
        evidence=evidence,
    )


def _with_filtered_all_chain_comparison(
    report: UCSCAssessmentReport,
    filtered_candidates: tuple[NormalizedCandidate, ...],
) -> UCSCAssessmentReport:
    comparison = build_filtered_all_chain_comparison(
        SOURCE,
        report.candidates,
        filtered_candidates,
        all_chain_provenance=CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )
    relationship = build_comparative_evidence_relationship(comparison)
    profile = build_result_profile(
        SOURCE,
        report.candidates,
        evidence_tier=report.evidence_tier,
        consumed_resource_roles=report.result_profile.consumed_resource_roles,
        filtered_all_chain_comparison=comparison,
        comparative_evidence_relationship=relationship,
    )
    comparison_resource = _resource(
        UCSCBundleResourceRole.CHAIN,
        consumed=True,
        provenance=FILTERED_CHAIN,
        digest_char="f",
    )
    return replace(
        report,
        result_profile=profile,
        filtered_all_chain_comparison=comparison,
        comparative_evidence_relationship=relationship,
        filtered_chain_comparison_resource=comparison_resource,
    )


def test_format_display_interval_makes_coordinate_convention_explicit() -> None:
    interval = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 12344, 12400)

    assert format_display_interval(interval) == "chr1:12345-12400 (1-based inclusive)"


def test_clean_default_summary_is_compact_facts_first_and_verdict_free() -> None:
    summary = render_assessment_summary(_report((_candidate(42),)))

    assert summary.splitlines()[:7] == [
        "* ONE COMPLETE CHAIN PROJECTION *",
        "Source:",
        "    chr1:101-200 (1-based inclusive)",
        "Source coverage:",
        "    100/100 source bases",
        "Target:",
        "    chrA:1001-1100 (1-based inclusive; same orientation)",
    ]
    assert "Assessment:" not in summary
    assert "WELL SUPPORTED" not in summary
    assert "Preferred candidate" not in summary
    assert "LIFTOVER-ONLY" in summary
    assert "named-variant and gene/transcript identity not assessed" in summary


def test_partial_fragmented_summary_expands_with_exact_coverage_and_gaps() -> None:
    report = _report(
        (
            _candidate(
                42,
                source_spans=((100, 150), (160, 190)),
                target_spans=((1000, 1050), (1060, 1090)),
                target_gaps=((1050, 1060),),
            ),
        )
    )

    summary = render_assessment_summary(report)

    assert summary.startswith("* PARTIAL AND FRAGMENTED PROJECTION *")
    assert "Source coverage:\n    80/100 source bases" in summary
    assert "Geometric mapped segments:\n    2" in summary
    assert "Uncovered source:\n    chr1:151-160" in summary
    assert "chr1:191-200" in summary
    assert "Target gaps:\n    chrA:1051-1060" in summary
    assert "bounding span of 2 geometric mapped segments" in summary


def test_multiple_projection_summary_leads_with_coverage_before_count() -> None:
    report = _report(
        (
            _candidate(
                1,
                source_spans=((100, 160),),
                target_spans=((1000, 1060),),
            ),
            _candidate(
                2,
                source_spans=((160, 200),),
                target_spans=((3000, 3040),),
            ),
        )
    )

    lines = render_assessment_summary(report).splitlines()

    assert lines[0] == "* SOURCE INTERVAL SPLITS ACROSS MULTIPLE PROJECTIONS *"
    assert lines.index("Maximum candidate source coverage:") < lines.index(
        "Chain projections:"
    )
    assert "    60/100 bases" in lines
    assert "Source bases represented across all projections:" in lines
    assert "    100/100" in lines
    assert "Projection order:" in lines
    assert "    reproducibility only; not rank." in lines


def test_large_multiple_projection_summary_is_bounded_without_candidate_sampling() -> (
    None
):
    report = _report(
        tuple(
            _candidate(
                chain_id,
                target_spans=((1000 + chain_id * 200, 1100 + chain_id * 200),),
            )
            for chain_id in range(1, 6)
        )
    )

    summary = render_assessment_summary(report)
    lines = summary.splitlines()

    assert "Chain projections:" in lines
    assert "    5" in lines
    assert "Target sequences represented:" in lines
    assert "    1" in lines
    assert "Projection orientations:" in lines
    assert "    SAME" in lines
    assert "Geometric mapped segments per projection:" in lines
    assert "Projections at maximum source coverage:" in lines
    assert not any(line.startswith("    - ") for line in lines)
    assert (
        "    omitted from default output for this candidate set; use --details or "
        "--json for every projection." in lines
    )


def test_summary_distinguishes_exact_blocks_from_geometric_fragmentation() -> None:
    candidate = _candidate(
        41,
        source_spans=((100, 150), (150, 200)),
        target_spans=((1000, 1050), (1050, 1100)),
    )
    report = _report((candidate,))

    summary = render_assessment_summary(report)
    details = render_assessment_details(report)

    candidate_profile = report.result_profile.candidate_profiles[0]
    assert (
        report.result_profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION
    )
    assert candidate_profile.exact_mapped_segment_count == 2
    assert candidate_profile.geometric_segment_count == 1
    assert "Geometric mapped segments: 2" not in summary
    assert "bounding span of 2" not in summary
    assert "Geometric mapped segments: 1" in details
    assert "Exact chain-derived mapped segments (2):" in details


def test_comparative_summary_names_consumed_resources_and_dependency_boundary() -> None:
    report = _report(
        (_candidate(42, reciprocal_best=ReciprocalBestMembershipStatus.FULL),),
        tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    summary = render_assessment_summary(report)

    assert "COMPARATIVE" in summary
    assert "CHAIN, NET, RECIPROCAL_BEST_CHAIN" in summary
    assert "not independent votes" in summary
    assert "exact shared processing-run provenance is not verified" in summary
    assert (
        report.result_profile.headline is FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION
    )


def test_details_exposes_profile_evidence_resources_and_scope() -> None:
    report = _report((_candidate(42),))

    details = render_assessment_details(report)

    assert "Detailed factual result dossier" in details
    assert "Headline: ONE COMPLETE CHAIN PROJECTION" in details
    assert "Actual reverse mapping: NOT_RUN" in details
    assert "Point/neighborhood context: NOT_RUN" in details
    assert "Comparative relationship synthesis: NOT_ASSESSED" in details
    assert "Chain 42" in details
    assert "MAPPING_COVERAGE: FULL; 100/100 source bases covered" in details
    assert "CHAIN_GAPS: 0 chain gap(s)" in details
    assert "CHAIN [consumed]" in details
    assert "Candidate order is preserved for reproducibility" in details
    assert "verdict" not in details.lower()
    assert details.endswith("This does not establish biological correctness.")


def test_json_schema_v2_uses_result_profile_and_removes_legacy_verdict_fields() -> None:
    report = _report((_candidate(42),))

    payload = json.loads(reporting.render_assessment_json(report))

    assert payload["schema_version"] == 2
    assert payload["report_type"] == "liftassess.ucsc_result"
    assert "aggregate_verdict" not in payload["semantics"]
    assert "assessment" not in payload
    assert "verdict" not in payload["result_profile"]
    assert "decision_reason" not in payload["result_profile"]
    assert "preferred_candidate_id" not in payload["result_profile"]
    assert payload["result_profile"]["headline"] == "ONE_COMPLETE_CHAIN_PROJECTION"
    assert payload["result_profile"]["source_coverage"] == {
        "maximum_candidate_covered_source_bases": 100,
        "maximum_coverage_candidate_ids": [chain_candidate_id(CHAIN.source_id, 42)],
        "source_bases": 100,
        "state": "COMPLETE",
        "union_covered_source_bases": 100,
    }
    assert payload["result_profile"]["scope"]["actual_reverse_mapping"] == "NOT_RUN"
    assert payload["result_profile"]["comparative_relationship"]["state"] == (
        "NOT_ASSESSED"
    )
    assert payload["filtered_all_chain_comparison"] == {"assessed": False}
    assert payload["semantics"]["comparative_relationships"] == (
        "categorical_not_scores_or_votes"
    )
    assert (
        payload["candidates"][0]["target_bounding_interval"]["coordinate_system"]
        == "0-based-half-open"
    )
    assert "assessment_role" not in payload["candidates"][0]["evidence"][0]


def test_json_and_details_preserve_reverse_orientation_coordinates() -> None:
    report = _report(
        (
            _candidate(
                42,
                orientation=MappingOrientation.REVERSE,
                target_spans=((1000, 1100),),
            ),
        )
    )

    payload = json.loads(reporting.render_assessment_json(report))
    candidate = payload["candidates"][0]

    assert candidate["orientation"] == "REVERSE"
    assert candidate["target_bounding_interval"]["start"] == 1000
    assert candidate["target_bounding_interval"]["end"] == 1100
    assert "reverse orientation" in render_assessment_details(report)


def test_comparative_json_keeps_exact_observations() -> None:
    net_observation = EvidenceObservation(
        "net",
        EvidenceKind.NET_HIERARCHY,
        NetHierarchySummary(
            depth=3,
            source_fill_interval=GenomicInterval(SOURCE_ASSEMBLY, "chr1", 90, 210),
        ),
        NET,
    )
    report = _report(
        (
            _candidate(
                42,
                reciprocal_best=ReciprocalBestMembershipStatus.FULL,
                extra_evidence=(net_observation,),
            ),
        ),
        tier=EvidenceAvailabilityTier.COMPARATIVE,
    )

    payload = json.loads(reporting.render_assessment_json(report))
    evidence_by_kind = {
        item["kind"]: item for item in payload["candidates"][0]["evidence"]
    }

    assert evidence_by_kind["NET_HIERARCHY"]["value"]["depth"] == 3
    assert evidence_by_kind["RECIPROCAL_BEST_MEMBERSHIP"]["value"]["status"] == "FULL"
    assert payload["provenance"]["sources"]


def test_comparative_summary_explains_why_one_placement_is_favored() -> None:
    favored = _with_depth1_top_net(
        _candidate(
            41,
            target_spans=((1000, 1100),),
            reciprocal_best=ReciprocalBestMembershipStatus.FULL,
        )
    )
    competitor = _candidate(
        42,
        target_spans=((2000, 2100),),
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )
    report = _with_filtered_all_chain_comparison(
        _report((favored, competitor), tier=EvidenceAvailabilityTier.COMPARATIVE),
        (_filtered_candidate(favored),),
    )

    summary = render_assessment_summary(report)

    assert "all-chain reveals 1 additional placement" in summary
    assert "available categorical evidence favors one placement" in summary
    assert "    Favored placement:\n        chrA:1001-1100" in summary
    assert (
        "only complete placement retained by the ordinary filtered liftOver chain"
        in summary
    )
    assert "depth-1 top-net support plus full reciprocal-best membership" in summary
    assert "no competing complete placement has that same joint support" in summary
    assert "not independent votes" in summary


def test_mixed_comparative_summary_names_the_conflicting_placements() -> None:
    filtered_only = _candidate(
        51,
        target_spans=((1000, 1100),),
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )
    top_rbest = _with_depth1_top_net(
        _candidate(
            52,
            target_spans=((3000, 3100),),
            reciprocal_best=ReciprocalBestMembershipStatus.FULL,
        )
    )
    report = _with_filtered_all_chain_comparison(
        _report((filtered_only, top_rbest), tier=EvidenceAvailabilityTier.COMPARATIVE),
        (_filtered_candidate(filtered_only),),
    )

    summary = render_assessment_summary(report)

    assert (
        report.comparative_evidence_relationship is not None
        and report.comparative_evidence_relationship.relationship
        is ComparativeEvidenceRelationship.MIXED_CONFLICTING
    )
    assert "available categorical evidence is mixed/conflicting" in summary
    assert (
        "    Complete placements retained by filtered chain:\n        chrA:1001-1100"
        in summary
    )
    assert (
        "    Complete placements with depth-1 top-net support:\n        chrA:3001-3100"
        in summary
    )
    assert (
        "    Complete placements with full reciprocal-best membership:\n"
        "        chrA:3001-3100" in summary
    )


def test_nonseparating_comparative_summary_exposes_missing_category_support() -> None:
    retained = _candidate(
        61,
        target_spans=((1000, 1100),),
        reciprocal_best=ReciprocalBestMembershipStatus.FULL,
    )
    competitor = _candidate(
        62,
        target_spans=((4000, 4100),),
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )
    report = _with_filtered_all_chain_comparison(
        _report((retained, competitor), tier=EvidenceAvailabilityTier.COMPARATIVE),
        (_filtered_candidate(retained),),
    )

    summary = render_assessment_summary(report)

    assert (
        report.comparative_evidence_relationship is not None
        and report.comparative_evidence_relationship.relationship
        is ComparativeEvidenceRelationship.DOES_NOT_SEPARATE_PLACEMENTS
    )
    assert "does not separate the complete placements" in summary
    assert (
        "    Complete placements retained by filtered chain:\n        chrA:1001-1100"
        in summary
    )
    assert (
        "    Complete placements with depth-1 top-net support:\n        none" in summary
    )
    assert (
        "    Complete placements with full reciprocal-best membership:\n"
        "        chrA:1001-1100" in summary
    )


def test_comparative_details_and_json_expose_inventory_support_and_provenance() -> None:
    favored = _with_depth1_top_net(
        _candidate(
            71,
            target_spans=((1000, 1100),),
            reciprocal_best=ReciprocalBestMembershipStatus.FULL,
        )
    )
    competitor = _candidate(
        72,
        target_spans=((5000, 5100),),
        reciprocal_best=ReciprocalBestMembershipStatus.NONE,
    )
    report = _with_filtered_all_chain_comparison(
        _report((favored, competitor), tier=EvidenceAvailabilityTier.COMPARATIVE),
        (_filtered_candidate(favored),),
    )

    details = render_assessment_details(report)
    payload = json.loads(reporting.render_assessment_json(report))

    assert "Filtered/all-chain comparative relationship" in details
    assert "Inventory state: ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS" in details
    assert "Categorical relationship: FAVORS_ONE_PLACEMENT" in details
    assert "retained by filtered chain=yes" in details
    assert "depth-1 top-net=yes" in details
    assert "full reciprocal-best=yes" in details
    assert "Filtered-chain comparison resource" in details
    assert "consumed for paired comparison" in details
    assert "UCSC pair dependency group: alignment" in details
    assert "Exact shared processing-run provenance: not verified" in details
    assert "Derived from: alignment" in details

    profile = payload["result_profile"]["comparative_relationship"]
    assert profile["state"] == "FAVORS_ONE_PLACEMENT"
    assert profile["favored_candidate_id"] == favored.candidate_id
    assert profile["placement_support"][0] == {
        "candidate_id": favored.candidate_id,
        "complete_source_coverage": True,
        "depth1_top_net": True,
        "full_reciprocal_best": True,
        "retained_by_filtered_chain": True,
    }

    assert payload["semantics"]["ucsc_pair_dependency_group"] == (
        "conservative_grouping_not_processing_run_proof"
    )

    comparison = payload["filtered_all_chain_comparison"]
    assert comparison["assessed"] is True
    assert comparison["inventory_state"] == "ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS"
    assert comparison["categorical_relationship"] == "FAVORS_ONE_PLACEMENT"
    assert comparison["candidate_matches"] == [
        {
            "all_chain_candidate_id": favored.candidate_id,
            "filtered_candidate_id": f"filtered:{favored.candidate_id}",
        }
    ]
    assert comparison["additional_all_chain_candidate_ids"] == [competitor.candidate_id]
    assert comparison["filtered_chain_resource"]["consumed_by_engine"] is True
    assert comparison["provenance"] == {
        "all_chain_source_id": CHAIN.source_id,
        "filtered_chain_source_id": FILTERED_CHAIN.source_id,
        "shared_ucsc_pair_dependency_source_ids": [ALIGNMENT.source_id],
        "shared_processing_run_provenance_verified": False,
    }
    provenance_ids = {item["source_id"] for item in payload["provenance"]["sources"]}
    assert FILTERED_CHAIN.source_id in provenance_ids


def test_provider_checksum_text_and_resource_json_preserve_transfer_metadata() -> None:
    checksum_url = "https://example.test/md5sum.txt"
    resource = CachedResource(
        path=Path("/cache/chain.gz"),
        source_url="https://example.test/chain.gz",
        retrieved_at="2026-08-19T00:00:00Z",
        sha256="sha256:" + "f" * 64,
        size_bytes=321,
        provider_checksum=ProviderChecksum(
            algorithm=ResourceChecksumAlgorithm.MD5,
            value="c" * 32,
            source_url=checksum_url,
        ),
        terms=ucsc_resource_terms(
            "https://hgdownload.soe.ucsc.edu/goldenPath/sourceAsm/liftOver/"
            "sourceAsmToTargetAsm.over.chain.gz"
        ),
        cache_hit=True,
    )
    provenance = ProvenanceSource(
        "checksum-resource",
        "checksum resource",
        identifiers=(
            ProvenanceIdentifier(ProvenanceIdentifierKind.SHA256, resource.sha256),
        ),
    )
    assessment_resource = UCSCAssessmentResource(
        role=UCSCBundleResourceRole.CHAIN,
        resource=resource,
        consumed_by_engine=True,
        file_provenance=provenance,
    )

    assert reporting._provider_checksum_text(resource) == (
        f"md5:{'c' * 32} (from {checksum_url})"
    )
    resource_json = reporting._assessment_resource_json(assessment_resource)
    assert resource_json["provider_checksum"] == {
        "algorithm": "md5",
        "source_url": checksum_url,
        "value": "c" * 32,
    }


def _target_role_catalog_for_reporting(
    *, provider_role: str
) -> AssemblySequenceCatalog:
    role_provenance = ProvenanceSource(
        "target-role-report", "version-matched NCBI sequence report"
    )
    return AssemblySequenceCatalog(
        assembly=TARGET_ASSEMBLY,
        sequences=(
            AssemblySequenceMetadata(
                sequence_name="chrA",
                length=10_000,
                role_context=AssemblySequenceRoleContext(
                    assembly_accession="GCA_000000001.1",
                    assembly_unit="Primary Assembly",
                    provider_role=provider_role,
                    length=10_000,
                    ucsc_style_name="chrA",
                ),
            ),
        ),
        sequence_provenance=ProvenanceSource(
            "target-chrom-info", "target chromInfo metadata"
        ),
        role_provenance=role_provenance,
    )


def _target_role_metadata_for_reporting() -> CachedTargetAssemblyRoleMetadata:
    accession = "GCA_000000001.1"
    description = CachedAssemblyRoleArtifact(
        path=Path("/cache/description.html"),
        source_url=(
            "https://hgdownload.soe.ucsc.edu/gbdb/targetAsm/html/description.html"
        ),
        retrieved_at="2026-08-28T00:00:00Z",
        sha256="sha256:" + "1" * 64,
        size_bytes=100,
        cache_hit=True,
    )
    sequence_report = CachedAssemblyRoleArtifact(
        path=Path("/cache/sequence_report.jsonl"),
        source_url=(
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{accession}/download?include_annotation_type=SEQUENCE_REPORT&"
            "hydrated=FULLY_HYDRATED"
        ),
        retrieved_at="2026-08-28T00:00:00Z",
        sha256="sha256:" + "2" * 64,
        size_bytes=200,
        cache_hit=True,
        archive_member=f"ncbi_dataset/data/{accession}/sequence_report.jsonl",
    )
    return CachedTargetAssemblyRoleMetadata(
        db="targetAsm",
        assembly_accession=accession,
        assembly_description=description,
        sequence_report=sequence_report,
    )


def test_reporting_marks_unavailable_target_role_without_name_inference() -> None:
    report = _report((_candidate(1),))
    profile = build_result_profile(
        SOURCE,
        report.candidates,
        evidence_tier=report.evidence_tier,
        consumed_resource_roles=report.result_profile.consumed_resource_roles,
        target_role_unavailable=True,
    )
    report = replace(report, result_profile=profile)

    summary = render_assessment_summary(report)
    payload = json.loads(reporting.render_assessment_json(report))

    assert "unavailable; no role was inferred from sequence naming" in summary
    assert payload["result_profile"]["target_role"]["state"] == "UNAVAILABLE"
    assert payload["target_role_metadata"]["assembly_accession"] is None
    assert payload["target_role_metadata"]["resources"] == []


def test_reporting_preserves_unusual_provider_target_role_and_provenance() -> None:
    report = _report((_candidate(1),))
    catalog = _target_role_catalog_for_reporting(provider_role="unplaced-scaffold")
    metadata = _target_role_metadata_for_reporting()
    profile = build_result_profile(
        SOURCE,
        report.candidates,
        evidence_tier=report.evidence_tier,
        consumed_resource_roles=report.result_profile.consumed_resource_roles,
        target_role_catalog=catalog,
    )
    report = replace(
        report,
        result_profile=profile,
        target_role_metadata=metadata,
        target_role_provenance=catalog.role_provenance,
    )

    summary = render_assessment_summary(report)
    payload = json.loads(reporting.render_assessment_json(report))

    assert "role=unplaced-scaffold" in summary
    assert profile.scope.target_role is TargetRoleState.ASSESSED
    role = payload["result_profile"]["target_role"]["sequences"][0]
    assert role["provider_role"] == "unplaced-scaffold"
    assert role["assembly_unit"] == "Primary Assembly"
    assert payload["target_role_metadata"]["assembly_accession"] == "GCA_000000001.1"
    provenance_ids = {item["source_id"] for item in payload["provenance"]["sources"]}
    assert catalog.role_provenance is not None
    assert catalog.role_provenance.source_id in provenance_ids
