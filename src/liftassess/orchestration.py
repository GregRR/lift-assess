"""Scientific report orchestration for already-acquired UCSC resources.

This module connects one exact cached UCSC resource bundle to the file-backed
candidate/evidence engine and the derived factual result profile.  It deliberately
does not discover or download resources; network policy, terms acknowledgement, and
transfer confirmation remain separate CLI boundaries.

The report retains every artifact in the acquired bundle as retrieval context while
marking actual engine consumption explicitly. Cached presence is never presented as
scientific evidence unless the engine consumed the resource.
"""

from dataclasses import dataclass, replace

from .chain_index import ChainIndex
from .comparative_inventory import (
    FilteredAllChainComparisonResult,
    build_filtered_all_chain_comparison,
)
from .comparative_relationship import (
    ComparativeEvidenceRelationshipResult,
    build_comparative_evidence_relationship,
)
from .models import (
    AssemblyIdentifier,
    EvidenceAvailabilityTier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceIdentifierKind,
    ProvenanceSource,
)
from .query_context import (
    DEFAULT_POINT_CONTEXT_BASES,
    PointQueryContextResult,
    QueryContextNotRunReason,
    QueryContextState,
    build_centered_point_context_interval,
    point_context_not_run,
)
from .resource_cache import (
    CachedResource,
    CachedUCSCChainResource,
    CachedUCSCResourceBundle,
    UCSCBundleResourceRole,
)
from .resource_files import (
    ResourceReadProgressCallback,
    _cached_bundle_resource_provenance,
    _cached_chain_resource_provenance,
    build_ucsc_candidates_from_cached_bundle,
    build_ucsc_chain_candidates_for_intervals_from_cached_chain,
)
from .result_profile import (
    ComparativeRelationshipState,
    ResultProfile,
    build_result_profile,
)
from .reverse_mapping import (
    CandidateReverseMappingResult,
    ReverseCheckState,
    build_reverse_mapping_results_from_cached_bundle,
    reverse_mapping_unavailable,
)


@dataclass(frozen=True)
class UCSCAssessmentResource:
    """One cached bundle artifact and its actual role in a result run."""

    role: UCSCBundleResourceRole
    resource: CachedResource
    consumed_by_engine: bool
    file_provenance: ProvenanceSource | None

    def __post_init__(self) -> None:
        if self.consumed_by_engine != (self.file_provenance is not None):
            raise ValueError(
                "consumed UCSC assessment resources must carry file provenance and "
                "unconsumed resources must not"
            )
        if self.file_provenance is None:
            return

        sha256_identifiers = tuple(
            identifier
            for identifier in self.file_provenance.identifiers
            if identifier.kind is ProvenanceIdentifierKind.SHA256
        )
        if len(sha256_identifiers) != 1:
            raise ValueError(
                "consumed UCSC assessment resource provenance must carry exactly one "
                "SHA256 identifier"
            )
        if sha256_identifiers[0].value != self.resource.sha256:
            raise ValueError(
                "UCSC assessment resource provenance must identify the cached bytes"
            )


@dataclass(frozen=True)
class UCSCAssessmentReport:
    """Scientific candidate/evidence report plus the derived factual result profile."""

    source_interval: GenomicInterval
    target_assembly: AssemblyIdentifier
    candidates: tuple[NormalizedCandidate, ...]
    evidence_tier: EvidenceAvailabilityTier
    result_profile: ResultProfile
    source_db: str
    target_db: str
    alignment_provenance: ProvenanceSource
    resources: tuple[UCSCAssessmentResource, ...]
    query_context_result: PointQueryContextResult | None = None
    reverse_mapping_results: tuple[CandidateReverseMappingResult, ...] | None = None
    reverse_alignment_provenance: ProvenanceSource | None = None
    reverse_mapping_resource: UCSCAssessmentResource | None = None
    filtered_all_chain_comparison: FilteredAllChainComparisonResult | None = None
    comparative_evidence_relationship: ComparativeEvidenceRelationshipResult | None = (
        None
    )
    filtered_chain_comparison_resource: UCSCAssessmentResource | None = None

    def __post_init__(self) -> None:
        expected_roles = _resource_roles_for_tier(self.evidence_tier)
        actual_roles = tuple(resource.role for resource in self.resources)
        if actual_roles != expected_roles:
            raise ValueError(
                "UCSC assessment report resources must preserve the complete ordered "
                "bundle roles for the evidence tier"
            )

        expected_consumed_roles = _consumed_resource_roles(
            self.evidence_tier,
            candidates_exist=bool(self.candidates),
        )
        actual_consumed_roles = {
            resource.role for resource in self.resources if resource.consumed_by_engine
        }
        if actual_consumed_roles != expected_consumed_roles:
            raise ValueError(
                "UCSC assessment report resource-consumption metadata does not match "
                "the engine path"
            )

        if self.result_profile.source_interval != self.source_interval:
            raise ValueError("result profile source interval must match the report")
        if self.result_profile.evidence_tier is not self.evidence_tier:
            raise ValueError("result profile evidence tier must match the report")
        expected_profile_consumed_roles = tuple(
            resource.role.value
            for resource in self.resources
            if resource.consumed_by_engine
        )
        if (
            self.result_profile.consumed_resource_roles
            != expected_profile_consumed_roles
        ):
            raise ValueError(
                "result profile consumed resource roles must match the report"
            )
        profile_candidate_ids = tuple(
            profile.candidate_id for profile in self.result_profile.candidate_profiles
        )
        report_candidate_ids = tuple(
            candidate.candidate_id for candidate in self.candidates
        )
        if profile_candidate_ids != report_candidate_ids:
            raise ValueError("result profile candidates must match report candidates")

        self._validate_query_context()
        self._validate_reverse_mapping_context()
        self._validate_filtered_all_chain_comparison()

        for resource in self.resources:
            if resource.file_provenance is not None and (
                resource.file_provenance.derived_from != (self.alignment_provenance,)
            ):
                raise ValueError(
                    "consumed UCSC file provenance must derive from the report "
                    "alignment provenance"
                )

    def _validate_query_context(self) -> None:
        result = self.query_context_result
        profile = self.result_profile.query_context
        if self.result_profile.scope.query_context is not profile.check_state:
            raise ValueError("result profile query-context scope is inconsistent")

        if result is None:
            if profile.check_state is not QueryContextState.NOT_RUN:
                raise ValueError("missing query context must remain NOT_RUN")
            return

        if self.source_interval.length != 1:
            raise ValueError("query context can only be attached to a one-base report")
        if profile.check_state is not result.check_state:
            raise ValueError("result profile query context must match report context")
        if profile.requested_window_bases != result.requested_window_bases:
            raise ValueError("query context requested window must match the report")

        if result.check_state is QueryContextState.NOT_RUN:
            if profile.not_run_reason is not result.not_run_reason:
                raise ValueError("query context not-run reason must match the report")
            return

        if profile.tested_source_interval != result.tested_source_interval:
            raise ValueError("query context tested interval must match the report")
        profile_candidate_ids = tuple(
            candidate.candidate_id for candidate in profile.candidate_profiles
        )
        result_candidate_ids = tuple(
            candidate.candidate_id for candidate in result.candidates
        )
        if profile_candidate_ids != result_candidate_ids:
            raise ValueError("query context candidates must match the result profile")

        chain_resource = next(
            resource
            for resource in self.resources
            if resource.role is UCSCBundleResourceRole.CHAIN
        )
        chain_provenance = chain_resource.file_provenance
        if chain_provenance is None:
            raise ValueError("query context requires consumed forward-chain provenance")
        for candidate in result.candidates:
            if candidate.mapping_provenance != chain_provenance:
                raise ValueError(
                    "query-context candidate provenance must identify the forward chain"
                )
            if any(
                observation.provenance != chain_provenance
                for observation in candidate.evidence
            ):
                raise ValueError(
                    "query-context evidence must remain chain-only forward evidence"
                )

    def _validate_reverse_mapping_context(self) -> None:
        results = self.reverse_mapping_results
        if results is None:
            if self.reverse_alignment_provenance is not None:
                raise ValueError(
                    "reverse alignment provenance requires reverse mapping results"
                )
            if self.reverse_mapping_resource is not None:
                raise ValueError("reverse resource requires reverse mapping results")
            if (
                self.result_profile.scope.reverse_result
                is not ReverseCheckState.NOT_RUN
            ):
                raise ValueError("result profile reverse state must be NOT_RUN")
            return

        result_candidate_ids = tuple(result.forward_candidate_id for result in results)
        report_candidate_ids = tuple(
            candidate.candidate_id for candidate in self.candidates
        )
        if result_candidate_ids != report_candidate_ids:
            raise ValueError("reverse mapping results must match report candidates")

        check_states = {result.check_state for result in results}
        if len(check_states) > 1:
            raise ValueError("reverse mapping results must share one check state")
        check_state = next(iter(check_states), ReverseCheckState.NOT_RUN)
        if self.result_profile.scope.reverse_result is not check_state:
            raise ValueError("result profile reverse state must match reverse results")

        if check_state is ReverseCheckState.RUN:
            if self.reverse_alignment_provenance is None:
                raise ValueError("completed reverse mapping requires provenance")
            resource = self.reverse_mapping_resource
            if resource is None or resource.role is not UCSCBundleResourceRole.CHAIN:
                raise ValueError(
                    "completed reverse mapping requires its chain resource"
                )
            if not resource.consumed_by_engine or resource.file_provenance is None:
                raise ValueError("reverse chain resource must be marked consumed")
            if resource.file_provenance.derived_from != (
                self.reverse_alignment_provenance,
            ):
                raise ValueError(
                    "reverse chain provenance must derive from reverse alignment "
                    "provenance"
                )
            for result in results:
                for segment_result in result.segment_results:
                    for candidate in segment_result.candidates:
                        if candidate.mapping_provenance != resource.file_provenance:
                            raise ValueError(
                                "reverse candidate mapping provenance must identify "
                                "the consumed reverse chain"
                            )
        else:
            if self.reverse_alignment_provenance is not None:
                raise ValueError("unperformed reverse mapping cannot carry provenance")
            if self.reverse_mapping_resource is not None:
                raise ValueError(
                    "unperformed reverse mapping cannot consume a resource"
                )

    def _validate_filtered_all_chain_comparison(self) -> None:
        comparison = self.filtered_all_chain_comparison
        resource = self.filtered_chain_comparison_resource
        relationship = self.comparative_evidence_relationship
        profile = self.result_profile.comparative_relationship
        if self.result_profile.scope.comparative_relationship is not profile.state:
            raise ValueError(
                "result profile comparative scope must match comparative profile"
            )

        if comparison is None:
            if resource is not None:
                raise ValueError(
                    "filtered-chain comparison resource requires comparison results"
                )
            if relationship is not None:
                raise ValueError(
                    "comparative evidence relationship requires paired inventory"
                )
            if profile.state is not ComparativeRelationshipState.NOT_ASSESSED:
                raise ValueError(
                    "missing paired inventory must remain NOT_ASSESSED in the profile"
                )
            return
        if relationship is None:
            raise ValueError(
                "paired filtered/all-chain inventory requires comparative evidence "
                "relationship"
            )

        if self.evidence_tier is not EvidenceAvailabilityTier.COMPARATIVE:
            raise ValueError(
                "filtered/all-chain comparison requires a COMPARATIVE forward report"
            )
        if resource is None or resource.file_provenance is None:
            raise ValueError(
                "filtered/all-chain comparison requires consumed filtered-chain "
                "resource provenance"
            )
        if resource.role is not UCSCBundleResourceRole.CHAIN:
            raise ValueError("filtered/all-chain comparison resource must be a chain")
        if not resource.consumed_by_engine:
            raise ValueError("filtered/all-chain comparison chain must be consumed")
        if comparison.source_interval != self.source_interval:
            raise ValueError(
                "filtered/all-chain comparison source interval must match the report"
            )
        if comparison.all_chain_candidates != self.candidates:
            raise ValueError(
                "filtered/all-chain comparison all-chain inventory must match "
                "the report"
            )

        all_chain_resource = next(
            item for item in self.resources if item.role is UCSCBundleResourceRole.CHAIN
        )
        if all_chain_resource.file_provenance is None:
            raise ValueError(
                "filtered/all-chain comparison requires consumed all-chain provenance"
            )
        if comparison.all_chain_provenance != all_chain_resource.file_provenance:
            raise ValueError(
                "filtered/all-chain comparison must identify the report all-chain"
            )
        if comparison.filtered_chain_provenance != resource.file_provenance:
            raise ValueError(
                "filtered/all-chain comparison must identify the consumed filtered "
                "chain"
            )
        if resource.file_provenance.derived_from != (self.alignment_provenance,):
            raise ValueError(
                "filtered-chain comparison provenance must preserve the report "
                "dependency provenance"
            )
        expected_relationship = build_comparative_evidence_relationship(comparison)
        if relationship != expected_relationship:
            raise ValueError(
                "comparative evidence relationship must match the paired inventory "
                "and candidate evidence"
            )
        if profile.inventory_state is not comparison.relationship:
            raise ValueError(
                "result profile comparative inventory state must match the report"
            )
        if profile.state.value != relationship.relationship.value:
            raise ValueError(
                "result profile comparative relationship must match the report"
            )
        if profile.favored_candidate_id != relationship.favored_candidate_id:
            raise ValueError(
                "result profile favored candidate must match the report relationship"
            )
        if (
            profile.additional_all_chain_candidate_ids
            != comparison.additional_all_chain_candidate_ids
        ):
            raise ValueError(
                "result profile additional all-chain placements must match the report"
            )
        expected_support = tuple(
            (
                item.candidate_id,
                item.complete_source_coverage,
                item.retained_by_filtered_chain,
                item.depth1_top_net,
                item.full_reciprocal_best,
            )
            for item in relationship.placement_support
        )
        actual_support = tuple(
            (
                item.candidate_id,
                item.complete_source_coverage,
                item.retained_by_filtered_chain,
                item.depth1_top_net,
                item.full_reciprocal_best,
            )
            for item in profile.placement_support
        )
        if actual_support != expected_support:
            raise ValueError(
                "result profile comparative placement support must match the report"
            )


def assess_ucsc_cached_bundle(
    source_interval: GenomicInterval,
    bundle: CachedUCSCResourceBundle,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
    progress_callback: ResourceReadProgressCallback | None = None,
    chain_index: ChainIndex | None = None,
) -> UCSCAssessmentReport:
    """Run candidate/evidence generation and derive the factual result profile."""

    candidates = build_ucsc_candidates_from_cached_bundle(
        source_interval,
        bundle,
        target_assembly=target_assembly,
        alignment_provenance=alignment_provenance,
        progress_callback=progress_callback,
        chain_index=chain_index,
    )
    resources = _assessment_resources(
        bundle,
        alignment_provenance=alignment_provenance,
        candidates_exist=bool(candidates),
    )
    consumed_resource_roles = tuple(
        resource.role.value for resource in resources if resource.consumed_by_engine
    )
    result_profile = build_result_profile(
        source_interval,
        candidates,
        evidence_tier=bundle.evidence_tier,
        consumed_resource_roles=consumed_resource_roles,
    )
    return UCSCAssessmentReport(
        source_interval=source_interval,
        target_assembly=target_assembly,
        candidates=candidates,
        evidence_tier=bundle.evidence_tier,
        result_profile=result_profile,
        source_db=bundle.source_db,
        target_db=bundle.target_db,
        alignment_provenance=alignment_provenance,
        resources=resources,
    )


def attach_filtered_all_chain_comparison(
    report: UCSCAssessmentReport,
    *,
    filtered_chain: CachedUCSCChainResource,
    filtered_chain_index: ChainIndex | None,
) -> UCSCAssessmentReport:
    """Attach one indexed ordinary-filtered versus all-chain candidate inventory.

    The comparison remains index-only. A missing prepared filtered-chain index is
    an explicit precondition failure and never falls back to exhaustive traversal.
    The attached categorical relationship uses only filtered retention, depth-1 top-net
    support, full reciprocal-best membership, and complete-placement geometry; human
    and machine rendering remain a later layer.
    """

    if report.filtered_all_chain_comparison is not None:
        raise ValueError("filtered/all-chain comparison is already attached")
    if report.evidence_tier is not EvidenceAvailabilityTier.COMPARATIVE:
        raise ValueError(
            "filtered/all-chain comparison requires a COMPARATIVE forward report"
        )
    if (
        filtered_chain.source_db != report.source_db
        or filtered_chain.target_db != report.target_db
    ):
        raise ValueError("filtered liftOver chain must match the forward database pair")
    if filtered_chain.evidence_tier is not EvidenceAvailabilityTier.LIFTOVER_ONLY:
        raise ValueError(
            "filtered/all-chain comparison requires the ordinary filtered liftOver "
            "chain publication class"
        )
    if filtered_chain_index is None:
        raise ValueError(
            "filtered/all-chain comparison requires a prepared filtered-chain index"
        )

    filtered_chain_provenance = _cached_chain_resource_provenance(
        filtered_chain,
        alignment_provenance=report.alignment_provenance,
    )
    filtered_candidates = build_ucsc_chain_candidates_for_intervals_from_cached_chain(
        (report.source_interval,),
        filtered_chain,
        target_assembly=report.target_assembly,
        alignment_provenance=report.alignment_provenance,
        chain_index=filtered_chain_index,
    )[0]

    all_chain_resource = next(
        item for item in report.resources if item.role is UCSCBundleResourceRole.CHAIN
    )
    if all_chain_resource.file_provenance is None:
        raise ValueError(
            "filtered/all-chain comparison requires consumed all-chain provenance"
        )
    comparison = build_filtered_all_chain_comparison(
        report.source_interval,
        report.candidates,
        filtered_candidates,
        all_chain_provenance=all_chain_resource.file_provenance,
        filtered_chain_provenance=filtered_chain_provenance,
    )
    relationship = build_comparative_evidence_relationship(comparison)
    comparison_resource = UCSCAssessmentResource(
        role=UCSCBundleResourceRole.CHAIN,
        resource=filtered_chain.chain,
        consumed_by_engine=True,
        file_provenance=filtered_chain_provenance,
    )
    profile = build_result_profile(
        report.source_interval,
        report.candidates,
        evidence_tier=report.evidence_tier,
        consumed_resource_roles=report.result_profile.consumed_resource_roles,
        reverse_mapping_results=report.reverse_mapping_results,
        query_context_result=report.query_context_result,
        filtered_all_chain_comparison=comparison,
        comparative_evidence_relationship=relationship,
    )
    return replace(
        report,
        result_profile=profile,
        filtered_all_chain_comparison=comparison,
        comparative_evidence_relationship=relationship,
        filtered_chain_comparison_resource=comparison_resource,
    )


def attach_query_context_result(
    report: UCSCAssessmentReport,
    query_context_result: PointQueryContextResult,
) -> UCSCAssessmentReport:
    """Attach one already-computed point-neighborhood chain-geometry result."""

    if report.query_context_result is not None:
        raise ValueError("query context is already attached")
    if report.source_interval.length != 1:
        raise ValueError("query context can only be attached to a one-base report")

    profile = build_result_profile(
        report.source_interval,
        report.candidates,
        evidence_tier=report.evidence_tier,
        consumed_resource_roles=report.result_profile.consumed_resource_roles,
        reverse_mapping_results=report.reverse_mapping_results,
        query_context_result=query_context_result,
        filtered_all_chain_comparison=report.filtered_all_chain_comparison,
        comparative_evidence_relationship=report.comparative_evidence_relationship,
    )
    return replace(
        report,
        result_profile=profile,
        query_context_result=query_context_result,
    )


def attach_point_query_context(
    report: UCSCAssessmentReport,
    *,
    chain_context: CachedUCSCChainResource,
    chain_index: ChainIndex | None,
    requested_window_bases: int = DEFAULT_POINT_CONTEXT_BASES,
) -> UCSCAssessmentReport:
    """Run indexed chain-only local context for a one-base source query.

    This deliberately does not re-run net or reciprocal-best resources.  The context
    dimension describes local chain geometry under the same chain publication class as
    the forward assessment.  When no usable index or source-sequence bound is available,
    the result remains explicitly ``NOT_RUN`` rather than falling back to a full scan.
    """

    if report.query_context_result is not None:
        raise ValueError("query context is already attached")
    if report.source_interval.length != 1:
        raise ValueError("point query context requires a one-base report")
    if (
        chain_context.source_db != report.source_db
        or chain_context.target_db != report.target_db
    ):
        raise ValueError("point-context chain must match the forward database pair")
    if chain_context.evidence_tier is not report.evidence_tier:
        raise ValueError(
            "point-context chain publication class must match the forward assessment"
        )

    if chain_index is None:
        return attach_query_context_result(
            report,
            point_context_not_run(
                requested_window_bases=requested_window_bases,
                reason=QueryContextNotRunReason.INDEX_UNAVAILABLE,
            ),
        )

    source_sequence_query_bound = chain_index.source_sequence_query_bound(
        report.source_interval.sequence_name
    )
    if source_sequence_query_bound is None:
        return attach_query_context_result(
            report,
            point_context_not_run(
                requested_window_bases=requested_window_bases,
                reason=QueryContextNotRunReason.SOURCE_BOUNDS_UNAVAILABLE,
            ),
        )

    context_interval = build_centered_point_context_interval(
        report.source_interval,
        requested_window_bases=requested_window_bases,
        source_sequence_query_bound=source_sequence_query_bound,
    )
    context_candidates = build_ucsc_chain_candidates_for_intervals_from_cached_chain(
        (context_interval,),
        chain_context,
        target_assembly=report.target_assembly,
        alignment_provenance=report.alignment_provenance,
        chain_index=chain_index,
    )[0]
    return attach_query_context_result(
        report,
        PointQueryContextResult(
            check_state=QueryContextState.RUN,
            requested_window_bases=requested_window_bases,
            tested_source_interval=context_interval,
            candidates=context_candidates,
        ),
    )


def attach_reverse_mapping_results(
    report: UCSCAssessmentReport,
    reverse_mapping_results: tuple[CandidateReverseMappingResult, ...],
    *,
    reverse_chain: CachedUCSCChainResource | None = None,
    reverse_alignment_provenance: ProvenanceSource | None = None,
) -> UCSCAssessmentReport:
    """Attach already-computed candidate-level actual reverse facts.

    Completed reverse runs require the exact chain context and provenance that produced
    them. ``NOT_RUN`` and ``UNAVAILABLE`` results must not claim a consumed resource.
    """

    if report.reverse_mapping_results is not None:
        raise ValueError("reverse mapping context is already attached")
    if not report.candidates:
        if reverse_mapping_results:
            raise ValueError(
                "reverse mapping results cannot be attached without forward candidates"
            )
        return report

    states = {result.check_state for result in reverse_mapping_results}
    if len(states) > 1:
        raise ValueError("reverse mapping results must share one check state")
    check_state = next(iter(states), ReverseCheckState.NOT_RUN)

    reverse_resource: UCSCAssessmentResource | None = None
    if check_state is ReverseCheckState.RUN:
        if reverse_chain is None or reverse_alignment_provenance is None:
            raise ValueError(
                "completed reverse mapping requires reverse chain context and "
                "provenance"
            )
        if (
            reverse_chain.source_db != report.target_db
            or reverse_chain.target_db != report.source_db
        ):
            raise ValueError("reverse chain must invert the report UCSC database pair")
        if reverse_chain.evidence_tier is not report.evidence_tier:
            raise ValueError(
                "reverse chain publication class must match the forward assessment"
            )
        reverse_file_provenance = _cached_chain_resource_provenance(
            reverse_chain,
            alignment_provenance=reverse_alignment_provenance,
        )
        reverse_resource = UCSCAssessmentResource(
            role=UCSCBundleResourceRole.CHAIN,
            resource=reverse_chain.chain,
            consumed_by_engine=True,
            file_provenance=reverse_file_provenance,
        )
    elif reverse_chain is not None or reverse_alignment_provenance is not None:
        raise ValueError(
            "unperformed reverse mapping must not claim a consumed reverse resource"
        )

    profile = build_result_profile(
        report.source_interval,
        report.candidates,
        evidence_tier=report.evidence_tier,
        consumed_resource_roles=report.result_profile.consumed_resource_roles,
        reverse_mapping_results=reverse_mapping_results,
        query_context_result=report.query_context_result,
        filtered_all_chain_comparison=report.filtered_all_chain_comparison,
        comparative_evidence_relationship=report.comparative_evidence_relationship,
    )
    return replace(
        report,
        result_profile=profile,
        reverse_mapping_results=reverse_mapping_results,
        reverse_alignment_provenance=reverse_alignment_provenance,
        reverse_mapping_resource=reverse_resource,
    )


def attach_reverse_mapping_context(
    report: UCSCAssessmentReport,
    *,
    reverse_bundle: CachedUCSCResourceBundle | None,
    reverse_alignment_provenance: ProvenanceSource | None = None,
    progress_callback: ResourceReadProgressCallback | None = None,
    chain_index: ChainIndex | None = None,
) -> UCSCAssessmentReport:
    """Compatibility helper that executes reverse mapping from a cached bundle.

    Only the bundle's chain is consumed. ``None`` means no usable reverse resource was
    available and is represented as ``UNAVAILABLE``. Automatic CLI execution uses the
    narrower chain-only cache boundary directly.
    """

    if report.reverse_mapping_results is not None:
        raise ValueError("reverse mapping context is already attached")
    if not report.candidates:
        return report
    if reverse_bundle is None:
        return attach_reverse_mapping_results(
            report,
            tuple(
                reverse_mapping_unavailable(candidate)
                for candidate in report.candidates
            ),
        )
    if reverse_alignment_provenance is None:
        raise ValueError("reverse bundle requires reverse alignment provenance")
    if (
        reverse_bundle.source_db != report.target_db
        or reverse_bundle.target_db != report.source_db
    ):
        raise ValueError("reverse bundle must invert the report UCSC database pair")
    if reverse_bundle.evidence_tier is not report.evidence_tier:
        raise ValueError(
            "reverse bundle publication class must match the forward assessment"
        )

    results = build_reverse_mapping_results_from_cached_bundle(
        report.candidates,
        reverse_bundle,
        reverse_alignment_provenance=reverse_alignment_provenance,
        progress_callback=progress_callback,
        chain_index=chain_index,
    )
    return attach_reverse_mapping_results(
        report,
        results,
        reverse_chain=CachedUCSCChainResource(
            source_db=reverse_bundle.source_db,
            target_db=reverse_bundle.target_db,
            evidence_tier=reverse_bundle.evidence_tier,
            chain=reverse_bundle.chain,
        ),
        reverse_alignment_provenance=reverse_alignment_provenance,
    )


def _assessment_resources(
    bundle: CachedUCSCResourceBundle,
    *,
    alignment_provenance: ProvenanceSource,
    candidates_exist: bool,
) -> tuple[UCSCAssessmentResource, ...]:
    consumed_roles = _consumed_resource_roles(
        bundle.evidence_tier, candidates_exist=candidates_exist
    )

    return tuple(
        UCSCAssessmentResource(
            role=role,
            resource=resource,
            consumed_by_engine=role in consumed_roles,
            file_provenance=(
                _cached_bundle_resource_provenance(
                    bundle,
                    role,
                    resource,
                    alignment_provenance=alignment_provenance,
                )
                if role in consumed_roles
                else None
            ),
        )
        for role, resource in _bundle_resources(bundle)
    )


def _bundle_resources(
    bundle: CachedUCSCResourceBundle,
) -> tuple[tuple[UCSCBundleResourceRole, CachedResource], ...]:
    resources: list[tuple[UCSCBundleResourceRole, CachedResource]] = [
        (UCSCBundleResourceRole.CHAIN, bundle.chain)
    ]
    if bundle.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        assert bundle.net is not None
        assert bundle.syntenic_net is not None
        assert bundle.reciprocal_best_chain is not None
        assert bundle.reciprocal_best_net is not None
        resources.extend(
            (
                (UCSCBundleResourceRole.NET, bundle.net),
                (UCSCBundleResourceRole.SYNTENIC_NET, bundle.syntenic_net),
                (
                    UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
                    bundle.reciprocal_best_chain,
                ),
                (
                    UCSCBundleResourceRole.RECIPROCAL_BEST_NET,
                    bundle.reciprocal_best_net,
                ),
            )
        )
    return tuple(resources)


def _consumed_resource_roles(
    evidence_tier: EvidenceAvailabilityTier,
    *,
    candidates_exist: bool,
) -> frozenset[UCSCBundleResourceRole]:
    roles = {UCSCBundleResourceRole.CHAIN}
    if candidates_exist and evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        roles.update(
            {
                UCSCBundleResourceRole.NET,
                UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
            }
        )
    return frozenset(roles)


def _resource_roles_for_tier(
    evidence_tier: EvidenceAvailabilityTier,
) -> tuple[UCSCBundleResourceRole, ...]:
    if evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return (UCSCBundleResourceRole.CHAIN,)
    if evidence_tier is not EvidenceAvailabilityTier.COMPARATIVE:
        raise ValueError(f"unsupported evidence tier: {evidence_tier!r}")
    return (
        UCSCBundleResourceRole.CHAIN,
        UCSCBundleResourceRole.NET,
        UCSCBundleResourceRole.SYNTENIC_NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
        UCSCBundleResourceRole.RECIPROCAL_BEST_NET,
    )
