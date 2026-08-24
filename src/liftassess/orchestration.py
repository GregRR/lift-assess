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
from .models import (
    AssemblyIdentifier,
    EvidenceAvailabilityTier,
    GenomicInterval,
    NormalizedCandidate,
    ProvenanceIdentifierKind,
    ProvenanceSource,
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
)
from .result_profile import ResultProfile, build_result_profile
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
    reverse_mapping_results: tuple[CandidateReverseMappingResult, ...] | None = None
    reverse_alignment_provenance: ProvenanceSource | None = None
    reverse_mapping_resource: UCSCAssessmentResource | None = None

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

        self._validate_reverse_mapping_context()

        for resource in self.resources:
            if resource.file_provenance is not None and (
                resource.file_provenance.derived_from != (self.alignment_provenance,)
            ):
                raise ValueError(
                    "consumed UCSC file provenance must derive from the report "
                    "alignment provenance"
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
