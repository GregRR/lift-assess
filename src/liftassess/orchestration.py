"""Scientific report orchestration for already-acquired UCSC resources.

This module connects one exact cached UCSC resource bundle to the file-backed
candidate/evidence engine and the derived factual result profile.  It deliberately
does not discover or download resources; network policy, terms acknowledgement, and
transfer confirmation remain separate CLI boundaries.

The report retains every artifact in the acquired bundle as retrieval context while
marking actual engine consumption explicitly. Cached presence is never presented as
scientific evidence unless the engine consumed the resource.
"""

from dataclasses import dataclass

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
    CachedUCSCResourceBundle,
    UCSCBundleResourceRole,
)
from .resource_files import (
    ResourceReadProgressCallback,
    _cached_bundle_resource_provenance,
    build_ucsc_candidates_from_cached_bundle,
)
from .result_profile import ResultProfile, build_result_profile


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

        for resource in self.resources:
            if resource.file_provenance is not None and (
                resource.file_provenance.derived_from != (self.alignment_provenance,)
            ):
                raise ValueError(
                    "consumed UCSC file provenance must derive from the report "
                    "alignment provenance"
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
