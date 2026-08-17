"""Assessment/report orchestration for already-acquired UCSC resources.

This module is the first end-to-end assessment boundary in v1. It connects one exact
cached UCSC resource bundle to the existing file-backed candidate engine and then to
the deterministic assessor. It deliberately does not discover or download resources;
network policy, terms acknowledgement, and transfer confirmation remain separate
boundaries that the future CLI can compose around this function.

The report retains every artifact in the acquired bundle as retrieval context, but it
marks resource consumption explicitly. This distinction matters scientifically:
``COMPARATIVE`` bundles contain syntenic-net and reciprocal-best-net files that the
current engine does not parse, and when chain projection produces no candidates the
engine returns before consuming the ordinary net or reciprocal-best chain. Cached
presence must therefore never be presented as evidence that a file contributed to a
verdict.
"""

from dataclasses import dataclass

from .assessor import assess_candidates
from .models import (
    AssemblyIdentifier,
    Assessment,
    EvidenceAvailabilityTier,
    GenomicInterval,
    ProvenanceIdentifierKind,
    ProvenanceSource,
)
from .resource_cache import (
    CachedResource,
    CachedUCSCResourceBundle,
    UCSCBundleResourceRole,
)
from .resource_files import (
    _cached_bundle_resource_provenance,
    build_ucsc_candidates_from_cached_bundle,
)


@dataclass(frozen=True)
class UCSCAssessmentResource:
    """One cached bundle artifact and its actual role in an assessment run.

    ``resource`` preserves retrieval metadata such as source URL, retrieval time,
    provider checksum information, terms references, and exact cache content digest.
    ``consumed_by_engine`` says whether the current run actually parsed that artifact.
    ``file_provenance`` is present only for consumed artifacts, preventing retrieval
    context from being mistaken for evidence provenance. Chain and reciprocal-best
    observations use their file provenance directly; net observations add a per-fill
    provenance node whose parent is the net file provenance recorded here.
    """

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
    """End-to-end v1 assessment plus auditable UCSC retrieval context."""

    assessment: Assessment
    source_db: str
    target_db: str
    alignment_provenance: ProvenanceSource
    resources: tuple[UCSCAssessmentResource, ...]

    def __post_init__(self) -> None:
        expected_roles = _resource_roles_for_tier(self.assessment.evidence_tier)
        actual_roles = tuple(resource.role for resource in self.resources)
        if actual_roles != expected_roles:
            raise ValueError(
                "UCSC assessment report resources must preserve the complete ordered "
                "bundle roles for the assessment evidence tier"
            )

        expected_consumed_roles = _consumed_resource_roles(
            self.assessment.evidence_tier,
            candidates_exist=bool(self.assessment.candidates),
        )
        actual_consumed_roles = {
            resource.role for resource in self.resources if resource.consumed_by_engine
        }
        if actual_consumed_roles != expected_consumed_roles:
            raise ValueError(
                "UCSC assessment report resource-consumption metadata does not match "
                "the v1 engine path"
            )

        for resource in self.resources:
            if resource.file_provenance is not None and (
                resource.file_provenance.derived_from != (self.alignment_provenance,)
            ):
                raise ValueError(
                    "consumed UCSC file provenance must derive from the report alignment "
                    "provenance"
                )

    @property
    def evidence_tier(self) -> EvidenceAvailabilityTier:
        """Expose evidence availability separately from the assessment verdict."""

        return self.assessment.evidence_tier


def assess_ucsc_cached_bundle(
    source_interval: GenomicInterval,
    bundle: CachedUCSCResourceBundle,
    *,
    target_assembly: AssemblyIdentifier,
    alignment_provenance: ProvenanceSource,
) -> UCSCAssessmentReport:
    """Run the v1 candidate engine and assessor over one cached UCSC bundle.

    ``target_assembly`` is intentionally forwarded to the existing cached-bundle
    bridge, which owns assembly-pair validation. The report's evidence tier comes from
    the verified bundle shape and is passed independently to the assessor; it is never
    inferred from the resulting verdict.
    """

    candidates = build_ucsc_candidates_from_cached_bundle(
        source_interval,
        bundle,
        target_assembly=target_assembly,
        alignment_provenance=alignment_provenance,
    )
    assessment = assess_candidates(
        source_interval,
        candidates,
        evidence_tier=bundle.evidence_tier,
    )
    resources = _assessment_resources(
        bundle,
        alignment_provenance=alignment_provenance,
        candidates_exist=bool(candidates),
    )
    return UCSCAssessmentReport(
        assessment=assessment,
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
