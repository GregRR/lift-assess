"""Categorical comparative relationships across UCSC alignment-derived evidence.

This module interprets the already-paired filtered/all-chain inventory using only
explicit categorical observations.  It does not use chain score, ``ali``, ``qDup``,
encounter order, hidden weights, or a numeric confidence score.

The first accepted positive rule is intentionally narrow and mirrors the B14-style
relationship from the validation corpus: multiple complete all-chain placements,
exactly one placement retained by the ordinary filtered liftOver chain, that same
placement represented by a depth-1 ``top`` net fill and full reciprocal-best
membership, and no competing complete placement with the same joint top-net + full
reciprocal-best support.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .comparative_inventory import FilteredAllChainComparisonResult
from .models import (
    EvidenceKind,
    MappingCoverageStatus,
    MappingCoverageSummary,
    NetHierarchySummary,
    NormalizedCandidate,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
)
from .net import NetClassification


class ComparativeEvidenceRelationship(str, Enum):
    """Categorical relationship among complete all-chain placements."""

    NO_COMPETING_FULL_PLACEMENTS = "NO_COMPETING_FULL_PLACEMENTS"
    FAVORS_ONE_PLACEMENT = "FAVORS_ONE_PLACEMENT"
    DOES_NOT_SEPARATE_PLACEMENTS = "DOES_NOT_SEPARATE_PLACEMENTS"
    MIXED_CONFLICTING = "MIXED_CONFLICTING"


@dataclass(frozen=True)
class ComparativePlacementSupport:
    """Categorical comparative observations for one all-chain placement."""

    candidate_id: str
    complete_source_coverage: bool
    retained_by_filtered_chain: bool
    depth1_top_net: bool
    full_reciprocal_best: bool

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("comparative placement candidate ID must not be empty")

    @property
    def top_net_and_full_reciprocal_best(self) -> bool:
        """Return whether both accepted categorical support observations are present."""

        return self.depth1_top_net and self.full_reciprocal_best


@dataclass(frozen=True)
class ComparativeEvidenceRelationshipResult:
    """Deterministic categorical synthesis over one paired chain inventory.

    Placement-support order follows the all-chain inventory only for reproducibility;
    it is not scientific rank.
    """

    relationship: ComparativeEvidenceRelationship
    placement_support: tuple[ComparativePlacementSupport, ...]
    favored_candidate_id: str | None = None

    def __post_init__(self) -> None:
        candidate_ids = tuple(item.candidate_id for item in self.placement_support)
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError("comparative placement support IDs must be unique")

        expected_relationship, expected_favored_id = _relationship_for_support(
            self.placement_support
        )
        if self.relationship is not expected_relationship:
            raise ValueError(
                "comparative evidence relationship does not match categorical "
                "placement support"
            )
        if self.favored_candidate_id != expected_favored_id:
            raise ValueError(
                "comparative favored candidate does not match categorical placement "
                "support"
            )

    @property
    def full_candidate_ids(self) -> tuple[str, ...]:
        return tuple(
            item.candidate_id
            for item in self.placement_support
            if item.complete_source_coverage
        )

    @property
    def filtered_retained_candidate_ids(self) -> tuple[str, ...]:
        return tuple(
            item.candidate_id
            for item in self.placement_support
            if item.retained_by_filtered_chain
        )

    @property
    def filtered_retained_full_candidate_ids(self) -> tuple[str, ...]:
        return tuple(
            item.candidate_id
            for item in self.placement_support
            if item.complete_source_coverage and item.retained_by_filtered_chain
        )

    @property
    def depth1_top_net_full_candidate_ids(self) -> tuple[str, ...]:
        return tuple(
            item.candidate_id
            for item in self.placement_support
            if item.complete_source_coverage and item.depth1_top_net
        )

    @property
    def full_rbest_full_candidate_ids(self) -> tuple[str, ...]:
        return tuple(
            item.candidate_id
            for item in self.placement_support
            if item.complete_source_coverage and item.full_reciprocal_best
        )

    @property
    def joint_top_net_full_rbest_candidate_ids(self) -> tuple[str, ...]:
        return tuple(
            item.candidate_id
            for item in self.placement_support
            if item.complete_source_coverage and item.top_net_and_full_reciprocal_best
        )


def build_comparative_evidence_relationship(
    comparison: FilteredAllChainComparisonResult,
) -> ComparativeEvidenceRelationshipResult:
    """Classify comparative support without numeric weighting or candidate ranking."""

    retained_ids = {
        match.all_chain_candidate_id for match in comparison.candidate_matches
    }
    support = tuple(
        _placement_support(candidate, retained_ids=retained_ids)
        for candidate in comparison.all_chain_candidates
    )
    relationship, favored_id = _relationship_for_support(support)
    return ComparativeEvidenceRelationshipResult(
        relationship=relationship,
        placement_support=support,
        favored_candidate_id=favored_id,
    )


def _relationship_for_support(
    support: tuple[ComparativePlacementSupport, ...],
) -> tuple[ComparativeEvidenceRelationship, str | None]:
    full_support = tuple(item for item in support if item.complete_source_coverage)
    if len(full_support) < 2:
        return ComparativeEvidenceRelationship.NO_COMPETING_FULL_PLACEMENTS, None

    retained = tuple(item for item in support if item.retained_by_filtered_chain)
    filtered_full = tuple(
        item for item in full_support if item.retained_by_filtered_chain
    )
    joint = tuple(
        item for item in full_support if item.top_net_and_full_reciprocal_best
    )
    if (
        len(retained) == 1
        and len(filtered_full) == 1
        and len(joint) == 1
        and filtered_full[0].candidate_id == joint[0].candidate_id
    ):
        favored_id = filtered_full[0].candidate_id
        return ComparativeEvidenceRelationship.FAVORS_ONE_PLACEMENT, favored_id

    if _has_explicit_categorical_conflict(full_support):
        return ComparativeEvidenceRelationship.MIXED_CONFLICTING, None
    return ComparativeEvidenceRelationship.DOES_NOT_SEPARATE_PLACEMENTS, None


def _placement_support(
    candidate: NormalizedCandidate,
    *,
    retained_ids: set[str],
) -> ComparativePlacementSupport:
    coverage = _mapping_coverage(candidate)
    return ComparativePlacementSupport(
        candidate_id=candidate.candidate_id,
        complete_source_coverage=coverage.status is MappingCoverageStatus.FULL,
        retained_by_filtered_chain=candidate.candidate_id in retained_ids,
        depth1_top_net=_has_depth1_top_net(candidate),
        full_reciprocal_best=_has_full_reciprocal_best(candidate),
    )


def _mapping_coverage(candidate: NormalizedCandidate) -> MappingCoverageSummary:
    observations = tuple(
        observation
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.MAPPING_COVERAGE
    )
    if len(observations) != 1 or not isinstance(
        observations[0].value, MappingCoverageSummary
    ):
        raise ValueError(
            f"candidate {candidate.candidate_id!r} must carry exactly one typed "
            "mapping-coverage observation"
        )
    return observations[0].value


def _has_depth1_top_net(candidate: NormalizedCandidate) -> bool:
    """Require ``top`` classification and depth 1 from the same net fill."""

    top_fill_ids = {
        observation.provenance.source_id
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.NET_CLASSIFICATION
        and observation.value == NetClassification.TOP.value
    }
    depth1_fill_ids = {
        observation.provenance.source_id
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.NET_HIERARCHY
        and isinstance(observation.value, NetHierarchySummary)
        and observation.value.depth == 1
    }
    return bool(top_fill_ids & depth1_fill_ids)


def _has_full_reciprocal_best(candidate: NormalizedCandidate) -> bool:
    observations = tuple(
        observation
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
    )
    if len(observations) != 1 or not isinstance(
        observations[0].value, ReciprocalBestMembershipSummary
    ):
        raise ValueError(
            f"candidate {candidate.candidate_id!r} must carry exactly one typed "
            "reciprocal-best observation"
        )
    return observations[0].value.status is ReciprocalBestMembershipStatus.FULL


def _has_explicit_categorical_conflict(
    full_support: tuple[ComparativePlacementSupport, ...],
) -> bool:
    """Return whether uniquely identifying categorical observations disagree.

    Each observation family is treated as a categorical relationship, not a vote.
    A conflict exists only when at least two families each identify exactly one
    complete placement and those placements differ.  Non-unique or absent support is
    non-separating rather than silently weighted.
    """

    families = (
        tuple(
            item.candidate_id
            for item in full_support
            if item.retained_by_filtered_chain
        ),
        tuple(item.candidate_id for item in full_support if item.depth1_top_net),
        tuple(item.candidate_id for item in full_support if item.full_reciprocal_best),
    )
    unique_ids = tuple(ids[0] for ids in families if len(ids) == 1)
    return len(set(unique_ids)) > 1
