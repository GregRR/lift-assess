"""Filtered liftOver versus all-chain candidate-inventory comparison.

This module records one factual cross-resource relationship: whether the ordinary
filtered liftOver chain yields the same local mapping inventory as the broader UCSC
all-chain resource, or whether the all-chain adds placements that the filtered chain
suppresses.  Correspondence is based on canonical normalized mapping geometry, never
candidate encounter order, chain score, or candidate identifiers.

The result is deliberately below the later comparative-interpretation layer.  It does
not synthesize net/reciprocal-best evidence, rank candidates, assign confidence, or
name a biologically correct locus.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from ._candidate_geometry import (
    canonical_mapping_geometry,
    validate_distinct_candidate_geometries,
)
from .models import GenomicInterval, NormalizedCandidate, ProvenanceSource


class FilteredAllChainCorrespondenceError(ValueError):
    """Raised when filtered and all-chain placements cannot be paired safely."""


class FilteredAllChainInventoryState(str, Enum):
    """Factual relationship between filtered and all-chain candidate inventories."""

    FILTERED_AND_ALL_CHAIN_AGREE = "FILTERED_AND_ALL_CHAIN_AGREE"
    ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS = "ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS"


@dataclass(frozen=True)
class FilteredAllChainCandidateMatch:
    """One geometry-safe correspondence between the two chain publications."""

    filtered_candidate_id: str
    all_chain_candidate_id: str

    def __post_init__(self) -> None:
        if not self.filtered_candidate_id:
            raise ValueError("filtered candidate ID must not be empty")
        if not self.all_chain_candidate_id:
            raise ValueError("all-chain candidate ID must not be empty")


@dataclass(frozen=True)
class FilteredAllChainComparisonResult:
    """Structured paired inventory for one exact source interval.

    Candidate ID order is retained only for reproducibility.  ``candidate_matches``
    follows filtered-chain encounter order and ``additional_all_chain_candidate_ids``
    follows all-chain encounter order; neither order is scientific rank.
    """

    source_interval: GenomicInterval
    relationship: FilteredAllChainInventoryState
    all_chain_candidates: tuple[NormalizedCandidate, ...]
    filtered_candidates: tuple[NormalizedCandidate, ...]
    candidate_matches: tuple[FilteredAllChainCandidateMatch, ...]
    additional_all_chain_candidate_ids: tuple[str, ...]
    all_chain_provenance: ProvenanceSource
    filtered_chain_provenance: ProvenanceSource

    def __post_init__(self) -> None:
        if self.source_interval.length <= 0:
            raise ValueError("filtered/all-chain comparison requires a non-empty query")

        _validate_candidate_tuple(self.all_chain_candidates, label="all-chain")
        _validate_candidate_tuple(self.filtered_candidates, label="filtered-chain")
        _validate_candidate_mapping_provenance(
            self.all_chain_candidates,
            self.all_chain_provenance,
            label="all-chain",
        )
        _validate_candidate_mapping_provenance(
            self.filtered_candidates,
            self.filtered_chain_provenance,
            label="filtered-chain",
        )
        _validate_shared_alignment_lineage(
            self.all_chain_provenance,
            self.filtered_chain_provenance,
        )

        all_chain_by_id = {
            candidate.candidate_id: candidate for candidate in self.all_chain_candidates
        }
        filtered_by_id = {
            candidate.candidate_id: candidate for candidate in self.filtered_candidates
        }
        filtered_ids = tuple(filtered_by_id)
        matched_filtered_ids = tuple(
            match.filtered_candidate_id for match in self.candidate_matches
        )
        if matched_filtered_ids != filtered_ids:
            raise ValueError(
                "filtered/all-chain matches must cover filtered candidates in "
                "filtered encounter order"
            )

        matched_all_chain_ids = tuple(
            match.all_chain_candidate_id for match in self.candidate_matches
        )
        if len(set(matched_all_chain_ids)) != len(matched_all_chain_ids):
            raise ValueError(
                "one all-chain candidate cannot correspond to multiple filtered "
                "candidates"
            )

        for match in self.candidate_matches:
            all_chain_candidate = all_chain_by_id.get(match.all_chain_candidate_id)
            if all_chain_candidate is None:
                raise ValueError(
                    "filtered/all-chain match references unknown all-chain ID"
                )
            filtered_candidate = filtered_by_id[match.filtered_candidate_id]
            if canonical_mapping_geometry(filtered_candidate) != (
                canonical_mapping_geometry(all_chain_candidate)
            ):
                raise ValueError(
                    "filtered/all-chain candidate match must preserve canonical "
                    "mapping geometry"
                )

        matched_all_chain_id_set = set(matched_all_chain_ids)
        expected_additional = tuple(
            candidate.candidate_id
            for candidate in self.all_chain_candidates
            if candidate.candidate_id not in matched_all_chain_id_set
        )
        if self.additional_all_chain_candidate_ids != expected_additional:
            raise ValueError(
                "additional all-chain candidate IDs must be the unmatched inventory"
            )

        agree = FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
        additional = (
            FilteredAllChainInventoryState.ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS
        )
        expected_relationship = agree if not expected_additional else additional
        if self.relationship is not expected_relationship:
            raise ValueError(
                "filtered/all-chain relationship does not match the paired inventory"
            )

        for candidate in self.filtered_candidates:
            if any(
                observation.provenance != self.filtered_chain_provenance
                for observation in candidate.evidence
            ):
                raise ValueError(
                    "filtered candidate evidence must remain filtered-chain-only"
                )

    @property
    def all_chain_candidate_ids(self) -> tuple[str, ...]:
        """All-chain candidate IDs in reproducible encounter order, not rank."""

        return tuple(candidate.candidate_id for candidate in self.all_chain_candidates)


def build_filtered_all_chain_comparison(
    source_interval: GenomicInterval,
    all_chain_candidates: tuple[NormalizedCandidate, ...],
    filtered_candidates: tuple[NormalizedCandidate, ...],
    *,
    all_chain_provenance: ProvenanceSource,
    filtered_chain_provenance: ProvenanceSource,
) -> FilteredAllChainComparisonResult:
    """Pair filtered candidates to all-chain candidates by canonical geometry.

    Every filtered placement must have exactly matching canonical geometry in the
    all-chain inventory before liftAssess will synthesize a comparative relationship.
    UCSC's ordinary liftOver-chain construction can clip an original chain to net-fill
    boundaries, so a mismatch is a correspondence limitation, not evidence that the
    primary assessment is invalid.
    """

    if source_interval.length <= 0:
        raise ValueError("filtered/all-chain comparison requires a non-empty query")

    _validate_candidate_tuple(all_chain_candidates, label="all-chain")
    _validate_candidate_tuple(filtered_candidates, label="filtered-chain")
    _validate_candidate_mapping_provenance(
        all_chain_candidates,
        all_chain_provenance,
        label="all-chain",
    )
    _validate_candidate_mapping_provenance(
        filtered_candidates,
        filtered_chain_provenance,
        label="filtered-chain",
    )

    _validate_shared_alignment_lineage(
        all_chain_provenance,
        filtered_chain_provenance,
    )

    all_chain_by_geometry = {
        canonical_mapping_geometry(candidate): candidate
        for candidate in all_chain_candidates
    }
    matches: list[FilteredAllChainCandidateMatch] = []
    matched_all_chain_ids: set[str] = set()
    for filtered_candidate in filtered_candidates:
        all_chain_candidate = all_chain_by_geometry.get(
            canonical_mapping_geometry(filtered_candidate)
        )
        if all_chain_candidate is None:
            raise FilteredAllChainCorrespondenceError(
                "filtered-chain placement cannot be paired to identical all-chain "
                "geometry: "
                f"{filtered_candidate.candidate_id!r}"
            )
        matches.append(
            FilteredAllChainCandidateMatch(
                filtered_candidate_id=filtered_candidate.candidate_id,
                all_chain_candidate_id=all_chain_candidate.candidate_id,
            )
        )
        matched_all_chain_ids.add(all_chain_candidate.candidate_id)

    additional_ids = tuple(
        candidate.candidate_id
        for candidate in all_chain_candidates
        if candidate.candidate_id not in matched_all_chain_ids
    )
    agree = FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    additional = FilteredAllChainInventoryState.ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS
    relationship = agree if not additional_ids else additional
    return FilteredAllChainComparisonResult(
        source_interval=source_interval,
        relationship=relationship,
        all_chain_candidates=all_chain_candidates,
        filtered_candidates=filtered_candidates,
        candidate_matches=tuple(matches),
        additional_all_chain_candidate_ids=additional_ids,
        all_chain_provenance=all_chain_provenance,
        filtered_chain_provenance=filtered_chain_provenance,
    )


def _validate_candidate_tuple(
    candidates: tuple[NormalizedCandidate, ...],
    *,
    label: str,
) -> None:
    _validate_candidate_ids(
        tuple(candidate.candidate_id for candidate in candidates),
        label=label,
    )
    validate_distinct_candidate_geometries(candidates)


def _validate_candidate_ids(candidate_ids: tuple[str, ...], *, label: str) -> None:
    if any(not candidate_id for candidate_id in candidate_ids):
        raise ValueError(f"{label} candidate IDs must not be empty")
    if len(set(candidate_ids)) != len(candidate_ids):
        raise ValueError(f"{label} candidate IDs must be unique")


def _validate_shared_alignment_lineage(
    all_chain_provenance: ProvenanceSource,
    filtered_chain_provenance: ProvenanceSource,
) -> None:
    if (
        not all_chain_provenance.derived_from
        or all_chain_provenance.derived_from != filtered_chain_provenance.derived_from
    ):
        raise ValueError(
            "filtered and all-chain provenance must preserve the same upstream "
            "alignment lineage"
        )


def _validate_candidate_mapping_provenance(
    candidates: tuple[NormalizedCandidate, ...],
    provenance: ProvenanceSource,
    *,
    label: str,
) -> None:
    for candidate in candidates:
        if candidate.mapping_provenance != provenance:
            raise ValueError(
                f"{label} candidate mapping provenance must identify its chain"
            )
