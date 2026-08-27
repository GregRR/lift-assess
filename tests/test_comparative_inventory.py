from __future__ import annotations

from dataclasses import replace

import pytest

from liftassess import (
    AssemblyIdentifier,
    FilteredAllChainCorrespondenceError,
    FilteredAllChainInventoryState,
    GenomicInterval,
    MappingOrientation,
    MappingSegment,
    NormalizedCandidate,
    ProvenanceSource,
    build_filtered_all_chain_comparison,
)

SOURCE_ASSEMBLY = AssemblyIdentifier("sourceAsm", "test")
TARGET_ASSEMBLY = AssemblyIdentifier("targetAsm", "test")
SOURCE_INTERVAL = GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 110)
ALIGNMENT = ProvenanceSource("alignment", "shared alignment lineage")
ALL_CHAIN = ProvenanceSource(
    "all-chain",
    "all-chain bytes",
    derived_from=(ALIGNMENT,),
)
FILTERED_CHAIN = ProvenanceSource(
    "filtered-chain",
    "filtered liftOver chain bytes",
    derived_from=(ALIGNMENT,),
)


def _candidate(
    candidate_id: str,
    *,
    provenance: ProvenanceSource,
    orientation: MappingOrientation = MappingOrientation.SAME,
    target_sequence: str = "chrA",
    target_start: int = 500,
    split: bool = False,
) -> NormalizedCandidate:
    segments: tuple[MappingSegment, ...]
    if split and orientation is MappingOrientation.SAME:
        segments = (
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 105),
                GenomicInterval(
                    TARGET_ASSEMBLY,
                    target_sequence,
                    target_start,
                    target_start + 5,
                ),
            ),
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 105, 110),
                GenomicInterval(
                    TARGET_ASSEMBLY,
                    target_sequence,
                    target_start + 5,
                    target_start + 10,
                ),
            ),
        )
    elif split:
        segments = (
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 100, 105),
                GenomicInterval(
                    TARGET_ASSEMBLY,
                    target_sequence,
                    target_start + 5,
                    target_start + 10,
                ),
            ),
            MappingSegment(
                GenomicInterval(SOURCE_ASSEMBLY, "chr1", 105, 110),
                GenomicInterval(
                    TARGET_ASSEMBLY,
                    target_sequence,
                    target_start,
                    target_start + 5,
                ),
            ),
        )
    else:
        segments = (
            MappingSegment(
                SOURCE_INTERVAL,
                GenomicInterval(
                    TARGET_ASSEMBLY,
                    target_sequence,
                    target_start,
                    target_start + 10,
                ),
            ),
        )

    return NormalizedCandidate(
        candidate_id=candidate_id,
        target_interval=GenomicInterval(
            TARGET_ASSEMBLY,
            target_sequence,
            target_start,
            target_start + 10,
        ),
        orientation=orientation,
        mapping_provenance=provenance,
        segments=segments,
    )


@pytest.mark.parametrize(
    "orientation",
    (MappingOrientation.SAME, MappingOrientation.REVERSE),
)
def test_filtered_candidate_corresponds_by_canonical_geometry_across_orientations(
    orientation: MappingOrientation,
) -> None:
    all_chain = _candidate(
        "all",
        provenance=ALL_CHAIN,
        orientation=orientation,
        split=True,
    )
    filtered = _candidate(
        "filtered",
        provenance=FILTERED_CHAIN,
        orientation=orientation,
    )

    comparison = build_filtered_all_chain_comparison(
        SOURCE_INTERVAL,
        (all_chain,),
        (filtered,),
        all_chain_provenance=ALL_CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )

    assert (
        comparison.relationship
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    )
    assert comparison.candidate_matches[0].filtered_candidate_id == "filtered"
    assert comparison.candidate_matches[0].all_chain_candidate_id == "all"
    assert comparison.additional_all_chain_candidate_ids == ()


def test_all_chain_additional_placements_remain_explicit_and_ordered() -> None:
    retained = _candidate("retained", provenance=ALL_CHAIN)
    extra_one = _candidate(
        "extra-one",
        provenance=ALL_CHAIN,
        target_sequence="chrB",
        target_start=700,
    )
    extra_two = _candidate(
        "extra-two",
        provenance=ALL_CHAIN,
        target_sequence="chrC",
        target_start=900,
    )
    filtered = _candidate("filtered", provenance=FILTERED_CHAIN)

    comparison = build_filtered_all_chain_comparison(
        SOURCE_INTERVAL,
        (extra_one, retained, extra_two),
        (filtered,),
        all_chain_provenance=ALL_CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )

    assert (
        comparison.relationship
        is FilteredAllChainInventoryState.ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS
    )
    assert comparison.additional_all_chain_candidate_ids == ("extra-one", "extra-two")


def test_zero_filtered_and_zero_all_chain_is_exact_inventory_agreement() -> None:
    comparison = build_filtered_all_chain_comparison(
        SOURCE_INTERVAL,
        (),
        (),
        all_chain_provenance=ALL_CHAIN,
        filtered_chain_provenance=FILTERED_CHAIN,
    )

    assert (
        comparison.relationship
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    )
    assert comparison.candidate_matches == ()
    assert comparison.additional_all_chain_candidate_ids == ()


def test_filtered_placement_missing_from_all_chain_is_invariant_failure() -> None:
    all_chain = _candidate("all", provenance=ALL_CHAIN)
    filtered = _candidate(
        "filtered",
        provenance=FILTERED_CHAIN,
        target_sequence="chrB",
        target_start=700,
    )

    with pytest.raises(
        FilteredAllChainCorrespondenceError,
        match="cannot be paired to identical all-chain geometry",
    ):
        build_filtered_all_chain_comparison(
            SOURCE_INTERVAL,
            (all_chain,),
            (filtered,),
            all_chain_provenance=ALL_CHAIN,
            filtered_chain_provenance=FILTERED_CHAIN,
        )


def test_orientation_mismatch_is_inconsistent_geometry() -> None:
    all_chain = _candidate("all", provenance=ALL_CHAIN)
    filtered = _candidate("filtered", provenance=FILTERED_CHAIN)
    filtered = replace(filtered, orientation=MappingOrientation.REVERSE)

    with pytest.raises(
        FilteredAllChainCorrespondenceError,
        match="cannot be paired to identical all-chain geometry",
    ):
        build_filtered_all_chain_comparison(
            SOURCE_INTERVAL,
            (all_chain,),
            (filtered,),
            all_chain_provenance=ALL_CHAIN,
            filtered_chain_provenance=FILTERED_CHAIN,
        )


def test_comparison_rejects_distinct_upstream_lineages() -> None:
    unrelated = ProvenanceSource("other", "different alignment lineage")
    filtered_provenance = ProvenanceSource(
        "filtered-other",
        "filtered chain from another lineage",
        derived_from=(unrelated,),
    )

    with pytest.raises(ValueError, match="same upstream dependency group"):
        build_filtered_all_chain_comparison(
            SOURCE_INTERVAL,
            (),
            (),
            all_chain_provenance=ALL_CHAIN,
            filtered_chain_provenance=filtered_provenance,
        )
