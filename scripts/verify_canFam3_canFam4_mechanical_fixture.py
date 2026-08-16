from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, cast

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCResourceBundle,
    ChainGapSummary,
    EvidenceAvailabilityTier,
    EvidenceKind,
    GenomicInterval,
    MappingCoverageStatus,
    MappingCoverageSummary,
    MappingOrientation,
    NetHierarchySummary,
    ProvenanceIdentifierKind,
    ProvenanceSource,
    ProviderChecksum,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
    ResourceChecksumAlgorithm,
    build_ucsc_candidates_from_cached_bundle,
    iter_chain_file,
    ucsc_resource_terms,
)
from liftassess.models import EvidenceObservation, NormalizedCandidate

SOURCE_DB = "canFam3"
TARGET_DB = "canFam4"
SOURCE_SEQUENCE = "chrUn_JH373233"
SOURCE_START = 1_845_735
SOURCE_END = 1_845_835

FORWARD_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
RBEST_BASE = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/"
    "canFam4/vsCanFam3/reciprocalBest/"
)

RESOURCE_EXPECTATIONS = {
    "chain": (
        f"{FORWARD_BASE}canFam3.canFam4.all.chain.gz",
        "f10a6b48b5461bb8378ffaff311fb7355b1910511131ce0f5df5402c4db67519",
        2_652_632_416,
    ),
    "net": (
        f"{FORWARD_BASE}canFam3.canFam4.net.gz",
        "c889134e95ff82741c0092b1673b3e5fe0125aa82f611ee97d91e468b56d51ac",
        10_518_762,
    ),
    "syntenic_net": (
        f"{FORWARD_BASE}canFam3.canFam4.syn.net.gz",
        "39ba8ca12f935755ced5eaa555b9b476460c615cc7d6c0f122ac43194a9fabce",
        9_511_838,
    ),
    "reciprocal_best_chain": (
        f"{RBEST_BASE}canFam3.canFam4.rbest.chain.gz",
        "34f4061fd29e7720c7eb2adc1ea8299e86f21f08e18f97f9ba468cf8b466690c",
        5_403_921,
    ),
    "reciprocal_best_net": (
        f"{RBEST_BASE}canFam3.canFam4.rbest.net.gz",
        "cd4891e80eaa8b8625620162306c50fc291eab6912227d78c566d9dae7fe716e",
        8_175_917,
    ),
}

# These three candidates were nominated by the earlier exhaustive mechanical scan.
# The direct rbest preflight below independently counts the chain records relevant to
# each source/target/orientation triple before the production bridge is exercised.
RBEST_PAIR_EXPECTATIONS = {
    573: ("chr35", MappingOrientation.REVERSE),
    5170: ("chrUn_MU018764v1", MappingOrientation.SAME),
    2692: ("chrUn_JAAHUQ010000602v1", MappingOrientation.SAME),
}


class VerificationError(RuntimeError):
    """Raised when a mechanical fixture expectation is not satisfied."""


def _check(condition: bool, message: str) -> None:
    """Enforce a verifier condition even when Python optimization is enabled."""

    if not condition:
        raise VerificationError(message)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify the selected real canFam3->canFam4 mechanical fixture through "
            "liftAssess's production cached-bundle bridge. No network access is used."
        )
    )
    parser.add_argument(
        "--cache-root",
        type=Path,
        required=True,
        help=(
            "Existing liftAssess cache root containing the acquired comparative bundle"
        ),
    )
    return parser.parse_args()


def _provider_checksum(value: object) -> ProviderChecksum | None:
    if value is None:
        return None
    _check(isinstance(value, dict), "cached provider_checksum is not an object")
    value = cast(dict[str, Any], value)
    algorithm = value.get("algorithm")
    checksum = value.get("value")
    source_url = value.get("source_url")
    _check(
        isinstance(algorithm, str),
        "cached provider checksum algorithm is missing",
    )
    _check(isinstance(checksum, str), "cached provider checksum value is missing")
    _check(
        isinstance(source_url, str),
        "cached provider checksum source URL is missing",
    )
    algorithm = cast(str, algorithm)
    checksum = cast(str, checksum)
    source_url = cast(str, source_url)
    try:
        parsed_algorithm = ResourceChecksumAlgorithm(algorithm)
    except ValueError as exc:
        raise VerificationError(
            f"unsupported cached provider checksum algorithm: {algorithm!r}"
        ) from exc
    return ProviderChecksum(
        algorithm=parsed_algorithm,
        value=checksum,
        source_url=source_url,
    )


def _load_cached_resource(
    cache_root: Path,
    *,
    url: str,
    expected_sha256: str,
    expected_size: int,
) -> CachedResource:
    artifact_path = (
        cache_root
        / "artifacts"
        / "sha256"
        / expected_sha256[:2]
        / expected_sha256
    )
    if not artifact_path.is_file():
        raise SystemExit(f"missing cached artifact: {artifact_path}")
    actual_size = artifact_path.stat().st_size
    if actual_size != expected_size:
        raise SystemExit(
            f"cached artifact size mismatch for {artifact_path}: "
            f"expected {expected_size}, got {actual_size}"
        )

    # Prefer the acquisition URL-index record when this cache copy has one. A
    # content-addressed artifact can legitimately be copied without its by-url
    # index. The production bridge consumes path + recorded SHA identity, and its
    # file adapter rehashes every consumed raw stream before candidates return.
    payload: dict[str, Any] | None = None
    index_key = hashlib.sha256(url.encode("utf-8")).hexdigest()
    index_path = cache_root / "by-url" / f"{index_key}.json"
    if index_path.is_file():
        loaded: Any = json.loads(index_path.read_text(encoding="utf-8"))
        _check(isinstance(loaded, dict), f"cache index is not an object: {index_path}")
        loaded = cast(dict[str, Any], loaded)
        _check(
            loaded.get("schema_version") == 1,
            f"unexpected cache schema_version: {index_path}",
        )
        _check(
            loaded.get("source_url") == url,
            f"cache source_url mismatch: {index_path}",
        )
        _check(
            loaded.get("sha256") == expected_sha256,
            f"cache sha256 mismatch: {index_path}",
        )
        _check(
            loaded.get("size_bytes") == expected_size,
            f"cache size_bytes mismatch: {index_path}",
        )
        payload = loaded

    if payload is None:
        print(
            "cache URL index unavailable; verifying from exact SHA-addressed "
            f"artifact: {expected_sha256[:12]}…",
            flush=True,
        )
        retrieved_at = "unavailable-in-copied-cache"
        provider_checksum = None
    else:
        retrieved_at_value = payload.get("retrieved_at")
        _check(
            isinstance(retrieved_at_value, str) and bool(retrieved_at_value),
            f"cache index lacks retrieved_at: {index_path}",
        )
        retrieved_at_value = cast(str, retrieved_at_value)
        retrieved_at = retrieved_at_value
        provider_checksum = _provider_checksum(payload.get("provider_checksum"))

    return CachedResource(
        path=artifact_path,
        source_url=url,
        retrieved_at=retrieved_at,
        sha256=f"sha256:{expected_sha256}",
        size_bytes=expected_size,
        provider_checksum=provider_checksum,
        terms=ucsc_resource_terms(url),
        cache_hit=True,
    )


def _load_bundle(cache_root: Path) -> CachedUCSCResourceBundle:
    loaded = {
        role: _load_cached_resource(
            cache_root,
            url=url,
            expected_sha256=sha256,
            expected_size=size,
        )
        for role, (url, sha256, size) in RESOURCE_EXPECTATIONS.items()
    }
    return CachedUCSCResourceBundle(
        source_db=SOURCE_DB,
        target_db=TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain=loaded["chain"],
        net=loaded["net"],
        syntenic_net=loaded["syntenic_net"],
        reciprocal_best_chain=loaded["reciprocal_best_chain"],
        reciprocal_best_net=loaded["reciprocal_best_net"],
    )


def _direct_rbest_pair_counts(path: Path) -> dict[int, int]:
    """Independently count rbest chains relevant to the three fixture candidates."""

    counts = {chain_id: 0 for chain_id in RBEST_PAIR_EXPECTATIONS}
    print(
        "Preflight: scanning the small reciprocal-best chain for audit counts...",
        flush=True,
    )
    for chain in iter_chain_file(path):
        if chain.target_name != SOURCE_SEQUENCE:
            continue
        for chain_id, (target_sequence, orientation) in RBEST_PAIR_EXPECTATIONS.items():
            if chain.query_name == target_sequence and chain.orientation is orientation:
                counts[chain_id] += 1
    for chain_id in sorted(counts):
        print(
            f"  direct_rbest_pair_count chain={chain_id}: {counts[chain_id]}",
            flush=True,
        )
    return counts


def _candidate_for_chain(
    candidates: tuple[NormalizedCandidate, ...], chain_id: int
) -> NormalizedCandidate:
    suffix = f":chain:{chain_id}"
    matches = [
        candidate
        for candidate in candidates
        if candidate.candidate_id.endswith(suffix)
    ]
    _check(
        len(matches) == 1,
        f"expected exactly one candidate for chain {chain_id}, got {len(matches)}",
    )
    return matches[0]


def _observations(
    candidate: NormalizedCandidate, kind: EvidenceKind
) -> list[EvidenceObservation]:
    return [
        observation
        for observation in candidate.evidence
        if observation.kind is kind
    ]


def _single_observation(
    candidate: NormalizedCandidate, kind: EvidenceKind
) -> EvidenceObservation:
    values = _observations(candidate, kind)
    _check(
        len(values) == 1,
        f"{candidate.candidate_id}: expected one {kind.value} observation, "
        f"got {len(values)}",
    )
    return values[0]


def _assert_file_sha256(provenance: ProvenanceSource, expected: str) -> None:
    identifiers = [
        item.value
        for item in provenance.identifiers
        if item.kind is ProvenanceIdentifierKind.SHA256
    ]
    _check(
        identifiers == [f"sha256:{expected}"],
        f"unexpected SHA-256 provenance identifiers: {identifiers!r}",
    )


def _coverage(candidate: NormalizedCandidate) -> MappingCoverageSummary:
    value = _single_observation(candidate, EvidenceKind.MAPPING_COVERAGE).value
    _check(
        isinstance(value, MappingCoverageSummary),
        f"{candidate.candidate_id}: mapping coverage has unexpected value type",
    )
    return cast(MappingCoverageSummary, value)


def _rbest(candidate: NormalizedCandidate) -> ReciprocalBestMembershipSummary:
    value = _single_observation(
        candidate, EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
    ).value
    _check(
        isinstance(value, ReciprocalBestMembershipSummary),
        f"{candidate.candidate_id}: reciprocal-best has unexpected value type",
    )
    return cast(ReciprocalBestMembershipSummary, value)


def main() -> None:
    args = parse_args()
    bundle = _load_bundle(args.cache_root)

    # COMPLETE_RESOURCE remains an explicit caller trust claim for a manually
    # constructed bundle. This verifier relies on the previously completed UCSC
    # acquisition for completeness; it independently checks exact byte size here,
    # while the production file-backed path checks SHA-256 over each consumed stream.
    direct_rbest_counts = _direct_rbest_pair_counts(bundle.reciprocal_best_chain.path)

    source = AssemblyIdentifier(name=SOURCE_DB, provider="UCSC")
    target = AssemblyIdentifier(name=TARGET_DB, provider="UCSC")
    alignment = ProvenanceSource(
        source_id="fixture:ucsc:canFam3-canFam4:shared-alignment",
        label=(
            "caller-declared shared upstream UCSC canFam3->canFam4 comparative "
            "alignment for mechanical fixture verification"
        ),
    )
    interval = GenomicInterval(
        assembly=source,
        sequence_name=SOURCE_SEQUENCE,
        start=SOURCE_START,
        end=SOURCE_END,
    )

    print(
        "Running selected locus through the full cached-bundle engine path...",
        flush=True,
    )
    candidates = build_ucsc_candidates_from_cached_bundle(
        interval,
        bundle,
        target_assembly=target,
        alignment_provenance=alignment,
    )

    _check(len(candidates) == 170, f"expected 170 candidates, got {len(candidates)}")
    target_sequence_count = len(
        {candidate.target_interval.sequence_name for candidate in candidates}
    )
    _check(
        target_sequence_count == 114,
        f"expected 114 distinct target sequences, got {target_sequence_count}",
    )

    primary = _candidate_for_chain(candidates, 573)
    _check(primary.target_interval.sequence_name == "chr35", "chain 573 target sequence")
    _check(
        (primary.target_interval.start, primary.target_interval.end) == (925_644, 925_938),
        "chain 573 target interval",
    )
    _check(primary.orientation is MappingOrientation.REVERSE, "chain 573 orientation")
    _check(len(primary.segments) == 2, "chain 573 segment count")
    _check(
        _single_observation(primary, EvidenceKind.CHAIN_SCORE).value == 16_617_372.0,
        "chain 573 score",
    )

    primary_coverage = _coverage(primary)
    _check(primary_coverage.status is MappingCoverageStatus.FULL, "chain 573 coverage status")
    _check(primary_coverage.covered_source_bases == 100, "chain 573 covered source bases")
    _check(primary_coverage.source_bases == 100, "chain 573 source bases")

    primary_gaps = _single_observation(primary, EvidenceKind.CHAIN_GAPS).value
    _check(isinstance(primary_gaps, ChainGapSummary), "chain 573 gap summary type")
    primary_gaps = cast(ChainGapSummary, primary_gaps)
    _check(len(primary_gaps.gaps) == 1, "chain 573 gap count")

    _check(
        [item.value for item in _observations(primary, EvidenceKind.ALIGNED_BASES)] == [3603],
        "chain 573 ali",
    )
    _check(
        [item.value for item in _observations(primary, EvidenceKind.DUPLICATED_QUERY_BASES)]
        == [4098],
        "chain 573 qDup",
    )
    _check(
        [item.value for item in _observations(primary, EvidenceKind.NET_CLASSIFICATION)]
        == ["nonSyn"],
        "chain 573 net classification",
    )

    hierarchy_observation = _single_observation(primary, EvidenceKind.NET_HIERARCHY)
    hierarchy = hierarchy_observation.value
    _check(isinstance(hierarchy, NetHierarchySummary), "chain 573 net hierarchy type")
    hierarchy = cast(NetHierarchySummary, hierarchy)
    _check(hierarchy.depth == 7, "chain 573 net hierarchy depth")
    _check(
        hierarchy.source_fill_interval.sequence_name == SOURCE_SEQUENCE,
        "chain 573 net source sequence",
    )
    _check(
        (hierarchy.source_fill_interval.start, hierarchy.source_fill_interval.end)
        == (1_843_971, 1_847_599),
        "chain 573 net source fill interval",
    )

    primary_rbest_observation = _single_observation(
        primary, EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
    )
    primary_rbest = primary_rbest_observation.value
    _check(
        isinstance(primary_rbest, ReciprocalBestMembershipSummary),
        "chain 573 reciprocal-best summary type",
    )
    primary_rbest = cast(ReciprocalBestMembershipSummary, primary_rbest)
    _check(primary_rbest.status is ReciprocalBestMembershipStatus.FULL, "chain 573 rbest status")
    _check(
        primary_rbest.resource_completeness
        is ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE,
        "chain 573 rbest completeness",
    )
    _check(primary_rbest.covered_source_bases == 100, "chain 573 rbest covered bases")
    _check(primary_rbest.candidate_source_bases == 100, "chain 573 rbest candidate bases")
    _check(
        primary_rbest.chains_examined == direct_rbest_counts[573],
        "chain 573 chains_examined disagrees with direct rbest pair count",
    )

    partial = _candidate_for_chain(candidates, 5170)
    _check(partial.target_interval.sequence_name == "chrUn_MU018764v1", "chain 5170 target sequence")
    _check(
        (partial.target_interval.start, partial.target_interval.end) == (171_661, 171_760),
        "chain 5170 target interval",
    )
    _check(partial.orientation is MappingOrientation.SAME, "chain 5170 orientation")
    _check(len(partial.segments) == 2, "chain 5170 segment count")
    partial_coverage = _coverage(partial)
    _check(partial_coverage.status is MappingCoverageStatus.PARTIAL, "chain 5170 coverage status")
    _check(partial_coverage.covered_source_bases == 99, "chain 5170 covered source bases")
    _check(partial_coverage.source_bases == 100, "chain 5170 source bases")
    _check(
        len(partial_coverage.uncovered_source_intervals) == 1,
        "chain 5170 uncovered interval count",
    )
    partial_gaps = _single_observation(partial, EvidenceKind.CHAIN_GAPS).value
    _check(isinstance(partial_gaps, ChainGapSummary), "chain 5170 gap summary type")
    partial_gaps = cast(ChainGapSummary, partial_gaps)
    _check(len(partial_gaps.gaps) == 1, "chain 5170 gap count")
    partial_rbest = _rbest(partial)
    _check(partial_rbest.status is ReciprocalBestMembershipStatus.NONE, "chain 5170 rbest status")
    _check(partial_rbest.covered_source_bases == 0, "chain 5170 rbest covered bases")
    _check(partial_rbest.candidate_source_bases == 99, "chain 5170 rbest candidate bases")
    _check(
        partial_rbest.chains_examined == direct_rbest_counts[5170],
        "chain 5170 chains_examined disagrees with direct rbest pair count",
    )

    alternative = _candidate_for_chain(candidates, 2692)
    _check(
        alternative.target_interval.sequence_name == "chrUn_JAAHUQ010000602v1",
        "chain 2692 target sequence",
    )
    _check(
        (alternative.target_interval.start, alternative.target_interval.end) == (62_326, 62_622),
        "chain 2692 target interval",
    )
    _check(alternative.orientation is MappingOrientation.SAME, "chain 2692 orientation")
    _check(len(alternative.segments) == 2, "chain 2692 segment count")
    alternative_coverage = _coverage(alternative)
    _check(
        alternative_coverage.status is MappingCoverageStatus.FULL,
        "chain 2692 coverage status",
    )
    alternative_gaps = _single_observation(alternative, EvidenceKind.CHAIN_GAPS).value
    _check(isinstance(alternative_gaps, ChainGapSummary), "chain 2692 gap summary type")
    alternative_gaps = cast(ChainGapSummary, alternative_gaps)
    _check(len(alternative_gaps.gaps) == 1, "chain 2692 gap count")
    alternative_rbest = _rbest(alternative)
    _check(alternative_rbest.status is ReciprocalBestMembershipStatus.NONE, "chain 2692 rbest status")
    _check(alternative_rbest.covered_source_bases == 0, "chain 2692 rbest covered bases")
    _check(alternative_rbest.candidate_source_bases == 100, "chain 2692 rbest candidate bases")
    _check(
        alternative_rbest.chains_examined == direct_rbest_counts[2692],
        "chain 2692 chains_examined disagrees with direct rbest pair count",
    )

    chain_sha = RESOURCE_EXPECTATIONS["chain"][1]
    net_sha = RESOURCE_EXPECTATIONS["net"][1]
    rbest_sha = RESOURCE_EXPECTATIONS["reciprocal_best_chain"][1]
    _assert_file_sha256(primary.mapping_provenance, chain_sha)
    _check(
        primary.mapping_provenance.derived_from == (alignment,),
        "chain file provenance did not preserve caller-declared alignment ancestor",
    )

    net_file_provenance = hierarchy_observation.provenance.derived_from
    _check(len(net_file_provenance) == 1, "net fill provenance should have one file parent")
    _assert_file_sha256(net_file_provenance[0], net_sha)
    _check(
        net_file_provenance[0].derived_from == (alignment,),
        "net file provenance did not preserve caller-declared alignment ancestor",
    )

    _assert_file_sha256(primary_rbest_observation.provenance, rbest_sha)
    _check(
        primary_rbest_observation.provenance.derived_from == (alignment,),
        "rbest file provenance did not preserve caller-declared alignment ancestor",
    )

    print()
    print("=== MECHANICAL FIXTURE VERIFICATION: PASS ===")
    print(f"source_0based_half_open={SOURCE_SEQUENCE}:{SOURCE_START}-{SOURCE_END}")
    print(f"candidate_count={len(candidates)}")
    print(f"distinct_target_sequences={target_sequence_count}")
    print(
        "primary_chain_573="
        f"{primary.target_interval.sequence_name}:{primary.target_interval.start}-"
        f"{primary.target_interval.end} orientation={primary.orientation.value} "
        f"segments={len(primary.segments)} coverage=FULL:100/100 gaps=1 "
        "score=16617372 ali=3603 qDup=4098 net=nonSyn depth=7 "
        f"rbest=FULL:100/100 chains_examined={primary_rbest.chains_examined}"
    )
    print(
        "partial_chain_5170="
        f"{partial.target_interval.sequence_name}:{partial.target_interval.start}-"
        f"{partial.target_interval.end} coverage=PARTIAL:99/100 gaps=1 "
        f"rbest=NONE:0/99 chains_examined={partial_rbest.chains_examined}"
    )
    print(
        "alternative_chain_2692="
        f"{alternative.target_interval.sequence_name}:"
        f"{alternative.target_interval.start}-{alternative.target_interval.end} "
        "coverage=FULL:100/100 gaps=1 "
        f"rbest=NONE:0/100 chains_examined={alternative_rbest.chains_examined}"
    )
    print(
        "provenance_wiring=chain/net/rbest file nodes preserve one "
        "caller-declared shared alignment ancestor"
    )
    print(
        "resource_completeness_basis=caller-declared COMPLETE_RESOURCE from prior "
        "complete acquisition; exact consumed bytes SHA-verified by production path"
    )
    print("verdict_computed=no")
    print("biological_ground_truth_claim=no")
    print()
    print("Exact primary segment/gap/rbest geometry for fixture freezing:")
    for index, segment in enumerate(primary.segments, start=1):
        print(
            f"  segment_{index}: source={segment.source_interval.start}-"
            f"{segment.source_interval.end} target={segment.target_interval.start}-"
            f"{segment.target_interval.end}"
        )
    for index, gap in enumerate(primary_gaps.gaps, start=1):
        source_gap = gap.source_gap_overlap
        target_gap = gap.target_gap_interval
        source_text = (
            "none"
            if source_gap is None
            else f"{source_gap.sequence_name}:{source_gap.start}-{source_gap.end}"
        )
        target_text = (
            "none"
            if target_gap is None
            else f"{target_gap.sequence_name}:{target_gap.start}-{target_gap.end}"
        )
        print(
            f"  gap_{index}: source_boundary={gap.source_boundary} "
            f"source_gap={source_text} target_gap={target_text}"
        )
    for index, covered in enumerate(primary_rbest.covered_source_intervals, start=1):
        print(
            f"  rbest_covered_{index}: {covered.sequence_name}:"
            f"{covered.start}-{covered.end}"
        )


if __name__ == "__main__":
    main()
