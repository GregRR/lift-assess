from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

SOURCE_DB = "canFam3"
TARGET_DB = "canFam4"

ALL_CHAIN_SHA256 = (
    "sha256:f10a6b48b5461bb8378ffaff311fb7355b1910511131ce0f5df5402c4db67519"
)
FILTERED_CHAIN_SHA256 = (
    "sha256:c79c9e7c2a3d546f7a9d7efe27cc8815da611d79adb0da4e4ff1556810f28f48"
)
NET_SHA256 = "sha256:c889134e95ff82741c0092b1673b3e5fe0125aa82f611ee97d91e468b56d51ac"
RBEST_CHAIN_SHA256 = (
    "sha256:34f4061fd29e7720c7eb2adc1ea8299e86f21f08e18f97f9ba468cf8b466690c"
)


class VerificationError(RuntimeError):
    """Raised when the recorded M21 real-data expectation is not satisfied."""


@dataclass(frozen=True)
class CaseExpectation:
    label: str
    locus: str
    headline: str
    projection_count: str
    inventory_state: str
    relationship: str
    all_chain_count: int
    filtered_chain_count: int
    additional_count: int
    favored_chain_id: int | None
    expected_chain_ids: tuple[int, ...]
    context_findings: frozenset[str]


CASES = (
    CaseExpectation(
        label="B12",
        locus="chrX:26956239-26956239",
        headline="ONE_COMPLETE_CHAIN_PROJECTION",
        projection_count="ONE",
        inventory_state="FILTERED_AND_ALL_CHAIN_AGREE",
        relationship="NO_COMPETING_FULL_PLACEMENTS",
        all_chain_count=1,
        filtered_chain_count=1,
        additional_count=0,
        favored_chain_id=None,
        expected_chain_ids=(2,),
        context_findings=frozenset({"AGREES_WITH_POINT"}),
    ),
    CaseExpectation(
        label="B13",
        locus="chr28:1484906-1484906",
        headline="ONE_COMPLETE_CHAIN_PROJECTION",
        projection_count="ONE",
        inventory_state="FILTERED_AND_ALL_CHAIN_AGREE",
        relationship="NO_COMPETING_FULL_PLACEMENTS",
        all_chain_count=1,
        filtered_chain_count=1,
        additional_count=0,
        favored_chain_id=None,
        expected_chain_ids=(30,),
        context_findings=frozenset({"AGREES_WITH_POINT"}),
    ),
    CaseExpectation(
        label="B14",
        locus="chr5:31705136-31705136",
        headline="MULTIPLE_CHAIN_PROJECTIONS",
        projection_count="MULTIPLE",
        inventory_state="ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS",
        relationship="FAVORS_ONE_PLACEMENT",
        all_chain_count=9,
        filtered_chain_count=1,
        additional_count=8,
        favored_chain_id=4,
        expected_chain_ids=(
            4,
            8_870_537,
            11_433_457,
            12_687_949,
            13_039_430,
            13_503_555,
            13_506_928,
            13_601_144,
            16_339_947,
        ),
        context_findings=frozenset(
            {
                "REVEALS_PARTIAL_COVERAGE",
                "REVEALS_FRAGMENTATION",
                "REVEALS_TARGET_DISCONTINUITY",
                "CHANGES_WITH_QUERY_SCALE",
            }
        ),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Verify the recorded M21 B12-B14 canFam3->canFam4 real-data outcomes "
            "through the production CLI. The verifier is offline and requires "
            "the comparative and filtered chain indexes to be prepared already."
        )
    )
    parser.add_argument(
        "--cache-root",
        type=Path,
        required=True,
        help="Existing liftAssess cache root containing the recorded UCSC resources.",
    )
    return parser.parse_args()


def _check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def _object(value: object, description: str) -> dict[str, Any]:
    _check(isinstance(value, dict), f"{description} is not an object")
    return cast(dict[str, Any], value)


def _array(value: object, description: str) -> list[Any]:
    _check(isinstance(value, list), f"{description} is not an array")
    return cast(list[Any], value)


def _chain_id(candidate_id: object) -> int:
    _check(isinstance(candidate_id, str), "candidate ID is not a string")
    candidate_id = cast(str, candidate_id)
    marker = ":chain:"
    _check(marker in candidate_id, f"candidate ID lacks chain suffix: {candidate_id}")
    _, chain_text = candidate_id.rsplit(marker, maxsplit=1)
    try:
        return int(chain_text)
    except ValueError as exc:
        raise VerificationError(
            f"candidate ID has non-integer chain suffix: {candidate_id}"
        ) from exc


def _run_case(cache_root: Path, case: CaseExpectation) -> dict[str, Any]:
    command = [
        sys.executable,
        "-m",
        "liftassess.cli",
        SOURCE_DB,
        TARGET_DB,
        case.locus,
        "--cache-dir",
        str(cache_root),
        "--evidence-tier",
        "COMPARATIVE",
        "--offline",
        "--json",
        "--quiet",
    ]
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise VerificationError(
            f"{case.label}: CLI failed with exit {completed.returncode}: "
            f"{completed.stderr.strip()}"
        )
    try:
        payload: object = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise VerificationError(f"{case.label}: CLI did not emit valid JSON") from exc
    return _object(payload, f"{case.label} JSON result")


def _resource_sha(payload: dict[str, Any], role: str) -> str:
    resources = _array(payload.get("resources"), "resources")
    matches = [
        _object(item, "resource")
        for item in resources
        if isinstance(item, dict) and item.get("role") == role
    ]
    _check(len(matches) == 1, f"expected exactly one {role} resource")
    sha = matches[0].get("sha256")
    _check(isinstance(sha, str), f"{role} resource lacks SHA-256")
    return cast(str, sha)


def _verify_resource_identity(payload: dict[str, Any], label: str) -> None:
    _check(
        _resource_sha(payload, "CHAIN") == ALL_CHAIN_SHA256,
        f"{label}: comparative all-chain fixture identity changed",
    )
    _check(
        _resource_sha(payload, "NET") == NET_SHA256,
        f"{label}: net fixture identity changed",
    )
    _check(
        _resource_sha(payload, "RECIPROCAL_BEST_CHAIN") == RBEST_CHAIN_SHA256,
        f"{label}: reciprocal-best fixture identity changed",
    )

    comparison = _object(
        payload.get("filtered_all_chain_comparison"),
        "filtered/all-chain comparison",
    )
    filtered_resource = _object(
        comparison.get("filtered_chain_resource"),
        "filtered-chain comparison resource",
    )
    _check(
        filtered_resource.get("sha256") == FILTERED_CHAIN_SHA256,
        f"{label}: ordinary filtered-chain fixture identity changed",
    )


def _verify_case(payload: dict[str, Any], case: CaseExpectation) -> None:
    profile = _object(payload.get("result_profile"), f"{case.label} result_profile")
    _check(profile.get("headline") == case.headline, f"{case.label}: headline changed")
    _check(
        profile.get("projection_count") == case.projection_count,
        f"{case.label}: projection-count state changed",
    )

    comparison = _object(
        payload.get("filtered_all_chain_comparison"),
        f"{case.label} filtered/all-chain comparison",
    )
    _check(comparison.get("assessed") is True, f"{case.label}: comparison not assessed")
    _check(
        comparison.get("inventory_state") == case.inventory_state,
        f"{case.label}: inventory state changed",
    )
    _check(
        comparison.get("categorical_relationship") == case.relationship,
        f"{case.label}: categorical relationship changed",
    )

    all_chain_ids = _array(
        comparison.get("all_chain_candidate_ids"),
        f"{case.label} all-chain candidate IDs",
    )
    filtered_ids = _array(
        comparison.get("filtered_candidate_ids"),
        f"{case.label} filtered candidate IDs",
    )
    additional_ids = _array(
        comparison.get("additional_all_chain_candidate_ids"),
        f"{case.label} additional all-chain candidate IDs",
    )
    _check(
        len(all_chain_ids) == case.all_chain_count,
        f"{case.label}: all-chain placement count changed",
    )
    _check(
        len(filtered_ids) == case.filtered_chain_count,
        f"{case.label}: filtered placement count changed",
    )
    _check(
        len(additional_ids) == case.additional_count,
        f"{case.label}: additional placement count changed",
    )

    actual_chain_ids = tuple(_chain_id(candidate_id) for candidate_id in all_chain_ids)
    _check(
        set(actual_chain_ids) == set(case.expected_chain_ids),
        f"{case.label}: all-chain candidate identities changed: {actual_chain_ids!r}",
    )

    favored_id = comparison.get("favored_candidate_id")
    if case.favored_chain_id is None:
        _check(favored_id is None, f"{case.label}: unexpectedly favors a placement")
    else:
        _check(
            _chain_id(favored_id) == case.favored_chain_id,
            f"{case.label}: favored placement changed",
        )

    relationship_profile = _object(
        profile.get("comparative_relationship"),
        f"{case.label} comparative relationship profile",
    )
    _check(
        relationship_profile.get("state") == case.relationship,
        f"{case.label}: profile relationship disagrees with comparison",
    )
    support = _array(
        relationship_profile.get("placement_support"),
        f"{case.label} placement support",
    )
    _check(
        len(support) == case.all_chain_count,
        f"{case.label}: placement-support count changed",
    )

    if case.label == "B14":
        support_by_chain = {
            _chain_id(_object(item, "placement support").get("candidate_id")): _object(
                item, "placement support"
            )
            for item in support
        }
        favored = support_by_chain[4]
        _check(
            favored.get("complete_source_coverage") is True,
            "B14: chain 4 not complete",
        )
        _check(
            favored.get("retained_by_filtered_chain") is True,
            "B14: chain 4 not filtered-retained",
        )
        _check(
            favored.get("depth1_top_net") is True,
            "B14: chain 4 lost depth-1 top-net support",
        )
        _check(
            favored.get("full_reciprocal_best") is True,
            "B14: chain 4 lost full rbest support",
        )
        for chain_id, item in support_by_chain.items():
            if chain_id == 4:
                continue
            _check(
                item.get("complete_source_coverage") is True,
                f"B14: chain {chain_id} not complete",
            )
            _check(
                item.get("retained_by_filtered_chain") is False,
                f"B14: chain {chain_id} unexpectedly filtered-retained",
            )
            _check(
                item.get("depth1_top_net") is False,
                f"B14: chain {chain_id} unexpectedly depth-1 top-net",
            )
            _check(
                item.get("full_reciprocal_best") is False,
                f"B14: chain {chain_id} unexpectedly full rbest",
            )

    query_context = _object(
        profile.get("query_context"), f"{case.label} query-context profile"
    )
    _check(
        query_context.get("check_state") == "RUN",
        f"{case.label}: context did not run",
    )
    findings = _array(query_context.get("findings"), f"{case.label} context findings")
    _check(
        frozenset(cast(list[str], findings)) == case.context_findings,
        f"{case.label}: context findings changed: {findings!r}",
    )

    provenance = _object(
        comparison.get("provenance"),
        f"{case.label} comparison provenance",
    )
    _check(
        provenance.get("shared_processing_run_provenance_verified") is False,
        f"{case.label}: comparison unexpectedly claims verified shared processing run",
    )

    _verify_resource_identity(payload, case.label)


def main() -> None:
    args = parse_args()
    cache_root = args.cache_root.expanduser().resolve()
    _check(cache_root.is_dir(), f"cache root does not exist: {cache_root}")

    print("M21 B12-B14 real-data verification")
    print(f"cache_root={cache_root}")
    print("network_access=no")
    print("index_builds=no")

    for case in CASES:
        print(f"verifying {case.label}: {case.locus} ...", flush=True)
        payload = _run_case(cache_root, case)
        _verify_case(payload, case)
        comparison = _object(payload["filtered_all_chain_comparison"], "comparison")
        print(
            f"  PASS {case.label}: "
            f"all_chain={len(comparison['all_chain_candidate_ids'])} "
            f"filtered={len(comparison['filtered_candidate_ids'])} "
            f"relationship={comparison['categorical_relationship']}"
        )

    print("M21 B12-B14 verification: PASS")


if __name__ == "__main__":
    try:
        main()
    except VerificationError as exc:
        print(f"verification failed: {exc}", file=sys.stderr)
        print(
            "Prerequisite: the recorded canFam3->canFam4 COMPARATIVE and "
            "LIFTOVER-ONLY resources must be cached, and both chain indexes must "
            "already be prepared. This verifier never downloads or builds indexes.",
            file=sys.stderr,
        )
        raise SystemExit(1) from exc
