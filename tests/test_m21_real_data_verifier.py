from __future__ import annotations

import runpy
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest


def _load_verifier() -> dict[str, Any]:
    return runpy.run_path("scripts/verify_m21_b12_b14.py")


def _resource(role: str, sha256: str) -> dict[str, object]:
    return {"role": role, "sha256": sha256}


def _candidate_id(chain_id: int) -> str:
    return f"canFam3:canFam4:chain:{chain_id}"


def _payload_for_case(namespace: dict[str, Any], case: Any) -> dict[str, object]:
    all_ids = [_candidate_id(chain_id) for chain_id in case.expected_chain_ids]
    filtered_ids = all_ids[: case.filtered_chain_count]
    additional_ids = all_ids[case.filtered_chain_count :]

    support: list[dict[str, object]] = []
    for chain_id in case.expected_chain_ids:
        item: dict[str, object] = {
            "candidate_id": _candidate_id(chain_id),
            "complete_source_coverage": True,
            "retained_by_filtered_chain": chain_id
            in case.expected_chain_ids[: case.filtered_chain_count],
            "depth1_top_net": False,
            "full_reciprocal_best": False,
        }
        if case.label == "B14" and chain_id == 4:
            item["depth1_top_net"] = True
            item["full_reciprocal_best"] = True
        support.append(item)

    favored_id = (
        None if case.favored_chain_id is None else _candidate_id(case.favored_chain_id)
    )

    return {
        "resources": [
            _resource("CHAIN", namespace["ALL_CHAIN_SHA256"]),
            _resource("NET", namespace["NET_SHA256"]),
            _resource("RECIPROCAL_BEST_CHAIN", namespace["RBEST_CHAIN_SHA256"]),
        ],
        "filtered_all_chain_comparison": {
            "assessed": True,
            "inventory_state": case.inventory_state,
            "categorical_relationship": case.relationship,
            "all_chain_candidate_ids": all_ids,
            "filtered_candidate_ids": filtered_ids,
            "additional_all_chain_candidate_ids": additional_ids,
            "favored_candidate_id": favored_id,
            "filtered_chain_resource": {"sha256": namespace["FILTERED_CHAIN_SHA256"]},
            "provenance": {
                "shared_processing_run_provenance_verified": False,
            },
        },
        "result_profile": {
            "headline": case.headline,
            "projection_count": case.projection_count,
            "comparative_relationship": {
                "state": case.relationship,
                "placement_support": support,
            },
            "query_context": {
                "check_state": "RUN",
                "findings": sorted(case.context_findings),
            },
        },
    }


def test_m21_real_data_verifier_loads_recorded_cases() -> None:
    namespace = _load_verifier()

    assert "main" in namespace
    cases = namespace["CASES"]
    assert tuple(case.label for case in cases) == ("B12", "B13", "B14")
    assert tuple(case.locus for case in cases) == (
        "chrX:26956239-26956239",
        "chr28:1484906-1484906",
        "chr5:31705136-31705136",
    )
    assert cases[2].all_chain_count == 9
    assert cases[2].filtered_chain_count == 1
    assert cases[2].relationship == "FAVORS_ONE_PLACEMENT"


def test_chain_id_parses_and_rejects_invalid_suffixes() -> None:
    namespace = _load_verifier()
    chain_id = namespace["_chain_id"]
    verification_error = namespace["VerificationError"]

    assert chain_id("canFam3:canFam4:chain:4") == 4
    with pytest.raises(verification_error, match="lacks chain suffix"):
        chain_id("candidate-4")
    with pytest.raises(verification_error, match="non-integer chain suffix"):
        chain_id("canFam3:canFam4:chain:not-an-int")


@pytest.mark.parametrize("case_index", (0, 1, 2))
def test_verify_case_accepts_recorded_synthetic_payload(case_index: int) -> None:
    namespace = _load_verifier()
    case = namespace["CASES"][case_index]
    payload = _payload_for_case(namespace, case)

    namespace["_verify_case"](payload, case)


def test_verify_resource_identity_rejects_changed_filtered_chain() -> None:
    namespace = _load_verifier()
    case = namespace["CASES"][0]
    payload = _payload_for_case(namespace, case)
    comparison = payload["filtered_all_chain_comparison"]
    assert isinstance(comparison, dict)
    filtered_resource = comparison["filtered_chain_resource"]
    assert isinstance(filtered_resource, dict)
    filtered_resource["sha256"] = "sha256:" + "0" * 64

    with pytest.raises(
        namespace["VerificationError"],
        match="ordinary filtered-chain fixture identity changed",
    ):
        namespace["_verify_resource_identity"](payload, case.label)


def test_verify_case_rejects_changed_relationship() -> None:
    namespace = _load_verifier()
    case = namespace["CASES"][2]
    payload = _payload_for_case(namespace, case)
    comparison = payload["filtered_all_chain_comparison"]
    assert isinstance(comparison, dict)
    comparison["categorical_relationship"] = "DOES_NOT_SEPARATE_PLACEMENTS"

    with pytest.raises(
        namespace["VerificationError"],
        match="categorical relationship changed",
    ):
        namespace["_verify_case"](payload, case)


def test_prepared_index_preflight_requires_both_exact_indexes(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    namespace = _load_verifier()
    verify_indexes = namespace["_verify_prepared_indexes"]
    function_globals = verify_indexes.__globals__
    evidence_tier = namespace["EvidenceAvailabilityTier"]
    all_sha = namespace["ALL_CHAIN_SHA256"]
    filtered_sha = namespace["FILTERED_CHAIN_SHA256"]
    resources = {
        evidence_tier.COMPARATIVE: SimpleNamespace(
            chain=SimpleNamespace(sha256=all_sha)
        ),
        evidence_tier.LIFTOVER_ONLY: SimpleNamespace(
            chain=SimpleNamespace(sha256=filtered_sha)
        ),
    }
    loaded: list[str] = []

    def fake_resolve(
        cache_root: Path,
        source_db: str,
        target_db: str,
        *,
        evidence_tier: Any,
    ) -> Any:
        assert cache_root == tmp_path
        assert source_db == "canFam3"
        assert target_db == "canFam4"
        return resources[evidence_tier]

    def fake_load(
        cache_root: Path,
        resource: Any,
        *,
        verify_database: bool,
    ) -> object:
        assert cache_root == tmp_path
        assert verify_database is False
        loaded.append(resource.sha256)
        return object()

    monkeypatch.setitem(
        function_globals,
        "resolve_cached_ucsc_chain_resource_metadata",
        fake_resolve,
    )
    monkeypatch.setitem(function_globals, "load_cached_chain_index", fake_load)

    verify_indexes(tmp_path)

    assert loaded == [all_sha, filtered_sha]


def test_prepared_index_preflight_rejects_missing_index(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    namespace = _load_verifier()
    verify_index = namespace["_verify_prepared_chain_index"]
    function_globals = verify_index.__globals__
    evidence_tier = namespace["EvidenceAvailabilityTier"]
    expected_sha = namespace["ALL_CHAIN_SHA256"]

    monkeypatch.setitem(
        function_globals,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: SimpleNamespace(
            chain=SimpleNamespace(sha256=expected_sha)
        ),
    )
    monkeypatch.setitem(
        function_globals,
        "load_cached_chain_index",
        lambda *args, **kwargs: None,
    )

    with pytest.raises(
        namespace["VerificationError"],
        match="COMPARATIVE: prepared chain index is missing",
    ):
        verify_index(tmp_path, evidence_tier.COMPARATIVE, expected_sha)
