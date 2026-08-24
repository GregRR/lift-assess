from __future__ import annotations

import argparse
import gzip
from io import StringIO
from pathlib import Path

import pytest

from liftassess import (
    CachedResource,
    CachedUCSCResourceBundle,
    EvidenceAvailabilityTier,
    index_cli,
    sha256_identifier_for_file,
    ucsc_resource_terms,
)
from liftassess.chain_index import load_cached_chain_index

_SOURCE_DB = "canFam3"
_TARGET_DB = "canFam4"
_CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
    "canFam3ToCanFam4.over.chain.gz"
)


def _cached_bundle(tmp_path: Path) -> CachedUCSCResourceBundle:
    chain_path = tmp_path / "chain.gz"
    chain = "chain 100 chr1 1000 + 100 120 chrA 2000 + 500 520 1\n20\n\n"
    with gzip.open(chain_path, mode="wt", encoding="ascii", newline="\n") as handle:
        handle.write(chain)
    resource = CachedResource(
        path=chain_path,
        source_url=_CHAIN_URL,
        retrieved_at="2026-08-20T00:00:00Z",
        sha256=sha256_identifier_for_file(chain_path).value,
        size_bytes=chain_path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(_CHAIN_URL),
        cache_hit=True,
    )
    return CachedUCSCResourceBundle(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=resource,
    )


def _args(cache_root: Path, *extra: str) -> argparse.Namespace:
    return index_cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "--cache-dir", str(cache_root), *extra]
    )


def test_prepare_index_requires_verified_cached_bundle(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        index_cli,
        "load_cached_ucsc_resource_bundle",
        lambda *args, **kwargs: None,
    )
    stderr = StringIO()

    exit_code = index_cli._run(
        _args(tmp_path / "cache"),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 1
    assert "requires a complete verified cached UCSC bundle" in stderr.getvalue()
    assert "assess-liftover" in stderr.getvalue()


def test_prepare_index_builds_and_reuses_exact_cached_chain(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cache_root = tmp_path / "cache"
    cache_root.mkdir()
    bundle = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        index_cli,
        "load_cached_ucsc_resource_bundle",
        lambda *args, **kwargs: bundle,
    )

    stdout = StringIO()
    exit_code = index_cli._run(
        _args(cache_root, "--quiet"),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert "Prepared: reusable chain index" in stdout.getvalue()
    index = load_cached_chain_index(cache_root, bundle.chain)
    assert index is not None
    assert index.manifest.source_chain_sha256_identifier == bundle.chain.sha256

    stdout = StringIO()
    exit_code = index_cli._run(
        _args(cache_root, "--quiet"),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert "Already prepared: reusable chain index" in stdout.getvalue()


def test_prepare_index_reports_corruption_and_requires_explicit_rebuild(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cache_root = tmp_path / "cache"
    cache_root.mkdir()
    bundle = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        index_cli,
        "load_cached_ucsc_resource_bundle",
        lambda *args, **kwargs: bundle,
    )
    assert (
        index_cli._run(
            _args(cache_root, "--quiet"),
            stdout=StringIO(),
            stderr=StringIO(),
        )
        == 0
    )
    index = load_cached_chain_index(cache_root, bundle.chain)
    assert index is not None
    with index.database_path.open("r+b") as handle:
        handle.seek(100)
        original = handle.read(1)
        handle.seek(100)
        handle.write(bytes([original[0] ^ 1]))

    stderr = StringIO()
    exit_code = index_cli._run(
        _args(cache_root, "--quiet"),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 1
    assert "--rebuild" in stderr.getvalue()

    stdout = StringIO()
    exit_code = index_cli._run(
        _args(cache_root, "--quiet", "--rebuild"),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert "Prepared: reusable chain index" in stdout.getvalue()
    assert load_cached_chain_index(cache_root, bundle.chain) is not None


def test_prepare_index_can_select_exact_chain_publication_class(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    from liftassess import CachedUCSCChainResource

    cache_root = tmp_path / "cache"
    cache_root.mkdir()
    bundle = _cached_bundle(tmp_path)
    chain_context = CachedUCSCChainResource(
        source_db=bundle.source_db,
        target_db=bundle.target_db,
        evidence_tier=bundle.evidence_tier,
        chain=bundle.chain,
    )
    seen_tiers: list[EvidenceAvailabilityTier] = []

    def load_chain(*args: object, **kwargs: object) -> CachedUCSCChainResource:
        del args
        tier = kwargs["evidence_tier"]
        assert isinstance(tier, EvidenceAvailabilityTier)
        seen_tiers.append(tier)
        return chain_context

    monkeypatch.setattr(index_cli, "load_cached_ucsc_chain_resource", load_chain)
    monkeypatch.setattr(
        index_cli,
        "load_cached_ucsc_resource_bundle",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("exact-tier preparation must not require a full bundle")
        ),
    )

    stdout = StringIO()
    exit_code = index_cli._run(
        _args(cache_root, "--evidence-tier", "LIFTOVER-ONLY", "--quiet"),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert seen_tiers == [EvidenceAvailabilityTier.LIFTOVER_ONLY]
    assert "Publication class: LIFTOVER-ONLY" in stdout.getvalue()
    assert load_cached_chain_index(cache_root, bundle.chain) is not None
