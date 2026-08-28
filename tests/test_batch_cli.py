from __future__ import annotations

import json
from io import StringIO
from pathlib import Path

import pytest

from liftassess import cli
from liftassess.chain_index import ChainIndex, build_chain_index
from liftassess.models import EvidenceAvailabilityTier
from liftassess.resource_cache import (
    CachedResource,
    CachedUCSCChainResource,
    ucsc_resource_terms,
)
from liftassess.resource_identity import sha256_identifier_for_file

SOURCE_DB = "hg38"
TARGET_DB = "hg19"


def _chain_context(tmp_path: Path) -> tuple[CachedUCSCChainResource, ChainIndex]:
    path = tmp_path / "hg38ToHg19.over.chain"
    path.write_text(
        """\
chain 100 chr1 1000 + 0 10 chrA 2000 + 100 110 1
10

chain 90 chr1 1000 + 20 30 chrA 2000 + 100 110 2
10

chain 80 chr1 1000 + 40 50 chrA 2000 + 105 115 3
10

""",
        encoding="ascii",
    )
    identifier = sha256_identifier_for_file(path).value
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
        "hg38ToHg19.over.chain.gz"
    )
    resource = CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-28T00:00:00Z",
        sha256=identifier,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=True,
    )
    context = CachedUCSCChainResource(
        source_db=SOURCE_DB,
        target_db=TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=resource,
    )
    index = build_chain_index(
        path,
        tmp_path / "index",
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=path.stat().st_size,
    ).index
    return context, index


def _bed_path(tmp_path: Path) -> Path:
    path = tmp_path / "batch.bed"
    path.write_text(
        "chr1\t0\t10\tfirst\n"
        "chr1\t20\t30\tsecond\n"
        "chr1\t40\t50\tthird\n"
        "chr1\t60\t70\tunmapped\n",
        encoding="utf-8",
    )
    return path


def _install_batch_cache(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> CachedUCSCChainResource:
    context, index = _chain_context(tmp_path)
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: context,
    )
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_chain_resource",
        lambda *args, **kwargs: context,
    )
    monkeypatch.setattr(cli, "load_cached_chain_index", lambda *args, **kwargs: index)

    def forbidden_discovery(*args: object, **kwargs: object) -> None:
        raise AssertionError("BED batch mode must not contact the provider")

    monkeypatch.setattr(cli, "discover_ucsc_resources", forbidden_discovery)
    return context


def test_bed_batch_cli_emits_indexed_chain_only_json_without_provider_access(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context = _install_batch_cache(monkeypatch, tmp_path)
    args = cli._build_parser().parse_args(
        [
            SOURCE_DB,
            TARGET_DB,
            "--bed",
            str(_bed_path(tmp_path)),
            "--cache-dir",
            str(tmp_path / "cache"),
            "--json",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["schema_version"] == 2
    assert payload["report_type"] == "liftassess.ucsc_batch_result"
    assert payload["semantics"]["evidence_scope"] == "indexed_chain_only"
    assert payload["evidence"]["assessment_scope"] == "CHAIN_ONLY"
    assert payload["evidence"]["comparative_net_reciprocal_best"] == "NOT_ASSESSED"
    assert payload["scope"]["authoritative_source_sequence_preflight"] == (
        "NOT_ASSESSED"
    )
    assert len(payload["records"]) == 4
    assert [len(record["candidates"]) for record in payload["records"]] == [1, 1, 1, 0]
    assert payload["records"][0]["source_interval"]["start"] == 0
    assert payload["records"][0]["source_interval"]["coordinate_system"] == (
        "0-based-half-open"
    )
    kinds = [relationship["kind"] for relationship in payload["relationships"]]
    assert kinds.count("EXACT_TARGET_COLLISION") == 1
    assert kinds.count("OVERLAPPING_TARGET_PROJECTIONS") == 2
    assert payload["resource"]["sha256"] == context.chain.sha256
    assert payload["resource"]["consumed_by_engine"] is True
    assert payload["provenance"]["alignment_source_id"] == (
        f"ucsc-pair:{SOURCE_DB}:{TARGET_DB}"
    )
    assert "whole-chain fallback are disabled" in stderr.getvalue()


def test_bed_batch_cli_summary_preserves_bed_coordinate_convention(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_batch_cache(monkeypatch, tmp_path)
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(_bed_path(tmp_path))]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    rendered = stdout.getvalue()
    assert "* BATCH CHAIN PROJECTIONS *" in rendered
    assert "row-1 [first]" in rendered
    assert "chr1:0-10 (0-based half-open)" in rendered
    assert "EXACT_TARGET_COLLISION=1" in rendered
    assert "OVERLAPPING_TARGET_PROJECTIONS=2" in rendered


def test_bed_batch_cli_rejects_details_before_batch_execution(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    def forbidden_execution(*args: object, **kwargs: object) -> None:
        raise AssertionError("batch execution must not start for unsupported --details")

    monkeypatch.setattr(cli, "run_indexed_chain_batch", forbidden_execution)
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(_bed_path(tmp_path)), "--details"]
    )

    with pytest.raises(ValueError, match="--details is not yet available with --bed"):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_requires_prepared_index_without_fallback(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context, _index = _chain_context(tmp_path)
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: context,
    )
    monkeypatch.setattr(cli, "load_cached_chain_index", lambda *args, **kwargs: None)
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(_bed_path(tmp_path))]
    )

    with pytest.raises(ValueError, match="prepare-liftassess-index"):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_rejects_single_locus_and_bed_together(tmp_path: Path) -> None:
    args = cli._build_parser().parse_args(
        [
            SOURCE_DB,
            TARGET_DB,
            "chr1:1-10",
            "--bed",
            str(_bed_path(tmp_path)),
        ]
    )

    with pytest.raises(ValueError, match="either a single locus or --bed"):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_rejects_refresh_and_point_context(tmp_path: Path) -> None:
    refresh_args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(_bed_path(tmp_path)), "--refresh"]
    )
    with pytest.raises(ValueError, match="--refresh is not available with --bed"):
        cli._run(
            refresh_args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )

    context_args = cli._build_parser().parse_args(
        [
            SOURCE_DB,
            TARGET_DB,
            "--bed",
            str(_bed_path(tmp_path)),
            "--context-bases",
            "101",
        ]
    )
    with pytest.raises(ValueError, match="not yet available with --bed"):
        cli._run(
            context_args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_accepts_stdin(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_batch_cache(monkeypatch, tmp_path)
    args = cli._build_parser().parse_args([SOURCE_DB, TARGET_DB, "--bed", "-"])
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO("chr1\t0\t10\tstdin-row\n"),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert "row-1 [stdin-row]" in stdout.getvalue()


def test_cli_requires_locus_or_bed() -> None:
    args = cli._build_parser().parse_args([SOURCE_DB, TARGET_DB])

    with pytest.raises(ValueError, match="source locus is required unless --bed"):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )
