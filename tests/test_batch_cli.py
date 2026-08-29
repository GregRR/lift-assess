from __future__ import annotations

import gzip
import json
from io import StringIO
from pathlib import Path

import pytest

from liftassess import cli
from liftassess.assembly_metadata_cache import CachedUCSCAssemblyMetadata
from liftassess.chain_index import ChainIndex, build_chain_index
from liftassess.models import EvidenceAvailabilityTier
from liftassess.resource_cache import (
    CachedResource,
    CachedUCSCChainResource,
    CachedUCSCResourceBundle,
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


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _cached_resource(path: Path, url: str) -> CachedResource:
    return CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-28T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=True,
    )


def _source_metadata(tmp_path: Path) -> CachedUCSCAssemblyMetadata:
    chrom_info_path = tmp_path / "chromInfo.txt.gz"
    chrom_alias_path = tmp_path / "chromAlias.txt.gz"
    _write_gzip(
        chrom_info_path,
        "chr1\t1000\t/gbdb/hg38/hg38.2bit\n"
        + "chrValidNoChain\t1000\t/gbdb/hg38/hg38.2bit\n",
    )
    _write_gzip(
        chrom_alias_path,
        "1\tchr1\tucscToEnsembl\n",
    )
    base = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/"
    return CachedUCSCAssemblyMetadata(
        db=SOURCE_DB,
        chrom_info=_cached_resource(chrom_info_path, f"{base}chromInfo.txt.gz"),
        chrom_alias=_cached_resource(chrom_alias_path, f"{base}chromAlias.txt.gz"),
    )


@pytest.fixture(autouse=True)
def _install_source_metadata_cache(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    metadata = _source_metadata(tmp_path)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_assembly_metadata",
        lambda *args, **kwargs: metadata,
    )


def _comparative_bundle_context(
    tmp_path: Path,
) -> tuple[CachedUCSCResourceBundle, ChainIndex]:
    chain_path = tmp_path / "hg38.hg19.all.chain.gz"
    net_path = tmp_path / "hg38.hg19.net.gz"
    syn_path = tmp_path / "hg38.hg19.syn.net.gz"
    rbest_chain_path = tmp_path / "hg38.hg19.rbest.chain.gz"
    rbest_net_path = tmp_path / "hg38.hg19.rbest.net.gz"
    _write_gzip(
        chain_path,
        "chain 100 chr1 1000 + 0 10 chrA 2000 + 100 110 1\n10\n\n",
    )
    _write_gzip(
        net_path,
        "net chr1 1000\n"
        " fill 0 10 chrA + 100 10 id 1 score 100 ali 10 qDup 0 type syn\n",
    )
    syn_path.write_bytes(b"not consumed")
    _write_gzip(
        rbest_chain_path,
        "chain 1 chr1 1000 + 0 10 chrA 2000 + 100 110 101\n10\n\n",
    )
    rbest_net_path.write_bytes(b"not consumed")
    forward = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsHg19/"
    reciprocal = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/vsHg38/reciprocalBest/"
    )
    bundle = CachedUCSCResourceBundle(
        source_db=SOURCE_DB,
        target_db=TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain=_cached_resource(chain_path, f"{forward}hg38.hg19.all.chain.gz"),
        net=_cached_resource(net_path, f"{forward}hg38.hg19.net.gz"),
        syntenic_net=_cached_resource(syn_path, f"{forward}hg38.hg19.syn.net.gz"),
        reciprocal_best_chain=_cached_resource(
            rbest_chain_path, f"{reciprocal}hg38.hg19.rbest.chain.gz"
        ),
        reciprocal_best_net=_cached_resource(
            rbest_net_path, f"{reciprocal}hg38.hg19.rbest.net.gz"
        ),
    )
    index = build_chain_index(
        chain_path,
        tmp_path / "comparative-index",
        source_chain_sha256_identifier=bundle.chain.sha256,
        source_chain_size_bytes=bundle.chain.size_bytes,
    ).index
    return bundle, index


def _point_context_chain_context(
    tmp_path: Path,
) -> tuple[CachedUCSCChainResource, ChainIndex]:
    path = tmp_path / "hg38ToHg19-point-context.over.chain"
    path.write_text(
        """\
chain 100 chr1 1000 + 100 201 chrA 2000 + 500 601 11
101

chain 90 chr1 1000 + 300 401 chrA 2000 + 500 601 12
101

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
        tmp_path / "point-context-index",
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=path.stat().st_size,
    ).index
    return context, index


def _point_bed_path(tmp_path: Path) -> Path:
    path = tmp_path / "points.bed"
    path.write_text(
        "chr1\t150\t151\tpoint-a\nchr1\t350\t351\tpoint-b\n",
        encoding="utf-8",
    )
    return path


def _install_specific_batch_cache(
    monkeypatch: pytest.MonkeyPatch,
    context: CachedUCSCChainResource,
    index: ChainIndex,
) -> None:
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


def _interval_table_path(tmp_path: Path) -> Path:
    path = tmp_path / "batch.tsv"
    path.write_text(
        "sequence\tstart\tend\tlabel\n"
        "chr1\t1\t10\tfirst\n"
        "chr1\t21\t30\tsecond\n"
        "chr1\t41\t50\tthird\n"
        "chr1\t61\t70\tunmapped\n",
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


def test_interval_table_batch_cli_normalizes_one_based_inclusive_input(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_batch_cache(monkeypatch, tmp_path)
    args = cli._build_parser().parse_args(
        [
            SOURCE_DB,
            TARGET_DB,
            "--interval-table",
            str(_interval_table_path(tmp_path)),
            "--json",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert [
        (
            record["source_interval"]["start"],
            record["source_interval"]["end"],
            record["label"],
        )
        for record in payload["records"]
    ] == [
        (0, 10, "first"),
        (20, 30, "second"),
        (40, 50, "third"),
        (60, 70, "unmapped"),
    ]
    assert payload["records"][0]["source_interval"]["coordinate_system"] == (
        "0-based-half-open"
    )
    assert "Input coordinates are 1-based inclusive" in stderr.getvalue()
    assert "Assessing interval-table batch" in stderr.getvalue()


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
    assert payload["scope"]["authoritative_source_sequence_preflight"] == "ASSESSED"
    assert [record["source_preflight"]["state"] for record in payload["records"]] == [
        "VALID",
        "VALID",
        "VALID",
        "VALID",
    ]
    assert len(payload["source_preflight_resources"]) == 1
    preflight_resource = payload["source_preflight_resources"][0]
    assert preflight_resource["source_url"].endswith("/chromInfo.txt.gz")
    assert preflight_resource["terms"]["restricted_liftover_chain"] is False
    provenance_ids = {item["source_id"] for item in payload["provenance"]["sources"]}
    assert f"file:{preflight_resource['sha256']}" in provenance_ids
    assert all(
        not item.get("source_url", "").endswith("/chromAlias.txt.gz")
        for item in payload["source_preflight_resources"]
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


def test_bed_batch_cli_uses_shared_comparative_evidence_when_bundle_is_complete(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    bundle, index = _comparative_bundle_context(tmp_path)
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_resource_bundle_metadata",
        lambda *args, **kwargs: bundle,
    )

    def forbidden_bundle_rehash(*args: object, **kwargs: object) -> None:
        raise AssertionError(
            "COMPARATIVE batch mode must verify net/rbest during shared parsing, "
            "not pre-hash them through the ordinary bundle loader"
        )

    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle_for_indexed_assessment",
        forbidden_bundle_rehash,
    )
    monkeypatch.setattr(cli, "load_cached_chain_index", lambda *args, **kwargs: index)
    bed = tmp_path / "comparative.bed"
    bed.write_text("chr1\t0\t10\tfirst\n", encoding="utf-8")
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(bed), "--json"]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["semantics"]["evidence_scope"] == (
        "indexed_chain_plus_shared_comparative"
    )
    assert payload["evidence"]["assessment_scope"] == "CHAIN_NET_RECIPROCAL_BEST"
    assert payload["evidence"]["comparative_net_reciprocal_best"] == (
        "ASSESSED_FOR_SUBMITTED_RECORDS"
    )
    assert payload["scope"]["filtered_all_chain_comparison"] == "NOT_ASSESSED"
    assert payload["scope"]["comparative_relationship_interpretation"] == (
        "NOT_ASSESSED"
    )
    assert [item["role"] for item in payload["comparative_resources"]] == [
        "NET",
        "RECIPROCAL_BEST_CHAIN",
    ]
    candidate_evidence = payload["records"][0]["candidates"][0]["evidence"]
    kinds = {item["kind"] for item in candidate_evidence}
    assert "NET_CLASSIFICATION" in kinds
    assert "RECIPROCAL_BEST_MEMBERSHIP" in kinds
    assert "shared net/reciprocal-best scans are enabled" in stderr.getvalue()


def test_bed_batch_cli_does_not_claim_unused_comparative_evidence(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    bundle, index = _comparative_bundle_context(tmp_path)
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_resource_bundle_metadata",
        lambda *args, **kwargs: bundle,
    )
    monkeypatch.setattr(cli, "load_cached_chain_index", lambda *args, **kwargs: index)
    bed = tmp_path / "comparative-no-candidate.bed"
    bed.write_text("chr1\t60\t70\tnone\n", encoding="utf-8")
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(bed), "--json"]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["semantics"]["evidence_scope"] == "indexed_chain_only"
    assert payload["evidence"]["assessment_scope"] == "CHAIN_ONLY"
    assert payload["evidence"]["comparative_net_reciprocal_best"] == (
        "NOT_USED_NO_SUBMITTED_CANDIDATES"
    )
    assert all(
        item["consumed_by_engine"] is False for item in payload["comparative_resources"]
    )
    provenance_ids = {item["source_id"] for item in payload["provenance"]["sources"]}
    assert bundle.net is not None
    assert bundle.reciprocal_best_chain is not None
    assert f"file:{bundle.net.sha256}" not in provenance_ids
    assert f"file:{bundle.reciprocal_best_chain.sha256}" not in provenance_ids


def test_bed_batch_cli_requires_cached_source_metadata_before_chain_lookup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_assembly_metadata",
        lambda *args, **kwargs: None,
    )

    def forbidden_chain_lookup(*args: object, **kwargs: object) -> None:
        raise AssertionError(
            "batch chain lookup must not start before source preflight"
        )

    monkeypatch.setattr(
        cli, "_resolve_preferred_cached_batch_chain", forbidden_chain_lookup
    )
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(_bed_path(tmp_path))]
    )

    with pytest.raises(ValueError, match="requires verified cached UCSC chromInfo"):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_rejects_verified_alias_before_chain_lookup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    bed = tmp_path / "alias.bed"
    bed.write_text("1\t0\t1\talias-row\n", encoding="utf-8")

    def forbidden_chain_lookup(*args: object, **kwargs: object) -> None:
        raise AssertionError(
            "batch chain lookup must not start for invalid source input"
        )

    monkeypatch.setattr(
        cli, "_resolve_preferred_cached_batch_chain", forbidden_chain_lookup
    )
    args = cli._build_parser().parse_args([SOURCE_DB, TARGET_DB, "--bed", str(bed)])

    with pytest.raises(
        ValueError,
        match=r"batch record row-1 .*chromAlias verifies 'chr1'.*not attempted",
    ):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_rejects_out_of_bounds_row_before_chain_lookup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    bed = tmp_path / "bounds.bed"
    bed.write_text("chr1\t999\t1001\tbounds-row\n", encoding="utf-8")

    def forbidden_chain_lookup(*args: object, **kwargs: object) -> None:
        raise AssertionError(
            "batch chain lookup must not start for invalid source input"
        )

    monkeypatch.setattr(
        cli, "_resolve_preferred_cached_batch_chain", forbidden_chain_lookup
    )
    args = cli._build_parser().parse_args([SOURCE_DB, TARGET_DB, "--bed", str(bed)])

    with pytest.raises(
        ValueError,
        match=r"batch record row-1 .*sequence length is 1000.*not attempted",
    ):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_valid_no_chain_sequence_uses_authoritative_point_bounds(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_batch_cache(monkeypatch, tmp_path)
    bed = tmp_path / "valid-no-chain.bed"
    bed.write_text(
        "chrValidNoChain\t100\t101\tvalid-no-chain\n",
        encoding="utf-8",
    )
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(bed), "--json"]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    record = payload["records"][0]
    assert record["source_preflight"]["state"] == "VALID"
    assert record["source_preflight"]["sequence_length"] == 1000
    assert record["candidates"] == []
    context = payload["point_context"]["records"][0]
    assert context["state"] == "RUN"
    assert context["tested_source_interval"]["start"] == 50
    assert context["tested_source_interval"]["end"] == 151
    assert context["candidates"] == []


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
    assert "Source preflight:" in rendered
    assert (
        "4/4 records valid against authoritative UCSC assembly-sequence metadata"
        in rendered
    )
    assert "Source metadata SHA-256:" in rendered
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

    with pytest.raises(
        ValueError, match="either a single locus or one batch input option"
    ):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_rejects_refresh_and_context_for_interval_only_batch(
    tmp_path: Path,
) -> None:
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
    with pytest.raises(ValueError, match="requires at least one 1-bp BED record"):
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


def test_interval_table_batch_cli_accepts_stdin(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_batch_cache(monkeypatch, tmp_path)
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--interval-table", "-"]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO("sequence\tstart\tend\tlabel\nchr1\t1\t10\tstdin-row\n"),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert "row-1 [stdin-row]" in stdout.getvalue()
    assert "chr1:0-10 (0-based half-open)" in stdout.getvalue()


def test_cli_requires_locus_or_batch_input() -> None:
    args = cli._build_parser().parse_args([SOURCE_DB, TARGET_DB])

    with pytest.raises(
        ValueError, match="source locus is required unless --bed or --interval-table"
    ):
        cli._run(
            args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )


def test_bed_batch_cli_automatically_reports_neighborhood_level_collision(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context, index = _point_context_chain_context(tmp_path)
    _install_specific_batch_cache(monkeypatch, context, index)
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(_point_bed_path(tmp_path))]
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
    assert "requested=101 bp; point records=2; run=2; not run=0" in rendered
    assert "NEIGHBORHOOD_LEVEL_TARGET_COLLISION=1" in rendered
    assert "source coverage=101/101" in rendered


def test_interval_table_points_receive_automatic_point_context(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context, index = _point_context_chain_context(tmp_path)
    _install_specific_batch_cache(monkeypatch, context, index)
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--interval-table", "-", "--json"]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(
            "sequence\tstart\tend\tlabel\n"
            "chr1\t151\t151\tpoint-a\n"
            "chr1\t351\t351\tpoint-b\n"
        ),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert [
        (record["source_interval"]["start"], record["source_interval"]["end"])
        for record in payload["records"]
    ] == [(150, 151), (350, 351)]
    assert [item["state"] for item in payload["point_context"]["records"]] == [
        "RUN",
        "RUN",
    ]
    assert [item["kind"] for item in payload["point_context"]["relationships"]] == [
        "NEIGHBORHOOD_LEVEL_TARGET_COLLISION"
    ]


def test_bed_batch_cli_json_exposes_point_context_as_separate_scale(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context, index = _point_context_chain_context(tmp_path)
    _install_specific_batch_cache(monkeypatch, context, index)
    args = cli._build_parser().parse_args(
        [
            SOURCE_DB,
            TARGET_DB,
            "--bed",
            str(_point_bed_path(tmp_path)),
            "--json",
        ]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["scope"]["point_context"] == "ASSESSED_FOR_ALL_POINT_RECORDS"
    assert payload["point_context"]["requested_window_bases"] == 101
    assert [item["state"] for item in payload["point_context"]["records"]] == [
        "RUN",
        "RUN",
    ]
    relationships = payload["point_context"]["relationships"]
    assert [item["kind"] for item in relationships] == [
        "NEIGHBORHOOD_LEVEL_TARGET_COLLISION"
    ]
    assert relationships[0]["overlap_intervals"][0]["start"] == 500
    assert relationships[0]["overlap_intervals"][0]["end"] == 601


def test_bed_batch_cli_json_keeps_point_context_record_shape_consistent(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context, index = _point_context_chain_context(tmp_path)
    _install_specific_batch_cache(monkeypatch, context, index)
    bed_path = tmp_path / "mixed-context.bed"
    bed_path.write_text(
        "chr1\t150\t151\tpoint\nchr1\t150\t160\tinterval\n",
        encoding="utf-8",
    )
    args = cli._build_parser().parse_args(
        [SOURCE_DB, TARGET_DB, "--bed", str(bed_path), "--json"]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    records = payload["point_context"]["records"]
    assert records[1] == {
        "record_id": "row-2",
        "state": "NOT_APPLICABLE",
        "requested_window_bases": 101,
        "not_run_reason": "SOURCE_INTERVAL_IS_NOT_ONE_BASE",
    }
    assert set(records[0]) >= {
        "record_id",
        "state",
        "requested_window_bases",
        "not_run_reason",
    }
    assert set(records[1]) == {
        "record_id",
        "state",
        "requested_window_bases",
        "not_run_reason",
    }


def test_bed_batch_cli_accepts_explicit_context_window_for_point_rows(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    context, index = _point_context_chain_context(tmp_path)
    _install_specific_batch_cache(monkeypatch, context, index)
    args = cli._build_parser().parse_args(
        [
            SOURCE_DB,
            TARGET_DB,
            "--bed",
            str(_point_bed_path(tmp_path)),
            "--context-bases",
            "51",
            "--json",
        ]
    )
    stdout = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=StringIO(),
    )

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["point_context"]["requested_window_bases"] == 51
    assert [
        item["tested_source_interval"]["end"] - item["tested_source_interval"]["start"]
        for item in payload["point_context"]["records"]
    ] == [51, 51]


def test_batch_cli_parser_rejects_bed_and_interval_table_together(
    tmp_path: Path,
) -> None:
    with pytest.raises(SystemExit):
        cli._build_parser().parse_args(
            [
                SOURCE_DB,
                TARGET_DB,
                "--bed",
                str(_bed_path(tmp_path)),
                "--interval-table",
                str(_interval_table_path(tmp_path)),
            ]
        )
