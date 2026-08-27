from __future__ import annotations

import gzip
import hashlib
import json
from io import StringIO
from pathlib import Path

import pytest

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCChainResource,
    CachedUCSCResourceBundle,
    ChainIndex,
    ChainIndexCorruptionError,
    ComparativeRelationshipState,
    EvidenceAvailabilityTier,
    GenomicInterval,
    ProvenanceSource,
    QueryContextState,
    ResourceReadProgressCallback,
    ReverseCheckState,
    UCSCAssessmentReport,
    UCSCBundleResourceRole,
    UCSCBundleTransferInspection,
    UCSCBundleTransferInspectionItem,
    UCSCRemoteResourceMetadata,
    UCSCResourceAcquisitionError,
    UCSCResourceBundle,
    UCSCResourceDiscoveryError,
    assess_ucsc_cached_bundle,
    build_cached_chain_index,
    cli,
    plan_ucsc_bundle_acquisition,
    sha256_identifier_for_file,
    ucsc_resource_terms,
)
from liftassess.resource_cache import _write_url_index

_SOURCE_DB = "canFam3"
_TARGET_DB = "canFam4"
_CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/"
    "canFam3ToCanFam4.over.chain.gz"
)


def _discovered_bundle() -> UCSCResourceBundle:
    return UCSCResourceBundle(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain_url=_CHAIN_URL,
    )


def _inspection(
    *,
    content_length_bytes: int | None = 2048,
    content_encoding: str | None = None,
) -> UCSCBundleTransferInspection:
    return UCSCBundleTransferInspection(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        items=(
            UCSCBundleTransferInspectionItem(
                role=UCSCBundleResourceRole.CHAIN,
                metadata=UCSCRemoteResourceMetadata(
                    url=_CHAIN_URL,
                    terms=ucsc_resource_terms(_CHAIN_URL),
                    content_length_bytes=content_length_bytes,
                    accept_ranges="bytes",
                    last_modified="Sun, 16 Aug 2026 00:00:00 GMT",
                    etag='"fixture"',
                    content_encoding=content_encoding,
                ),
            ),
        ),
    )


def _cached_bundle(tmp_path: Path) -> CachedUCSCResourceBundle:
    path = tmp_path / "chain.gz"
    chain = "chain 100 chr1 1000 + 100 120 chrA 2000 + 500 520 1\n20\n\n"
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(chain)
    resource = CachedResource(
        path=path,
        source_url=_CHAIN_URL,
        retrieved_at="2026-08-16T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(_CHAIN_URL),
        cache_hit=False,
    )
    return CachedUCSCResourceBundle(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=resource,
    )


def _cached_point_context_bundle(tmp_path: Path) -> CachedUCSCResourceBundle:
    path = tmp_path / "point-context-chain.gz"
    chain = "chain 100 chr1 2000 + 0 1500 chrA 3000 + 500 2000 1\n1500\n\n"
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(chain)
    resource = CachedResource(
        path=path,
        source_url=_CHAIN_URL,
        retrieved_at="2026-08-16T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(_CHAIN_URL),
        cache_hit=False,
    )
    return CachedUCSCResourceBundle(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=resource,
    )


def _cached_reverse_bundle(tmp_path: Path) -> CachedUCSCResourceBundle:
    path = tmp_path / "reverse-chain.gz"
    chain = "chain 100 chrA 2000 + 500 520 chr1 1000 + 100 120 2\n20\n\n"
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(chain)
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/"
        "canFam4ToCanFam3.over.chain.gz"
    )
    resource = CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-16T00:00:00Z",
        sha256=sha256_identifier_for_file(path).value,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=False,
    )
    return CachedUCSCResourceBundle(
        source_db=_TARGET_DB,
        target_db=_SOURCE_DB,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=resource,
    )


def _cached_reverse_chain(tmp_path: Path) -> CachedUCSCChainResource:
    bundle = _cached_reverse_bundle(tmp_path)
    return CachedUCSCChainResource(
        source_db=bundle.source_db,
        target_db=bundle.target_db,
        evidence_tier=bundle.evidence_tier,
        chain=bundle.chain,
    )


def _publish_cached_chain_for_resolution(
    cache_root: Path,
    *,
    source_url: str,
    data: bytes,
) -> None:
    digest = hashlib.sha256(data).hexdigest()
    artifact = cache_root / "artifacts" / "sha256" / digest[:2] / digest
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_bytes(data)
    index_key = hashlib.sha256(source_url.encode("utf-8")).hexdigest()
    _write_url_index(
        cache_root / "by-url" / f"{index_key}.json",
        source_url=source_url,
        retrieved_at="2026-08-23T00:00:00Z",
        sha256=digest,
        size_bytes=len(data),
        provider_checksum=None,
        terms=ucsc_resource_terms(source_url),
    )


def _install_successful_resource_flow(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: _discovered_bundle(),
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _inspection(),
    )
    cached = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        cli,
        "acquire_ucsc_resource_bundle",
        lambda plan, cache_root, **kwargs: cached,
    )


def test_cli_runs_end_to_end_with_interactive_acknowledgements(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_successful_resource_flow(monkeypatch, tmp_path)
    stdin = StringIO("yes\nyes\n")
    stdout = StringIO()
    stderr = StringIO()
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--cache-dir", str(tmp_path / "cache")]
    )

    exit_code = cli._run(args, stdin=stdin, stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert "Source:\n    chr1:101-120 (1-based inclusive)" in stdout.getvalue()
    assert "Evidence:\n    LIFTOVER-ONLY" in stdout.getvalue()
    assert "This does not establish biological correctness." in stdout.getvalue()
    assert "UCSC terms to review" in stderr.getvalue()
    assert "Transfer plan: LIFTOVER_ONLY (1 resource(s))" in stderr.getvalue()
    assert (
        "Provider-advertised total identity resource size: 2.0 KiB" in stderr.getvalue()
    )
    assert "Verified cache hits may avoid resource-body transfer." in stderr.getvalue()


def test_run_uses_cached_indexed_reverse_chain_without_provider_access(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_bundle(tmp_path)
    reverse = _cached_reverse_chain(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, reverse.chain)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target: forward,
    )

    seen_tiers: list[EvidenceAvailabilityTier] = []

    def resolve_reverse(*args: object, **kwargs: object) -> CachedUCSCChainResource:
        del args
        tier = kwargs["evidence_tier"]
        assert isinstance(tier, EvidenceAvailabilityTier)
        seen_tiers.append(tier)
        return reverse

    monkeypatch.setattr(
        cli, "resolve_cached_ucsc_chain_resource_metadata", resolve_reverse
    )
    monkeypatch.setattr(
        cli, "load_cached_ucsc_chain_resource", lambda *args, **kwargs: reverse
    )
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: (_ for _ in ()).throw(
            AssertionError("cached forward/reverse assessment must not contact UCSC")
        ),
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(cache_root),
            "--offline",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert seen_tiers == [EvidenceAvailabilityTier.LIFTOVER_ONLY]
    assert (
        "Reverse mapping:\n    exactly reconstructs the original aligned source geometry"
        in stdout.getvalue()
    )
    assert (
        "Assessing actual reverse mapping from cached indexed canFam4→canFam3 chain"
        in stderr.getvalue()
    )


def test_refresh_leaves_reverse_not_run_without_reverse_cache_access(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_bundle(tmp_path)
    report = assess_ucsc_cached_bundle(
        GenomicInterval(
            AssemblyIdentifier(name=_SOURCE_DB, provider="UCSC"),
            "chr1",
            100,
            120,
        ),
        forward,
        target_assembly=AssemblyIdentifier(name=_TARGET_DB, provider="UCSC"),
        alignment_provenance=ProvenanceSource("forward", "forward lineage"),
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("--refresh must not inspect automatic reverse resources")
        ),
    )
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--refresh"]
    )
    stderr = StringIO()

    enriched = cli._attach_cached_reverse_mapping_context(
        report,
        args=args,
        cache_root=tmp_path / "cache",
        stderr=stderr,
    )

    assert enriched.result_profile.scope.reverse_result is ReverseCheckState.NOT_RUN
    assert "not run during --refresh" in stderr.getvalue()


def test_run_reports_reverse_unavailable_without_matching_cached_chain(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target: forward,
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: None,
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--offline",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert (
        "Reverse mapping:\n    unavailable from the current prepared reverse resources"
        in stdout.getvalue()
    )
    assert "UCSC was not contacted" in stderr.getvalue()


def test_reverse_cache_resolution_requires_matching_publication_class(
    tmp_path: Path,
) -> None:
    forward = _comparative_cached_bundle(tmp_path)
    report = assess_ucsc_cached_bundle(
        GenomicInterval(
            AssemblyIdentifier(name=_SOURCE_DB, provider="UCSC"),
            "chr1",
            100,
            120,
        ),
        forward,
        target_assembly=AssemblyIdentifier(name=_TARGET_DB, provider="UCSC"),
        alignment_provenance=ProvenanceSource("forward", "forward lineage"),
    )
    assert report.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE

    cache_root = tmp_path / "reverse-cache"
    reverse_liftover_url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/"
        "canFam4ToCanFam3.over.chain.gz"
    )
    _publish_cached_chain_for_resolution(
        cache_root,
        source_url=reverse_liftover_url,
        data=b"filtered reverse chain",
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(cache_root),
            "--offline",
        ]
    )
    stderr = StringIO()

    enriched = cli._attach_cached_reverse_mapping_context(
        report,
        args=args,
        cache_root=cache_root,
        stderr=stderr,
    )

    assert enriched.result_profile.scope.reverse_result is ReverseCheckState.UNAVAILABLE
    assert "matching COMPARATIVE publication class" in stderr.getvalue()


def test_reverse_index_load_corruption_marks_not_run_before_chain_load(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_bundle(tmp_path)
    report = assess_ucsc_cached_bundle(
        GenomicInterval(
            AssemblyIdentifier(name=_SOURCE_DB, provider="UCSC"),
            "chr1",
            100,
            120,
        ),
        forward,
        target_assembly=AssemblyIdentifier(name=_TARGET_DB, provider="UCSC"),
        alignment_provenance=ProvenanceSource("forward", "forward lineage"),
    )
    reverse = _cached_reverse_chain(tmp_path)
    monkeypatch.setattr(
        cli, "resolve_cached_ucsc_chain_resource_metadata", lambda *a, **k: reverse
    )
    monkeypatch.setattr(
        cli,
        "load_cached_chain_index",
        lambda *a, **k: (_ for _ in ()).throw(
            ChainIndexCorruptionError("fixture reverse index load corruption")
        ),
    )
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_chain_resource",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("reverse chain must not load after index-load corruption")
        ),
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--offline",
        ]
    )
    stderr = StringIO()

    enriched = cli._attach_cached_reverse_mapping_context(
        report,
        args=args,
        cache_root=tmp_path / "cache",
        stderr=stderr,
    )

    assert enriched.result_profile.scope.reverse_result is ReverseCheckState.NOT_RUN
    assert "cached reverse chain index is unusable" in stderr.getvalue()


def test_run_marks_reverse_not_run_after_index_lookup_corruption(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_bundle(tmp_path)
    reverse = _cached_reverse_chain(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, reverse.chain)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli, "resolve_cached_ucsc_chain_resource_metadata", lambda *a, **k: reverse
    )
    monkeypatch.setattr(cli, "load_cached_ucsc_chain_resource", lambda *a, **k: reverse)
    monkeypatch.setattr(
        cli,
        "build_reverse_mapping_results_from_cached_chain",
        lambda *a, **k: (_ for _ in ()).throw(
            ChainIndexCorruptionError("fixture reverse query corruption")
        ),
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(cache_root),
            "--offline",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert "failed during lookup" in stderr.getvalue()
    assert "no full reverse-chain scan was started" in stderr.getvalue()
    assert "Reverse mapping:\n    not run" in stdout.getvalue()


def test_run_marks_reverse_not_run_when_matching_chain_is_not_indexed(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_bundle(tmp_path)
    reverse = _cached_reverse_chain(tmp_path)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli, "resolve_cached_ucsc_chain_resource_metadata", lambda *a, **k: reverse
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--offline",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert "no prepared index is available" in stderr.getvalue()
    assert "no full reverse-chain scan was started" in stderr.getvalue()
    assert "Reverse mapping:\n    not run" in stdout.getvalue()


def test_cli_explicit_acknowledgement_flags_skip_prompts(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _install_successful_resource_flow(monkeypatch, tmp_path)
    stdout = StringIO()
    stderr = StringIO()
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
            "--quiet",
        ]
    )

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert "Discovering UCSC resources" not in stderr.getvalue()
    assert "[y/N]" not in stderr.getvalue()
    assert "UCSC terms to review" in stderr.getvalue()
    assert "Transfer plan:" in stderr.getvalue()


def test_cli_declining_terms_stops_before_provider_inspection(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: _discovered_bundle(),
    )

    def unexpected_inspection(*args: object, **kwargs: object) -> None:
        raise AssertionError(
            "provider inspection must not run after terms are declined"
        )

    monkeypatch.setattr(cli, "inspect_ucsc_bundle_transfer_plan", unexpected_inspection)
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--cache-dir", str(tmp_path / "cache")]
    )
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO("no\n"), stdout=StringIO(), stderr=stderr)

    assert exit_code == 1
    assert (
        "Cancelled before UCSC resource inspection or acquisition." in stderr.getvalue()
    )


def test_cli_declining_transfer_stops_before_acquisition(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: _discovered_bundle(),
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _inspection(),
    )

    def unexpected_acquisition(*args: object, **kwargs: object) -> None:
        raise AssertionError("acquisition must not run after transfer is declined")

    monkeypatch.setattr(cli, "acquire_ucsc_resource_bundle", unexpected_acquisition)
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--cache-dir", str(tmp_path / "cache")]
    )
    stderr = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO("yes\nno\n"),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 1
    assert "Cancelled before UCSC resource acquisition." in stderr.getvalue()


def test_cli_reports_missing_resource_pair_without_planning(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(cli, "discover_ucsc_resources", lambda source, target: None)
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--cache-dir", str(tmp_path / "cache")]
    )
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=StringIO(), stderr=stderr)

    assert exit_code == 1
    assert (
        f"no supported UCSC resources found for {_SOURCE_DB}→{_TARGET_DB}"
        in stderr.getvalue()
    )


def test_default_cache_root_uses_macos_user_cache() -> None:
    assert cli._default_user_cache_root(
        platform_name="darwin",
        os_name="posix",
        environ={},
        home=Path("/Users/researcher"),
    ) == Path("/Users/researcher/Library/Caches/liftassess")


def test_default_cache_root_uses_xdg_cache_on_linux() -> None:
    assert cli._default_user_cache_root(
        platform_name="linux",
        os_name="posix",
        environ={"XDG_CACHE_HOME": "/scratch/cache"},
        home=Path("/home/researcher"),
    ) == Path("/scratch/cache/liftassess")


def test_default_cache_root_uses_localappdata_on_windows() -> None:
    assert (
        cli._default_user_cache_root(
            platform_name="win32",
            os_name="nt",
            environ={"LOCALAPPDATA": r"C:\Users\researcher\AppData\Local"},
            home=Path(r"C:\Users\researcher"),
        )
        == Path(r"C:\Users\researcher\AppData\Local") / "liftassess" / "Cache"
    )


def test_format_bytes_is_human_readable() -> None:
    assert cli._format_bytes(512) == "512 B"
    assert cli._format_bytes(2048) == "2.0 KiB"
    assert cli._format_bytes(3 * 1024**3) == "3.0 GiB"


def _comparative_discovered_bundle() -> UCSCResourceBundle:
    forward = "https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/vsCanFam4/"
    reciprocal = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/vsCanFam3/reciprocalBest/"
    )
    return UCSCResourceBundle(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain_url=f"{forward}canFam3.canFam4.all.chain.gz",
        net_url=f"{forward}canFam3.canFam4.net.gz",
        syntenic_net_url=f"{forward}canFam3.canFam4.syn.net.gz",
        reciprocal_best_chain_url=f"{reciprocal}canFam3.canFam4.rbest.chain.gz",
        reciprocal_best_net_url=f"{reciprocal}canFam3.canFam4.rbest.net.gz",
    )


def _comparative_inspection() -> UCSCBundleTransferInspection:
    plan = plan_ucsc_bundle_acquisition(_comparative_discovered_bundle())
    return UCSCBundleTransferInspection(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        items=tuple(
            UCSCBundleTransferInspectionItem(
                role=item.role,
                metadata=UCSCRemoteResourceMetadata(
                    url=item.url,
                    terms=item.terms,
                    content_length_bytes=1024,
                    accept_ranges="bytes",
                    last_modified="Sun, 16 Aug 2026 00:00:00 GMT",
                    etag='"fixture"',
                    content_encoding=None,
                ),
            )
            for item in plan.items
        ),
    )


def _comparative_cached_bundle(tmp_path: Path) -> CachedUCSCResourceBundle:
    bundle = _comparative_discovered_bundle()
    assert bundle.net_url is not None
    assert bundle.syntenic_net_url is not None
    assert bundle.reciprocal_best_chain_url is not None
    assert bundle.reciprocal_best_net_url is not None

    def cached_resource(name: str, url: str, content: bytes) -> CachedResource:
        path = tmp_path / name
        path.write_bytes(content)
        return CachedResource(
            path=path,
            source_url=url,
            retrieved_at="2026-08-16T00:00:00Z",
            sha256=sha256_identifier_for_file(path).value,
            size_bytes=path.stat().st_size,
            provider_checksum=None,
            terms=ucsc_resource_terms(url),
            cache_hit=False,
        )

    def gzip_bytes(text: str) -> bytes:
        return gzip.compress(text.encode("utf-8"))

    chain = "chain 100 chr1 1000 + 100 120 chrA 2000 + 500 520 1\n20\n\n"
    rbest_chain = "chain 100 chr1 1000 + 100 120 chrA 2000 + 500 520 101\n20\n\n"
    net = (
        "net chr1 1000\n"
        " fill 100 20 chrA + 500 20 id 1 score 100 ali 20 qDup 0 type syn\n"
    )

    return CachedUCSCResourceBundle(
        source_db=_SOURCE_DB,
        target_db=_TARGET_DB,
        evidence_tier=EvidenceAvailabilityTier.COMPARATIVE,
        chain=cached_resource("chain.gz", bundle.chain_url, gzip_bytes(chain)),
        net=cached_resource("net.gz", bundle.net_url, gzip_bytes(net)),
        syntenic_net=cached_resource(
            "syn.net.gz", bundle.syntenic_net_url, b"retrieval context only"
        ),
        reciprocal_best_chain=cached_resource(
            "rbest.chain.gz",
            bundle.reciprocal_best_chain_url,
            gzip_bytes(rbest_chain),
        ),
        reciprocal_best_net=cached_resource(
            "rbest.net.gz",
            bundle.reciprocal_best_net_url,
            b"retrieval context only",
        ),
    )


def test_main_runs_success_path_through_console_boundary(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _install_successful_resource_flow(monkeypatch, tmp_path)

    exit_code = cli.main(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "Evidence:\n    LIFTOVER-ONLY" in captured.out
    assert "This does not establish biological correctness." in captured.out
    assert "UCSC terms to review" in captured.err


def test_main_translates_invalid_locus_to_error_exit(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = cli.main([_SOURCE_DB, _TARGET_DB, "chr1:120-101"])

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "error: locus end must be greater than or equal to start" in captured.err


def test_main_preserves_discovery_failure_as_error(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def fail_discovery(source: str, target: str) -> UCSCResourceBundle | None:
        raise UCSCResourceDiscoveryError("fixture discovery failure")

    monkeypatch.setattr(cli, "discover_ucsc_resources", fail_discovery)

    exit_code = cli.main(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "empty-cache"),
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "error: fixture discovery failure" in captured.err
    assert "no supported UCSC resources found" not in captured.err


def test_main_surfaces_acquisition_failure_without_report(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: _discovered_bundle(),
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _inspection(),
    )

    def fail_acquisition(*args: object, **kwargs: object) -> CachedUCSCResourceBundle:
        raise UCSCResourceAcquisitionError("fixture partial bundle failure")

    monkeypatch.setattr(cli, "acquire_ucsc_resource_bundle", fail_acquisition)

    exit_code = cli.main(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "error: fixture partial bundle failure" in captured.err
    assert captured.out == ""


def test_comparative_cli_run_shares_pair_provenance_across_consumed_resources(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    discovered = _comparative_discovered_bundle()
    cached = _comparative_cached_bundle(tmp_path)
    reports: list[UCSCAssessmentReport] = []

    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: discovered,
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _comparative_inspection(),
    )
    monkeypatch.setattr(
        cli,
        "acquire_ucsc_resource_bundle",
        lambda plan, cache_root, **kwargs: cached,
    )

    def capture_assessment(
        source_interval: GenomicInterval,
        bundle: CachedUCSCResourceBundle,
        *,
        target_assembly: AssemblyIdentifier,
        alignment_provenance: ProvenanceSource,
        progress_callback: object = None,
        chain_index: object = None,
    ) -> UCSCAssessmentReport:
        assert progress_callback is None
        assert chain_index is None
        report = assess_ucsc_cached_bundle(
            source_interval,
            bundle,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
        )
        reports.append(report)
        return report

    monkeypatch.setattr(cli, "assess_ucsc_cached_bundle", capture_assessment)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
        ]
    )
    stdout = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=StringIO())

    assert exit_code == 0
    assert "Evidence:\n    COMPARATIVE" in stdout.getvalue()
    assert len(reports) == 1
    report = reports[0]
    assert report.alignment_provenance.source_id == "ucsc-pair:canFam3:canFam4"
    consumed = tuple(
        resource for resource in report.resources if resource.consumed_by_engine
    )
    assert {resource.role for resource in consumed} == {
        UCSCBundleResourceRole.CHAIN,
        UCSCBundleResourceRole.NET,
        UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
    }
    assert all(
        resource.file_provenance is not None
        and resource.file_provenance.derived_from == (report.alignment_provenance,)
        for resource in consumed
    )


def test_comparative_cli_attaches_filtered_all_chain_from_prepared_index(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    comparative_dir = tmp_path / "comparative"
    comparative_dir.mkdir()
    filtered_dir = tmp_path / "filtered"
    filtered_dir.mkdir()
    report = assess_ucsc_cached_bundle(
        GenomicInterval(
            AssemblyIdentifier(name=_SOURCE_DB, provider="UCSC"),
            "chr1",
            100,
            120,
        ),
        _comparative_cached_bundle(comparative_dir),
        target_assembly=AssemblyIdentifier(name=_TARGET_DB, provider="UCSC"),
        alignment_provenance=cli._ucsc_pair_lineage_provenance(_SOURCE_DB, _TARGET_DB),
    )
    filtered_bundle = _cached_bundle(filtered_dir)
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    filtered_index = build_cached_chain_index(
        tmp_path / "index-cache",
        filtered_bundle.chain,
    ).index

    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: filtered_chain,
    )
    monkeypatch.setattr(
        cli,
        "load_cached_chain_index",
        lambda *args, **kwargs: filtered_index,
    )
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_chain_resource",
        lambda *args, **kwargs: filtered_chain,
    )
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--offline"]
    )
    stderr = StringIO()

    enriched = cli._attach_cached_filtered_all_chain_comparison(
        report,
        args=args,
        cache_root=tmp_path / "cache",
        stderr=stderr,
    )

    assert enriched.filtered_all_chain_comparison is not None
    assert (
        enriched.result_profile.scope.comparative_relationship
        is ComparativeRelationshipState.NO_COMPETING_FULL_PLACEMENTS
    )
    assert "Comparing ordinary filtered liftOver and all-chain placements" in (
        stderr.getvalue()
    )


def test_comparative_cli_run_renders_paired_filtered_all_chain_result(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    discovered = _comparative_discovered_bundle()
    comparative_dir = tmp_path / "comparative-cli"
    comparative_dir.mkdir()
    cached = _comparative_cached_bundle(comparative_dir)
    filtered_dir = tmp_path / "filtered-cli"
    filtered_dir.mkdir()
    filtered_bundle = _cached_bundle(filtered_dir)
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    filtered_index = build_cached_chain_index(
        tmp_path / "filtered-cli-index-cache",
        filtered_bundle.chain,
    ).index

    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: discovered,
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _comparative_inspection(),
    )
    monkeypatch.setattr(
        cli,
        "acquire_ucsc_resource_bundle",
        lambda plan, cache_root, **kwargs: cached,
    )

    def resolve_chain(
        cache_root: object,
        source_db: str,
        target_db: str,
        *,
        evidence_tier: EvidenceAvailabilityTier,
    ) -> CachedUCSCChainResource | None:
        del cache_root
        if (
            source_db == _SOURCE_DB
            and target_db == _TARGET_DB
            and evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
        ):
            return filtered_chain
        return None

    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        resolve_chain,
    )
    monkeypatch.setattr(
        cli,
        "load_cached_chain_index",
        lambda cache_root, resource: (
            filtered_index if resource.sha256 == filtered_bundle.chain.sha256 else None
        ),
    )
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_chain_resource",
        lambda *args, **kwargs: filtered_chain,
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == 0
    assert "Filtered/all-chain comparison:\n    inventories agree" in stdout.getvalue()
    assert (
        "\n    Comparing ordinary filtered liftOver and all-chain placements"
        in stderr.getvalue()
    )


def test_comparative_cli_skips_filtered_comparison_without_prepared_index(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    comparative_dir = tmp_path / "comparative-no-filtered-index"
    comparative_dir.mkdir()
    filtered_dir = tmp_path / "filtered-no-index"
    filtered_dir.mkdir()
    report = assess_ucsc_cached_bundle(
        GenomicInterval(
            AssemblyIdentifier(name=_SOURCE_DB, provider="UCSC"),
            "chr1",
            100,
            120,
        ),
        _comparative_cached_bundle(comparative_dir),
        target_assembly=AssemblyIdentifier(name=_TARGET_DB, provider="UCSC"),
        alignment_provenance=cli._ucsc_pair_lineage_provenance(_SOURCE_DB, _TARGET_DB),
    )
    filtered_bundle = _cached_bundle(filtered_dir)
    filtered_chain = CachedUCSCChainResource(
        source_db=filtered_bundle.source_db,
        target_db=filtered_bundle.target_db,
        evidence_tier=filtered_bundle.evidence_tier,
        chain=filtered_bundle.chain,
    )
    executed = False

    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: filtered_chain,
    )
    monkeypatch.setattr(
        cli,
        "load_cached_chain_index",
        lambda *args, **kwargs: None,
    )

    def fail_if_loaded(*args: object, **kwargs: object) -> CachedUCSCChainResource:
        nonlocal executed
        del args, kwargs
        executed = True
        raise AssertionError("filtered chain must not be queried without an index")

    monkeypatch.setattr(cli, "load_cached_ucsc_chain_resource", fail_if_loaded)
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--offline"]
    )
    stderr = StringIO()

    unchanged = cli._attach_cached_filtered_all_chain_comparison(
        report,
        args=args,
        cache_root=tmp_path / "cache",
        stderr=stderr,
    )

    assert unchanged is report
    assert executed is False
    assert (
        report.result_profile.scope.comparative_relationship
        is ComparativeRelationshipState.NOT_ASSESSED
    )
    assert "no prepared index is available" in stderr.getvalue()
    assert "no full filtered-chain scan was started" in stderr.getvalue()


def test_comparative_cli_refresh_does_not_mix_unrefreshed_filtered_chain(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    comparative_dir = tmp_path / "comparative-refresh"
    comparative_dir.mkdir()
    report = assess_ucsc_cached_bundle(
        GenomicInterval(
            AssemblyIdentifier(name=_SOURCE_DB, provider="UCSC"),
            "chr1",
            100,
            120,
        ),
        _comparative_cached_bundle(comparative_dir),
        target_assembly=AssemblyIdentifier(name=_TARGET_DB, provider="UCSC"),
        alignment_provenance=cli._ucsc_pair_lineage_provenance(_SOURCE_DB, _TARGET_DB),
    )
    called = False

    def fail_if_called(*args: object, **kwargs: object) -> None:
        nonlocal called
        del args, kwargs
        called = True
        raise AssertionError("refresh must not consume an unrefreshed filtered chain")

    monkeypatch.setattr(
        cli, "resolve_cached_ucsc_chain_resource_metadata", fail_if_called
    )
    args = cli._build_parser().parse_args(
        [_SOURCE_DB, _TARGET_DB, "chr1:101-120", "--refresh"]
    )
    stderr = StringIO()

    unchanged = cli._attach_cached_filtered_all_chain_comparison(
        report,
        args=args,
        cache_root=tmp_path / "cache",
        stderr=stderr,
    )

    assert unchanged is report
    assert called is False
    assert "not run during --refresh" in stderr.getvalue()


def test_pair_lineage_is_stable_per_direction_and_distinct_across_pairs() -> None:
    first = cli._ucsc_pair_lineage_provenance("canFam3", "canFam4")
    repeated = cli._ucsc_pair_lineage_provenance("canFam3", "canFam4")
    reverse = cli._ucsc_pair_lineage_provenance("canFam4", "canFam3")

    assert first == repeated
    assert first.source_id == "ucsc-pair:canFam3:canFam4"
    assert reverse.source_id == "ucsc-pair:canFam4:canFam3"
    assert first != reverse


def test_refresh_propagates_to_acquisition_and_transfer_display(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: _discovered_bundle(),
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _inspection(),
    )
    captured_kwargs: dict[str, object] = {}
    cached = _cached_bundle(tmp_path)

    def capture_acquisition(
        plan: object,
        cache_root: object,
        **kwargs: object,
    ) -> CachedUCSCResourceBundle:
        captured_kwargs.update(kwargs)
        return cached

    monkeypatch.setattr(cli, "acquire_ucsc_resource_bundle", capture_acquisition)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--refresh",
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
        ]
    )
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=StringIO(), stderr=stderr)

    assert exit_code == 0
    assert captured_kwargs["refresh"] is True
    assert captured_kwargs["terms_acknowledged"] is True
    assert captured_kwargs["transfer_plan_acknowledged"] is True
    assert "Refresh cached resources: yes" in stderr.getvalue()
    assert (
        "Verified cache hits may avoid resource-body transfer." not in stderr.getvalue()
    )


def test_transfer_plan_reports_unavailable_content_length(tmp_path: Path) -> None:
    plan = plan_ucsc_bundle_acquisition(_discovered_bundle())
    stderr = StringIO()

    cli._print_transfer_plan(
        plan,
        _inspection(content_length_bytes=None),
        cache_root=tmp_path,
        refresh=False,
        stderr=stderr,
    )

    text = stderr.getvalue()
    assert "HTTP Content-Length unavailable" in text
    assert "Provider-advertised total identity resource size: unknown" in text


def test_transfer_plan_labels_nonidentity_content_encoding(tmp_path: Path) -> None:
    plan = plan_ucsc_bundle_acquisition(_discovered_bundle())
    stderr = StringIO()

    cli._print_transfer_plan(
        plan,
        _inspection(content_length_bytes=2048, content_encoding="gzip"),
        cache_root=tmp_path,
        refresh=False,
        stderr=stderr,
    )

    text = stderr.getvalue()
    assert "HTTP Content-Length 2.0 KiB; Content-Encoding gzip" in text
    assert "Provider-advertised total identity resource size: unknown" in text


def test_confirm_fails_closed_on_eof() -> None:
    stderr = StringIO()

    accepted = cli._confirm("Continue?", stdin=StringIO(""), stderr=stderr)

    assert accepted is False
    assert stderr.getvalue() == "Continue? [y/N] \n"


class _TTYStringIO(StringIO):
    def isatty(self) -> bool:
        return True


def test_transfer_progress_display_is_resume_aware_and_marks_cache_hits() -> None:
    plan = plan_ucsc_bundle_acquisition(_discovered_bundle())
    stderr = _TTYStringIO()
    display = cli._TransferProgressDisplay(plan, _inspection(), stderr=stderr)

    display.start()
    display.update(UCSCBundleResourceRole.CHAIN, 1024, 2048, False)
    resumed = stderr.getvalue()
    display.update(UCSCBundleResourceRole.CHAIN, 2048, 2048, True)

    assert "50%" in resumed
    assert "1.00 KiB / 2.00 KiB" in resumed
    assert "cached (2.00 KiB)" in stderr.getvalue()


def test_transfer_progress_display_handles_unknown_size_without_fake_percentage() -> (
    None
):
    plan = plan_ucsc_bundle_acquisition(_discovered_bundle())
    stderr = _TTYStringIO()
    display = cli._TransferProgressDisplay(
        plan,
        _inspection(content_length_bytes=None),
        stderr=stderr,
    )

    display.start()
    display.update(UCSCBundleResourceRole.CHAIN, 512, None, False)

    text = stderr.getvalue()
    assert "512 B complete" in text
    assert "100%" not in text


def test_run_wires_measured_transfer_progress_on_tty(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: _discovered_bundle(),
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _inspection(),
    )
    cached = _cached_bundle(tmp_path)

    def acquire(
        plan: object,
        cache_root: object,
        **kwargs: object,
    ) -> CachedUCSCResourceBundle:
        del plan, cache_root
        progress_callback = kwargs["progress_callback"]
        assert callable(progress_callback)
        progress_callback(UCSCBundleResourceRole.CHAIN, 0, 2048, False)
        progress_callback(UCSCBundleResourceRole.CHAIN, 1024, 2048, False)
        progress_callback(UCSCBundleResourceRole.CHAIN, 2048, 2048, False)
        return cached

    monkeypatch.setattr(cli, "acquire_ucsc_resource_bundle", acquire)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--refresh",
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
        ]
    )
    stderr = _TTYStringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 0
    text = stderr.getvalue()
    assert "Acquiring/verifying UCSC resources..." in text
    assert "50%" in text
    assert "1.00 KiB / 2.00 KiB" in text


def test_transfer_progress_callback_is_suppressed_for_quiet_and_non_tty(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "discover_ucsc_resources",
        lambda source, target: _discovered_bundle(),
    )
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _inspection(),
    )
    cached = _cached_bundle(tmp_path)
    callbacks: list[object] = []

    def acquire(
        plan: object,
        cache_root: object,
        **kwargs: object,
    ) -> CachedUCSCResourceBundle:
        del plan, cache_root
        callbacks.append(kwargs["progress_callback"])
        return cached

    monkeypatch.setattr(cli, "acquire_ucsc_resource_bundle", acquire)
    base = [
        _SOURCE_DB,
        _TARGET_DB,
        "chr1:101-120",
        "--cache-dir",
        str(tmp_path / "cache"),
        "--refresh",
        "--acknowledge-ucsc-terms",
        "--accept-transfer-plan",
    ]

    quiet_args = cli._build_parser().parse_args([*base, "--quiet"])
    assert (
        cli._run(
            quiet_args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=_TTYStringIO(),
        )
        == 0
    )
    non_tty_args = cli._build_parser().parse_args(base)
    assert (
        cli._run(
            non_tty_args,
            stdin=StringIO(""),
            stdout=StringIO(),
            stderr=StringIO(),
        )
        == 0
    )

    assert callbacks == [None, None]


def test_run_uses_matching_cached_chain_index_when_present(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    built = build_cached_chain_index(cache_root, cached.chain)

    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target: cached,
    )
    seen_indexes: list[object] = []

    def capture_assessment(
        source_interval: GenomicInterval,
        bundle: CachedUCSCResourceBundle,
        *,
        target_assembly: AssemblyIdentifier,
        alignment_provenance: ProvenanceSource,
        progress_callback: object = None,
        chain_index: object = None,
    ) -> UCSCAssessmentReport:
        seen_indexes.append(chain_index)
        return assess_ucsc_cached_bundle(
            source_interval,
            bundle,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
            chain_index=built.index,
        )

    monkeypatch.setattr(cli, "assess_ucsc_cached_bundle", capture_assessment)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(cache_root),
            "--offline",
        ]
    )
    stderr = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 0
    assert seen_indexes == [built.index]
    assert "Using verified cached chain index" in stderr.getvalue()


def test_run_retries_full_traversal_after_mid_query_index_corruption(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    built = build_cached_chain_index(cache_root, cached.chain)

    def load_cached(
        cache_root: Path,
        source: str,
        target: str,
        **kwargs: object,
    ) -> CachedUCSCResourceBundle:
        del cache_root, source, target, kwargs
        return cached

    monkeypatch.setattr(cli, "load_cached_ucsc_resource_bundle", load_cached)

    seen_indexes: list[object] = []

    def assess_with_corrupt_index_once(
        source_interval: GenomicInterval,
        bundle: CachedUCSCResourceBundle,
        *,
        target_assembly: AssemblyIdentifier,
        alignment_provenance: ProvenanceSource,
        progress_callback: object = None,
        chain_index: object = None,
    ) -> UCSCAssessmentReport:
        del progress_callback
        seen_indexes.append(chain_index)
        if chain_index is not None:
            raise ChainIndexCorruptionError("fixture query corruption")
        return assess_ucsc_cached_bundle(
            source_interval,
            bundle,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
            chain_index=None,
        )

    monkeypatch.setattr(
        cli, "assess_ucsc_cached_bundle", assess_with_corrupt_index_once
    )

    progress_modes: list[bool] = []

    class RecordingProgressDisplay:
        def __init__(
            self,
            bundle: CachedUCSCResourceBundle,
            *,
            stderr: StringIO,
            indexed_chain: bool = False,
        ) -> None:
            del bundle, stderr
            progress_modes.append(indexed_chain)

        def start(self) -> None:
            pass

        def update(
            self,
            role: UCSCBundleResourceRole,
            bytes_read: int,
            total_bytes: int,
        ) -> None:
            del role, bytes_read, total_bytes

        def finish(self, *, candidates_exist: bool) -> None:
            del candidates_exist

    monkeypatch.setattr(cli, "_AssessmentProgressDisplay", RecordingProgressDisplay)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(cache_root),
            "--offline",
        ]
    )
    stdout = StringIO()
    stderr = _TTYStringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == 0
    assert seen_indexes == [built.index, None]
    assert progress_modes == [True, False]
    assert "retrying with full traversal" in stderr.getvalue()
    assert "Source:\n    chr1:101-120 (1-based inclusive)" in stdout.getvalue()


def test_run_uses_validated_index_identity_to_skip_redundant_chain_rehash(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, cached.chain)
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_resource_bundle_metadata",
        lambda cache_root, source, target: cached,
    )
    trusted_seen: list[frozenset[str]] = []

    def load_cached(
        cache_root: Path,
        source: str,
        target: str,
        *,
        trusted_artifact_sha256_identifiers: frozenset[str],
    ) -> CachedUCSCResourceBundle:
        del cache_root, source, target
        trusted_seen.append(trusted_artifact_sha256_identifiers)
        return cached

    monkeypatch.setattr(
        cli, "load_cached_ucsc_resource_bundle_for_indexed_assessment", load_cached
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(cache_root),
            "--offline",
        ]
    )

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=StringIO(),
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert trusted_seen == [frozenset({cached.chain.sha256})]


def test_assessment_progress_marks_indexed_chain_without_fake_byte_progress(
    tmp_path: Path,
) -> None:
    stderr = _TTYStringIO()
    display = cli._AssessmentProgressDisplay(
        _cached_bundle(tmp_path),
        stderr=stderr,
        indexed_chain=True,
    )
    display.start()

    assert "indexed" in stderr.getvalue()


def test_cache_verification_progress_display_waits_for_integrity_success() -> None:
    stderr = _TTYStringIO()
    display = cli._CacheVerificationProgressDisplay(stderr=stderr)
    total = 4096

    display.update(0, total, False)
    display.update(total // 2, total, False)
    display.update(total, total, False)
    before_success = stderr.getvalue()
    display.update(total, total, True)

    assert "Cache verification" in before_success
    assert "50%" in before_success
    assert "99%" in before_success
    assert "100%" not in before_success
    assert "100%" in stderr.getvalue()
    assert cli._format_progress_bytes(total) in stderr.getvalue()


def test_run_wires_measured_cache_verification_progress_on_tty(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)

    def load_cached(
        cache_root: Path,
        source: str,
        target: str,
        *,
        progress_callback: object,
    ) -> CachedUCSCResourceBundle:
        del cache_root, source, target
        assert callable(progress_callback)
        progress_callback(0, 4096, False)
        progress_callback(2048, 4096, False)
        progress_callback(4096, 4096, False)
        progress_callback(4096, 4096, True)
        return cached

    monkeypatch.setattr(cli, "load_cached_ucsc_resource_bundle", load_cached)

    def forbidden_discovery(source: str, target: str) -> UCSCResourceBundle | None:
        raise AssertionError("provider discovery must not run for a complete cache hit")

    monkeypatch.setattr(cli, "discover_ucsc_resources", forbidden_discovery)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
        ]
    )
    stderr = _TTYStringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 0
    progress = stderr.getvalue()
    assert "Checking/verifying local UCSC cache..." in progress
    assert "Cache verification" in progress
    assert "50%" in progress
    assert "99%" in progress
    assert "100%" in progress


def test_quiet_suppresses_cache_verification_progress_callback(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)

    def load_cached(
        cache_root: Path, source: str, target: str
    ) -> CachedUCSCResourceBundle:
        del cache_root, source, target
        return cached

    monkeypatch.setattr(cli, "load_cached_ucsc_resource_bundle", load_cached)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--quiet",
        ]
    )
    stderr = _TTYStringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 0
    assert "Cache verification" not in stderr.getvalue()


def test_cli_reuses_complete_verified_cache_without_provider_access(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target: cached,
    )

    def forbidden_discovery(source: str, target: str) -> UCSCResourceBundle | None:
        raise AssertionError("provider discovery must not run for a complete cache hit")

    monkeypatch.setattr(cli, "discover_ucsc_resources", forbidden_discovery)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert "Evidence:\n    LIFTOVER-ONLY" in stdout.getvalue()
    assert (
        "Checking/verifying local UCSC cache...\n    Using verified cached"
        in stderr.getvalue()
    )
    assert "UCSC was not contacted" in stderr.getvalue()
    assert "UCSC terms to review" not in stderr.getvalue()


def test_details_flag_emits_full_dossier_from_cached_assessment(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target: cached,
    )

    def forbidden_discovery(source: str, target: str) -> UCSCResourceBundle | None:
        raise AssertionError("provider discovery must not run for a complete cache hit")

    monkeypatch.setattr(cli, "discover_ucsc_resources", forbidden_discovery)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--details",
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
    output = stdout.getvalue()
    assert "Detailed factual result dossier" in output
    assert "Chain 1" in output
    assert "Resources" in output
    assert "Provenance dependency graph" in output
    assert "Source locus: chr1:101-120 (1-based inclusive)" in output


def test_json_flag_emits_machine_readable_cached_assessment(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target: cached,
    )

    def forbidden_discovery(source: str, target: str) -> UCSCResourceBundle | None:
        raise AssertionError("provider discovery must not run for a complete cache hit")

    monkeypatch.setattr(cli, "discover_ucsc_resources", forbidden_discovery)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--json",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == 0
    assert "Checking/verifying local UCSC cache" in stderr.getvalue()
    payload = json.loads(stdout.getvalue())
    assert payload["schema_version"] == 2
    assert payload["source_interval"]["start"] == 100
    assert payload["source_interval"]["end"] == 120
    assert payload["result_profile"]["headline"] == "ONE_COMPLETE_CHAIN_PROJECTION"
    assert "aggregate_verdict" not in payload["semantics"]
    assert "assessment" not in payload
    assert payload["resources"][0]["role"] == "CHAIN"
    assert payload["caveat"] == "This does not establish biological correctness."


def test_details_and_json_output_modes_are_mutually_exclusive() -> None:
    parser = cli._build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                _SOURCE_DB,
                _TARGET_DB,
                "chr1:101-120",
                "--details",
                "--json",
            ]
        )


def test_offline_requires_complete_cached_bundle_without_provider_access(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target: None,
    )

    def forbidden_discovery(source: str, target: str) -> UCSCResourceBundle | None:
        raise AssertionError("--offline must guarantee zero provider access")

    monkeypatch.setattr(cli, "discover_ucsc_resources", forbidden_discovery)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "empty-cache"),
            "--offline",
        ]
    )
    stderr = StringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=StringIO(),
        stderr=stderr,
    )

    assert exit_code == 1
    assert (
        "--offline requires a complete verified cached UCSC bundle" in stderr.getvalue()
    )


def test_refresh_and_offline_are_mutually_exclusive() -> None:
    parser = cli._build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                _SOURCE_DB,
                _TARGET_DB,
                "chr1:101-120",
                "--refresh",
                "--offline",
            ]
        )


def test_run_integrates_zero_candidate_progress_without_consuming_comparative_evidence(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cached = _comparative_cached_bundle(tmp_path)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: cached,
    )

    def forbidden_discovery(source: str, target: str) -> UCSCResourceBundle | None:
        raise AssertionError("complete offline cache must prevent provider discovery")

    monkeypatch.setattr(cli, "discover_ucsc_resources", forbidden_discovery)
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:201-220",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--offline",
        ]
    )
    stdout = StringIO()
    stderr = _TTYStringIO()

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=stdout,
        stderr=stderr,
    )

    assert exit_code == 0
    assert "Chain projections:\n    0" in stdout.getvalue()
    progress = stderr.getvalue()
    assert "Chain" in progress
    assert "100%" in progress
    assert "Net" in progress
    assert "Reciprocal-best" in progress
    assert progress.count("not used") == 2


def test_assessment_progress_display_reports_measured_byte_percentage(
    tmp_path: Path,
) -> None:
    bundle = _cached_bundle(tmp_path)
    stderr = _TTYStringIO()
    display = cli._AssessmentProgressDisplay(bundle, stderr=stderr)
    total = bundle.chain.size_bytes

    display.start()
    display.update(UCSCBundleResourceRole.CHAIN, total // 2, total)
    display.update(UCSCBundleResourceRole.CHAIN, total, total)
    display.finish(candidates_exist=True)

    text = stderr.getvalue()
    assert "Chain" in text
    assert "50%" in text or "49%" in text
    assert "100%" in text
    assert "████████████████████" in text
    assert cli._format_progress_bytes(total) in text


def test_comparative_progress_display_starts_later_resources_as_pending(
    tmp_path: Path,
) -> None:
    bundle = _comparative_cached_bundle(tmp_path)
    stderr = _TTYStringIO()
    display = cli._AssessmentProgressDisplay(bundle, stderr=stderr)

    display.start()

    text = stderr.getvalue()
    assert "Chain" in text
    assert "Net" in text
    assert "Reciprocal-best" in text
    assert text.count("pending") == 3


def test_run_automatically_assesses_101bp_point_context_from_forward_index(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_point_context_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, forward.chain)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: None,
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-101",
            "--cache-dir",
            str(cache_root),
            "--offline",
            "--json",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    context = payload["query_context"]
    assert context["check_state"] == QueryContextState.RUN.value
    assert context["requested_window_bases"] == 101
    assert context["actual_window_bases"] == 101
    assert context["tested_source_interval"]["start"] == 50
    assert context["tested_source_interval"]["end"] == 151
    assert context["evidence_scope"] == "forward_chain_only"
    assert len(context["candidates"]) == 1
    assert "Assessing 101-bp point context" in stderr.getvalue()


def test_run_point_context_without_forward_index_is_not_run_without_extra_scan(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_point_context_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: None,
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-101",
            "--cache-dir",
            str(cache_root),
            "--offline",
            "--json",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    payload = json.loads(stdout.getvalue())
    assert payload["query_context"]["check_state"] == QueryContextState.NOT_RUN.value
    assert payload["query_context"]["not_run_reason"] == "INDEX_UNAVAILABLE"
    assert "no additional full chain scan was started" in stderr.getvalue()


def test_run_marks_point_context_not_run_after_index_lookup_corruption(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_point_context_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, forward.chain)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: None,
    )

    real_assess = assess_ucsc_cached_bundle
    assessment_calls = 0

    def count_assessment(
        source_interval: GenomicInterval,
        bundle: CachedUCSCResourceBundle,
        *,
        target_assembly: AssemblyIdentifier,
        alignment_provenance: ProvenanceSource,
        progress_callback: ResourceReadProgressCallback | None = None,
        chain_index: ChainIndex | None = None,
    ) -> UCSCAssessmentReport:
        nonlocal assessment_calls
        assessment_calls += 1
        return real_assess(
            source_interval,
            bundle,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
            progress_callback=progress_callback,
            chain_index=chain_index,
        )

    monkeypatch.setattr(cli, "assess_ucsc_cached_bundle", count_assessment)
    monkeypatch.setattr(
        cli,
        "attach_point_query_context",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            ChainIndexCorruptionError("fixture point-context query corruption")
        ),
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-101",
            "--cache-dir",
            str(cache_root),
            "--offline",
            "--json",
        ]
    )
    stdout = StringIO()
    stderr = StringIO()

    exit_code = cli._run(args, stdin=StringIO(""), stdout=stdout, stderr=stderr)

    assert exit_code == 0
    assert assessment_calls == 1
    context = json.loads(stdout.getvalue())["query_context"]
    assert context["check_state"] == QueryContextState.NOT_RUN.value
    assert context["not_run_reason"] == "INDEX_UNUSABLE"
    assert "no full chain fallback was started" in stderr.getvalue()


def test_run_accepts_explicit_larger_odd_point_context_window(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_point_context_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, forward.chain)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: None,
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:751-751",
            "--context-bases",
            "1001",
            "--cache-dir",
            str(cache_root),
            "--offline",
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
    context = json.loads(stdout.getvalue())["query_context"]
    assert context["requested_window_bases"] == 1001
    assert context["actual_window_bases"] == 1001
    assert context["tested_source_interval"]["start"] == 250
    assert context["tested_source_interval"]["end"] == 1251


def test_main_rejects_explicit_point_context_for_interval_query(
    capsys: pytest.CaptureFixture[str],
) -> None:
    exit_code = cli.main(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--context-bases",
            "1001",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "--context-bases currently requires a 1-bp point query" in captured.err


def test_run_point_context_details_state_chain_only_scope(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_point_context_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, forward.chain)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: None,
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-101",
            "--cache-dir",
            str(cache_root),
            "--offline",
            "--details",
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
    rendered = stdout.getvalue()
    assert "Point neighborhood context" in rendered
    assert (
        "Evidence scope: forward chain only; net/reciprocal-best not re-run" in rendered
    )
    assert "Tested source window: chr1:51-151 (1-based inclusive)" in rendered


def test_run_point_context_summary_reports_exact_tested_window(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    forward = _cached_point_context_bundle(tmp_path)
    cache_root = tmp_path / "cache"
    build_cached_chain_index(cache_root, forward.chain)
    monkeypatch.setattr(
        cli,
        "load_cached_ucsc_resource_bundle",
        lambda cache_root, source, target, **kwargs: forward,
    )
    monkeypatch.setattr(
        cli,
        "resolve_cached_ucsc_chain_resource_metadata",
        lambda *args, **kwargs: None,
    )
    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-101",
            "--cache-dir",
            str(cache_root),
            "--offline",
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
    rendered = stdout.getvalue()
    assert "Local context (forward chain only)" in rendered
    assert "chr1:51-151 (1-based inclusive); 101 bp tested" in rendered
    assert "point and local context map together" in rendered


def test_cli_can_explicitly_acquire_liftover_only_chain(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    requested: list[EvidenceAvailabilityTier | None] = []

    def discover(
        source: str,
        target: str,
        *,
        evidence_tier: EvidenceAvailabilityTier | None = None,
    ) -> UCSCResourceBundle:
        assert (source, target) == (_SOURCE_DB, _TARGET_DB)
        requested.append(evidence_tier)
        return _discovered_bundle()

    monkeypatch.setattr(cli, "discover_ucsc_resources", discover)
    monkeypatch.setattr(
        cli,
        "inspect_ucsc_bundle_transfer_plan",
        lambda plan, *, terms_acknowledged: _inspection(),
    )
    cached = _cached_bundle(tmp_path)
    monkeypatch.setattr(
        cli,
        "acquire_ucsc_resource_bundle",
        lambda plan, cache_root, **kwargs: cached,
    )

    args = cli._build_parser().parse_args(
        [
            _SOURCE_DB,
            _TARGET_DB,
            "chr1:101-120",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--evidence-tier",
            "LIFTOVER-ONLY",
            "--acknowledge-ucsc-terms",
            "--accept-transfer-plan",
        ]
    )

    exit_code = cli._run(
        args,
        stdin=StringIO(""),
        stdout=StringIO(),
        stderr=StringIO(),
    )

    assert exit_code == 0
    assert requested == [EvidenceAvailabilityTier.LIFTOVER_ONLY]
