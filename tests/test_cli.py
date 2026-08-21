from __future__ import annotations

import gzip
import json
from io import StringIO
from pathlib import Path

import pytest

from liftassess import (
    AssemblyIdentifier,
    CachedResource,
    CachedUCSCResourceBundle,
    ChainIndexCorruptionError,
    EvidenceAvailabilityTier,
    GenomicInterval,
    ProvenanceSource,
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
    assert "Source: chr1:101-120 (1-based inclusive)" in stdout.getvalue()
    assert "Evidence: LIFTOVER-ONLY" in stdout.getvalue()
    assert "This does not establish biological correctness." in stdout.getvalue()
    assert "UCSC terms to review" in stderr.getvalue()
    assert "Transfer plan: LIFTOVER_ONLY (1 resource(s))" in stderr.getvalue()
    assert (
        "Provider-advertised total identity resource size: 2.0 KiB" in stderr.getvalue()
    )
    assert "Verified cache hits may avoid resource-body transfer." in stderr.getvalue()


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
    assert "Evidence: LIFTOVER-ONLY" in captured.out
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
    assert "Evidence: COMPARATIVE" in stdout.getvalue()
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
    assert "Source: chr1:101-120 (1-based inclusive)" in stdout.getvalue()


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
    assert "Evidence: LIFTOVER-ONLY" in stdout.getvalue()
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
    assert "Chain projections: 0" in stdout.getvalue()
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
