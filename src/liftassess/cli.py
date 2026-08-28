"""Command-line entry point for the common liftAssess UCSC workflow.

The CLI is deliberately a thin composition layer over the existing library boundaries:
UCSC-style input parsing, resource discovery/planning, explicit provider-terms and
transfer acknowledgements, cached acquisition, assessment orchestration, and summary
rendering.  It does not duplicate candidate interpretation or verdict logic.

The automatic UCSC path uses one conservative pair-level provenance node as the shared
upstream dependency for consumed UCSC resources.  This is a dependency-grouping choice,
not a claim that liftAssess reconstructed the provider's exact alignment pipeline from
file bytes.  Exact resource identity remains represented by the child file SHA-256
provenance nodes created at the cached-resource boundary.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TextIO

from .batch_execution import run_indexed_chain_batch
from .batch_input import parse_bed_batch
from .batch_reporting import (
    render_indexed_chain_batch_json,
    render_indexed_chain_batch_summary,
)
from .chain_index import ChainIndex, ChainIndexCorruptionError, load_cached_chain_index
from .cli_input import parse_ucsc_locus, ucsc_assembly_identifier
from .comparative_inventory import FilteredAllChainCorrespondenceError
from .models import AssemblyIdentifier, EvidenceAvailabilityTier, ProvenanceSource
from .orchestration import (
    UCSCAssessmentReport,
    assess_ucsc_cached_bundle,
    attach_filtered_all_chain_comparison,
    attach_point_query_context,
    attach_query_context_result,
    attach_reverse_mapping_results,
)
from .query_context import (
    DEFAULT_POINT_CONTEXT_BASES,
    QueryContextNotRunReason,
    point_context_not_run,
)
from .reporting import (
    render_assessment_details,
    render_assessment_json,
    render_assessment_summary,
)
from .resource_cache import (
    CachedUCSCChainResource,
    CachedUCSCResourceBundle,
    CacheVerificationProgressCallback,
    UCSCBundleAcquisitionPlan,
    UCSCBundleResourceRole,
    UCSCBundleTransferInspection,
    UCSCBundleTransferProgressCallback,
    UCSCResourceAcquisitionError,
    acquire_ucsc_resource_bundle,
    inspect_ucsc_bundle_transfer_plan,
    load_cached_ucsc_chain_resource,
    load_cached_ucsc_resource_bundle,
    load_cached_ucsc_resource_bundle_for_indexed_assessment,
    plan_ucsc_bundle_acquisition,
    resolve_cached_ucsc_chain_resource_metadata,
    resolve_cached_ucsc_resource_bundle_metadata,
)
from .resource_files import ResourceReadProgressCallback
from .resources import UCSCResourceDiscoveryError, discover_ucsc_resources
from .reverse_mapping import (
    build_reverse_mapping_results_from_cached_chain,
    reverse_mapping_not_run,
    reverse_mapping_unavailable,
)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the ``assess-liftover`` command and return its process exit code."""

    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return _run(args, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr)
    except (
        OSError,
        UCSCResourceAcquisitionError,
        UCSCResourceDiscoveryError,
        ValueError,
    ) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="assess-liftover",
        description=(
            "Assess evidence supporting UCSC genomic coordinate liftover candidates."
        ),
    )
    parser.add_argument(
        "source_db", help="source UCSC database identifier, e.g. canFam3"
    )
    parser.add_argument(
        "target_db", help="target UCSC database identifier, e.g. canFam4"
    )
    parser.add_argument(
        "locus",
        nargs="?",
        help=(
            "source locus in UCSC-style 1-based inclusive coordinates, "
            "e.g. chr1:100-200; omit when using --bed"
        ),
    )
    parser.add_argument(
        "--bed",
        metavar="PATH",
        help=(
            "assess BED3-or-later batch input with native 0-based half-open "
            "coordinates; use '-' to read BED from stdin"
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="resource cache directory (default: platform user cache)",
    )
    parser.add_argument(
        "--evidence-tier",
        choices=("COMPARATIVE", "LIFTOVER-ONLY"),
        help=(
            "request one exact UCSC resource publication class instead of the "
            "default COMPARATIVE-preferred discovery/cache selection"
        ),
    )
    network_mode = parser.add_mutually_exclusive_group()
    network_mode.add_argument(
        "--refresh",
        action="store_true",
        help=(
            "contact UCSC and reacquire current resource bytes instead of reusing "
            "cache entries"
        ),
    )
    network_mode.add_argument(
        "--offline",
        action="store_true",
        help=(
            "guarantee zero provider access and require a complete verified local "
            "cache bundle"
        ),
    )
    parser.add_argument(
        "--acknowledge-ucsc-terms",
        action="store_true",
        help=(
            "confirm that the displayed UCSC terms have been reviewed and acknowledged"
        ),
    )
    parser.add_argument(
        "--accept-transfer-plan",
        action="store_true",
        help=(
            "confirm the displayed resource transfer plan without an interactive prompt"
        ),
    )
    output_mode = parser.add_mutually_exclusive_group()
    output_mode.add_argument(
        "--details",
        action="store_true",
        help="emit the full human-readable evidence, resource, and provenance dossier",
    )
    output_mode.add_argument(
        "--json",
        dest="json_output",
        action="store_true",
        help="emit the full machine-readable assessment report as JSON",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help=(
            "suppress nonessential progress messages "
            "(terms and transfer confirmations remain)"
        ),
    )
    parser.add_argument(
        "--context-bases",
        type=_context_window_bases_arg,
        help=(
            "override the automatic 101-bp local-context window for a 1-bp point "
            "query with an odd number of bases, e.g. 1001"
        ),
    )
    return parser


def _requested_evidence_tier(
    args: argparse.Namespace,
) -> EvidenceAvailabilityTier | None:
    if args.evidence_tier is None:
        return None
    return EvidenceAvailabilityTier(args.evidence_tier.replace("-", "_"))


def _run(
    args: argparse.Namespace,
    *,
    stdin: TextIO,
    stdout: TextIO,
    stderr: TextIO,
) -> int:
    source_assembly = ucsc_assembly_identifier(args.source_db)
    target_assembly = ucsc_assembly_identifier(args.target_db)
    cache_root = args.cache_dir or default_user_cache_root()
    requested_evidence_tier = _requested_evidence_tier(args)

    if args.bed is not None:
        if args.locus is not None:
            raise ValueError("provide either a single locus or --bed, not both")
        if args.context_bases is not None:
            raise ValueError(
                "--context-bases is not yet available with --bed batch assessment"
            )
        return _run_indexed_bed_batch(
            args,
            source_assembly=source_assembly,
            target_assembly=target_assembly,
            cache_root=cache_root,
            requested_evidence_tier=requested_evidence_tier,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
        )

    if args.locus is None:
        raise ValueError("a source locus is required unless --bed is provided")
    source_interval = parse_ucsc_locus(args.locus, assembly=source_assembly)
    if args.context_bases is not None and source_interval.length != 1:
        raise ValueError("--context-bases currently requires a 1-bp point query")

    cached_bundle = None
    chain_index = None
    unusable_index_sha256: str | None = None
    if not args.refresh:
        _status(
            "Checking/verifying local UCSC cache...",
            quiet=args.quiet,
            stderr=stderr,
        )

        if requested_evidence_tier is None:
            structural_bundle = resolve_cached_ucsc_resource_bundle_metadata(
                cache_root,
                args.source_db,
                args.target_db,
            )
        else:
            structural_bundle = resolve_cached_ucsc_resource_bundle_metadata(
                cache_root,
                args.source_db,
                args.target_db,
                evidence_tier=requested_evidence_tier,
            )
        if structural_bundle is not None:
            try:
                chain_index = load_cached_chain_index(
                    cache_root, structural_bundle.chain
                )
            except ChainIndexCorruptionError as exc:
                unusable_index_sha256 = structural_bundle.chain.sha256
                _status(
                    "Ignoring unusable cached chain index; "
                    f"using full chain verification/traversal ({exc}).",
                    quiet=args.quiet,
                    stderr=stderr,
                    indent=4,
                )

        cache_progress_display = _CacheVerificationProgressDisplay(stderr=stderr)
        cache_progress_callback: CacheVerificationProgressCallback | None = None
        if not args.quiet and _is_interactive_terminal(stderr):
            cache_progress_callback = cache_progress_display.update

        trusted_identifiers = (
            frozenset({structural_bundle.chain.sha256})
            if structural_bundle is not None and chain_index is not None
            else frozenset()
        )
        if cache_progress_callback is None:
            if trusted_identifiers and requested_evidence_tier is None:
                cached_bundle = load_cached_ucsc_resource_bundle_for_indexed_assessment(
                    cache_root,
                    args.source_db,
                    args.target_db,
                    trusted_artifact_sha256_identifiers=trusted_identifiers,
                )
            elif trusted_identifiers:
                cached_bundle = load_cached_ucsc_resource_bundle_for_indexed_assessment(
                    cache_root,
                    args.source_db,
                    args.target_db,
                    evidence_tier=requested_evidence_tier,
                    trusted_artifact_sha256_identifiers=trusted_identifiers,
                )
            elif requested_evidence_tier is None:
                cached_bundle = load_cached_ucsc_resource_bundle(
                    cache_root,
                    args.source_db,
                    args.target_db,
                )
            else:
                cached_bundle = load_cached_ucsc_resource_bundle(
                    cache_root,
                    args.source_db,
                    args.target_db,
                    evidence_tier=requested_evidence_tier,
                )
        elif trusted_identifiers and requested_evidence_tier is None:
            cached_bundle = load_cached_ucsc_resource_bundle_for_indexed_assessment(
                cache_root,
                args.source_db,
                args.target_db,
                progress_callback=cache_progress_callback,
                trusted_artifact_sha256_identifiers=trusted_identifiers,
            )
        elif trusted_identifiers:
            cached_bundle = load_cached_ucsc_resource_bundle_for_indexed_assessment(
                cache_root,
                args.source_db,
                args.target_db,
                evidence_tier=requested_evidence_tier,
                progress_callback=cache_progress_callback,
                trusted_artifact_sha256_identifiers=trusted_identifiers,
            )
        elif requested_evidence_tier is None:
            cached_bundle = load_cached_ucsc_resource_bundle(
                cache_root,
                args.source_db,
                args.target_db,
                progress_callback=cache_progress_callback,
            )
        else:
            cached_bundle = load_cached_ucsc_resource_bundle(
                cache_root,
                args.source_db,
                args.target_db,
                evidence_tier=requested_evidence_tier,
                progress_callback=cache_progress_callback,
            )
        if cached_bundle is not None:
            _status(
                "Using verified cached "
                f"{cached_bundle.evidence_tier.value} bundle; UCSC was not contacted. "
                "Use --refresh to check current provider resources.",
                quiet=args.quiet,
                stderr=stderr,
                indent=4,
            )

    if cached_bundle is None:
        if args.offline:
            requested = (
                ""
                if requested_evidence_tier is None
                else " " + requested_evidence_tier.value.replace("_", "-")
            )
            print(
                "error: --offline requires a complete verified cached"
                f"{requested} UCSC bundle for "
                f"{args.source_db}→{args.target_db} under {cache_root}",
                file=stderr,
            )
            return 1
        cached_bundle = _discover_and_acquire_bundle(
            args,
            cache_root=cache_root,
            stdin=stdin,
            stderr=stderr,
        )
        if cached_bundle is None:
            return 1

    if (
        chain_index is not None
        and chain_index.manifest.source_chain_sha256_identifier
        != cached_bundle.chain.sha256
    ):
        chain_index = None

    if chain_index is None and cached_bundle.chain.sha256 != unusable_index_sha256:
        try:
            chain_index = load_cached_chain_index(cache_root, cached_bundle.chain)
        except ChainIndexCorruptionError as exc:
            _status(
                "Ignoring unusable cached chain index; "
                f"using full chain traversal ({exc}).",
                quiet=args.quiet,
                stderr=stderr,
                indent=4,
            )
    if chain_index is not None:
        _status(
            "Using verified cached chain index for candidate lookup.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )

    _status("Assessing locus...", quiet=args.quiet, stderr=stderr)
    progress_display = _AssessmentProgressDisplay(
        cached_bundle,
        stderr=stderr,
        indexed_chain=chain_index is not None,
    )
    progress_callback: ResourceReadProgressCallback | None = None
    if not args.quiet and _is_interactive_terminal(stderr):
        progress_display.start()
        progress_callback = progress_display.update

    alignment_provenance = _ucsc_pair_dependency_provenance(
        args.source_db,
        args.target_db,
    )
    try:
        report = assess_ucsc_cached_bundle(
            source_interval,
            cached_bundle,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
            progress_callback=progress_callback,
            chain_index=chain_index,
        )
    except ChainIndexCorruptionError as exc:
        if chain_index is None:
            raise
        _status(
            f"Cached chain index query failed; retrying with full traversal ({exc}).",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        if progress_callback is not None:
            progress_display = _AssessmentProgressDisplay(
                cached_bundle, stderr=stderr, indexed_chain=False
            )
            progress_display.start()
            progress_callback = progress_display.update
        report = assess_ucsc_cached_bundle(
            source_interval,
            cached_bundle,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
            progress_callback=progress_callback,
            chain_index=None,
        )
        chain_index = None
    if progress_callback is not None:
        progress_display.finish(
            candidates_exist=bool(report.candidates),
        )

    report = _attach_cached_filtered_all_chain_comparison(
        report,
        args=args,
        cache_root=cache_root,
        stderr=stderr,
    )

    report = _attach_cached_point_query_context(
        report,
        args=args,
        cached_bundle=cached_bundle,
        chain_index=chain_index,
        stderr=stderr,
    )

    report = _attach_cached_reverse_mapping_context(
        report,
        args=args,
        cache_root=cache_root,
        stderr=stderr,
    )

    if args.json_output:
        rendered = render_assessment_json(report)
    elif args.details:
        rendered = render_assessment_details(report)
    else:
        rendered = render_assessment_summary(report)
    print(rendered, file=stdout)
    return 0


def _run_indexed_bed_batch(
    args: argparse.Namespace,
    *,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    cache_root: Path,
    requested_evidence_tier: EvidenceAvailabilityTier | None,
    stdin: TextIO,
    stdout: TextIO,
    stderr: TextIO,
) -> int:
    """Run cache-only, index-only BED batch assessment."""

    if args.details:
        raise ValueError(
            "--details is not yet available with --bed batch assessment; "
            "use default summary or --json"
        )
    if args.refresh:
        raise ValueError(
            "--refresh is not available with --bed batch assessment; batch mode is "
            "cache-only and requires a prepared exact-resource chain index"
        )

    if args.bed == "-":
        records = parse_bed_batch(stdin, assembly=source_assembly)
    else:
        bed_path = Path(args.bed)
        with bed_path.open("r", encoding="utf-8", newline="") as handle:
            records = parse_bed_batch(handle, assembly=source_assembly)

    _status(
        f"Loaded {len(records)} BED record(s).",
        quiet=args.quiet,
        stderr=stderr,
    )
    structural_chain = _resolve_preferred_cached_batch_chain(
        cache_root,
        args.source_db,
        args.target_db,
        requested_evidence_tier=requested_evidence_tier,
    )
    if structural_chain is None:
        requested = (
            "preferred COMPARATIVE/LIFTOVER-ONLY"
            if requested_evidence_tier is None
            else requested_evidence_tier.value.replace("_", "-")
        )
        raise ValueError(
            "--bed batch assessment requires a cached "
            f"{requested} chain for {args.source_db}→{args.target_db}; "
            "batch mode does not download resources implicitly"
        )

    try:
        chain_index = load_cached_chain_index(cache_root, structural_chain.chain)
    except ChainIndexCorruptionError as exc:
        tier = structural_chain.evidence_tier.value.replace("_", "-")
        raise ValueError(
            "--bed batch assessment requires a usable prepared chain index for "
            f"{tier}; rerun prepare-liftassess-index {args.source_db} "
            f"{args.target_db} --evidence-tier {tier} --rebuild ({exc})"
        ) from exc
    if chain_index is None:
        tier = structural_chain.evidence_tier.value.replace("_", "-")
        raise ValueError(
            "--bed batch assessment requires a prepared chain index for "
            f"{tier}; run prepare-liftassess-index {args.source_db} "
            f"{args.target_db} --evidence-tier {tier}"
        )

    chain_context = load_cached_ucsc_chain_resource(
        cache_root,
        args.source_db,
        args.target_db,
        evidence_tier=structural_chain.evidence_tier,
        trusted_artifact_sha256_identifiers=frozenset({structural_chain.chain.sha256}),
    )
    if chain_context is None:
        raise ValueError(
            "cached chain metadata changed or failed validation after prepared-index "
            "selection; batch assessment was not started"
        )

    _status(
        "Assessing BED batch with prepared chain index; provider access and "
        "whole-chain fallback are disabled.",
        quiet=args.quiet,
        stderr=stderr,
    )
    alignment_provenance = _ucsc_pair_dependency_provenance(
        args.source_db, args.target_db
    )
    try:
        result = run_indexed_chain_batch(
            records,
            chain_context,
            target_assembly=target_assembly,
            alignment_provenance=alignment_provenance,
            chain_index=chain_index,
        )
    except ChainIndexCorruptionError as exc:
        raise ValueError(
            "prepared chain index failed during batch lookup; no whole-chain fallback "
            f"was started ({exc})"
        ) from exc

    if args.json_output:
        rendered = render_indexed_chain_batch_json(result)
    else:
        rendered = render_indexed_chain_batch_summary(result)
    print(rendered, file=stdout)
    return 0


def _resolve_preferred_cached_batch_chain(
    cache_root: Path,
    source_db: str,
    target_db: str,
    *,
    requested_evidence_tier: EvidenceAvailabilityTier | None,
) -> CachedUCSCChainResource | None:
    tiers: tuple[EvidenceAvailabilityTier, ...]
    if requested_evidence_tier is None:
        tiers = (
            EvidenceAvailabilityTier.COMPARATIVE,
            EvidenceAvailabilityTier.LIFTOVER_ONLY,
        )
    else:
        tiers = (requested_evidence_tier,)
    for tier in tiers:
        resource = resolve_cached_ucsc_chain_resource_metadata(
            cache_root,
            source_db,
            target_db,
            evidence_tier=tier,
        )
        if resource is not None:
            return resource
    return None


def _context_window_bases_arg(value: str) -> int:
    try:
        window_bases = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "context window must be an integer number of bases"
        ) from exc
    if window_bases < 3:
        raise argparse.ArgumentTypeError("context window must be at least 3 bases")
    if window_bases % 2 == 0:
        raise argparse.ArgumentTypeError(
            "context window must contain an odd number of bases"
        )
    return window_bases


def _discover_and_acquire_bundle(
    args: argparse.Namespace,
    *,
    cache_root: Path,
    stdin: TextIO,
    stderr: TextIO,
) -> CachedUCSCResourceBundle | None:
    _status("Discovering UCSC resources...", quiet=args.quiet, stderr=stderr)
    requested_evidence_tier = _requested_evidence_tier(args)
    if requested_evidence_tier is None:
        discovered = discover_ucsc_resources(args.source_db, args.target_db)
    else:
        discovered = discover_ucsc_resources(
            args.source_db,
            args.target_db,
            evidence_tier=requested_evidence_tier,
        )
    if discovered is None:
        resource_label = (
            "UCSC resources"
            if requested_evidence_tier is None
            else requested_evidence_tier.value.replace("_", "-") + " UCSC resources"
        )
        print(
            f"error: no supported {resource_label} found for "
            f"{args.source_db}→{args.target_db}",
            file=stderr,
        )
        return None

    plan = plan_ucsc_bundle_acquisition(discovered)
    _print_terms(plan, stderr=stderr)
    terms_acknowledged = args.acknowledge_ucsc_terms or _confirm(
        "Acknowledge the UCSC terms shown above and continue to provider "
        "metadata inspection?",
        stdin=stdin,
        stderr=stderr,
    )
    if not terms_acknowledged:
        print("Cancelled before UCSC resource inspection or acquisition.", file=stderr)
        return None

    _status("Inspecting UCSC transfer metadata...", quiet=args.quiet, stderr=stderr)
    inspection = inspect_ucsc_bundle_transfer_plan(
        plan,
        terms_acknowledged=True,
    )
    _print_transfer_plan(
        plan,
        inspection,
        cache_root=cache_root,
        refresh=args.refresh,
        stderr=stderr,
    )
    transfer_acknowledged = args.accept_transfer_plan or _confirm(
        "Accept this transfer plan and continue?",
        stdin=stdin,
        stderr=stderr,
    )
    if not transfer_acknowledged:
        print("Cancelled before UCSC resource acquisition.", file=stderr)
        return None

    _status("Acquiring/verifying UCSC resources...", quiet=args.quiet, stderr=stderr)
    transfer_progress_display = _TransferProgressDisplay(
        plan,
        inspection,
        stderr=stderr,
    )
    transfer_progress_callback: UCSCBundleTransferProgressCallback | None = None
    if not args.quiet and _is_interactive_terminal(stderr):
        transfer_progress_display.start()
        transfer_progress_callback = transfer_progress_display.update

    return acquire_ucsc_resource_bundle(
        plan,
        cache_root,
        transfer_plan_acknowledged=True,
        terms_acknowledged=True,
        refresh=args.refresh,
        progress_callback=transfer_progress_callback,
    )


def default_user_cache_root() -> Path:
    """Return the platform-appropriate per-user liftAssess cache directory."""

    return _default_user_cache_root(
        platform_name=sys.platform,
        os_name=os.name,
        environ=os.environ,
        home=Path.home(),
    )


def _default_user_cache_root(
    *,
    platform_name: str,
    os_name: str,
    environ: Mapping[str, str],
    home: Path,
) -> Path:
    if platform_name == "darwin":
        return home / "Library" / "Caches" / "liftassess"
    if os_name == "nt":
        local_app_data = environ.get("LOCALAPPDATA")
        base = Path(local_app_data) if local_app_data else home / "AppData" / "Local"
        return base / "liftassess" / "Cache"

    xdg_cache_home = environ.get("XDG_CACHE_HOME")
    base = Path(xdg_cache_home) if xdg_cache_home else home / ".cache"
    return base / "liftassess"


def _human_evidence_tier(tier: EvidenceAvailabilityTier) -> str:
    return tier.value.replace("_", "-")


def _ucsc_pair_dependency_provenance(
    source_db: str, target_db: str
) -> ProvenanceSource:
    """Return the conservative pair-level dependency node used by the automatic CLI.

    Exact per-file bytes are content-addressed independently.  The CLI cannot infer from
    those bytes whether separately published UCSC resources came from one exact provider
    processing run, so this node only groups resources for the same UCSC assembly
    direction to prevent them from being presented as independent confirmation.
    """

    return ProvenanceSource(
        source_id=f"ucsc-pair:{source_db}:{target_db}",
        label=(
            f"Conservative UCSC {source_db}→{target_db} pair-level dependency group "
            "(liftAssess CLI; exact provider processing-run provenance not inferred)"
        ),
    )


def _attach_cached_filtered_all_chain_comparison(
    report: UCSCAssessmentReport,
    *,
    args: argparse.Namespace,
    cache_root: Path,
    stderr: TextIO,
) -> UCSCAssessmentReport:
    """Attach paired filtered/all-chain facts without provider access or full scans."""

    if report.evidence_tier is not EvidenceAvailabilityTier.COMPARATIVE:
        return report

    if args.refresh:
        _status(
            "Filtered/all-chain comparison not run during --refresh: the ordinary "
            "filtered liftOver chain was not refreshed automatically.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return report

    structural = resolve_cached_ucsc_chain_resource_metadata(
        cache_root,
        report.source_db,
        report.target_db,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
    )
    if structural is None:
        _status(
            "Filtered/all-chain comparison not run: no cached ordinary filtered "
            "liftOver chain is available; UCSC was not contacted.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return report

    try:
        filtered_index = load_cached_chain_index(cache_root, structural.chain)
    except ChainIndexCorruptionError as exc:
        _status(
            "Filtered/all-chain comparison not run: cached filtered-chain index is "
            f"unusable ({exc}). Rebuild it with prepare-liftassess-index "
            f"{report.source_db} {report.target_db} --evidence-tier LIFTOVER-ONLY "
            "--rebuild.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return report

    if filtered_index is None:
        _status(
            "Filtered/all-chain comparison not run: the ordinary filtered liftOver "
            "chain is cached but no prepared index is available. Run "
            f"prepare-liftassess-index {report.source_db} {report.target_db} "
            "--evidence-tier LIFTOVER-ONLY; no full filtered-chain scan was started.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return report

    filtered_chain = load_cached_ucsc_chain_resource(
        cache_root,
        report.source_db,
        report.target_db,
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        trusted_artifact_sha256_identifiers=frozenset({structural.chain.sha256}),
    )
    if filtered_chain is None:
        _status(
            "Filtered/all-chain comparison not run: cached filtered chain does not "
            "match the validated index identity.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return report

    _status(
        "Comparing ordinary filtered liftOver and all-chain placements using the "
        "prepared filtered-chain index...",
        quiet=args.quiet,
        stderr=stderr,
        indent=4,
    )
    try:
        return attach_filtered_all_chain_comparison(
            report,
            filtered_chain=filtered_chain,
            filtered_chain_index=filtered_index,
        )
    except ChainIndexCorruptionError as exc:
        _status(
            "Filtered/all-chain comparison not run: cached filtered-chain index "
            f"failed during lookup ({exc}). Rebuild it with prepare-liftassess-index "
            f"{report.source_db} {report.target_db} --evidence-tier LIFTOVER-ONLY "
            "--rebuild; no full filtered-chain scan was started.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return report
    except FilteredAllChainCorrespondenceError as exc:
        _status(
            "Filtered/all-chain comparison not run: filtered-chain geometry could "
            f"not be paired safely to the all-chain inventory ({exc}). The primary "
            "assessment remains valid; no comparative relationship was synthesized.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return report


def _attach_cached_point_query_context(
    report: UCSCAssessmentReport,
    *,
    args: argparse.Namespace,
    cached_bundle: CachedUCSCResourceBundle,
    chain_index: ChainIndex | None,
    stderr: TextIO,
) -> UCSCAssessmentReport:
    """Attach automatic point context without another whole-resource traversal."""

    if report.source_interval.length != 1:
        return report

    requested_window_bases = args.context_bases or DEFAULT_POINT_CONTEXT_BASES
    if chain_index is None:
        _status(
            "Point context not run: no prepared forward chain index is available; "
            "no additional full chain scan was started. Run prepare-liftassess-index "
            f"{report.source_db} {report.target_db} --evidence-tier "
            f"{_human_evidence_tier(report.evidence_tier)}.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_point_query_context(
            report,
            chain_context=CachedUCSCChainResource(
                source_db=cached_bundle.source_db,
                target_db=cached_bundle.target_db,
                evidence_tier=cached_bundle.evidence_tier,
                chain=cached_bundle.chain,
            ),
            chain_index=None,
            requested_window_bases=requested_window_bases,
        )

    _status(
        f"Assessing {requested_window_bases}-bp point context from the cached "
        "forward chain index...",
        quiet=args.quiet,
        stderr=stderr,
        indent=4,
    )
    chain_context = CachedUCSCChainResource(
        source_db=cached_bundle.source_db,
        target_db=cached_bundle.target_db,
        evidence_tier=cached_bundle.evidence_tier,
        chain=cached_bundle.chain,
    )
    try:
        enriched = attach_point_query_context(
            report,
            chain_context=chain_context,
            chain_index=chain_index,
            requested_window_bases=requested_window_bases,
        )
    except ChainIndexCorruptionError as exc:
        _status(
            "Point context not run: cached forward chain index failed during "
            f"context lookup ({exc}); no full chain fallback was started. Rebuild "
            "it with prepare-liftassess-index "
            f"{report.source_db} {report.target_db} --evidence-tier "
            f"{_human_evidence_tier(report.evidence_tier)} --rebuild.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_query_context_result(
            report,
            point_context_not_run(
                requested_window_bases=requested_window_bases,
                reason=QueryContextNotRunReason.INDEX_UNUSABLE,
            ),
        )

    context_profile = enriched.result_profile.query_context
    if (
        context_profile.not_run_reason
        is QueryContextNotRunReason.SOURCE_BOUNDS_UNAVAILABLE
    ):
        _status(
            "Point context not run: the prepared chain index does not provide a "
            f"source-sequence bound for {report.source_interval.sequence_name!r}; "
            "no full chain fallback was started.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
    return enriched


def _attach_cached_reverse_mapping_context(
    report: UCSCAssessmentReport,
    *,
    args: argparse.Namespace,
    cache_root: Path,
    stderr: TextIO,
) -> UCSCAssessmentReport:
    """Attach automatic reverse facts without provider access or full scans."""

    if not report.candidates:
        return report

    unavailable = tuple(
        reverse_mapping_unavailable(candidate) for candidate in report.candidates
    )
    not_run = tuple(
        reverse_mapping_not_run(candidate) for candidate in report.candidates
    )
    reverse_source_db = report.target_db
    reverse_target_db = report.source_db

    if args.refresh:
        _status(
            "Reverse mapping not run during --refresh: reverse-direction resources "
            "were not refreshed automatically.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_reverse_mapping_results(report, not_run)

    structural = resolve_cached_ucsc_chain_resource_metadata(
        cache_root,
        reverse_source_db,
        reverse_target_db,
        evidence_tier=report.evidence_tier,
    )
    if structural is None:
        _status(
            "Reverse mapping unavailable: no cached reverse-direction chain with "
            f"matching {report.evidence_tier.value.replace('_', '-')} publication "
            "class; UCSC was not contacted.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_reverse_mapping_results(report, unavailable)

    try:
        reverse_index = load_cached_chain_index(cache_root, structural.chain)
    except ChainIndexCorruptionError as exc:
        _status(
            "Reverse mapping not run: cached reverse chain index is unusable "
            f"({exc}). Rebuild it with prepare-liftassess-index "
            f"{reverse_source_db} {reverse_target_db} --evidence-tier "
            f"{_human_evidence_tier(report.evidence_tier)} --rebuild.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_reverse_mapping_results(report, not_run)

    if reverse_index is None:
        _status(
            "Reverse mapping not run: the matching reverse chain is cached but no "
            "prepared index is available. Run prepare-liftassess-index "
            f"{reverse_source_db} {reverse_target_db} --evidence-tier "
            f"{_human_evidence_tier(report.evidence_tier)}; no full reverse-chain "
            "scan was started.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_reverse_mapping_results(report, not_run)

    reverse_chain = load_cached_ucsc_chain_resource(
        cache_root,
        reverse_source_db,
        reverse_target_db,
        evidence_tier=report.evidence_tier,
        trusted_artifact_sha256_identifiers=frozenset({structural.chain.sha256}),
    )
    if reverse_chain is None:
        _status(
            "Reverse mapping unavailable: cached reverse chain does not match the "
            "validated index identity.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_reverse_mapping_results(report, unavailable)

    _status(
        "Assessing actual reverse mapping from cached indexed "
        f"{reverse_source_db}→{reverse_target_db} chain...",
        quiet=args.quiet,
        stderr=stderr,
        indent=4,
    )
    reverse_alignment = _ucsc_pair_dependency_provenance(
        reverse_source_db, reverse_target_db
    )
    try:
        reverse_results = build_reverse_mapping_results_from_cached_chain(
            report.candidates,
            reverse_chain,
            reverse_alignment_provenance=reverse_alignment,
            chain_index=reverse_index,
        )
    except ChainIndexCorruptionError as exc:
        _status(
            "Reverse mapping not run: cached reverse chain index failed during "
            f"lookup ({exc}). Rebuild it with prepare-liftassess-index "
            f"{reverse_source_db} {reverse_target_db} --evidence-tier "
            f"{_human_evidence_tier(report.evidence_tier)} --rebuild; no full "
            "reverse-chain scan "
            "was started.",
            quiet=args.quiet,
            stderr=stderr,
            indent=4,
        )
        return attach_reverse_mapping_results(report, not_run)

    return attach_reverse_mapping_results(
        report,
        reverse_results,
        reverse_chain=reverse_chain,
        reverse_alignment_provenance=reverse_alignment,
    )


_PROGRESS_BAR_WIDTH = 20
_ASSESSMENT_PROGRESS_LABELS = {
    UCSCBundleResourceRole.CHAIN: "Chain",
    UCSCBundleResourceRole.NET: "Net",
    UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN: "Reciprocal-best",
}
_TRANSFER_PROGRESS_LABELS = {
    UCSCBundleResourceRole.CHAIN: "Chain",
    UCSCBundleResourceRole.NET: "Net",
    UCSCBundleResourceRole.SYNTENIC_NET: "Syntenic net",
    UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN: "Rbest chain",
    UCSCBundleResourceRole.RECIPROCAL_BEST_NET: "Rbest net",
}


class _TransferProgressDisplay:
    """Render measured, resume-aware acquisition progress for each bundle resource."""

    def __init__(
        self,
        plan: UCSCBundleAcquisitionPlan,
        inspection: UCSCBundleTransferInspection,
        *,
        stderr: TextIO,
    ) -> None:
        self._stderr = stderr
        self._roles = tuple(item.role for item in plan.items)
        self._bytes_complete = {role: 0 for role in self._roles}
        self._total_bytes = {
            item.role: item.identity_content_length_bytes for item in inspection.items
        }
        self._cache_hit = {role: False for role in self._roles}
        self._last_rows = {role: "" for role in self._roles}
        self._started = False

    def start(self) -> None:
        self._started = True
        self._render(initial=True)

    def update(
        self,
        role: UCSCBundleResourceRole,
        bytes_complete: int,
        total_bytes: int | None,
        cache_hit: bool,
    ) -> None:
        if role not in self._bytes_complete:
            return
        bounded = max(bytes_complete, 0)
        if total_bytes is not None:
            total_bytes = max(total_bytes, 0)
            bounded = min(bounded, total_bytes)
            self._total_bytes[role] = total_bytes
        self._bytes_complete[role] = bounded
        self._cache_hit[role] = cache_hit
        self._render()

    def _row(self, role: UCSCBundleResourceRole) -> str:
        return _progress_row(
            _TRANSFER_PROGRESS_LABELS[role],
            bytes_read=self._bytes_complete[role],
            total_bytes=self._total_bytes[role],
            cached=self._cache_hit[role],
        )

    def _render(self, *, initial: bool = False) -> None:
        if not self._started:
            return
        rows = {role: self._row(role) for role in self._roles}
        if not initial and rows == self._last_rows:
            return
        if not initial:
            self._stderr.write(f"\x1b[{len(self._roles)}A")
        for role in self._roles:
            self._stderr.write("\x1b[2K")
            self._stderr.write(rows[role])
            self._stderr.write("\n")
        self._stderr.flush()
        self._last_rows = rows


class _CacheVerificationProgressDisplay:
    """Render one aggregate row for measured cached-artifact SHA-256 work."""

    def __init__(self, *, stderr: TextIO) -> None:
        self._stderr = stderr
        self._started = False
        self._last_percent = -1
        self._last_complete = False

    def update(
        self,
        bytes_hashed: int,
        total_bytes: int,
        complete: bool,
    ) -> None:
        bounded_total = max(total_bytes, 0)
        bounded_hashed = min(max(bytes_hashed, 0), bounded_total)
        percent = _progress_percent(bounded_hashed, bounded_total)
        if not complete and percent == 100:
            # Reading the final byte is not the same as validating the expected digest.
            # Hold the visual completion state until the cache loader confirms that
            # every required artifact passed its SHA-256 identity check.
            percent = 99
        if percent == self._last_percent and complete == self._last_complete:
            return

        if self._started:
            self._stderr.write("\x1b[1A")
        self._stderr.write("\x1b[2K")
        self._stderr.write(
            _progress_row(
                "Cache verification",
                bytes_read=bounded_hashed,
                total_bytes=bounded_total,
                percent_override=percent,
            )
        )
        self._stderr.write("\n")
        self._stderr.flush()
        self._started = True
        self._last_percent = percent
        self._last_complete = complete


class _AssessmentProgressDisplay:
    """Render measured raw-byte progress for resources consumed by assessment."""

    def __init__(
        self,
        bundle: CachedUCSCResourceBundle,
        *,
        stderr: TextIO,
        indexed_chain: bool = False,
    ) -> None:
        self._stderr = stderr
        self._indexed_chain = indexed_chain
        self._resources = {
            UCSCBundleResourceRole.CHAIN: bundle.chain,
        }
        if bundle.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
            assert bundle.net is not None
            assert bundle.reciprocal_best_chain is not None
            self._resources[UCSCBundleResourceRole.NET] = bundle.net
            self._resources[UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN] = (
                bundle.reciprocal_best_chain
            )
        self._bytes_read = {role: 0 for role in self._resources}
        self._last_percent = {role: -1 for role in self._resources}
        self._started = False
        self._finished = False

    def start(self) -> None:
        self._started = True
        self._render(initial=True)

    def update(
        self,
        role: UCSCBundleResourceRole,
        bytes_read: int,
        total_bytes: int,
    ) -> None:
        resource = self._resources.get(role)
        if resource is None:
            return
        if total_bytes != resource.size_bytes:
            raise ValueError(
                "assessment progress total does not match cached resource size"
            )
        bounded = min(max(bytes_read, 0), total_bytes)
        percent = _progress_percent(bounded, total_bytes)
        self._bytes_read[role] = bounded
        if percent == self._last_percent[role] and bounded != total_bytes:
            return
        self._last_percent[role] = percent
        self._render()

    def finish(self, *, candidates_exist: bool) -> None:
        self._finished = True
        if not candidates_exist:
            for role in (
                UCSCBundleResourceRole.NET,
                UCSCBundleResourceRole.RECIPROCAL_BEST_CHAIN,
            ):
                if role in self._resources and self._bytes_read[role] == 0:
                    self._last_percent[role] = -2
        self._render()

    def _render(self, *, initial: bool = False) -> None:
        if not self._started:
            return
        row_count = len(self._resources)
        if not initial:
            self._stderr.write(f"\x1b[{row_count}A")
        for role, resource in self._resources.items():
            self._stderr.write("\x1b[2K")
            self._stderr.write(
                _progress_row(
                    _ASSESSMENT_PROGRESS_LABELS[role],
                    bytes_read=self._bytes_read[role],
                    total_bytes=resource.size_bytes,
                    not_used=(self._finished and self._last_percent[role] == -2),
                    indexed=(
                        role is UCSCBundleResourceRole.CHAIN and self._indexed_chain
                    ),
                )
            )
            self._stderr.write("\n")
        self._stderr.flush()


def _progress_row(
    label: str,
    *,
    bytes_read: int,
    total_bytes: int | None,
    not_used: bool = False,
    cached: bool = False,
    indexed: bool = False,
    percent_override: int | None = None,
) -> str:
    if not_used:
        return f"  {label:<18} [{'—' * _PROGRESS_BAR_WIDTH}]  --   not used"
    if indexed:
        return f"  {label:<18} [{'█' * _PROGRESS_BAR_WIDTH}]  --   indexed"
    if cached:
        amount = _format_progress_bytes(bytes_read)
        return f"  {label:<18} [{'█' * _PROGRESS_BAR_WIDTH}]  --   cached ({amount})"
    if total_bytes is None:
        amount = (
            "pending"
            if bytes_read == 0
            else f"{_format_progress_bytes(bytes_read)} complete"
        )
        return f"  {label:<18} [{'—' * _PROGRESS_BAR_WIDTH}]  --   {amount}"
    percent = (
        _progress_percent(bytes_read, total_bytes)
        if percent_override is None
        else min(100, max(0, percent_override))
    )
    filled = min(_PROGRESS_BAR_WIDTH, percent * _PROGRESS_BAR_WIDTH // 100)
    bar = "█" * filled + "-" * (_PROGRESS_BAR_WIDTH - filled)
    if bytes_read == 0:
        amount = "pending"
    else:
        amount = (
            f"{_format_progress_bytes(bytes_read)} / "
            f"{_format_progress_bytes(total_bytes)}"
        )
    return f"  {label:<18} [{bar}] {percent:3d}%  {amount}"


def _progress_percent(bytes_read: int, total_bytes: int) -> int:
    if total_bytes <= 0:
        return 100
    return min(100, int(bytes_read * 100 / total_bytes))


def _format_progress_bytes(size: int) -> str:
    if size < 1024:
        return f"{size} B"
    units = ("KiB", "MiB", "GiB", "TiB")
    value = float(size)
    for unit in units:
        value /= 1024.0
        if value < 1024.0 or unit == units[-1]:
            return f"{value:.2f} {unit}"
    raise AssertionError("unreachable")


def _is_interactive_terminal(stream: TextIO) -> bool:
    try:
        return stream.isatty()
    except (AttributeError, OSError):
        return False


def _print_terms(plan: UCSCBundleAcquisitionPlan, *, stderr: TextIO) -> None:
    print(
        "UCSC terms to review before provider metadata access/acquisition:",
        file=stderr,
    )
    print(f"  General terms: {plan.items[0].terms.general_terms_url}", file=stderr)

    seen_directories: set[str] = set()
    for item in plan.items:
        directory_terms = item.terms.directory_terms_url
        if directory_terms not in seen_directories:
            print(f"  Directory terms: {directory_terms}", file=stderr)
            seen_directories.add(directory_terms)
        if item.terms.restricted_liftover_chain:
            print(
                "  Restricted liftOver-chain terms apply to the planned chain "
                "resource.",
                file=stderr,
            )


def _print_transfer_plan(
    plan: UCSCBundleAcquisitionPlan,
    inspection: UCSCBundleTransferInspection,
    *,
    cache_root: Path,
    refresh: bool,
    stderr: TextIO,
) -> None:
    print(
        f"Transfer plan: {plan.evidence_tier.value} ({len(plan.items)} resource(s))",
        file=stderr,
    )
    for plan_item, inspected_item in zip(plan.items, inspection.items, strict=True):
        metadata = inspected_item.metadata
        if metadata.content_length_bytes is None:
            size_text = "HTTP Content-Length unavailable"
        else:
            size_text = (
                f"HTTP Content-Length {_format_bytes(metadata.content_length_bytes)}"
            )
        if (
            metadata.content_encoding is not None
            and metadata.content_encoding.casefold() != "identity"
        ):
            size_text += f"; Content-Encoding {metadata.content_encoding}"
        print(f"  {plan_item.role.value}: {plan_item.url} ({size_text})", file=stderr)

    total = inspection.total_content_length_bytes
    total_text = _format_bytes(total) if total is not None else "unknown"
    print(
        f"  Provider-advertised total identity resource size: {total_text}",
        file=stderr,
    )
    print(f"  Cache: {cache_root}", file=stderr)
    print(f"  Refresh cached resources: {'yes' if refresh else 'no'}", file=stderr)
    if not refresh:
        print(
            "  Verified cache hits may avoid resource-body transfer.",
            file=stderr,
        )


def _format_bytes(size_bytes: int) -> str:
    if size_bytes < 1024:
        return f"{size_bytes} B"

    value = float(size_bytes)
    units = ("KiB", "MiB", "GiB", "TiB")
    for unit in units:
        value /= 1024.0
        if value < 1024.0 or unit == units[-1]:
            return f"{value:.1f} {unit}"
    raise AssertionError("unreachable byte-size formatting branch")


def _confirm(prompt: str, *, stdin: TextIO, stderr: TextIO) -> bool:
    print(f"{prompt} [y/N] ", end="", file=stderr, flush=True)
    response = stdin.readline()
    if response == "":
        print(file=stderr)
        return False
    return response.strip().casefold() in {"y", "yes"}


def _status(
    message: str,
    *,
    quiet: bool,
    stderr: TextIO,
    indent: int = 0,
) -> None:
    if not quiet:
        print(" " * indent + message, file=stderr)


if __name__ == "__main__":  # pragma: no cover - console-script path is primary.
    raise SystemExit(main())
