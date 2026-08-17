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

from .cli_input import parse_ucsc_locus, ucsc_assembly_identifier
from .models import ProvenanceSource
from .orchestration import assess_ucsc_cached_bundle
from .reporting import render_assessment_summary
from .resource_cache import (
    UCSCBundleAcquisitionPlan,
    UCSCBundleTransferInspection,
    UCSCResourceAcquisitionError,
    acquire_ucsc_resource_bundle,
    inspect_ucsc_bundle_transfer_plan,
    plan_ucsc_bundle_acquisition,
)
from .resources import UCSCResourceDiscoveryError, discover_ucsc_resources


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
        help=(
            "source locus in UCSC-style 1-based inclusive coordinates, "
            "e.g. chr1:100-200"
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="resource cache directory (default: platform user cache)",
    )
    parser.add_argument(
        "--refresh",
        action="store_true",
        help=(
            "contact UCSC and reacquire current resource bytes instead of reusing "
            "cache entries"
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
    parser.add_argument(
        "--quiet",
        action="store_true",
        help=(
            "suppress nonessential progress messages "
            "(terms and transfer confirmations remain)"
        ),
    )
    return parser


def _run(
    args: argparse.Namespace,
    *,
    stdin: TextIO,
    stdout: TextIO,
    stderr: TextIO,
) -> int:
    source_assembly = ucsc_assembly_identifier(args.source_db)
    target_assembly = ucsc_assembly_identifier(args.target_db)
    source_interval = parse_ucsc_locus(args.locus, assembly=source_assembly)
    cache_root = args.cache_dir or default_user_cache_root()

    _status("Discovering UCSC resources...", quiet=args.quiet, stderr=stderr)
    discovered = discover_ucsc_resources(args.source_db, args.target_db)
    if discovered is None:
        print(
            "error: no supported UCSC resources found for "
            f"{args.source_db}→{args.target_db}",
            file=stderr,
        )
        return 1

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
        return 1

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
        return 1

    _status("Acquiring/verifying UCSC resources...", quiet=args.quiet, stderr=stderr)
    cached_bundle = acquire_ucsc_resource_bundle(
        plan,
        cache_root,
        transfer_plan_acknowledged=True,
        terms_acknowledged=True,
        refresh=args.refresh,
    )

    _status("Assessing locus...", quiet=args.quiet, stderr=stderr)
    report = assess_ucsc_cached_bundle(
        source_interval,
        cached_bundle,
        target_assembly=target_assembly,
        alignment_provenance=_ucsc_pair_lineage_provenance(
            args.source_db,
            args.target_db,
        ),
    )
    print(render_assessment_summary(report.assessment), file=stdout)
    return 0


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


def _ucsc_pair_lineage_provenance(source_db: str, target_db: str) -> ProvenanceSource:
    """Return the conservative shared-dependency node used by the automatic CLI.

    The CLI cannot infer exact upstream process provenance from downloaded bytes.  It
    therefore groups consumed resources for one UCSC assembly direction under a single
    pair-level lineage node so related observations are not presented as independent
    confirmation.  This deliberately conservative dependency statement does not replace
    the exact content-addressed provenance recorded for each consumed file.
    """

    return ProvenanceSource(
        source_id=f"ucsc-pair:{source_db}:{target_db}",
        label=(
            f"Conservative shared UCSC {source_db}→{target_db} resource lineage "
            "(liftAssess CLI dependency grouping)"
        ),
    )


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


def _status(message: str, *, quiet: bool, stderr: TextIO) -> None:
    if not quiet:
        print(message, file=stderr)


if __name__ == "__main__":  # pragma: no cover - console-script path is primary.
    raise SystemExit(main())
