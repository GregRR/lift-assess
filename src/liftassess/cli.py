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
from .models import EvidenceAvailabilityTier, ProvenanceSource
from .orchestration import assess_ucsc_cached_bundle
from .reporting import (
    render_assessment_details,
    render_assessment_json,
    render_assessment_summary,
)
from .resource_cache import (
    CachedUCSCResourceBundle,
    CacheVerificationProgressCallback,
    UCSCBundleAcquisitionPlan,
    UCSCBundleResourceRole,
    UCSCBundleTransferInspection,
    UCSCBundleTransferProgressCallback,
    UCSCResourceAcquisitionError,
    acquire_ucsc_resource_bundle,
    inspect_ucsc_bundle_transfer_plan,
    load_cached_ucsc_resource_bundle,
    plan_ucsc_bundle_acquisition,
)
from .resource_files import ResourceReadProgressCallback
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

    cached_bundle = None
    if not args.refresh:
        _status(
            "Checking/verifying local UCSC cache...",
            quiet=args.quiet,
            stderr=stderr,
        )
        cache_progress_display = _CacheVerificationProgressDisplay(stderr=stderr)
        cache_progress_callback: CacheVerificationProgressCallback | None = None
        if not args.quiet and _is_interactive_terminal(stderr):
            cache_progress_callback = cache_progress_display.update

        if cache_progress_callback is None:
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
                progress_callback=cache_progress_callback,
            )
        if cached_bundle is not None:
            _status(
                "Using verified cached "
                f"{cached_bundle.evidence_tier.value} bundle; UCSC was not contacted. "
                "Use --refresh to check current provider resources.",
                quiet=args.quiet,
                stderr=stderr,
            )

    if cached_bundle is None:
        if args.offline:
            print(
                "error: --offline requires a complete verified cached UCSC bundle for "
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

    _status("Assessing locus...", quiet=args.quiet, stderr=stderr)
    progress_display = _AssessmentProgressDisplay(cached_bundle, stderr=stderr)
    progress_callback: ResourceReadProgressCallback | None = None
    if not args.quiet and _is_interactive_terminal(stderr):
        progress_display.start()
        progress_callback = progress_display.update

    report = assess_ucsc_cached_bundle(
        source_interval,
        cached_bundle,
        target_assembly=target_assembly,
        alignment_provenance=_ucsc_pair_lineage_provenance(
            args.source_db,
            args.target_db,
        ),
        progress_callback=progress_callback,
    )
    if progress_callback is not None:
        progress_display.finish(
            candidates_exist=bool(report.candidates),
        )
    if args.json_output:
        rendered = render_assessment_json(report)
    elif args.details:
        rendered = render_assessment_details(report)
    else:
        rendered = render_assessment_summary(report)
    print(rendered, file=stdout)
    return 0


def _discover_and_acquire_bundle(
    args: argparse.Namespace,
    *,
    cache_root: Path,
    stdin: TextIO,
    stderr: TextIO,
) -> CachedUCSCResourceBundle | None:
    _status("Discovering UCSC resources...", quiet=args.quiet, stderr=stderr)
    discovered = discover_ucsc_resources(args.source_db, args.target_db)
    if discovered is None:
        print(
            "error: no supported UCSC resources found for "
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

    def __init__(self, bundle: CachedUCSCResourceBundle, *, stderr: TextIO) -> None:
        self._stderr = stderr
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
    percent_override: int | None = None,
) -> str:
    if not_used:
        return f"  {label:<18} [{'—' * _PROGRESS_BAR_WIDTH}]  --   not used"
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


def _status(message: str, *, quiet: bool, stderr: TextIO) -> None:
    if not quiet:
        print(message, file=stderr)


if __name__ == "__main__":  # pragma: no cover - console-script path is primary.
    raise SystemExit(main())
