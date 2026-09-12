"""Human-readable and machine-readable factual result reporting.

Both renderers consume the same derived result profile. Exact candidate geometry,
evidence, resource identity, and provenance remain available alongside the profile;
no renderer assigns an aggregate verdict, confidence score, candidate rank, or
biological correctness claim.
"""

import json
import re
from textwrap import wrap
from urllib.parse import urlencode

from .chain import chain_id_from_candidate_id
from .comparative_inventory import FilteredAllChainInventoryState
from .models import (
    AssemblyIdentifier,
    ChainGapSummary,
    EvidenceAvailabilityTier,
    EvidenceObservation,
    EvidenceValue,
    GenomicInterval,
    MappingCoverageSummary,
    MappingOrientation,
    NetHierarchySummary,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestMembershipSummary,
)
from .orchestration import UCSCAssessmentReport, UCSCAssessmentResource
from .query_context import QueryContextState
from .resource_cache import CachedResource
from .result_profile import (
    CandidateResultProfile,
    CandidateReverseMappingProfile,
    ComparativePlacementProfile,
    ComparativeRelationshipProfile,
    ComparativeRelationshipState,
    ExternalContextProfile,
    ExternalContextState,
    FactualHeadline,
    QueryContextFinding,
    QueryContextProfile,
    ResultProfile,
    TargetRoleState,
    TargetSequenceRoleProfile,
)
from .reverse_mapping import (
    CandidateReverseMappingResult,
    ReverseCheckState,
    ReverseRelationshipState,
)
from .segmental_duplication import (
    CandidateSegmentalDuplicationOverlap,
    UCSCSegmentalDuplicationContextResult,
    UCSCSegmentalDuplicationOverlap,
    UCSCSegmentalDuplicationRecord,
)

_BIOLOGICAL_CORRECTNESS_CAVEAT = "This does not establish biological correctness."
_JSON_SCHEMA_VERSION = 2
_JSON_REPORT_TYPE = "liftassess.ucsc_result"
_JSON_INTERVAL_COORDINATE_SYSTEM = "0-based-half-open"
_DEFAULT_INLINE_PROJECTION_LIMIT = 4


def format_display_interval(interval: GenomicInterval) -> str:
    """Format a canonical interval as explicit 1-based inclusive coordinates."""

    if interval.length <= 0:
        raise ValueError("displayed genomic intervals must span at least one base")
    return (
        f"{interval.sequence_name}:{interval.start + 1}-{interval.end} "
        "(1-based inclusive)"
    )


def render_invalid_source_coordinate(
    database: str,
    source_interval: GenomicInterval,
    *,
    sequence_length: int,
) -> str:
    """Render an out-of-bounds source-coordinate error for terminal users."""

    if sequence_length <= 0:
        raise ValueError("source sequence length must be positive")
    if source_interval.end <= sequence_length:
        raise ValueError(
            "invalid-coordinate rendering requires an out-of-bounds interval"
        )

    excess = source_interval.end - sequence_length
    sequence = source_interval.sequence_name
    boundary_phrase = (
        "coordinate is" if source_interval.length == 1 else "interval ends"
    )
    return "\n".join(
        (
            "* INVALID SOURCE COORDINATE *",
            "",
            "Source:",
            f"    {_human_interval_text(database, source_interval)}",
            "",
            f"{database} {sequence} length:",
            f"    {sequence_length:,} bp",
            "",
            "= INPUT ERROR =",
            "",
            (
                f"The requested {boundary_phrase} {excess:,} bp beyond the end of "
                f"{database} {sequence}."
            ),
            "",
            "liftOver was not attempted.",
            "",
            "= CHECK PERFORMED =",
            "",
            f"    Source coordinate checked against UCSC {database} sequence size",
        )
    )


def _wrap_summary_paragraphs(lines: list[str], *, width: int = 88) -> list[str]:
    """Wrap explanatory summary prose while preserving labels, indentation, and URLs."""

    wrapped: list[str] = []
    for line in lines:
        if (
            not line
            or line.startswith("    ")
            or line.endswith(":")
            or "https://" in line
        ):
            wrapped.append(line)
            continue
        wrapped.extend(
            wrap(
                line,
                width=width,
                break_long_words=False,
                break_on_hyphens=False,
            )
            or [line]
        )
    return wrapped


def _local_why_lines(paragraphs: list[str]) -> list[str]:
    """Render finding-specific interpretation immediately below that finding."""

    if not paragraphs:
        return []
    lines = ["", "    Why this matters:"]
    for paragraph in paragraphs:
        lines.extend(
            "        " + line
            for line in wrap(
                paragraph,
                width=80,
                break_long_words=False,
                break_on_hyphens=False,
            )
        )
    return lines


def render_assessment_summary(report: UCSCAssessmentReport) -> str:
    """Render the progressive-disclosure default factual result summary."""

    profile = report.result_profile
    if len(profile.candidate_profiles) <= 1:
        return _render_single_mapping_summary(report)
    return _render_multiple_mapping_summary(report)


def _render_multiple_mapping_summary(report: UCSCAssessmentReport) -> str:
    """Render multiple mappings using standard liftOver terminology."""

    profile = report.result_profile
    lines = [
        f"* {_summary_headline_text(report)} *",
        "",
        "Source:",
        f"    {_human_interval_text(report.source_db, profile.source_interval)}",
        "",
        "= KEY FINDINGS =",
        "",
        "Mapping:",
        *_multiple_mapping_overview_lines(report),
    ]

    reverse_lines = _multiple_reverse_liftover_lines(report)
    if reverse_lines:
        lines.extend(("", *reverse_lines))

    comparative_lines = _multiple_comparative_finding_lines(report)
    if comparative_lines:
        lines.extend(("", *comparative_lines))

    context_lines = _single_flanking_interval_lines(report)
    if context_lines:
        lines.extend(("", *context_lines))

    duplication_lines = _multiple_segmental_duplication_lines(report)
    if duplication_lines:
        lines.extend(("", *duplication_lines))

    target_metadata_lines = _single_target_sequence_metadata_lines(report)
    if target_metadata_lines:
        lines.extend(("", *target_metadata_lines))

    why_lines, next_step_lines = _multiple_mapping_guidance_lines(report)
    if why_lines:
        lines.extend(
            ("", "= WHY THIS MATTERS =", "", *_wrap_summary_paragraphs(why_lines))
        )
    if next_step_lines:
        lines.extend(
            ("", "= NEXT STEP =", "", *_wrap_summary_paragraphs(next_step_lines))
        )
    lines.extend(
        (
            "",
            "Details:",
            "    Use --details for every mapping, alignment blocks, evidence,",
            "    and provenance, or --json for machine-readable output.",
        )
    )
    return "\n".join(lines)


def _multiple_mapping_overview_lines(report: UCSCAssessmentReport) -> list[str]:
    profile = report.result_profile
    count = len(profile.candidate_profiles)
    complete_count = sum(
        candidate.covered_source_bases == candidate.source_bases
        for candidate in profile.candidate_profiles
    )
    source_bases = profile.source_bases

    lines: list[str] = []
    if complete_count == count:
        if report.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
            lines.append(
                f"    The UCSC all-chain alignments contain {count} complete mappings."
            )
        else:
            lines.append(f"    {count} complete liftOver mappings were found.")
        lines.append(f"    Each maps {source_bases}/{source_bases} input bases.")
    else:
        lines.append(
            "    Best single mapping: "
            f"{profile.maximum_candidate_covered_source_bases}/{source_bases} "
            "input bases mapped."
        )
        lines.append(f"    {count} liftOver mappings were found.")
        if complete_count:
            noun = "mapping" if complete_count == 1 else "mappings"
            lines.append(
                f"    {complete_count} {noun} cover the entire source interval."
            )
        if (
            profile.union_covered_source_bases
            != profile.maximum_candidate_covered_source_bases
        ):
            lines.append(
                "    Across all mappings, "
                f"{profile.union_covered_source_bases}/{source_bases} source bases "
                "are represented."
            )

    target_sequences = {
        candidate.target_interval.sequence_name for candidate in report.candidates
    }
    if len(target_sequences) > 1:
        if all(_is_standard_ucsc_chromosome_name(name) for name in target_sequences):
            noun = "chromosomes"
        else:
            noun = "target sequences"
        lines.append(
            f"    The mappings span {len(target_sequences)} {report.target_db} {noun}."
        )

    if count <= _DEFAULT_INLINE_PROJECTION_LIMIT:
        lines.extend(("", "Mappings:"))
        for candidate, candidate_profile in zip(
            report.candidates,
            profile.candidate_profiles,
            strict=True,
        ):
            lines.extend(
                (
                    (
                        "    "
                        + _human_interval_text(
                            report.target_db, candidate.target_interval
                        )
                    ),
                    (
                        "        Input bases mapped: "
                        f"{candidate_profile.covered_source_bases}/"
                        f"{candidate_profile.source_bases}"
                    ),
                    "        Orientation:",
                    f"            {candidate.orientation.value.lower()}",
                )
            )
    else:
        lines.append(f"    Use --details to view all {count} mappings.")
    return lines


def _multiple_reverse_liftover_lines(
    report: UCSCAssessmentReport,
) -> list[str]:
    state = report.result_profile.scope.reverse_result
    if state is ReverseCheckState.UNAVAILABLE:
        return [
            "Reverse liftOver:",
            "    Reverse liftOver was unavailable for this assessment.",
        ]
    return []


def _multiple_comparative_finding_lines(
    report: UCSCAssessmentReport,
) -> list[str]:
    profile = report.result_profile.comparative_relationship
    if profile.state is ComparativeRelationshipState.NOT_ASSESSED:
        return []
    comparison = report.filtered_all_chain_comparison
    if comparison is None:
        raise ValueError("comparative summary requires filtered/all-chain comparison")

    all_chain_count = len(comparison.all_chain_candidates)
    lines = ["Comparative UCSC evidence:"]

    if profile.state is ComparativeRelationshipState.FAVORS_ONE_PLACEMENT:
        favored_id = profile.favored_candidate_id
        if favored_id is None:
            raise ValueError("favored comparative relationship requires candidate ID")
        favored = _report_candidate_for_id(report, favored_id)
        additional = len(profile.additional_all_chain_candidate_ids)
        lines.extend(
            (
                (
                    f"    The standard {report.source_db}→{report.target_db} liftOver "
                    f"chain retains one of the {all_chain_count} mappings:"
                ),
                (
                    "        "
                    + _human_interval_text(report.target_db, favored.target_interval)
                ),
                "        Orientation:",
                f"            {favored.orientation.value.lower()}",
            )
        )
        if additional:
            lines.extend(
                (
                    (
                        f"    The other {additional} complete "
                        f"{'mapping is' if additional == 1 else 'mappings are'} present in"
                    ),
                    (
                        "    the all-chain alignments but not in the standard "
                        "liftOver chain."
                    ),
                )
            )
        lines.extend(
            (
                "    The retained mapping is represented by a top-level net fill.",
                (
                    f"    All {report.source_interval.length}/"
                    f"{report.source_interval.length} input bases are present in the "
                    "reciprocal-best chain."
                ),
                (
                    f"    None of the other {all_chain_count - 1} complete "
                    f"{'mapping has' if all_chain_count - 1 == 1 else 'mappings have'} "
                    "that same combination."
                ),
            )
        )
        return lines

    if (
        profile.inventory_state
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    ):
        lines.append(
            "    The standard liftOver chain and all-chain alignments contain the "
            "same mappings."
        )
    else:
        additional = len(profile.additional_all_chain_candidate_ids)
        lines.extend(
            (
                (
                    f"    The all-chain alignments contain {additional} additional "
                    f"{'mapping' if additional == 1 else 'mappings'} not retained by"
                ),
                "    the standard liftOver chain.",
            )
        )

    if profile.state is ComparativeRelationshipState.NO_COMPETING_FULL_PLACEMENTS:
        lines.append(
            "    There are not two competing complete mappings for the comparative "
            "resources to distinguish."
        )
        return lines

    if profile.state is ComparativeRelationshipState.DOES_NOT_SEPARATE_PLACEMENTS:
        lines.append(
            "    The assessed UCSC chain/net relationships do not distinguish among "
            "the complete mappings."
        )
    elif profile.state is ComparativeRelationshipState.MIXED_CONFLICTING:
        lines.append(
            "    The standard liftOver chain, net, and reciprocal-best chain "
            "distinguish different complete mappings."
        )
    else:
        raise ValueError(
            f"unsupported comparative relationship state: {profile.state!r}"
        )

    lines.extend(_multiple_comparative_support_lines(report, profile))
    return lines


def _multiple_comparative_support_lines(
    report: UCSCAssessmentReport,
    profile: ComparativeRelationshipProfile,
) -> list[str]:
    complete = tuple(
        item for item in profile.placement_support if item.complete_source_coverage
    )
    return [
        "    Standard liftOver chain:",
        "        "
        + _summary_comparative_support_set_text(
            report, tuple(item for item in complete if item.retained_by_filtered_chain)
        ),
        "    Top-level net fill:",
        "        "
        + _summary_comparative_support_set_text(
            report, tuple(item for item in complete if item.depth1_top_net)
        ),
        "    Reciprocal-best chain:",
        "        "
        + _summary_comparative_support_set_text(
            report, tuple(item for item in complete if item.full_reciprocal_best)
        ),
    ]


def _summary_comparative_support_set_text(
    report: UCSCAssessmentReport,
    support: tuple[ComparativePlacementProfile, ...],
) -> str:
    if not support:
        return "none"
    if len(support) > _DEFAULT_INLINE_PROJECTION_LIMIT:
        return f"{len(support)} complete mappings; use --details for coordinates"
    return ", ".join(
        _human_interval_text(
            report.target_db,
            _report_candidate_for_id(report, item.candidate_id).target_interval,
        )
        for item in support
    )


def _report_candidate_for_id(
    report: UCSCAssessmentReport,
    candidate_id: str,
) -> NormalizedCandidate:
    matches = tuple(
        candidate
        for candidate in report.candidates
        if candidate.candidate_id == candidate_id
    )
    if len(matches) != 1:
        raise ValueError(
            "comparative placement support must identify exactly one report candidate"
        )
    return matches[0]


def _multiple_segmental_duplication_lines(
    report: UCSCAssessmentReport,
) -> list[str]:
    context = report.segmental_duplication_context_result
    if context is None:
        return []
    source_overlap = bool(context.source_overlaps)
    target_overlap = bool(context.target_overlaps)
    if not source_overlap and not target_overlap:
        return []

    lines = ["Segmental Duplications:"]
    if source_overlap:
        lines.append(
            "    The source interval overlaps the UCSC Segmental Duplications track."
        )
    if target_overlap:
        candidate_count = len({item.candidate_id for item in context.target_overlaps})
        noun = "mapping" if candidate_count == 1 else "mappings"
        lines.append(
            f"    {candidate_count} target {noun} overlap the UCSC Segmental "
            "Duplications track."
        )
    return lines


def _multiple_mapping_guidance_lines(
    report: UCSCAssessmentReport,
) -> tuple[list[str], list[str]]:
    """Explain why multiple mappings matter and what a user can do next."""

    profile = report.result_profile
    comparative = profile.comparative_relationship
    why: list[str] = []
    next_steps: list[str] = []

    complete_count = sum(
        candidate.covered_source_bases == candidate.source_bases
        for candidate in profile.candidate_profiles
    )
    if complete_count > 1:
        why.append(
            "More than one complete coordinate mapping exists, so liftOver alone does "
            "not identify which mapped locus, if any, represents the same biological "
            "feature."
        )
    elif (
        profile.union_covered_source_bases
        > profile.maximum_candidate_covered_source_bases
    ):
        why.append(
            "Different mappings cover different parts of the source interval. Their "
            "combined coverage is not one continuous mapping and should not be joined "
            "into a single target interval."
        )
    else:
        why.append(
            "Multiple liftOver mappings exist, so mapping order does not identify a "
            "preferred biological locus."
        )

    if comparative.state is ComparativeRelationshipState.FAVORS_ONE_PLACEMENT:
        why.append(
            "Within the assessed UCSC alignment resources, one mapping is distinguished "
            "because it is retained by the standard liftOver chain and also has "
            "top-level net and full reciprocal-best support."
        )
        why.append(
            "Those UCSC resources share alignment lineage and are not independent "
            "confirmations; the distinction does not establish biological correctness."
        )
    elif comparative.state is ComparativeRelationshipState.DOES_NOT_SEPARATE_PLACEMENTS:
        why.append(
            "The assessed UCSC alignment relationships do not distinguish among the "
            "complete mappings, so they do not provide a basis for choosing one by "
            "mapping order."
        )
        why.append(
            "The standard liftOver chain, net, and reciprocal-best chain are related "
            "UCSC alignment evidence rather than independent confirmations."
        )
    elif comparative.state is ComparativeRelationshipState.MIXED_CONFLICTING:
        why.append(
            "The assessed UCSC alignment relationships distinguish different complete "
            "mappings, so these resources do not provide one consistent mapping choice."
        )
        why.append(
            "The standard liftOver chain, net, and reciprocal-best chain are related "
            "UCSC alignment evidence rather than independent confirmations."
        )
    elif comparative.state is ComparativeRelationshipState.NO_COMPETING_FULL_PLACEMENTS:
        why.append(
            "The additional mappings are not competing complete mappings, so source "
            "coverage remains the important distinction among them."
        )

    reverse_state = profile.scope.reverse_result
    if reverse_state is ReverseCheckState.UNAVAILABLE:
        why.append(
            "Reverse liftOver was unavailable, so reciprocity was not assessed for "
            "these mappings."
        )

    duplication = report.segmental_duplication_context_result
    if duplication is not None and (
        duplication.source_overlaps or duplication.target_overlaps
    ):
        why.append(
            "Segmental-duplication overlap is relevant duplicated-sequence context, "
            "but it does not by itself identify the biologically corresponding mapping."
        )

    next_steps.append(
        "Use --details to inspect every mapping, its source coverage, alignment blocks, "
        "and the comparative evidence attached to it."
    )

    if comparative.state is ComparativeRelationshipState.FAVORS_ONE_PLACEMENT:
        favored_id = comparative.favored_candidate_id
        if favored_id is not None:
            favored = _report_candidate_for_id(report, favored_id)
            next_steps.extend(
                (
                    (
                        "Review the source and distinguished mapping in the "
                        "UCSC Genome Browser:"
                    ),
                    (
                        "    Source: "
                        + _ucsc_browser_url(report.source_db, report.source_interval)
                    ),
                    (
                        "    Distinguished mapping: "
                        + _ucsc_browser_url(report.target_db, favored.target_interval)
                    ),
                )
            )
    elif len(report.candidates) <= _DEFAULT_INLINE_PROJECTION_LIMIT:
        next_steps.extend(
            (
                "Review the source and mapped loci in the UCSC Genome Browser:",
                f"    Source: {_ucsc_browser_url(report.source_db, report.source_interval)}",
            )
        )
        next_steps.extend(
            f"    Mapping {index}: "
            + _ucsc_browser_url(report.target_db, candidate.target_interval)
            for index, candidate in enumerate(report.candidates, start=1)
        )

    if (
        comparative.state is ComparativeRelationshipState.NOT_ASSESSED
        and report.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
    ):
        next_steps.append(
            "If comparative resources are available, rerun with --evidence-tier "
            "COMPARATIVE to compare the standard liftOver chain with broader UCSC "
            "alignment relationships."
        )

    next_steps.append(
        "If the choice of locus affects a named variant, gene, transcript, or other "
        "biological feature, use target-assembly-specific feature evidence before "
        "selecting a mapping."
    )
    return why, next_steps


def _render_single_mapping_summary(report: UCSCAssessmentReport) -> str:
    """Render a compact standard-terminology summary for zero/one mappings."""

    profile = report.result_profile
    lines = [
        f"* {_summary_headline_text(report)} *",
        "",
        "Source:",
        f"    {_human_interval_text(report.source_db, profile.source_interval)}",
    ]

    if not profile.candidate_profiles:
        lines.extend(
            (
                "",
                "= KEY FINDINGS =",
                "",
                "Mapping:",
                "    No liftOver mapping was found for the requested source interval.",
            )
        )
        context_lines = _single_flanking_interval_lines(report)
        if context_lines:
            lines.extend(("", *context_lines))
        duplication_lines = _single_segmental_duplication_lines(report)
        if duplication_lines:
            lines.extend(("", *duplication_lines))
        why_lines, next_step_lines = _no_mapping_guidance_lines(report)
        lines.extend(
            ("", "= WHY THIS MATTERS =", "", *_wrap_summary_paragraphs(why_lines))
        )
        lines.extend(
            ("", "= NEXT STEP =", "", *_wrap_summary_paragraphs(next_step_lines))
        )
        evidence_scope_lines = _single_evidence_scope_lines(report)
        if evidence_scope_lines:
            lines.extend(("", *evidence_scope_lines))
        lines.extend(
            (
                "",
                "Details:",
                "    Use --details for alignment evidence, resources, and provenance;",
                "    use --json for machine-readable output.",
            )
        )
        return "\n".join(lines)

    candidate = report.candidates[0]
    candidate_profile = profile.candidate_profiles[0]
    target_label = _single_mapping_target_label(profile, candidate_profile)
    lines.extend(
        (
            "",
            f"{target_label}:",
            f"    {_human_interval_text(report.target_db, candidate.target_interval)}",
        )
    )
    if candidate.orientation.value == "REVERSE":
        lines.extend(("    Orientation:", "        reverse"))
    lines.extend(
        (
            "",
            "= KEY FINDINGS =",
            "",
        )
    )
    lines.extend(_single_mapping_finding_lines(report, candidate_profile))

    why_lines, next_step_lines = _single_mapping_guidance_lines(
        report, candidate_profile
    )
    if why_lines:
        lines.extend(
            ("", "= WHY THIS MATTERS =", "", *_wrap_summary_paragraphs(why_lines))
        )
    if next_step_lines:
        lines.extend(
            ("", "= NEXT STEP =", "", *_wrap_summary_paragraphs(next_step_lines))
        )
    if not why_lines and not next_step_lines:
        limitation_lines = _single_mapping_limitation_lines(report)
        if limitation_lines:
            lines.extend(("", "= LIMITATIONS =", "", *limitation_lines))

    evidence_scope_lines = _single_evidence_scope_lines(report)
    if evidence_scope_lines:
        lines.extend(("", *evidence_scope_lines))

    lines.extend(
        (
            "",
            "Details:",
            "    Use --details for alignment blocks, gaps, evidence, and provenance;",
            "    use --json for machine-readable output.",
        )
    )
    return "\n".join(lines)


def _single_evidence_scope_lines(report: UCSCAssessmentReport) -> list[str]:
    lines: list[str] = []
    consumed_roles = set(report.result_profile.consumed_resource_roles)

    if report.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        resources: list[str] = []
        if "CHAIN" in consumed_roles:
            resources.append("UCSC all-chain alignments")
        if "NET" in consumed_roles:
            resources.append("UCSC net alignment")
        if "RECIPROCAL_BEST_CHAIN" in consumed_roles:
            resources.append("UCSC reciprocal-best chain")
        if resources:
            lines.extend(("Evidence:", "    " + ", ".join(resources) + "."))

    if report.result_profile.scope.target_role is TargetRoleState.UNAVAILABLE:
        lines.extend(
            (
                "Target sequence role:",
                (
                    "    Not assessed because version-matched NCBI assembly sequence "
                    "metadata was unavailable."
                ),
            )
        )
    return lines


def _ucsc_browser_url(database: str, interval: GenomicInterval) -> str:
    """Return a direct UCSC Genome Browser link for one assembly interval."""

    position = f"{interval.sequence_name}:{interval.start + 1}-{interval.end}"
    query = urlencode({"db": database, "position": position})
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?{query}"


def _no_mapping_guidance_lines(
    report: UCSCAssessmentReport,
) -> tuple[list[str], list[str]]:
    """Explain a valid source interval with no mapping in the consumed chain."""

    why = [
        (
            "No liftOver mapping means the consumed chain does not provide a coordinate "
            "conversion for this source interval; it does not establish biological "
            "deletion or absence of homologous sequence from the target assembly."
        )
    ]
    next_steps = [
        "Review the source locus in the UCSC Genome Browser:",
        f"    {_ucsc_browser_url(report.source_db, report.source_interval)}",
    ]
    if report.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        next_steps.append(
            "If comparative resources are available, rerun with --evidence-tier "
            "COMPARATIVE to look for additional UCSC chain alignments."
        )
    next_steps.append(
        "If a named variant, gene, or transcript is the real target of the conversion, "
        "check that feature directly in a target-assembly-specific source."
    )
    return why, next_steps


def _single_browser_review_lines(report: UCSCAssessmentReport) -> list[str]:
    """Render direct Browser navigation for an unusual single-mapping result."""

    if len(report.candidates) != 1:
        return []
    return [
        "Review the source and mapped loci in the UCSC Genome Browser:",
        f"    Source: {_ucsc_browser_url(report.source_db, report.source_interval)}",
        (
            "    Mapped: "
            + _ucsc_browser_url(report.target_db, report.candidates[0].target_interval)
        ),
    ]


def _single_mapping_guidance_lines(
    report: UCSCAssessmentReport,
    profile: CandidateResultProfile,
) -> tuple[list[str], list[str]]:
    """Explain why unusual single-mapping observations matter and what to do next."""

    why: list[str] = []
    next_steps: list[str] = []
    browser_review = False
    local_why_present = False

    if profile.covered_source_bases < profile.source_bases:
        why.append(
            "Not all requested source bases map, so a feature spanning the unmapped "
            "bases cannot be transferred as one complete interval."
        )
        next_steps.append(
            "Use --details to inspect the exact alignment blocks, gaps, and unmapped "
            "source intervals before transferring a larger feature."
        )
        browser_review = True

    if profile.geometric_segment_count > 1 or profile.target_discontinuous:
        why.append(
            "The mapped target span is not one continuous alignment; treating the "
            "bounding span as continuous would include unaligned target sequence."
        )
        if not any("alignment blocks" in line for line in next_steps):
            next_steps.append(
                "Use --details to inspect the exact alignment blocks and gaps before "
                "using the target span as one interval."
            )
        browser_review = True

    if profile.orientation is MappingOrientation.REVERSE:
        why.append(
            "Reverse orientation means the source and target align on opposite "
            "strands; that strand relationship is not by itself a mapping error."
        )
        next_steps.append(
            "Account for strand when interpreting strand-sensitive features, alleles, "
            "or transcript structure at the mapped locus."
        )
        browser_review = True

    reverse = profile.reverse_mapping
    if reverse.check_state is ReverseCheckState.UNAVAILABLE:
        why.append(
            "Reverse liftOver was not available, so this run does not establish "
            "whether the mapping is reciprocal."
        )
        next_steps.append(
            "Prepare the matching reverse liftOver chain index/resources and rerun "
            "if reciprocity matters for the intended use."
        )
    elif reverse.check_state is ReverseCheckState.RUN:
        if (
            reverse.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY
            and not reverse.exact_original_geometry_return
        ):
            why.append(
                "Reverse liftOver returns only to the source locus but does not "
                "reconstruct the complete original aligned geometry."
            )
            next_steps.append(
                "Use --details to inspect the recovered source coverage before "
                "treating the mapping as fully reciprocal."
            )
            browser_review = True
        elif reverse.relationship in {
            ReverseRelationshipState.ELSEWHERE_ONLY,
            ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE,
            ReverseRelationshipState.NO_PROJECTION,
        }:
            local_why_present = True
            browser_review = True

    context_profile = report.result_profile.query_context
    context_findings = set(context_profile.findings)
    context_disagrees = bool(
        context_findings
        & {
            QueryContextFinding.REVEALS_PARTIAL_COVERAGE,
            QueryContextFinding.REVEALS_FRAGMENTATION,
            QueryContextFinding.REVEALS_TARGET_DISCONTINUITY,
            QueryContextFinding.CHANGES_WITH_QUERY_SCALE,
        }
    )
    if context_disagrees:
        why.append(
            "The point and its flanking interval do not show the same alignment "
            "behavior, so the local sequence context is more complex than the point "
            "alone suggests."
        )
        next_steps.append(
            "Use --details to inspect the flanking interval geometry; request a wider "
            "--context-bases window only when a larger biological context is justified."
        )
        browser_review = True

    duplication = report.segmental_duplication_context_result
    if duplication is not None and (
        duplication.source_overlaps or duplication.target_overlaps
    ):
        local_why_present = True
        browser_review = True

    target_roles = report.result_profile.target_sequence_roles
    if any(item.context is None for item in target_roles):
        why.append(
            "The version-matched assembly metadata did not identify the target "
            "sequence role, so primary/alternate/unplaced status is unresolved."
        )
        next_steps.append(
            "Check the target assembly sequence metadata before applying a workflow "
            "that filters by primary, alternate, or unplaced sequence role."
        )
        browser_review = True
    elif any(
        item.context is not None
        and (
            item.context.provider_role != "assembled-molecule"
            or item.context.assembly_unit != "Primary Assembly"
        )
        for item in target_roles
    ):
        why.append(
            "The mapped sequence is not an assembled molecule on the Primary Assembly in "
            "the version-matched assembly metadata, which can matter for downstream "
            "tools that restrict analyses to primary-assembly sequences."
        )
        next_steps.append(
            "Confirm that the downstream workflow accepts this target sequence role "
            "before filtering or discarding the mapping."
        )
        browser_review = True

    comparative = report.result_profile.comparative_relationship
    if (
        comparative.state is not ComparativeRelationshipState.NOT_ASSESSED
        and comparative.inventory_state
        is FilteredAllChainInventoryState.ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS
    ):
        why.append(
            "The all-chain alignments contain additional alignment relationships that "
            "are not retained by the standard liftOver chain."
        )
        next_steps.append(
            "Use --details to inspect the additional all-chain mappings and their "
            "coverage before treating the standard liftOver result as unique."
        )

    if not why and not local_why_present:
        return [], []

    if browser_review:
        next_steps[:0] = _single_browser_review_lines(report)

    if (
        report.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
        and reverse.check_state is ReverseCheckState.RUN
        and reverse.relationship
        in {
            ReverseRelationshipState.ELSEWHERE_ONLY,
            ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE,
            ReverseRelationshipState.NO_PROJECTION,
        }
    ):
        next_steps.append(
            "If comparative resources are available, rerun with --evidence-tier "
            "COMPARATIVE to look for additional UCSC chain alignments."
        )

    next_steps.append(
        "If this coordinate will stand in for a named variant, gene, or transcript, "
        "verify that identity separately in a target-assembly-specific source."
    )
    return why, next_steps


def _single_mapping_target_label(
    profile: ResultProfile,
    candidate: CandidateResultProfile,
) -> str:
    if profile.source_interval.length == 1:
        return "Mapped coordinate"
    if candidate.geometric_segment_count > 1:
        return "Target span"
    return "Mapped interval"


def _single_mapping_finding_lines(
    report: UCSCAssessmentReport,
    profile: CandidateResultProfile,
) -> list[str]:
    source_bases = profile.source_bases
    lines = [
        "Mapping:",
        (
            f"    One liftOver chain maps {profile.covered_source_bases}/"
            f"{source_bases} bp"
            + (
                "; no chain gap at the query."
                if source_bases == 1 and not profile.source_gap_intervals
                else "."
            )
        ),
    ]

    if profile.geometric_segment_count > 1:
        lines.extend(
            (
                (
                    "    The mapping contains "
                    f"{profile.geometric_segment_count} alignment blocks."
                ),
                (
                    "    The target span above is a bounding span, not one "
                    "continuous alignment."
                ),
            )
        )
    if profile.uncovered_source_intervals:
        lines.extend(
            (
                "    Unmapped source interval(s):",
                *(
                    "        " + _human_interval_text(report.source_db, interval)
                    for interval in profile.uncovered_source_intervals
                ),
            )
        )
    if profile.target_gap_intervals:
        lines.extend(
            (
                "    Target gap(s) within the mapping:",
                *(
                    "        " + _human_interval_text(report.target_db, interval)
                    for interval in profile.target_gap_intervals
                ),
            )
        )

    reverse_lines = _single_reverse_liftover_lines(report, profile)
    if reverse_lines:
        lines.extend(("", *reverse_lines))

    context_lines = _single_flanking_interval_lines(report)
    if context_lines:
        lines.extend(("", *context_lines))

    comparative_lines = _single_comparative_lines(report)
    if comparative_lines:
        lines.extend(("", *comparative_lines))

    duplication_lines = _single_segmental_duplication_lines(report)
    if duplication_lines:
        lines.extend(("", *duplication_lines))

    target_role_lines = _single_target_sequence_metadata_lines(report)
    if target_role_lines:
        lines.extend(("", *target_role_lines))
    return lines


def _single_reverse_liftover_lines(
    report: UCSCAssessmentReport,
    profile: CandidateResultProfile,
) -> list[str]:
    reverse = profile.reverse_mapping
    if reverse.check_state is ReverseCheckState.NOT_RUN:
        return []
    if reverse.check_state is ReverseCheckState.UNAVAILABLE:
        return [
            "Reverse liftOver:",
            "    Reverse liftOver was unavailable from the prepared resources.",
        ]

    queried = reverse.queried_target_segments
    if len(queried) == 1:
        query_text = _human_interval_text(report.target_db, queried[0])
    else:
        query_text = f"{len(queried)} mapped target blocks"

    if reverse.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_ONLY:
        if reverse.exact_original_geometry_return:
            return [
                "Reverse liftOver:",
                "    Returns exactly to the source locus.",
            ]
        assert reverse.original_source_covered_bases is not None
        return [
            "Reverse liftOver:",
            "    Returns only to the source locus.",
            (
                "    Recovered original source bases: "
                f"{reverse.original_source_covered_bases}/"
                f"{reverse.original_source_bases}."
            ),
        ]

    if reverse.relationship is ReverseRelationshipState.NO_PROJECTION:
        lines = [
            "Reverse liftOver:",
            (
                f"    {query_text} has no mapping in the "
                f"{report.target_db}→{report.source_db} chain."
            ),
        ]
        lines.extend(
            _local_why_lines(
                [
                    (
                        "Because the mapped target does not map back through the available "
                        "reverse chain, this check does not demonstrate a one-to-one "
                        "reciprocal coordinate relationship."
                    ),
                    (
                        "The forward mapping is not necessarily wrong, but coordinate "
                        "conversion alone is not sufficient evidence of biological feature "
                        "identity in this case."
                    ),
                ]
            )
        )
        return lines

    mapped_back = _reverse_target_intervals(report)
    lines = ["Reverse liftOver:"]
    if reverse.relationship is ReverseRelationshipState.ELSEWHERE_ONLY:
        lines.append("    Does not return to the source locus.")
    elif reverse.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE:
        lines.append("    Returns to the source locus and to at least one other locus.")

    if mapped_back:
        if len(mapped_back) == 1:
            mapped_text = _human_interval_text(report.source_db, mapped_back[0])
            if reverse.relationship is ReverseRelationshipState.ELSEWHERE_ONLY:
                lines[-1] = (
                    "    Does not return to the source locus; "
                    f"maps instead to {mapped_text}."
                )
            else:
                lines.append(f"    Also maps to {mapped_text}.")
        else:
            prefix = (
                "    Maps instead to"
                if reverse.relationship is ReverseRelationshipState.ELSEWHERE_ONLY
                else "    Also maps to"
            )
            lines.append(f"{prefix} {len(mapped_back)} source-assembly locations:")
            lines.extend(
                "        " + _human_interval_text(report.source_db, interval)
                for interval in mapped_back[:_DEFAULT_INLINE_PROJECTION_LIMIT]
            )

    if reverse.relationship is ReverseRelationshipState.ELSEWHERE_ONLY:
        lines.extend(
            _local_why_lines(
                [
                    (
                        "The mapped coordinate does not return to the original locus when "
                        "lifted back. The source and mapped loci therefore do not have a "
                        "one-to-one reciprocal correspondence under the available forward "
                        "and reverse chain mappings."
                    ),
                    (
                        "The forward mapping is not necessarily wrong, but it should not by "
                        "itself be treated as sufficient evidence that the target coordinate "
                        "represents the same variant, gene, transcript, or other biological "
                        "feature."
                    ),
                ]
            )
        )
    elif reverse.relationship is ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE:
        lines.extend(
            _local_why_lines(
                [
                    (
                        "Reverse liftOver returns to the source locus and to other loci. "
                        "That one-to-many reverse relationship means the source and mapped "
                        "loci do not have a one-to-one reciprocal correspondence under the "
                        "available chain mappings."
                    ),
                    (
                        "The forward mapping is not necessarily wrong, but it should not by "
                        "itself be treated as sufficient evidence that the target coordinate "
                        "represents the same biological feature."
                    ),
                ]
            )
        )
    return lines


def _reverse_target_intervals(
    report: UCSCAssessmentReport,
) -> tuple[GenomicInterval, ...]:
    if not report.reverse_mapping_results:
        return ()
    intervals: list[GenomicInterval] = []
    for segment_result in report.reverse_mapping_results[0].segment_results:
        intervals.extend(
            candidate.target_interval for candidate in segment_result.candidates
        )
    return tuple(dict.fromkeys(intervals))


def _single_flanking_interval_lines(report: UCSCAssessmentReport) -> list[str]:
    result = report.query_context_result
    profile = report.result_profile.query_context
    if result is None or profile.check_state is QueryContextState.NOT_RUN:
        return []
    tested = result.tested_source_interval
    if tested is None:
        raise ValueError("completed flanking-interval check requires a source interval")

    lines = ["Flanking interval:"]
    if profile.point_and_local_context_map_together and len(result.candidates) == 1:
        context_candidate = result.candidates[0]
        covered = profile.maximum_candidate_covered_source_bases
        assert covered is not None
        lines.extend(
            (
                (
                    "    "
                    + _human_interval_text(report.source_db, tested)
                    + " → "
                    + _human_interval_text(
                        report.target_db, context_candidate.target_interval
                    )
                ),
                f"    {covered}/{tested.length} bp map through the same liftOver chain.",
            )
        )
        return lines

    lines.append(
        f"    A {tested.length}-bp interval centered on the input coordinate was "
        "also assessed."
    )
    findings = set(profile.findings)
    if QueryContextFinding.REVEALS_PARTIAL_COVERAGE in findings:
        lines.append("    The flanking interval has partial source coverage.")
    if QueryContextFinding.REVEALS_FRAGMENTATION in findings:
        lines.append("    The flanking interval maps in multiple alignment blocks.")
    if QueryContextFinding.REVEALS_TARGET_DISCONTINUITY in findings:
        lines.append("    The flanking interval contains a target-side alignment gap.")
    if QueryContextFinding.CHANGES_WITH_QUERY_SCALE in findings:
        lines.append(
            "    The liftOver result changes when the larger interval is assessed."
        )
    if QueryContextFinding.NO_PROJECTION_AT_EITHER_SCALE in findings:
        lines.append(
            "    No liftOver mapping was found for the point or the flanking interval."
        )
    return lines


def _single_comparative_lines(report: UCSCAssessmentReport) -> list[str]:
    profile = report.result_profile.comparative_relationship
    if profile.state is ComparativeRelationshipState.NOT_ASSESSED:
        return []
    comparison = report.filtered_all_chain_comparison
    if comparison is None:
        raise ValueError("comparative summary requires filtered/all-chain comparison")

    lines = ["Comparative UCSC evidence:"]
    if (
        profile.inventory_state
        is FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE
    ):
        lines.append(
            "    The standard liftOver chain and the UCSC all-chain alignments "
            "contain the same mapping."
        )
    else:
        additional = len(profile.additional_all_chain_candidate_ids)
        lines.append(
            "    The UCSC all-chain alignments contain "
            f"{additional} additional {'mapping' if additional == 1 else 'mappings'} "
            "not retained by the standard liftOver chain."
        )

    if profile.state is ComparativeRelationshipState.NO_COMPETING_FULL_PLACEMENTS:
        lines.append(
            "    No additional complete mapping is present in the all-chain alignments."
        )
    return lines


def _single_segmental_duplication_lines(report: UCSCAssessmentReport) -> list[str]:
    context = report.segmental_duplication_context_result
    if context is None:
        return []
    source_overlap = bool(context.source_overlaps)
    target_overlap = bool(context.target_overlaps)
    if not source_overlap and not target_overlap:
        return []
    if source_overlap and target_overlap:
        text = (
            "Both the source and mapped coordinates overlap the UCSC "
            "Segmental Duplications track."
        )
    elif source_overlap:
        text = "The source coordinate overlaps the UCSC Segmental Duplications track."
    else:
        text = "The mapped coordinate overlaps the UCSC Segmental Duplications track."

    lines = ["Segmental Duplications:", f"    {text}"]
    profile = report.result_profile
    if len(profile.candidate_profiles) == 1 and context.source_overlaps:
        source_sequence = profile.source_interval.sequence_name
        target_sequence = report.candidates[0].target_interval.sequence_name
        if source_sequence != target_sequence and any(
            overlap.record.paired_interval.sequence_name == target_sequence
            for overlap in context.source_overlaps
        ):
            lines.append(
                f"    One overlapping {report.source_db} Segmental Duplications record "
                f"pairs the {source_sequence} source region with {report.source_db} "
                f"{target_sequence}."
            )
    lines.extend(
        _local_why_lines(
            [
                (
                    "Duplicated sequence can complicate interpretation of mappings between "
                    "loci, so this context is reported when present. The overlap does not "
                    "by itself establish paralogy, mapping error, non-uniqueness, or the "
                    "cause of this mapping."
                )
            ]
        )
    )
    return lines


def _single_target_sequence_metadata_lines(report: UCSCAssessmentReport) -> list[str]:
    profile = report.result_profile
    state = profile.scope.target_role
    if state is not TargetRoleState.ASSESSED:
        return []
    unusual = [
        item
        for item in profile.target_sequence_roles
        if item.context is None
        or item.context.provider_role != "assembled-molecule"
        or item.context.assembly_unit != "Primary Assembly"
    ]
    if not unusual:
        return []
    lines = ["Target assembly sequence:"]
    for item in unusual[:_DEFAULT_INLINE_PROJECTION_LIMIT]:
        if item.context is None:
            lines.append(
                f"    {item.sequence_name}: no matching row in the version-matched "
                "NCBI sequence report."
            )
            continue
        lines.extend(
            (
                f"    {item.sequence_name}",
                f"        Assembly unit: {item.context.assembly_unit}",
                f"        Sequence role: {item.context.provider_role}",
            )
        )
    return lines


def _single_mapping_limitation_lines(report: UCSCAssessmentReport) -> list[str]:
    mapped_object = (
        "mapped coordinate"
        if report.result_profile.source_interval.length == 1
        else "mapped interval"
    )
    lines = [
        (
            f"This result does not establish that the {mapped_object} is unique or "
            "represents the same variant, gene, transcript, or other feature."
        ),
    ]
    context = report.segmental_duplication_context_result
    if context is not None and (context.source_overlaps or context.target_overlaps):
        candidate_profiles = report.result_profile.candidate_profiles
        reverse_returns_elsewhere = bool(
            len(candidate_profiles) == 1
            and candidate_profiles[0].reverse_mapping.relationship
            in {
                ReverseRelationshipState.ELSEWHERE_ONLY,
                ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE,
            }
        )
        lines.append(
            "Overlap with UCSC Segmental Duplications does not by itself establish "
            "paralogy, mapping error, or non-uniqueness."
        )
        if reverse_returns_elsewhere:
            lines.append("It does not explain the non-reciprocal mapping.")
    return lines


def _summary_check_lines(
    report: UCSCAssessmentReport,
    *,
    include_forward: bool = False,
) -> list[str]:
    lines: list[str] = []
    consumed_roles = set(report.result_profile.consumed_resource_roles)

    if include_forward:
        if report.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
            lines.append(f"    {report.source_db} → {report.target_db} liftOver")
        elif "CHAIN" in consumed_roles:
            lines.append("    UCSC all-chain alignments")

    reverse_state = report.result_profile.scope.reverse_result
    if reverse_state is ReverseCheckState.RUN:
        lines.append(f"    {report.target_db} → {report.source_db} reverse liftOver")
    elif reverse_state is ReverseCheckState.UNAVAILABLE:
        lines.append("    Reverse liftOver (resource unavailable)")

    context = report.result_profile.query_context
    if (
        context.check_state is QueryContextState.RUN
        and context.actual_window_bases is not None
    ):
        lines.append(
            f"    {context.actual_window_bases}-bp flanking interval centered on "
            "the input coordinate"
        )

    target_role = report.result_profile.scope.target_role
    if target_role is TargetRoleState.ASSESSED:
        lines.append("    NCBI assembly sequence metadata")
    elif target_role is TargetRoleState.UNAVAILABLE:
        lines.append("    NCBI assembly sequence metadata (unavailable)")

    if report.evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
        if not include_forward and "CHAIN" in consumed_roles:
            lines.append("    UCSC all-chain alignments")
        if "NET" in consumed_roles:
            lines.append("    UCSC net alignment")
        if "RECIPROCAL_BEST_CHAIN" in consumed_roles:
            lines.append("    UCSC reciprocal-best chain")
        if report.filtered_all_chain_comparison is not None:
            lines.append(
                "    Standard liftOver chain compared with UCSC all-chain alignments"
            )

    external = report.result_profile.scope.external_context
    if external is ExternalContextState.ASSESSED:
        lines.append("    UCSC Segmental Duplications track")
    elif external is ExternalContextState.PARTIALLY_ASSESSED:
        lines.append("    UCSC Segmental Duplications track (partially available)")
    elif external is ExternalContextState.UNAVAILABLE:
        lines.append("    UCSC Segmental Duplications track (unavailable)")
    return lines


def _summary_headline_text(report: UCSCAssessmentReport) -> str:
    profile = report.result_profile
    if len(profile.candidate_profiles) == 1:
        candidate = report.candidates[0]
        source_sequence = profile.source_interval.sequence_name
        target_sequence = candidate.target_interval.sequence_name
        if (
            _is_standard_ucsc_chromosome_name(source_sequence)
            and _is_standard_ucsc_chromosome_name(target_sequence)
            and source_sequence != target_sequence
        ):
            return "INTERCHROMOSOMAL LIFTOVER MAPPING"

    return _headline_text(profile.headline)


def _is_standard_ucsc_chromosome_name(sequence_name: str) -> bool:
    """Return whether a UCSC sequence label plainly names a chromosome."""

    return re.fullmatch(r"chr(?:[1-9][0-9]*|X|Y)", sequence_name) is not None


def _human_interval_text(database: str, interval: GenomicInterval) -> str:
    if interval.length == 1:
        coordinate = f"{interval.sequence_name}:{interval.start + 1}"
    else:
        coordinate = f"{interval.sequence_name}:{interval.start + 1}-{interval.end}"
    return f"{database} {coordinate}"


def _comparative_candidate_label(
    report: UCSCAssessmentReport,
    candidate_id: str,
) -> str:
    candidates = tuple(
        candidate
        for candidate in report.candidates
        if candidate.candidate_id == candidate_id
    )
    if len(candidates) != 1:
        raise ValueError(
            "comparative placement support must identify exactly one report candidate"
        )
    candidate = candidates[0]
    interval = candidate.target_interval
    return (
        f"{interval.sequence_name}:{interval.start + 1}-{interval.end} "
        f"({candidate.orientation.value.lower()} orientation; {candidate_id})"
    )


def _comparative_detail_lines(report: UCSCAssessmentReport) -> list[str]:
    comparison = report.filtered_all_chain_comparison
    relationship = report.comparative_evidence_relationship
    profile = report.result_profile.comparative_relationship
    if comparison is None or relationship is None:
        raise ValueError(
            "comparative details require paired inventory and relationship"
        )

    lines = [
        f"  Mapping inventory: {_detail_inventory_state_text(comparison.relationship)}",
        f"  All-chain mappings: {len(comparison.all_chain_candidates)}",
        f"  Standard liftOver mappings: {len(comparison.filtered_candidates)}",
        (
            "  Additional all-chain mappings: "
            + (
                ", ".join(comparison.additional_all_chain_candidate_ids)
                if comparison.additional_all_chain_candidate_ids
                else "none"
            )
        ),
        f"  Comparative result: {_detail_comparative_relationship_text(profile.state)}",
        (
            "  Favored mapping: "
            + (
                relationship.favored_candidate_id
                if relationship.favored_candidate_id is not None
                else "none"
            )
        ),
        (
            "  UCSC pair dependency group: "
            + ", ".join(
                parent.source_id
                for parent in comparison.all_chain_provenance.derived_from
            )
        ),
        "  Exact shared processing-run provenance: not verified",
        (
            "  Dependency note: filtered chain, net, and reciprocal-best observations "
            "are conservatively grouped as dependent UCSC-derived evidence, not "
            "independent votes; the pair group does not prove that the files came from one "
            "processing run."
        ),
        "  Mapping support (all-chain order is reproducibility only, not rank):",
    ]
    if not profile.placement_support:
        lines.append("    none")
        return lines
    for item in profile.placement_support:
        lines.append(
            "    - "
            + _comparative_candidate_label(report, item.candidate_id)
            + "; complete source coverage="
            + _yes_no(item.complete_source_coverage)
            + "; retained by filtered chain="
            + _yes_no(item.retained_by_filtered_chain)
            + "; depth-1 top-net="
            + _yes_no(item.depth1_top_net)
            + "; full reciprocal-best="
            + _yes_no(item.full_reciprocal_best)
        )
    return lines


def _yes_no(value: bool) -> str:
    return "yes" if value else "no"


def _headline_text(headline: FactualHeadline) -> str:
    replacements = {
        FactualHeadline.NO_CHAIN_PROJECTION: "NO LIFTOVER MAPPING",
        FactualHeadline.ONE_COMPLETE_CHAIN_PROJECTION: "ONE LIFTOVER MAPPING",
        FactualHeadline.PARTIAL_SOURCE_COVERAGE: "PARTIAL LIFTOVER MAPPING",
        FactualHeadline.PARTIAL_AND_FRAGMENTED_PROJECTION: (
            "PARTIAL MAPPING ACROSS MULTIPLE ALIGNMENT BLOCKS"
        ),
        FactualHeadline.COMPLETE_BUT_DISCONTINUOUS_PROJECTION: (
            "LIFTOVER MAPPING WITH A TARGET GAP"
        ),
        FactualHeadline.MULTIPLE_CHAIN_PROJECTIONS: "MULTIPLE LIFTOVER MAPPINGS",
        FactualHeadline.SOURCE_INTERVAL_SPLITS_ACROSS_MULTIPLE_PROJECTIONS: (
            "SOURCE INTERVAL MAPS TO MULTIPLE LOCATIONS"
        ),
    }
    return replacements[headline]


def _candidate_text(
    candidate: NormalizedCandidate,
    profile: CandidateResultProfile,
) -> str:
    interval = candidate.target_interval
    coordinate_text = f"{interval.sequence_name}:{interval.start + 1}-{interval.end}"
    details = [
        "1-based inclusive",
        f"{candidate.orientation.value.lower()} orientation",
    ]
    if profile.geometric_segment_count > 1:
        details.append(
            f"bounding span of {profile.geometric_segment_count} mapped segments"
        )
    return f"{coordinate_text} ({'; '.join(details)})"


def _detail_state_text(value: str) -> str:
    replacements = {
        "VALID": "valid",
        "NOT_ASSESSED": "not assessed",
        "UNAVAILABLE": "unavailable",
        "NO_TARGET_PROJECTIONS": "no target mappings",
        "ASSESSED": "assessed",
        "PARTIALLY_ASSESSED": "partially assessed",
        "NOT_CHECKED": "not checked",
        "NOT_RUN": "not performed",
        "RUN": "performed",
        "NONE": "none",
        "COMPLETE": "complete",
        "FULL": "complete",
        "PARTIAL": "partial",
        "SAME": "same",
        "REVERSE": "reverse",
        "MIXED": "mixed",
    }
    return replacements.get(value, value.lower().replace("_", " "))


def _detail_evidence_tier_text(tier: EvidenceAvailabilityTier) -> str:
    if tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
        return "standard liftOver chain only"
    return "comparative UCSC alignments"


def _detail_resource_role_text(
    role: str,
    *,
    evidence_tier: EvidenceAvailabilityTier | None = None,
) -> str:
    if role == "CHAIN":
        if evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY:
            return "standard liftOver chain"
        if evidence_tier is EvidenceAvailabilityTier.COMPARATIVE:
            return "all-chain alignments"
    replacements = {
        "CHAIN": "chain alignment",
        "NET": "net alignment",
        "SYNTENIC_NET": "syntenic net alignment",
        "RECIPROCAL_BEST_CHAIN": "reciprocal-best chain",
        "RECIPROCAL_BEST_NET": "reciprocal-best net",
    }
    return replacements.get(role, role.lower().replace("_", " "))


def _capitalize_first(text: str) -> str:
    return text[:1].upper() + text[1:]


def _detail_reverse_relationship_text(relationship: ReverseRelationshipState) -> str:
    replacements = {
        ReverseRelationshipState.NO_PROJECTION: "no reverse mapping",
        ReverseRelationshipState.ORIGINAL_SOURCE_ONLY: (
            "returns only to the source locus"
        ),
        ReverseRelationshipState.ELSEWHERE_ONLY: (
            "does not return to the source locus"
        ),
        ReverseRelationshipState.ORIGINAL_SOURCE_AND_ELSEWHERE: (
            "returns to the source locus and elsewhere"
        ),
    }
    return replacements[relationship]


def _detail_query_context_finding_text(finding: QueryContextFinding) -> str:
    replacements = {
        QueryContextFinding.AGREES_WITH_POINT: "agrees with point mapping",
        QueryContextFinding.NO_PROJECTION_AT_EITHER_SCALE: (
            "no mapping at either scale"
        ),
        QueryContextFinding.REVEALS_PARTIAL_COVERAGE: "reveals partial coverage",
        QueryContextFinding.REVEALS_FRAGMENTATION: (
            "reveals multiple alignment blocks"
        ),
        QueryContextFinding.REVEALS_TARGET_DISCONTINUITY: ("reveals a target gap"),
        QueryContextFinding.CHANGES_WITH_QUERY_SCALE: "changes with query scale",
    }
    return replacements[finding]


def _detail_evidence_kind_text(kind: str) -> str:
    replacements = {
        "MAPPING_COVERAGE": "Source coverage",
        "CHAIN_GAPS": "Chain gaps",
        "CANDIDATE_RANK": "Mapping rank",
        "TARGET_PLACEMENT": "Target placement",
        "CHAIN_SCORE": "Chain score",
        "ALIGNED_BASES": "Aligned bases",
        "DUPLICATED_QUERY_BASES": "Duplicated query bases",
        "NET_CLASSIFICATION": "Net classification",
        "NET_HIERARCHY": "Net hierarchy",
        "RECIPROCAL_BEST_MEMBERSHIP": "Reciprocal-best chain coverage",
        "FLANKING_GENE_SYNTENY": "Flanking-gene synteny",
    }
    return replacements.get(kind, kind.lower().replace("_", " ").capitalize())


def _detail_inventory_state_text(state: FilteredAllChainInventoryState) -> str:
    replacements = {
        FilteredAllChainInventoryState.FILTERED_AND_ALL_CHAIN_AGREE: (
            "standard liftOver and all-chain mappings agree"
        ),
        FilteredAllChainInventoryState.ALL_CHAIN_REVEALS_ADDITIONAL_PLACEMENTS: (
            "all-chain alignments contain additional mappings"
        ),
    }
    return replacements[state]


def _detail_comparative_relationship_text(
    state: ComparativeRelationshipState,
) -> str:
    replacements = {
        ComparativeRelationshipState.NOT_ASSESSED: "not assessed",
        ComparativeRelationshipState.NO_COMPETING_FULL_PLACEMENTS: (
            "no competing complete mappings"
        ),
        ComparativeRelationshipState.FAVORS_ONE_PLACEMENT: "favors one mapping",
        ComparativeRelationshipState.DOES_NOT_SEPARATE_PLACEMENTS: (
            "does not separate mappings"
        ),
        ComparativeRelationshipState.MIXED_CONFLICTING: "mixed or conflicting",
    }
    return replacements[state]


def render_assessment_details(report: UCSCAssessmentReport) -> str:
    """Render the complete factual profile, evidence, resources, and provenance."""

    profile = report.result_profile
    lines = [
        "Detailed liftOver assessment",
        f"UCSC database pair: {report.source_db} -> {report.target_db}",
        f"Source locus: {format_display_interval(report.source_interval)}",
        f"Headline: {_headline_text(profile.headline)}",
        f"Interpretation: {profile.interpretation}",
        f"Source validation: {_detail_state_text(profile.input_validity.value)}",
        f"Mapping count: {len(report.candidates)}",
        f"Mapping orientation: {_detail_state_text(profile.orientation.value)}",
        (
            "Maximum source coverage for one mapping: "
            f"{profile.maximum_candidate_covered_source_bases}/{profile.source_bases}"
        ),
        (
            "Union source coverage across mappings: "
            f"{profile.union_covered_source_bases}/{profile.source_bases}"
        ),
        f"Evidence resources: {_detail_evidence_tier_text(profile.evidence_tier)}",
        "Consumed UCSC resources: "
        + (
            ", ".join(
                _detail_resource_role_text(role, evidence_tier=profile.evidence_tier)
                for role in profile.consumed_resource_roles
            )
            or "none"
        ),
        "",
        "Current scope boundaries",
        f"  Target sequence role: {_detail_state_text(profile.scope.target_role.value)}",
        f"  Reverse liftOver: {_detail_state_text(profile.scope.reverse_result.value)}",
        f"  Flanking-interval assessment: {_detail_state_text(profile.scope.query_context.value)}",
        (
            "  Standard liftOver/all-chain comparison: "
            f"{_detail_comparative_relationship_text(profile.scope.comparative_relationship)}"
        ),
        f"  Batch relationships: {_detail_state_text(profile.scope.batch_relationship.value)}",
        f"  Segmental Duplications context: {_detail_state_text(profile.scope.external_context.value)}",
        "  Named variant / rsID identity: not assessed",
        "  Gene / transcript identity: not assessed",
        "  File / downstream workflow: not assessed",
    ]

    if report.source_preflight is not None:
        preflight = report.source_preflight
        lines.extend(
            (
                "",
                "Source validation details",
                f"  State: {_detail_state_text(preflight.state.value)}",
                f"  Canonical sequence: {preflight.canonical_sequence_name}",
                f"  Authoritative sequence length: {preflight.sequence_length}",
                "  Metadata provenance: "
                + ", ".join(
                    source.source_id for source in preflight.provenance_sources
                ),
            )
        )
        for resource in report.source_preflight_resources:
            lines.extend(
                (
                    f"  Metadata resource: {resource.source_url}",
                    f"    Cache path: {resource.path}",
                    f"    Retrieved at: {resource.retrieved_at}",
                    f"    Size: {resource.size_bytes} bytes",
                    f"    SHA-256: {resource.sha256}",
                )
            )

    if profile.scope.target_role is not TargetRoleState.NOT_ASSESSED:
        lines.extend(("", "Target sequence role"))
        lines.append(f"  State: {_detail_state_text(profile.scope.target_role.value)}")
        for item in profile.target_sequence_roles:
            lines.extend(_target_role_detail_lines(item))
        if report.target_role_metadata is not None:
            metadata = report.target_role_metadata
            lines.extend(
                (
                    f"  UCSC database: {metadata.db}",
                    f"  Version-matched NCBI assembly: {metadata.assembly_accession}",
                    "  UCSC assembly description:",
                    f"    Source URL: {metadata.assembly_description.source_url}",
                    f"    Cache path: {metadata.assembly_description.path}",
                    f"    SHA-256: {metadata.assembly_description.sha256}",
                    "  NCBI sequence report:",
                    f"    Source URL: {metadata.sequence_report.source_url}",
                    f"    Archive member: {metadata.sequence_report.archive_member}",
                    f"    Cache path: {metadata.sequence_report.path}",
                    f"    SHA-256: {metadata.sequence_report.sha256}",
                )
            )

    if profile.scope.external_context is not ExternalContextState.NOT_ASSESSED:
        lines.extend(("", "UCSC Segmental Duplications context"))
        lines.extend(_segmental_duplication_detail_lines(report))

    lines.extend(
        (
            "",
            "Mappings",
            (
                "Mapping order is preserved for reproducibility and does not "
                "indicate rank or preference."
            ),
        )
    )

    if not report.candidates:
        lines.append("  none")
    for candidate, candidate_profile in zip(
        report.candidates,
        profile.candidate_profiles,
        strict=True,
    ):
        lines.extend(_candidate_detail_lines(candidate, candidate_profile))

    if report.query_context_result is not None:
        lines.extend(("", "Flanking-interval assessment"))
        lines.extend(_query_context_detail_lines(report))

    if report.filtered_all_chain_comparison is not None:
        lines.extend(("", "Standard liftOver/all-chain comparison"))
        lines.extend(_comparative_detail_lines(report))

    if report.reverse_mapping_results is not None:
        lines.extend(("", "Reverse liftOver results"))
        if not report.reverse_mapping_results:
            lines.append("  none")
        for result in report.reverse_mapping_results:
            lines.extend(_reverse_mapping_detail_lines(result))

    lines.extend(("", "Resources"))
    for assessment_resource in report.resources:
        resource = assessment_resource.resource
        consumption = (
            "consumed" if assessment_resource.consumed_by_engine else "not consumed"
        )
        lines.extend(
            (
                f"{_capitalize_first(_detail_resource_role_text(assessment_resource.role.value, evidence_tier=profile.evidence_tier))} [{consumption}]",
                f"  Source URL: {resource.source_url}",
                f"  Cache path: {resource.path}",
                f"  Retrieved at: {resource.retrieved_at}",
                f"  Size: {resource.size_bytes} bytes",
                f"  SHA-256: {resource.sha256}",
                f"  Cache hit at acquisition: {'yes' if resource.cache_hit else 'no'}",
                f"  Provider checksum: {_provider_checksum_text(resource)}",
                (
                    "  File provenance: "
                    + (
                        assessment_resource.file_provenance.source_id
                        if assessment_resource.file_provenance is not None
                        else "none"
                    )
                ),
            )
        )

    if report.reverse_mapping_resource is not None:
        lines.extend(("", "Reverse mapping resource"))
        reverse_resource = report.reverse_mapping_resource
        resource = reverse_resource.resource
        lines.extend(
            (
                (
                    f"{report.target_db}->{report.source_db} "
                    f"{_capitalize_first(_detail_resource_role_text(reverse_resource.role.value, evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY))} [consumed]"
                ),
                f"  Source URL: {resource.source_url}",
                f"  Cache path: {resource.path}",
                f"  Retrieved at: {resource.retrieved_at}",
                f"  Size: {resource.size_bytes} bytes",
                f"  SHA-256: {resource.sha256}",
                f"  Provider checksum: {_provider_checksum_text(resource)}",
                (
                    "  File provenance: "
                    + (
                        reverse_resource.file_provenance.source_id
                        if reverse_resource.file_provenance is not None
                        else "none"
                    )
                ),
            )
        )

    if report.filtered_chain_comparison_resource is not None:
        lines.extend(("", "Filtered-chain comparison resource"))
        comparison_resource = report.filtered_chain_comparison_resource
        resource = comparison_resource.resource
        lines.extend(
            (
                (
                    f"{report.source_db}->{report.target_db} ordinary filtered "
                    "liftOver chain [consumed for paired comparison]"
                ),
                f"  Source URL: {resource.source_url}",
                f"  Cache path: {resource.path}",
                f"  Retrieved at: {resource.retrieved_at}",
                f"  Size: {resource.size_bytes} bytes",
                f"  SHA-256: {resource.sha256}",
                f"  Provider checksum: {_provider_checksum_text(resource)}",
                (
                    "  File provenance: "
                    + (
                        comparison_resource.file_provenance.source_id
                        if comparison_resource.file_provenance is not None
                        else "none"
                    )
                ),
            )
        )

    lines.extend(("", "Provenance dependency graph"))
    for source in _report_provenance_sources(report):
        lines.append(source.source_id)
        lines.append(f"  Label: {source.label}")
        identifiers = ", ".join(
            f"{identifier.kind.value}={identifier.value}"
            for identifier in source.identifiers
        )
        lines.append(f"  Identifiers: {identifiers or 'none'}")
        parents = ", ".join(parent.source_id for parent in source.derived_from)
        lines.append(f"  Derived from: {parents or 'none'}")

    lines.extend(
        (
            "",
            (
                "Dependency note: provenance edges record shared upstream dependence; "
                "they do not establish independent confirmation."
            ),
            "",
            _BIOLOGICAL_CORRECTNESS_CAVEAT,
        )
    )
    return "\n".join(lines)


def render_assessment_json(report: UCSCAssessmentReport) -> str:
    """Render the schema-v2 factual result report."""

    payload: dict[str, object] = {
        "schema_version": _JSON_SCHEMA_VERSION,
        "report_type": _JSON_REPORT_TYPE,
        "semantics": {
            "interval_coordinates": _JSON_INTERVAL_COORDINATE_SYSTEM,
            "candidate_order": "reproducibility_only_not_rank",
            "result_dimensions": "orthogonal_not_votes",
            "comparative_relationships": "categorical_not_scores_or_votes",
            "provenance_edges": "dependence_not_independent_confirmation",
            "ucsc_pair_dependency_group": (
                "conservative_grouping_not_processing_run_proof"
            ),
        },
        "ucsc_database_pair": {
            "source_db": report.source_db,
            "target_db": report.target_db,
        },
        "source_assembly": _assembly_json(report.source_interval.assembly),
        "target_assembly": _assembly_json(report.target_assembly),
        "source_interval": _interval_json(report.source_interval),
        "source_preflight": _source_preflight_json(report),
        "target_role_metadata": _target_role_metadata_json(report),
        "typed_external_context": _external_context_json(report),
        "result_profile": _result_profile_json(report.result_profile),
        "candidates": [_candidate_json(candidate) for candidate in report.candidates],
        "query_context": _query_context_json(report),
        "reverse_mapping": _reverse_mapping_json(report),
        "filtered_all_chain_comparison": _filtered_all_chain_comparison_json(report),
        "resources": [
            _assessment_resource_json(assessment_resource)
            for assessment_resource in report.resources
        ],
        "provenance": {
            "alignment_source_id": report.alignment_provenance.source_id,
            "sources": [
                _provenance_source_json(source)
                for source in _report_provenance_sources(report)
            ],
        },
        "caveat": _BIOLOGICAL_CORRECTNESS_CAVEAT,
    }
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False)


def _source_preflight_json(report: UCSCAssessmentReport) -> dict[str, object]:
    preflight = report.source_preflight
    if preflight is None:
        return {
            "state": "NOT_ASSESSED",
            "canonical_sequence_name": None,
            "sequence_length": None,
            "suggested_sequence_name": None,
            "alias_sources": [],
            "provenance_source_ids": [],
            "resources": [],
        }
    return {
        "state": preflight.state.value,
        "canonical_sequence_name": preflight.canonical_sequence_name,
        "sequence_length": preflight.sequence_length,
        "suggested_sequence_name": preflight.suggested_sequence_name,
        "alias_sources": list(preflight.alias_sources),
        "provenance_source_ids": [
            source.source_id for source in preflight.provenance_sources
        ],
        "resources": [
            _preflight_resource_json(resource)
            for resource in report.source_preflight_resources
        ],
    }


def _preflight_resource_json(resource: CachedResource) -> dict[str, object]:
    checksum = resource.provider_checksum
    return {
        "source_url": resource.source_url,
        "cache_path": str(resource.path),
        "retrieved_at": resource.retrieved_at,
        "size_bytes": resource.size_bytes,
        "sha256": resource.sha256,
        "cache_hit_at_acquisition": resource.cache_hit,
        "provider_checksum": (
            {
                "algorithm": checksum.algorithm.value,
                "value": checksum.value,
                "source_url": checksum.source_url,
            }
            if checksum is not None
            else None
        ),
        "terms": {
            "resource_class": resource.terms.resource_class.value,
            "general_terms_url": resource.terms.general_terms_url,
            "directory_terms_url": resource.terms.directory_terms_url,
            "restricted_liftover_chain": resource.terms.restricted_liftover_chain,
        },
    }


def _target_role_detail_lines(item: TargetSequenceRoleProfile) -> list[str]:
    lines = [f"  {item.sequence_name}"]
    if item.context is None:
        lines.append("    Provider role: unavailable for this UCSC target sequence")
    else:
        context = item.context
        lines.extend(
            (
                f"    Assembly accession: {context.assembly_accession}",
                f"    Assembly unit: {context.assembly_unit}",
                f"    Provider role: {context.provider_role}",
                f"    Length: {context.length}",
                f"    Chromosome name: {context.chromosome_name or 'none'}",
                f"    GenBank accession: {context.genbank_accession or 'none'}",
                f"    RefSeq accession: {context.refseq_accession or 'none'}",
            )
        )
    lines.append(f"    Role provenance: {item.provenance_source_id or 'none'}")
    return lines


def _target_sequence_role_json(
    item: TargetSequenceRoleProfile,
) -> dict[str, object]:
    context = item.context
    return {
        "sequence_name": item.sequence_name,
        "metadata_available": context is not None,
        "assembly_accession": context.assembly_accession if context else None,
        "assembly_unit": context.assembly_unit if context else None,
        "provider_role": context.provider_role if context else None,
        "length": context.length if context else None,
        "chromosome_name": context.chromosome_name if context else None,
        "ucsc_style_name": context.ucsc_style_name if context else None,
        "genbank_accession": context.genbank_accession if context else None,
        "refseq_accession": context.refseq_accession if context else None,
        "provenance_source_id": item.provenance_source_id,
    }


def _target_role_metadata_json(report: UCSCAssessmentReport) -> dict[str, object]:
    metadata = report.target_role_metadata
    if metadata is None:
        return {
            "state": report.result_profile.scope.target_role.value,
            "ucsc_database": report.target_db,
            "assembly_accession": None,
            "resources": [],
        }
    return {
        "state": report.result_profile.scope.target_role.value,
        "ucsc_database": metadata.db,
        "assembly_accession": metadata.assembly_accession,
        "resources": [
            {
                "kind": "UCSC_ASSEMBLY_DESCRIPTION",
                "source_url": metadata.assembly_description.source_url,
                "archive_member": metadata.assembly_description.archive_member,
                "cache_path": str(metadata.assembly_description.path),
                "retrieved_at": metadata.assembly_description.retrieved_at,
                "size_bytes": metadata.assembly_description.size_bytes,
                "sha256": metadata.assembly_description.sha256,
                "cache_hit_at_acquisition": metadata.assembly_description.cache_hit,
            },
            {
                "kind": "NCBI_SEQUENCE_REPORT",
                "source_url": metadata.sequence_report.source_url,
                "archive_member": metadata.sequence_report.archive_member,
                "cache_path": str(metadata.sequence_report.path),
                "retrieved_at": metadata.sequence_report.retrieved_at,
                "size_bytes": metadata.sequence_report.size_bytes,
                "sha256": metadata.sequence_report.sha256,
                "cache_hit_at_acquisition": metadata.sequence_report.cache_hit,
            },
        ],
    }


def _segmental_duplication_detail_lines(report: UCSCAssessmentReport) -> list[str]:
    result = report.segmental_duplication_context_result
    if result is None:
        return ["  State: NOT_ASSESSED"]

    lines = [
        f"  Overall state: {_detail_state_text(report.result_profile.scope.external_context.value)}",
        f"  Source query state: {_detail_state_text(result.source_state.value)}",
        f"  Target mapping state: {_detail_state_text(result.target_state.value)}",
    ]
    if report.source_segmental_duplication_resource is not None:
        resource = report.source_segmental_duplication_resource
        lines.extend(
            (
                "  Source track resource:",
                f"    Source URL: {resource.source_url}",
                f"    Cache path: {resource.path}",
                f"    SHA-256: {resource.sha256}",
            )
        )
    if report.target_segmental_duplication_resource is not None:
        resource = report.target_segmental_duplication_resource
        lines.extend(
            (
                "  Target track resource:",
                f"    Source URL: {resource.source_url}",
                f"    Cache path: {resource.path}",
                f"    SHA-256: {resource.sha256}",
            )
        )

    lines.append(f"  Source overlap rows: {len(result.source_overlaps)}")
    for source_overlap in result.source_overlaps:
        lines.extend(_source_segmental_duplication_detail(source_overlap))
    lines.append(
        f"  Target mapping/row overlap observations: {len(result.target_overlaps)}"
    )
    for target_overlap in result.target_overlaps:
        lines.extend(_target_segmental_duplication_detail(target_overlap))
    lines.append(
        "  Interpretation boundary: overlap is descriptive context only; it does not "
        "penalize a mapping or establish biological correctness."
    )
    return lines


def _source_segmental_duplication_detail(
    item: UCSCSegmentalDuplicationOverlap,
) -> list[str]:
    return [
        f"    Overlap: {format_display_interval(item.overlap_interval)}",
        f"      Track interval: {format_display_interval(item.record.interval)}",
        (
            "      Paired interval: "
            f"{format_display_interval(item.record.paired_interval)}"
        ),
        f"      Pair strand: {item.record.strand}",
        f"      Provider UID: {item.record.uid}",
        f"      Aligned bases: {item.record.aligned_bases}",
        f"      Fraction matching bases: {item.record.fraction_matching_bases}",
    ]


def _target_segmental_duplication_detail(
    item: CandidateSegmentalDuplicationOverlap,
) -> list[str]:
    return [
        f"    Mapping: {item.candidate_id}",
        "      Exact mapped overlap: "
        + ", ".join(
            format_display_interval(interval) for interval in item.overlap_intervals
        ),
        f"      Track interval: {format_display_interval(item.record.interval)}",
        (
            "      Paired interval: "
            f"{format_display_interval(item.record.paired_interval)}"
        ),
        f"      Pair strand: {item.record.strand}",
        f"      Provider UID: {item.record.uid}",
        f"      Aligned bases: {item.record.aligned_bases}",
        f"      Fraction matching bases: {item.record.fraction_matching_bases}",
    ]


def _external_context_json(report: UCSCAssessmentReport) -> dict[str, object]:
    result = report.segmental_duplication_context_result
    resources: list[dict[str, object]] = []
    if report.source_segmental_duplication_resource is not None:
        resources.append(
            _segmental_duplication_resource_json(
                "SOURCE", report.source_segmental_duplication_resource
            )
        )
    if report.target_segmental_duplication_resource is not None:
        resources.append(
            _segmental_duplication_resource_json(
                "TARGET", report.target_segmental_duplication_resource
            )
        )
    if result is None:
        return {
            "state": ExternalContextState.NOT_ASSESSED.value,
            "ucsc_segmental_duplications": None,
            "resources": resources,
        }
    return {
        "state": report.result_profile.scope.external_context.value,
        "ucsc_segmental_duplications": _segmental_duplication_context_json(result),
        "resources": resources,
    }


def _external_context_profile_json(
    profile: ExternalContextProfile,
) -> dict[str, object]:
    return {
        "state": profile.state.value,
        "ucsc_segmental_duplications": (
            _segmental_duplication_context_json(profile.ucsc_segmental_duplication)
            if profile.ucsc_segmental_duplication is not None
            else None
        ),
    }


def _segmental_duplication_context_json(
    result: UCSCSegmentalDuplicationContextResult,
) -> dict[str, object]:
    return {
        "source_query_state": result.source_state.value,
        "target_projection_state": result.target_state.value,
        "source_provenance_source_id": (
            result.source_provenance.source_id if result.source_provenance else None
        ),
        "target_provenance_source_id": (
            result.target_provenance.source_id if result.target_provenance else None
        ),
        "source_overlaps": [
            _source_segmental_duplication_json(item) for item in result.source_overlaps
        ],
        "target_overlaps": [
            _target_segmental_duplication_json(item) for item in result.target_overlaps
        ],
    }


def _source_segmental_duplication_json(
    item: UCSCSegmentalDuplicationOverlap,
) -> dict[str, object]:
    return {
        "overlap_interval": _interval_json(item.overlap_interval),
        "record": _segmental_duplication_record_json(item.record),
    }


def _target_segmental_duplication_json(
    item: CandidateSegmentalDuplicationOverlap,
) -> dict[str, object]:
    return {
        "candidate_id": item.candidate_id,
        "overlap_intervals": [
            _interval_json(interval) for interval in item.overlap_intervals
        ],
        "record": _segmental_duplication_record_json(item.record),
    }


def _segmental_duplication_record_json(
    record: UCSCSegmentalDuplicationRecord,
) -> dict[str, object]:
    return {
        "interval": _interval_json(record.interval),
        "paired_interval": _interval_json(record.paired_interval),
        "strand": record.strand,
        "uid": record.uid,
        "aligned_bases": record.aligned_bases,
        "fraction_matching_bases": record.fraction_matching_bases,
    }


def _segmental_duplication_resource_json(
    side: str, resource: CachedResource
) -> dict[str, object]:
    return {
        "side": side,
        "source_url": resource.source_url,
        "cache_path": str(resource.path),
        "retrieved_at": resource.retrieved_at,
        "size_bytes": resource.size_bytes,
        "sha256": resource.sha256,
        "cache_hit_at_acquisition": resource.cache_hit,
    }


def _result_profile_json(profile: ResultProfile) -> dict[str, object]:
    return {
        "input_validity": profile.input_validity.value,
        "headline": profile.headline.value,
        "interpretation": profile.interpretation,
        "projection_count": profile.projection_count.value,
        "orientation": profile.orientation.value,
        "source_coverage": {
            "state": profile.source_coverage.value,
            "maximum_candidate_covered_source_bases": (
                profile.maximum_candidate_covered_source_bases
            ),
            "source_bases": profile.source_bases,
            "maximum_coverage_candidate_ids": list(
                profile.maximum_coverage_candidate_ids
            ),
            "union_covered_source_bases": profile.union_covered_source_bases,
        },
        "evidence": {
            "tier": profile.evidence_tier.value,
            "consumed_resource_roles": list(profile.consumed_resource_roles),
        },
        "candidate_profiles": [
            _candidate_profile_json(candidate)
            for candidate in profile.candidate_profiles
        ],
        "query_context": _query_context_profile_json(profile.query_context),
        "comparative_relationship": _comparative_relationship_profile_json(
            profile.comparative_relationship
        ),
        "target_role": {
            "state": profile.scope.target_role.value,
            "sequences": [
                _target_sequence_role_json(item)
                for item in profile.target_sequence_roles
            ],
        },
        "external_context": _external_context_profile_json(profile.external_context),
        "scope": {
            "target_role": profile.scope.target_role.value,
            "actual_reverse_mapping": profile.scope.reverse_result.value,
            "query_context": profile.scope.query_context.value,
            "comparative_relationship": profile.scope.comparative_relationship.value,
            "batch_relationship": profile.scope.batch_relationship.value,
            "external_context": profile.scope.external_context.value,
            "named_variant_identity_assessed": (
                profile.scope.named_variant_identity_assessed
            ),
            "gene_transcript_identity_assessed": (
                profile.scope.gene_transcript_identity_assessed
            ),
            "downstream_workflow_assessed": profile.scope.downstream_workflow_assessed,
        },
    }


def _comparative_relationship_profile_json(
    profile: ComparativeRelationshipProfile,
) -> dict[str, object]:
    return {
        "state": profile.state.value,
        "inventory_state": (
            profile.inventory_state.value
            if profile.inventory_state is not None
            else None
        ),
        "favored_candidate_id": profile.favored_candidate_id,
        "additional_all_chain_candidate_ids": list(
            profile.additional_all_chain_candidate_ids
        ),
        "placement_support": [
            {
                "candidate_id": item.candidate_id,
                "complete_source_coverage": item.complete_source_coverage,
                "retained_by_filtered_chain": item.retained_by_filtered_chain,
                "depth1_top_net": item.depth1_top_net,
                "full_reciprocal_best": item.full_reciprocal_best,
            }
            for item in profile.placement_support
        ],
    }


def _filtered_all_chain_comparison_json(
    report: UCSCAssessmentReport,
) -> dict[str, object]:
    comparison = report.filtered_all_chain_comparison
    relationship = report.comparative_evidence_relationship
    resource = report.filtered_chain_comparison_resource
    if comparison is None:
        if relationship is not None or resource is not None:
            raise ValueError(
                "unassessed filtered/all-chain comparison cannot carry relationship "
                "or resource state"
            )
        return {"assessed": False}
    if relationship is None or resource is None:
        raise ValueError(
            "assessed filtered/all-chain comparison requires relationship and resource"
        )

    return {
        "assessed": True,
        "inventory_state": comparison.relationship.value,
        "categorical_relationship": relationship.relationship.value,
        "favored_candidate_id": relationship.favored_candidate_id,
        "all_chain_candidate_ids": list(comparison.all_chain_candidate_ids),
        "filtered_candidate_ids": [
            candidate.candidate_id for candidate in comparison.filtered_candidates
        ],
        "candidate_matches": [
            {
                "filtered_candidate_id": match.filtered_candidate_id,
                "all_chain_candidate_id": match.all_chain_candidate_id,
            }
            for match in comparison.candidate_matches
        ],
        "additional_all_chain_candidate_ids": list(
            comparison.additional_all_chain_candidate_ids
        ),
        "filtered_candidates": [
            _candidate_json(candidate) for candidate in comparison.filtered_candidates
        ],
        "filtered_chain_resource": _assessment_resource_json(resource),
        "provenance": {
            "all_chain_source_id": comparison.all_chain_provenance.source_id,
            "filtered_chain_source_id": comparison.filtered_chain_provenance.source_id,
            "shared_ucsc_pair_dependency_source_ids": [
                parent.source_id
                for parent in comparison.all_chain_provenance.derived_from
            ],
            "shared_processing_run_provenance_verified": False,
        },
    }


def _query_context_profile_json(profile: QueryContextProfile) -> dict[str, object]:
    return {
        "check_state": profile.check_state.value,
        "findings": [finding.value for finding in profile.findings],
        "requested_window_bases": profile.requested_window_bases,
        "tested_source_interval": (
            _interval_json(profile.tested_source_interval)
            if profile.tested_source_interval is not None
            else None
        ),
        "actual_window_bases": profile.actual_window_bases,
        "not_run_reason": (
            profile.not_run_reason.value if profile.not_run_reason is not None else None
        ),
        "projection_count": (
            profile.projection_count.value
            if profile.projection_count is not None
            else None
        ),
        "source_coverage": (
            profile.source_coverage.value
            if profile.source_coverage is not None
            else None
        ),
        "maximum_candidate_covered_source_bases": (
            profile.maximum_candidate_covered_source_bases
        ),
        "union_covered_source_bases": profile.union_covered_source_bases,
        "headline": profile.headline.value if profile.headline is not None else None,
        "point_and_local_context_map_together": (
            profile.point_and_local_context_map_together
        ),
        "candidate_profiles": [
            _candidate_profile_json(candidate)
            for candidate in profile.candidate_profiles
        ],
    }


def _query_context_json(report: UCSCAssessmentReport) -> dict[str, object]:
    result = report.query_context_result
    profile = report.result_profile.query_context
    return {
        "check_state": profile.check_state.value,
        "evidence_scope": "forward_chain_only",
        "requested_window_bases": profile.requested_window_bases,
        "tested_source_interval": (
            _interval_json(profile.tested_source_interval)
            if profile.tested_source_interval is not None
            else None
        ),
        "actual_window_bases": profile.actual_window_bases,
        "not_run_reason": (
            profile.not_run_reason.value if profile.not_run_reason is not None else None
        ),
        "findings": [finding.value for finding in profile.findings],
        "candidates": (
            [_candidate_json(candidate) for candidate in result.candidates]
            if result is not None and result.check_state is QueryContextState.RUN
            else []
        ),
    }


def _candidate_profile_json(profile: CandidateResultProfile) -> dict[str, object]:
    return {
        "candidate_id": profile.candidate_id,
        "source_coverage": {
            "state": profile.coverage_state.value,
            "covered_source_bases": profile.covered_source_bases,
            "source_bases": profile.source_bases,
            "uncovered_source_intervals": [
                _interval_json(interval)
                for interval in profile.uncovered_source_intervals
            ],
            "largest_uncovered_source_span_bases": (
                profile.largest_uncovered_source_span_bases
            ),
        },
        "geometry": {
            "exact_mapped_segment_count": profile.exact_mapped_segment_count,
            "geometric_segment_count": profile.geometric_segment_count,
            "fragmented": profile.fragmented,
            "target_discontinuous": profile.target_discontinuous,
            "target_bounding_span": _interval_json(profile.target_bounding_span),
            "source_gap_intervals": [
                _interval_json(interval) for interval in profile.source_gap_intervals
            ],
            "target_gap_intervals": [
                _interval_json(interval) for interval in profile.target_gap_intervals
            ],
            "largest_source_gap_bases": profile.largest_source_gap_bases,
            "largest_target_gap_bases": profile.largest_target_gap_bases,
        },
        "orientation": profile.orientation.value,
        "reverse_mapping": _reverse_profile_json(profile.reverse_mapping),
    }


def _reverse_profile_json(
    profile: CandidateReverseMappingProfile,
) -> dict[str, object]:
    return {
        "check_state": profile.check_state.value,
        "relationship": (
            profile.relationship.value if profile.relationship is not None else None
        ),
        "original_source_bases": profile.original_source_bases,
        "original_source_covered_bases": profile.original_source_covered_bases,
        "original_source_coverage": (
            profile.original_source_coverage.value
            if profile.original_source_coverage is not None
            else None
        ),
        "exact_original_geometry_return": profile.exact_original_geometry_return,
        "reverse_projection_count": profile.reverse_projection_count,
        "segments_with_reverse_projection": profile.segments_with_reverse_projection,
        "queried_target_segments": [
            _interval_json(interval) for interval in profile.queried_target_segments
        ],
    }


def _reverse_mapping_json(report: UCSCAssessmentReport) -> dict[str, object]:
    results = report.reverse_mapping_results
    return {
        "check_state": report.result_profile.scope.reverse_result.value,
        "reverse_database_pair": (
            {
                "source_db": report.target_db,
                "target_db": report.source_db,
            }
            if results is not None
            and any(result.check_state is ReverseCheckState.RUN for result in results)
            else None
        ),
        "resource": (
            _assessment_resource_json(report.reverse_mapping_resource)
            if report.reverse_mapping_resource is not None
            else None
        ),
        "candidate_results": (
            [_candidate_reverse_mapping_json(result) for result in results]
            if results is not None
            else []
        ),
    }


def _candidate_reverse_mapping_json(
    result: CandidateReverseMappingResult,
) -> dict[str, object]:
    return {
        "forward_candidate_id": result.forward_candidate_id,
        "check_state": result.check_state.value,
        "relationship": (
            result.relationship.value if result.relationship is not None else None
        ),
        "original_source_bases": result.original_source_bases,
        "original_source_covered_bases": (
            result.original_source_covered_bases
            if result.check_state is ReverseCheckState.RUN
            else None
        ),
        "original_source_coverage": (
            result.original_source_coverage.value
            if result.check_state is ReverseCheckState.RUN
            else None
        ),
        "exact_original_geometry_return": (
            result.exact_original_geometry_return
            if result.check_state is ReverseCheckState.RUN
            else None
        ),
        "queried_target_segments": [
            _interval_json(interval) for interval in result.queried_target_segments
        ],
        "segment_results": [
            {
                "queried_target_segment": _interval_json(
                    segment_result.queried_target_segment
                ),
                "expected_original_source_segment": _interval_json(
                    segment_result.expected_original_source_segment
                ),
                "reverse_candidates": [
                    _candidate_json(candidate)
                    for candidate in segment_result.candidates
                ],
            }
            for segment_result in result.segment_results
        ],
    }


def assembly_json_payload(assembly: AssemblyIdentifier) -> dict[str, object]:
    """Return the canonical schema-v2 assembly payload for package reporters."""

    return _assembly_json(assembly)


def interval_json_payload(interval: GenomicInterval) -> dict[str, object]:
    """Return the canonical schema-v2 interval payload for package reporters."""

    return _interval_json(interval)


def candidate_json_payload(candidate: NormalizedCandidate) -> dict[str, object]:
    """Return the canonical schema-v2 candidate payload for package reporters."""

    return _candidate_json(candidate)


def provenance_source_json_payload(source: ProvenanceSource) -> dict[str, object]:
    """Return the canonical schema-v2 provenance payload for package reporters."""

    return _provenance_source_json(source)


def _assembly_json(assembly: AssemblyIdentifier) -> dict[str, object]:
    # Kept as a helper so assembly serialization is identical inside/outside intervals.
    return {
        "name": assembly.name,
        "provider": assembly.provider,
        "accession": assembly.accession,
        "aliases": list(assembly.aliases),
    }


def _interval_json(interval: GenomicInterval) -> dict[str, object]:
    return {
        "assembly": _assembly_json(interval.assembly),
        "sequence_name": interval.sequence_name,
        "start": interval.start,
        "end": interval.end,
        "coordinate_system": _JSON_INTERVAL_COORDINATE_SYSTEM,
    }


def _candidate_json(candidate: NormalizedCandidate) -> dict[str, object]:
    return {
        "candidate_id": candidate.candidate_id,
        "ucsc_chain_id": _chain_id_from_candidate(candidate),
        "orientation": candidate.orientation.value,
        "target_bounding_interval": _interval_json(candidate.target_interval),
        "mapping_provenance_source_id": candidate.mapping_provenance.source_id,
        "segments": [
            {
                "source_interval": _interval_json(segment.source_interval),
                "target_interval": _interval_json(segment.target_interval),
            }
            for segment in candidate.segments
        ],
        "evidence": [
            _evidence_observation_json(observation)
            for observation in candidate.evidence
        ],
    }


def _evidence_observation_json(observation: EvidenceObservation) -> dict[str, object]:
    return {
        "observation_id": observation.observation_id,
        "kind": observation.kind.value,
        "value": _evidence_value_json(observation.value),
        "provenance_source_id": observation.provenance.source_id,
    }


def _evidence_value_json(value: EvidenceValue) -> dict[str, object]:
    if isinstance(value, MappingCoverageSummary):
        return {
            "type": "MAPPING_COVERAGE_SUMMARY",
            "status": value.status.value,
            "covered_source_bases": value.covered_source_bases,
            "source_bases": value.source_bases,
            "uncovered_source_intervals": [
                _interval_json(interval)
                for interval in value.uncovered_source_intervals
            ],
        }

    if isinstance(value, ChainGapSummary):
        return {
            "type": "CHAIN_GAP_SUMMARY",
            "gaps": [
                {
                    "source_boundary_0_based": gap.source_boundary,
                    "source_gap_overlap": (
                        _interval_json(gap.source_gap_overlap)
                        if gap.source_gap_overlap is not None
                        else None
                    ),
                    "target_gap_interval": (
                        _interval_json(gap.target_gap_interval)
                        if gap.target_gap_interval is not None
                        else None
                    ),
                }
                for gap in value.gaps
            ],
        }

    if isinstance(value, NetHierarchySummary):
        return {
            "type": "NET_HIERARCHY_SUMMARY",
            "depth": value.depth,
            "source_fill_interval": _interval_json(value.source_fill_interval),
        }

    if isinstance(value, ReciprocalBestMembershipSummary):
        return {
            "type": "RECIPROCAL_BEST_MEMBERSHIP_SUMMARY",
            "status": value.status.value,
            "resource_completeness": value.resource_completeness.value,
            "chains_examined": value.chains_examined,
            "covered_source_bases": value.covered_source_bases,
            "candidate_source_bases": value.candidate_source_bases,
            "covered_source_intervals": [
                _interval_json(interval) for interval in value.covered_source_intervals
            ],
        }

    if isinstance(value, (str, int, float, bool)):
        return {"type": "SCALAR", "value": value}

    raise TypeError(f"unsupported evidence value for JSON reporting: {type(value)!r}")


def _assessment_resource_json(
    assessment_resource: UCSCAssessmentResource,
) -> dict[str, object]:
    role = assessment_resource.role
    resource = assessment_resource.resource
    file_provenance = assessment_resource.file_provenance
    checksum = resource.provider_checksum
    return {
        "role": role.value,
        "consumed_by_engine": assessment_resource.consumed_by_engine,
        "file_provenance_source_id": (
            file_provenance.source_id if file_provenance is not None else None
        ),
        "source_url": resource.source_url,
        "cache_path": str(resource.path),
        "retrieved_at": resource.retrieved_at,
        "size_bytes": resource.size_bytes,
        "sha256": resource.sha256,
        "cache_hit_at_acquisition": resource.cache_hit,
        "provider_checksum": (
            {
                "algorithm": checksum.algorithm.value,
                "value": checksum.value,
                "source_url": checksum.source_url,
            }
            if checksum is not None
            else None
        ),
        "terms": {
            "resource_class": resource.terms.resource_class.value,
            "general_terms_url": resource.terms.general_terms_url,
            "directory_terms_url": resource.terms.directory_terms_url,
            "restricted_liftover_chain": resource.terms.restricted_liftover_chain,
        },
    }


def _provenance_source_json(source: ProvenanceSource) -> dict[str, object]:
    return {
        "source_id": source.source_id,
        "label": source.label,
        "identifiers": [
            {"kind": identifier.kind.value, "value": identifier.value}
            for identifier in source.identifiers
        ],
        "derived_from_source_ids": [parent.source_id for parent in source.derived_from],
    }


def _candidate_detail_lines(
    candidate: NormalizedCandidate,
    profile: CandidateResultProfile,
    *,
    include_reverse_mapping: bool = True,
) -> list[str]:
    lines = [
        _candidate_heading(candidate),
        f"  Mapping ID: {candidate.candidate_id}",
        f"  Target: {_candidate_text(candidate, profile)}",
        (
            f"  Source coverage: {profile.covered_source_bases}/{profile.source_bases} "
            f"({_detail_state_text(profile.coverage_state.value)})"
        ),
        f"  Mapped segments: {profile.geometric_segment_count}",
        f"  Fragmented: {'yes' if profile.fragmented else 'no'}",
        f"  Target discontinuous: {'yes' if profile.target_discontinuous else 'no'}",
        (
            "  Largest uncovered source span: "
            f"{profile.largest_uncovered_source_span_bases} bases"
        ),
        f"  Largest source chain gap: {profile.largest_source_gap_bases} bases",
        f"  Largest target gap: {profile.largest_target_gap_bases} bases",
        f"  Mapping provenance: {candidate.mapping_provenance.source_id}",
        (f"  Mapped chain alignment blocks ({profile.exact_mapped_segment_count}):"),
    ]
    for segment in candidate.segments:
        lines.append(
            "    "
            f"{format_display_interval(segment.source_interval)} -> "
            f"{format_display_interval(segment.target_interval)}"
        )

    if profile.uncovered_source_intervals:
        lines.append("  Uncovered source intervals:")
        lines.extend(
            f"    {format_display_interval(interval)}"
            for interval in profile.uncovered_source_intervals
        )
    else:
        lines.append("  Uncovered source intervals: none")

    if profile.source_gap_intervals:
        lines.append("  Source chain-gap intervals:")
        lines.extend(
            f"    {format_display_interval(interval)}"
            for interval in profile.source_gap_intervals
        )
    else:
        lines.append("  Source chain-gap intervals: none")

    if profile.target_gap_intervals:
        lines.append("  Target gap intervals:")
        lines.extend(
            f"    {format_display_interval(interval)}"
            for interval in profile.target_gap_intervals
        )
    else:
        lines.append("  Target gap intervals: none")

    if include_reverse_mapping:
        reverse = profile.reverse_mapping
        lines.append(
            f"  Reverse mapping check: {_detail_state_text(reverse.check_state.value)}"
        )
        if reverse.check_state is ReverseCheckState.RUN:
            assert reverse.relationship is not None
            assert reverse.original_source_covered_bases is not None
            assert reverse.original_source_coverage is not None
            assert reverse.exact_original_geometry_return is not None
            assert reverse.reverse_projection_count is not None
            assert reverse.segments_with_reverse_projection is not None
            lines.extend(
                (
                    f"  Reverse result: {_detail_reverse_relationship_text(reverse.relationship)}",
                    (
                        "  Reverse original-source coverage: "
                        f"{reverse.original_source_covered_bases}/"
                        f"{reverse.original_source_bases} "
                        f"({_detail_state_text(reverse.original_source_coverage.value)})"
                    ),
                    (
                        "  Exact original aligned geometry reconstructed: "
                        f"{'yes' if reverse.exact_original_geometry_return else 'no'}"
                    ),
                    f"  Reverse mappings: {reverse.reverse_projection_count}",
                    (
                        "  Forward target segments with reverse mapping: "
                        f"{reverse.segments_with_reverse_projection}/"
                        f"{len(reverse.queried_target_segments)}"
                    ),
                )
            )

    lines.append(f"  Evidence ({len(candidate.evidence)}):")
    for observation in candidate.evidence:
        value_lines = _evidence_value_lines(observation)
        lines.append(
            f"    {_detail_evidence_kind_text(observation.kind.value)}: {value_lines[0]}"
        )
        lines.extend(f"      {line}" for line in value_lines[1:])
        lines.append(f"      provenance: {observation.provenance.source_id}")
    return lines


def _query_context_detail_lines(report: UCSCAssessmentReport) -> list[str]:
    result = report.query_context_result
    if result is None:
        return ["  Check state: not performed"]
    profile = report.result_profile.query_context
    lines = [
        f"  Check state: {_detail_state_text(profile.check_state.value)}",
        f"  Requested window: {profile.requested_window_bases} bases",
        "  Evidence scope: forward liftOver chain only; net and reciprocal-best resources were not reassessed",
    ]
    if profile.check_state is QueryContextState.NOT_RUN:
        if profile.not_run_reason is None:
            raise ValueError("unperformed query context requires a not-run reason")
        lines.append(
            f"  Not-run reason: {_detail_state_text(profile.not_run_reason.value)}"
        )
        return lines

    assert profile.tested_source_interval is not None
    assert profile.actual_window_bases is not None
    assert profile.projection_count is not None
    assert profile.source_coverage is not None
    assert profile.headline is not None
    lines.extend(
        (
            (
                "  Tested source window: "
                f"{format_display_interval(profile.tested_source_interval)}"
            ),
            f"  Actual tested width: {profile.actual_window_bases} bases",
            f"  Flanking-interval headline: {_headline_text(profile.headline)}",
            f"  Flanking-interval mapping count: {len(result.candidates)}",
            f"  Flanking-interval source coverage: {_detail_state_text(profile.source_coverage.value)}",
            (
                "  Maximum flanking-interval source coverage for one mapping: "
                f"{profile.maximum_candidate_covered_source_bases}/"
                f"{profile.actual_window_bases}"
            ),
            (
                "  Union flanking-interval source coverage: "
                f"{profile.union_covered_source_bases}/{profile.actual_window_bases}"
            ),
            "  Findings: "
            + (
                ", ".join(
                    _detail_query_context_finding_text(finding)
                    for finding in profile.findings
                )
                or "none"
            ),
            (
                "  Point and flanking interval map together: "
                f"{'yes' if profile.point_and_local_context_map_together else 'no'}"
            ),
            "  Flanking-interval mappings:",
        )
    )
    if not result.candidates:
        lines.append("    none")
        return lines
    for candidate, candidate_profile in zip(
        result.candidates,
        profile.candidate_profiles,
        strict=True,
    ):
        candidate_lines = _candidate_detail_lines(
            candidate,
            candidate_profile,
            include_reverse_mapping=False,
        )
        lines.extend(f"  {line}" for line in candidate_lines)
    return lines


def _reverse_mapping_detail_lines(
    result: CandidateReverseMappingResult,
) -> list[str]:
    lines = [
        f"Mapping {result.forward_candidate_id}",
        f"  Check state: {_detail_state_text(result.check_state.value)}",
    ]
    if result.check_state is not ReverseCheckState.RUN:
        return lines

    assert result.relationship is not None
    lines.extend(
        (
            f"  Result: {_detail_reverse_relationship_text(result.relationship)}",
            (
                "  Original aligned source coverage: "
                f"{result.original_source_covered_bases}/"
                f"{result.original_source_bases} "
                f"({_detail_state_text(result.original_source_coverage.value)})"
            ),
            (
                "  Exact original aligned geometry reconstructed: "
                f"{'yes' if result.exact_original_geometry_return else 'no'}"
            ),
        )
    )
    for index, segment_result in enumerate(result.segment_results, start=1):
        lines.append(
            f"  Segment {index} reverse query: "
            f"{format_display_interval(segment_result.queried_target_segment)}"
        )
        expected_source = format_display_interval(
            segment_result.expected_original_source_segment
        )
        lines.append(f"    Expected original source: {expected_source}")
        if not segment_result.candidates:
            lines.append("    Reverse mappings: none")
            continue
        lines.append(f"    Reverse mappings: {len(segment_result.candidates)}")
        for candidate in segment_result.candidates:
            lines.append(
                "      "
                f"{candidate.candidate_id}: "
                f"{format_display_interval(candidate.target_interval)}"
            )
    return lines


def _candidate_heading(candidate: NormalizedCandidate) -> str:
    chain_id = _chain_id_from_candidate(candidate)
    if chain_id is None:
        return f"Mapping {candidate.candidate_id}"
    return f"Chain {chain_id}"


def _chain_id_from_candidate(candidate: NormalizedCandidate) -> int | None:
    return chain_id_from_candidate_id(candidate.candidate_id)


def _evidence_value_lines(observation: EvidenceObservation) -> list[str]:
    value = observation.value
    if isinstance(value, MappingCoverageSummary):
        lines = [
            (
                f"{_detail_state_text(value.status.value)}; {value.covered_source_bases}/"
                f"{value.source_bases} source bases covered"
            )
        ]
        if value.uncovered_source_intervals:
            lines.append(
                "uncovered source intervals: "
                + ", ".join(
                    format_display_interval(interval)
                    for interval in value.uncovered_source_intervals
                )
            )
        else:
            lines.append("uncovered source intervals: none")
        return lines

    if isinstance(value, ChainGapSummary):
        lines = [f"{len(value.gaps)} chain gap(s) through the requested locus"]
        for gap in value.gaps:
            source_gap = (
                format_display_interval(gap.source_gap_overlap)
                if gap.source_gap_overlap is not None
                else "none"
            )
            target_gap = (
                format_display_interval(gap.target_gap_interval)
                if gap.target_gap_interval is not None
                else "none"
            )
            lines.append(
                f"source boundary={gap.source_boundary} (0-based boundary); "
                f"source gap={source_gap}; target gap={target_gap}"
            )
        return lines

    if isinstance(value, NetHierarchySummary):
        return [
            (
                f"depth={value.depth}; "
                f"fill span={format_display_interval(value.source_fill_interval)}"
            )
        ]

    if isinstance(value, ReciprocalBestMembershipSummary):
        lines = [
            (
                f"{_detail_state_text(value.status.value)}; {value.covered_source_bases}/"
                f"{value.candidate_source_bases} mapped source bases covered; "
                f"resource completeness={_detail_state_text(value.resource_completeness.value)}; "
                f"chains examined={value.chains_examined}"
            )
        ]
        if value.covered_source_intervals:
            lines.append(
                "covered source intervals: "
                + ", ".join(
                    format_display_interval(interval)
                    for interval in value.covered_source_intervals
                )
            )
        else:
            lines.append("covered source intervals: none")
        return lines

    return [str(value)]


def _provider_checksum_text(resource: CachedResource) -> str:
    checksum = resource.provider_checksum
    if checksum is None:
        return "none"
    return f"{checksum.algorithm.value}:{checksum.value} (from {checksum.source_url})"


def _report_provenance_sources(
    report: UCSCAssessmentReport,
) -> tuple[ProvenanceSource, ...]:
    roots: list[ProvenanceSource] = [report.alignment_provenance]
    if report.source_preflight is not None:
        roots.extend(report.source_preflight.provenance_sources)
    if report.target_role_provenance is not None:
        roots.append(report.target_role_provenance)
    if report.segmental_duplication_context_result is not None:
        context = report.segmental_duplication_context_result
        if context.source_provenance is not None:
            roots.append(context.source_provenance)
        if context.target_provenance is not None:
            roots.append(context.target_provenance)
    for assessment_resource in report.resources:
        if assessment_resource.file_provenance is not None:
            roots.append(assessment_resource.file_provenance)
    for candidate in report.candidates:
        roots.append(candidate.mapping_provenance)
        roots.extend(observation.provenance for observation in candidate.evidence)
    if report.reverse_alignment_provenance is not None:
        roots.append(report.reverse_alignment_provenance)
    if (
        report.reverse_mapping_resource is not None
        and report.reverse_mapping_resource.file_provenance is not None
    ):
        roots.append(report.reverse_mapping_resource.file_provenance)
    if (
        report.filtered_chain_comparison_resource is not None
        and report.filtered_chain_comparison_resource.file_provenance is not None
    ):
        roots.append(report.filtered_chain_comparison_resource.file_provenance)
    if report.filtered_all_chain_comparison is not None:
        roots.append(report.filtered_all_chain_comparison.all_chain_provenance)
        roots.append(report.filtered_all_chain_comparison.filtered_chain_provenance)
        for candidate in report.filtered_all_chain_comparison.filtered_candidates:
            roots.append(candidate.mapping_provenance)
            roots.extend(observation.provenance for observation in candidate.evidence)
    if report.reverse_mapping_results is not None:
        for result in report.reverse_mapping_results:
            for segment_result in result.segment_results:
                for candidate in segment_result.candidates:
                    roots.append(candidate.mapping_provenance)
                    roots.extend(
                        observation.provenance for observation in candidate.evidence
                    )

    by_id: dict[str, ProvenanceSource] = {}
    pending = list(roots)
    while pending:
        source = pending.pop()
        existing = by_id.get(source.source_id)
        if existing is not None:
            if _provenance_definition(existing) != _provenance_definition(source):
                raise ValueError(
                    "provenance source ID refers to conflicting source definitions"
                )
            continue
        by_id[source.source_id] = source
        pending.extend(source.derived_from)
    return tuple(by_id[source_id] for source_id in sorted(by_id))


def _provenance_definition(source: ProvenanceSource) -> tuple[object, ...]:
    return (
        source.label,
        source.identifiers,
        tuple(parent.source_id for parent in source.derived_from),
    )
