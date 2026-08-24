"""Single-pass orchestration for the v1 UCSC candidate-generation engine.

This module connects the already-validated mechanical layers without assigning an
aggregate result verdict.  It consumes chain/net/reciprocal-best record streams once,
generates normalized candidates from chains, and attaches the comparative evidence
that can be established from the supplied records.

The single-pass behavior is deliberate.  UCSC comparative resources can be large,
so orchestration must not require materializing a whole net/chain file or rescanning
one external resource once per candidate.  Only candidate-relevant net fills and
reciprocal-best chains are retained while the corresponding streams are consumed.
"""

from collections.abc import Iterable

from .chain import ChainRecord
from .models import (
    AssemblyIdentifier,
    GenomicInterval,
    MappingOrientation,
    NormalizedCandidate,
    ProvenanceSource,
    ReciprocalBestResourceCompleteness,
)
from .net import NetRecord
from .net_evidence import _annotate_candidate_with_net_records
from .projection import project_interval_through_chain
from .reciprocal_best import _annotate_candidate_with_reciprocal_best_chains


def build_ucsc_candidates(
    source_interval: GenomicInterval,
    chains: Iterable[ChainRecord],
    *,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
    net_records: Iterable[NetRecord] | None = None,
    net_provenance: ProvenanceSource | None = None,
    reciprocal_best_chains: Iterable[ChainRecord] | None = None,
    reciprocal_best_provenance: ProvenanceSource | None = None,
    reciprocal_best_completeness: ReciprocalBestResourceCompleteness | None = None,
) -> tuple[NormalizedCandidate, ...]:
    """Build chain-backed candidates and attach available UCSC comparative evidence.

    ``chains`` is consumed once to generate every candidate overlapping
    ``source_interval``.  Optional net and reciprocal-best streams are likewise
    consumed once across the whole candidate set.  This is an engine-boundary
    operation only: it performs no candidate ranking and assigns no assessment
    verdict.

    Net records and their provenance must be supplied together. Reciprocal-best
    records require both provenance and an explicit completeness claim. A complete
    reciprocal-best stream is filtered during the single pass to the sequence-pair
    and orientation combinations relevant to the generated candidates. The evidence
    preserves the caller's completeness basis for the stream that was exhaustively
    consumed; ``COMPLETE_CANDIDATE_SUBSET`` therefore means complete for *all*
    generated candidates in this engine call.
    """

    _validate_optional_inputs(
        net_records=net_records,
        net_provenance=net_provenance,
        reciprocal_best_chains=reciprocal_best_chains,
        reciprocal_best_provenance=reciprocal_best_provenance,
        reciprocal_best_completeness=reciprocal_best_completeness,
    )

    candidate_entries: list[tuple[int, NormalizedCandidate]] = []
    candidate_ids: set[str] = set()

    for chain in chains:
        candidate = project_interval_through_chain(
            source_interval,
            chain,
            target_assembly=target_assembly,
            mapping_provenance=chain_provenance,
        )
        if candidate is None:
            continue
        if candidate.candidate_id in candidate_ids:
            raise ValueError(
                "chain stream produced duplicate candidate identity: "
                f"{candidate.candidate_id}"
            )
        candidate_ids.add(candidate.candidate_id)
        candidate_entries.append((chain.chain_id, candidate))

    if not candidate_entries:
        return ()

    if net_records is not None:
        assert net_provenance is not None
        candidate_entries = _attach_net_evidence(
            candidate_entries,
            net_records=net_records,
            net_provenance=net_provenance,
        )

    if reciprocal_best_chains is not None:
        assert reciprocal_best_provenance is not None
        assert reciprocal_best_completeness is not None
        candidate_entries = _attach_reciprocal_best_evidence(
            candidate_entries,
            reciprocal_best_chains=reciprocal_best_chains,
            reciprocal_best_provenance=reciprocal_best_provenance,
            reciprocal_best_completeness=reciprocal_best_completeness,
        )

    return tuple(candidate for _, candidate in candidate_entries)


def build_ucsc_chain_candidates_for_intervals(
    source_intervals: Iterable[GenomicInterval],
    chains: Iterable[ChainRecord],
    *,
    target_assembly: AssemblyIdentifier,
    chain_provenance: ProvenanceSource,
) -> tuple[tuple[NormalizedCandidate, ...], ...]:
    """Project many source intervals through one shared chain traversal.

    This is the chain-only shared-traversal primitive used by reverse and later
    multi-query features. ``chains`` is consumed exactly once. Candidate ordering for
    each interval follows the original chain encounter order, matching
    :func:`build_ucsc_candidates` without attaching net or reciprocal-best evidence.
    """

    intervals = tuple(source_intervals)
    if not intervals:
        return ()
    if any(interval.length == 0 for interval in intervals):
        raise ValueError(
            "zero-length source interval projection is not defined for liftAssess v1"
        )

    interval_indices_by_sequence: dict[str, list[int]] = {}
    maximum_end_by_sequence: dict[str, int] = {}
    for index, interval in enumerate(intervals):
        interval_indices_by_sequence.setdefault(interval.sequence_name, []).append(
            index
        )
        maximum_end_by_sequence[interval.sequence_name] = max(
            maximum_end_by_sequence.get(interval.sequence_name, 0), interval.end
        )

    candidates_by_interval: list[list[NormalizedCandidate]] = [[] for _ in intervals]
    candidate_ids_by_interval: list[set[str]] = [set() for _ in intervals]

    for chain in chains:
        interval_indices = interval_indices_by_sequence.get(chain.target_name)
        if interval_indices is None:
            continue
        if maximum_end_by_sequence[chain.target_name] > chain.target_size:
            raise ValueError("source interval exceeds chain target sequence bounds")

        for index in interval_indices:
            interval = intervals[index]
            if interval.end <= chain.target_start or interval.start >= chain.target_end:
                continue
            candidate = project_interval_through_chain(
                interval,
                chain,
                target_assembly=target_assembly,
                mapping_provenance=chain_provenance,
            )
            if candidate is None:
                continue
            candidate_ids = candidate_ids_by_interval[index]
            if candidate.candidate_id in candidate_ids:
                raise ValueError(
                    "chain stream produced duplicate candidate identity: "
                    f"{candidate.candidate_id}"
                )
            candidate_ids.add(candidate.candidate_id)
            candidates_by_interval[index].append(candidate)

    return tuple(tuple(candidates) for candidates in candidates_by_interval)


def _validate_optional_inputs(
    *,
    net_records: Iterable[NetRecord] | None,
    net_provenance: ProvenanceSource | None,
    reciprocal_best_chains: Iterable[ChainRecord] | None,
    reciprocal_best_provenance: ProvenanceSource | None,
    reciprocal_best_completeness: ReciprocalBestResourceCompleteness | None,
) -> None:
    if (net_records is None) != (net_provenance is None):
        raise ValueError("net records and net provenance must be supplied together")

    reciprocal_values = (
        reciprocal_best_chains,
        reciprocal_best_provenance,
        reciprocal_best_completeness,
    )
    if any(value is None for value in reciprocal_values) and any(
        value is not None for value in reciprocal_values
    ):
        raise ValueError(
            "reciprocal-best chains, provenance, and completeness must be supplied together"
        )


def _attach_net_evidence(
    candidate_entries: list[tuple[int, NormalizedCandidate]],
    *,
    net_records: Iterable[NetRecord],
    net_provenance: ProvenanceSource,
) -> list[tuple[int, NormalizedCandidate]]:
    """Attach net evidence after one pass through ``net_records``.

    Chain IDs identify which fills *might* be relevant, but they are not enough to
    establish a match: the candidate-level annotator still checks sequence,
    orientation, and exact source-segment overlap.  Retaining all fills for a
    candidate chain ID is necessary because UCSC nets may contain the same chain ID
    at multiple hierarchy positions.
    """

    candidate_chain_ids = {chain_id for chain_id, _ in candidate_entries}
    records_by_chain_id: dict[int, list[NetRecord]] = {
        chain_id: [] for chain_id in candidate_chain_ids
    }

    for record in net_records:
        chain_id = record.chain_id
        if chain_id is not None and chain_id in records_by_chain_id:
            records_by_chain_id[chain_id].append(record)

    return [
        (
            chain_id,
            _annotate_candidate_with_net_records(
                candidate,
                chain_id=chain_id,
                net_records=records_by_chain_id[chain_id],
                net_provenance=net_provenance,
            ),
        )
        for chain_id, candidate in candidate_entries
    ]


def _attach_reciprocal_best_evidence(
    candidate_entries: list[tuple[int, NormalizedCandidate]],
    *,
    reciprocal_best_chains: Iterable[ChainRecord],
    reciprocal_best_provenance: ProvenanceSource,
    reciprocal_best_completeness: ReciprocalBestResourceCompleteness,
) -> list[tuple[int, NormalizedCandidate]]:
    """Attach reciprocal-best evidence after one exhaustive input-stream pass.

    The caller's completeness claim applies to the input stream and is preserved in
    the resulting evidence; this internal filter does not redefine that scope. During
    the one exhaustive pass, retain every chain whose source/target sequence pair and
    orientation could apply to at least one generated candidate. These are the same
    necessary predicates enforced again by the candidate-level reciprocal-best
    annotator, so a chain excluded here could not contribute membership evidence to
    any generated candidate. The filter is therefore lossless with respect to every
    candidate in this call: ``COMPLETE_RESOURCE`` remains a truthful description when
    the caller supplied the complete resource, while ``COMPLETE_CANDIDATE_SUBSET``
    remains truthful when that narrower scope was supplied.

    As everywhere in liftAssess provenance/completeness handling, the library records
    and propagates the caller's claim; it cannot prove that an external file was not
    truncated before being supplied.
    """

    relevant_keys = {
        _candidate_pair_key(candidate) for _, candidate in candidate_entries
    }
    relevant_chains = [
        chain
        for chain in reciprocal_best_chains
        if _chain_pair_key(chain) in relevant_keys
    ]

    # The evidence basis records what the caller asserts was exhaustively consumed by
    # this engine call. Even though only candidate-relevant chains are retained after
    # that pass, COMPLETE_RESOURCE remains truthful when the full resource was scanned;
    # COMPLETE_CANDIDATE_SUBSET remains truthful when that was the supplied scope.

    return [
        (
            chain_id,
            _annotate_candidate_with_reciprocal_best_chains(
                candidate,
                reciprocal_best_chains=relevant_chains,
                resource_completeness=reciprocal_best_completeness,
                reciprocal_best_provenance=reciprocal_best_provenance,
            ),
        )
        for chain_id, candidate in candidate_entries
    ]


def _candidate_pair_key(
    candidate: NormalizedCandidate,
) -> tuple[str, str, MappingOrientation]:
    return (
        candidate.segments[0].source_interval.sequence_name,
        candidate.target_interval.sequence_name,
        candidate.orientation,
    )


def _chain_pair_key(chain: ChainRecord) -> tuple[str, str, MappingOrientation]:
    return (chain.target_name, chain.query_name, chain.orientation)
