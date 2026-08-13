"""Core data structures for liftAssess.

The model deliberately contains no candidate-generation logic, evidence scoring,
or verdict heuristics. Coordinates are represented internally as 0-based,
half-open intervals.
"""

from dataclasses import dataclass
from enum import Enum
from typing import TypeAlias


class Verdict(str, Enum):
    """The three v1 assessment verdicts."""

    WELL_SUPPORTED = "WELL_SUPPORTED"
    CONTESTED = "CONTESTED"
    INDETERMINATE = "INDETERMINATE"


class EvidenceAvailabilityTier(str, Enum):
    """How much source evidence is available, separate from verdict strength."""

    COMPARATIVE = "COMPARATIVE"
    LIFTOVER_ONLY = "LIFTOVER_ONLY"


class MappingOrientation(str, Enum):
    """Relative alignment orientation between source and target spans."""

    SAME = "SAME"
    REVERSE = "REVERSE"


class ProvenanceIdentifierKind(str, Enum):
    """Identifier schemes established by the v1 provenance design."""

    REFGET_SEQUENCE = "REFGET_SEQUENCE"
    SEQCOL = "SEQCOL"
    SHA256 = "SHA256"


class EvidenceKind(str, Enum):
    """Evidence categories explicitly in scope for v1."""

    MAPPING_COVERAGE = "MAPPING_COVERAGE"
    CHAIN_GAPS = "CHAIN_GAPS"
    CANDIDATE_RANK = "CANDIDATE_RANK"
    TARGET_PLACEMENT = "TARGET_PLACEMENT"
    CHAIN_SCORE = "CHAIN_SCORE"
    ALIGNED_BASES = "ALIGNED_BASES"
    DUPLICATED_QUERY_BASES = "DUPLICATED_QUERY_BASES"
    NET_CLASSIFICATION = "NET_CLASSIFICATION"
    NET_HIERARCHY = "NET_HIERARCHY"
    RECIPROCAL_BEST_MEMBERSHIP = "RECIPROCAL_BEST_MEMBERSHIP"
    FLANKING_GENE_SYNTENY = "FLANKING_GENE_SYNTENY"


@dataclass(frozen=True)
class AssemblyIdentifier:
    """A structured assembly identity without attempting alias resolution."""

    name: str
    provider: str
    accession: str | None = None
    aliases: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("assembly name must not be empty")
        if not self.provider:
            raise ValueError("assembly provider must not be empty")
        if self.accession == "":
            raise ValueError("assembly accession must not be empty when provided")
        if any(not alias for alias in self.aliases):
            raise ValueError("assembly aliases must not contain empty values")


@dataclass(frozen=True)
class GenomicInterval:
    """Canonical 0-based, half-open genomic interval."""

    assembly: AssemblyIdentifier
    sequence_name: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if not self.sequence_name:
            raise ValueError("sequence name must not be empty")
        if self.start < 0:
            raise ValueError("interval start must be non-negative")
        if self.end < self.start:
            raise ValueError("interval end must be greater than or equal to start")

    @property
    def length(self) -> int:
        return self.end - self.start

    def contains(self, position: int) -> bool:
        """Return whether a 0-based position lies inside this half-open interval."""

        return self.start <= position < self.end


class MappingCoverageStatus(str, Enum):
    """Whether all bases in the requested source locus are aligned."""

    FULL = "FULL"
    PARTIAL = "PARTIAL"


class ReciprocalBestMembershipStatus(str, Enum):
    """Evidence-detail coverage state, not an assessment verdict."""

    FULL = "FULL"
    PARTIAL = "PARTIAL"
    NONE = "NONE"


class ReciprocalBestResourceCompleteness(str, Enum):
    """Explicit basis for treating a reciprocal-best scan as exhaustive.

    There is deliberately no incomplete value: an incomplete scan must not produce
    exhaustive reciprocal-best membership evidence.
    """

    COMPLETE_RESOURCE = "COMPLETE_RESOURCE"
    COMPLETE_CANDIDATE_SUBSET = "COMPLETE_CANDIDATE_SUBSET"


@dataclass(frozen=True)
class MappingCoverageSummary:
    """Mechanical source-locus coverage for one candidate mapping."""

    status: MappingCoverageStatus
    covered_source_bases: int
    source_bases: int
    uncovered_source_intervals: tuple[GenomicInterval, ...] = ()

    def __post_init__(self) -> None:
        if self.source_bases <= 0:
            raise ValueError("mapping coverage source_bases must be positive")
        if not 0 < self.covered_source_bases <= self.source_bases:
            raise ValueError(
                "mapping coverage covered_source_bases must be positive and no greater "
                "than source_bases"
            )

        uncovered_bases = sum(
            interval.length for interval in self.uncovered_source_intervals
        )
        if uncovered_bases != self.source_bases - self.covered_source_bases:
            raise ValueError(
                "mapping coverage uncovered intervals must account for every "
                "unaligned source base"
            )

        if self.status is MappingCoverageStatus.FULL:
            if self.covered_source_bases != self.source_bases:
                raise ValueError("full mapping coverage must align every source base")
            if self.uncovered_source_intervals:
                raise ValueError("full mapping coverage cannot contain uncovered intervals")
        elif self.covered_source_bases == self.source_bases:
            raise ValueError("partial mapping coverage must leave source bases uncovered")


@dataclass(frozen=True)
class ReciprocalBestMembershipSummary:
    """Exact reciprocal-best coverage of one candidate's aligned source bases.

    The denominator is the candidate's aligned source mapping, not the full
    requested source locus. A candidate may already be partial because of chain
    gaps; reciprocal-best evidence answers a different question: how much of the
    mapping that *does* exist survives UCSC's reciprocal-best netting pipeline.

    ``resource_completeness`` is the caller's explicit completeness claim.
    ``chains_examined`` is audit context only; a chain count cannot itself prove that
    the external resource was complete.
    """

    status: ReciprocalBestMembershipStatus
    resource_completeness: ReciprocalBestResourceCompleteness
    chains_examined: int
    covered_source_bases: int
    candidate_source_bases: int
    covered_source_intervals: tuple[GenomicInterval, ...] = ()

    def __post_init__(self) -> None:
        if self.chains_examined < 0:
            raise ValueError("reciprocal-best chains_examined must be non-negative")
        if self.candidate_source_bases <= 0:
            raise ValueError("reciprocal-best candidate_source_bases must be positive")
        if not 0 <= self.covered_source_bases <= self.candidate_source_bases:
            raise ValueError(
                "reciprocal-best covered_source_bases must be between zero and "
                "candidate_source_bases"
            )

        covered_from_intervals = sum(
            interval.length for interval in self.covered_source_intervals
        )
        if covered_from_intervals != self.covered_source_bases:
            raise ValueError(
                "reciprocal-best covered intervals must account for every covered "
                "source base"
            )

        for previous, current in zip(
            self.covered_source_intervals, self.covered_source_intervals[1:]
        ):
            if current.assembly != previous.assembly or (
                current.sequence_name != previous.sequence_name
            ):
                raise ValueError(
                    "reciprocal-best covered intervals must share one source sequence"
                )
            if current.start < previous.end:
                raise ValueError(
                    "reciprocal-best covered intervals must be ordered and non-overlapping"
                )

        if self.status is ReciprocalBestMembershipStatus.FULL:
            if self.covered_source_bases != self.candidate_source_bases:
                raise ValueError(
                    "full reciprocal-best membership must cover every candidate source base"
                )
        elif self.status is ReciprocalBestMembershipStatus.NONE:
            if self.covered_source_bases != 0 or self.covered_source_intervals:
                raise ValueError(
                    "no reciprocal-best membership cannot contain covered source bases"
                )
        elif not 0 < self.covered_source_bases < self.candidate_source_bases:
            raise ValueError(
                "partial reciprocal-best membership must cover some but not all "
                "candidate source bases"
            )


@dataclass(frozen=True)
class ChainGap:
    """One chain block gap that intersects or lies within the requested locus.

    ``source_gap_overlap`` is the portion of the source-side chain gap that
    overlaps the requested source locus. It is ``None`` for a destination-only
    gap (UCSC ``dt == 0``). ``target_gap_interval`` is the corresponding full
    destination-side gap in forward-reference coordinates and is ``None`` when
    UCSC ``dq == 0``.
    """

    source_boundary: int
    source_gap_overlap: GenomicInterval | None = None
    target_gap_interval: GenomicInterval | None = None

    def __post_init__(self) -> None:
        if self.source_boundary < 0:
            raise ValueError("chain gap source boundary must be non-negative")
        if self.source_gap_overlap is None and self.target_gap_interval is None:
            raise ValueError("chain gap must contain a source or target gap")
        if (
            self.source_gap_overlap is not None
            and self.source_gap_overlap.length <= 0
        ):
            raise ValueError("chain source gap overlap must span at least one base")
        if (
            self.target_gap_interval is not None
            and self.target_gap_interval.length <= 0
        ):
            raise ValueError("chain target gap interval must span at least one base")


@dataclass(frozen=True)
class ChainGapSummary:
    """Exact chain block gaps observed through one requested source locus."""

    gaps: tuple[ChainGap, ...] = ()


@dataclass(frozen=True)
class NetHierarchySummary:
    """Mechanical hierarchy context for one matched UCSC net fill.

    ``depth`` is the fill's indentation-derived hierarchy level in the net.
    ``source_fill_interval`` identifies the exact target-side fill span that
    contributed the observation. Depth is context only; it is not a support
    score or confidence measure.
    """

    depth: int
    source_fill_interval: GenomicInterval

    def __post_init__(self) -> None:
        if self.depth < 1:
            raise ValueError("net hierarchy depth must be at least 1")
        if self.source_fill_interval.length <= 0:
            raise ValueError("net fill interval must span at least one base")


EvidenceValue: TypeAlias = (
    str
    | int
    | float
    | bool
    | MappingCoverageSummary
    | ReciprocalBestMembershipSummary
    | ChainGapSummary
    | NetHierarchySummary
)


@dataclass(frozen=True)
class ProvenanceIdentifier:
    """A typed identifier or digest attached to a provenance source."""

    kind: ProvenanceIdentifierKind
    value: str

    def __post_init__(self) -> None:
        if not self.value:
            raise ValueError("provenance identifier value must not be empty")
        if self.kind is ProvenanceIdentifierKind.SHA256:
            prefix = "sha256:"
            digest = self.value.removeprefix(prefix)
            if (
                not self.value.startswith(prefix)
                or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
            ):
                raise ValueError(
                    "SHA256 provenance identifier must use canonical "
                    "sha256:<64 lowercase hexadecimal characters> form"
                )


@dataclass(frozen=True)
class ProvenanceSource:
    """A source node in the evidence provenance graph.

    ``derived_from`` records upstream source dependence. Two observations may
    point to the same source, or to distinct source nodes that trace to the same
    upstream source. The model deliberately does not store an ``independent``
    flag because independence is established from provenance case by case.
    """

    source_id: str
    label: str
    identifiers: tuple[ProvenanceIdentifier, ...] = ()
    derived_from: tuple["ProvenanceSource", ...] = ()

    def __post_init__(self) -> None:
        if not self.source_id:
            raise ValueError("provenance source_id must not be empty")
        if not self.label:
            raise ValueError("provenance label must not be empty")


@dataclass(frozen=True)
class EvidenceObservation:
    """One evidence observation, distinct from the source that produced it."""

    observation_id: str
    kind: EvidenceKind
    value: EvidenceValue
    provenance: ProvenanceSource

    def __post_init__(self) -> None:
        if not self.observation_id:
            raise ValueError("observation_id must not be empty")


@dataclass(frozen=True)
class MappingSegment:
    """One exact ungapped source-to-target portion of a candidate mapping."""

    source_interval: GenomicInterval
    target_interval: GenomicInterval

    def __post_init__(self) -> None:
        if self.source_interval.length <= 0 or self.target_interval.length <= 0:
            raise ValueError("mapping segments must span at least one base")
        if self.source_interval.length != self.target_interval.length:
            raise ValueError("mapping segment source and target lengths must match")


@dataclass(frozen=True)
class NormalizedCandidate:
    """A candidate mapping consumed by the assessor core.

    ``target_interval`` is the smallest forward-reference bounding span that
    contains every exact ``segment`` target interval. It is a summary span,
    not a claim that bases between segments are aligned.
    """

    candidate_id: str
    target_interval: GenomicInterval
    orientation: MappingOrientation
    mapping_provenance: ProvenanceSource
    segments: tuple[MappingSegment, ...]
    evidence: tuple[EvidenceObservation, ...] = ()

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("candidate_id must not be empty")
        if not self.segments:
            raise ValueError("candidate must contain at least one mapping segment")

        source_assembly = self.segments[0].source_interval.assembly
        source_sequence = self.segments[0].source_interval.sequence_name
        for segment in self.segments:
            if (
                segment.source_interval.assembly != source_assembly
                or segment.source_interval.sequence_name != source_sequence
            ):
                raise ValueError(
                    "candidate mapping segments must share one source sequence"
                )
            if (
                segment.target_interval.assembly != self.target_interval.assembly
                or segment.target_interval.sequence_name
                != self.target_interval.sequence_name
            ):
                raise ValueError(
                    "candidate mapping segments must share the candidate target sequence"
                )

        for previous, current in zip(self.segments, self.segments[1:]):
            if current.source_interval.start < previous.source_interval.end:
                raise ValueError(
                    "candidate mapping segments must be ordered and non-overlapping "
                    "on the source sequence"
                )

            if self.orientation is MappingOrientation.SAME:
                if current.target_interval.start < previous.target_interval.end:
                    raise ValueError(
                        "same-orientation candidate segments must be ordered and "
                        "non-overlapping on the target sequence"
                    )
            elif current.target_interval.end > previous.target_interval.start:
                raise ValueError(
                    "reverse-orientation candidate segments must be ordered and "
                    "non-overlapping on the target sequence"
                )

        expected_target_start = min(
            segment.target_interval.start for segment in self.segments
        )
        expected_target_end = max(
            segment.target_interval.end for segment in self.segments
        )
        if (
            self.target_interval.start != expected_target_start
            or self.target_interval.end != expected_target_end
        ):
            raise ValueError(
                "candidate target_interval must exactly bound its mapping segments"
            )

        observation_ids = [observation.observation_id for observation in self.evidence]
        if len(observation_ids) != len(set(observation_ids)):
            raise ValueError("evidence observation IDs must be unique within a candidate")


@dataclass(frozen=True)
class EvidenceReference:
    """A typed reference from an assessment to candidate evidence."""

    candidate_id: str
    observation_id: str

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("candidate_id must not be empty")
        if not self.observation_id:
            raise ValueError("observation_id must not be empty")


@dataclass(frozen=True)
class Assessment:
    """Assessment result supporting summary and detailed reporting.

    The model stores the verdict and references to the evidence used in the
    assessment. It does not compute a verdict, confidence score, or biological
    truth claim. Verdict-specific semantic constraints belong to the assessor
    logic; this container enforces only referential integrity.
    """

    source_interval: GenomicInterval
    verdict: Verdict
    evidence_tier: EvidenceAvailabilityTier
    candidates: tuple[NormalizedCandidate, ...]
    preferred_candidate_id: str | None = None
    supporting_evidence: tuple[EvidenceReference, ...] = ()
    contradicting_evidence: tuple[EvidenceReference, ...] = ()

    def __post_init__(self) -> None:
        candidate_ids = [candidate.candidate_id for candidate in self.candidates]
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError("candidate IDs must be unique within an assessment")

        candidate_by_id = {
            candidate.candidate_id: candidate for candidate in self.candidates
        }
        if (
            self.preferred_candidate_id is not None
            and self.preferred_candidate_id not in candidate_by_id
        ):
            raise ValueError("preferred candidate must reference an assessment candidate")

        for reference in (*self.supporting_evidence, *self.contradicting_evidence):
            candidate = candidate_by_id.get(reference.candidate_id)
            if candidate is None:
                raise ValueError("evidence reference must name an assessment candidate")
            observation_ids = {
                observation.observation_id for observation in candidate.evidence
            }
            if reference.observation_id not in observation_ids:
                raise ValueError(
                    "evidence reference must name an observation on its candidate"
                )
