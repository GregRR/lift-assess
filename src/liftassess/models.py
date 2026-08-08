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


EvidenceValue: TypeAlias = str | int | float | bool


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


@dataclass(frozen=True)
class ProvenanceIdentifier:
    """A typed identifier or digest attached to a provenance source."""

    kind: ProvenanceIdentifierKind
    value: str

    def __post_init__(self) -> None:
        if not self.value:
            raise ValueError("provenance identifier value must not be empty")


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
class NormalizedCandidate:
    """A candidate mapping consumed by the assessor core."""

    candidate_id: str
    target_interval: GenomicInterval
    orientation: MappingOrientation
    mapping_provenance: ProvenanceSource
    evidence: tuple[EvidenceObservation, ...] = ()

    def __post_init__(self) -> None:
        if not self.candidate_id:
            raise ValueError("candidate_id must not be empty")
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
    truth claim.
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
