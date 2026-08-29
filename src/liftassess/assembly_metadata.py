"""Assembly-sequence metadata and source-interval preflight.

This module deliberately does not infer assembly membership from chain resources.
Callers provide verified assembly metadata, and preflight answers only whether a
submitted source interval uses a canonical assembly sequence name and lies within the
metadata-defined sequence bounds.  Exact assembly-aware aliases may be retained as
suggestions, but they are not silently substituted for the submitted sequence name.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType

from .models import AssemblyIdentifier, GenomicInterval, ProvenanceSource


class SourceIntervalPreflightState(str, Enum):
    """Authoritative source-input states established before scientific assessment."""

    VALID = "VALID"
    UNRECOGNIZED_SOURCE_SEQUENCE_NAME = "UNRECOGNIZED_SOURCE_SEQUENCE_NAME"
    INVALID_SOURCE_COORDINATE = "INVALID_SOURCE_COORDINATE"


@dataclass(frozen=True)
class AssemblySequenceAlias:
    """One exact assembly-aware alias for a canonical sequence name."""

    name: str
    sources: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("assembly sequence alias must not be empty")
        if any(not source for source in self.sources):
            raise ValueError(
                "assembly sequence alias sources must not contain empty values"
            )
        if len(set(self.sources)) != len(self.sources):
            raise ValueError("assembly sequence alias sources must be unique")


@dataclass(frozen=True)
class AssemblySequenceMetadata:
    """Canonical sequence identity and authoritative assembly length."""

    sequence_name: str
    length: int
    aliases: tuple[AssemblySequenceAlias, ...] = ()

    def __post_init__(self) -> None:
        if not self.sequence_name:
            raise ValueError("assembly sequence name must not be empty")
        if self.length <= 0:
            raise ValueError("assembly sequence length must be positive")

        alias_names = tuple(alias.name for alias in self.aliases)
        if len(set(alias_names)) != len(alias_names):
            raise ValueError("assembly sequence aliases must be unique")
        if self.sequence_name in alias_names:
            raise ValueError("canonical sequence name must not be repeated as an alias")


@dataclass(frozen=True)
class AssemblySequenceCatalog:
    """One assembly's canonical sequences, bounds, exact aliases, and provenance."""

    assembly: AssemblyIdentifier
    sequences: tuple[AssemblySequenceMetadata, ...]
    sequence_provenance: ProvenanceSource
    alias_provenance: ProvenanceSource | None = None
    _sequences_by_name: Mapping[str, AssemblySequenceMetadata] = field(
        init=False,
        repr=False,
        compare=False,
    )
    _aliases_by_name: Mapping[str, tuple[str, AssemblySequenceAlias]] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not self.sequences:
            raise ValueError(
                "assembly sequence catalog must contain at least one sequence"
            )

        sequences_by_name: dict[str, AssemblySequenceMetadata] = {}
        for sequence in self.sequences:
            if sequence.sequence_name in sequences_by_name:
                raise ValueError(
                    "assembly sequence catalog contains duplicate canonical sequence "
                    f"name: {sequence.sequence_name}"
                )
            sequences_by_name[sequence.sequence_name] = sequence

        aliases_by_name: dict[str, tuple[str, AssemblySequenceAlias]] = {}
        for sequence in self.sequences:
            for alias in sequence.aliases:
                if alias.name in sequences_by_name:
                    raise ValueError(
                        "assembly sequence alias conflicts with a canonical sequence "
                        f"name: {alias.name}"
                    )
                previous = aliases_by_name.get(alias.name)
                if previous is not None:
                    raise ValueError(
                        "assembly sequence alias maps to more than one canonical "
                        f"sequence: {alias.name}"
                    )
                aliases_by_name[alias.name] = (sequence.sequence_name, alias)

        if aliases_by_name and self.alias_provenance is None:
            raise ValueError("assembly sequence aliases require alias provenance")

        object.__setattr__(
            self,
            "_sequences_by_name",
            MappingProxyType(sequences_by_name),
        )
        object.__setattr__(
            self,
            "_aliases_by_name",
            MappingProxyType(aliases_by_name),
        )

    def sequence(self, sequence_name: str) -> AssemblySequenceMetadata | None:
        """Return canonical metadata for an exact assembly sequence name."""

        return self._sequences_by_name.get(sequence_name)

    def alias_match(self, alias_name: str) -> tuple[str, AssemblySequenceAlias] | None:
        """Return the canonical sequence and metadata for an exact unique alias."""

        return self._aliases_by_name.get(alias_name)


@dataclass(frozen=True)
class SourceIntervalPreflightResult:
    """Authoritative source-input facts established before mapping is attempted."""

    state: SourceIntervalPreflightState
    source_interval: GenomicInterval
    canonical_sequence_name: str | None
    sequence_length: int | None
    suggested_sequence_name: str | None = None
    alias_sources: tuple[str, ...] = ()
    provenance_sources: tuple[ProvenanceSource, ...] = ()

    def __post_init__(self) -> None:
        if not self.provenance_sources:
            raise ValueError(
                "source preflight requires authoritative metadata provenance"
            )
        source_ids = tuple(source.source_id for source in self.provenance_sources)
        if len(set(source_ids)) != len(source_ids):
            raise ValueError("source preflight provenance sources must be unique")
        if self.sequence_length is not None and self.sequence_length <= 0:
            raise ValueError("preflight sequence length must be positive")
        if self.state is SourceIntervalPreflightState.UNRECOGNIZED_SOURCE_SEQUENCE_NAME:
            if (
                self.canonical_sequence_name is not None
                or self.sequence_length is not None
            ):
                raise ValueError(
                    "unrecognized source sequence cannot carry canonical bounds"
                )
            if self.suggested_sequence_name is None and self.alias_sources:
                raise ValueError(
                    "alias sources require a suggested canonical sequence name"
                )
            return

        if self.canonical_sequence_name is None or self.sequence_length is None:
            raise ValueError(
                "recognized source sequence preflight requires canonical name "
                "and bounds"
            )
        if self.suggested_sequence_name is not None or self.alias_sources:
            raise ValueError(
                "recognized source sequence preflight cannot carry an alias suggestion"
            )

    @property
    def mapping_may_proceed(self) -> bool:
        """Return whether scientific mapping may begin for this source interval."""

        return self.state is SourceIntervalPreflightState.VALID


def preflight_source_interval(
    source_interval: GenomicInterval,
    catalog: AssemblySequenceCatalog,
) -> SourceIntervalPreflightResult:
    """Validate a source interval against assembly metadata, never chain membership."""

    if source_interval.assembly != catalog.assembly:
        raise ValueError(
            "source interval assembly does not match assembly sequence catalog"
        )
    if source_interval.length <= 0:
        raise ValueError("source interval preflight requires a non-empty interval")

    canonical_provenance = (catalog.sequence_provenance,)
    alias_lookup_provenance = canonical_provenance + (
        (catalog.alias_provenance,) if catalog.alias_provenance is not None else ()
    )

    sequence = catalog.sequence(source_interval.sequence_name)
    if sequence is None:
        alias_match = catalog.alias_match(source_interval.sequence_name)
        if alias_match is None:
            return SourceIntervalPreflightResult(
                state=(SourceIntervalPreflightState.UNRECOGNIZED_SOURCE_SEQUENCE_NAME),
                source_interval=source_interval,
                canonical_sequence_name=None,
                sequence_length=None,
                provenance_sources=alias_lookup_provenance,
            )

        canonical_name, alias = alias_match
        return SourceIntervalPreflightResult(
            state=SourceIntervalPreflightState.UNRECOGNIZED_SOURCE_SEQUENCE_NAME,
            source_interval=source_interval,
            canonical_sequence_name=None,
            sequence_length=None,
            suggested_sequence_name=canonical_name,
            alias_sources=alias.sources,
            provenance_sources=alias_lookup_provenance,
        )

    if source_interval.end > sequence.length:
        return SourceIntervalPreflightResult(
            state=SourceIntervalPreflightState.INVALID_SOURCE_COORDINATE,
            source_interval=source_interval,
            canonical_sequence_name=sequence.sequence_name,
            sequence_length=sequence.length,
            provenance_sources=canonical_provenance,
        )

    return SourceIntervalPreflightResult(
        state=SourceIntervalPreflightState.VALID,
        source_interval=source_interval,
        canonical_sequence_name=sequence.sequence_name,
        sequence_length=sequence.length,
        provenance_sources=canonical_provenance,
    )


def parse_ucsc_chrom_info(
    lines: Iterable[str],
) -> tuple[AssemblySequenceMetadata, ...]:
    """Parse UCSC ``chromInfo`` table rows into canonical sequence metadata.

    UCSC's published table schema is ``chrom, size, fileName``.  The file path is
    deliberately ignored here because source validity and bounds depend only on the
    canonical browser sequence name and its assembly length.
    """

    sequences: list[AssemblySequenceMetadata] = []
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(
                "UCSC chromInfo row must contain exactly 3 tab-separated fields "
                f"at line {line_number}"
            )
        sequence_name, length_text, _file_name = fields
        if not sequence_name:
            raise ValueError(
                f"UCSC chromInfo sequence name is empty at line {line_number}"
            )
        try:
            length = int(length_text)
        except ValueError as exc:
            raise ValueError(
                "UCSC chromInfo sequence length is not an integer at line "
                f"{line_number}"
            ) from exc
        sequences.append(
            AssemblySequenceMetadata(sequence_name=sequence_name, length=length)
        )

    if not sequences:
        raise ValueError("UCSC chromInfo data contains no sequence rows")
    return tuple(sequences)


def build_ucsc_assembly_sequence_catalog(
    assembly: AssemblyIdentifier,
    chrom_info_lines: Iterable[str],
    *,
    sequence_provenance: ProvenanceSource,
    chrom_alias_lines: Iterable[str] | None = None,
    alias_provenance: ProvenanceSource | None = None,
) -> AssemblySequenceCatalog:
    """Build one UCSC assembly catalog from verified database table rows.

    ``chromAlias`` rows use UCSC's published ``alias, chrom, source`` schema. Alias
    names are retained only when their canonical ``chrom`` exists in the supplied
    ``chromInfo`` metadata. They are suggestions for user input; this function does
    not redefine the assembly's canonical sequence namespace.
    """

    sequences = parse_ucsc_chrom_info(chrom_info_lines)
    if chrom_alias_lines is None:
        if alias_provenance is not None:
            raise ValueError("alias provenance requires UCSC chromAlias data")
        return AssemblySequenceCatalog(
            assembly=assembly,
            sequences=sequences,
            sequence_provenance=sequence_provenance,
        )
    if alias_provenance is None:
        raise ValueError("UCSC chromAlias data requires alias provenance")

    aliases_by_sequence: dict[str, list[AssemblySequenceAlias]] = {
        sequence.sequence_name: [] for sequence in sequences
    }
    seen_aliases: set[str] = set()
    for line_number, raw_line in enumerate(chrom_alias_lines, start=1):
        line = raw_line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(
                "UCSC chromAlias row must contain exactly 3 tab-separated fields "
                f"at line {line_number}"
            )
        alias_name, canonical_name, source_text = fields
        if not alias_name or not canonical_name:
            raise ValueError(
                "UCSC chromAlias alias and canonical sequence names must not be empty "
                f"at line {line_number}"
            )
        if canonical_name not in aliases_by_sequence:
            raise ValueError(
                "UCSC chromAlias row references a sequence absent from chromInfo at "
                f"line {line_number}: {canonical_name}"
            )
        if alias_name in seen_aliases:
            raise ValueError(
                f"UCSC chromAlias contains duplicate alias at line {line_number}: "
                f"{alias_name}"
            )
        seen_aliases.add(alias_name)
        sources = tuple(source for source in source_text.split(",") if source)
        aliases_by_sequence[canonical_name].append(
            AssemblySequenceAlias(name=alias_name, sources=sources)
        )

    sequences_with_aliases = tuple(
        AssemblySequenceMetadata(
            sequence_name=sequence.sequence_name,
            length=sequence.length,
            aliases=tuple(aliases_by_sequence[sequence.sequence_name]),
        )
        for sequence in sequences
    )
    return AssemblySequenceCatalog(
        assembly=assembly,
        sequences=sequences_with_aliases,
        sequence_provenance=sequence_provenance,
        alias_provenance=alias_provenance,
    )
