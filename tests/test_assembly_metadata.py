from __future__ import annotations

import pytest

from liftassess import (
    AssemblyIdentifier,
    GenomicInterval,
    ProvenanceSource,
)
from liftassess.assembly_metadata import (
    AssemblySequenceAlias,
    AssemblySequenceCatalog,
    AssemblySequenceMetadata,
    SourceIntervalPreflightState,
    build_ucsc_assembly_sequence_catalog,
    parse_ucsc_chrom_info,
    preflight_source_interval,
)

ASSEMBLY = AssemblyIdentifier("testDb", "UCSC")
CHROM_INFO = ProvenanceSource("chrom-info", "UCSC testDb chromInfo table")
CHROM_ALIAS = ProvenanceSource("chrom-alias", "UCSC testDb chromAlias table")


def _catalog() -> AssemblySequenceCatalog:
    return build_ucsc_assembly_sequence_catalog(
        ASSEMBLY,
        (
            "chr1\t1000\t/gbdb/testDb/testDb.2bit\n",
            "chrNoChain\t200\t/gbdb/testDb/testDb.2bit\n",
        ),
        sequence_provenance=CHROM_INFO,
        chrom_alias_lines=(
            "1\tchr1\tassembly\n",
            "NC_000001.11\tchr1\trefseq\n",
        ),
        alias_provenance=CHROM_ALIAS,
    )


def test_valid_source_interval_uses_authoritative_sequence_bound() -> None:
    interval = GenomicInterval(ASSEMBLY, "chr1", 999, 1000)

    result = preflight_source_interval(interval, _catalog())

    assert result.state is SourceIntervalPreflightState.VALID
    assert result.mapping_may_proceed
    assert result.canonical_sequence_name == "chr1"
    assert result.sequence_length == 1000
    assert result.suggested_sequence_name is None


def test_valid_assembly_sequence_does_not_require_chain_membership() -> None:
    interval = GenomicInterval(ASSEMBLY, "chrNoChain", 10, 20)

    result = preflight_source_interval(interval, _catalog())

    assert result.state is SourceIntervalPreflightState.VALID
    assert result.mapping_may_proceed
    assert result.canonical_sequence_name == "chrNoChain"
    assert result.sequence_length == 200


def test_exact_alias_is_rejected_with_authoritative_suggestion() -> None:
    interval = GenomicInterval(ASSEMBLY, "NC_000001.11", 10, 20)

    result = preflight_source_interval(interval, _catalog())

    assert (
        result.state is SourceIntervalPreflightState.UNRECOGNIZED_SOURCE_SEQUENCE_NAME
    )
    assert not result.mapping_may_proceed
    assert result.canonical_sequence_name is None
    assert result.sequence_length is None
    assert result.suggested_sequence_name == "chr1"
    assert result.alias_sources == ("refseq",)


def test_unknown_sequence_has_no_invented_alias_suggestion() -> None:
    interval = GenomicInterval(ASSEMBLY, "chrMissing", 10, 20)

    result = preflight_source_interval(interval, _catalog())

    assert (
        result.state is SourceIntervalPreflightState.UNRECOGNIZED_SOURCE_SEQUENCE_NAME
    )
    assert result.suggested_sequence_name is None
    assert result.alias_sources == ()


def test_out_of_bounds_interval_reports_authoritative_maximum() -> None:
    interval = GenomicInterval(ASSEMBLY, "chr1", 999, 1001)

    result = preflight_source_interval(interval, _catalog())

    assert result.state is SourceIntervalPreflightState.INVALID_SOURCE_COORDINATE
    assert not result.mapping_may_proceed
    assert result.canonical_sequence_name == "chr1"
    assert result.sequence_length == 1000


def test_preflight_rejects_catalog_for_another_assembly() -> None:
    interval = GenomicInterval(AssemblyIdentifier("otherDb", "UCSC"), "chr1", 0, 1)

    with pytest.raises(ValueError, match="does not match assembly sequence catalog"):
        preflight_source_interval(interval, _catalog())


def test_parse_ucsc_chrom_info_requires_published_three_field_shape() -> None:
    with pytest.raises(ValueError, match="exactly 3 tab-separated fields"):
        parse_ucsc_chrom_info(("chr1\t1000\n",))


def test_parse_ucsc_chrom_info_rejects_nonpositive_length() -> None:
    with pytest.raises(ValueError, match="length must be positive"):
        parse_ucsc_chrom_info(("chr1\t0\t/gbdb/testDb/testDb.2bit\n",))


def test_ucsc_catalog_rejects_alias_for_sequence_absent_from_chrom_info() -> None:
    with pytest.raises(ValueError, match="absent from chromInfo"):
        build_ucsc_assembly_sequence_catalog(
            ASSEMBLY,
            ("chr1\t1000\t/gbdb/testDb/testDb.2bit\n",),
            sequence_provenance=CHROM_INFO,
            chrom_alias_lines=("2\tchr2\tassembly\n",),
            alias_provenance=CHROM_ALIAS,
        )


def test_catalog_rejects_alias_collision_with_canonical_sequence() -> None:
    with pytest.raises(ValueError, match="conflicts with a canonical sequence"):
        AssemblySequenceCatalog(
            assembly=ASSEMBLY,
            sequences=(
                AssemblySequenceMetadata(
                    "chr1",
                    1000,
                    aliases=(AssemblySequenceAlias("chr2", ("assembly",)),),
                ),
                AssemblySequenceMetadata("chr2", 500),
            ),
            sequence_provenance=CHROM_INFO,
            alias_provenance=CHROM_ALIAS,
        )
