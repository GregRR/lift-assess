from __future__ import annotations

import json

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
    AssemblySequenceRoleContext,
    SourceIntervalPreflightState,
    attach_ncbi_sequence_role_context,
    build_ucsc_assembly_sequence_catalog,
    parse_ncbi_genome_sequence_report,
    parse_ucsc_assembly_description_accession,
    parse_ucsc_chrom_info,
    preflight_source_interval,
)

ASSEMBLY = AssemblyIdentifier("testDb", "UCSC")
CHROM_INFO = ProvenanceSource("chrom-info", "UCSC testDb chromInfo table")
CHROM_ALIAS = ProvenanceSource("chrom-alias", "UCSC testDb chromAlias table")
ROLE_REPORT = ProvenanceSource("sequence-report", "NCBI genome sequence report")


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


def test_parse_ucsc_description_accession_supports_accession_id_label() -> None:
    html = "<html><body><b>Accession ID:</b> GCA_000002285.2</body></html>"

    assert parse_ucsc_assembly_description_accession(html) == "GCA_000002285.2"


def test_parse_ucsc_assembly_description_accession_supports_refseq_accession() -> None:
    html = "<html><body><b>Accession ID:</b> GCF_000002285.5</body></html>"

    assert parse_ucsc_assembly_description_accession(html) == "GCF_000002285.5"


def test_parse_ucsc_description_accession_supports_assembly_accession_label() -> None:
    html = (
        "<html><body>Assembly accession: "
        '<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.29/">'
        "GCA_000001405.29</a></body></html>"
    )

    assert parse_ucsc_assembly_description_accession(html) == "GCA_000001405.29"


def test_parse_ncbi_sequence_report_preserves_provider_native_role_and_unit() -> None:
    rows = (
        json.dumps(
            {
                "assemblyAccession": "GCA_000001405.29",
                "assemblyUnit": "Primary Assembly",
                "genbankAccession": "KI270752.1",
                "length": 27745,
                "role": "unplaced-scaffold",
                "sequenceName": "KI270752.1",
                "ucscStyleName": "chrUn_KI270752v1",
            }
        ),
    )

    contexts = parse_ncbi_genome_sequence_report(
        rows, expected_assembly_accession="GCA_000001405.29"
    )

    assert contexts == (
        AssemblySequenceRoleContext(
            assembly_accession="GCA_000001405.29",
            assembly_unit="Primary Assembly",
            provider_role="unplaced-scaffold",
            length=27745,
            sequence_name="KI270752.1",
            ucsc_style_name="chrUn_KI270752v1",
            genbank_accession="KI270752.1",
        ),
    )


def test_parse_ncbi_sequence_report_rejects_non_object_row() -> None:
    rows = (json.dumps(["not", "an", "object"]),)

    with pytest.raises(TypeError, match="must be a JSON object"):
        parse_ncbi_genome_sequence_report(
            rows, expected_assembly_accession="GCA_000001405.29"
        )


def test_parse_ncbi_sequence_report_rejects_assembly_version_mismatch() -> None:
    rows = (
        json.dumps(
            {
                "assemblyAccession": "GCA_000001405.28",
                "assemblyUnit": "Primary Assembly",
                "length": 1000,
                "role": "assembled-molecule",
                "ucscStyleName": "chr1",
            }
        ),
    )

    with pytest.raises(ValueError, match="assembly accession mismatch"):
        parse_ncbi_genome_sequence_report(
            rows, expected_assembly_accession="GCA_000001405.29"
        )


def test_attach_ncbi_role_context_by_exact_ucsc_style_name() -> None:
    catalog = _catalog()
    rows = (
        json.dumps(
            {
                "assemblyAccession": "GCA_test.1",
                "assemblyUnit": "Primary Assembly",
                "chrName": "1",
                "genbankAccession": "CM000001.1",
                "length": 1000,
                "refseqAccession": "NC_000001.11",
                "role": "assembled-molecule",
                "sequenceName": "1",
                "ucscStyleName": "chr1",
            }
        ),
    )

    enriched = attach_ncbi_sequence_role_context(
        catalog,
        rows,
        expected_assembly_accession="GCA_test.1",
        role_provenance=ROLE_REPORT,
    )

    sequence = enriched.sequence("chr1")
    assert sequence is not None
    assert sequence.role_context is not None
    assert sequence.role_context.provider_role == "assembled-molecule"
    assert sequence.role_context.assembly_unit == "Primary Assembly"
    assert enriched.role_provenance == ROLE_REPORT


def test_attach_ncbi_role_context_can_join_through_verified_sequence_alias() -> None:
    catalog = _catalog()
    rows = (
        json.dumps(
            {
                "assemblyAccession": "GCA_test.1",
                "assemblyUnit": "Primary Assembly",
                "length": 1000,
                "refseqAccession": "NC_000001.11",
                "role": "assembled-molecule",
            }
        ),
    )

    enriched = attach_ncbi_sequence_role_context(
        catalog,
        rows,
        expected_assembly_accession="GCA_test.1",
        role_provenance=ROLE_REPORT,
    )

    sequence = enriched.sequence("chr1")
    assert sequence is not None
    assert sequence.role_context is not None
    assert sequence.role_context.refseq_accession == "NC_000001.11"


def test_attach_ncbi_role_context_rejects_sequence_length_mismatch() -> None:
    rows = (
        json.dumps(
            {
                "assemblyAccession": "GCA_test.1",
                "assemblyUnit": "Primary Assembly",
                "length": 999,
                "role": "assembled-molecule",
                "ucscStyleName": "chr1",
            }
        ),
    )

    with pytest.raises(ValueError, match="does not match authoritative UCSC chromInfo"):
        attach_ncbi_sequence_role_context(
            _catalog(),
            rows,
            expected_assembly_accession="GCA_test.1",
            role_provenance=ROLE_REPORT,
        )


def test_attach_ncbi_role_context_rejects_conflicting_exact_identifiers() -> None:
    catalog = build_ucsc_assembly_sequence_catalog(
        ASSEMBLY,
        (
            "chr1\t1000\t/gbdb/testDb/testDb.2bit\n",
            "chr2\t1000\t/gbdb/testDb/testDb.2bit\n",
        ),
        sequence_provenance=CHROM_INFO,
        chrom_alias_lines=("NC_000001.11\tchr2\trefseq\n",),
        alias_provenance=CHROM_ALIAS,
    )
    rows = (
        json.dumps(
            {
                "assemblyAccession": "GCA_test.1",
                "assemblyUnit": "Primary Assembly",
                "length": 1000,
                "refseqAccession": "NC_000001.11",
                "role": "assembled-molecule",
                "ucscStyleName": "chr1",
            }
        ),
    )

    with pytest.raises(ValueError, match="conflicting UCSC canonical names"):
        attach_ncbi_sequence_role_context(
            catalog,
            rows,
            expected_assembly_accession="GCA_test.1",
            role_provenance=ROLE_REPORT,
        )


def test_attach_ncbi_role_context_rejects_zero_sequence_matches() -> None:
    rows = (
        json.dumps(
            {
                "assemblyAccession": "GCA_test.1",
                "assemblyUnit": "Primary Assembly",
                "length": 1000,
                "role": "assembled-molecule",
                "ucscStyleName": "chrMissing",
            }
        ),
    )

    with pytest.raises(ValueError, match="did not resolve any sequence"):
        attach_ncbi_sequence_role_context(
            _catalog(),
            rows,
            expected_assembly_accession="GCA_test.1",
            role_provenance=ROLE_REPORT,
        )
