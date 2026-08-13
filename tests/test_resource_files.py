from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from liftassess import (
    AssemblyIdentifier,
    EvidenceKind,
    GenomicInterval,
    MappingOrientation,
    NormalizedCandidate,
    ProvenanceIdentifierKind,
    ProvenanceSource,
    ReciprocalBestMembershipStatus,
    ReciprocalBestMembershipSummary,
    ReciprocalBestResourceCompleteness,
    ResourceIdentityMismatchError,
    build_ucsc_candidates_from_files,
    iter_chain_file,
    iter_net_file,
    provenance_source_for_file,
)
from liftassess.models import EvidenceValue
from liftassess.net import NetClassification


@pytest.fixture
def source_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier(name="canFam3", provider="UCSC")


@pytest.fixture
def target_assembly() -> AssemblyIdentifier:
    return AssemblyIdentifier(name="canFam4", provider="UCSC")


@pytest.fixture
def provenance() -> tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource]:
    alignment = ProvenanceSource(
        source_id="ucsc-canFam3-canFam4-alignment",
        label="UCSC canFam3 to canFam4 comparative alignment",
    )
    return (
        ProvenanceSource(
            source_id="chain-file",
            label="test chain resource",
            derived_from=(alignment,),
        ),
        ProvenanceSource(
            source_id="net-file",
            label="test net resource",
            derived_from=(alignment,),
        ),
        ProvenanceSource(
            source_id="rbest-chain-file",
            label="test reciprocal-best chain resource",
            derived_from=(alignment,),
        ),
    )


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _chain_text(*, chain_id: int = 1) -> str:
    return (
        f"chain 100 chr1 1000 + 100 120 chrA 2000 + 500 520 {chain_id}\n"
        "20\n\n"
    )


def _net_text() -> str:
    return (
        "net chr1 1000\n"
        " fill 100 20 chrA + 500 20 id 1 score 100 ali 20 qDup 0 type syn\n"
    )


def _observation_value(
    candidate: NormalizedCandidate, kind: EvidenceKind
) -> EvidenceValue:
    return next(item.value for item in candidate.evidence if item.kind is kind)


def test_iter_chain_file_reads_plain_text(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain"
    _write_text(chain_path, _chain_text(chain_id=7))

    (record,) = tuple(iter_chain_file(chain_path))

    assert record.chain_id == 7
    assert record.target_name == "chr1"
    assert record.query_name == "chrA"


def test_iter_chain_and_net_files_read_gzip(tmp_path: Path) -> None:
    chain_path = tmp_path / "example.chain.gz"
    net_path = tmp_path / "example.net.gz"
    _write_gzip(chain_path, _chain_text())
    _write_gzip(net_path, _net_text())

    (chain,) = tuple(iter_chain_file(chain_path))
    (net,) = tuple(iter_net_file(net_path))

    assert chain.chain_id == 1
    assert net.chain_id == 1
    assert net.classification is NetClassification.SYNTENIC



def test_resource_compression_is_detected_from_magic_bytes_not_suffix(
    tmp_path: Path,
) -> None:
    gzip_without_suffix = tmp_path / "renamed.chain"
    plain_with_gzip_suffix = tmp_path / "plain.chain.gz"
    _write_gzip(gzip_without_suffix, _chain_text(chain_id=8))
    _write_text(plain_with_gzip_suffix, _chain_text(chain_id=9))

    (compressed_record,) = tuple(iter_chain_file(gzip_without_suffix))
    (plain_record,) = tuple(iter_chain_file(plain_with_gzip_suffix))

    assert compressed_record.chain_id == 8
    assert plain_record.chain_id == 9


def test_file_adapter_runs_parsers_and_engine_end_to_end(
    tmp_path: Path,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    chain_path = tmp_path / "canFam3.canFam4.all.chain.gz"
    net_path = tmp_path / "canFam3.canFam4.net.gz"
    rbest_path = tmp_path / "canFam3.canFam4.rbest.chain.gz"
    _write_gzip(chain_path, _chain_text(chain_id=1))
    _write_gzip(net_path, _net_text())
    _write_gzip(rbest_path, _chain_text(chain_id=101))

    alignment = ProvenanceSource(
        source_id="ucsc-canFam3-canFam4-alignment",
        label="UCSC canFam3 to canFam4 comparative alignment",
    )
    chain_provenance = provenance_source_for_file(
        chain_path, label="test chain resource", derived_from=(alignment,)
    )
    net_provenance = provenance_source_for_file(
        net_path, label="test net resource", derived_from=(alignment,)
    )
    rbest_provenance = provenance_source_for_file(
        rbest_path, label="test reciprocal-best resource", derived_from=(alignment,)
    )

    (candidate,) = build_ucsc_candidates_from_files(
        GenomicInterval(source_assembly, "chr1", 105, 115),
        chain_path,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_path=net_path,
        net_provenance=net_provenance,
        reciprocal_best_chain_path=rbest_path,
        reciprocal_best_provenance=rbest_provenance,
        reciprocal_best_completeness=(
            ReciprocalBestResourceCompleteness.COMPLETE_RESOURCE
        ),
    )

    assert candidate.orientation is MappingOrientation.SAME
    assert (candidate.target_interval.start, candidate.target_interval.end) == (
        505,
        515,
    )
    assert _observation_value(candidate, EvidenceKind.CHAIN_SCORE) == 100.0
    assert _observation_value(candidate, EvidenceKind.NET_CLASSIFICATION) == "syn"

    reciprocal = _observation_value(
        candidate, EvidenceKind.RECIPROCAL_BEST_MEMBERSHIP
    )
    assert isinstance(reciprocal, ReciprocalBestMembershipSummary)
    assert reciprocal.status is ReciprocalBestMembershipStatus.FULL
    assert reciprocal.covered_source_bases == 10
    assert reciprocal.chains_examined == 1


def test_file_adapter_preserves_content_addressed_file_provenance(
    tmp_path: Path,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    chain_path = tmp_path / "example.chain"
    _write_text(chain_path, _chain_text(chain_id=17))
    alignment = ProvenanceSource("alignment", "upstream alignment")
    chain_provenance = provenance_source_for_file(
        chain_path,
        label="local chain resource",
        derived_from=(alignment,),
    )

    (candidate,) = build_ucsc_candidates_from_files(
        GenomicInterval(source_assembly, "chr1", 105, 115),
        chain_path,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
    )

    assert candidate.mapping_provenance == chain_provenance
    assert candidate.candidate_id == f"{chain_provenance.source_id}:chain:17"
    assert chain_provenance.identifiers[0].kind is ProvenanceIdentifierKind.SHA256


def test_file_adapter_requires_sha256_file_provenance(
    tmp_path: Path,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    chain_path = tmp_path / "example.chain"
    _write_text(chain_path, _chain_text())
    chain_provenance = ProvenanceSource("chain-file", "unhashed chain file")

    with pytest.raises(ValueError, match="exactly one canonical SHA256 identifier"):
        build_ucsc_candidates_from_files(
            GenomicInterval(source_assembly, "chr1", 105, 115),
            chain_path,
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
        )


def test_file_adapter_rejects_bytes_changed_after_provenance_was_created(
    tmp_path: Path,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    chain_path = tmp_path / "example.chain"
    _write_text(chain_path, _chain_text(chain_id=17))
    chain_provenance = provenance_source_for_file(
        chain_path, label="local chain resource", derived_from=()
    )

    _write_text(chain_path, _chain_text(chain_id=18))

    with pytest.raises(ResourceIdentityMismatchError, match="provenance mismatch"):
        build_ucsc_candidates_from_files(
            GenomicInterval(source_assembly, "chr1", 105, 115),
            chain_path,
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
        )


def test_file_backed_chain_and_net_provenance_preserve_shared_upstream_source(
    tmp_path: Path,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
) -> None:
    chain_path = tmp_path / "example.chain"
    net_path = tmp_path / "example.net"
    _write_text(chain_path, _chain_text())
    _write_text(net_path, _net_text())
    alignment = ProvenanceSource("alignment", "shared alignment")
    chain_provenance = provenance_source_for_file(
        chain_path, label="chain resource", derived_from=(alignment,)
    )
    net_provenance = provenance_source_for_file(
        net_path, label="net resource", derived_from=(alignment,)
    )

    (candidate,) = build_ucsc_candidates_from_files(
        GenomicInterval(source_assembly, "chr1", 105, 115),
        chain_path,
        target_assembly=target_assembly,
        chain_provenance=chain_provenance,
        net_path=net_path,
        net_provenance=net_provenance,
    )

    net_evidence = next(
        observation
        for observation in candidate.evidence
        if observation.kind is EvidenceKind.NET_CLASSIFICATION
    )
    assert candidate.mapping_provenance.derived_from == (alignment,)
    assert net_evidence.provenance.derived_from == (net_provenance,)
    assert net_evidence.provenance.derived_from[0].derived_from == (alignment,)


def test_file_adapter_uses_engine_validation_for_optional_groups(
    tmp_path: Path,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, _ = provenance
    chain_path = tmp_path / "example.chain"
    net_path = tmp_path / "example.net"
    _write_text(chain_path, _chain_text())
    _write_text(net_path, _net_text())

    with pytest.raises(ValueError, match="net records and net provenance"):
        build_ucsc_candidates_from_files(
            GenomicInterval(source_assembly, "chr1", 105, 115),
            chain_path,
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
            net_path=net_path,
        )


def test_file_adapter_uses_engine_validation_for_reciprocal_best_group(
    tmp_path: Path,
    source_assembly: AssemblyIdentifier,
    target_assembly: AssemblyIdentifier,
    provenance: tuple[ProvenanceSource, ProvenanceSource, ProvenanceSource],
) -> None:
    chain_provenance, _, rbest_provenance = provenance
    chain_path = tmp_path / "example.chain"
    rbest_path = tmp_path / "example.rbest.chain"
    _write_text(chain_path, _chain_text())
    _write_text(rbest_path, _chain_text(chain_id=101))

    with pytest.raises(
        ValueError,
        match="reciprocal-best chains, provenance, and completeness",
    ):
        build_ucsc_candidates_from_files(
            GenomicInterval(source_assembly, "chr1", 105, 115),
            chain_path,
            target_assembly=target_assembly,
            chain_provenance=chain_provenance,
            reciprocal_best_chain_path=rbest_path,
            reciprocal_best_provenance=rbest_provenance,
        )
