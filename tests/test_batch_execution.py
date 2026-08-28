from pathlib import Path

import pytest

from liftassess.batch import BatchInputRecord, BatchTargetRelationshipKind
from liftassess.batch_execution import run_indexed_chain_batch
from liftassess.chain_index import ChainIndex, build_chain_index
from liftassess.models import (
    AssemblyIdentifier,
    EvidenceAvailabilityTier,
    GenomicInterval,
    ProvenanceSource,
)
from liftassess.query_context import QueryContextNotRunReason, QueryContextState
from liftassess.resource_cache import (
    CachedResource,
    CachedUCSCChainResource,
    ucsc_resource_terms,
)
from liftassess.resource_identity import sha256_identifier_for_file

SOURCE = AssemblyIdentifier(name="hg38", provider="UCSC")
TARGET = AssemblyIdentifier(name="hg19", provider="UCSC")
ALIGNMENT = ProvenanceSource(source_id="alignment", label="shared alignment")


def _chain_text() -> str:
    return """\
chain 100 chr1 1000 + 0 10 chrA 2000 + 100 110 1
10

chain 90 chr1 1000 + 20 30 chrA 2000 + 100 110 2
10

chain 80 chr1 1000 + 40 50 chrA 2000 + 105 115 3
10

"""


def _chain_context(
    tmp_path: Path,
    *,
    chain_text: str | None = None,
) -> tuple[CachedUCSCChainResource, ChainIndex]:
    path = tmp_path / "hg38ToHg19.over.chain"
    path.write_text(chain_text or _chain_text(), encoding="ascii")
    identifier = sha256_identifier_for_file(path).value
    url = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
        "hg38ToHg19.over.chain.gz"
    )
    resource = CachedResource(
        path=path,
        source_url=url,
        retrieved_at="2026-08-28T00:00:00Z",
        sha256=identifier,
        size_bytes=path.stat().st_size,
        provider_checksum=None,
        terms=ucsc_resource_terms(url),
        cache_hit=True,
    )
    context = CachedUCSCChainResource(
        source_db="hg38",
        target_db="hg19",
        evidence_tier=EvidenceAvailabilityTier.LIFTOVER_ONLY,
        chain=resource,
    )
    index = build_chain_index(
        path,
        tmp_path / "index",
        source_chain_sha256_identifier=identifier,
        source_chain_size_bytes=path.stat().st_size,
    ).index
    return context, index


def _record(record_id: str, start: int, end: int) -> BatchInputRecord:
    return BatchInputRecord(
        record_id=record_id,
        source_interval=GenomicInterval(SOURCE, "chr1", start, end),
        source_line_number=int(record_id.removeprefix("row-")),
    )


def test_indexed_chain_batch_projects_all_rows_and_derives_relationships(
    tmp_path: Path,
) -> None:
    chain_context, index = _chain_context(tmp_path)
    records = (
        _record("row-1", 0, 10),
        _record("row-2", 20, 30),
        _record("row-3", 40, 50),
        _record("row-4", 60, 70),
    )

    result = run_indexed_chain_batch(
        records,
        chain_context,
        target_assembly=TARGET,
        alignment_provenance=ALIGNMENT,
        chain_index=index,
    )

    assert result.chain_sha256_identifier == chain_context.chain.sha256
    assert result.chain_resource == chain_context.chain
    assert result.alignment_provenance == ALIGNMENT
    assert result.evidence_tier is EvidenceAvailabilityTier.LIFTOVER_ONLY
    assert [len(item.candidates) for item in result.record_assessments] == [1, 1, 1, 0]
    kinds = [relationship.kind for relationship in result.relationships.relationships]
    assert kinds.count(BatchTargetRelationshipKind.EXACT_TARGET_COLLISION) == 1
    assert kinds.count(BatchTargetRelationshipKind.OVERLAPPING_TARGET_PROJECTIONS) == 2


def test_indexed_chain_batch_preserves_chain_provenance_on_each_candidate(
    tmp_path: Path,
) -> None:
    chain_context, index = _chain_context(tmp_path)

    result = run_indexed_chain_batch(
        (_record("row-1", 0, 10),),
        chain_context,
        target_assembly=TARGET,
        alignment_provenance=ALIGNMENT,
        chain_index=index,
    )

    candidate = result.record_assessments[0].candidates[0]
    identifiers = {
        identifier.value for identifier in candidate.mapping_provenance.identifiers
    }
    assert chain_context.chain.sha256 in identifiers


def test_indexed_chain_batch_refuses_missing_index(tmp_path: Path) -> None:
    chain_context, _index = _chain_context(tmp_path)

    with pytest.raises(ValueError, match="prepared chain index"):
        run_indexed_chain_batch(
            (_record("row-1", 0, 10),),
            chain_context,
            target_assembly=TARGET,
            alignment_provenance=ALIGNMENT,
            chain_index=None,
        )


def test_indexed_chain_batch_rejects_empty_record_set(tmp_path: Path) -> None:
    chain_context, index = _chain_context(tmp_path)

    with pytest.raises(ValueError, match="at least one"):
        run_indexed_chain_batch(
            (),
            chain_context,
            target_assembly=TARGET,
            alignment_provenance=ALIGNMENT,
            chain_index=index,
        )


def test_indexed_chain_batch_automatically_detects_neighborhood_level_collision(
    tmp_path: Path,
) -> None:
    chain_context, index = _chain_context(
        tmp_path,
        chain_text="""\
chain 100 chr1 1000 + 100 201 chrA 2000 + 500 601 11
101

chain 90 chr1 1000 + 300 401 chrA 2000 + 500 601 12
101

""",
    )
    records = (
        _record("row-1", 150, 151),
        _record("row-2", 350, 351),
    )

    result = run_indexed_chain_batch(
        records,
        chain_context,
        target_assembly=TARGET,
        alignment_provenance=ALIGNMENT,
        chain_index=index,
    )

    contexts = tuple(item.context_result for item in result.point_context_records)
    assert all(context is not None for context in contexts)
    assert [context.check_state for context in contexts if context is not None] == [
        QueryContextState.RUN,
        QueryContextState.RUN,
    ]
    assert [
        context.tested_source_interval for context in contexts if context is not None
    ] == [
        GenomicInterval(SOURCE, "chr1", 100, 201),
        GenomicInterval(SOURCE, "chr1", 300, 401),
    ]
    assert [
        relationship.kind for relationship in result.relationships.relationships
    ] == [BatchTargetRelationshipKind.EXACT_TARGET_COLLISION]
    assert [
        relationship.kind
        for relationship in result.point_context_relationships.relationships
    ] == [BatchTargetRelationshipKind.EXACT_TARGET_COLLISION]


def test_indexed_chain_batch_keeps_offset_point_context_overlap_distinct(
    tmp_path: Path,
) -> None:
    chain_context, index = _chain_context(
        tmp_path,
        chain_text="""\
chain 100 chr1 1000 + 100 201 chrA 2000 + 500 601 11
101

chain 90 chr1 1000 + 300 401 chrA 2000 + 502 603 12
101

""",
    )
    records = (
        _record("row-1", 150, 151),
        _record("row-2", 350, 351),
    )

    result = run_indexed_chain_batch(
        records,
        chain_context,
        target_assembly=TARGET,
        alignment_provenance=ALIGNMENT,
        chain_index=index,
    )

    assert result.relationships.relationships == ()
    context_relationship = result.point_context_relationships.relationships[0]
    assert (
        context_relationship.kind
        is BatchTargetRelationshipKind.OVERLAPPING_TARGET_PROJECTIONS
    )
    assert context_relationship.overlap_intervals == (
        GenomicInterval(TARGET, "chrA", 502, 601),
    )


def test_indexed_chain_batch_reports_point_context_not_run_without_indexed_bound(
    tmp_path: Path,
) -> None:
    chain_context, index = _chain_context(tmp_path)
    missing_sequence_record = BatchInputRecord(
        record_id="row-1",
        source_interval=GenomicInterval(SOURCE, "chrMissing", 10, 11),
        source_line_number=1,
    )

    result = run_indexed_chain_batch(
        (missing_sequence_record,),
        chain_context,
        target_assembly=TARGET,
        alignment_provenance=ALIGNMENT,
        chain_index=index,
    )

    assert result.record_assessments[0].candidates == ()
    context = result.point_context_records[0].context_result
    assert context is not None
    assert context.check_state is QueryContextState.NOT_RUN
    assert context.not_run_reason is QueryContextNotRunReason.SOURCE_BOUNDS_UNAVAILABLE
    assert result.point_context_relationships.relationships == ()


def test_indexed_chain_batch_does_not_widen_non_point_rows(tmp_path: Path) -> None:
    chain_context, index = _chain_context(tmp_path)

    result = run_indexed_chain_batch(
        (_record("row-1", 0, 10),),
        chain_context,
        target_assembly=TARGET,
        alignment_provenance=ALIGNMENT,
        chain_index=index,
    )

    assert result.point_context_records[0].context_result is None
    assert result.point_context_relationships.relationships == ()
