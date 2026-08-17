from io import StringIO

import pytest

from liftassess import MappingOrientation
from liftassess.chain import (
    ChainFormatError,
    ChainStrand,
    chain_candidate_id,
    chain_id_from_candidate_id,
    iter_chain_records,
)


def test_chain_candidate_id_round_trips_through_shared_decoder() -> None:
    candidate_id = chain_candidate_id("file:sha256:abc", 42)

    assert candidate_id == "file:sha256:abc:chain:42"
    assert chain_id_from_candidate_id(candidate_id) == 42
    assert chain_id_from_candidate_id("not-a-chain-candidate") is None
    assert chain_id_from_candidate_id(":chain:42") is None
    assert chain_id_from_candidate_id("source:chain:-1") is None
    assert chain_id_from_candidate_id("source:chain:+1") is None
    assert chain_id_from_candidate_id("source:chain:01") is None
    assert chain_id_from_candidate_id("source:chain:not-an-int") is None


def test_parser_reads_chain_header_blocks_and_gaps() -> None:
    text = """\
chain 100 chrSource 1000 + 100 135 chrTarget 2000 + 500 538 7
10 5 8
20

"""

    (record,) = tuple(iter_chain_records(StringIO(text)))

    assert record.score == 100
    assert record.target_name == "chrSource"
    assert record.target_start == 100
    assert record.target_end == 135
    assert record.query_name == "chrTarget"
    assert record.query_start == 500
    assert record.query_end == 538
    assert record.chain_id == 7
    assert len(record.blocks) == 2
    assert record.blocks[0].size == 10
    assert record.blocks[0].target_gap == 5
    assert record.blocks[0].query_gap == 8
    assert record.blocks[1].size == 20
    assert record.blocks[1].is_terminal


def test_negative_query_coordinates_convert_to_forward_reference() -> None:
    text = """\
chain 50 chrSource 1000 + 100 120 chrTarget 2000 - 300 320 8
20

"""

    (record,) = tuple(iter_chain_records(StringIO(text)))

    assert record.query_strand is ChainStrand.MINUS
    assert record.query_forward_start == 1680
    assert record.query_forward_end == 1700
    assert record.orientation is MappingOrientation.REVERSE


def test_positive_query_strand_has_same_orientation() -> None:
    text = """\
chain 50 chrSource 1000 + 100 120 chrTarget 2000 + 300 320 8
20

"""

    (record,) = tuple(iter_chain_records(StringIO(text)))

    assert record.orientation is MappingOrientation.SAME


def test_parser_rejects_negative_target_strand() -> None:
    text = """\
chain 50 chrSource 1000 - 100 120 chrTarget 2000 - 300 320 8
20

"""

    with pytest.raises(ChainFormatError, match="target strand must be"):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_streams_multiple_records() -> None:
    text = """\
chain 10 chr1 100 + 0 5 chrA 100 + 10 15 1
5

chain 20 chr2 100 + 20 25 chrB 100 + 30 35 2
5

"""

    records = tuple(iter_chain_records(StringIO(text)))

    assert [record.chain_id for record in records] == [1, 2]


def test_parser_rejects_blocks_that_do_not_match_header_span() -> None:
    text = """\
chain 100 chrSource 1000 + 100 136 chrTarget 2000 + 500 538 7
10 5 8
20

"""

    with pytest.raises(ChainFormatError, match="target header interval"):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_rejects_blank_line_before_terminal_block() -> None:
    text = """\
chain 100 chrSource 1000 + 100 115 chrTarget 2000 + 500 518 7
10 5 8

"""

    with pytest.raises(ChainFormatError, match="before terminal chain block"):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_rejects_non_chain_header() -> None:
    with pytest.raises(ChainFormatError, match="expected chain header"):
        tuple(iter_chain_records(StringIO("not-chain data\n")))


def test_parser_rejects_non_chain_text_between_records() -> None:
    text = """\
chain 10 chr1 100 + 0 5 chrA 100 + 10 15 1
5

not-chain data
"""

    with pytest.raises(ChainFormatError, match=r"line 4: expected chain header"):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_rejects_negative_block_gap_with_line_context() -> None:
    text = """\
chain 100 chrSource 1000 + 100 130 chrTarget 2000 + 500 530 7
10 -1 0
20

"""

    with pytest.raises(ChainFormatError, match=r"line 2: chain target gap"):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_rejects_eof_before_terminal_block() -> None:
    text = "chain 100 chrSource 1000 + 100 120 chrTarget 2000 + 500 520 7\n10 5 8\n"

    with pytest.raises(ChainFormatError, match="ended before terminal block"):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_rejects_malformed_block_field_count() -> None:
    text = "chain 100 chrSource 1000 + 100 120 chrTarget 2000 + 500 520 7\n10 5\n20\n\n"

    with pytest.raises(ChainFormatError, match="must contain 1 or 3 fields"):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_accepts_terminal_block_at_eof_without_blank_line() -> None:
    text = "chain 100 chrSource 1000 + 100 120 chrTarget 2000 + 500 520 7\n20"

    (record,) = tuple(iter_chain_records(StringIO(text)))

    assert record.chain_id == 7
    assert len(record.blocks) == 1
    assert record.blocks[0].size == 20
    assert record.blocks[0].is_terminal


def test_parser_skips_ucsc_metadata_while_looking_for_chain_headers() -> None:
    text = """\
##matrix=axtChain 16 91,-114,-31,-123
# lastz.v1.04.03 H=2000 M=254 --format=axt+
#
chain 10 chr1 100 + 0 5 chrA 100 + 10 15 1
5

##gapPenalties=axtChain O=400 E=30
chain 20 chr2 100 + 20 25 chrB 100 + 30 35 2
5

"""

    records = tuple(iter_chain_records(StringIO(text)))

    assert [record.chain_id for record in records] == [1, 2]


def test_parser_rejects_metadata_inside_chain_record() -> None:
    text = """\
# allowed before the record
chain 100 chrSource 1000 + 100 120 chrTarget 2000 + 500 520 7
# metadata is not valid inside the block list
20

"""

    with pytest.raises(
        ChainFormatError,
        match=r"line 3: metadata/comment line inside chain record",
    ):
        tuple(iter_chain_records(StringIO(text)))


def test_parser_skips_trailing_metadata_after_last_record() -> None:
    text = """\
chain 10 chr1 100 + 0 5 chrA 100 + 10 15 1
5

# trailing metadata
##gapPenalties=axtChain O=400 E=30
"""

    records = tuple(iter_chain_records(StringIO(text)))

    assert [record.chain_id for record in records] == [1]


def test_aligned_block_iterator_normalizes_reverse_query_geometry_once() -> None:
    text = """\
chain 50 chrSource 1000 + 100 125 chrTarget 2000 - 300 328 8
10 5 8
10

"""

    (record,) = tuple(iter_chain_records(StringIO(text)))
    blocks = tuple(record.iter_aligned_blocks())

    assert [(block.target_start, block.target_end) for block in blocks] == [
        (100, 110),
        (115, 125),
    ]
    assert [
        (block.query_forward_start, block.query_forward_end) for block in blocks
    ] == [(1690, 1700), (1672, 1682)]
    assert blocks[1].query_interval_for_target_interval(117, 121) == (1676, 1680)
