from io import StringIO

import pytest

from liftassess import MappingOrientation
from liftassess.chain import ChainFormatError, ChainStrand, iter_chain_records


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


def test_parser_rejects_negative_block_gap_with_line_context() -> None:
    text = """\
chain 100 chrSource 1000 + 100 130 chrTarget 2000 + 500 530 7
10 -1 0
20

"""

    with pytest.raises(ChainFormatError, match=r"line 2: chain target gap"):
        tuple(iter_chain_records(StringIO(text)))
