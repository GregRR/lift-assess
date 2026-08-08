from io import StringIO

import pytest

from liftassess import MappingOrientation
from liftassess.net import (
    NetClassification,
    NetFormatError,
    NetRecordKind,
    iter_net_records,
)


def test_parser_reads_documented_fill_gap_hierarchy_and_v1_fields() -> None:
    text = """\
net chr2L 23011544
 fill 6004 3278 chrXR_group3a - 1396397 2164 id 25606 score 23114 ali 782 qDup 576 type top tN 0
  gap 6096 1485 chrXR_group3a - 1397572 897 tN 0 qN 0
   fill 6096 513 chrU - 5570675 533 id 48675 score 4435 ali 465 qDup 533 type nonSyn
"""

    records = tuple(iter_net_records(StringIO(text)))

    assert [record.depth for record in records] == [1, 2, 3]
    top, gap, nested = records

    assert top.kind is NetRecordKind.FILL
    assert top.target_name == "chr2L"
    assert top.target_sequence_size == 23011544
    assert top.target_start == 6004
    assert top.target_end == 9282
    assert top.query_name == "chrXR_group3a"
    assert top.orientation is MappingOrientation.REVERSE
    assert top.query_start == 1396397
    assert top.query_end == 1398561
    assert top.chain_id == 25606
    assert top.score == 23114
    assert top.aligned_bases == 782
    assert top.duplicated_query_bases == 576
    assert top.classification is NetClassification.TOP
    assert ("tN", "0") in top.attributes

    assert gap.kind is NetRecordKind.GAP
    assert gap.chain_id is None
    assert gap.classification is None

    assert nested.kind is NetRecordKind.FILL
    assert nested.chain_id == 48675
    assert nested.classification is NetClassification.NON_SYNTENIC


def test_parser_streams_multiple_net_sections() -> None:
    text = """\
net chr1 100
 fill 0 10 chrA + 20 10 id 1
net chr2 200
 fill 50 20 chrB - 70 20 id 2
"""

    records = tuple(iter_net_records(StringIO(text)))

    assert [(record.target_name, record.chain_id) for record in records] == [
        ("chr1", 1),
        ("chr2", 2),
    ]
    assert [record.depth for record in records] == [1, 1]


def test_parser_preserves_unknown_optional_fields_opaquely() -> None:
    text = """\
net chr1 100
 fill 0 10 chrA + 20 10 id 1 futureField someValue
"""

    (record,) = tuple(iter_net_records(StringIO(text)))

    assert record.attributes == (("id", "1"), ("futureField", "someValue"))


def test_parser_ignores_blank_and_comment_lines() -> None:
    text = """\
# metadata comment

net chr1 100
 fill 0 10 chrA + 20 10 id 1

"""

    (record,) = tuple(iter_net_records(StringIO(text)))

    assert record.chain_id == 1


def test_parser_rejects_record_before_section_header() -> None:
    with pytest.raises(NetFormatError, match="expected net section header"):
        tuple(iter_net_records(StringIO(" fill 0 10 chrA + 20 10\n")))


def test_parser_rejects_hierarchy_depth_jump() -> None:
    text = """\
net chr1 100
 fill 0 10 chrA + 20 10 id 1
   gap 2 2 chrA + 22 0
"""

    with pytest.raises(NetFormatError, match="increase by only one level"):
        tuple(iter_net_records(StringIO(text)))


def test_parser_rejects_tab_indentation() -> None:
    text = "net chr1 100\n\tfill 0 10 chrA + 20 10 id 1\n"

    with pytest.raises(NetFormatError, match="indentation must use spaces"):
        tuple(iter_net_records(StringIO(text)))


def test_parser_rejects_unpaired_optional_field() -> None:
    text = """\
net chr1 100
 fill 0 10 chrA + 20 10 id
"""

    with pytest.raises(NetFormatError, match="name/value pairs"):
        tuple(iter_net_records(StringIO(text)))


def test_parser_rejects_duplicate_optional_field_names() -> None:
    text = """\
net chr1 100
 fill 0 10 chrA + 20 10 id 1 id 2
"""

    with pytest.raises(NetFormatError, match="duplicate net optional field"):
        tuple(iter_net_records(StringIO(text)))


def test_parser_rejects_unknown_net_classification() -> None:
    text = """\
net chr1 100
 fill 0 10 chrA + 20 10 id 1 type unexpected
"""

    with pytest.raises(NetFormatError, match="top, syn, inv, or nonSyn"):
        tuple(iter_net_records(StringIO(text)))


def test_parser_rejects_target_span_past_sequence_bounds() -> None:
    text = """\
net chr1 100
 fill 95 10 chrA + 20 10 id 1
"""

    with pytest.raises(NetFormatError, match="exceeds target sequence bounds"):
        tuple(iter_net_records(StringIO(text)))


def test_parser_allows_same_chain_id_on_multiple_fills() -> None:
    # UCSC's published net-format example uses chain 48675 in more than one fill.
    text = """\
net chr2L 23011544
 fill 6004 3278 chrXR_group3a - 1396397 2164 id 25606 type top
  gap 6096 1485 chrXR_group3a - 1397572 897
   fill 6096 513 chrU - 5570675 533 id 48675 type nonSyn
 fill 9833 1251 chrU - 5562980 1239 id 48675 type top
"""

    records = tuple(iter_net_records(StringIO(text)))
    repeated = [record for record in records if record.chain_id == 48675]

    assert len(repeated) == 2
    assert [record.depth for record in repeated] == [3, 1]
    assert [record.classification for record in repeated] == [
        NetClassification.NON_SYNTENIC,
        NetClassification.TOP,
    ]
