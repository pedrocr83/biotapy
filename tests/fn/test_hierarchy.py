import gzip
from pathlib import Path

import pytest

import biotapy as bt

DATA = Path(__file__).parents[1] / "data" / "humann"


def test_parent_first_lines_become_child_parent_rows():
    edges = bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group")
    assert edges[["child", "parent"]].values.tolist() == [
        ["UniRef90_A", "G1"],
        ["UniRef90_B", "G1"],
        ["UniRef90_B", "G2"],
        ["UniRef90_C", "G2"],
        ["UniRef90_E", "G2"],
        ["UniRef90_E", "G3"],
    ]
    assert list(edges.columns) == ["child", "parent", "level", "parent_name"]
    assert set(edges["level"]) == {"group"} and edges["parent_name"].isna().all()


def test_child_first_reads_two_column_tables(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("ko:K00001\tpath:map00010\nko:K00001\tpath:map00071\nko:K00002\tpath:map00010\n")
    edges = bt.fn.load_hierarchy(path, "pathway", layout="child_first")
    assert edges[["child", "parent"]].values.tolist() == [
        ["ko:K00001", "path:map00010"],
        ["ko:K00001", "path:map00071"],
        ["ko:K00002", "path:map00010"],
    ]


def test_repeated_lines_and_blank_lines_give_distinct_pairs(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text("P1\tK1\n\nP1\tK1\tK2\t\n")
    assert len(bt.fn.load_hierarchy(path, "pathway")) == 2


def test_reads_gzip(tmp_path):
    path = tmp_path / "map.tsv.gz"
    path.write_bytes(gzip.compress(b"P1\tK1\tK2\n"))
    assert len(bt.fn.load_hierarchy(path, "pathway")) == 2


def test_attrs_name_the_file():
    assert bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group").attrs["source"].endswith("regroup_map.tsv")


def test_a_line_with_one_id_is_named(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text("P1\tK1\nP2\n")
    with pytest.raises(ValueError, match=r"lines 2 "):
        bt.fn.load_hierarchy(path, "pathway")


def test_unknown_layout_raises():
    with pytest.raises(ValueError, match="layout="):
        bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group", layout="wide")


@pytest.mark.parametrize("text", ["", "\n\n", "# only a comment\n"], ids=["empty", "blank", "comment"])
def test_a_file_with_no_edges_raises(tmp_path, text):
    path = tmp_path / "map.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match="path=.*no edges"):
        bt.fn.load_hierarchy(path, "pathway")


@pytest.mark.parametrize("line", ["\tK1\tK2", "  \tK1", "P1\t\tK2"], ids=["leading tab", "blank first", "middle gap"])
def test_an_empty_cell_before_the_last_id_is_named(tmp_path, line):
    path = tmp_path / "map.tsv"
    path.write_text(f"P0\tK0\n{line}\n")
    with pytest.raises(ValueError, match=r"path=.*lines 2 have an empty cell"):
        bt.fn.load_hierarchy(path, "pathway")


def test_comment_lines_are_skipped(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text("# pathway\tchild\nP1\tK1\n")
    assert bt.fn.load_hierarchy(path, "pathway")[["child", "parent"]].values.tolist() == [["K1", "P1"]]


def test_a_utf8_bom_is_not_part_of_the_first_id(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_bytes(b"\xef\xbb\xbfP1\tK1\n")
    assert bt.fn.load_hierarchy(path, "pathway")["parent"].tolist() == ["P1"]


def test_whitespace_only_lines_are_blank(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text("P1\tK1\n   \n \t \nP2\tK2\n")
    assert len(bt.fn.load_hierarchy(path, "pathway")) == 2


def test_long_line_lists_say_they_are_truncated(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text("P1\tK1\n" + "P\n" * 7)
    with pytest.raises(ValueError, match=r"lines 2, 3, 4 \(and 4 more\)"):
        bt.fn.load_hierarchy(path, "pathway")
