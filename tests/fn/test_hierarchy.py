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
    with pytest.raises(ValueError, match=r"lines \[2\]"):
        bt.fn.load_hierarchy(path, "pathway")


def test_unknown_layout_raises():
    with pytest.raises(ValueError, match="layout="):
        bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group", layout="wide")
