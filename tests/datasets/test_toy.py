import treedata as td

import biotapy as bt
from biotapy._core import get_tree


def test_toy_shape_taxonomy_and_tree():
    tdata = bt.datasets.toy()
    assert tdata.shape == (6, 8)
    assert tdata.var["genus"].isna().tolist() == [False] * 7 + [True]
    assert set(get_tree(tdata).successors("n5")) == {"f4", "f5"}


def test_toy_is_counts_with_one_provenance_entry():
    meta = bt.datasets.toy().uns["biotapy"]
    assert meta["x_kind"] == "counts" and len(meta["provenance"]) == 1


def test_toy_round_trips_through_h5td(tmp_path):
    path = tmp_path / "toy.h5td"
    bt.datasets.toy().write_h5td(path)
    back = td.read_h5td(path)
    assert (back.X != bt.datasets.toy().X).nnz == 0
    assert back.vart["phylo"].edges["n4", "f1"]["length"] == 0.1
    assert list(back.uns["biotapy"]["provenance"]) == list(bt.datasets.toy().uns["biotapy"]["provenance"])
