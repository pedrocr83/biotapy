import numpy as np
import pandas as pd
import pytest

from biotapy._core import get_tree, make_treedata, tree_from_edges


def _make(tree):
    return make_treedata(
        np.ones((1, 2)),
        obs=pd.DataFrame(index=["s1"]),
        var=pd.DataFrame(index=["a", "b"]),
        tree=tree,
        x_kind="counts",
        source="test",
    )


def test_tree_from_edges_stores_branch_length():
    assert tree_from_edges([("r", "a", 0.5)]).edges["r", "a"]["length"] == 0.5


def test_make_treedata_adds_no_label_column():
    assert "tree" not in _make(tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)])).var.columns


def test_get_tree_without_phylogeny_names_the_key():
    with pytest.raises(KeyError, match="phylo"):
        get_tree(_make(None))
