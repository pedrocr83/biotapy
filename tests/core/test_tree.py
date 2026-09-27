import math

import numpy as np
import pandas as pd
import pytest

from biotapy._core import get_tree, make_treedata, relabel_tips, tree_from_edges, tree_from_newick, tree_tips


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


def test_relabel_tips_renames_listed_tips_only():
    tree = relabel_tips(tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)]), {"a": "x"})
    assert _leaves(tree) == {"x", "b"} and tree.edges["r", "x"]["length"] == 1.0


def test_relabel_tips_rejects_a_name_already_in_the_tree():
    tree = tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)])
    with pytest.raises(ValueError, match="names=.*'b'"):
        relabel_tips(tree, {"a": "b"})


def test_relabel_tips_allows_swapping_names():
    tree = relabel_tips(tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)]), {"a": "b", "b": "a"})
    assert tree.edges["r", "b"]["length"] == 1.0 and tree.edges["r", "a"]["length"] == 2.0


def test_tree_tips_lists_nodes_without_children():
    assert sorted(tree_tips(tree_from_edges([("r", "a", 1.0), ("r", "n", 1.0), ("n", "b", 1.0)]))) == ["a", "b"]


def test_make_treedata_adds_no_label_column():
    assert "tree" not in _make(tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)])).var.columns


def test_get_tree_without_phylogeny_names_the_key():
    with pytest.raises(KeyError, match="phylo"):
        get_tree(_make(None))


def _leaves(tree) -> set[str]:
    return {n for n in tree.nodes if tree.out_degree(n) == 0}


def _frames(samples, features):
    return pd.DataFrame(index=samples), pd.DataFrame(index=features)


def test_tree_from_newick_keeps_underscores_in_tip_names():
    assert _leaves(tree_from_newick("((ASV_1:0.1,ASV_2:0.2):0.05,ASV_3:0.3);")) == {"ASV_1", "ASV_2", "ASV_3"}


def test_tree_from_newick_stores_branch_lengths():
    tree = tree_from_newick("((a:0.1,b:0.2):0.05,c:0.3);")
    parent = next(iter(tree.predecessors("a")))
    assert tree.edges[parent, "a"]["length"] == 0.1


def test_tree_from_newick_gives_internal_nodes_unique_names():
    tree = tree_from_newick("((a:1,b:1)0.95:1,(c:1,d:1)0.95:1);")
    internal = set(tree.nodes) - {"a", "b", "c", "d"}
    assert len(internal) == 3 and "0.95" not in internal


def test_tree_from_newick_internal_names_skip_tip_names():
    tree = tree_from_newick("((n0:1,n1:1):1,n2:1);")
    assert tree.number_of_nodes() == 5 and _leaves(tree) == {"n0", "n1", "n2"}


def test_tree_from_newick_missing_length_is_nan():
    tree = tree_from_newick("(a,b:2);")
    assert math.isnan(tree.edges[next(iter(tree.predecessors("a"))), "a"]["length"])


def test_tree_from_newick_repeated_tip_names_raise():
    with pytest.raises(ValueError, match="'a'"):
        tree_from_newick("(a:1,a:1);")


def test_tree_from_newick_single_tip_is_a_one_node_tree():
    tree = tree_from_newick("a;")
    assert set(tree.nodes) == {"a"}
    assert tree.number_of_edges() == 0


def test_make_treedata_casts_ids_to_str():
    obs, var = _frames([1, 2], [10, 20])
    tdata = make_treedata(np.ones((2, 2)), obs=obs, var=var, tree=None, x_kind="counts", source="test")
    assert list(tdata.obs_names) == ["1", "2"] and list(tdata.var_names) == ["10", "20"]


def test_make_treedata_names_duplicate_ids():
    obs, var = _frames(["s1", "s1"], ["a", "b"])
    with pytest.raises(ValueError, match="duplicate obs ids.*s1"):
        make_treedata(np.ones((2, 2)), obs=obs, var=var, tree=None, x_kind="counts", source="test")


@pytest.mark.parametrize("missing", [np.nan, None])
@pytest.mark.parametrize("axis", ["obs", "var"])
def test_make_treedata_rejects_missing_ids(axis, missing):
    obs, var = _frames(["s1", "s2"], ["a", "b"])
    frames = {"obs": obs, "var": var}
    frames[axis] = pd.DataFrame(index=pd.Index(["x", missing], dtype=object))
    with pytest.raises(ValueError, match=f"{axis} ids"):
        make_treedata(np.ones((2, 2)), **frames, tree=None, x_kind="counts", source="test")


def test_make_treedata_keeps_shared_features_and_warns_once():
    obs, var = _frames(["s1"], ["a", "b", "c"])
    tree = tree_from_edges([("r", "a", 1.0), ("r", "b", 1.0), ("r", "d", 1.0)])
    with pytest.warns(UserWarning, match=r"1 feature\(s\) not in the tree and 1 tree tip\(s\)") as record:
        tdata = make_treedata(np.array([[1, 2, 3]]), obs=obs, var=var, tree=tree, x_kind="counts", source="test")
    assert len(record) == 1
    assert list(tdata.var_names) == ["a", "b"]
    assert tdata.X.toarray().tolist() == [[1, 2]]
    assert _leaves(get_tree(tdata)) == {"a", "b"}


def test_make_treedata_without_shared_features_raises():
    obs, var = _frames(["s1"], ["a"])
    with pytest.raises(ValueError, match="no feature"):
        make_treedata(
            np.ones((1, 1)), obs=obs, var=var, tree=tree_from_edges([("r", "z", 1.0)]), x_kind="counts", source="test"
        )


def test_make_treedata_leaves_input_frames_alone():
    obs, var = _frames([1], ["a"])
    make_treedata(np.ones((1, 1)), obs=obs, var=var, tree=None, x_kind="counts", source="test")
    assert list(obs.index) == [1]
