import tempfile
import warnings
from pathlib import Path

import h5py
import mudata
import networkx as nx
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from hypothesis import given, settings
from hypothesis import strategies as st
from mudata import MuData

import biotapy as bt
from biotapy._core import TreeData, get_tree


def _tree_with_attributes(leaves: list[str]) -> nx.DiGraph:
    tree = nx.DiGraph()
    tree.add_edge("root", "left", length=0.5)
    tree.add_edge("root", "right", length=0.25)
    for index, leaf in enumerate(leaves):
        tree.add_edge("left" if index % 2 else "right", leaf, length=1.0 + index)
    nx.set_node_attributes(tree, {"left": 0.9}, "support")
    return tree


def _rich_tdata(*, alignment: str = "leaves", label: str | None = "tree") -> TreeData:
    base = bt.datasets.toy()
    obst = {"samples": _tree_with_attributes(base.obs_names.tolist())}
    tdata = TreeData(
        base.X.copy(),
        obs=base.obs.copy(),
        var=base.var.copy(),
        vart={"phylo": get_tree(base)},
        obst=obst,
        label=label,
        allow_overlap=False,
        alignment=alignment,
    )
    tdata.layers["rel"] = tdata.X.copy() / 2
    tdata.obsm["pc"] = np.arange(tdata.n_obs * 2, dtype=float).reshape(tdata.n_obs, 2)
    tdata.obsp["dist"] = np.ones((tdata.n_obs, tdata.n_obs))
    tdata.uns["note"] = "kept"
    return tdata


def _metabolites(like: AnnData) -> AnnData:
    values = np.arange(like.n_obs * 3, dtype=np.float64).reshape(like.n_obs, 3)
    return AnnData(X=values, obs=like.obs[[]].copy(), var={"name": ["m1", "m2", "m3"]})


def _assert_same_trees(left: nx.DiGraph, right: nx.DiGraph) -> None:
    assert set(left.nodes) == set(right.nodes)
    assert {(a, b, d["length"]) for a, b, d in left.edges(data=True)} == {
        (a, b, d["length"]) for a, b, d in right.edges(data=True)
    }
    assert dict(left.nodes(data=True)) == dict(right.nodes(data=True))


def test_round_trip_keeps_a_treedata_modality_and_a_plain_one(tmp_path):
    tdata = _rich_tdata()
    mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata)})
    path = tmp_path / "study.h5mu"
    bt.io.write_h5mu(mdata, path)
    back = bt.io.read_h5mu(path)
    assert isinstance(back, MuData) and list(back.mod) == ["taxa", "metabolites"]
    taxa = back["taxa"]
    assert type(taxa) is TreeData
    assert type(back["metabolites"]) is AnnData
    _assert_same_trees(get_tree(taxa), get_tree(tdata))
    _assert_same_trees(taxa.obst["samples"], tdata.obst["samples"])
    assert (taxa.label, taxa.allow_overlap, taxa.alignment) == ("tree", False, "leaves")
    assert taxa.obs_names.tolist() == tdata.obs_names.tolist()
    assert taxa.var_names.tolist() == tdata.var_names.tolist()
    np.testing.assert_array_equal(taxa.X.toarray(), tdata.X.toarray())
    np.testing.assert_array_equal(taxa.layers["rel"].toarray(), tdata.layers["rel"].toarray())
    np.testing.assert_array_equal(taxa.obsm["pc"], tdata.obsm["pc"])
    np.testing.assert_array_equal(taxa.obsp["dist"], tdata.obsp["dist"])
    assert taxa.uns["note"] == "kept"


def test_the_toy_dataset_keeps_its_tree(tmp_path):
    toy = bt.datasets.toy()
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": toy}), tmp_path / "toy.h5mu")
    taxa = bt.io.read_h5mu(tmp_path / "toy.h5mu")["taxa"]
    _assert_same_trees(get_tree(taxa), get_tree(toy))
    # mudata writes string columns as categoricals; the values are what must survive.
    pd.testing.assert_frame_equal(taxa.var, toy.var, check_dtype=False, check_categorical=False)


def test_plain_mudata_still_reads_the_file_without_a_warning(tmp_path):
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": _rich_tdata()}), tmp_path / "study.h5mu")
    # Only warnings about the file matter here (UserWarning and its subclasses, such as anndata's); dependencies'
    # deprecation warnings, such as pre-release pandas' on mudata's own `drop(inplace=True)`, are not this test's concern.
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        plain = mudata.read_h5mu(tmp_path / "study.h5mu")
    assert type(plain["taxa"]) is AnnData


def test_a_plain_mudata_file_reads_back_as_anndata(tmp_path):
    mdata = bt.io.to_mudata({"metabolites": _metabolites(bt.datasets.toy())})
    mdata.write_h5mu(tmp_path / "plain.h5mu")
    assert type(bt.io.read_h5mu(tmp_path / "plain.h5mu")["metabolites"]) is AnnData


def test_a_treedata_without_trees_round_trips(tmp_path):
    base = bt.datasets.toy()
    bare = TreeData(base.X.copy(), obs=base.obs.copy(), var=base.var.copy())
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": bare}), tmp_path / "bare.h5mu")
    taxa = bt.io.read_h5mu(tmp_path / "bare.h5mu")["taxa"]
    assert type(taxa) is TreeData and len(taxa.vart) == 0 and len(taxa.obst) == 0


def test_label_none_round_trips(tmp_path):
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": _rich_tdata(label=None)}), tmp_path / "l.h5mu")
    assert bt.io.read_h5mu(tmp_path / "l.h5mu")["taxa"].label is None


def test_subset_alignment_round_trips(tmp_path):
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": _rich_tdata(alignment="subset")}), tmp_path / "s.h5mu")
    assert bt.io.read_h5mu(tmp_path / "s.h5mu")["taxa"].alignment == "subset"


def test_writing_does_not_change_the_input(tmp_path, assert_unchanged):
    tdata = _rich_tdata()
    mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata)})
    before = mdata["taxa"].copy()
    before_obs_dtypes = mdata["taxa"].obs.dtypes.copy()
    bt.io.write_h5mu(mdata, tmp_path / "study.h5mu")
    assert_unchanged(before, mdata["taxa"])
    assert mdata["taxa"].obs.dtypes.equals(before_obs_dtypes)
    assert type(mdata["taxa"]) is TreeData
    _assert_same_trees(get_tree(mdata["taxa"]), get_tree(before))


def _caterpillar(leaves: list[str], lengths: list[float]) -> nx.DiGraph:
    tree = nx.DiGraph()
    parent = "root"
    for index, leaf in enumerate(leaves[:-1]):
        tree.add_edge(parent, leaf, length=lengths[index])
        tree.add_edge(parent, f"n{index}", length=lengths[index])
        parent = f"n{index}"
    tree.add_edge(parent, leaves[-1], length=lengths[-1])
    return tree


@st.composite
def _random_trees(draw, leaves: list[str]) -> nx.DiGraph:
    """A binary tree over ``leaves`` from a drawn merge order, with drawn positive branch lengths."""
    tree = nx.DiGraph()
    pending = list(leaves)
    for step in range(len(leaves) - 1):
        first = draw(st.integers(0, len(pending) - 1))
        second = draw(st.integers(0, len(pending) - 2))
        second += second >= first
        parent = f"n{step}"
        for child in (pending[first], pending[second]):
            tree.add_edge(parent, child, length=draw(st.floats(0.01, 10.0)))
        pending = [node for index, node in enumerate(pending) if index not in (first, second)] + [parent]
    return tree


def _small_tdata(
    counts: np.ndarray,
    *,
    var: pd.DataFrame | None = None,
    vart: bool = True,
    tree: nx.DiGraph | None = None,
) -> TreeData:
    features = [f"f_{i}" for i in range(counts.shape[1])]
    samples = [f"s_{i}" for i in range(counts.shape[0])]
    tree = tree if tree is not None else _caterpillar(features, [0.5 + i for i in range(len(features))])
    return TreeData(
        sp.csr_matrix(counts.astype(np.float64)),
        obs=pd.DataFrame(index=samples),
        var=var if var is not None else pd.DataFrame(index=features),
        vart={"phylo": tree} if vart else None,
    )


def _round_trip(tdata: TreeData, path) -> TreeData:
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": tdata}), path)
    return bt.io.read_h5mu(path)["taxa"]


@st.composite
def _tables_with_trees(draw) -> tuple[np.ndarray, nx.DiGraph]:
    n_samples, n_features = draw(st.integers(1, 4)), draw(st.integers(1, 5))
    counts = draw(
        st.lists(
            st.lists(st.integers(0, 20), min_size=n_features, max_size=n_features),
            min_size=n_samples,
            max_size=n_samples,
        )
    )
    return np.array(counts), draw(_random_trees([f"f_{i}" for i in range(n_features)]))


@settings(max_examples=25, deadline=None)
@given(_tables_with_trees())
def test_round_trip_keeps_the_table_the_names_and_the_tree(table):
    counts, tree = table
    tdata = _small_tdata(counts, tree=tree)
    with tempfile.TemporaryDirectory() as directory:
        back = _round_trip(tdata, Path(directory) / "p.h5mu")
    np.testing.assert_array_equal(back.X.toarray(), tdata.X.toarray())
    assert back.obs_names.tolist() == tdata.obs_names.tolist()
    assert back.var_names.tolist() == tdata.var_names.tolist()
    _assert_same_trees(get_tree(back), get_tree(tdata))


@pytest.mark.parametrize(
    "counts",
    [
        np.array([[0, 0, 0], [1, 2, 3]]),
        np.array([[0, 1, 2], [0, 3, 4]]),
        np.array([[4, 5, 6]]),
    ],
    ids=["all-zero sample", "all-zero feature", "single sample"],
)
def test_degenerate_tables_round_trip(counts, tmp_path):
    tdata = _small_tdata(counts)
    back = _round_trip(tdata, tmp_path / "d.h5mu")
    np.testing.assert_array_equal(back.X.toarray(), tdata.X.toarray())
    _assert_same_trees(get_tree(back), get_tree(tdata))


def test_a_missing_taxonomy_rank_round_trips(tmp_path):
    var = pd.DataFrame({"genus": ["Alpha", None, "Beta"]}, index=["f_0", "f_1", "f_2"])
    back = _round_trip(_small_tdata(np.ones((2, 3)), var=var), tmp_path / "n.h5mu")
    assert back.var["genus"].isna().tolist() == [False, True, False]


def test_trees_on_the_sample_axis_only_round_trip(tmp_path):
    base = _small_tdata(np.ones((3, 2)), vart=False)
    base.obst["samples"] = _tree_with_attributes(base.obs_names.tolist())
    back = _round_trip(base, tmp_path / "o.h5mu")
    assert len(back.vart) == 0
    _assert_same_trees(back.obst["samples"], base.obst["samples"])


def test_trees_on_the_feature_axis_only_round_trip(tmp_path):
    base = _small_tdata(np.ones((3, 2)))
    back = _round_trip(base, tmp_path / "v.h5mu")
    assert len(back.obst) == 0
    _assert_same_trees(get_tree(back), get_tree(base))


def test_a_newer_layout_version_raises(tmp_path):
    path = tmp_path / "study.h5mu"
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": bt.datasets.toy()}), path)
    with h5py.File(path, "a") as handle:
        handle["mod"]["taxa"].attrs["biotapy-treedata-encoding"] = "99"
    with pytest.raises(
        ValueError, match=r"unknown tree-slot layout version '99' \(this biotapy reads '1'\).*newer biotapy"
    ):
        bt.io.read_h5mu(path)


def test_a_non_mudata_raises_naming_mdata(tmp_path):
    with pytest.raises(TypeError, match="mdata must be a MuData, got TreeData"):
        bt.io.write_h5mu(bt.datasets.toy(), tmp_path / "x.h5mu")


def test_a_backed_mudata_raises_naming_mdata(tmp_path):
    bt.io.write_h5mu(bt.io.to_mudata({"taxa": bt.datasets.toy()}), tmp_path / "in.h5mu")
    backed = mudata.read_h5mu(tmp_path / "in.h5mu", backed=True)
    try:
        with pytest.raises(ValueError, match=r"mdata is backed.*backed=True"):
            bt.io.write_h5mu(backed, tmp_path / "out.h5mu")
    finally:
        backed.file.close()
