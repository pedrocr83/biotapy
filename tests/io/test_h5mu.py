import warnings

import mudata
import networkx as nx
import numpy as np
import pandas as pd
from anndata import AnnData
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
    with warnings.catch_warnings():
        warnings.simplefilter("error")
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
