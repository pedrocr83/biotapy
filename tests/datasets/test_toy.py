import numpy as np
import treedata as td

import biotapy as bt


def test_toy_shape_taxonomy_and_tree():
    tdata = bt.datasets.toy()
    assert tdata.shape == (6, 8)
    assert tdata.var["genus"].isna().tolist() == [False] * 7 + [True]
    assert set(tdata.vart["phylo"].successors("n5")) == {"f4", "f5"}


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


def test_toy_humann_has_both_modalities_over_the_toy_samples():
    mdata = bt.datasets.toy_humann()
    assert mdata["function"].shape == (6, 6) and mdata["function_by_taxon"].shape == (6, 7)
    assert mdata["function"].obs_names.tolist() == bt.datasets.toy().obs_names.tolist()
    assert mdata["function"].obs["group"].tolist() == ["A"] * 3 + ["B"] * 3
    assert mdata["function"].uns["biotapy"]["x_kind"] == "rpk"


def test_toy_humann_community_rows_are_the_sums_of_their_strata():
    mdata = bt.datasets.toy_humann()
    function, by_taxon = mdata["function"], mdata["function_by_taxon"]
    for name in ["UNGROUPED", "1.1.1.1", "2.7.1.1", "2.7.1.2", "3.2.1.4"]:
        strata = by_taxon[:, (by_taxon.var["function"] == name).to_numpy()].X.sum(axis=1)
        np.testing.assert_array_equal(np.asarray(strata).ravel(), function[:, name].X.toarray().ravel())
