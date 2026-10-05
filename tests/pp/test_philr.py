import json

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt
from biotapy._core import TreeData, get_tree, tree_from_edges


def _toy6() -> TreeData:
    # toy()'s root has three children; its first six features sit under a binary root.
    return bt.datasets.toy()[:, :6].copy()


def _gmean(values: np.ndarray) -> float:
    return float(np.exp(np.log(values).mean()))


def test_balances_are_scaled_log_ratios_of_geometric_means():
    tdata = _toy6()
    x = tdata.X.toarray()[0] + 0.5
    balances = bt.pp.philr(tdata).obsm["X_philr"].loc["s1"]
    # n4 = (f1, f2); n1 = (n4, f3); root = (n1, n2): the first child is the numerator, as in philr::philr.
    assert balances["n4"] == pytest.approx(np.sqrt(1 / 2) * np.log(x[0] / x[1]))
    assert balances["n1"] == pytest.approx(np.sqrt(2 * 1 / 3) * np.log(_gmean(x[:2]) / x[2]))
    assert balances["root"] == pytest.approx(np.sqrt(3 * 3 / 6) * np.log(_gmean(x[:3]) / _gmean(x[3:])))


def test_one_column_per_internal_node_in_preorder():
    out = bt.pp.philr(_toy6()).obsm["X_philr"]
    assert isinstance(out, pd.DataFrame) and out.index.tolist() == [f"s{i}" for i in range(1, 7)]
    assert out.columns.tolist() == ["root", "n1", "n4", "n2", "n5"]


def test_keeps_x_tree_and_input(assert_unchanged):
    tdata = _toy6()
    before = tdata.copy()
    edges = [(u, v, dict(data)) for u, v, data in get_tree(tdata).edges(data=True)]
    out = bt.pp.philr(tdata)
    assert_unchanged(before, tdata)
    assert list(get_tree(tdata).edges(data=True)) == edges
    assert isinstance(out, TreeData) and "phylo" in out.vart and (out.X != tdata.X).nnz == 0


def test_one_child_nodes_are_skipped_keeping_child_order():
    # Without f2, n4 has one child: n1 becomes (f1, f3) with f1 still first, so still the numerator.
    tdata = _toy6()[:, ["f1", "f3", "f4", "f5", "f6"]].copy()
    x = tdata.X.toarray()[0] + 0.5
    out = bt.pp.philr(tdata).obsm["X_philr"]
    assert out.columns.tolist() == ["root", "n1", "n2", "n5"]
    assert out.loc["s1", "n1"] == pytest.approx(np.sqrt(1 / 2) * np.log(x[0] / x[1]))


def test_three_child_root_raises():
    with pytest.raises(ValueError, match=r"the node\(s\) \['root'\] have more than two children"):
        bt.pp.philr(bt.datasets.toy())


def test_internal_polytomy_raises():
    edges = [("r", "a", 1.0), ("r", "p", 1.0), ("p", "b", 1.0), ("p", "c", 1.0), ("p", "d", 1.0)]
    tdata = TreeData(
        X=sp.csr_matrix(np.ones((2, 4))),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=list("abcd")),
        vart={"phylo": tree_from_edges(edges)},
        label=None,
    )
    with pytest.raises(ValueError, match=r"\['p'\] have more than two children"):
        bt.pp.philr(tdata)


def test_feature_outside_the_tree_raises():
    tdata = TreeData(
        X=sp.csr_matrix(np.ones((2, 3))),
        obs=pd.DataFrame(index=["s1", "s2"]),
        var=pd.DataFrame(index=["a", "b", "extra"]),
        vart={"phylo": tree_from_edges([("r", "a", 1.0), ("r", "b", 1.0)])},
        label=None,
    )
    with pytest.raises(ValueError, match=r"every feature to be a tip of the tree; 1 feature\(s\) are not"):
        bt.pp.philr(tdata)


def test_fewer_than_two_features_raise():
    with pytest.raises(ValueError, match="at least two features, got 1"):
        bt.pp.philr(_toy6()[:, :1].copy())


def test_plain_anndata_raises():
    with pytest.raises(TypeError, match="needs a TreeData"):
        bt.pp.philr(_toy6().to_adata())


def test_all_zero_sample_has_zero_balances():
    tdata = _toy6()
    dense = tdata.X.toarray()
    dense[0] = 0
    tdata.X = sp.csr_matrix(dense)
    np.testing.assert_array_equal(bt.pp.philr(tdata).obsm["X_philr"].loc["s1"], 0.0)


def test_all_zero_feature_is_finite():
    tdata = _toy6()
    dense = tdata.X.toarray()
    dense[:, 0] = 0
    tdata.X = sp.csr_matrix(dense)
    assert np.isfinite(bt.pp.philr(tdata).obsm["X_philr"].to_numpy()).all()


def test_single_sample():
    assert bt.pp.philr(_toy6()[:1].copy()).obsm["X_philr"].shape == (1, 5)


def test_philr_pseudocount_above_the_smallest_value_warns():
    relative = _toy6()
    relative.X = bt.pp.relative(relative).layers["relative"]
    with pytest.warns(UserWarning, match=r"pseudocount=0.5 is larger than the smallest non-zero value in X \(0.0"):
        bt.pp.philr(relative)


def test_zero_without_pseudocount_raises():
    with pytest.raises(ValueError, match="pass pseudocount > 0 to pp.philr"):
        bt.pp.philr(_toy6(), pseudocount=0)


@pytest.mark.parametrize("pseudocount", [True, "0.5", None])
def test_non_numeric_pseudocount_raises(pseudocount):
    with pytest.raises(TypeError, match="pseudocount must be a real number"):
        bt.pp.philr(_toy6(), pseudocount=pseudocount)


def test_records_provenance():
    entries = bt.pp.philr(_toy6()).uns["biotapy"]["provenance"]
    assert json.loads(entries[-1])["step"] == "pp.philr" and json.loads(entries[-1])["params"] == {"pseudocount": 0.5}


def test_round_trips_through_h5td(tmp_path):
    import treedata

    out = bt.pp.philr(_toy6())
    out.write_h5td(tmp_path / "philr.h5td")
    back = treedata.read_h5td(tmp_path / "philr.h5td")
    pd.testing.assert_frame_equal(back.obsm["X_philr"], out.obsm["X_philr"])


def test_feature_filter_drops_the_balances():
    out = bt.pp.filter_features(bt.pp.philr(_toy6()), min_prevalence=0.5)
    assert "X_philr" not in out.obsm


@settings(deadline=None)
@given(arrays(np.int64, st.tuples(st.integers(1, 5), st.just(6)), elements=st.integers(0, 1000)), st.floats(0.01, 1))
def test_balances_keep_the_clr_distance(dense, pseudocount):
    # An isometric log-ratio basis: each sample's balances have the length of its CLR vector.
    tdata = _toy6()[: dense.shape[0]].copy()
    tdata.X = sp.csr_matrix(dense)
    balances = bt.pp.philr(tdata, pseudocount=pseudocount).obsm["X_philr"].to_numpy()
    clr = bt.pp.clr(tdata, pseudocount=pseudocount).layers["clr"]
    np.testing.assert_allclose(np.linalg.norm(balances, axis=1), np.linalg.norm(clr, axis=1), rtol=1e-9, atol=1e-9)
