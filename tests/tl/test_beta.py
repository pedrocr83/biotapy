import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt
from biotapy._core import get_tree, make_treedata, tree_from_edges


def _adata(dense) -> ad.AnnData:
    return ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
    )


def test_braycurtis_by_hand():
    out = bt.tl.beta(_adata(np.array([[1, 2, 3], [3, 2, 1]])))
    assert out.loc["s0", "s1"] == pytest.approx(4 / 12)
    assert list(out.index) == list(out.columns) == ["s0", "s1"]


def test_jaccard_is_presence_absence():
    out = bt.tl.beta(_adata(np.array([[1, 0, 30], [5, 5, 0]])), metric="jaccard")
    assert out.loc["s0", "s1"] == pytest.approx(2 / 3)


def test_all_zero_sample_and_feature():
    out = bt.tl.beta(_adata(np.array([[0, 0, 0], [1, 2, 0], [0, 0, 0]])))
    assert out.loc["s0", "s1"] == 1.0 and np.isnan(out.loc["s0", "s2"])


def test_single_sample():
    assert bt.tl.beta(_adata(np.array([[1, 2]]))).to_numpy().tolist() == [[0.0]]


def test_missing_rank_is_irrelevant():
    tdata = bt.datasets.toy()
    tdata.var["genus"] = np.nan
    assert bt.tl.beta(tdata).shape == (6, 6)


def test_beta_inplace_writes_obsp(assert_unchanged):
    tdata = bt.datasets.toy()
    expected = bt.tl.beta(tdata, metric="jaccard")
    assert bt.tl.beta(tdata, metric="jaccard", inplace=True) is None
    np.testing.assert_array_equal(tdata.obsp["jaccard"], expected.to_numpy())


def test_beta_input_unchanged(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.tl.beta(tdata)
    assert_unchanged(before, tdata)


def test_unknown_metric_raises():
    with pytest.raises(ValueError, match="metric must be one of"):
        bt.tl.beta(bt.datasets.toy(), metric="euclidean")


@given(
    arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 20)),
    st.sampled_from(["braycurtis", "jaccard"]),
)
def test_beta_is_symmetric_with_zero_diagonal(dense, metric):
    out = bt.tl.beta(_adata(dense), metric=metric).to_numpy()
    np.testing.assert_array_equal(out, out.T)
    np.testing.assert_array_equal(np.diag(out), 0.0)


def test_unifrac_by_hand():
    # root -> a (1), root -> b (3); s0 holds a, s1 holds b: unweighted 4/4, weighted normalized 1.
    tdata = make_treedata(
        np.array([[2, 0], [0, 5]]),
        obs=pd.DataFrame(index=["s0", "s1"]),
        var=pd.DataFrame(index=["a", "b"]),
        tree=tree_from_edges([("r", "a", 1.0), ("r", "b", 3.0)]),
        x_kind="counts",
        source="test",
    )
    assert bt.tl.unifrac(tdata).loc["s0", "s1"] == pytest.approx(1.0)
    assert bt.tl.unifrac(tdata, weighted=True).loc["s0", "s1"] == pytest.approx(1.0)
    assert bt.tl.unifrac(tdata, weighted=True, normalized=False).loc["s0", "s1"] == pytest.approx(4.0)


def test_unifrac_multifurcating_root_is_accepted():
    # toy's root has three children; scikit-bio alone raises "The tree must be rooted."
    out = bt.tl.unifrac(bt.datasets.toy())
    assert out.shape == (6, 6) and out.loc["s1", "s4"] == pytest.approx(0.0916030534)


def test_unifrac_after_filtering_keeps_path_lengths():
    # Filtering leaves unary nodes in the tree; distances between the kept features' samples must not change.
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[:, [4, 7]] = 0  # f5, f8 absent everywhere
    tdata.X = sp.csr_matrix(dense)
    filtered = bt.pp.filter_features(tdata, min_total=1)
    assert filtered.n_vars == 6
    pd.testing.assert_frame_equal(bt.tl.unifrac(filtered), bt.tl.unifrac(tdata))
    pd.testing.assert_frame_equal(bt.tl.unifrac(filtered, weighted=True), bt.tl.unifrac(tdata, weighted=True))


def _toy_with_f1_length(length: float) -> ad.AnnData:
    tdata = bt.datasets.toy()
    get_tree(tdata).edges["n4", "f1"]["length"] = length
    return tdata


@pytest.mark.parametrize("weighted", [False, True])
def test_unifrac_counts_a_nan_branch_length_as_zero(weighted):
    out = bt.tl.unifrac(_toy_with_f1_length(np.nan), weighted=weighted)
    assert np.isfinite(out.to_numpy()).all()
    pd.testing.assert_frame_equal(out, bt.tl.unifrac(_toy_with_f1_length(0.0), weighted=weighted))


def test_unifrac_all_zero_sample_single_sample_and_missing_rank():
    tdata = bt.datasets.toy()
    tdata.var["genus"] = np.nan
    assert bt.tl.unifrac(tdata[:1].copy()).shape == (1, 1)
    dense = tdata.X.toarray()
    dense[0] = 0
    tdata.X = sp.csr_matrix(dense)
    assert bt.tl.unifrac(tdata, weighted=True).shape == (6, 6)


def test_two_all_zero_samples_pin_scikit_bios_convention():
    # scikit-bio's own values, kept as is with no custom NaN mapping (rules.md R2.1).
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[[0, 1]] = 0  # s1, s2 both empty
    tdata.X = sp.csr_matrix(dense)
    assert np.isnan(bt.tl.beta(tdata).loc["s1", "s2"])
    assert bt.tl.beta(tdata, metric="jaccard").loc["s1", "s2"] == 0.0
    assert bt.tl.unifrac(tdata).loc["s1", "s2"] == 0.0
    assert bt.tl.unifrac(tdata, weighted=True).loc["s1", "s2"] == 0.0


def _toy_proportions(x_kind: str) -> ad.AnnData:
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    tdata.X = sp.csr_matrix(dense / dense.sum(axis=1, keepdims=True))
    tdata.uns["biotapy"]["x_kind"] = x_kind
    return tdata


@pytest.mark.parametrize("x_kind", ["relative", "counts"])
def test_weighted_unifrac_needs_counts(x_kind):
    # scikit-bio's tree code casts abundances to int64, which would put every pair of proportions 0 apart.
    with pytest.raises(ValueError, match=r"tl.unifrac\(weighted=True\) needs raw counts"):
        bt.tl.unifrac(_toy_proportions(x_kind), weighted=True)


def test_unweighted_unifrac_on_proportions_equals_counts():
    pd.testing.assert_frame_equal(bt.tl.unifrac(_toy_proportions("relative")), bt.tl.unifrac(bt.datasets.toy()))


def test_unifrac_inplace_writes_the_contract_keys():
    tdata = bt.datasets.toy()
    assert bt.tl.unifrac(tdata, inplace=True) is None
    assert bt.tl.unifrac(tdata, weighted=True, inplace=True) is None
    assert {"unweighted_unifrac", "weighted_unifrac"} <= set(tdata.obsp.keys())


def test_unifrac_input_unchanged(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.tl.unifrac(tdata, weighted=True)
    assert_unchanged(before, tdata)


def test_unifrac_without_a_tree_names_it():
    with pytest.raises(TypeError, match="needs a TreeData"):
        bt.tl.unifrac(_adata(np.array([[1, 2]])))
    tdata = bt.datasets.toy()
    del tdata.vart["phylo"]
    with pytest.raises(KeyError, match="phylo"):
        bt.tl.unifrac(tdata)


@given(arrays(np.int64, st.tuples(st.integers(1, 5), st.just(8)), elements=st.integers(0, 20)), st.booleans())
def test_unifrac_is_symmetric_with_zero_diagonal(dense, weighted):
    tdata = bt.datasets.toy()[: dense.shape[0]].copy()
    tdata.X = sp.csr_matrix(dense)
    out = bt.tl.unifrac(tdata, weighted=weighted).to_numpy()
    np.testing.assert_array_equal(out, out.T)
    np.testing.assert_array_equal(np.diag(out), 0.0)
