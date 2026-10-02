import tracemalloc

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt
from biotapy._core import get_tree

ALL = ["observed_features", "shannon", "simpson", "chao1", "faith_pd"]


def test_default_metrics_by_hand(make_adata):
    out = bt.tl.alpha(make_adata(np.array([[1, 2, 3, 0], [4, 4, 0, 0]])))
    assert list(out.columns) == ["observed_features", "shannon", "simpson", "chao1"]
    p = np.array([1, 2, 3]) / 6
    np.testing.assert_allclose(out.loc["s0"].to_numpy(), [3, -(p * np.log(p)).sum(), 1 - (p**2).sum(), 3], rtol=1e-12)
    np.testing.assert_allclose(out.loc["s1"].to_numpy(), [2, np.log(2), 0.5, 2], rtol=1e-12)


def test_chao1_is_bias_corrected(make_adata):
    # 2 singletons, 1 doubleton: 4 + 2 * 1 / (2 * (1 + 1)) = 4.5
    assert bt.tl.alpha(make_adata(np.array([[1, 1, 2, 5]])), metrics=["chao1"]).loc["s0", "chao1"] == 4.5


def test_faith_pd_sums_branches_to_the_root():
    # s1 lacks f5 and f8: every branch except n5-f5 (0.03) and n3-f8 (0.12); toy's branches total 1.34.
    out = bt.tl.alpha(bt.datasets.toy(), metrics=["faith_pd"])
    assert out.loc["s1", "faith_pd"] == pytest.approx(1.34 - 0.03 - 0.12)


def _toy_with_f1_length(length: float) -> ad.AnnData:
    tdata = bt.datasets.toy()
    get_tree(tdata).edges["n4", "f1"]["length"] = length
    return tdata


def test_faith_pd_counts_a_nan_branch_length_as_zero():
    out = bt.tl.alpha(_toy_with_f1_length(np.nan), metrics=["faith_pd"])
    assert np.isfinite(out["faith_pd"]).all()
    pd.testing.assert_frame_equal(out, bt.tl.alpha(_toy_with_f1_length(0.0), metrics=["faith_pd"]))


def test_all_zero_sample_is_nan_or_zero_never_raises():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[0] = 0
    tdata.X = sp.csr_matrix(dense)
    out = bt.tl.alpha(tdata, metrics=ALL)
    assert out.loc["s1", ["observed_features", "chao1", "faith_pd"]].tolist() == [0, 0, 0]
    assert out.loc["s1", ["shannon", "simpson"]].isna().all()


def test_all_zero_feature_changes_nothing(make_adata):
    dense = np.array([[1, 2, 0], [3, 1, 0]])
    pd.testing.assert_frame_equal(bt.tl.alpha(make_adata(dense)), bt.tl.alpha(make_adata(dense[:, :2])))


def test_single_sample(make_adata):
    assert bt.tl.alpha(make_adata(np.array([[4, 0, 1]])), metrics=["observed_features"]).shape == (1, 1)


def test_missing_rank_is_irrelevant():
    tdata = bt.datasets.toy()
    tdata.var["genus"] = np.nan
    assert bt.tl.alpha(tdata).shape == (6, 4)


def test_inplace_writes_obs_and_returns_none():
    tdata = bt.datasets.toy()
    expected = bt.tl.alpha(tdata, metrics=["shannon", "faith_pd"])
    assert bt.tl.alpha(tdata, metrics=["shannon", "faith_pd"], inplace=True) is None
    np.testing.assert_array_equal(tdata.obs["alpha_shannon"], expected["shannon"])
    np.testing.assert_array_equal(tdata.obs["alpha_faith_pd"], expected["faith_pd"])


def test_input_unchanged(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.tl.alpha(tdata, metrics=ALL)
    assert_unchanged(before, tdata)


def test_densifies_in_bounded_chunks():
    # 2**19 + 1 features: each chunk of at most 2**20 values holds one sample (4 MiB);
    # densifying all 8 samples at once would take 32 MiB.
    n_obs, n_vars = 8, 2**19 + 1
    rng = np.random.default_rng(0)
    rows = np.repeat(np.arange(n_obs), 20)
    cols = rng.choice(n_vars, size=rows.size)
    X = sp.csr_matrix((rng.integers(1, 50, size=rows.size).astype(np.float64), (rows, cols)), shape=(n_obs, n_vars))
    wide = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"s{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(n_vars)]),
    )
    tracemalloc.start()
    out = bt.tl.alpha(wide)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak < n_obs * n_vars * 8 / 2
    pd.testing.assert_frame_equal(out, bt.tl.alpha(wide[:, np.unique(X.indices)].copy()))


def test_string_metrics_raise():
    with pytest.raises(TypeError, match="not the string 'shannon'"):
        bt.tl.alpha(bt.datasets.toy(), metrics="shannon")


def test_unknown_or_empty_metrics_raise():
    with pytest.raises(ValueError, match="metrics must name"):
        bt.tl.alpha(bt.datasets.toy(), metrics=["pielou"])
    with pytest.raises(ValueError, match="metrics must name"):
        bt.tl.alpha(bt.datasets.toy(), metrics=[])


def test_repeated_metrics_raise():
    with pytest.raises(ValueError, match=r"metrics repeats \['shannon'\]"):
        bt.tl.alpha(bt.datasets.toy(), metrics=["shannon", "simpson", "shannon"])


def test_count_metrics_need_counts_others_do_not():
    rel = bt.datasets.toy()
    rel.uns["biotapy"]["x_kind"] = "relative"
    with pytest.raises(ValueError, match=r"tl.alpha with \['chao1'\] needs raw counts"):
        bt.tl.alpha(rel, metrics=["shannon", "chao1"])
    assert bt.tl.alpha(rel, metrics=["shannon", "simpson", "faith_pd"]).shape == (6, 3)


def _toy_proportions() -> ad.AnnData:
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    tdata.X = sp.csr_matrix(dense / dense.sum(axis=1, keepdims=True))
    tdata.uns["biotapy"]["x_kind"] = "relative"
    return tdata


def test_faith_pd_on_proportions_equals_counts():
    # scikit-bio's tree code casts abundances to int64, which would turn every proportion into 0.
    out = bt.tl.alpha(_toy_proportions(), metrics=["faith_pd"])
    pd.testing.assert_frame_equal(out, bt.tl.alpha(bt.datasets.toy(), metrics=["faith_pd"]))


def test_count_metrics_reject_fractional_values_labelled_counts():
    tdata = _toy_proportions()
    tdata.uns["biotapy"]["x_kind"] = "counts"
    with pytest.raises(ValueError, match=r"tl.alpha with \['chao1'\] needs raw counts in X, but X holds non-integer"):
        bt.tl.alpha(tdata, metrics=["chao1"])


def test_faith_pd_needs_a_tree(make_adata):
    with pytest.raises(TypeError, match="needs a TreeData"):
        bt.tl.alpha(make_adata(np.array([[1, 2]])), metrics=["faith_pd"])
    tdata = bt.datasets.toy()
    del tdata.vart["phylo"]
    with pytest.raises(KeyError, match="phylo"):
        bt.tl.alpha(tdata, metrics=["faith_pd"])


@given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 8)), elements=st.integers(0, 20)))
def test_observed_features_counts_nonzero_features(make_adata, dense):
    out = bt.tl.alpha(make_adata(dense), metrics=["observed_features", "simpson"])
    np.testing.assert_array_equal(out["observed_features"], (dense > 0).sum(axis=1))
    assert (out["observed_features"] <= dense.shape[1]).all()
    simpson = out["simpson"].dropna()
    assert ((simpson >= 0) & (simpson < 1)).all()
