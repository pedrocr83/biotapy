import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.stats import t as t_dist

import biotapy as bt


def _counts_adata(dense, **obs):
    return AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(obs, index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
    )


def _ols(counts, design):
    """Group coefficient and its standard error for each feature, as R's lm on log2 CLR values."""
    values = counts + 0.5 if (counts == 0).any() else counts
    ratios = np.log2(values) - np.log2(values).mean(axis=1, keepdims=True)
    coef, *_ = np.linalg.lstsq(design, ratios, rcond=None)
    dof = design.shape[0] - design.shape[1]
    sigma2 = ((ratios - design @ coef) ** 2).sum(axis=0) / dof
    return coef[1], np.sqrt(sigma2 * np.linalg.inv(design.T @ design)[1, 1]), dof


def test_linda_is_least_squares_on_log2_ratios_shifted_by_one_bias():
    tdata = bt.datasets.toy()
    out = bt.da.linda(tdata, "group")
    group = (tdata.obs["group"] == "B").to_numpy(np.float64)
    beta, se, dof = _ols(tdata.X.toarray().astype(float), np.c_[np.ones(6), group])
    np.testing.assert_allclose(out["se"], se, rtol=1e-12)
    # Every feature's effect is its coefficient minus the same bias, the mode of all coefficients.
    np.testing.assert_allclose(np.ptp(beta - out["effect"]), 0, atol=1e-12)
    np.testing.assert_allclose(out["pvalue"], 2 * t_dist.sf(np.abs(out["effect"] / out["se"]), dof), rtol=1e-12)


def test_bias_correction_recovers_a_fold_change():
    # B repeats A with f5-f7 eight times as abundant: the CLR moves every other feature too, the bias removes that.
    rng = np.random.default_rng(0)
    a = rng.integers(20, 200, size=(4, 8)).astype(float)
    b = a * np.array([1, 1, 1, 1, 1, 8, 8, 8])
    out = bt.da.linda(_counts_adata(np.r_[a, b], g=["a"] * 4 + ["b"] * 4), "g")
    # Measured: the bias-corrected effects are 0.0093 from 0 and 3 (the mode estimate is not exact on 8 features).
    np.testing.assert_allclose(out["effect"].iloc[:5], 0, atol=0.01)
    np.testing.assert_allclose(out["effect"].iloc[5:], 3, atol=0.01)


def test_pseudocount_is_added_only_when_x_has_a_zero():
    counts = np.array([[5.0, 9, 2], [4, 8, 3], [7, 2, 6], [6, 3, 9]])
    out = bt.da.linda(_counts_adata(counts, g=["a", "a", "b", "b"]), "g")
    _, se, _ = _ols(counts, np.c_[np.ones(4), [0.0, 0, 1, 1]])
    np.testing.assert_allclose(out["se"], se, rtol=1e-12)


def test_covariates_enter_the_model():
    tdata = bt.datasets.toy()
    tdata.obs["batch"] = pd.Categorical(["x", "y", "y", "x", "y", "x"])
    tdata.obs["age"] = [30.0, 41, 52, 38, 45, 60]
    out = bt.da.linda(tdata, "group", covariates=["batch", "age"])
    group = (tdata.obs["group"] == "B").to_numpy(np.float64)
    batch = (tdata.obs["batch"] == "y").to_numpy(np.float64)
    age = tdata.obs["age"].to_numpy()
    design = np.c_[np.ones(6), group, batch, (age - age.mean()) / age.std(ddof=1)]
    _, se, _ = _ols(tdata.X.toarray().astype(float), design)
    np.testing.assert_allclose(out["se"], se, rtol=1e-12)


def test_numeric_group_is_per_standard_deviation():
    tdata = bt.datasets.toy()
    tdata.obs["ph"] = [5.1, 5.6, 6.0, 6.8, 7.1, 7.4]
    out = bt.da.linda(tdata, "ph")
    tdata.obs["ph"] *= 10
    pd.testing.assert_frame_equal(bt.da.linda(tdata, "ph"), out, rtol=1e-9)
    assert (out["contrast"] == "ph").all()


def test_no_feature_is_dropped():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[:, 4] = 0
    tdata.X = sp.csr_matrix(dense)
    out = bt.da.linda(tdata, "group")
    assert out.index.tolist() == tdata.var_names.tolist()
    # As in MicrobiomeStat: the all-zero feature is 0.5 everywhere, so it is tested like the others.
    assert np.isfinite(out[["effect", "se", "pvalue", "qvalue"]].to_numpy()).all()


def test_linda_keeps_input(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.da.linda(tdata, "group")
    assert_unchanged(before, tdata)
