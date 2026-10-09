import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.stats import false_discovery_control
from skbio.stats.composition import ancombc2

import biotapy as bt


def _toy_without(feature, samples):
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[samples, tdata.var_names.get_loc(feature)] = 0
    tdata.X = sp.csr_matrix(dense)
    return tdata


def test_ancombc2_is_scikit_bio_in_log2():
    tdata = bt.datasets.toy()
    out = bt.da.ancombc2(tdata, "group")
    counts = pd.DataFrame(tdata.X.toarray(), index=tdata.obs_names, columns=tdata.var_names)
    expected = ancombc2(counts, tdata.obs, "group").result.xs("group[T.B]", level="Covariate")
    np.testing.assert_allclose(out["effect"], expected["Log(FC)"] / np.log(2), rtol=1e-12)
    np.testing.assert_allclose(out["se"], expected["SE"] / np.log(2), rtol=1e-12)
    np.testing.assert_allclose(out["pvalue"], expected["pvalue"], rtol=1e-12)


def test_feature_absent_from_a_group_is_not_tested():
    out = bt.da.ancombc2(_toy_without("f5", [0, 1, 2]), "group")  # no f5 read in group A
    assert out.loc["f5", ["effect", "se", "pvalue", "qvalue"]].isna().all() and out.loc["f5", "direction"] == 0
    tested = out.drop(index="f5")
    assert np.isfinite(tested[["effect", "se", "pvalue", "qvalue"]].to_numpy()).all()
    # scikit-bio and R count it with p = 1; biotapy leaves it out of the correction.
    np.testing.assert_allclose(tested["qvalue"], false_discovery_control(tested["pvalue"]), rtol=1e-12)


def test_all_zero_feature_is_not_tested():
    out = bt.da.ancombc2(_toy_without("f7", slice(None)), "group")
    assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
    assert out.loc["f7", ["effect", "pvalue"]].isna().all() and out["effect"].notna().sum() == 7


def test_ancombc2_keeps_input(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.da.ancombc2(tdata, "group")
    assert_unchanged(before, tdata)


def test_feature_with_no_residual_degrees_of_freedom_is_not_tested():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[:, tdata.var_names.get_loc("f8")] = [0, 4, 0, 0, 0, 3]  # two reads in all: as many as model terms
    tdata.X = sp.csr_matrix(dense)
    out = bt.da.ancombc2(tdata, "group")
    assert out.loc["f8", ["effect", "se", "pvalue", "qvalue"]].isna().all() and out.loc["f8", "direction"] == 0
    assert out.drop(index="f8")["effect"].notna().all()


def test_scikit_bio_failure_names_the_function():
    tdata = bt.datasets.toy()[:, ["f1", "f2"]].copy()
    tdata.X = sp.csr_matrix(
        np.array([[5, 0], [6, 0], [7, 0], [0, 5], [0, 6], [0, 7]])
    )  # each feature in one group only
    with pytest.raises(ValueError, match=r"da\.ancombc2: scikit-bio could not fit the model: .*estimable"):
        bt.da.ancombc2(tdata, "group")


# Two reference samples whose centred log abundances of f5 nearly agree (variance 3.3e-6): scikit-bio 0.7.4's bias
# E-M divides 0 by 0 where R's .bias_em sets the responsibility to 0. R fits this table (B vs A, log2): f0 2.002,
# f1 0.744, f2 0.088, f3 -1.315, f4 -0.218, f5 -0.325, f6 -0.605. When scikit-bio fixes it, this test fails:
# turn it into a check against those values.
EM_UNDERFLOW = [
    [11, 39, 115, 65, 122, 90, 12],
    [81, 41, 194, 22, 174, 192, 186],
    [86, 68, 125, 1, 103, 14, 42],
    [21, 24, 43, 20, 38, 80, 6],
    [134, 61, 155, 8, 84, 169, 31],
    [71, 17, 64, 28, 63, 0, 10],
]


def test_two_reference_samples_that_underflow_scikit_bio_raise_naming_it():
    counts = np.array(EM_UNDERFLOW)
    adata = ad.AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame({"g": ["a"] * 2 + ["b"] * 4}, index=[f"s{i}" for i in range(6)]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(7)]),
    )
    with pytest.raises(ValueError, match=r"da\.ancombc2: scikit-bio could not fit the model: Quantiles must be in"):
        bt.da.ancombc2(adata, "g")
