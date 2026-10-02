import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import biotapy as bt


def _toy_with(*metrics):
    tdata = bt.datasets.toy()
    bt.tl.alpha(tdata, metrics=list(metrics), inplace=True)
    return tdata


def _points(ax):
    return np.concatenate([collection.get_offsets() for collection in ax.collections])


def test_richness_one_point_per_sample(ax):
    tdata = _toy_with("shannon")
    bt.pl.richness(tdata, "shannon", ax=ax)
    np.testing.assert_array_equal(_points(ax), np.c_[np.arange(6), tdata.obs["alpha_shannon"]])
    assert [label.get_text() for label in ax.get_xticklabels()] == list(tdata.obs_names)
    assert ax.get_ylabel() == "shannon"


def test_richness_x_places_samples_by_group(ax):
    tdata = _toy_with("observed_features")
    bt.pl.richness(tdata, "observed_features", x="group", ax=ax)
    np.testing.assert_array_equal(_points(ax)[:, 0], [0, 0, 0, 1, 1, 1])
    assert [label.get_text() for label in ax.get_xticklabels()] == ["A", "B"]


def test_richness_color_draws_one_scatter_per_group(ax):
    tdata = _toy_with("shannon")
    tdata.obs["site"] = pd.Categorical(["u", "v", None, "u", "v", "u"])
    bt.pl.richness(tdata, "shannon", x="group", color="site", ax=ax)
    assert [text.get_text() for text in ax.get_legend().get_texts()] == ["u", "v", "NA"]
    assert [len(collection.get_offsets()) for collection in ax.collections] == [3, 2, 1]
    np.testing.assert_array_equal(np.sort(_points(ax)[:, 1]), np.sort(tdata.obs["alpha_shannon"]))


def test_richness_nan_values_are_left_out(ax):
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[0] = 0
    tdata.X = sp.csr_matrix(dense)
    bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)  # s1 is all-zero: Shannon NaN
    bt.pl.richness(tdata, "shannon", color="group", ax=ax)
    assert len(_points(ax)) == 5 and np.isfinite(_points(ax)).all()


def test_richness_single_sample(ax):
    tdata = bt.datasets.toy()[:1].copy()
    bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)
    bt.pl.richness(tdata, "shannon", ax=ax)
    np.testing.assert_array_equal(_points(ax), [[0, tdata.obs["alpha_shannon"].iloc[0]]])


def test_richness_input_unchanged(assert_unchanged, ax):
    tdata = _toy_with("shannon")
    before = tdata.copy()
    bt.pl.richness(tdata, "shannon", x="group", color="group", ax=ax)
    assert_unchanged(before, tdata)


def test_richness_missing_metric_names_the_call(ax):
    with pytest.raises(KeyError, match=r"bt.tl.alpha\(adata, metrics=\['chao1'\], inplace=True\)"):
        bt.pl.richness(_toy_with("shannon"), "chao1", ax=ax)


def test_richness_bad_columns_raise(ax):
    tdata = _toy_with("shannon")
    with pytest.raises(KeyError, match="color='site' is not a column of obs"):
        bt.pl.richness(tdata, "shannon", color="site", ax=ax)
    with pytest.raises(TypeError, match="x='alpha_shannon' is a numeric column"):
        bt.pl.richness(tdata, "shannon", x="alpha_shannon", ax=ax)
