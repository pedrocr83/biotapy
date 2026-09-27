import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import treedata as td
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist, squareform

import biotapy as bt


def _toy_with(metric="braycurtis"):
    tdata = bt.datasets.toy()
    bt.tl.beta(tdata, metric=metric, inplace=True)
    return tdata


def test_pcoa_reproduces_euclidean_distances():
    # PCoA of Euclidean distances recovers them exactly from all positive axes.
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray().astype(float)
    tdata.obsp["euclidean"] = squareform(pdist(dense))
    coords, axes = bt.tl.pcoa(tdata, distance="euclidean")
    assert list(coords.columns) == [f"PC{i}" for i in range(1, 6)]  # at most n_obs - 1 axes
    np.testing.assert_allclose(squareform(pdist(coords.to_numpy())), tdata.obsp["euclidean"], atol=1e-9)
    assert axes["proportion_explained"].sum() == pytest.approx(1.0)


def test_pcoa_proportion_divides_by_all_eigenvalues():
    coords, axes = bt.tl.pcoa(_toy_with(), n_components=2)
    assert coords.shape == (6, 2) and list(axes.index) == ["PC1", "PC2"]
    assert (axes["proportion_explained"] <= 1).all() and axes["eigenvalue"].is_monotonic_decreasing


def test_pcoa_inplace_writes_obsm_and_uns():
    tdata = _toy_with()
    coords, axes = bt.tl.pcoa(tdata, n_components=3)
    assert bt.tl.pcoa(tdata, n_components=3, inplace=True) is None
    np.testing.assert_array_equal(tdata.obsm["X_pcoa"], coords.to_numpy())
    np.testing.assert_array_equal(tdata.uns["biotapy"]["pcoa"]["eigenvalues"], axes["eigenvalue"].to_numpy())
    np.testing.assert_array_equal(
        tdata.uns["biotapy"]["pcoa"]["proportion_explained"], axes["proportion_explained"].to_numpy()
    )


def test_pcoa_input_unchanged(assert_unchanged):
    tdata = _toy_with()
    before = tdata.copy()
    bt.tl.pcoa(tdata)
    assert_unchanged(before, tdata)


def test_pcoa_missing_distance_names_the_call():
    with pytest.raises(KeyError, match=r"bt.tl.unifrac\(tdata, weighted=True, inplace=True\)"):
        bt.tl.pcoa(bt.datasets.toy(), distance="weighted_unifrac")


def test_pcoa_nan_distances_raise():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[[0, 1]] = 0  # two all-zero samples are NaN apart under Bray-Curtis
    tdata.X = sp.csr_matrix(dense)
    bt.tl.beta(tdata, inplace=True)
    with pytest.raises(ValueError, match="NaN"):
        bt.tl.pcoa(tdata)


def test_pcoa_single_sample_and_bad_components_raise():
    with pytest.raises(ValueError, match="at least 2 samples"):
        bt.tl.pcoa(_toy_with()[:1].copy())
    with pytest.raises(ValueError, match="n_components"):
        bt.tl.pcoa(_toy_with(), n_components=0)


def test_pcoa_after_filter_samples_uses_the_subset():
    tdata = bt.pp.filter_samples(_toy_with(), 70)
    coords, _ = bt.tl.pcoa(tdata)
    assert list(coords.index) == ["s3", "s4", "s5", "s6"] and coords.shape[1] == 3


def test_nmds_same_seed_same_result():
    first, stress = bt.tl.nmds(_toy_with(), seed=3)
    again, stress_again = bt.tl.nmds(_toy_with(), seed=np.random.default_rng(3))
    np.testing.assert_array_equal(first.to_numpy(), again.to_numpy())
    assert stress == stress_again and 0 <= stress < 0.2
    assert list(first.columns) == ["NMDS1", "NMDS2"]


def test_nmds_recovers_a_planar_configuration():
    points = np.random.default_rng(0).uniform(size=(10, 2))
    adata = ad.AnnData(X=sp.csr_matrix(np.ones((10, 1))), obs=pd.DataFrame(index=[f"s{i}" for i in range(10)]))
    adata.obsp["euclidean"] = squareform(pdist(points))
    coords, stress = bt.tl.nmds(adata, distance="euclidean", seed=0)
    _, _, disparity = procrustes(points, coords.to_numpy())
    assert stress < 0.05 and np.sqrt(1 - disparity) > 0.99


def test_nmds_inplace_writes_obsm_and_stress():
    tdata = _toy_with()
    coords, stress = bt.tl.nmds(tdata, seed=0)
    assert bt.tl.nmds(tdata, seed=0, inplace=True) is None
    np.testing.assert_array_equal(tdata.obsm["X_nmds"], coords.to_numpy())
    assert tdata.uns["biotapy"]["nmds"] == {"stress": stress}


def test_nmds_input_unchanged(assert_unchanged):
    tdata = _toy_with()
    before = tdata.copy()
    bt.tl.nmds(tdata, seed=0)
    assert_unchanged(before, tdata)


def test_nmds_nan_distances_raise():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[[0, 1]] = 0  # two all-zero samples are NaN apart under Bray-Curtis
    tdata.X = sp.csr_matrix(dense)
    bt.tl.beta(tdata, inplace=True)
    with pytest.raises(ValueError, match=r"distance='braycurtis': obsp\['braycurtis'\] holds NaN"):
        bt.tl.nmds(tdata, seed=0)


def test_nmds_too_few_samples_raise():
    with pytest.raises(ValueError, match="more than n_components"):
        bt.tl.nmds(_toy_with()[:3].copy(), seed=0)


def test_feature_changes_drop_stored_ordinations():
    tdata = _toy_with()
    bt.tl.pcoa(tdata, inplace=True)
    bt.tl.nmds(tdata, seed=0, inplace=True)
    out = bt.pp.filter_features(tdata, min_total=19)
    assert not {"X_pcoa", "X_nmds"} & set(out.obsm.keys())
    assert set(out.uns["biotapy"]) == {"x_kind", "provenance"}


def test_stored_results_survive_h5td(tmp_path):
    tdata = _toy_with()
    bt.tl.alpha(tdata, inplace=True)
    bt.tl.unifrac(tdata, weighted=True, inplace=True)
    bt.tl.pcoa(tdata, inplace=True)
    bt.tl.nmds(tdata, seed=0, inplace=True)
    tdata.write_h5td(tmp_path / "toy.h5td")
    back = td.read_h5td(tmp_path / "toy.h5td")
    np.testing.assert_array_equal(back.obsp["weighted_unifrac"], tdata.obsp["weighted_unifrac"])
    np.testing.assert_array_equal(back.obsm["X_pcoa"], tdata.obsm["X_pcoa"])
    np.testing.assert_array_equal(
        back.uns["biotapy"]["pcoa"]["eigenvalues"], tdata.uns["biotapy"]["pcoa"]["eigenvalues"]
    )
    assert back.uns["biotapy"]["nmds"]["stress"] == tdata.uns["biotapy"]["nmds"]["stress"]
    np.testing.assert_array_equal(back.obs["alpha_shannon"], tdata.obs["alpha_shannon"])
