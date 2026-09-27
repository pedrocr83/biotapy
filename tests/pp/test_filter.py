import json

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


def _adata(dense) -> ad.AnnData:
    return ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
    )


def _tips(tdata) -> set[str]:
    tree = get_tree(tdata)
    return {n for n in tree.nodes if tree.out_degree(n) == 0}


# toy(): f1..f8 are non-zero in 5, 5, 6, 6, 4, 5, 5, 4 of 6 samples and total 33, 19, 75, 117, 7, 130, 50, 11.


def test_min_prevalence_is_inclusive():
    # f5 and f8 are in exactly 4 of 6 samples.
    assert bt.pp.filter_features(bt.datasets.toy(), min_prevalence=4 / 6).n_vars == 8
    assert list(bt.pp.filter_features(bt.datasets.toy(), min_prevalence=5 / 6).var_names) == [
        "f1",
        "f2",
        "f3",
        "f4",
        "f6",
        "f7",
    ]


def test_min_total_is_inclusive():
    out = bt.pp.filter_features(bt.datasets.toy(), min_total=19)
    assert list(out.var_names) == ["f1", "f2", "f3", "f4", "f6", "f7"]


def test_both_thresholds_must_pass():
    out = bt.pp.filter_features(bt.datasets.toy(), min_prevalence=1.0, min_total=100)
    assert list(out.var_names) == ["f4"]


def test_decimal_prevalence_boundary_is_kept():
    # 7 / 25 == 0.28 in floating point, while 0.28 * 25 is just above 7.
    dense = np.ones((25, 2), dtype=np.int64)
    dense[7:, 0] = 0
    assert list(bt.pp.filter_features(_adata(dense), min_prevalence=0.28).var_names) == ["f0", "f1"]


def test_explicit_zeros_do_not_count_as_present():
    adata = _adata(np.array([[1, 3], [2, 4]]))
    adata.X.data[0] = 0  # a stored zero: still 4 stored entries
    assert adata.X.nnz == 4
    assert list(bt.pp.filter_features(adata, min_prevalence=1.0).var_names) == ["f1"]


def test_all_zero_feature_and_sample():
    dense = np.array([[0, 0, 0], [1, 0, 2], [3, 0, 0]])
    out = bt.pp.filter_features(_adata(dense), min_prevalence=0.1)
    assert list(out.var_names) == ["f0", "f2"] and out.n_obs == 3


def test_single_sample():
    assert list(bt.pp.filter_features(_adata(np.array([[0, 4, 1]])), min_total=1).var_names) == ["f1", "f2"]


def test_keeps_missing_ranks_and_prunes_the_tree():
    out = bt.pp.filter_features(bt.datasets.toy(), min_total=11)
    assert "f8" in out.var_names and pd.isna(out.var.loc["f8", "genus"])
    assert _tips(out) == set(out.var_names) == {"f1", "f2", "f3", "f4", "f6", "f7", "f8"}


def test_drops_derived_slots_and_records_provenance():
    tdata = bt.pp.relative(bt.datasets.toy())
    tdata.obsp["braycurtis"] = np.zeros((6, 6))
    out = bt.pp.filter_features(tdata, min_total=19)
    assert "relative" not in out.layers and "braycurtis" not in out.obsp
    entry = json.loads(out.uns["biotapy"]["provenance"][-1])
    assert entry["step"] == "pp.filter_features" and entry["params"] == {"min_prevalence": None, "min_total": 19}


def test_numpy_thresholds_are_recorded_as_python_numbers():
    out = bt.pp.filter_features(bt.datasets.toy(), min_prevalence=np.float32(0.5), min_total=np.int64(19))
    assert list(out.var_names) == ["f1", "f2", "f3", "f4", "f6", "f7"]
    assert json.loads(out.uns["biotapy"]["provenance"][-1])["params"] == {"min_prevalence": 0.5, "min_total": 19}


def test_input_unchanged(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.pp.filter_features(tdata, min_prevalence=0.9)
    assert_unchanged(before, tdata)


def test_neither_threshold_raises():
    with pytest.raises(ValueError, match="min_prevalence=, min_total="):
        bt.pp.filter_features(bt.datasets.toy())


def test_prevalence_outside_zero_to_one_raises():
    with pytest.raises(ValueError, match="min_prevalence"):
        bt.pp.filter_features(bt.datasets.toy(), min_prevalence=5)


def test_nothing_passing_raises():
    with pytest.raises(ValueError, match="no feature passes"):
        bt.pp.filter_features(bt.datasets.toy(), min_total=10_000)


@given(
    arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 8)), elements=st.integers(0, 20)),
    st.sampled_from([0.0, 0.25, 0.5, 1.0]),
    st.integers(0, 30),
)
def test_keeps_exactly_the_features_passing_both(dense, prevalence, total):
    passing = ((dense > 0).mean(axis=0) >= prevalence) & (dense.sum(axis=0) >= total)
    if not passing.any():
        with pytest.raises(ValueError, match="no feature passes"):
            bt.pp.filter_features(_adata(dense), min_prevalence=prevalence, min_total=total)
        return
    out = bt.pp.filter_features(_adata(dense), min_prevalence=prevalence, min_total=total)
    assert list(out.var_names) == [f"f{i}" for i in np.flatnonzero(passing)]
    np.testing.assert_array_equal(out.X.toarray(), dense[:, passing])


# toy() sample depths: s1..s6 = 68, 66, 77, 76, 82, 73.


def test_min_depth_is_inclusive():
    assert list(bt.pp.filter_samples(bt.datasets.toy(), 73).obs_names) == ["s3", "s4", "s5", "s6"]


def test_filter_samples_keeps_every_slot():
    tdata = bt.pp.relative(bt.datasets.toy())
    tdata.obsp["braycurtis"] = np.arange(36.0).reshape(6, 6)
    tdata.obsm["X_pcoa"] = np.arange(12.0).reshape(6, 2)
    out = bt.pp.filter_samples(tdata, 73)
    np.testing.assert_array_equal(out.obsp["braycurtis"], np.arange(36.0).reshape(6, 6)[2:, 2:])
    np.testing.assert_array_equal(out.obsm["X_pcoa"], np.arange(12.0).reshape(6, 2)[2:])
    assert "relative" in out.layers and out.n_vars == 8 and _tips(out) == set(out.var_names)
    assert json.loads(out.uns["biotapy"]["provenance"][-1])["step"] == "pp.filter_samples"


def test_filter_samples_drops_an_all_zero_sample_and_keeps_all_zero_features():
    out = bt.pp.filter_samples(_adata(np.array([[0, 0], [0, 3]])), 1)
    assert list(out.obs_names) == ["s1"] and out.n_vars == 2


def test_filter_samples_accepts_a_numpy_threshold():
    # rarefy's X is int64, so a depth taken from it is an np.int64, not an int.
    tdata = bt.pp.rarefy(bt.datasets.toy(), depth=60, seed=0)
    depths = np.asarray(tdata.X.sum(axis=1)).ravel()
    out = bt.pp.filter_samples(tdata, depths.min())
    assert out.n_obs == 6
    assert json.loads(out.uns["biotapy"]["provenance"][-1])["params"] == {"min_depth": 60}


def test_filter_samples_single_sample():
    assert bt.pp.filter_samples(_adata(np.array([[2, 3]])), 5).n_obs == 1


def test_filter_samples_input_unchanged(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.pp.filter_samples(tdata, 73)
    assert_unchanged(before, tdata)


def test_filter_samples_nothing_passing_raises():
    with pytest.raises(ValueError, match="min_depth=1000"):
        bt.pp.filter_samples(bt.datasets.toy(), 1000)


@given(
    arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 5)), elements=st.integers(0, 10)), st.integers(0, 20)
)
def test_filter_samples_keeps_exactly_the_deep_samples(dense, depth):
    deep = dense.sum(axis=1) >= depth
    if not deep.any():
        with pytest.raises(ValueError, match="min_depth"):
            bt.pp.filter_samples(_adata(dense), depth)
        return
    assert list(bt.pp.filter_samples(_adata(dense), depth).obs_names) == [f"s{i}" for i in np.flatnonzero(deep)]
