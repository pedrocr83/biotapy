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


def _row_sums(M) -> np.ndarray:
    return np.asarray(M.sum(axis=1)).ravel()


def test_relative_rows_sum_to_one():
    np.testing.assert_allclose(_row_sums(bt.pp.relative(bt.datasets.toy()).layers["relative"]), 1.0)


def test_relative_keeps_x_and_input(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    out = bt.pp.relative(tdata)
    assert_unchanged(before, tdata)
    assert (out.X != tdata.X).nnz == 0


def test_relative_all_zero_sample_stays_zero():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[0] = 0
    tdata.X = sp.csr_matrix(dense)
    assert bt.pp.relative(tdata).layers["relative"][0].nnz == 0


def test_relative_single_sample():
    tdata = bt.datasets.toy()[:1].copy()
    np.testing.assert_allclose(_row_sums(bt.pp.relative(tdata).layers["relative"]), 1.0)


def test_relative_all_zero_feature_stays_zero():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[:, 0] = 0
    tdata.X = sp.csr_matrix(dense)
    relative = bt.pp.relative(tdata).layers["relative"]
    assert relative[:, 0].nnz == 0
    np.testing.assert_allclose(_row_sums(relative), 1.0)


def test_relative_records_provenance():
    entries = bt.pp.relative(bt.datasets.toy()).uns["biotapy"]["provenance"]
    assert json.loads(entries[-1])["step"] == "pp.relative"


def test_relative_subnormal_total_stays_finite():
    # 1 / 5e-324 overflows to inf; dividing each value by the total does not.
    tiny = np.nextafter(0.0, 1.0)
    adata = ad.AnnData(
        X=sp.csr_matrix(np.array([[tiny, 0.0], [1.0, 3.0]])),
        obs=pd.DataFrame(index=["s0", "s1"]),
        var=pd.DataFrame(index=["f0", "f1"]),
    )
    rel = bt.pp.relative(adata).layers["relative"].toarray()
    np.testing.assert_array_equal(rel, [[1.0, 0.0], [0.25, 0.75]])


@given(arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 8)), elements=st.integers(0, 1000)))
def test_relative_rows_sum_to_one_or_zero(dense):
    adata = ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
    )
    expected = np.where(dense.sum(axis=1) > 0, 1.0, 0.0)
    np.testing.assert_allclose(_row_sums(bt.pp.relative(adata).layers["relative"]), expected)


def test_relative_float32_rows_sum_to_one_in_float64():
    adata = ad.AnnData(
        X=sp.csr_matrix(np.array([[0.1, 0.2, 0.3, 0.7], [1.1, 2.3, 0.0, 4.9]]), dtype=np.float32),
        obs=pd.DataFrame(index=["s0", "s1"]),
        var=pd.DataFrame(index=list("abcd")),
    )
    np.testing.assert_allclose(_row_sums(bt.pp.relative(adata).layers["relative"]), 1.0, rtol=0, atol=1e-12)


def _clr_by_hand(dense: np.ndarray, pseudocount: float) -> np.ndarray:
    logs = np.log(dense + pseudocount)
    return logs - logs.mean(axis=1, keepdims=True)


def test_clr_is_the_log_minus_the_mean_log(make_adata):
    dense = np.array([[1.0, 3.0, 0.0], [0.0, 2.0, 6.0]])
    np.testing.assert_allclose(bt.pp.clr(make_adata(dense)).layers["clr"], _clr_by_hand(dense, 0.5), rtol=1e-12)


def test_clr_is_dense_float64_and_rows_sum_to_zero():
    out = bt.pp.clr(bt.datasets.toy())
    assert isinstance(out.layers["clr"], np.ndarray) and out.layers["clr"].dtype == np.float64
    np.testing.assert_allclose(out.layers["clr"].sum(axis=1), 0.0, atol=1e-12)


def test_clr_keeps_x_and_input(assert_unchanged):
    tdata = bt.datasets.toy()
    before = tdata.copy()
    out = bt.pp.clr(tdata)
    assert_unchanged(before, tdata)
    assert (out.X != tdata.X).nnz == 0 and type(out) is type(tdata)


def test_clr_all_zero_sample_is_all_zero(make_adata):
    out = bt.pp.clr(make_adata(np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])))
    np.testing.assert_array_equal(out.layers["clr"][0], 0.0)


def test_clr_all_zero_feature_is_finite(make_adata):
    out = bt.pp.clr(make_adata(np.array([[0.0, 4.0, 1.0], [0.0, 2.0, 3.0]])))
    assert np.isfinite(out.layers["clr"]).all()
    np.testing.assert_allclose(out.layers["clr"].sum(axis=1), 0.0, atol=1e-12)


def test_clr_single_sample():
    out = bt.pp.clr(bt.datasets.toy()[:1].copy())
    assert out.layers["clr"].shape == (1, 8)


def test_clr_without_pseudocount_needs_no_zeros(make_adata):
    dense = np.array([[1.0, 2.0], [3.0, 5.0]])
    np.testing.assert_allclose(
        bt.pp.clr(make_adata(dense), pseudocount=0).layers["clr"], _clr_by_hand(dense, 0.0), rtol=1e-12
    )


def test_clr_zero_without_pseudocount_raises(make_adata):
    with pytest.raises(ValueError, match="pass pseudocount > 0 to pp.clr"):
        bt.pp.clr(make_adata(np.array([[0.0, 2.0], [3.0, 5.0]])), pseudocount=0)


@pytest.mark.parametrize("pseudocount", [-0.5, np.nan, np.inf])
def test_clr_bad_pseudocount_raises(pseudocount):
    with pytest.raises(ValueError, match="pseudocount must be a finite number >= 0"):
        bt.pp.clr(bt.datasets.toy(), pseudocount=pseudocount)


@pytest.mark.parametrize("bad", [-1.0, np.nan])
def test_clr_negative_or_missing_value_raises(make_adata, bad):
    with pytest.raises(ValueError, match="pp.clr needs finite, non-negative values in X"):
        bt.pp.clr(make_adata(np.array([[bad, 2.0], [3.0, 5.0]])))


def test_clr_pseudocount_above_the_smallest_value_warns():
    relative = bt.pp.relative(bt.datasets.toy())
    relative.X = relative.layers["relative"]
    with pytest.warns(UserWarning, match=r"pseudocount=0.5 is larger than the smallest non-zero value in X \(0.013\)"):
        bt.pp.clr(relative)


def test_clr_on_counts_does_not_warn(recwarn):
    bt.pp.clr(bt.datasets.toy())
    assert not [w for w in recwarn if issubclass(w.category, UserWarning)]


def test_clr_records_provenance():
    entries = bt.pp.clr(bt.datasets.toy(), pseudocount=1).uns["biotapy"]["provenance"]
    assert json.loads(entries[-1]) == {"step": "pp.clr", "version": bt.__version__, "params": {"pseudocount": 1}}


@given(
    arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 1000)),
    # At most 1, the smallest non-zero count, so no call warns.
    st.floats(0.01, 1),
    st.floats(0.01, 100),
)
def test_clr_rows_sum_to_zero_and_ignore_scale(dense, pseudocount, scale):
    adata = ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
    )
    out = bt.pp.clr(adata, pseudocount=pseudocount).layers["clr"]
    np.testing.assert_allclose(out.sum(axis=1), 0.0, atol=1e-9)
    # CLR is scale invariant: scaling X and the pseudocount together changes nothing.
    scaled = adata.copy()
    scaled.X = sp.csr_matrix(dense * scale)
    np.testing.assert_allclose(
        bt.pp.clr(scaled, pseudocount=pseudocount * scale).layers["clr"], out, rtol=1e-9, atol=1e-9
    )


@pytest.mark.parametrize("pseudocount", [True, False, "0.5", None])
def test_clr_non_numeric_pseudocount_raises(pseudocount):
    with pytest.raises(TypeError, match="pseudocount must be a real number"):
        bt.pp.clr(bt.datasets.toy(), pseudocount=pseudocount)
