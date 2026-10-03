import json

import anndata as ad
import numpy as np
import pandas as pd
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
