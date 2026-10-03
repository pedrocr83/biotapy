import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt
from biotapy._core import make_function_mudata


def _row_totals(adata):
    return np.asarray(adata.X.sum(axis=1)).ravel()


def test_community_rows_sum_to_one_per_sample():
    out = bt.fn.renorm(bt.datasets.toy_humann(), "relab")
    np.testing.assert_allclose(_row_totals(out["function"]), 1.0)
    assert out["function"].uns["biotapy"]["x_kind"] == "relative"


def test_strata_are_divided_by_the_community_total():
    mdata = bt.datasets.toy_humann()
    out = bt.fn.renorm(mdata, "cpm")
    totals = _row_totals(mdata["function"])
    np.testing.assert_allclose(
        out["function_by_taxon"].X.toarray(), mdata["function_by_taxon"].X.toarray() / totals[:, None] * 1e6
    )
    assert out["function_by_taxon"].uns["biotapy"]["x_kind"] == "cpm"


def test_special_false_drops_specials_before_totalling():
    out = bt.fn.renorm(bt.datasets.toy_humann(), "relab", special=False)
    assert not out["function"].var["special"].any() and not out["function_by_taxon"].var["special"].any()
    np.testing.assert_allclose(_row_totals(out["function"]), 1.0)


def test_keeps_global_obs_and_other_modalities():
    mdata = bt.datasets.toy_humann()
    mdata.obs["subject"] = ["a", "b", "c", "d", "e", "f"]
    out = bt.fn.renorm(mdata, "relab")
    assert out.obs["subject"].tolist() == ["a", "b", "c", "d", "e", "f"]
    assert '"step": "fn.renorm"' in out["function"].uns["biotapy"]["provenance"][-1]


def test_input_unchanged(assert_unchanged):
    mdata = bt.datasets.toy_humann()
    before = mdata.copy()
    bt.fn.renorm(mdata, "relab", special=False)
    for key in ("function", "function_by_taxon"):
        assert_unchanged(before[key], mdata[key])


def test_all_zero_sample_stays_zero():
    mdata = bt.datasets.toy_humann()
    for key in ("function", "function_by_taxon"):
        dense = mdata[key].X.toarray()
        dense[0] = 0
        mdata.mod[key].X = sp.csr_matrix(dense)
    out = bt.fn.renorm(mdata, "relab")
    assert out["function"].X[0].nnz == 0 and np.isfinite(out["function_by_taxon"].X.toarray()).all()


def test_single_sample():
    out = bt.fn.renorm(bt.datasets.toy_humann()[:1].copy(), "relab")
    np.testing.assert_allclose(_row_totals(out["function"]), [1.0])


def test_unknown_units_raise():
    with pytest.raises(ValueError, match="units="):
        bt.fn.renorm(bt.datasets.toy_humann(), "tpm")


def test_missing_modality_names_it():
    mdata = bt.datasets.toy_humann()
    del mdata.mod["function_by_taxon"]
    with pytest.raises(KeyError, match="function_by_taxon"):
        bt.fn.renorm(mdata, "relab")


def test_stratified_only_table_raises():
    mdata = bt.datasets.toy_humann()
    mdata.mod["function"] = mdata["function"][:, []].copy()
    with pytest.raises(ValueError, match="community"):
        bt.fn.renorm(mdata, "relab")


@given(arrays(np.float64, st.tuples(st.integers(1, 4), st.integers(1, 5)), elements=st.floats(0, 1e6)))
def test_community_totals_are_one_or_zero(dense):
    obs = pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])])
    ids = pd.Index([f"K{j}" for j in range(dense.shape[1])] + [f"K{j}|unclassified" for j in range(dense.shape[1])])
    mdata = make_function_mudata(np.hstack([dense, dense]), obs=obs, row_ids=ids, x_kind="rpk", source="test")
    totals = _row_totals(bt.fn.renorm(mdata, "relab")["function"])
    np.testing.assert_allclose(totals, np.where(dense.sum(axis=1) > 0, 1.0, 0.0), rtol=1e-12)
