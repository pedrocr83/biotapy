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
    with pytest.warns(UserWarning, match="no community abundance"):
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


@st.composite
def _community_and_strata(draw):
    n_samples, n_features = draw(st.integers(1, 4)), draw(st.integers(1, 5))
    community = draw(arrays(np.float64, (n_samples, n_features), elements=st.floats(1, 1e6)))
    strata = draw(arrays(np.float64, (n_samples, n_features), elements=st.floats(0, 1e6)))
    return community, strata


@given(_community_and_strata())
def test_strata_are_raw_values_over_the_community_total(tables):
    community, strata = tables
    obs = pd.DataFrame(index=[f"s{i}" for i in range(community.shape[0])])
    n = community.shape[1]
    ids = pd.Index([f"K{j}" for j in range(n)] + [f"K{j}|unclassified" for j in range(n)])
    mdata = make_function_mudata(np.hstack([community, strata]), obs=obs, row_ids=ids, x_kind="rpk", source="test")
    out = bt.fn.renorm(mdata, "relab")
    np.testing.assert_allclose(_row_totals(out["function"]), 1.0, rtol=1e-12)
    totals = community.sum(axis=1, keepdims=True)
    np.testing.assert_allclose(out["function_by_taxon"].X.toarray(), strata / totals, rtol=1e-12)


def test_zero_total_sample_warns_once_naming_it():
    mdata = bt.datasets.toy_humann()
    for key in ("function", "function_by_taxon"):
        dense = mdata[key].X.toarray()
        dense[0] = 0
        mdata.mod[key].X = sp.csr_matrix(dense)
    name = mdata.obs_names[0]
    with pytest.warns(UserWarning, match=rf"1 sample.*{name}") as record:
        bt.fn.renorm(mdata, "relab")
    assert len(record) == 1 and record[0].filename == __file__


@pytest.mark.parametrize(
    "samples", [slice(None, None, -1), slice(0, 4)], ids=["same samples reordered", "fewer samples"]
)
def test_modalities_with_different_samples_raise_naming_mdata(samples):
    # Totals come from "function" and are applied by row to "function_by_taxon".
    mdata = bt.datasets.toy_humann()
    mdata.mod["function_by_taxon"] = mdata["function_by_taxon"][samples].copy()
    with pytest.raises(ValueError, match=r"mdata\['function'\] and mdata\['function_by_taxon'\].*same samples"):
        bt.fn.renorm(mdata, "relab")
