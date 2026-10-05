import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays
from scipy.stats import false_discovery_control

import biotapy as bt

METHODS = [bt.da.linda]
COLUMNS = ["effect", "se", "pvalue", "qvalue", "direction", "method", "contrast"]


@pytest.mark.parametrize("method", METHODS)
def test_result_has_the_schema_columns_and_dtypes(method):
    out = method(bt.datasets.toy(), "group")
    assert out.columns.tolist() == COLUMNS and out.index.name == "feature"
    assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
    assert (out.dtypes.iloc[:4] == np.float64).all() and out["direction"].dtype == np.int8
    assert pd.api.types.is_string_dtype(out["method"]) and pd.api.types.is_string_dtype(out["contrast"])
    assert (out["method"] == method.__name__).all()


@pytest.mark.parametrize("method", METHODS)
def test_qvalue_is_benjamini_hochberg(method):
    out = method(bt.datasets.toy(), "group")
    tested = out["pvalue"].notna()
    np.testing.assert_allclose(
        out.loc[tested, "qvalue"], false_discovery_control(out.loc[tested, "pvalue"]), rtol=1e-12
    )


@pytest.mark.parametrize("method", METHODS)
@settings(max_examples=25, deadline=None)
@given(
    counts=arrays(np.int64, (6, 7), elements=st.integers(0, 50)),
    split=st.integers(2, 4),
)
def test_direction_is_the_sign_and_qvalue_bounds_pvalue(method, counts, split):
    counts[:, 0] += 1  # no empty sample
    adata = ad.AnnData(
        X=sp.csr_matrix(counts),
        obs=pd.DataFrame({"g": ["a"] * split + ["b"] * (6 - split)}, index=[f"s{i}" for i in range(6)]),
        var=pd.DataFrame(index=[f"f{i}" for i in range(7)]),
    )
    out = method(adata, "g")
    tested = out["pvalue"].notna()
    assert out.index.tolist() == adata.var_names.tolist()
    assert (out["direction"] == np.sign(out["effect"].fillna(0))).all()
    assert (out.loc[tested, "qvalue"] >= out.loc[tested, "pvalue"] - 1e-15).all()
    assert out.loc[tested, ["pvalue", "qvalue"]].stack().between(0, 1).all()
    assert out.loc[~tested, ["effect", "qvalue"]].isna().all().all()
