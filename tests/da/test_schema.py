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

METHODS = [bt.da.ancombc2, bt.da.linda]
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
    # Distinct counts: no empty sample and no feature fitted exactly, where ANCOM-BC2's variances are 0 and it fails.
    counts=arrays(np.int64, (6, 7), elements=st.integers(0, 200), unique=True),
    split=st.integers(2, 4),
)
def test_direction_is_the_sign_and_qvalue_bounds_pvalue(method, counts, split):
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


def _broken(change):
    table = bt.da.linda(bt.datasets.toy(), "group")
    change(table)
    return table


def _drop_effect(table, feature):
    table.loc[feature, ["effect", "direction"]] = [np.nan, 0]


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda t: t.drop(columns="se", inplace=True), r"lacks the result columns \['se'\]"),
        (lambda t: t.__setitem__("direction", t["direction"].astype(float)), "direction an integer"),
        (lambda t: t.__setitem__("pvalue", t["pvalue"] * 10), "must lie between 0 and 1"),
        (lambda t: t.__setitem__("qvalue", t["qvalue"].where(t.index != "f1")), "NaN exactly where pvalue is"),
        (lambda t: t.__setitem__("direction", -t["direction"]), "direction must be the sign of effect"),
        (lambda t: t.__setitem__("contrast", ["B vs A"] * 7 + ["A vs B"]), "must hold one contrast"),
        (lambda t: t.__setitem__("contrast", [np.nan] + ["B vs A"] * 7), "must hold one contrast"),
        (lambda t: _drop_effect(t, "f1"), "effect must be NaN exactly where pvalue is"),
    ],
    ids=["column", "dtype", "range", "missing q", "direction", "contrast", "nan contrast", "missing effect"],
)
def test_consensus_refuses_tables_that_break_the_schema(change, message):
    table = _broken(change)
    with pytest.raises(ValueError, match=message):
        bt.da.consensus([table, bt.da.ancombc2(bt.datasets.toy(), "group")])


def test_consensus_refuses_repeated_features():
    table = bt.da.linda(bt.datasets.toy(), "group")
    with pytest.raises(ValueError, match=r"results\[0\] repeats features \['f1'\]"):
        bt.da.consensus([table.set_axis(["f1"] * 8), bt.da.ancombc2(bt.datasets.toy(), "group")])


def test_consensus_refuses_what_is_not_a_table():
    with pytest.raises(TypeError, match=r"results\[1\] must be a result table of a bt.da method, got str"):
        bt.da.consensus([bt.da.linda(bt.datasets.toy(), "group"), "ancombc2"])
