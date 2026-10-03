import json
import re
from pathlib import Path

import mudata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given, settings
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


# No deadline: under --cov, building a MuData per example can exceed Hypothesis's 200 ms default.
@settings(deadline=None)
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
    message = rf"1 sample.*{name}"
    with pytest.warns(UserWarning, match=message) as record:
        bt.fn.renorm(mdata, "relab")
    # Count only biotapy's warning: pre-release pandas makes anndata emit its own deprecation warnings.
    own = [w for w in record if issubclass(w.category, UserWarning) and re.search(message, str(w.message))]
    assert len(own) == 1 and own[0].filename == __file__


@pytest.mark.parametrize(
    "samples", [slice(None, None, -1), slice(0, 4)], ids=["same samples reordered", "fewer samples"]
)
def test_modalities_with_different_samples_raise_naming_mdata(samples):
    # Totals come from "function" and are applied by row to "function_by_taxon".
    mdata = bt.datasets.toy_humann()
    mdata.mod["function_by_taxon"] = mdata["function_by_taxon"][samples].copy()
    with pytest.raises(ValueError, match=r"mdata\['function'\] and mdata\['function_by_taxon'\].*same samples"):
        bt.fn.renorm(mdata, "relab")


def test_special_false_on_a_table_of_specials_says_so():
    # HUMAnN's demo pathway table holds only UNMAPPED and UNINTEGRATED community rows.
    mdata = bt.io.read_humann(Path(__file__).parents[1] / "data" / "humann" / "demo_pathabundance_with_names.tsv")
    with pytest.raises(ValueError, match=r"special=False dropped every community row"):
        bt.fn.renorm(mdata, "relab", special=False)


def _stage(adata):
    """(var columns and dtypes, x_kind, provenance steps) of one modality."""
    meta = adata.uns["biotapy"]
    steps = [json.loads(entry)["step"] for entry in meta["provenance"]]
    return dict(adata.var.dtypes.astype(str)), meta["x_kind"], steps


def test_toy_humann_composes_through_func_glom_renorm_and_func_glom_again():
    # The review's chain on HMP2, here on toy_humann: every stage keeps the var contract and adds one step.
    community = {"name": "str", "special": "bool"}
    by_taxon = dict.fromkeys(["function", "name", "taxon", "genus", "species"], "str") | {"special": "bool"}
    classes = pd.DataFrame(
        {"child": ["1.1.1.1", "2.7.1.1", "2.7.1.2"], "parent": ["1.-.-.-", "2.-.-.-", "2.-.-.-"], "level": "class"}
    )
    everything = pd.DataFrame({"child": ["1.-.-.-", "2.-.-.-"], "parent": "all", "level": "all"})
    mdata = bt.datasets.toy_humann()
    grouped = mudata.MuData({key: bt.fn.func_glom(mod, "class", hierarchy=classes) for key, mod in mdata.mod.items()})
    renormed = bt.fn.renorm(grouped, "relab")
    again = {key: bt.fn.func_glom(mod, "all", hierarchy=everything) for key, mod in renormed.mod.items()}
    stages = [
        (mdata.mod, "rpk", ["datasets.toy_humann"]),
        (grouped.mod, "rpk", ["datasets.toy_humann", "fn.func_glom"]),
        (renormed.mod, "relative", ["datasets.toy_humann", "fn.func_glom", "fn.renorm"]),
        (again, "relative", ["datasets.toy_humann", "fn.func_glom", "fn.renorm", "fn.func_glom"]),
    ]
    for mods, kind, steps in stages:
        assert _stage(mods["function"]) == (community, kind, steps)
        assert _stage(mods["function_by_taxon"]) == (by_taxon, kind, steps)
    assert again["function"].var_names.tolist() == ["UNGROUPED", "UNMAPPED", "all"]
