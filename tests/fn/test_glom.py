import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from hypothesis import given
from hypothesis import strategies as st

import biotapy as bt

EC = pd.DataFrame(
    {
        "child": ["1.1.1.1", "2.7.1.1", "2.7.1.2", "2.7.1.1", "2.7.1.2"],
        "parent": ["1.-.-.-", "2.-.-.-", "2.-.-.-", "kinase", "kinase"],
        "level": ["class", "class", "class", "role", "role"],
        "parent_name": ["Oxidoreductases", "Transferases", "Transferases", np.nan, np.nan],
    }
)


def _function():
    return bt.datasets.toy_humann()["function"]


def _by_taxon():
    return bt.datasets.toy_humann()["function_by_taxon"]


def _column(adata, name):
    return adata[:, name].X.toarray().ravel()


def test_sums_children_into_parents():
    out = bt.fn.func_glom(_function(), "class", hierarchy=EC)
    assert out.var_names.tolist() == ["1.-.-.-", "2.-.-.-", "UNGROUPED", "UNMAPPED"]
    np.testing.assert_array_equal(
        _column(out, "2.-.-.-"), _column(_function(), "2.7.1.1") + _column(_function(), "2.7.1.2")
    )


def test_many_to_many_counts_a_child_in_every_parent():
    hierarchy = pd.DataFrame({"child": ["2.7.1.1", "2.7.1.1"], "parent": ["P1", "P2"], "level": "pathway"})
    out = bt.fn.func_glom(_function(), "pathway", hierarchy=hierarchy)
    np.testing.assert_array_equal(_column(out, "P1"), _column(_function(), "2.7.1.1"))
    np.testing.assert_array_equal(_column(out, "P2"), _column(_function(), "2.7.1.1"))


def test_unmapped_features_go_to_ungrouped_and_protected_specials_pass_through():
    out = bt.fn.func_glom(_function(), "class", hierarchy=EC)
    # 3.2.1.4 has no parent here; the input's own UNGROUPED row joins it.
    np.testing.assert_array_equal(
        _column(out, "UNGROUPED"), _column(_function(), "UNGROUPED") + _column(_function(), "3.2.1.4")
    )
    np.testing.assert_array_equal(_column(out, "UNMAPPED"), _column(_function(), "UNMAPPED"))
    assert out.var["special"].tolist() == [False, False, True, True]


def test_parent_names_become_var_name():
    out = bt.fn.func_glom(_function(), "class", hierarchy=EC)
    assert out.var["name"].tolist()[:2] == ["Oxidoreductases", "Transferases"] and out.var["name"].isna().tolist()[
        2:
    ] == [True, True]


def test_stratified_input_is_grouped_per_taxon():
    out = bt.fn.func_glom(_by_taxon(), "role", hierarchy=EC)
    assert out.var_names.tolist() == [
        "UNGROUPED|g__Bacteroides.s__Bacteroides_ovatus",
        "UNGROUPED|g__Faecalibacterium.s__Faecalibacterium_prausnitzii",
        "UNGROUPED|unclassified",
        "kinase|g__Bacteroides.s__Bacteroides_ovatus",
        "kinase|g__Blautia.s__Blautia_obeum",
    ]
    assert out.var.loc["kinase|g__Blautia.s__Blautia_obeum", ["function", "taxon", "genus"]].tolist() == [
        "kinase",
        "g__Blautia.s__Blautia_obeum",
        "Blautia",
    ]
    by_taxon = _by_taxon()
    np.testing.assert_array_equal(
        _column(out, "kinase|g__Blautia.s__Blautia_obeum"),
        _column(by_taxon, "2.7.1.1|g__Blautia.s__Blautia_obeum")
        + _column(by_taxon, "2.7.1.2|g__Blautia.s__Blautia_obeum"),
    )


def test_mean_divides_by_members_present():
    hierarchy = pd.DataFrame({"child": ["2.7.1.1", "2.7.1.2", "9.9.9.9"], "parent": "P", "level": "pathway"})
    out = bt.fn.func_glom(_function(), "pathway", hierarchy=hierarchy, agg="mean")
    np.testing.assert_allclose(
        _column(out, "P"), (_column(_function(), "2.7.1.1") + _column(_function(), "2.7.1.2")) / 2
    )


def test_keeps_obs_and_x_kind_and_drops_derived_slots():
    function = bt.pp.relative(_function())
    out = bt.fn.func_glom(function, "class", hierarchy=EC)
    assert out.obs["group"].tolist() == function.obs["group"].tolist()
    assert out.uns["biotapy"]["x_kind"] == "rpk" and "relative" not in out.layers
    assert '"step": "fn.func_glom"' in out.uns["biotapy"]["provenance"][-1]


def test_output_round_trips_through_h5ad(tmp_path):
    out = bt.fn.func_glom(_by_taxon(), "role", hierarchy=EC)
    out.write_h5ad(tmp_path / "glom.h5ad")
    back = ad.read_h5ad(tmp_path / "glom.h5ad")
    assert back.var["taxon"].tolist() == out.var["taxon"].tolist() and back.var["name"].isna().all()


def test_input_unchanged(assert_unchanged):
    function = _function()
    before = function.copy()
    bt.fn.func_glom(function, "class", hierarchy=EC)
    assert_unchanged(before, function)


def test_all_zero_sample_stays_zero():
    function = _function()
    dense = function.X.toarray()
    dense[0] = 0
    function.X = sp.csr_matrix(dense)
    assert bt.fn.func_glom(function, "class", hierarchy=EC).X[0].nnz == 0


def test_all_zero_feature_keeps_its_group():
    function = _function()
    dense = function.X.toarray()
    dense[:, function.var_names.get_loc("1.1.1.1")] = 0
    function.X = sp.csr_matrix(dense)
    out = bt.fn.func_glom(function, "class", hierarchy=EC)
    assert "1.-.-.-" in out.var_names and _column(out, "1.-.-.-").sum() == 0


def test_single_sample():
    out = bt.fn.func_glom(_function()[:1].copy(), "class", hierarchy=EC)
    assert out.shape == (1, 4)


def test_empty_modality_gives_no_groups():
    empty = _by_taxon()[:, []].copy()
    assert bt.fn.func_glom(empty, "class", hierarchy=EC).shape == (6, 0)


def test_no_feature_in_the_hierarchy_raises_with_examples():
    hierarchy = pd.DataFrame({"child": ["EC:2.7.1.1"], "parent": ["P"], "level": "pathway"})
    with pytest.raises(
        ValueError, match=r"features: \['1.1.1.1', '2.7.1.1', '2.7.1.2'\], children: \['EC:2.7.1.1'\]"
    ):
        bt.fn.func_glom(_function(), "pathway", hierarchy=hierarchy)


def test_error_examples_are_distinct_ids():
    hierarchy = pd.DataFrame({"child": ["x", "x", "y"], "parent": ["P", "Q", "P"], "level": "pathway"})
    with pytest.raises(ValueError, match=r"children: \['x', 'y'\]"):
        bt.fn.func_glom(_function(), "pathway", hierarchy=hierarchy)


def test_only_specials_are_summed_or_passed_through():
    # humann_regroup_table: UNGROUPED is summed into UNGROUPED, UNMAPPED maps to itself.
    function = _function()[:, ["UNMAPPED", "UNGROUPED"]].copy()
    out = bt.fn.func_glom(function, "class", hierarchy=EC)
    assert out.var_names.tolist() == ["UNGROUPED", "UNMAPPED"]
    np.testing.assert_array_equal(_column(out, "UNGROUPED"), _column(function, "UNGROUPED"))
    np.testing.assert_array_equal(_column(out, "UNMAPPED"), _column(function, "UNMAPPED"))


def test_missing_child_or_parent_raises_naming_hierarchy():
    for column in ("child", "parent"):
        broken = EC.copy()
        broken.loc[0, column] = np.nan
        with pytest.raises(ValueError, match="hierarchy has a missing value"):
            bt.fn.func_glom(_function(), "class", hierarchy=broken)


def test_generic_input_gets_no_taxon_columns():
    function = _function()
    function.var["genus"] = "Foo"
    out = bt.fn.func_glom(function, "class", hierarchy=EC)
    assert "genus" not in out.var.columns


def test_unknown_level_names_the_levels():
    with pytest.raises(KeyError, match=r"its levels: \['class', 'role'\]"):
        bt.fn.func_glom(_function(), "pathway", hierarchy=EC)


def test_missing_hierarchy_column_raises():
    with pytest.raises(KeyError, match="missing"):
        bt.fn.func_glom(_function(), "class", hierarchy=EC.drop(columns="level"))


def test_unknown_agg_raises():
    with pytest.raises(ValueError, match="agg="):
        bt.fn.func_glom(_function(), "class", hierarchy=EC, agg="median")


@st.composite
def _cases(draw):
    n_obs, n_var = draw(st.integers(1, 4)), draw(st.integers(1, 6))
    x = np.array(draw(st.lists(st.integers(0, 20), min_size=n_obs * n_var, max_size=n_obs * n_var))).reshape(
        n_obs, n_var
    )
    pairs = draw(st.sets(st.tuples(st.integers(0, n_var - 1), st.integers(0, 2)), min_size=1))
    return x, sorted(pairs)


@given(_cases())
def test_total_equals_each_feature_times_its_parent_count(case):
    # The roadmap's invariant: mapped abundance x membership count is preserved; unmapped goes to UNGROUPED.
    x, pairs = case
    adata = AnnData(
        X=sp.csr_matrix(x),
        obs=pd.DataFrame(index=[f"s{i}" for i in range(x.shape[0])]),
        var=pd.DataFrame(index=[f"f{j}" for j in range(x.shape[1])]),
    )
    hierarchy = pd.DataFrame(
        {"child": [f"f{j}" for j, _ in pairs], "parent": [f"P{g}" for _, g in pairs], "level": "l"}
    )
    out = bt.fn.func_glom(adata, "l", hierarchy=hierarchy)
    parents = np.bincount([j for j, _ in pairs], minlength=x.shape[1])
    weights = np.where(parents > 0, parents, 1)
    np.testing.assert_array_equal(np.asarray(out.X.sum(axis=1)).ravel(), x @ weights)
