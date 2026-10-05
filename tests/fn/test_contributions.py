import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt
from biotapy._core import make_function_mudata

EC = pd.DataFrame({"child": ["2.7.1.1", "2.7.1.2"], "parent": ["kinase", "kinase"], "level": "role"})


def _by_taxon():
    return bt.datasets.toy_humann()["function_by_taxon"]


def _strata(adata, function):
    """The function's stratified columns, samples x var_names, straight from X."""
    columns = adata.var_names[adata.var["function"] == function]
    return adata[:, columns].X.toarray()


def test_one_column_per_taxon_ordered_by_total():
    out = bt.fn.contributions(_by_taxon(), "2.7.1.2")
    # Totals over the six samples: Blautia obeum 25, Bacteroides ovatus 21.
    assert out.columns.tolist() == ["g__Blautia.s__Blautia_obeum", "g__Bacteroides.s__Bacteroides_ovatus"]
    assert out.columns.name == "taxon" and out.index.tolist() == ["s1", "s2", "s3", "s4", "s5", "s6"]
    np.testing.assert_array_equal(out.loc["s4"].to_numpy(), [7.0, 1.0])
    assert out.dtypes.eq(np.float64).all()


def test_unclassified_is_an_ordinary_taxon():
    out = bt.fn.contributions(_by_taxon(), "1.1.1.1")
    assert out.columns.tolist() == ["g__Bacteroides.s__Bacteroides_ovatus", "unclassified"]
    np.testing.assert_array_equal(out["unclassified"].to_numpy(), [3, 2, 4, 1, 0, 1])


def test_special_functions_are_queried_like_any_other():
    out = bt.fn.contributions(_by_taxon(), "UNGROUPED")
    assert out.columns.tolist() == ["unclassified"]


def test_top_keeps_the_largest_taxa_and_sums_the_rest_into_other():
    out = bt.fn.contributions(_by_taxon(), "2.7.1.2", top=1)
    assert out.columns.tolist() == ["g__Blautia.s__Blautia_obeum", "other"] and out.columns.name == "taxon"
    np.testing.assert_array_equal(out["other"].to_numpy(), [6, 5, 7, 1, 0, 2])


def test_top_at_or_above_the_taxon_count_adds_no_other():
    assert "other" not in bt.fn.contributions(_by_taxon(), "2.7.1.2", top=2).columns


def test_ties_are_broken_by_taxon_name():
    obs = pd.DataFrame(index=["s1", "s2"])
    ids = pd.Index(["K1", "K1|c", "K1|a", "K1|b"])
    mdata = make_function_mudata(np.array([[3, 1, 1, 1], [3, 1, 1, 1]]), obs=obs, row_ids=ids, x_kind="rpk", source="t")
    out = bt.fn.contributions(mdata["function_by_taxon"], "K1", top=2)
    assert out.columns.tolist() == ["a", "b", "other"]


def test_reads_regrouped_and_renormalised_tables():
    by_role = bt.fn.func_glom(_by_taxon(), "role", hierarchy=EC)
    out = bt.fn.contributions(by_role, "kinase")
    np.testing.assert_array_equal(out.sum(axis=1).to_numpy(), _strata(by_role, "kinase").sum(axis=1))
    renormed = bt.fn.renorm(bt.datasets.toy_humann(), "relab")["function_by_taxon"]
    shares = bt.fn.contributions(renormed, "2.7.1.2")
    np.testing.assert_allclose(shares.sum(axis=1).to_numpy(), _strata(renormed, "2.7.1.2").sum(axis=1))


def test_reads_a_picrust2_table(tmp_path):
    (tmp_path / "unstrat.tsv").write_text("function\tS1\tS2\nEC:1.1.1.1\t6\t4\n")
    header = "sample\tfunction\ttaxon\ttaxon_abun\ttaxon_rel_abun\tgenome_function_count"
    header += "\ttaxon_function_abun\ttaxon_rel_function_abun\tnorm_taxon_function_contrib\n"
    rows = "S1\tEC:1.1.1.1\tASV1\t3\t50\t2\t6\t100\t1\nS2\tEC:1.1.1.1\tRARE\t4\t100\t1\t4\t100\t1\n"
    (tmp_path / "contrib.tsv").write_text(header + rows)
    mdata = bt.io.read_picrust2(tmp_path / "unstrat.tsv", contrib=tmp_path / "contrib.tsv")
    out = bt.fn.contributions(mdata["function_by_taxon"], "1.1.1.1")
    assert out.columns.tolist() == ["ASV1", "RARE"] and out.loc["S2"].tolist() == [0.0, 4.0]


def test_input_unchanged(assert_unchanged):
    by_taxon = _by_taxon()
    before = by_taxon.copy()
    bt.fn.contributions(by_taxon, "2.7.1.2", top=1)
    assert_unchanged(before, by_taxon)


def test_all_zero_sample_is_a_row_of_zeros():
    by_taxon = _by_taxon()
    dense = by_taxon.X.toarray()
    dense[0] = 0
    by_taxon.X = sp.csr_matrix(dense)
    assert bt.fn.contributions(by_taxon, "2.7.1.2", top=1).loc["s1"].tolist() == [0.0, 0.0]


def test_all_zero_taxon_is_kept_last():
    by_taxon = _by_taxon()
    dense = by_taxon.X.toarray()
    dense[:, by_taxon.var_names.get_loc("2.7.1.2|g__Blautia.s__Blautia_obeum")] = 0
    by_taxon.X = sp.csr_matrix(dense)
    out = bt.fn.contributions(by_taxon, "2.7.1.2")
    assert out.columns[-1] == "g__Blautia.s__Blautia_obeum" and out.iloc[:, -1].eq(0).all()


def test_single_sample():
    out = bt.fn.contributions(_by_taxon()[:1].copy(), "2.7.1.2")
    assert out.shape == (1, 2)


def test_unknown_function_names_close_ids():
    match = r"function='2\.7\.1\.3' has no rows in adata; close function ids: \['2\.7\.1\.2', '2\.7\.1\.1'\]"
    with pytest.raises(KeyError, match=match):
        bt.fn.contributions(_by_taxon(), "2.7.1.3")


def test_community_modality_names_the_one_to_pass():
    with pytest.raises(KeyError, match=r"adata needs var columns \['function', 'taxon'\].*function_by_taxon"):
        bt.fn.contributions(bt.datasets.toy_humann()["function"], "2.7.1.2")


def test_mudata_is_refused_naming_the_modality():
    with pytest.raises(TypeError, match=r"mdata\['function_by_taxon'\]"):
        bt.fn.contributions(bt.datasets.toy_humann(), "2.7.1.2")


@pytest.mark.parametrize("top", [0, -1, 1.5, True, "3"])
def test_top_must_be_a_positive_integer(top):
    with pytest.raises(ValueError, match="top="):
        bt.fn.contributions(_by_taxon(), "2.7.1.2", top=top)


def test_a_taxon_named_other_cannot_be_hidden_by_top():
    obs = pd.DataFrame(index=["s1"])
    ids = pd.Index(["K1", "K1|other", "K1|a", "K1|b"])
    mdata = make_function_mudata(np.array([[3, 1, 1, 1]]), obs=obs, row_ids=ids, x_kind="rpk", source="t")
    with pytest.raises(ValueError, match="top=1.*'other'"):
        bt.fn.contributions(mdata["function_by_taxon"], "K1", top=1)
    assert "other" in bt.fn.contributions(mdata["function_by_taxon"], "K1").columns


@settings(deadline=None)
@given(
    arrays(np.float64, st.tuples(st.integers(1, 4), st.integers(1, 6)), elements=st.floats(0, 1e6)),
    st.none() | st.integers(1, 7),
)
def test_rows_sum_to_the_functions_strata(strata, top):
    n_obs, n_taxa = strata.shape
    ids = pd.Index(["K1", "K2", *(f"K1|t{j}" for j in range(n_taxa)), "K2|t0"])
    X = np.hstack([strata.sum(axis=1, keepdims=True), np.ones((n_obs, 1)), strata, np.ones((n_obs, 1))])
    obs = pd.DataFrame(index=[f"s{i}" for i in range(n_obs)])
    by_taxon = make_function_mudata(X, obs=obs, row_ids=ids, x_kind="rpk", source="t")["function_by_taxon"]
    out = bt.fn.contributions(by_taxon, "K1", top=top)
    np.testing.assert_allclose(out.sum(axis=1).to_numpy(), strata.sum(axis=1), rtol=1e-12)
    assert out.shape[1] == (n_taxa if top is None or top >= n_taxa else top + 1)
