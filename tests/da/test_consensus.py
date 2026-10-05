import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

import biotapy as bt


def _table(method, effect, qvalue, *, features=None, contrast="B vs A"):
    """A result table as a bt.da method writes it; NaN qvalue marks an untested feature."""
    effect, qvalue = np.asarray(effect, dtype=float), np.asarray(qvalue, dtype=float)
    features = features or [f"f{i}" for i in range(len(effect))]
    effect = np.where(np.isnan(qvalue), np.nan, effect)
    return pd.DataFrame(
        {
            "effect": effect,
            "se": np.where(np.isnan(qvalue), np.nan, 0.5),
            "pvalue": qvalue / 2,
            "qvalue": qvalue,
            "direction": np.sign(np.nan_to_num(effect)).astype(np.int8),
            "method": method,
            "contrast": contrast,
        },
        index=pd.Index(features, name="feature"),
    )


def test_consensus_counts_the_methods_that_call_each_feature():
    a = _table("a", [2.0, -1.0, 0.5, 3.0], [0.01, 0.01, 0.30, 0.02])
    b = _table("b", [1.5, -2.0, 0.4, 2.0], [0.02, 0.20, 0.01, 0.03])
    out = bt.da.consensus([a, b])
    assert out["n_tested"].tolist() == [2, 2, 2, 2]
    assert out["n_significant"].tolist() == [2, 1, 1, 2]
    assert out["direction"].tolist() == [1, -1, 1, 1]
    assert out["consensus"].tolist() == [True, False, False, True]
    assert out["significant_a"].tolist() == [True, True, False, True]


def test_columns_follow_the_results_order():
    out = bt.da.consensus([_table("b", [1.0], [0.01]), _table("a", [1.0], [0.01])])
    assert out.columns.tolist() == [
        *["effect_b", "qvalue_b", "significant_b", "effect_a", "qvalue_a", "significant_a"],
        *["n_tested", "n_significant", "direction", "consensus", "conflict"],
    ]
    assert out.index.name == "feature" and out["direction"].dtype == np.int8


def test_q_equal_to_alpha_is_not_called():
    out = bt.da.consensus([_table("a", [1.0], [0.05]), _table("b", [1.0], [0.0499])])
    assert out["significant_a"].tolist() == [False] and out["n_significant"].tolist() == [1]


def test_min_methods_sets_how_many_must_agree():
    results = [_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01]), _table("c", [1.0], [0.4])]
    assert bt.da.consensus(results, min_methods=2)["consensus"].tolist() == [True]
    assert not bt.da.consensus(results, min_methods=3)["consensus"].any()


def test_conflict_is_flagged_and_never_consensus():
    out = bt.da.consensus([_table("a", [2.0], [0.01]), _table("b", [-2.0], [0.01]), _table("c", [2.0], [0.01])])
    assert out["conflict"].tolist() == [True] and out["direction"].tolist() == [0]
    assert out["consensus"].tolist() == [False]


def test_untested_is_not_counted_as_not_significant():
    out = bt.da.consensus([_table("a", [1.0, 1.0], [0.01, np.nan]), _table("b", [1.0, 1.0], [0.01, 0.01])])
    assert out["n_tested"].tolist() == [2, 1] and out["n_significant"].tolist() == [2, 1]
    assert not out["significant_a"].iloc[1] and np.isnan(out["qvalue_a"].iloc[1])


def test_features_are_the_union_of_the_tables():
    a = _table("a", [1.0, 1.0], [0.01, 0.01], features=["f1", "f2"])
    b = _table("b", [1.0, 1.0], [0.01, 0.01], features=["f2", "f3"])
    out = bt.da.consensus([a, b])
    assert out.index.tolist() == ["f1", "f2", "f3"]
    assert out["n_tested"].tolist() == [1, 2, 1] and out["consensus"].tolist() == [False, True, False]


def test_results_with_different_contrasts_raise():
    with pytest.raises(ValueError, match=r"different contrasts \['A vs B', 'B vs A'\]"):
        bt.da.consensus([_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01], contrast="A vs B")])


def test_repeated_method_raises():
    with pytest.raises(ValueError, match=r"repeat the method\(s\) \['a'\]"):
        bt.da.consensus([_table("a", [1.0], [0.01]), _table("a", [1.0], [0.02])])


@pytest.mark.parametrize("min_methods", [0, 3])
def test_min_methods_out_of_range_raises(min_methods):
    with pytest.raises(ValueError, match=r"min_methods must be between 1 and the number of results \(2\)"):
        bt.da.consensus([_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01])], min_methods=min_methods)


@pytest.mark.parametrize("alpha", [0, 1, 1.5])
def test_alpha_out_of_range_raises(alpha):
    with pytest.raises(ValueError, match="alpha must be between 0 and 1"):
        bt.da.consensus([_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01])], alpha=alpha)


def test_one_table_instead_of_a_list_raises():
    with pytest.raises(TypeError, match=r"a list of da result tables, such as \[table_a, table_b\]"):
        bt.da.consensus(_table("a", [1.0], [0.01]))


def test_consensus_of_the_native_methods():
    tdata = bt.datasets.toy()
    out = bt.da.consensus([bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")])
    assert out.index[out["consensus"]].tolist() == ["f6", "f7"] and (out["direction"][out["consensus"]] == 1).all()


def test_consensus_keeps_its_inputs():
    results = [_table("a", [1.0, -2.0], [0.01, 0.2]), _table("b", [1.0, -1.0], [0.01, np.nan])]
    before = [table.copy() for table in results]
    bt.da.consensus(results)
    for table, copy in zip(results, before, strict=True):
        pd.testing.assert_frame_equal(table, copy)


@settings(deadline=None)
@given(
    st.lists(
        st.lists(
            st.tuples(st.floats(-3, 3), st.sampled_from([0.001, 0.04, 0.05, 0.3, np.nan])), min_size=4, max_size=4
        ),
        min_size=1,
        max_size=4,
    ),
    st.integers(1, 4),
)
def test_consensus_needs_enough_calls_and_no_conflict(tables, min_methods):
    results = [_table(f"m{i}", [effect for effect, _ in rows], [q for _, q in rows]) for i, rows in enumerate(tables)]
    out = bt.da.consensus(results, min_methods=min(min_methods, len(results)))
    assert (out["n_significant"] <= out["n_tested"]).all()
    hits = out[out["consensus"]]
    assert (hits["n_significant"] >= min(min_methods, len(results))).all() and not hits["conflict"].any()
    assert (hits["direction"] != 0).all() and (out.loc[out["conflict"], "direction"] == 0).all()
