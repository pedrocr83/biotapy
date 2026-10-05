import numpy as np
import pandas as pd
import pytest

import biotapy as bt


def _consensus_table():
    """da.consensus's layout for three features and two methods, hand-built."""
    return pd.DataFrame(
        {
            "effect_a": [2.0, -1.0, 0.1, 0.5],
            "qvalue_a": [0.01, 0.02, 0.5, np.nan],
            "significant_a": [True, True, False, False],
            "effect_b": [1.0, -3.0, 4.0, 0.2],
            "qvalue_b": [0.01, 0.30, 0.01, 0.9],
            "significant_b": [True, False, True, False],
            "n_tested": [2, 2, 2, 1],
            "n_significant": [2, 1, 1, 0],
            "direction": np.array([1, -1, 1, 0], dtype=np.int8),
            "consensus": [True, False, False, False],
            "conflict": [False] * 4,
        },
        index=pd.Index(["x", "y", "z", "w"], name="feature"),
    )


def _dots(ax):
    return {collection.get_label(): len(collection.get_offsets()) for collection in ax.collections}


def test_consensus_draws_one_dot_per_call_and_per_tested_feature(ax):
    assert bt.pl.consensus(_consensus_table(), ax=ax) is ax
    # x is called by both, y by a only (b tested it), z by b only (a tested it); w is called by none and not drawn.
    assert _dots(ax) == {"effect > 0": 3, "effect < 0": 1, "not significant": 2}
    assert [label.get_text() for label in ax.get_xticklabels()] == ["a", "b"]


def test_rows_are_sorted_by_calls_then_mean_absolute_effect(ax):
    bt.pl.consensus(_consensus_table(), ax=ax)
    # x has two calls; z (mean |effect| 2.05) comes before y (2.0).
    assert [label.get_text() for label in ax.get_yticklabels()] == ["x", "z", "y"]
    assert [label.get_fontweight() for label in ax.get_yticklabels()] == ["bold", "normal", "normal"]
    assert ax.get_ylim() == (2.5, -0.5)


def test_top_keeps_the_first_rows(ax):
    bt.pl.consensus(_consensus_table(), top=1, ax=ax)
    assert [label.get_text() for label in ax.get_yticklabels()] == ["x"]
    assert _dots(ax) == {"effect > 0": 2}


def test_nothing_called_raises(ax):
    table = _consensus_table()
    table["n_significant"] = 0
    with pytest.raises(ValueError, match="no method calls any feature significant: nothing to draw"):
        bt.pl.consensus(table, ax=ax)


def test_top_below_one_raises(ax):
    with pytest.raises(ValueError, match="top must be at least 1, got 0"):
        bt.pl.consensus(_consensus_table(), top=0, ax=ax)


def test_a_result_table_instead_of_a_consensus_table_raises(ax):
    with pytest.raises(KeyError, match="pl.consensus draws the table bt.da.consensus returns"):
        bt.pl.consensus(bt.da.linda(bt.datasets.toy(), "group"), ax=ax)


def test_consensus_keeps_the_table(ax):
    table = _consensus_table()
    before = table.copy()
    bt.pl.consensus(table, ax=ax)
    pd.testing.assert_frame_equal(table, before)
