import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba

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


def test_dots_sit_at_their_row_and_method_with_the_colour_of_their_sign(ax):
    bt.pl.consensus(_consensus_table(), ax=ax)
    by_label = {collection.get_label(): collection for collection in ax.collections}
    # Rows are x, z, y (top first); columns are a, b. x is called up by both, y down by a only.
    up, down = by_label["effect > 0"], by_label["effect < 0"]
    assert sorted(map(tuple, up.get_offsets().tolist())) == [(0, 0), (1, 0), (1, 1)]
    assert down.get_offsets().tolist() == [[0, 2]]
    assert up.get_facecolor().tolist() == [list(to_rgba("#d62728"))]
    assert down.get_facecolor().tolist() == [list(to_rgba("#1f77b4"))]


def test_a_non_table_raises(ax):
    with pytest.raises(TypeError, match="table must be a pandas DataFrame, got list"):
        bt.pl.consensus([1, 2], ax=ax)


@pytest.mark.parametrize("top", [2.5, True, "3"])
def test_top_that_is_not_an_integer_raises(ax, top):
    with pytest.raises(TypeError, match="top must be an integer"):
        bt.pl.consensus(_consensus_table(), top=top, ax=ax)


def test_numpy_integer_top_is_accepted(ax):
    bt.pl.consensus(_consensus_table(), top=np.int64(1), ax=ax)
    assert [label.get_text() for label in ax.get_yticklabels()] == ["x"]


def test_an_empty_table_raises(ax):
    with pytest.raises(ValueError, match="no method calls any feature significant: nothing to draw"):
        bt.pl.consensus(_consensus_table().iloc[0:0], ax=ax)


def test_the_legend_is_outside_the_axes(ax):
    bt.pl.consensus(_consensus_table(), ax=ax)
    assert ax.get_legend().get_bbox_to_anchor().transformed(ax.transAxes.inverted()).x0 > 1
