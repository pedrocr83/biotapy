import numpy as np
import pandas as pd
import pytest

import biotapy as bt


def _toy_with():
    tdata = bt.datasets.toy()
    bt.tl.beta(tdata, inplace=True)
    return tdata


def test_permanova_pseudo_f_by_hand():
    tdata = _toy_with()
    d2 = tdata.obsp["braycurtis"] ** 2
    total = d2[np.triu_indices(6, 1)].sum() / 6
    within = sum(d2[np.ix_(idx, idx)][np.triu_indices(3, 1)].sum() / 3 for idx in ([0, 1, 2], [3, 4, 5]))
    expected = (total - within) / (within / 4)
    result = bt.tl.permanova(tdata, "group", seed=0)
    assert result["test statistic"] == pytest.approx(expected)
    assert result["number of groups"] == 2 and result["sample size"] == 6


def test_permanova_same_seed_same_p_value():
    first = bt.tl.permanova(_toy_with(), "group", permutations=99, seed=1)
    again = bt.tl.permanova(_toy_with(), "group", permutations=99, seed=np.random.default_rng(1))
    assert first["p-value"] == again["p-value"] and first["number of permutations"] == 99


def test_permanova_input_unchanged(assert_unchanged):
    tdata = _toy_with()
    before = tdata.copy()
    bt.tl.permanova(tdata, "group", seed=0)
    assert_unchanged(before, tdata)


def test_permanova_unknown_grouping_is_named():
    with pytest.raises(KeyError, match="grouping='site'"):
        bt.tl.permanova(_toy_with(), "site")


def test_permanova_numeric_grouping_raises():
    # adonis2 would fit a numeric column as one continuous term, not one group per value.
    tdata = _toy_with()
    tdata.obs["dose"] = [1, 1, 2, 2, 3, 3]
    with pytest.raises(TypeError, match=r"grouping='dose' is a numeric column.*\.astype\(\"category\"\)"):
        bt.tl.permanova(tdata, "dose", seed=0)


@pytest.mark.parametrize(
    ("values", "n_groups"), [(pd.Categorical([1, 1, 2, 2, 3, 3]), 3), ([True, True, True, False, False, False], 2)]
)
def test_permanova_category_or_bool_grouping_works(values, n_groups):
    tdata = _toy_with()
    tdata.obs["dose"] = values
    assert bt.tl.permanova(tdata, "dose", seed=0)["number of groups"] == n_groups


def test_permanova_missing_group_raises():
    tdata = _toy_with()
    tdata.obs["group"] = tdata.obs["group"].cat.add_categories("C")
    tdata.obs.loc["s1", "group"] = np.nan
    with pytest.raises(ValueError, match=r"missing for 1 sample"):
        bt.tl.permanova(tdata, "group")


def test_permanova_missing_distance_names_the_call():
    with pytest.raises(KeyError, match=r"bt.tl.beta\(adata, metric='jaccard', inplace=True\)"):
        bt.tl.permanova(bt.datasets.toy(), "group", distance="jaccard")
