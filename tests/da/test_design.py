import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import biotapy as bt

METHODS = [bt.da.ancombc2, bt.da.linda]
# ANCOM-BC2's bias E-M is fitted against the reference level, so swapping it flips the effects only approximately:
# 9e-4 log2 on toy(), which 2e-3 covers, but 0.08-0.40 log2 on GlobalPatterns genera (in R too), where this would fail.
SWAP_ATOL = {"ancombc2": 2e-3, "linda": 1e-12}


def _toy_with(**columns):
    tdata = bt.datasets.toy()
    for name, values in columns.items():
        tdata.obs[name] = values
    return tdata


@pytest.mark.parametrize("method", METHODS)
def test_reference_sets_the_sign(method):
    tdata = bt.datasets.toy()
    a_first, b_first = method(tdata, "group"), method(tdata, "group", reference="B")
    assert (a_first["contrast"] == "B vs A").all() and (b_first["contrast"] == "A vs B").all()
    atol = SWAP_ATOL[method.__name__]
    np.testing.assert_allclose(b_first["effect"], -a_first["effect"], rtol=1e-9, atol=atol)
    clear = a_first["effect"].abs() > atol
    assert (b_first.loc[clear, "direction"] == -a_first.loc[clear, "direction"]).all() and clear.sum() >= 7
    np.testing.assert_allclose(b_first["pvalue"], a_first["pvalue"], rtol=1e-9, atol=atol)


@pytest.mark.parametrize("method", METHODS)
def test_string_group_takes_the_first_sorted_level_as_reference(method):
    tdata = _toy_with(site=["gut", "gut", "gut", "air", "air", "air"])
    assert (method(tdata, "site")["contrast"] == "gut vs air").all()


@pytest.mark.parametrize("method", METHODS)
def test_bool_group_compares_true_with_false(method):
    tdata = _toy_with(treated=[False] * 3 + [True] * 3)
    np.testing.assert_allclose(method(tdata, "treated")["effect"], method(tdata, "group")["effect"], rtol=1e-12)
    assert (method(tdata, "treated")["contrast"] == "True vs False").all()


@pytest.mark.parametrize("method", METHODS)
def test_missing_group_value_raises(method):
    tdata = _toy_with(group=pd.Categorical(["A", None, "A", "B", "B", "B"]))
    with pytest.raises(ValueError, match=r"obs\['group'\] is missing for 1 sample\(s\); drop them first"):
        method(tdata, "group")


@pytest.mark.parametrize("method", METHODS)
def test_missing_covariate_value_raises(method):
    tdata = _toy_with(age=[30.0, np.nan, 52, 38, 45, 60])
    with pytest.raises(ValueError, match=r"obs\['age'\] is missing for 1 sample"):
        method(tdata, "group", covariates=["age"])


@pytest.mark.parametrize("method", METHODS)
def test_three_levels_raise(method):
    tdata = _toy_with(site=["a", "b", "c", "a", "b", "c"])
    with pytest.raises(ValueError, match=r"compares two levels, but obs\['site'\] has 3: \['a', 'b', 'c'\]"):
        method(tdata, "site")


@pytest.mark.parametrize("method", METHODS)
def test_unused_categories_are_not_levels(method):
    tdata = _toy_with(group=pd.Categorical(["A"] * 3 + ["B"] * 3, categories=["A", "B", "C"]))
    assert (method(tdata, "group")["contrast"] == "B vs A").all()


@pytest.mark.parametrize("method", METHODS)
def test_unknown_reference_raises(method):
    with pytest.raises(ValueError, match=r"reference='C' is not a level of obs\['group'\], which has \['A', 'B'\]"):
        method(bt.datasets.toy(), "group", reference="C")


@pytest.mark.parametrize("method", METHODS)
def test_reference_for_a_numeric_group_raises(method):
    tdata = _toy_with(ph=[5.1, 5.6, 6.0, 6.8, 7.1, 7.4])
    with pytest.raises(ValueError, match=r"reference='A' needs a categorical group, but obs\['ph'\] is numeric"):
        method(tdata, "ph", reference="A")


@pytest.mark.parametrize("method", METHODS)
def test_absent_column_raises(method):
    with pytest.raises(KeyError, match=r"\['batch'\] not in obs"):
        method(bt.datasets.toy(), "group", covariates=["batch"])


@pytest.mark.parametrize("method", METHODS)
def test_covariates_as_a_string_raises(method):
    with pytest.raises(TypeError, match=r"covariates must be a list of obs columns, such as \['group'\]"):
        method(bt.datasets.toy(), "group", covariates="group")


@pytest.mark.parametrize("method", METHODS)
def test_group_among_covariates_raises(method):
    with pytest.raises(ValueError, match="repeat a column"):
        method(bt.datasets.toy(), "group", covariates=["group"])


@pytest.mark.parametrize("method", METHODS)
def test_collinear_covariate_raises(method):
    tdata = _toy_with(arm=["x"] * 3 + ["y"] * 3)
    with pytest.raises(ValueError, match="are collinear; drop the covariate that repeats another"):
        method(tdata, "group", covariates=["arm"])


@pytest.mark.parametrize("method", METHODS)
def test_single_sample_raises(method):
    tdata = _toy_with(ph=[5.1, 5.6, 6.0, 6.8, 7.1, 7.4])[:1].copy()
    with pytest.raises(ValueError, match=r"needs more samples \(1\) than model terms \(2, with the intercept\)"):
        method(tdata, "ph")


@pytest.mark.parametrize("method", METHODS)
def test_all_zero_sample_raises(method):
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[1] = 0
    tdata.X = sp.csr_matrix(dense)
    with pytest.raises(ValueError, match=r"1 sample\(s\) have none \(\['s2'\]\).*bt.pp.filter_samples"):
        method(tdata, "group")


@pytest.mark.parametrize("method", METHODS)
def test_relative_abundances_raise(method):
    relative = bt.pp.relative(bt.datasets.toy())
    relative.X = relative.layers["relative"]
    relative.uns["biotapy"]["x_kind"] = "relative"
    with pytest.raises(ValueError, match="needs raw counts in X"):
        method(relative, "group")


@pytest.mark.parametrize("method", METHODS)
def test_single_feature_raises(method):
    with pytest.raises(ValueError, match=r"needs at least two features"):
        method(bt.datasets.toy()[:, :1].copy(), "group")


@pytest.mark.parametrize("method", METHODS)
def test_constant_numeric_group_raises(method):
    tdata = _toy_with(ph=[5.0] * 6)
    with pytest.raises(ValueError, match=r"group column 'ph' is constant across samples"):
        method(tdata, "ph")


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("values", [[3.0] * 6, ["x"] * 6])
def test_constant_covariate_raises(method, values):
    tdata = _toy_with(c=values)
    with pytest.raises(ValueError, match=r"covariates column 'c' is constant across samples; drop it"):
        method(tdata, "group", covariates=["c"])


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("covariates", [None, [1], ("age", 2)])
def test_covariates_must_be_a_sequence_of_names(method, covariates):
    with pytest.raises(TypeError, match=r"covariates must be a list of obs columns"):
        method(bt.datasets.toy(), "group", covariates=covariates)


@pytest.mark.parametrize("method", METHODS)
def test_reference_must_be_a_string(method):
    tdata = _toy_with(treated=[False] * 3 + [True] * 3)
    with pytest.raises(TypeError, match=r"reference must be the level's name as a string, such as 'False'"):
        method(tdata, "treated", reference=True)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize(("axis", "fix"), [("var", "var_names_make_unique"), ("obs", "obs_names_make_unique")])
def test_repeated_names_raise(method, axis, fix):
    tdata = bt.datasets.toy()
    names = getattr(tdata, f"{axis}_names").tolist()
    setattr(tdata, f"{axis}_names", ["x", "x", *names[2:]])
    with pytest.raises(ValueError, match=rf"unique {axis} names.*\['x'\].*adata\.{fix}\(\)"):
        method(tdata, "group")
