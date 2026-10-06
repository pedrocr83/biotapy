import sys
import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.stats import false_discovery_control

import biotapy as bt

# What biotapy's R call returns for bt.da.maaslin3(toy, "group", seed=0): maaslin3's abundance results (sorted by its
# q-value) for the group term x0, 6 significant digits; test_maaslin3_on_toy_is_the_canned_table checks it in R.
CANNED = pd.DataFrame(
    {
        "feature": ["f6", "f7", "f8", "f1", "f3", "f2", "f4", "f5"],
        "metadata": "x0",
        "coef": [4.64115, 4.38614, 2.45111, -1.85211, -1.50564, -1.0101, -0.625606, 0.625606],
        "stderr": [0.452717, 0.372327, 0.667219, 0.378469, 0.320919, 0.507362, 0.187502, 0.994888],
        "pval": [0.0337446, 0.036925, 0.206379, 0.226817, 0.27898, 0.474965, 0.619162, 0.657504],
        "failed": False,
    }
)
EXPECTED = CANNED.set_index("feature").reindex(bt.datasets.toy().var_names)


def test_maaslin3_maps_the_group_rows_onto_the_schema(fake_rpy2):
    covariate = CANNED.assign(metadata="x1", coef=9.0, pval=1e-9)  # rows of a covariate are not the result
    fake_rpy2.output = pd.concat([covariate, CANNED], ignore_index=True)
    out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
    assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
    np.testing.assert_array_equal(out["effect"], EXPECTED["coef"])
    np.testing.assert_array_equal(out["se"], EXPECTED["stderr"])
    np.testing.assert_array_equal(out["pvalue"], EXPECTED["pval"])
    # BH over the group's p-values only; maaslin3's qval_individual would pool them with the covariate's.
    np.testing.assert_allclose(out["qvalue"], false_discovery_control(EXPECTED["pval"]), rtol=1e-12)
    assert (out["method"] == "maaslin3").all() and (out["contrast"] == "B vs A").all()


def test_maaslin3_gets_counts_plain_column_names_a_formula_and_one_seed(fake_rpy2):
    fake_rpy2.output = CANNED
    tdata = bt.datasets.toy()
    tdata.obs["body site"] = pd.Categorical(["x", "y", "x", "y", "x", "y"])
    bt.da.maaslin3(tdata, "group", covariates=["body site"], reference="B", seed=5)
    ((code, (counts, metadata, formula, seed)),) = fake_rpy2.calls
    assert 'evaluate_only = "abundance"' in code and "subtract_median = TRUE" in code and "set.seed(seed)" in code
    np.testing.assert_array_equal(counts.to_numpy(), tdata.X.toarray())
    assert counts.index.tolist() == tdata.obs_names.tolist() and counts.columns.tolist() == tdata.var_names.tolist()
    assert formula == "~ x0 + x1" and metadata.columns.tolist() == ["x0", "x1"]
    assert metadata["x0"].cat.categories.tolist() == ["B", "A"]  # maaslin3 keeps a factor's first level as reference
    assert seed == int(np.random.default_rng(5).integers(2**31 - 1))


def test_bool_levels_reach_r_as_strings_in_reference_order(fake_rpy2):
    fake_rpy2.output = CANNED
    tdata = bt.datasets.toy()
    tdata.obs["treated"] = [False] * 3 + [True] * 3
    out = bt.da.maaslin3(tdata, "treated", reference="True", seed=0)
    assert fake_rpy2.calls[0][1][1]["x0"].cat.categories.tolist() == ["True", "False"]
    assert (out["contrast"] == "False vs True").all()


def test_failed_fit_is_not_tested(fake_rpy2):
    failed = CANNED["feature"] == "f3"
    fake_rpy2.output = CANNED.assign(failed=failed)
    out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
    assert out.loc["f3", ["effect", "se", "pvalue", "qvalue"]].isna().all() and out.loc["f3", "direction"] == 0
    np.testing.assert_allclose(out["qvalue"].drop("f3"), false_discovery_control(EXPECTED["pval"].drop("f3")))


def test_feature_without_a_row_is_not_tested(fake_rpy2):
    fake_rpy2.output = CANNED[CANNED["feature"] != "f7"]
    out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
    assert out.loc["f7", ["effect", "pvalue"]].isna().all() and out["effect"].notna().sum() == 7


def test_numeric_group_reports_its_slope(fake_rpy2):
    fake_rpy2.output = CANNED
    tdata = bt.datasets.toy()
    tdata.obs["ph"] = [5.1, 5.6, 6.0, 6.8, 7.1, 7.4]
    out = bt.da.maaslin3(tdata, "ph", seed=0)
    assert (out["contrast"] == "ph").all()
    assert fake_rpy2.calls[0][1][1]["x0"].dtype == np.float64


def test_missing_rpy2_names_the_extra(monkeypatch):
    for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(ImportError, match=r"pip install 'biotapy\[r\]'"):
        bt.da.maaslin3(bt.datasets.toy(), "group")


def test_missing_r_package_names_the_install_line(fake_rpy2):
    fake_rpy2.installed = {"ALDEx2"}
    with pytest.raises(ImportError, match=r'BiocManager::install\("maaslin3"\)'):
        bt.da.maaslin3(bt.datasets.toy(), "group")


def test_shared_model_checks_apply(fake_rpy2):
    tdata = bt.datasets.toy()
    tdata.obs["arm"] = ["x"] * 3 + ["y"] * 3
    with pytest.raises(ValueError, match="are collinear"):
        bt.da.maaslin3(tdata, "group", covariates=["arm"])
    with pytest.raises(ValueError, match=r"obs\['age'\] is missing for 1 sample"):
        tdata.obs["age"] = [30.0, np.nan, 52, 38, 45, 60]
        bt.da.maaslin3(tdata, "group", covariates=["age"])
    assert fake_rpy2.calls == []


def test_maaslin3_keeps_input(fake_rpy2, assert_unchanged):
    fake_rpy2.output = CANNED
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.da.maaslin3(tdata, "group", seed=0)
    assert_unchanged(before, tdata)


def test_seed_accepts_a_generator(fake_rpy2):
    fake_rpy2.output = CANNED
    bt.da.maaslin3(bt.datasets.toy(), "group", seed=np.random.default_rng(3))
    bt.da.maaslin3(bt.datasets.toy(), "group", seed=3)
    assert fake_rpy2.calls[0][1][3] == fake_rpy2.calls[1][1][3]
    generator = np.random.default_rng(3)
    bt.da.maaslin3(bt.datasets.toy(), "group", seed=generator)
    bt.da.maaslin3(bt.datasets.toy(), "group", seed=generator)
    assert fake_rpy2.calls[2][1][3] != fake_rpy2.calls[3][1][3]  # the caller's generator advances


@pytest.mark.parametrize(("axis", "fix"), [("var", "var_names_make_unique"), ("obs", "obs_names_make_unique")])
def test_repeated_names_raise(fake_rpy2, axis, fix):
    tdata = bt.datasets.toy()
    names = getattr(tdata, f"{axis}_names").tolist()
    setattr(tdata, f"{axis}_names", ["x", "x", *names[2:]])
    with pytest.raises(ValueError, match=rf"unique {axis} names.*\['x'\].*adata\.{fix}\(\)"):
        bt.da.maaslin3(tdata, "group")
    assert fake_rpy2.calls == []


def test_seed_is_checked_before_rpy2_is_imported(monkeypatch):
    for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(TypeError, match="seed must be"):
        bt.da.maaslin3(bt.datasets.toy(), "group", seed="abc")


def test_each_r_warning_is_re_emitted_at_the_callers_line(fake_rpy2):
    fake_rpy2.output = CANNED
    fake_rpy2.warnings = ["values are unreliable"]
    with pytest.warns(UserWarning) as record:
        bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
    assert [str(w.message) for w in record] == ["da.maaslin3: R warned: values are unreliable"]
    assert record[0].filename == __file__


def test_r_errors_name_the_function(fake_rpy2):
    fake_rpy2.error = "boom\n"
    with pytest.raises(RuntimeError, match=r"^da\.maaslin3: R stopped: boom$"):
        bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)


def test_r_warnings_before_an_error_are_in_its_message(fake_rpy2):
    fake_rpy2.warnings = ["first", "second"]
    fake_rpy2.error = "boom"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(RuntimeError, match=r"^da\.maaslin3: R stopped: boom \(R warned first: first; second\)$"):
            bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)


def test_ordered_categoricals_reach_r_unordered(fake_rpy2):
    # R codes an ordered factor by polynomial contrasts, which would shrink a two-level effect by 1/sqrt(2).
    fake_rpy2.output = CANNED
    tdata = bt.datasets.toy()
    tdata.obs["group"] = pd.Categorical(tdata.obs["group"], categories=["A", "B"], ordered=True)
    tdata.obs["site"] = pd.Categorical(["x", "y", "z", "x", "y", "z"], categories=["x", "y", "z"], ordered=True)
    bt.da.maaslin3(tdata, "group", covariates=["site"], seed=0)
    metadata = fake_rpy2.calls[0][1][1]
    assert not metadata["x0"].cat.ordered and not metadata["x1"].cat.ordered


@pytest.mark.r
def test_maaslin3_on_toy_is_the_canned_table():
    out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
    # CANNED holds 6 significant digits of R's output, hence rtol=1e-5.
    for ours, theirs in {"effect": "coef", "se": "stderr", "pvalue": "pval"}.items():
        np.testing.assert_allclose(out[ours], EXPECTED[theirs], rtol=1e-5, err_msg=ours)


@pytest.mark.r
def test_seed_reproduces_the_table_in_r():
    tdata = bt.datasets.toy()
    first = bt.da.maaslin3(tdata, "group", seed=0)
    pd.testing.assert_frame_equal(bt.da.maaslin3(tdata, "group", seed=0), first)
    assert not bt.da.maaslin3(tdata, "group", seed=1)["pvalue"].equals(first["pvalue"])


@pytest.mark.r
def test_reference_sets_the_sign_in_r():
    tdata = bt.datasets.toy()
    a_first, b_first = bt.da.maaslin3(tdata, "group", seed=0), bt.da.maaslin3(tdata, "group", reference="B", seed=0)
    assert (b_first["contrast"] == "A vs B").all()
    np.testing.assert_allclose(b_first["effect"], -a_first["effect"], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(b_first["se"], a_first["se"], rtol=1e-9)


@pytest.mark.r
def test_covariates_and_odd_column_names_reach_r():
    tdata = bt.datasets.toy()
    tdata.obs["body site"] = pd.Categorical(["x", "y", "x", "y", "x", "y"])
    tdata.obs["treated"] = [False] * 3 + [True] * 3
    adjusted = bt.da.maaslin3(tdata, "treated", covariates=["body site"], seed=0)
    assert (adjusted["contrast"] == "True vs False").all() and adjusted["effect"].notna().all()
    assert not np.allclose(adjusted["se"], bt.da.maaslin3(tdata, "treated", seed=0)["se"])


@pytest.mark.r
def test_bool_reference_sets_the_sign_in_r():
    tdata = bt.datasets.toy()
    tdata.obs["treated"] = [False] * 3 + [True] * 3
    default, swapped = (
        bt.da.maaslin3(tdata, "treated", seed=0),
        bt.da.maaslin3(tdata, "treated", reference="True", seed=0),
    )
    np.testing.assert_allclose(swapped["effect"], -default["effect"], rtol=1e-9, atol=1e-12)


@pytest.mark.r
def test_all_zero_feature_is_not_tested_in_r():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[:, tdata.var_names.get_loc("f7")] = 0
    tdata.X = sp.csr_matrix(dense)
    out = bt.da.maaslin3(tdata, "group", seed=0)
    assert out.loc["f7", ["effect", "pvalue", "qvalue"]].isna().all()
    assert np.isfinite(out.drop(index="f7")[["effect", "pvalue", "qvalue"]].to_numpy()).all()


@pytest.mark.r
def test_maaslin3_table_passes_the_consensus_checks():
    out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
    assert bt.da.consensus([out], min_methods=1)["n_tested"].tolist() == [1] * 8


@pytest.mark.r
def test_ordered_group_gives_the_unordered_table_in_r():
    tdata = bt.datasets.toy()
    plain = bt.da.maaslin3(tdata, "group", seed=0)
    tdata.obs["group"] = pd.Categorical(tdata.obs["group"], categories=["A", "B"], ordered=True)
    pd.testing.assert_frame_equal(bt.da.maaslin3(tdata, "group", seed=0), plain)


def test_all_zero_sample_raises_before_r(fake_rpy2):
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[1] = 0
    tdata.X = sp.csr_matrix(dense)
    with pytest.raises(ValueError, match=r"1 sample\(s\) have none"):
        bt.da.maaslin3(tdata, "group")
    assert fake_rpy2.calls == []


def test_single_sample_raises(fake_rpy2):
    with pytest.raises(ValueError, match=r"compares two levels, but obs\['group'\] has 1"):
        bt.da.maaslin3(bt.datasets.toy()[:1].copy(), "group")
    assert fake_rpy2.calls == []


@pytest.mark.r
@pytest.mark.parametrize("seeded", [True, False])
def test_r_random_state_is_left_as_found(seeded):
    from rpy2.robjects import r

    drop = 'if (exists(".Random.seed", globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())'
    r("set.seed(1)" if seeded else drop)
    before = list(r(".Random.seed")) if seeded else None
    bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
    assert r('exists(".Random.seed", globalenv(), inherits = FALSE)')[0] == seeded
    if seeded:
        assert list(r(".Random.seed")) == before
