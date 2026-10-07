import sys
import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy.stats import false_discovery_control

import biotapy as bt

# ALDEx2::aldex(t(toy counts), rep(c("0", "1"), each = 3), mc.samples = 128, ...) after set.seed(1826701614), the
# integer bt.da.aldex2(toy, "group", seed=0) passes to R (6 significant digits; test_aldex2_on_toy_is_the_canned_table checks it in R).
CANNED = pd.DataFrame(
    {
        "diff.btw": [-3.70251, -2.8971, -2.66143, -1.72769, -0.198703, 4.06776, 3.66243, 2.69778],
        "effect": [-1.69089, -1.24091, -2.09973, -1.91244, -0.0639279, 2.14194, 1.99029, 1.03629],
        "we.ep": [0.08828, 0.125246, 0.0447655, 0.0520091, 0.736825, 0.074296, 0.0783649, 0.172146],
        "we.eBH": [0.219275, 0.269506, 0.154564, 0.161972, 0.880257, 0.28763, 0.296191, 0.46741],
    },
    index=[f"f{i}" for i in range(1, 9)],
)


def test_aldex2_maps_aldex_onto_the_schema(fake_rpy2):
    fake_rpy2.output = CANNED.iloc[::-1]  # R's row order does not matter
    out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
    assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
    np.testing.assert_array_equal(out["effect"], CANNED["diff.btw"])
    np.testing.assert_array_equal(out["pvalue"], CANNED["we.ep"])
    assert out["se"].isna().all()
    # BH of the expected p-values, as for every method; ALDEx2's we.eBH averages per-draw corrections instead.
    np.testing.assert_allclose(out["qvalue"], false_discovery_control(CANNED["we.ep"]), rtol=1e-12)
    assert (out["method"] == "aldex2").all() and (out["contrast"] == "B vs A").all()


def test_aldex2_sends_features_as_rows_the_reference_as_0_and_one_seed(fake_rpy2):
    fake_rpy2.output = CANNED
    bt.da.aldex2(bt.datasets.toy(), "group", mc_samples=64, reference="B", seed=7)
    ((code, (reads, conditions, mc_samples, seed)),) = fake_rpy2.calls
    assert "set.seed(seed)" in code and "ALDEx2::aldex(" in code
    assert reads.shape == (8, 6) and reads.index.tolist() == bt.datasets.toy().var_names.tolist()
    np.testing.assert_array_equal(reads.to_numpy(), bt.datasets.toy().X.toarray().T)
    assert conditions.tolist() == ["1"] * 3 + ["0"] * 3  # reference B is "0", which ALDEx2 sorts first
    assert mc_samples == 64 and seed == int(np.random.default_rng(7).integers(2**31 - 1))


def test_seed_accepts_a_generator(fake_rpy2):
    fake_rpy2.output = CANNED
    bt.da.aldex2(bt.datasets.toy(), "group", seed=np.random.default_rng(3))
    bt.da.aldex2(bt.datasets.toy(), "group", seed=3)
    assert fake_rpy2.calls[0][1][3] == fake_rpy2.calls[1][1][3]
    generator = np.random.default_rng(3)
    bt.da.aldex2(bt.datasets.toy(), "group", seed=generator)
    bt.da.aldex2(bt.datasets.toy(), "group", seed=generator)
    assert fake_rpy2.calls[2][1][3] != fake_rpy2.calls[3][1][3]  # the caller's generator advances


def test_feature_aldex2_drops_is_not_tested(fake_rpy2):
    fake_rpy2.output = CANNED.drop(index="f7")  # ALDEx2 removes a feature with no read before its draws
    out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
    assert out.loc["f7", ["effect", "pvalue", "qvalue"]].isna().all() and out.loc["f7", "direction"] == 0
    np.testing.assert_allclose(out["qvalue"].drop("f7"), false_discovery_control(CANNED["we.ep"].drop("f7")))


def test_missing_rpy2_names_the_extra(monkeypatch):
    for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(ImportError, match=r"pip install 'biotapy\[r\]'"):
        bt.da.aldex2(bt.datasets.toy(), "group")


def test_missing_r_package_names_the_install_line(fake_rpy2):
    fake_rpy2.installed = set()
    with pytest.raises(
        ImportError, match=r'da\.aldex2 needs the R package ALDEx2\. .*BiocManager::install\("ALDEx2"\)'
    ):
        bt.da.aldex2(bt.datasets.toy(), "group")


@pytest.mark.parametrize(("axis", "fix"), [("var", "var_names_make_unique"), ("obs", "obs_names_make_unique")])
def test_repeated_names_raise(fake_rpy2, axis, fix):
    tdata = bt.datasets.toy()
    names = getattr(tdata, f"{axis}_names").tolist()
    setattr(tdata, f"{axis}_names", ["x", "x", *names[2:]])
    with pytest.raises(ValueError, match=rf"unique {axis} names.*\['x'\].*adata\.{fix}\(\)"):
        bt.da.aldex2(tdata, "group")
    assert fake_rpy2.calls == []


def test_seed_is_checked_before_rpy2_is_imported(monkeypatch):
    for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
        monkeypatch.setitem(sys.modules, name, None)
    with pytest.raises(TypeError, match="seed must be"):
        bt.da.aldex2(bt.datasets.toy(), "group", seed="abc")


def test_each_r_warning_is_re_emitted_at_the_callers_line(fake_rpy2):
    fake_rpy2.output = CANNED
    fake_rpy2.warnings = ["values are unreliable", "second warning"]
    with pytest.warns(UserWarning) as record:
        bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
    assert [str(w.message) for w in record] == [
        "da.aldex2: R warned: values are unreliable",
        "da.aldex2: R warned: second warning",
    ]
    assert record[0].filename == __file__


def test_r_errors_name_the_function(fake_rpy2):
    fake_rpy2.error = "boom\n"
    with pytest.raises(RuntimeError, match=r"^da\.aldex2: R stopped: boom$"):
        bt.da.aldex2(bt.datasets.toy(), "group", seed=0)


def test_r_warnings_before_an_error_are_in_its_message(fake_rpy2):
    fake_rpy2.warnings = ["first", "second"]
    fake_rpy2.error = "boom"
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(RuntimeError, match=r"^da\.aldex2: R stopped: boom \(R warned first: first; second\)$"):
            bt.da.aldex2(bt.datasets.toy(), "group", seed=0)


def test_numeric_group_raises(fake_rpy2):
    tdata = bt.datasets.toy()
    tdata.obs["ph"] = [5.1, 5.6, 6.0, 6.8, 7.1, 7.4]
    with pytest.raises(ValueError, match=r"compares two levels, but obs\['ph'\] is numeric"):
        bt.da.aldex2(tdata, "ph")


def test_level_in_one_sample_raises(fake_rpy2):
    tdata = bt.datasets.toy()
    tdata.obs["arm"] = ["x"] * 5 + ["y"]
    with pytest.raises(ValueError, match=r"needs two samples in each level of obs\['arm'\]"):
        bt.da.aldex2(tdata, "arm")


@pytest.mark.parametrize(("mc_samples", "error"), [(1.5, TypeError), (True, TypeError), (0, ValueError)])
def test_mc_samples_must_be_a_positive_int(fake_rpy2, mc_samples, error):
    with pytest.raises(error, match="mc_samples"):
        bt.da.aldex2(bt.datasets.toy(), "group", mc_samples=mc_samples)


def test_shared_model_checks_apply(fake_rpy2):
    tdata = bt.datasets.toy()
    tdata.obs["site"] = ["a", "b", "c", "a", "b", "c"]
    with pytest.raises(ValueError, match=r"compares two levels, but obs\['site'\] has 3"):
        bt.da.aldex2(tdata, "site")
    dense = tdata.X.toarray()
    dense[1] = 0
    tdata.X = sp.csr_matrix(dense)
    with pytest.raises(ValueError, match=r"1 sample\(s\) have none"):
        bt.da.aldex2(tdata, "group")
    assert fake_rpy2.calls == []


def test_aldex2_keeps_input(fake_rpy2, assert_unchanged):
    fake_rpy2.output = CANNED
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.da.aldex2(tdata, "group", seed=0)
    assert_unchanged(before, tdata)


@pytest.mark.r
def test_aldex2_on_toy_is_the_canned_table():
    out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
    # CANNED holds 6 significant digits of R's output, hence rtol=1e-5.
    np.testing.assert_allclose(out["effect"], CANNED["diff.btw"], rtol=1e-5)
    np.testing.assert_allclose(out["pvalue"], CANNED["we.ep"], rtol=1e-5)


@pytest.mark.r
def test_seed_reproduces_the_table_in_r():
    tdata = bt.datasets.toy()
    first = bt.da.aldex2(tdata, "group", seed=0)
    pd.testing.assert_frame_equal(bt.da.aldex2(tdata, "group", seed=0), first)
    assert not bt.da.aldex2(tdata, "group", seed=1)["effect"].equals(first["effect"])


@pytest.mark.r
def test_reference_sets_the_sign_in_r():
    tdata = bt.datasets.toy()
    a_first, b_first = bt.da.aldex2(tdata, "group", seed=0), bt.da.aldex2(tdata, "group", reference="B", seed=0)
    assert (b_first["contrast"] == "A vs B").all()
    # Swapping the reference flips the sign, but diff.btw's resampling follows the labels, so effects are not exactly
    # opposite; the p-values do not change.
    assert (a_first["direction"] == -b_first["direction"]).all()
    np.testing.assert_allclose(b_first["effect"], -a_first["effect"], atol=0.7)


@pytest.mark.r
def test_all_zero_feature_is_not_tested_in_r():
    tdata = bt.datasets.toy()
    dense = tdata.X.toarray()
    dense[:, tdata.var_names.get_loc("f7")] = 0
    tdata.X = sp.csr_matrix(dense)
    out = bt.da.aldex2(tdata, "group", seed=0)
    assert out.loc["f7", ["effect", "pvalue", "qvalue"]].isna().all()
    assert np.isfinite(out.drop(index="f7")[["effect", "pvalue", "qvalue"]].to_numpy()).all()


@pytest.mark.r
def test_few_mc_samples_warn_in_python():
    with pytest.warns(UserWarning, match=r"da\.aldex2: R warned: .*unreliable"):
        bt.da.aldex2(bt.datasets.toy(), "group", mc_samples=16, seed=0)


@pytest.mark.r
def test_r_messages_are_not_warnings():
    # A direct test of the bridge: R's message() is progress text, not a warning condition.
    from biotapy.da._r import call_r, r_function

    function = r_function('function() { message("hi"); data.frame(a = 1) }', package="base", func="da.x")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert call_r(function)["a"].tolist() == [1.0]


@pytest.mark.r
def test_aldex2_table_passes_the_consensus_checks():
    out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
    assert bt.da.consensus([out], min_methods=1)["n_tested"].tolist() == [1] * 8


@pytest.mark.r
def test_r_warnings_before_an_error_are_in_its_message_in_r():
    from biotapy.da._r import call_r, r_function

    function = r_function('function() { warning("before"); stop("boom") }', package="base", func="da.x")
    with pytest.raises(RuntimeError, match=r"^da\.x: R stopped: boom \(R warned first: before\)$"):
        call_r(function)


def test_single_sample_raises(fake_rpy2):
    with pytest.raises(ValueError, match=r"compares two levels, but obs\['group'\] has 1"):
        bt.da.aldex2(bt.datasets.toy()[:1].copy(), "group")
    assert fake_rpy2.calls == []


@pytest.mark.r
@pytest.mark.parametrize("seeded", [True, False])
def test_r_random_state_is_left_as_found(seeded):
    from rpy2.robjects import r

    drop = 'if (exists(".Random.seed", globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())'
    r("set.seed(1)" if seeded else drop)
    before = list(r(".Random.seed")) if seeded else None
    bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
    assert r('exists(".Random.seed", globalenv(), inherits = FALSE)')[0] == seeded
    if seeded:
        assert list(r(".Random.seed")) == before
