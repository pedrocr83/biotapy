import sys
import types

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

import biotapy as bt

GROUP = "biotapy.embeddings"


def _installer(patch, root):
    """Install embedding plugins as a package would: a distribution's entry points under ``root`` on sys.path."""
    module = types.ModuleType(f"biotapy_test_plugins_{root.name}")
    patch.setitem(sys.modules, module.__name__, module)
    patch.syspath_prepend(str(root))
    registered = {}

    def install(name, plugin, *, distribution="biotapy-test-plugins"):
        attribute = f"plugin_{len(vars(module))}"
        setattr(module, attribute, plugin)
        registered.setdefault(distribution, []).append(f"{name} = {module.__name__}:{attribute}")
        info = root / f"{distribution.replace('-', '_')}-1.0.dist-info"
        info.mkdir(exist_ok=True)
        (info / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {distribution}\nVersion: 1.0\n", encoding="utf-8")
        lines = "\n".join(registered[distribution])
        (info / "entry_points.txt").write_text(f"[{GROUP}]\n{lines}\n", encoding="utf-8")

    return install


@pytest.fixture
def install(tmp_path, monkeypatch):
    return _installer(monkeypatch, tmp_path)


@pytest.fixture(scope="module")
def echo(tmp_path_factory):
    """A plugin that returns ``uns["echo"]``, installed once per module: Hypothesis rejects function-scoped fixtures."""
    with pytest.MonkeyPatch.context() as patch:
        _installer(patch, tmp_path_factory.mktemp("echo"))("echo", lambda adata: adata.uns["echo"])
        yield


def _ones(adata):
    return np.ones((adata.n_obs, 3))


def test_returns_the_plugin_s_embedding(install):
    install("fake", lambda adata: np.arange(adata.n_obs * 2, dtype=float).reshape(-1, 2))
    tdata = bt.datasets.toy()
    result = bt.ml.embed(tdata, "fake")
    np.testing.assert_array_equal(result, np.arange(12, dtype=float).reshape(6, 2))
    assert "X_fake" not in tdata.obsm


def test_inplace_writes_obsm_x_model_and_returns_none(install):
    install("fake", _ones)
    tdata = bt.datasets.toy()
    assert bt.ml.embed(tdata, "fake", inplace=True) is None
    np.testing.assert_array_equal(tdata.obsm["X_fake"], np.ones((6, 3)))


def test_the_plugin_gets_the_data_itself(install):
    calls = []
    install("fake", lambda adata: calls.append(adata) or np.ones((adata.n_obs, 1)))
    tdata = bt.datasets.toy()
    bt.ml.embed(tdata, "fake")
    assert len(calls) == 1 and calls[0] is tdata


def test_unknown_model_lists_the_installed_ones(install):
    install("fake", _ones)
    install("other", _ones)
    with pytest.raises(
        KeyError, match=r"model='nope' is not an installed embedding plugin; installed: \[.*'fake', .*'other'"
    ):
        bt.ml.embed(bt.datasets.toy(), "nope")


@pytest.mark.parametrize("model", ["", "a/b", "with space", "x=y", "../up"])
def test_a_model_name_that_is_not_an_obsm_key_raises(install, model):
    install("fake", _ones)
    with pytest.raises(ValueError, match=r"model=.* must be letters, digits, '_', '-' or '\.'"):
        bt.ml.embed(bt.datasets.toy(), model)


@pytest.mark.parametrize("model", ["fake", "My.Model_2", "mgm-v2"])
def test_a_model_name_of_letters_digits_and_separators_is_accepted(install, model):
    install(model, _ones)
    assert bt.ml.embed(bt.datasets.toy(), model).shape == (6, 3)


def test_a_package_registering_one_name_twice_is_named_once(install):
    install("fake", _ones, distribution="only-package")
    install("fake", _ones, distribution="only-package")
    with pytest.raises(ValueError, match=r"several packages \(only-package\); uninstall all but one"):
        bt.ml.embed(bt.datasets.toy(), "fake")


def test_two_packages_registering_one_name_raise(install):
    install("fake", _ones, distribution="first-package")
    install("fake", _ones, distribution="second-package")
    with pytest.raises(
        ValueError, match=r"model='fake' is registered by several packages \(first-package, second-package\)"
    ):
        bt.ml.embed(bt.datasets.toy(), "fake")


def test_plugin_with_wrong_row_count_raises(install):
    install("fake", lambda adata: np.ones((adata.n_obs - 1, 3)))
    with pytest.raises(ValueError, match=r"plugin 'fake' returned shape \(5, 3\); expected \(6, dimensions\)"):
        bt.ml.embed(bt.datasets.toy(), "fake")


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_plugin_with_nan_raises(install, bad):
    def plugin(adata):
        out = np.ones((adata.n_obs, 3))
        out[2, 1] = bad
        return out

    install("fake", plugin)
    with pytest.raises(ValueError, match="plugin 'fake' returned NaN or infinite values"):
        bt.ml.embed(bt.datasets.toy(), "fake")


@pytest.mark.parametrize("shape", [(6,), (6, 0), (6, 2, 2)])
def test_a_result_that_is_not_samples_by_dimensions_raises(install, shape):
    install("fake", lambda adata: np.ones(shape))
    with pytest.raises(ValueError, match=r"plugin 'fake' returned shape"):
        bt.ml.embed(bt.datasets.toy(), "fake")


def test_an_integer_result_raises(install):
    install("fake", lambda adata: np.ones((adata.n_obs, 3), dtype=np.int64))
    with pytest.raises(ValueError, match="plugin 'fake' returned int64 values; expected floats"):
        bt.ml.embed(bt.datasets.toy(), "fake")


def test_a_result_that_is_not_an_array_raises(install):
    install("fake", lambda adata: [[1.0, 2.0]] * adata.n_obs)
    with pytest.raises(TypeError, match="plugin 'fake' returned a list, not a NumPy array"):
        bt.ml.embed(bt.datasets.toy(), "fake")


@pytest.mark.parametrize(
    ("make", "kind"),
    [
        (lambda n: np.ma.masked_invalid(np.full((n, 3), np.nan)), "MaskedArray"),
        (lambda n: np.ones((n, 3)).view(np.matrix), "matrix"),
    ],
)
def test_an_array_subclass_raises(install, make, kind):
    install("fake", lambda adata: make(adata.n_obs))
    with pytest.raises(TypeError, match=f"plugin 'fake' returned a {kind}, not a NumPy array"):
        bt.ml.embed(bt.datasets.toy(), "fake")


@pytest.mark.parametrize("where", ["X", "layer", "obsm", "uns"])
def test_a_result_that_shares_memory_with_the_input_is_copied(install, where):
    tdata = bt.datasets.toy()
    tdata.X = np.ones((tdata.n_obs, tdata.n_vars))
    tdata.layers["own"] = np.ones((tdata.n_obs, tdata.n_vars))
    tdata.obsm["own"] = np.ones((tdata.n_obs, 3))
    tdata.uns["own"] = np.ones((tdata.n_obs, 3))
    held = {
        "X": lambda a: a.X,
        "layer": lambda a: a.layers["own"],
        "obsm": lambda a: a.obsm["own"],
        "uns": lambda a: a.uns["own"],
    }
    install("fake", held[where])
    result = bt.ml.embed(tdata, "fake")
    assert not np.shares_memory(result, held[where](tdata))
    np.testing.assert_array_equal(result, held[where](tdata))


def test_a_plugin_that_raises_keeps_its_error_and_is_named(install):
    def plugin(adata):
        msg = "model exploded"
        raise RuntimeError(msg)

    install("fake", plugin)
    with pytest.raises(RuntimeError, match="model exploded") as info:
        bt.ml.embed(bt.datasets.toy(), "fake")
    assert any("embedding plugin 'fake'" in note for note in info.value.__notes__)


def test_a_plugin_that_is_not_callable_is_named(install):
    install("fake", 3)
    with pytest.raises(TypeError, match="not callable") as info:
        bt.ml.embed(bt.datasets.toy(), "fake")
    assert any("embedding plugin 'fake'" in note for note in info.value.__notes__)


def test_a_plugin_that_cannot_be_loaded_is_named(install, tmp_path):
    install("fake", _ones)
    module = sys.modules[f"biotapy_test_plugins_{tmp_path.name}"]
    for attribute in [name for name in vars(module) if name.startswith("plugin_")]:
        delattr(module, attribute)
    with pytest.raises(AttributeError) as info:
        bt.ml.embed(bt.datasets.toy(), "fake")
    assert any("embedding plugin 'fake'" in note for note in info.value.__notes__)


def test_a_refused_result_is_not_written(install):
    install("fake", lambda adata: np.full((adata.n_obs, 3), np.nan))
    tdata = bt.datasets.toy()
    with pytest.raises(ValueError, match="NaN"):
        bt.ml.embed(tdata, "fake", inplace=True)
    assert "X_fake" not in tdata.obsm


def test_keeps_the_input(install, assert_unchanged):
    install("fake", _ones)
    tdata = bt.datasets.toy()
    before = tdata.copy()
    bt.ml.embed(tdata, "fake")
    assert_unchanged(before, tdata)


def test_a_feature_change_drops_the_embedding(install):
    install("fake", _ones)
    tdata = bt.datasets.toy()
    bt.ml.embed(tdata, "fake", inplace=True)
    assert "X_fake" not in bt.pp.filter_features(tdata, min_prevalence=0.5).obsm


def test_single_sample_and_all_zero_sample(install, make_adata):
    install("fake", lambda adata: np.asarray(adata.X.sum(axis=1), dtype=float))
    np.testing.assert_array_equal(bt.ml.embed(make_adata(np.array([[0, 0, 0]])), "fake"), [[0.0]])
    np.testing.assert_array_equal(bt.ml.embed(make_adata(np.array([[0, 0], [2, 3]])), "fake"), [[0.0], [5.0]])


def test_options_are_keyword_only(install):
    install("fake", _ones)
    with pytest.raises(TypeError):
        bt.ml.embed(bt.datasets.toy(), "fake", True)


@given(
    arrays(
        np.float64,
        st.tuples(st.just(6), st.integers(1, 5)),
        elements=st.floats(allow_nan=True, allow_infinity=True) | st.floats(-1e6, 1e6),
    )
)
def test_a_finite_result_comes_back_unchanged_and_any_other_raises(echo, result):
    tdata = bt.datasets.toy()
    tdata.uns["echo"] = result
    if np.isfinite(result).all():
        np.testing.assert_array_equal(bt.ml.embed(tdata, "echo"), result)
    else:
        with pytest.raises(ValueError, match="plugin 'echo' returned NaN or infinite values"):
            bt.ml.embed(tdata, "echo")
