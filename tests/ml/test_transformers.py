import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import make_pipeline
from sklearn.utils.estimator_checks import parametrize_with_checks

import biotapy as bt


# The checks feed random floats below CLR's default pseudocount of 0.5, which CLR warns about
# (tested below); the warning is not what they check.
@pytest.mark.filterwarnings("ignore:pseudocount=")
@parametrize_with_checks([bt.ml.PrevalenceFilter(), bt.ml.CLR()])
def test_scikit_learn_estimator_checks(estimator, check):
    check(estimator)


def _toy_x() -> sp.csr_matrix:
    return bt.datasets.toy().X


def test_prevalence_filter_keeps_what_pp_filter_features_keeps():
    tdata = bt.datasets.toy()
    kept = bt.pp.filter_features(tdata, min_prevalence=0.8).var_names
    selector = bt.ml.PrevalenceFilter(min_prevalence=0.8).fit(tdata.X)
    assert tdata.var_names[selector.get_support()].tolist() == kept.tolist()


def test_prevalence_filter_learns_only_from_the_samples_it_is_fitted_on():
    X = np.array([[1, 0, 3], [2, 0, 1], [0, 5, 2], [1, 7, 0]])
    selector = bt.ml.PrevalenceFilter(min_prevalence=0.5).fit(X[:2])
    np.testing.assert_array_equal(selector.prevalence_, [1.0, 0.0, 1.0])
    np.testing.assert_array_equal(selector.transform(X[2:]), [[0, 2], [1, 0]])


def test_prevalence_filter_boundary_is_inclusive_and_exact():
    X = np.zeros((25, 2))
    X[:7, 0] = 1.0
    X[:, 1] = 1.0
    # 7 / 25 >= 0.28 holds, while 7 >= 0.28 * 25 does not.
    assert bt.ml.PrevalenceFilter(min_prevalence=0.28).fit(X).get_support().tolist() == [True, True]


def test_prevalence_filter_keeps_csr_sparse():
    out = bt.ml.PrevalenceFilter(min_prevalence=1.0).fit_transform(_toy_x())
    assert sp.issparse(out) and out.format == "csr" and out.shape == (6, 2)


def test_prevalence_filter_keeps_feature_names():
    tdata = bt.datasets.toy()
    frame = pd.DataFrame(tdata.X.toarray(), index=tdata.obs_names, columns=tdata.var_names)
    selector = bt.ml.PrevalenceFilter(min_prevalence=1.0).set_output(transform="pandas").fit(frame)
    assert selector.transform(frame).columns.tolist() == ["f3", "f4"]


def test_prevalence_filter_all_zero_feature_is_dropped_and_all_zero_sample_counts():
    X = np.array([[0, 0], [0, 1], [0, 2], [0, 0]])
    selector = bt.ml.PrevalenceFilter(min_prevalence=0.5).fit(X)
    np.testing.assert_array_equal(selector.prevalence_, [0.0, 0.5])
    assert selector.get_support().tolist() == [False, True]


def test_prevalence_filter_single_sample():
    selector = bt.ml.PrevalenceFilter(min_prevalence=1.0).fit(np.array([[0, 3, 1]]))
    assert selector.get_support().tolist() == [False, True, True]


@pytest.mark.parametrize("value", [-0.1, 1.5])
def test_prevalence_filter_out_of_range_raises(value):
    with pytest.raises(ValueError, match=f"min_prevalence must be between 0 and 1, got {value}"):
        bt.ml.PrevalenceFilter(min_prevalence=value).fit(_toy_x())


def test_prevalence_filter_with_nothing_kept_raises():
    with pytest.raises(ValueError, match=r"no feature is non-zero in at least min_prevalence=0\.6 of the 2 samples"):
        bt.ml.PrevalenceFilter(min_prevalence=0.6).fit(np.array([[1, 0], [0, 0]]))


def test_prevalence_filter_stays_unfitted_when_nothing_is_kept():
    selector = bt.ml.PrevalenceFilter(min_prevalence=0.6)
    with pytest.raises(ValueError, match="no feature is non-zero"):
        selector.fit(np.array([[1, 0], [0, 0]]))
    # validate_data has already set n_features_in_, so check_is_fitted would pass; prevalence_ is the learned state.
    assert not hasattr(selector, "prevalence_")


def test_prevalence_filter_with_nothing_kept_prints_the_threshold_unrounded():
    with pytest.raises(ValueError, match=r"min_prevalence=0\.999 "):
        bt.ml.PrevalenceFilter(min_prevalence=0.999).fit(np.array([[1, 0], [0, 1]]))


@pytest.mark.parametrize("value", [True, "0.5", None])
def test_prevalence_filter_wrong_type_raises(value):
    with pytest.raises(TypeError, match="min_prevalence must be a real number"):
        bt.ml.PrevalenceFilter(min_prevalence=value).fit(_toy_x())


def test_prevalence_filter_keeps_its_input():
    X = sp.csr_matrix((np.array([0.0, 2.0, 3.0, 1.0]), np.array([0, 1, 1, 2]), np.array([0, 2, 4])), shape=(2, 3))
    before = (X.data.copy(), X.indices.copy(), X.indptr.copy())
    bt.ml.PrevalenceFilter(min_prevalence=0.5).fit_transform(X)
    for kept, original in zip((X.data, X.indices, X.indptr), before, strict=True):
        np.testing.assert_array_equal(kept, original)


def test_clr_equals_pp_clr():
    tdata = bt.datasets.toy()
    out = bt.ml.CLR().fit_transform(tdata.X)
    np.testing.assert_allclose(out, bt.pp.clr(tdata).layers["clr"], rtol=1e-12)


def test_clr_dense_and_sparse_input_agree():
    X = _toy_x()
    np.testing.assert_allclose(bt.ml.CLR().fit_transform(X.toarray()), bt.ml.CLR().fit_transform(X), rtol=1e-12)


def test_clr_all_zero_sample_is_all_zero():
    out = bt.ml.CLR().fit_transform(np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]))
    np.testing.assert_array_equal(out[0], 0.0)


def test_clr_keeps_its_input():
    X = _toy_x()
    before = X.copy()
    bt.ml.CLR(pseudocount=1).fit_transform(X)
    assert (X != before).nnz == 0


def test_clr_negative_value_raises_in_fit():
    with pytest.raises(ValueError, match="Negative values in data passed to CLR"):
        bt.ml.CLR().fit(np.array([[-1.0, 2.0]]))


@pytest.mark.parametrize(("pseudocount", "error"), [(-1, ValueError), (True, TypeError)])
def test_clr_bad_pseudocount_raises_in_fit(pseudocount, error):
    with pytest.raises(error, match="pseudocount must be"):
        bt.ml.CLR(pseudocount=pseudocount).fit(_toy_x())


def test_clr_pseudocount_above_the_smallest_value_warns():
    relative = np.array([[0.2, 0.8], [0.5, 0.5]])
    with pytest.warns(UserWarning, match=r"pseudocount=0.5 is larger than the smallest non-zero value in X \(0.2\)"):
        bt.ml.CLR().fit_transform(relative)


def test_clr_keeps_feature_names():
    tdata = bt.datasets.toy()
    frame = pd.DataFrame(tdata.X.toarray(), index=tdata.obs_names, columns=tdata.var_names)
    out = bt.ml.CLR().set_output(transform="pandas").fit_transform(frame)
    assert out.columns.tolist() == tdata.var_names.tolist() and out.index.tolist() == tdata.obs_names.tolist()


def test_a_pipeline_refits_the_filter_in_every_fold():
    tdata = bt.datasets.toy()
    pipeline = make_pipeline(bt.ml.PrevalenceFilter(min_prevalence=0.9), bt.ml.CLR(), LogisticRegression())
    result = cross_validate(
        pipeline, tdata.X, tdata.obs["group"], cv=StratifiedKFold(3), return_estimator=True, return_indices=True
    )
    for estimator, train in zip(result["estimator"], result["indices"]["train"], strict=True):
        own = bt.ml.PrevalenceFilter(min_prevalence=0.9).fit(tdata.X[train])
        np.testing.assert_array_equal(estimator[0].prevalence_, own.prevalence_)
    supports = {tuple(estimator[0].get_support()) for estimator in result["estimator"]}
    assert len(supports) > 1


@given(
    arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 6)), elements=st.integers(0, 3)),
    st.floats(0, 1),
)
def test_prevalence_filter_keeps_exactly_the_features_at_or_above_the_threshold(X, threshold):
    present = (X != 0).sum(axis=0) / X.shape[0]
    if not (present >= threshold).any():
        with pytest.raises(ValueError, match="no feature is non-zero"):
            bt.ml.PrevalenceFilter(min_prevalence=threshold).fit(X)
        return
    support = bt.ml.PrevalenceFilter(min_prevalence=threshold).fit(X).get_support()
    np.testing.assert_array_equal(support, present >= threshold)


@given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 1000)))
def test_clr_rows_sum_to_zero(X):
    np.testing.assert_allclose(bt.ml.CLR().fit_transform(X).sum(axis=1), 0.0, atol=1e-9)
