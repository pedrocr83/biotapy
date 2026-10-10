import numpy as np
import pytest
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from biotapy._core import argmax_by, as_csr, divide_rows, finite_non_negative, sum_by, sum_pairs

X = sp.csr_matrix(np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int64))


def test_as_csr_converts_dense():
    assert isinstance(as_csr(np.eye(2)), sp.csr_matrix)


def test_as_csr_does_not_copy_csr():
    assert as_csr(X) is X


def test_sum_by_sums_columns_per_group():
    np.testing.assert_array_equal(sum_by(X, np.array([0, 1, 0]), 2).toarray(), [[4, 2], [10, 5]])


def test_sum_by_drops_negative_codes():
    np.testing.assert_array_equal(sum_by(X, np.array([0, -1, 0]), 1).toarray(), [[4], [10]])


def test_sum_by_keeps_integer_dtype():
    assert sum_by(X, np.array([0, 0, 0]), 1).dtype == np.int64


@pytest.mark.parametrize(
    ("dtype", "value", "summed"),
    [
        (np.int8, 100, np.int64),
        (np.uint8, 200, np.uint64),
        (np.int16, 20_000, np.int64),
        (np.bool_, True, np.int64),
        (np.float32, 0.5, np.float32),
    ],
)
def test_sum_by_sums_in_the_dtype_numpy_sums_in(dtype, value, summed):
    # Three columns into one group: 3 * value overflows int8, uint8 and int16, and bool would saturate at True.
    out = sum_by(sp.csr_matrix(np.full((2, 3), value, dtype=dtype)), np.array([0, 0, 0]), 1)
    np.testing.assert_array_equal(out.toarray(), np.full((2, 1), 3 * value))
    assert out.dtype == summed


def test_argmax_by_one_index_per_group_in_code_order():
    np.testing.assert_array_equal(argmax_by(np.array([1.0, 9.0, 3.0, 7.0]), np.array([1, 0, 1, 0])), [1, 2])


def test_argmax_by_takes_first_on_ties():
    np.testing.assert_array_equal(argmax_by(np.array([5.0, 5.0]), np.array([0, 0])), [0])


def test_argmax_by_skips_negative_codes():
    np.testing.assert_array_equal(argmax_by(np.array([9.0, 1.0, 2.0]), np.array([-1, 0, 0])), [2])


def test_argmax_by_with_no_valid_codes_is_empty():
    assert argmax_by(np.array([1.0]), np.array([-1])).size == 0


@given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 50)), st.data())
def test_sum_by_preserves_sample_totals(dense, data):
    codes = np.array(data.draw(st.lists(st.integers(0, 2), min_size=dense.shape[1], max_size=dense.shape[1])))
    out = sum_by(sp.csr_matrix(dense), codes, 3)
    np.testing.assert_array_equal(np.asarray(out.sum(axis=1)).ravel(), dense.sum(axis=1))


def test_sum_pairs_counts_a_feature_in_every_group():
    # f1 belongs to groups 0 and 1, so it counts in full toward both.
    out = sum_pairs(X, np.array([0, 1, 1, 2]), np.array([0, 0, 1, 1]), n_groups=2)
    np.testing.assert_array_equal(out.toarray(), [[3, 5], [9, 11]])


def test_sum_pairs_counts_a_repeated_pair_once():
    out = sum_pairs(X, np.array([0, 0]), np.array([0, 0]), n_groups=1)
    np.testing.assert_array_equal(out.toarray(), [[1], [4]])


def test_sum_pairs_with_no_pairs_is_empty_per_group():
    out = sum_pairs(X, np.array([], dtype=np.intp), np.array([], dtype=np.intp), n_groups=2)
    assert out.shape == (2, 2) and out.nnz == 0


def test_sum_pairs_keeps_integer_dtype():
    assert sum_pairs(X, np.array([0]), np.array([0]), n_groups=1).dtype == np.int64


@given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 50)), st.data())
def test_sum_pairs_with_one_group_per_feature_equals_sum_by(dense, data):
    codes = np.array(data.draw(st.lists(st.integers(-1, 2), min_size=dense.shape[1], max_size=dense.shape[1])))
    kept = np.flatnonzero(codes >= 0)
    expected = sum_by(sp.csr_matrix(dense), codes, 3)
    np.testing.assert_array_equal(
        sum_pairs(sp.csr_matrix(dense), kept, codes[kept], n_groups=3).toarray(), expected.toarray()
    )


@given(
    arrays(np.int64, st.tuples(st.integers(1, 5), st.integers(1, 5)), elements=st.integers(0, 50)),
    st.lists(st.tuples(st.integers(0, 4), st.integers(0, 2)), max_size=12),
)
def test_sum_pairs_many_to_many_matches_a_dense_reference(dense, pairs):
    pairs = [(feature, group) for feature, group in pairs if feature < dense.shape[1]]
    membership = np.zeros((dense.shape[1], 3), dtype=np.int64)
    for feature, group in pairs:
        membership[feature, group] = 1  # a set: a repeated pair stays 1
    features = np.array([feature for feature, _ in pairs], dtype=np.intp)
    groups = np.array([group for _, group in pairs], dtype=np.intp)
    out = sum_pairs(sp.csr_matrix(dense), features, groups, n_groups=3)
    np.testing.assert_array_equal(out.toarray(), dense @ membership)


def test_sum_pairs_keeps_float_dtype():
    out = sum_pairs(sp.csr_matrix(X, dtype=np.float32), np.array([0, 2]), np.array([0, 0]), n_groups=1)
    assert out.dtype == np.float32
    np.testing.assert_array_equal(out.toarray(), [[4], [10]])


def test_sum_pairs_with_no_samples_gives_no_rows():
    out = sum_pairs(sp.csr_matrix((0, 3), dtype=np.int64), np.array([0, 1]), np.array([0, 1]), n_groups=2)
    assert out.shape == (0, 2) and out.nnz == 0


def test_divide_rows_divides_each_value_by_its_row_total():
    X = sp.csr_matrix(np.array([[1, 3], [2, 0]]))
    out = divide_rows(X, np.array([4.0, 2.0]))
    np.testing.assert_array_equal(out.toarray(), [[0.25, 0.75], [1.0, 0.0]])
    assert out.dtype == np.float64


def test_divide_rows_leaves_zero_total_rows_zero():
    out = divide_rows(sp.csr_matrix(np.array([[0.0, 0.0], [1.0, 1.0]])), np.array([0.0, 2.0]))
    np.testing.assert_array_equal(out.toarray(), [[0.0, 0.0], [0.5, 0.5]])


def test_divide_rows_with_a_subnormal_total_stays_finite():
    tiny = np.nextafter(0.0, 1.0)
    out = divide_rows(sp.csr_matrix(np.array([[tiny, 0.0]])), np.array([tiny]))
    np.testing.assert_array_equal(out.toarray(), [[1.0, 0.0]])


def test_divide_rows_does_not_change_its_input():
    X = sp.csr_matrix(np.array([[1.0, 3.0]]))
    divide_rows(X, np.array([4.0]))
    np.testing.assert_array_equal(X.toarray(), [[1.0, 3.0]])


@pytest.mark.parametrize(
    ("data", "expected"),
    [([0.0, 1.5], True), ([], True), ([1.0, -0.1], False), ([1.0, np.nan], False), ([np.inf], False)],
)
def test_finite_non_negative(data, expected):
    X = sp.csr_matrix(
        (np.array(data, dtype=float), (np.zeros(len(data), int), np.arange(len(data)))),
        shape=(1, max(len(data), 1)),
    )
    assert finite_non_negative(X) is expected
