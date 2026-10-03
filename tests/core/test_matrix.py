import numpy as np
import scipy.sparse as sp
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from biotapy._core import argmax_by, as_csr, sum_by, sum_pairs

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
