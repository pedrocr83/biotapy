import numpy as np
import pytest
import scipy.sparse as sp

from biotapy._core import check_pseudocount, pseudocounted


def test_dense_and_sparse_x_give_the_same_values():
    dense = np.array([[0.0, 2.0], [3.0, 1.0]])
    expected = dense + 0.5
    np.testing.assert_array_equal(pseudocounted(dense, 0.5, func="f"), expected)
    np.testing.assert_array_equal(pseudocounted(sp.csr_matrix(dense), 0.5, func="f"), expected)


def test_columns_are_taken_in_order():
    dense = np.array([[1.0, 2.0, 3.0]])
    np.testing.assert_array_equal(pseudocounted(dense, 1, func="f", columns=np.array([2, 0])), [[4.0, 2.0]])


def test_does_not_change_its_input():
    dense = np.array([[1.0, 2.0]])
    pseudocounted(dense, 0.5, func="f")
    np.testing.assert_array_equal(dense, [[1.0, 2.0]])


@pytest.mark.parametrize(
    ("value", "error"), [(True, TypeError), ("1", TypeError), (-1, ValueError), (np.inf, ValueError)]
)
def test_check_pseudocount_refuses(value, error):
    with pytest.raises(error, match="pseudocount must be"):
        check_pseudocount(value)


def test_messages_name_the_caller():
    with pytest.raises(ValueError, match="ml.CLR needs finite, non-negative values in X"):
        pseudocounted(np.array([[-1.0, 2.0]]), 0.5, func="ml.CLR")


def test_pseudocount_above_the_smallest_value_warns_naming_it():
    values = np.array([[0.0, 0.2], [0.8, 0.5]])
    with pytest.warns(UserWarning, match=r"pseudocount=0\.5 is larger than the smallest non-zero value in X \(0\.2\)"):
        pseudocounted(values, 0.5, func="ml.CLR")


def test_pseudocount_warning_ends_with_the_calling_step():
    with pytest.warns(UserWarning, match=r"on their scale \(ml\.CLR\)$"):
        pseudocounted(np.array([[0.0, 0.2], [0.8, 0.5]]), 0.5, func="ml.CLR")


def test_zero_pseudocount_with_a_zero_in_x_raises():
    with pytest.raises(ValueError, match="X holds zeros, whose logarithm is undefined; pass pseudocount > 0 to ml.CLR"):
        pseudocounted(np.array([[0.0, 2.0]]), 0, func="ml.CLR")
