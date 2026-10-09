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
