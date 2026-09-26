"""Sparse kernels shared by pp, fn and tl."""

from typing import Any, cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp


def as_csr(X: object) -> sp.csr_matrix:
    """Return ``X`` as a CSR matrix, without copying one that already is.

    ``X`` may be anything ``scipy.sparse.csr_matrix`` accepts, including
    ``AnnData.X`` (typed by anndata as a private union of array/sparse types).
    """
    if isinstance(X, sp.csr_matrix):
        return X
    # scipy-stubs' csr_matrix overloads cannot resolve an arbitrary object argument;
    # the cast is annotation-only, the runtime call is unchanged.
    return sp.csr_matrix(cast(Any, X))


def sum_by(X: sp.csr_matrix, codes: npt.NDArray[np.intp], n_groups: int) -> sp.csr_matrix:
    """Sum the columns of ``X`` that share a group code; negative codes are dropped."""
    rows = np.flatnonzero(codes >= 0)
    indicator = sp.csr_matrix(
        (np.ones(rows.size, dtype=X.dtype), (rows, codes[rows])),
        shape=(codes.size, n_groups),
    )
    return sp.csr_matrix(X @ indicator)


def argmax_by(values: npt.NDArray[np.float64], codes: npt.NDArray[np.intp]) -> npt.NDArray[np.intp]:
    """Index of the largest value per group, first on ties, ordered by group code."""
    valid = np.flatnonzero(codes >= 0)
    if valid.size == 0:
        return valid
    order = valid[np.lexsort((-values[valid], codes[valid]))]
    first = np.r_[True, np.diff(codes[order]) != 0]
    return order[first]
