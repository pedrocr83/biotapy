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


def finite_non_negative(X: sp.csr_matrix) -> bool:
    """True when every stored value of ``X`` is finite and >= 0."""
    return bool(np.isfinite(X.data).all() and not (X.data < 0).any())


def divide_rows(X: sp.csr_matrix, totals: npt.NDArray[np.float64]) -> sp.csr_matrix:
    """A float64 copy of ``X`` with each stored value divided by its row's total; zero-total rows stay zero."""
    out = X.astype(np.float64)
    row_totals = np.repeat(totals, np.diff(out.indptr))
    # Divide each value, never multiply by 1 / total: the reciprocal of a subnormal total overflows to inf.
    out.data = np.divide(out.data, row_totals, out=np.zeros_like(out.data), where=row_totals > 0)
    return out


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


def sum_pairs(
    X: sp.csr_matrix, features: npt.NDArray[np.intp], groups: npt.NDArray[np.intp], *, n_groups: int
) -> sp.csr_matrix:
    """Sum the columns of ``X`` into groups given ``(feature, group)`` membership pairs.

    A feature listed with several groups counts in full toward each of them
    (many-to-many, as ``humann_regroup_table`` does). The pairs are a set: a
    repeated pair counts once. A feature in no pair is left out.
    """
    indicator = sp.csr_matrix(
        (np.ones(features.size, dtype=X.dtype), (features, groups)),
        shape=(X.shape[1], n_groups),
    )
    # The constructor sums repeated (feature, group) entries; membership is a set, so reset them to 1.
    indicator.data[:] = 1
    return sp.csr_matrix(X @ indicator)
