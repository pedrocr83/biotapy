"""scikit-learn transformers, so preprocessing is fitted inside each cross-validation fold."""

from typing import Self

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from skbio.stats.composition import clr as skbio_clr
from sklearn.base import BaseEstimator, OneToOneFeatureMixin, TransformerMixin
from sklearn.feature_selection import SelectorMixin
from sklearn.utils import Tags
from sklearn.utils.validation import check_is_fitted, check_non_negative, validate_data

from biotapy._core import as_csr, check_pseudocount, pseudocounted

# What fit and transform take: a samples x features array, sparse matrix or DataFrame.
Table = npt.ArrayLike | sp.spmatrix


class PrevalenceFilter(SelectorMixin, BaseEstimator):
    """Keep the features that are non-zero in enough of the training samples.

    Parameters
    ----------
    min_prevalence
        Keep features that are non-zero in at least this fraction of the samples
        given to ``fit``, from 0 to 1.

    Attributes
    ----------
    prevalence_
        Each feature's fraction of training samples in which it is non-zero.
    n_features_in_
        Number of features seen in ``fit``.
    feature_names_in_
        Feature names seen in ``fit``, when it was given a DataFrame.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/machine_learning`

    The scikit-learn form of :func:`biotapy.pp.filter_features` with
    ``min_prevalence``, keeping the same features. Which features pass depends
    on the samples, so filtering before splitting lets the test samples choose
    the features: inside a :class:`~sklearn.pipeline.Pipeline` the filter is
    fitted on each training fold only. A sparse ``X`` stays sparse.

    Examples
    --------
    >>> import biotapy as bt
    >>> X = bt.datasets.toy().X
    >>> bt.ml.PrevalenceFilter(min_prevalence=1.0).fit_transform(X).shape
    (6, 2)
    """

    def __init__(self, min_prevalence: float = 0.1) -> None:
        self.min_prevalence = min_prevalence

    def fit(self, X: Table, y: object = None) -> Self:
        """Learn each feature's prevalence in ``X``, samples x features."""
        X = validate_data(self, X, accept_sparse="csr")
        if not 0 <= self.min_prevalence <= 1:
            msg = f"min_prevalence must be between 0 and 1, got {self.min_prevalence}"
            raise ValueError(msg)
        matrix = as_csr(X)
        present = np.bincount(matrix.indices[matrix.data != 0], minlength=matrix.shape[1])
        self.prevalence_ = present / matrix.shape[0]
        if not self._get_support_mask().any():
            msg = (
                f"no feature is non-zero in at least {self.min_prevalence:.0%} of the "
                f"{matrix.shape[0]} samples; lower min_prevalence"
            )
            raise ValueError(msg)
        return self

    def _get_support_mask(self) -> npt.NDArray[np.bool_]:
        check_is_fitted(self)
        # Divide rather than multiply, as pp.filter_features does: 7 / 25 >= 0.28 holds, 7 >= 0.28 * 25 does not.
        return np.asarray(self.prevalence_ >= self.min_prevalence)

    def __sklearn_tags__(self) -> Tags:
        # scikit-learn's tag methods are unannotated, and mypy's untyped_calls_exclude does not reach super().
        tags: Tags = super().__sklearn_tags__()  # type: ignore[no-untyped-call]
        tags.input_tags.sparse = True
        return tags


class CLR(OneToOneFeatureMixin, TransformerMixin, BaseEstimator):
    """Centred log-ratio transform of each sample, as a scikit-learn transformer.

    Parameters
    ----------
    pseudocount
        Added to every value before the logarithm, so zeros have one.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/machine_learning`

    Gives the values of :func:`biotapy.pp.clr`. Each sample is transformed on
    its own, so ``fit`` learns nothing (it only checks ``X``); the class exists
    so the transform can sit in a :class:`~sklearn.pipeline.Pipeline`. The
    output is a dense float64 array: 8 bytes x samples x features.

    Examples
    --------
    >>> import biotapy as bt
    >>> out = bt.ml.CLR().fit_transform(bt.datasets.toy().X)
    >>> round(float(out[0].sum()), 6)
    0.0
    """

    def __init__(self, pseudocount: float = 0.5) -> None:
        self.pseudocount = pseudocount

    def fit(self, X: Table, y: object = None) -> Self:
        """Check ``X`` (samples x features, non-negative) and ``pseudocount``; nothing is learned."""
        X = validate_data(self, X, accept_sparse="csr")
        check_pseudocount(self.pseudocount)
        check_non_negative(X, "CLR")
        return self

    def transform(self, X: Table) -> npt.NDArray[np.float64]:
        """Return the CLR of each sample of ``X`` as a dense float64 array."""
        check_is_fitted(self)
        X = validate_data(self, X, accept_sparse="csr", reset=False)
        return np.asarray(skbio_clr(pseudocounted(X, self.pseudocount, func="ml.CLR")), dtype=np.float64)

    def __sklearn_tags__(self) -> Tags:
        # scikit-learn's tag methods are unannotated, and mypy's untyped_calls_exclude does not reach super().
        tags: Tags = super().__sklearn_tags__()  # type: ignore[no-untyped-call]
        tags.input_tags.sparse = True
        tags.input_tags.positive_only = True
        tags.requires_fit = False
        return tags
