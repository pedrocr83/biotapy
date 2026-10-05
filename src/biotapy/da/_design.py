"""The model every da method fits: ``group``, ``covariates`` and ``reference`` checked once, without formulas."""

from collections.abc import Sequence
from typing import cast

import numpy as np
import numpy.typing as npt
import pandas as pd
from anndata import AnnData

from biotapy._core import as_csr, require_counts


def model(
    adata: AnnData, group: str, *, covariates: Sequence[str], reference: str | None, func: str
) -> tuple[pd.DataFrame, str]:
    """``obs[[group, *covariates]]`` ready to fit, and the ``contrast`` text.

    Numeric columns become float64; any other column a ``Categorical`` without
    unused levels, the group's ``reference`` first. Raises when ``covariates`` is
    not a list of column names, a column is missing, has missing values or is
    constant, the group does not have two levels, or the design has no residual
    degrees of freedom or collinear columns.
    """
    if (
        not isinstance(covariates, Sequence)
        or isinstance(covariates, str)
        or not all(isinstance(c, str) for c in covariates)
    ):
        example = covariates if isinstance(covariates, str) else "age"
        msg = f"{func}: covariates must be a list of obs columns, such as [{example!r}]"
        raise TypeError(msg)
    names = [group, *covariates]
    absent = [name for name in names if name not in adata.obs.columns]
    if absent:
        msg = f"{func}: {absent} not in obs; group and covariates name obs columns"
        raise KeyError(msg)
    if len(set(names)) < len(names):
        msg = f"{func}: group={group!r} and covariates={list(covariates)} repeat a column"
        raise ValueError(msg)
    # anndata types obs as DataFrame | Dataset2D (its lazy variant); the data model guarantees a DataFrame.
    obs = cast("pd.DataFrame", adata.obs)
    frame = pd.DataFrame({name: _column(obs[name], func=func) for name in names}, index=obs.index)
    frame[group], contrast = _group(frame[group], reference, func=func)
    design = design_matrix(frame, scale=False)
    if design.shape[0] <= design.shape[1]:
        msg = f"{func} needs more samples ({design.shape[0]}) than model terms ({design.shape[1]}, with the intercept)"
        raise ValueError(msg)
    _check_varies(frame, group, func=func)
    if np.linalg.matrix_rank(design) < design.shape[1]:
        msg = f"{func}: group={group!r} and covariates={list(covariates)} are collinear; drop the covariate that repeats another"
        raise ValueError(msg)
    return frame, contrast


def design_matrix(frame: pd.DataFrame, *, scale: bool) -> npt.NDArray[np.float64]:
    """Intercept, then each column of ``frame``: indicators against the first category, or the values.

    With ``scale=True`` numeric columns are centred and divided by their standard
    deviation, as R's ``scale()``. The group's column is always column 1.
    """
    blocks = [np.ones((len(frame), 1))]
    for name in frame.columns:
        values = frame[name]
        if isinstance(values.dtype, pd.CategoricalDtype):
            codes = values.cat.codes.to_numpy()
            blocks.append((codes[:, None] == np.arange(1, len(values.cat.categories))).astype(np.float64))
        elif scale:
            blocks.append(((values - values.mean()) / values.std()).to_numpy(np.float64)[:, None])
        else:
            blocks.append(values.to_numpy(np.float64)[:, None])
    return np.hstack(blocks)


def dense_counts(adata: AnnData, *, func: str) -> npt.NDArray[np.float64]:
    """``X`` as a dense float64 array of raw counts with no empty sample."""
    require_counts(adata, func=func)
    if adata.n_vars < 2:
        msg = f"{func} needs at least two features: log-ratio methods compare each feature with the others"
        raise ValueError(msg)
    X = as_csr(adata.X)
    empty = adata.obs_names[np.asarray(X.sum(axis=1)).ravel() == 0]
    if len(empty):
        msg = (
            f"{func} needs reads in every sample; {len(empty)} sample(s) have none ({empty[:3].tolist()}): "
            "drop them first, for example with bt.pp.filter_samples(adata, min_depth=1)"
        )
        raise ValueError(msg)
    # LinDA's log-ratios and scikit-bio's ancombc2 need a dense table (rules.md R6.2): the one dense copy of X.
    return X.toarray().astype(np.float64, copy=False)


def _check_varies(frame: pd.DataFrame, group: str, *, func: str) -> None:
    """Raise for a constant column; a one-level categorical group is left to ``_group``'s two-level message."""
    for name in frame.columns:
        if frame[name].nunique() > 1 or (name == group and isinstance(frame[name].dtype, pd.CategoricalDtype)):
            continue
        if name == group:
            msg = f"{func}: group column {name!r} is constant across samples; there is nothing to compare"
        else:
            msg = f"{func}: covariates column {name!r} is constant across samples; drop it"
        raise ValueError(msg)


def _column(values: pd.Series, *, func: str) -> pd.Series:
    """Float64 for a numeric column, else a ``Categorical`` of the levels present; missing values raise."""
    if values.isna().any():
        msg = f"{func}: obs[{values.name!r}] is missing for {int(values.isna().sum())} sample(s); drop them first"
        raise ValueError(msg)
    if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
        return values.astype(np.float64)
    return pd.Series(pd.Categorical(values).remove_unused_categories(), index=values.index, name=values.name)


def _group(values: pd.Series, reference: str | None, *, func: str) -> tuple[pd.Series, str]:
    """The group column with ``reference`` as its first level, and the contrast text."""
    if not isinstance(values.dtype, pd.CategoricalDtype):
        if reference is not None:
            msg = f"{func}: reference={reference!r} needs a categorical group, but obs[{values.name!r}] is numeric"
            raise ValueError(msg)
        return values, str(values.name)
    levels = [str(level) for level in values.cat.categories]
    if len(levels) != 2:
        msg = f"{func} compares two levels, but obs[{values.name!r}] has {len(levels)}: {levels}; subset the samples first"
        raise ValueError(msg)
    if reference is None:
        reference = levels[0]
    if not isinstance(reference, str):
        msg = f"{func}: reference must be the level's name as a string, such as {levels[0]!r}, not {reference!r}"
        raise TypeError(msg)
    if reference not in levels:
        msg = f"{func}: reference={reference!r} is not a level of obs[{values.name!r}], which has {levels}"
        raise ValueError(msg)
    first = values.cat.categories[levels.index(reference)]
    ordered = values.cat.reorder_categories([first, *[level for level in values.cat.categories if level != first]])
    return ordered, f"{levels[1 - levels.index(reference)]} vs {reference}"
