"""LinDA: linear models on log2 centred log-ratios, corrected for compositional bias."""

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd
from anndata import AnnData
from scipy.stats import t as t_dist

from ._design import dense_counts, design_matrix, model
from ._schema import result


def linda(adata: AnnData, group: str, *, covariates: Sequence[str] = (), reference: str | None = None) -> pd.DataFrame:
    """Differential abundance of each feature between two groups by LinDA.

    Parameters
    ----------
    adata
        Samples x features; ``X`` holds raw counts.
    group
        The ``obs`` column whose effect is reported: a categorical, string or bool
        column with two levels, or a numeric column.
    covariates
        ``obs`` columns to adjust for: numeric ones scaled to unit variance, others as
        one indicator per level against their first level.
    reference
        The level of a categorical ``group`` that the other level is compared
        with; by default its first category (sorted values for a string column).

    Returns
    -------
    pandas.DataFrame
        One row per feature, in ``var_names`` order, indexed by ``feature``:
        ``effect`` (log2 fold change, bias-corrected), ``se``, ``pvalue``,
        ``qvalue`` (Benjamini-Hochberg), ``direction`` (sign of ``effect``),
        ``method`` (``"linda"``) and ``contrast`` (``"<level> vs <reference>"``, or
        the column name for a numeric ``group``).

    Raises
    ------
    KeyError
        ``group`` or a covariate is not an ``obs`` column.
    TypeError
        ``covariates`` is not a list of column names, or ``reference`` is not a string.
    ValueError
        ``X`` does not hold raw counts, has an empty sample or fewer than two
        features; a used ``obs`` column has missing values or is constant;
        ``group`` has other than two levels; ``reference`` is not one of them or is
        given for a numeric ``group``; the model has at least as many terms as
        samples, or collinear columns.

    Notes
    -----
    R equivalent: ``MicrobiomeStat::linda``
    Guide: :doc:`/guide/differential_abundance`

    Equals ``MicrobiomeStat::linda(t(counts), meta, "~group + covariates",
    feature.dat.type = "count", is.winsor = FALSE)`` with fixed effects. When ``X``
    holds a zero, 0.5 is added to every count; the log2 counts are centred per
    sample; each feature gets one least-squares fit; the effect's bias is the mode
    of all features' effects, found by a Gaussian mean shift started at the
    shortest half of the data with ``bw.nrd0``'s bandwidth (``modeest::mlv``), and
    is subtracted; p-values are two-sided t-tests on ``n - p`` degrees of freedom.
    Numeric columns are scaled to unit variance first, so a numeric ``group``'s
    effect is per standard deviation.

    MicrobiomeStat 1.4 announces an imputation approach for zeros when library
    size depends on the model, but adds the pseudocount in every case (its switch
    compares ``"Imputation"`` with ``"imputation"``); biotapy computes what it
    returns. Not ported: winsorisation (MicrobiomeStat's default ``is.winsor =
    TRUE``), its prevalence and abundance filters (filter once with
    :func:`biotapy.pp.filter_features` before any method) and random effects.

    The log-ratios have no zeros, so ``X`` is densified once: 8 bytes x samples x
    features. Peak memory is about five such arrays (5.0x measured on 400 x 1,000
    and 200 x 2,000 tables), so budget for that on large tables.

    References
    ----------
    Zhou H, He K, Chen J, Zhang X (2022) LinDA: linear models for differential abundance
    analysis of microbiome compositional data. Genome Biology 23:95.

    Examples
    --------
    >>> import biotapy as bt
    >>> table = bt.da.linda(bt.datasets.toy(), "group")
    >>> round(float(table.loc["f6", "effect"]), 2), table.loc["f6", "contrast"]
    (6.3, 'B vs A')
    """
    frame, contrast = model(adata, group, covariates=covariates, reference=reference, func="da.linda")
    design = design_matrix(frame, scale=True)
    counts = dense_counts(adata, func="da.linda")
    if (counts == 0).any():
        counts += 0.5
    logs = np.log2(counts)
    ratios = logs - logs.mean(axis=1, keepdims=True)
    beta, se, dof = _fit(design, ratios)
    root_n = np.sqrt(adata.n_obs)
    effect = beta - _mode(root_n * beta) / root_n
    pvalue = 2 * t_dist.sf(np.abs(effect / se), dof)
    return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="linda", contrast=contrast)


def _fit(
    design: npt.NDArray[np.float64], ratios: npt.NDArray[np.float64]
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], int]:
    """Least squares of every feature on ``design``: the group coefficient, its standard error, the residual df."""
    coef, *_ = np.linalg.lstsq(design, ratios, rcond=None)
    dof = design.shape[0] - design.shape[1]
    sigma2 = ((ratios - design @ coef) ** 2).sum(axis=0) / dof
    se = np.sqrt(sigma2 * np.linalg.inv(design.T @ design)[1, 1])
    return coef[1], se, dof


def _mode(values: npt.NDArray[np.float64]) -> float:
    """``modeest::mlv(values, method = "meanshift", kernel = "gaussian")``, as modeest 2.4.0 computes it."""
    sd = values.std(ddof=1)
    iqr = np.subtract(*np.percentile(values, [75, 25]))
    # stats::bw.nrd0, including its fallbacks for a zero spread.
    spread = min(sd, iqr / 1.34) or sd or abs(values[0]) or 1.0
    bandwidth = 0.9 * spread * values.size**-0.2
    mode = _shorth(values)
    for _ in range(1000):
        weights = np.exp(-0.5 * ((values - mode) / bandwidth) ** 2)
        shifted = float(values @ weights / weights.sum())
        # modeest stops on a relative step below sqrt(eps) and returns the previous value, not the shifted one.
        # (It divides by a mode of exactly 0 and fails; an unchanged mode stops it either way.)
        if shifted == mode or (mode != 0 and abs(shifted / mode - 1) < np.sqrt(np.finfo(np.float64).eps)):
            break
        mode = shifted
    return mode


def _shorth(values: npt.NDArray[np.float64]) -> float:
    """``modeest::shorth``: the mean of the shortest window holding half of the sorted values."""
    ordered = np.sort(values)
    k = int(np.ceil(ordered.size / 2)) - 1
    widths = ordered[k:] - ordered[: ordered.size - k]
    ties = np.flatnonzero(widths == widths.min()) + 1
    # Tied windows: modeest takes the mean of their 1-based starts, and R's indexing truncates the fraction.
    start = int(ties.mean()) - 1
    return float(ordered[start : start + k + 1].mean())
