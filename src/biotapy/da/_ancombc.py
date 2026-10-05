"""ANCOM-BC2 through scikit-bio: bias-corrected log abundances, one linear model per feature."""

from collections.abc import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData
from skbio.stats.composition import ancombc2 as skbio_ancombc2

from ._design import dense_counts, model
from ._schema import result


def ancombc2(
    adata: AnnData, group: str, *, covariates: Sequence[str] = (), reference: str | None = None
) -> pd.DataFrame:
    """Differential abundance of each feature between two groups by ANCOM-BC2.

    Parameters
    ----------
    adata
        Samples x features; ``X`` holds raw counts.
    group
        The ``obs`` column whose effect is reported: a categorical, string or bool
        column with two levels, or a numeric column.
    covariates
        ``obs`` columns to adjust for: numeric ones as they are, others as one
        indicator per level against their first level.
    reference
        The level of a categorical ``group`` that the other level is compared
        with; by default its first category (sorted values for a string column).

    Returns
    -------
    pandas.DataFrame
        One row per feature, in ``var_names`` order, indexed by ``feature``:
        ``effect`` (log2 fold change, bias-corrected), ``se``, ``pvalue``,
        ``qvalue`` (Benjamini-Hochberg over the tested features), ``direction``,
        ``method`` (``"ancombc2"``) and ``contrast``. A feature the model cannot
        test has NaN ``effect``, ``se``, ``pvalue`` and ``qvalue`` and
        ``direction`` 0.

    Raises
    ------
    KeyError
        ``group`` or a covariate is not an ``obs`` column.
    TypeError
        ``covariates`` is a string rather than a list of column names.
    ValueError
        ``X`` does not hold raw counts or has an empty sample; a used ``obs`` column
        has missing values; ``group`` has other than two levels; ``reference`` is
        not one of them or is given for a numeric ``group``; the model has as many
        terms as samples, or collinear columns.

    Notes
    -----
    R equivalent: ``ANCOMBC::ancombc2``
    Guide: :doc:`/guide/differential_abundance`

    Runs scikit-bio's :func:`~skbio.stats.composition.ancombc2` (Lin and Peddada
    2024), which matches ``ANCOMBC::ancombc2(..., fix_formula, p_adj_method =
    "BH", prv_cut = 0, pseudo_sens = FALSE)`` with its other defaults: zeros are
    treated as missing (``pseudo = 0``), ``s0_perc = 0.05``, no structural-zero
    test. ANCOM-BC2 reports natural logs; ``effect`` and ``se`` are divided by
    ln 2. A feature whose zeros leave one level of ``group`` (or of a categorical
    covariate) without an observed value cannot be fitted: R and scikit-bio
    report it with p = 1 and count it in the correction, biotapy reports NaN and
    leaves it out, so its q-values are BH over the features actually tested.
    Not run: the pseudocount sensitivity analysis (R's default ``pseudo_sens =
    TRUE``), whose ``passed_ss`` flag the result table does not carry, and R's
    prevalence filter (``prv_cut = 0.10``; filter once with
    :func:`biotapy.pp.filter_features` before any method).

    The bias is estimated by an E-M algorithm capped at 100 iterations, R's
    ``em_control`` default, which biotapy keeps. On some data it has not
    converged by then, and the estimate depends on the cap: on GlobalPatterns'
    genera, human hosts against the rest, R moves every effect by about -0.28
    log2 with 1,000 iterations, and calls 220 genera instead of 208.

    scikit-bio needs a dense table, so ``X`` is densified once: 8 bytes x
    samples x features. Peak memory is about seven such arrays (7.2x to 7.5x
    measured on 400 x 1,000 and 200 x 2,000 tables).

    References
    ----------
    Lin H, Peddada SD (2024) Multigroup analysis of compositions of microbiomes with
    covariate adjustments and repeated measures. Nature Methods 21:83-91.

    Examples
    --------
    >>> import biotapy as bt
    >>> table = bt.da.ancombc2(bt.datasets.toy(), "group")
    >>> round(float(table.loc["f6", "effect"]), 2), table.loc["f6", "contrast"]
    (3.76, 'B vs A')
    """
    frame, contrast = model(adata, group, covariates=covariates, reference=reference, func="da.ancombc2")
    counts = pd.DataFrame(dense_counts(adata, func="da.ancombc2"), index=adata.obs_names, columns=adata.var_names)
    # Plain names, so patsy needs no quoting; categorical columns keep their levels, reference first.
    names = [f"x{i}" for i in range(frame.shape[1])]
    fit = skbio_ancombc2(counts, frame.set_axis(names, axis=1), " + ".join(names), p_adjust=None).result
    covariate = next(name for name in fit.index.unique("Covariate") if name == "x0" or name.startswith("x0["))
    table = fit.xs(covariate, level="Covariate").reindex(adata.var_names)
    effect = table["Log(FC)"].to_numpy(np.float64) / np.log(2)
    # scikit-bio reports a feature it could not fit with p = 1; the schema says "not tested" with NaN.
    pvalue = np.where(np.isnan(effect), np.nan, table["pvalue"].to_numpy(np.float64))
    se = table["SE"].to_numpy(np.float64) / np.log(2)
    return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="ancombc2", contrast=contrast)
