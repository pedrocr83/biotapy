"""MaAsLin 3 through the R bridge: its abundance model, linear models of log2 relative abundance where present."""

from collections.abc import Sequence

import numpy as np
import pandas as pd
from anndata import AnnData

from ._design import dense_counts, model
from ._r import call_r, r_function, r_seed
from ._schema import result

# evaluate_only = "abundance" needs warn_prevalence = FALSE (maaslin3 stops otherwise). subtract_median = TRUE reports
# each coefficient minus the median its p-value is tested against. The output folder maaslin3 insists on is deleted,
# and pbapply's progress bar, which verbosity does not reach, is off for the call.
_MAASLIN = """function(counts, metadata, formula, seed) {
  output <- tempfile("biotapy-maaslin3-")
  progress <- pbapply::pboptions(type = "none")
  on.exit({unlink(output, recursive = TRUE); pbapply::pboptions(progress)})
  set.seed(seed)
  fit <- maaslin3::maaslin3(
    counts, metadata, output, formula = formula, min_abundance = 0, min_prevalence = 0, evaluate_only = "abundance",
    warn_prevalence = FALSE, subtract_median = TRUE, plot_summary_plot = FALSE, plot_associations = FALSE,
    cores = 1, verbosity = "ERROR"
  )
  out <- fit$fit_data_abundance$results
  data.frame(feature = out$feature, metadata = out$metadata, coef = out$coef, stderr = out$stderr,
             pval = out$pval_individual, failed = !is.na(out$error))
}"""


def maaslin3(
    adata: AnnData,
    group: str,
    *,
    covariates: Sequence[str] = (),
    reference: str | None = None,
    seed: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """Differential abundance of each feature between two groups by MaAsLin 3's abundance model, run in R.

    Parameters
    ----------
    adata
        Samples x features; ``X`` holds raw counts.
    group
        The ``obs`` column whose effect is reported: a categorical, string or bool
        column with two levels, or a numeric column, whose effect is per standard
        deviation (MaAsLin 3 standardises numeric columns).
    covariates
        ``obs`` columns to adjust for: numeric ones standardised, others as one
        indicator per level against their first level.
    reference
        The level of a categorical ``group`` that the other level is compared
        with; by default its first category (sorted values for a string column).
    seed
        Seeds R's random number generator through ``set.seed`` for the call; MaAsLin 3's test
        against the median simulates, so the same seed gives the same table.

    Returns
    -------
    pandas.DataFrame
        One row per feature, in ``var_names`` order, indexed by ``feature``:
        ``effect`` (log2 fold change of the relative abundance where the feature
        is present, minus the median described in Notes), ``se``, ``pvalue``,
        ``qvalue`` (Benjamini-Hochberg over the tested features), ``direction``,
        ``method`` (``"maaslin3"``) and ``contrast``. A feature MaAsLin 3 could
        not fit has NaN values and ``direction`` 0.

    Raises
    ------
    ImportError
        rpy2 (the extra ``r``) or the R package maaslin3 is not installed.
    KeyError
        ``group`` or a covariate is not an ``obs`` column.
    TypeError
        ``covariates`` is not a list of column names, or ``reference`` is not a string.
    RuntimeError
        R stops with an error.
    ValueError
        ``X`` does not hold raw counts, has an empty sample, repeated
        ``var_names`` or ``obs_names``, or fewer than two
        features; a used ``obs`` column has missing values, is constant or is
        repeated; ``group`` has other than two levels; ``reference`` is not one of
        them or is given for a numeric ``group``; the model has at least as many
        terms as samples, or collinear columns.

    Notes
    -----
    R equivalent: ``maaslin3::maaslin3``
    Guide: :doc:`/guide/differential_abundance`

    Calls ``maaslin3::maaslin3(counts, metadata, output, formula, min_prevalence
    = 0, evaluate_only = "abundance", warn_prevalence = FALSE, subtract_median =
    TRUE)`` (Nickols et al.) through rpy2, with its other defaults: counts
    are divided by each sample's total (TSS), zeros are left out and the rest
    log2-transformed, one linear model per feature is fitted on the samples where
    it is present, and each coefficient is tested against the median coefficient
    (``median_comparison_abundance = TRUE``), MaAsLin 3's correction for
    compositionality, which draws 10,000 normal samples. That median is taken, as
    MaAsLin 3 computes it, over the features without a fit error whose own
    p-value is below 0.95. ``effect`` is that coefficient minus the median, so
    its sign is the side of the median the test is about; MaAsLin 3's default
    output reports the coefficient itself. Only the abundance model runs: the
    prevalence model's log-odds cannot share a column with fold changes. ``qvalue`` is the
    Benjamini-Hochberg correction of the group's p-values; MaAsLin 3's
    ``qval_individual`` corrects them together with every covariate's. A feature
    whose fit reports an error is not tested, as MaAsLin 3 leaves it out of its
    own correction. Swapping ``reference`` negates ``effect`` exactly, but the
    simulation draws around the coefficients, not their negatives, so the same
    ``seed`` moves ``pvalue`` by up to 1e-3 and ``qvalue`` by up to 2e-3 on the
    GlobalPatterns genera (the calls are the same there). Each R warning raised
    during the call is re-emitted as a ``UserWarning``, and an R error is raised
    as a ``RuntimeError``.

    Needs R with maaslin3 (``BiocManager::install("maaslin3")``) and ``pip
    install 'biotapy[r]'``, which builds rpy2 (GPL-2.0-or-later) against that R.
    ``seed`` becomes one integer for R's ``set.seed``; the call puts the
    embedded R session's random state back as it found it. Plots are off and MaAsLin 3's
    output folder is a temporary one, deleted afterwards.

    rpy2 converts dense tables only, so ``X`` is densified once: 8 bytes x
    samples x features, plus R's copies of the table.

    References
    ----------
    Nickols WA, et al. MaAsLin 3: refining and extending generalized multivariable
    linear models for meta-omic association discovery. Nature Methods,
    doi:10.1038/s41592-025-02923-9.

    Examples
    --------
    >>> import biotapy as bt
    >>> table = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)  # doctest: +SKIP
    >>> table.loc["f6", "contrast"]  # doctest: +SKIP
    'B vs A'
    """
    frame, contrast = model(adata, group, covariates=covariates, reference=reference, func="da.maaslin3")
    seed_for_r = r_seed(seed)
    # rpy2 has no sparse converter (rules.md R6.2): the one dense copy of X.
    counts = pd.DataFrame(dense_counts(adata, func="da.maaslin3"), index=adata.obs_names, columns=adata.var_names)
    # Plain names in the formula, so an obs column such as "body site" needs no quoting. rpy2 makes an R factor, whose
    # first level maaslin3 keeps as reference, only from string categories; bool ones would arrive as sorted strings.
    metadata = pd.DataFrame(
        {f"x{i}": _string_levels(frame[name]) for i, name in enumerate(frame.columns)}, index=frame.index
    )
    fit = r_function(_MAASLIN, package="maaslin3", func="da.maaslin3")
    table = call_r(fit, counts, metadata, "~ " + " + ".join(metadata.columns), seed_for_r)
    rows = table[table["metadata"] == "x0"].set_index("feature").reindex(adata.var_names)
    untested = rows["failed"].fillna(True).to_numpy(bool) | ~np.isfinite(rows["pval"].to_numpy(np.float64))
    effect, se, pvalue = (np.where(untested, np.nan, rows[c].to_numpy(np.float64)) for c in ("coef", "stderr", "pval"))
    return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="maaslin3", contrast=contrast)


def _string_levels(values: pd.Series) -> pd.Series:
    """A categorical column with its levels as strings, in the same order; a numeric one as it is."""
    if isinstance(values.dtype, pd.CategoricalDtype):
        return values.cat.rename_categories([str(level) for level in values.cat.categories])
    return values
