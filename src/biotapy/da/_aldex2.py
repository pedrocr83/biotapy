"""ALDEx2 through the R bridge: Monte Carlo CLR values of two groups compared by Welch's t-test."""

import numpy as np
import pandas as pd
from anndata import AnnData

from ._design import dense_counts, model
from ._r import call_r, r_function, r_seed
from ._schema import result

# aldex() is ALDEx2's documented entry point and the golden file's call; its progress notes are R messages, which
# rpy2 would log as warnings.
_ALDEX = """function(reads, conditions, mc_samples, seed) {
  set.seed(seed)
  suppressMessages(ALDEx2::aldex(reads, conditions, mc.samples = mc_samples, test = "t", effect = TRUE, denom = "all"))
}"""


def aldex2(
    adata: AnnData,
    group: str,
    *,
    mc_samples: int = 128,
    reference: str | None = None,
    seed: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """Differential abundance of each feature between two groups by ALDEx2, run in R.

    Parameters
    ----------
    adata
        Samples x features; ``X`` holds raw counts.
    group
        The ``obs`` column whose two levels are compared: a categorical, string or
        bool column with two levels, each in at least two samples.
    mc_samples
        Monte Carlo draws from each sample's Dirichlet posterior (``mc.samples``).
    reference
        The level of ``group`` that the other level is compared with; by default
        its first category (sorted values for a string column). It changes the
        Monte Carlo draws as well as the sign (Notes).
    seed
        Seeds R's random number generator through ``set.seed`` for the call; the same seed
        gives the same table.

    Returns
    -------
    pandas.DataFrame
        One row per feature, in ``var_names`` order, indexed by ``feature``:
        ``effect`` (ALDEx2's ``diff.btw``, the median log2 difference between the
        groups), ``se`` (NaN: ALDEx2 has none), ``pvalue`` (``we.ep``),
        ``qvalue`` (Benjamini-Hochberg over the tested features), ``direction``,
        ``method`` (``"aldex2"``) and ``contrast``. A feature with no read in any
        sample is not tested: NaN values and ``direction`` 0.

    Raises
    ------
    ImportError
        rpy2 (the extra ``r``) or the R package ALDEx2 is not installed.
    KeyError
        ``group`` is not an ``obs`` column.
    TypeError
        ``reference`` is not a string, or ``mc_samples`` not an int.
    RuntimeError
        R stops with an error.
    ValueError
        ``X`` does not hold raw counts, has an empty sample, repeated
        ``var_names`` or ``obs_names``, or fewer than two features; ``group`` is numeric, has missing values or other than two
        levels, or a level in fewer than two samples; ``reference`` is not one of
        them; ``mc_samples`` is below 1.

    Notes
    -----
    R equivalent: ``ALDEx2::aldex``
    Guide: :doc:`/guide/differential_abundance`

    Calls ``ALDEx2::aldex(reads, conditions, mc.samples, test = "t", effect =
    TRUE, denom = "all")`` (Fernandes et al. 2014) through rpy2: Dirichlet Monte
    Carlo draws of each sample's proportions (0.5 added to every count), their
    log2 centred log-ratios, and per draw a Welch t-test. ``we.ep`` doubles the
    one-sided p-value of each draw in each direction (capped at 1), averages each
    direction over the draws and keeps the smaller average. ``effect`` is
    ``diff.btw``, already log2; ALDEx2's own ``effect`` column is a standardised size, not a fold
    change, and is not carried. ``qvalue`` is the Benjamini-Hochberg correction
    of ``we.ep``, as for every method; ALDEx2's ``we.eBH`` averages the
    corrected values of the draws instead and is not carried. ALDEx2 compares
    two groups without covariates (its ``glm`` test is not wrapped).

    Swapping ``reference`` does more than flip the sign. ALDEx2 draws the same
    Monte Carlo instances, so ``we.ep`` and the calls do not change, but
    ``diff.btw`` comes from a random resampling of each group's values done in
    label order, so the same ``seed`` gives different effects, and the two runs
    are not exact mirror images: they differ from exact antisymmetry by 0.1 to
    0.3 log2 on ``toy()``, depending on the seed (its effects are about 4 log2
    wide), and by up to 0.65 log2 on the GlobalPatterns genera, 15 of 636 of which
    then do not change direction. The p-values and the calls are the same. Each R
    warning raised during the call, such as the one for fewer than 128
    ``mc_samples``, is re-emitted as a ``UserWarning``, and an R error is raised
    as a ``RuntimeError``.

    Needs R with ALDEx2 (``BiocManager::install("ALDEx2")``) and ``pip install
    'biotapy[r]'``, which builds rpy2 (GPL-2.0-or-later) against that R.
    ``seed`` becomes one integer for R's ``set.seed``; the call puts the
    embedded R session's random state back as it found it.

    rpy2 converts dense tables only, so ``X`` is densified once: 8 bytes x
    samples x features, plus R's copy and its Monte Carlo draws (about
    ``mc_samples`` times the table).

    References
    ----------
    Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB (2014)
    Unifying the analysis of high-throughput sequencing datasets: characterizing
    RNA-seq, 16S rRNA gene sequencing and selective growth experiments by
    compositional data analysis. Microbiome 2:15.

    Fernandes AD, Macklaim JM, Linn TG, Reid G, Gloor GB (2013) ANOVA-like
    differential gene expression analysis of single-organism and meta-RNA-seq.
    PLoS ONE 8:e67019.

    Examples
    --------
    >>> import biotapy as bt
    >>> table = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)  # doctest: +SKIP
    >>> table.loc["f6", "contrast"]  # doctest: +SKIP
    'B vs A'
    """
    frame, contrast = model(adata, group, covariates=(), reference=reference, func="da.aldex2")
    conditions = _conditions(frame[group], func="da.aldex2")
    if isinstance(mc_samples, bool) or not isinstance(mc_samples, int | np.integer):
        msg = f"da.aldex2: mc_samples must be an int, got {type(mc_samples).__name__}"
        raise TypeError(msg)
    if mc_samples < 1:
        msg = f"da.aldex2: mc_samples must be at least 1, got {mc_samples}"
        raise ValueError(msg)
    seed_for_r = r_seed(seed)
    # rpy2 has no sparse converter (rules.md R6.2): the one dense copy of X, features as rows as ALDEx2 wants them.
    reads = pd.DataFrame(dense_counts(adata, func="da.aldex2").T, index=adata.var_names, columns=adata.obs_names)
    aldex = r_function(_ALDEX, package="ALDEx2", func="da.aldex2")
    table = call_r(aldex, reads, conditions, int(mc_samples), seed_for_r).reindex(adata.var_names)
    effect = table["diff.btw"].to_numpy(np.float64)
    pvalue = table["we.ep"].to_numpy(np.float64)
    se = np.full(adata.n_vars, np.nan)
    return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="aldex2", contrast=contrast)


def _conditions(values: pd.Series, *, func: str) -> pd.Series:
    """The group as "0" (reference) and "1" labels; raises for a numeric group or a level in fewer than two samples."""
    if not isinstance(values.dtype, pd.CategoricalDtype):
        msg = f"{func} compares two levels, but obs[{values.name!r}] is numeric; ALDEx2's t-test has no slope"
        raise ValueError(msg)
    sizes = values.value_counts()
    if (sizes < 2).any():
        msg = f"{func} needs two samples in each level of obs[{values.name!r}], which has {sizes.to_dict()}"
        raise ValueError(msg)
    # ALDEx2 sorts the labels and reports the second minus the first; "0" < "1" in every R locale, unlike level names.
    return values.cat.codes.astype(str)
