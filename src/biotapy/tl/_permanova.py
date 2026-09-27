"""PERMANOVA: do groups of samples differ in a stored distance matrix?"""

import numpy as np
import pandas as pd
from anndata import AnnData
from skbio.stats.distance import permanova as skbio_permanova

from biotapy._core import as_generator

from ._beta import stored_distances


def permanova(
    adata: AnnData,
    grouping: str,
    *,
    distance: str = "braycurtis",
    permutations: int = 999,
    seed: int | np.random.Generator | None = None,
) -> pd.Series:
    """Permutational multivariate analysis of variance of a distance matrix in ``obsp``.

    Parameters
    ----------
    adata
        Samples x features with ``obsp[distance]``, written by :func:`biotapy.tl.beta`
        or :func:`biotapy.tl.unifrac` with ``inplace=True``.
    grouping
        The ``obs`` column holding each sample's group: categorical, string or bool.
    distance
        The ``obsp`` key to test.
    permutations
        Permutations for the p-value.
    seed
        Seed or generator for the permutations.

    Returns
    -------
    pandas.Series
        scikit-bio's result: ``test statistic`` (pseudo-F), ``p-value``,
        ``sample size``, ``number of groups``, ``number of permutations`` and the
        method and statistic names.

    Raises
    ------
    KeyError
        ``grouping`` is not an ``obs`` column, or ``obsp[distance]`` is missing.
    TypeError
        ``obs[grouping]`` is numeric (and not bool); convert it with ``.astype("category")`` for groups.
    ValueError
        ``grouping`` has missing values, or the distances hold NaN.

    Notes
    -----
    R equivalent: ``vegan::adonis2``
    Guide: :doc:`/guide/ordination`

    Matches ``adonis2(distance ~ grouping, data, permutations)`` with one term,
    whose ``F`` is the test statistic. P-values agree only up to permutation noise:
    R and NumPy random generators differ. ``grouping`` is categorical, one group
    per distinct value, like an R factor. adonis2 fits a numeric column as one
    continuous term instead, so a numeric column raises rather than silently
    becoming one group per value.

    References
    ----------
    Anderson MJ (2001) A new method for non-parametric multivariate analysis of variance.
    Austral Ecology 26:32-46.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> bt.tl.beta(tdata, inplace=True)
    >>> result = bt.tl.permanova(tdata, "group", seed=0)
    >>> int(result["number of groups"])
    2
    """
    if grouping not in adata.obs.columns:
        msg = f"grouping={grouping!r} is not a column of obs"
        raise KeyError(msg)
    groups = adata.obs[grouping]
    if pd.api.types.is_numeric_dtype(groups) and not pd.api.types.is_bool_dtype(groups):
        msg = (
            f"grouping={grouping!r} is a numeric column ({groups.dtype}); tl.permanova compares groups, "
            f'so convert it with .astype("category") for one group per value'
        )
        raise TypeError(msg)
    if groups.isna().any():
        msg = f"grouping={grouping!r} is missing for {int(groups.isna().sum())} sample(s); drop them first"
        raise ValueError(msg)
    distances = stored_distances(adata, distance)
    result: pd.Series = skbio_permanova(
        distances, groups.to_numpy(), permutations=permutations, seed=as_generator(seed)
    )
    return result
