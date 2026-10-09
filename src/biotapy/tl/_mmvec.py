"""mmvec: which metabolites co-occur with which microbes, over samples of both."""

from typing import cast

import numpy as np
import pandas as pd
from mudata import MuData
from skbio.stats.ordination import mmvec as skbio_mmvec

from biotapy._core import as_csr, as_generator


def mmvec(
    mdata: MuData,
    *,
    microbes: str = "taxa",
    metabolites: str = "metabolites",
    seed: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """Learn how likely each metabolite is given each microbe (Morton et al. 2019).

    Parameters
    ----------
    mdata
        Modalities over the same samples, in the same order, such as
        :func:`biotapy.io.to_mudata` returns; ``X`` holds counts or another
        non-negative abundance.
    microbes
        The modality whose features condition the model.
    metabolites
        The modality whose features are predicted.
    seed
        Seed or generator for the starting values of the fit.

    Returns
    -------
    pandas.DataFrame
        Microbes (rows, ``mdata[microbes].var_names``) x metabolites (columns):
        the log conditional probability of each metabolite given each microbe,
        centred so that every row sums to 0. Larger means more likely to
        co-occur.

    Raises
    ------
    KeyError
        ``microbes`` or ``metabolites`` is not a modality.
    ValueError
        The two modalities hold different samples or a different order; a value
        is negative or not finite; or a feature or sample is all zero.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/multiomics`

    Wraps :func:`skbio.stats.ordination.mmvec` with its defaults: three
    latent dimensions and up to 1,000 L-BFGS iterations. For the embeddings,
    the fitted probabilities or a prediction on new samples, call scikit-bio
    with the same tables.

    scikit-bio needs dense tables, so both ``X`` are densified once:
    8 bytes x samples x features of each modality.

    References
    ----------
    Morton JT et al. (2019) Learning representations of microbe-metabolite
    interactions. Nature Methods 16:1306-1314.

    Examples
    --------
    >>> import anndata as ad
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> metabolites = ad.AnnData(tdata.X[:, [5, 2]].toarray() + 1.0, obs=tdata.obs[[]])
    >>> mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
    >>> ranks = bt.tl.mmvec(mdata, seed=0)
    >>> ranks.shape
    (8, 2)
    """
    rng = as_generator(seed)
    x_table = _table(mdata, microbes, argument="microbes")
    y_table = _table(mdata, metabolites, argument="metabolites")
    if not x_table.index.equals(y_table.index):
        msg = (
            f"microbes={microbes!r} and metabolites={metabolites!r} hold different samples or a different order; "
            "align them with bt.io.to_mudata"
        )
        raise ValueError(msg)
    return cast("pd.DataFrame", skbio_mmvec(x_table, y_table, seed=rng).ranks)


def _table(mdata: MuData, key: str, *, argument: str) -> pd.DataFrame:
    """Modality ``key`` as a dense samples x features DataFrame, checked as mmvec needs it."""
    if key not in mdata.mod:
        msg = f"{argument}={key!r} is not a modality; found {list(mdata.mod)}"
        raise KeyError(msg)
    mod = mdata.mod[key]
    X = as_csr(mod.X).astype(np.float64)
    if not np.all(np.isfinite(X.data)) or np.any(X.data < 0):
        msg = f"tl.mmvec needs finite, non-negative values in {argument}={key!r}"
        raise ValueError(msg)
    for axis, names, what in ((0, mod.var_names, "features"), (1, mod.obs_names, "samples")):
        empty = names[np.asarray(X.sum(axis=axis)).ravel() == 0]
        if len(empty):
            msg = f"{argument}={key!r} has all-zero {what}: {empty[:5].tolist()}"
            raise ValueError(msg)
    # scikit-bio's mmvec needs dense input (rules.md R6.2): one dense copy per modality.
    return pd.DataFrame(X.toarray(), index=mod.obs_names, columns=mod.var_names)
