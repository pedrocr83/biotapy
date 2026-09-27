"""Ordination of a stored distance matrix: PCoA and non-metric MDS."""

import numpy as np
import pandas as pd
from anndata import AnnData
from skbio.stats.ordination import pcoa as skbio_pcoa
from sklearn.manifold import MDS

from biotapy._core import as_generator

from ._beta import stored_distances

# vegan::metaMDS's try = 20: the best of 20 random starts.
_NMDS_STARTS = 20


def pcoa(
    adata: AnnData, *, distance: str = "braycurtis", n_components: int = 10, inplace: bool = False
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Principal coordinates analysis of a distance matrix in ``obsp``.

    Parameters
    ----------
    adata
        Samples x features with ``obsp[distance]``, written by :func:`biotapy.tl.beta`
        or :func:`biotapy.tl.unifrac` with ``inplace=True``.
    distance
        The ``obsp`` key to ordinate.
    n_components
        Axes to keep; at most ``n_obs - 1`` are returned.
    inplace
        Write the coordinates to ``obsm['X_pcoa']`` and the axes to
        ``uns['biotapy']['pcoa']`` (``eigenvalues``, ``proportion_explained``), and return ``None``.

    Returns
    -------
    tuple of pandas.DataFrame, or None
        Coordinates (samples x ``PC1``..), and per axis its ``eigenvalue`` and
        ``proportion_explained`` (eigenvalue / sum of all eigenvalues, as ``ape::pcoa``'s ``Relative_eig``).

    Raises
    ------
    KeyError
        ``obsp[distance]`` is missing; the message names the call that writes it.
    ValueError
        ``n_components`` is below 1, there are fewer than 2 samples, or the distances hold NaN.

    Notes
    -----
    R equivalent: ``phyloseq::ordinate``, ``ape::pcoa``
    Guide: :doc:`/guide/ordination`

    Like ``ordinate(physeq, "PCoA", ...)`` (``ape::pcoa`` with ``correction = "none"``).
    Axis signs are arbitrary, as in R. scikit-bio sets negative eigenvalues, and their
    coordinates, to 0, whereas ``ape::pcoa`` reports those eigenvalues as negative.
    Proportions divide by the trace, the sum of all eigenvalues with the negative ones
    included, so they still match ape's ``Relative_eig`` on the positive axes.

    References
    ----------
    Gower JC (1966) Some distance properties of latent root and vector methods used in
    multivariate analysis. Biometrika 53:325-338.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> bt.tl.beta(tdata, inplace=True)
    >>> coords, axes = bt.tl.pcoa(tdata, n_components=2)
    >>> coords.shape
    (6, 2)
    """
    distances = stored_distances(adata, distance)
    if n_components < 1 or adata.n_obs < 2:
        msg = f"tl.pcoa needs n_components >= 1 and at least 2 samples, got {n_components} and {adata.n_obs}"
        raise ValueError(msg)
    # With fewer axes than samples, scikit-bio divides by the trace of the centred matrix, as ape::pcoa does.
    dimensions = min(n_components, adata.n_obs - 1)
    result = skbio_pcoa(distances, dimensions=dimensions)
    names = [f"PC{i}" for i in range(1, dimensions + 1)]
    coords = pd.DataFrame(result.samples.to_numpy(), index=adata.obs_names, columns=names)
    axes = pd.DataFrame(
        {"eigenvalue": result.eigvals.to_numpy(), "proportion_explained": result.proportion_explained.to_numpy()},
        index=names,
    )
    if not inplace:
        return coords, axes
    adata.obsm["X_pcoa"] = coords.to_numpy()
    adata.uns.setdefault("biotapy", {})["pcoa"] = {
        "eigenvalues": axes["eigenvalue"].to_numpy(),
        "proportion_explained": axes["proportion_explained"].to_numpy(),
    }
    return None


def nmds(
    adata: AnnData,
    *,
    distance: str = "braycurtis",
    n_components: int = 2,
    seed: int | np.random.Generator | None = None,
    inplace: bool = False,
) -> tuple[pd.DataFrame, float] | None:
    """Non-metric multidimensional scaling of a distance matrix in ``obsp``.

    Parameters
    ----------
    adata
        Samples x features with ``obsp[distance]``, written by :func:`biotapy.tl.beta`
        or :func:`biotapy.tl.unifrac` with ``inplace=True``.
    distance
        The ``obsp`` key to ordinate.
    n_components
        Dimensions of the configuration.
    seed
        Seed or generator for the random starts.
    inplace
        Write the configuration to ``obsm['X_nmds']`` and the stress to
        ``uns['biotapy']['nmds']['stress']``, and return ``None``.

    Returns
    -------
    tuple of pandas.DataFrame and float, or None
        The configuration (samples x ``NMDS1``..) and its Kruskal stress-1.

    Raises
    ------
    KeyError
        ``obsp[distance]`` is missing; the message names the call that writes it.
    ValueError
        ``n_components`` is below 1, there are not more samples than
        ``n_components + 1``, or the distances hold NaN.

    Notes
    -----
    R equivalent: ``phyloseq::ordinate``, ``vegan::metaMDS``
    Guide: :doc:`/guide/ordination`

    scikit-learn's SMACOF, non-metric, keeps the best of 20 random starts, as
    ``vegan::metaMDS`` does by default (``try = 20``). The configuration is not
    rotated or scaled afterwards, unlike metaMDS's ``postMDS``, so compare
    configurations up to rotation and scale (Procrustes).

    References
    ----------
    Kruskal JB (1964) Nonmetric multidimensional scaling: a numerical method. Psychometrika
    29:115-129.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> bt.tl.beta(tdata, inplace=True)
    >>> coords, stress = bt.tl.nmds(tdata, seed=0)
    >>> coords.shape
    (6, 2)
    """
    distances = stored_distances(adata, distance)
    if n_components < 1 or adata.n_obs <= n_components + 1:
        msg = f"tl.nmds needs n_components >= 1 and more than n_components + 1 samples, got {n_components} and {adata.n_obs}"
        raise ValueError(msg)
    # scikit-learn takes an int or a RandomState, not a Generator: draw its seed from ours.
    random_state = int(as_generator(seed).integers(2**31 - 1))
    mds = MDS(
        n_components=n_components,
        metric_mds=False,
        metric="precomputed",
        n_init=_NMDS_STARTS,
        init="random",
        normalized_stress="auto",
        random_state=random_state,
    )
    names = [f"NMDS{i}" for i in range(1, n_components + 1)]
    coords = pd.DataFrame(mds.fit_transform(distances.data), index=adata.obs_names, columns=names)
    stress = float(mds.stress_)
    if not inplace:
        return coords, stress
    adata.obsm["X_nmds"] = coords.to_numpy()
    adata.uns.setdefault("biotapy", {})["nmds"] = {"stress": stress}
    return None
