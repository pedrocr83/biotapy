"""Alpha diversity: one value per sample and metric."""

import math
from collections.abc import Sequence
from typing import Any, Literal, get_args

import numpy as np
import pandas as pd
from anndata import AnnData
from skbio.diversity import alpha_diversity

from biotapy._core import as_csr, get_skbio_tree, require_counts

AlphaMetric = Literal["observed_features", "shannon", "simpson", "chao1", "faith_pd"]
# phyloseq's estimate_richness: natural-log Shannon (vegan::diversity) and bias-corrected Chao1 (vegan::estimateR).
_PARAMS: dict[str, dict[str, Any]] = {"shannon": {"base": math.e}, "chao1": {"bias_corrected": True}}
# scikit-bio needs dense rows (rules.md R6.2); densify at most 2**20 values (8 MiB of float64) at a time.
_CHUNK_VALUES = 2**20


def alpha(
    adata: AnnData,
    *,
    metrics: Sequence[AlphaMetric] = ("observed_features", "shannon", "simpson", "chao1"),
    inplace: bool = False,
) -> pd.DataFrame | None:
    """Alpha diversity of every sample.

    Parameters
    ----------
    adata
        Samples x features. ``"faith_pd"`` needs a TreeData with ``vart['phylo']``.
    metrics
        Any of ``"observed_features"``, ``"shannon"`` (natural log), ``"simpson"``
        (Gini-Simpson, ``1 - sum(p**2)``), ``"chao1"`` (bias-corrected) and ``"faith_pd"``.
    inplace
        Write ``obs['alpha_<metric>']`` for every metric and return ``None``.

    Returns
    -------
    pandas.DataFrame or None
        One row per sample (index ``obs_names``) and one column per metric, in the
        order given. An all-zero sample gets 0 for ``observed_features``, ``chao1``
        and ``faith_pd`` and NaN for ``shannon`` and ``simpson``.

    Raises
    ------
    TypeError
        ``metrics`` is a string, or ``"faith_pd"`` is asked of an AnnData that is not a TreeData.
    ValueError
        ``metrics`` is empty, names an unknown metric or repeats one, or
        ``"observed_features"`` or ``"chao1"`` is asked of data that is not raw
        counts (``x_kind`` ``"counts"`` and whole numbers).
    KeyError
        ``"faith_pd"`` is asked of a TreeData without ``vart['phylo']``.

    Notes
    -----
    R equivalent: ``phyloseq::estimate_richness``, ``picante::pd``
    Guide: :doc:`/guide/diversity`

    scikit-bio needs dense input, so rows are densified in chunks of at most 2**20
    values (8 MiB of float64). Faith PD includes the root, as
    ``picante::pd(include.root = TRUE)``; a root with more than two children first
    gets a zero-length split, which changes no root-to-tip distance. Faith PD
    depends on presence only, so it is computed on presence/absence and runs on
    any abundance: scikit-bio's tree code casts abundances to integers, which
    truncates proportions to 0. For an all-zero sample phyloseq reports Shannon 0
    and Simpson 1; biotapy returns NaN.

    References
    ----------
    Shannon CE (1948) A mathematical theory of communication. Bell System Technical Journal
    27:379-423.

    Simpson EH (1949) Measurement of diversity. Nature 163:688.

    Chao A (1984) Nonparametric estimation of the number of classes in a population.
    Scandinavian Journal of Statistics 11:265-270.

    Faith DP (1992) Conservation evaluation and phylogenetic diversity. Biological
    Conservation 61:1-10.

    Examples
    --------
    >>> import biotapy as bt
    >>> out = bt.tl.alpha(bt.datasets.toy(), metrics=["observed_features", "shannon"])
    >>> out.loc["s1"].round(3).tolist()
    [6.0, 1.361]
    """
    _check_metrics(adata, metrics)
    params = {metric: _PARAMS.get(metric, {}) for metric in metrics}
    if "faith_pd" in metrics:
        params["faith_pd"] = {"taxa": adata.var_names.tolist(), "tree": get_skbio_tree(adata)}
    X = as_csr(adata.X)
    step = max(1, _CHUNK_VALUES // max(adata.n_vars, 1))
    frames = []
    for start in range(0, adata.n_obs, step):
        dense, ids = X[start : start + step].toarray(), adata.obs_names[start : start + step].tolist()
        # Faith PD depends on presence only, and scikit-bio's tree code casts abundances to int64,
        # which would turn every proportion into 0.
        inputs = {m: (dense > 0).astype(np.int64) if m == "faith_pd" else dense for m in metrics}
        frames.append(pd.DataFrame({m: alpha_diversity(m, inputs[m], ids=ids, **params[m]) for m in metrics}))
    result = pd.concat(frames)
    if not inplace:
        return result
    for metric in metrics:
        adata.obs[f"alpha_{metric}"] = result[metric].to_numpy()
    return None


def _check_metrics(adata: AnnData, metrics: Sequence[str]) -> None:
    if isinstance(metrics, str):
        msg = f"metrics must be a list of metric names, not the string {metrics!r}"
        raise TypeError(msg)
    known = get_args(AlphaMetric)
    if not metrics or any(metric not in known for metric in metrics):
        msg = f"metrics must name one or more of {list(known)}, got {list(metrics)}"
        raise ValueError(msg)
    repeated = [metric for metric in known if list(metrics).count(metric) > 1]
    if repeated:
        msg = f"metrics repeats {repeated}; name each metric once"
        raise ValueError(msg)
    needs_counts = [metric for metric in metrics if metric in ("observed_features", "chao1")]
    if needs_counts:
        require_counts(adata, func=f"tl.alpha with {needs_counts}")
