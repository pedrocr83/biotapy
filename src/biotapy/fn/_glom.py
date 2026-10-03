"""Aggregation along a function hierarchy, with humann_regroup_table's semantics."""

from typing import Literal, cast

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from biotapy._core import (
    PROTECTED_FEATURES,
    RANKS,
    SPECIAL_FEATURES,
    UNGROUPED,
    add_provenance,
    as_csr,
    replace_features,
    sum_pairs,
)

_EDGE_COLUMNS = ("child", "parent", "level")
# var columns that describe a stratum's taxon; they are the same for every member of a group.
_TAXON_COLUMNS = ("taxon", *RANKS)


def func_glom(adata: AnnData, level: str, *, hierarchy: pd.DataFrame, agg: Literal["sum", "mean"] = "sum") -> AnnData:
    """Aggregate function features to one level of a hierarchy.

    Each feature's abundance counts in full toward every parent it has at
    ``level`` (one KO in two pathways counts in both), as
    ``humann_regroup_table`` does.

    Parameters
    ----------
    adata
        Samples x functions: one modality of ``bt.io.read_humann``'s result,
        or any AnnData whose ``var_names`` are function ids. A
        ``"function_by_taxon"`` modality (``var`` columns ``function`` and
        ``taxon``) is grouped per taxon.
    level
        The ``level`` value of the hierarchy rows to use, for example
        ``"pathway"`` or ``"class"``.
    hierarchy
        Edge table with columns ``child``, ``parent`` and ``level`` (and
        optionally ``parent_name``), as ``bt.fn.load_hierarchy`` and
        ``bt.datasets.enzyme`` return.
    agg
        ``"sum"``, or ``"mean"`` over the members present in ``adata``.

    Returns
    -------
    AnnData
        Samples x groups, sorted by name. Features with no parent at
        ``level`` are summed into ``UNGROUPED`` (per taxon when stratified);
        ``UNMAPPED``, ``READS_UNMAPPED`` and ``UNINTEGRATED`` pass through
        unchanged. ``var`` holds ``name`` (from ``parent_name``) and
        ``special``; a stratified input also keeps ``function``, ``taxon``
        and its rank columns. ``obs`` and ``uns['biotapy']['x_kind']`` are
        kept; ``layers``, ``obsm``, ``obsp``, ``varm`` and ``varp`` are
        dropped because they described the old features.

    Raises
    ------
    KeyError
        ``hierarchy`` lacks a required column, or has no row at ``level``.
    ValueError
        ``agg`` is not ``"sum"`` or ``"mean"``; no feature of ``adata``
        is a child at ``level``.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/function`

    Matches ``humann_regroup_table`` (HUMAnN 3.9) with its defaults
    ``--ungrouped Y --protected Y``; the golden tests compare the two. One
    difference: ``READS_UNMAPPED`` passes through as in HUMAnN master, where
    3.9 sums it into ``UNGROUPED``. As
    there, a mean divides by the number of members present in the table,
    not by the group's size in ``hierarchy``.

    References
    ----------
    Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of
    diverse microbial communities with bioBakery 3. eLife 10:e65088.

    Examples
    --------
    >>> import pandas as pd
    >>> import biotapy as bt
    >>> edges = pd.DataFrame(
    ...     {"child": ["2.7.1.1", "2.7.1.2"], "parent": ["2.7.1.-", "2.7.1.-"], "level": "subsubclass"}
    ... )
    >>> out = bt.fn.func_glom(bt.datasets.toy_humann()["function"], "subsubclass", hierarchy=edges)
    >>> out.var_names.tolist()
    ['2.7.1.-', 'UNGROUPED', 'UNMAPPED']
    """
    if agg not in ("sum", "mean"):
        msg = f"agg={agg!r} must be 'sum' or 'mean'"
        raise ValueError(msg)
    edges = _edges_at(hierarchy, level)
    var = cast("pd.DataFrame", adata.var)
    function = var["function"] if "function" in var.columns else var.index.to_series()
    pairs = _pairs(function.astype(str).to_numpy(), edges, level=level)
    key = pairs["group"]
    if "taxon" in var.columns:
        key = key + "|" + var["taxon"].astype(str).to_numpy()[pairs["feature"]]
    codes, labels = pd.factorize(key, sort=True)
    X = sum_pairs(as_csr(adata.X), pairs["feature"].to_numpy(dtype=np.intp), codes, n_groups=labels.size)
    if agg == "mean":
        X = sp.csr_matrix(X @ sp.diags(1.0 / np.bincount(codes, minlength=labels.size)))
    out = replace_features(adata, X, _group_var(var, pairs.assign(code=codes), labels, edges=edges))
    add_provenance(out, "fn.func_glom", level=level, agg=agg, hierarchy=hierarchy.attrs.get("source"))
    return out


def _edges_at(hierarchy: pd.DataFrame, level: str) -> pd.DataFrame:
    missing = [column for column in _EDGE_COLUMNS if column not in hierarchy.columns]
    if missing:
        msg = f"hierarchy needs columns {list(_EDGE_COLUMNS)}; missing {missing}"
        raise KeyError(msg)
    if hierarchy[["child", "parent"]].isna().any().any():
        msg = "hierarchy has a missing value in column 'child' or 'parent'"
        raise ValueError(msg)
    edges = hierarchy[hierarchy["level"] == level]
    if edges.empty:
        msg = f"hierarchy has no row at level={level!r}; its levels: {sorted(hierarchy['level'].unique().tolist())}"
        raise KeyError(msg)
    return edges.drop_duplicates(["child", "parent"])


def _pairs(function: np.ndarray, edges: pd.DataFrame, *, level: str) -> pd.DataFrame:
    """(feature, group) rows: each parent of a feature, else itself if protected, else UNGROUPED."""
    features = pd.DataFrame({"feature": np.arange(function.size), "child": function})
    protected = features["child"].isin(PROTECTED_FEATURES)
    mapped = features[~protected].merge(edges[["child", "parent"]], on="child")
    plain = ~features["child"].isin(SPECIAL_FEATURES)
    if plain.any() and not mapped["feature"].isin(features["feature"][plain]).any():
        msg = (
            f"no feature of adata is a child at level={level!r}; "
            f"features: {features['child'][plain].unique()[:3].tolist()}, children: {edges['child'].unique()[:3].tolist()}"
        )
        raise ValueError(msg)
    alone = features[protected | ~features["feature"].isin(mapped["feature"])]
    alone = alone.assign(parent=alone["child"].where(alone["child"].isin(PROTECTED_FEATURES), UNGROUPED))
    rows = pd.concat([mapped, alone], ignore_index=True)
    return rows[["feature", "parent"]].rename(columns={"parent": "group"})


def _group_var(var: pd.DataFrame, pairs: pd.DataFrame, labels: pd.Index, *, edges: pd.DataFrame) -> pd.DataFrame:
    """One var row per group (``pairs["code"]``): function, name, special, and a stratified input's taxon columns."""
    group = pairs.groupby("code")["group"].first()
    first = pairs.groupby("code")["feature"].min().to_numpy()
    names = edges.drop_duplicates("parent").set_index("parent").get("parent_name")
    # The pandas str dtype, never object: anndata's writer rejects an all-NaN object column (contracts/data-model-slots).
    text = pd.StringDtype(na_value=np.nan)
    out = pd.DataFrame(index=labels)
    if "taxon" in var.columns:
        out["function"] = pd.array(group.to_numpy(), dtype=text)
    out["name"] = pd.array(group.map(names).to_numpy() if names is not None else [np.nan] * len(group), dtype=text)
    if "taxon" in var.columns:
        for column in [column for column in _TAXON_COLUMNS if column in var.columns]:
            out[column] = var[column].array[first]
    out["special"] = group.isin(SPECIAL_FEATURES).to_numpy()
    return out
