"""The taxa behind one function: its stratified rows as a samples x taxa table."""

import difflib
from numbers import Integral
from typing import cast

import numpy as np
import pandas as pd
from anndata import AnnData

from biotapy._core import BY_TAXON_KEY, as_csr

OTHER = "other"
_STRATIFIED_COLUMNS = ("function", "taxon")


def contributions(adata: AnnData, function: str, *, top: int | None = None) -> pd.DataFrame:
    """Per-taxon abundance of one function in every sample.

    Parameters
    ----------
    adata
        A stratified function table: the ``"function_by_taxon"`` modality of
        ``bt.io.read_humann`` or ``bt.io.read_picrust2``, or
        ``bt.fn.func_glom``'s or ``bt.fn.renorm``'s output for it (``var``
        columns ``function`` and ``taxon``).
    function
        A function id as it appears in ``var["function"]``, for example
        ``"PWY-5100"``, ``"2.7.1.1"`` or ``"UNINTEGRATED"``.
    top
        Keep the ``top`` taxa with the largest total over all samples and sum
        the rest into one column, ``"other"``. By default every taxon is a
        column.

    Returns
    -------
    pandas.DataFrame
        Samples x taxa, ``float64``, indexed by ``obs_names``; the column
        index is named ``"taxon"``. Values are the stratified rows as stored
        in ``X``, so each row sums to the function's strata in that sample.
        Columns are ordered by total over all samples, largest first (ties by
        taxon id), and ``"other"`` comes last. ``unclassified`` and ``RARE``
        are ordinary taxa.

    Raises
    ------
    TypeError
        ``adata`` is not an AnnData (for example the whole MuData).
    KeyError
        ``adata`` lacks the ``function`` or ``taxon`` column (for example the
        ``"function"`` modality); ``function`` has no stratified row in
        ``adata``. The message lists up to three close function ids.
    ValueError
        ``top`` is not a positive integer or ``None``; ``top`` would sum
        taxa into ``"other"`` while a taxon is itself named ``"other"``.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/function`

    For HUMAnN and PICRUSt2 pathways the strata need not sum to the
    community value, so the rows are not shares of it. To read them as
    shares of each sample's community total, pass the stratified modality
    of ``bt.fn.renorm(mdata, "relab")``. Only the function's columns of
    ``X`` are densified.

    Examples
    --------
    >>> import biotapy as bt
    >>> by_taxon = bt.datasets.toy_humann()["function_by_taxon"]
    >>> table = bt.fn.contributions(by_taxon, "2.7.1.2")
    >>> table.columns.tolist()
    ['g__Blautia.s__Blautia_obeum', 'g__Bacteroides.s__Bacteroides_ovatus']
    >>> table.loc["s4"].tolist()
    [7.0, 1.0]
    """
    var = _stratified_var(adata)
    if top is not None and (isinstance(top, bool) or not isinstance(top, Integral) or top < 1):
        msg = f"top={top!r} must be a positive integer or None"
        raise ValueError(msg)
    columns = np.flatnonzero((var["function"] == function).to_numpy())
    if columns.size == 0:
        close = difflib.get_close_matches(function, var["function"].unique().tolist(), n=3)
        msg = f"function={function!r} has no rows in adata; close function ids: {close}"
        raise KeyError(msg)
    # One function's columns only: samples x its taxa, small next to X.
    values = as_csr(adata.X)[:, columns].toarray().astype(np.float64)
    taxa = pd.Index(var["taxon"].to_numpy()[columns], name="taxon")
    table = pd.DataFrame(values, index=adata.obs_names.copy(), columns=taxa)
    table = table.iloc[:, np.lexsort((taxa.to_numpy(), -values.sum(axis=0)))]
    if top is None or top >= table.shape[1]:
        return table
    if OTHER in table.columns:
        msg = f"top={top} would sum taxa into a column named {OTHER!r}, but adata has a taxon named {OTHER!r}; pass top=None"
        raise ValueError(msg)
    kept = table.iloc[:, :top].copy()
    kept[OTHER] = table.iloc[:, top:].sum(axis=1)
    return kept


def _stratified_var(adata: AnnData) -> pd.DataFrame:
    """``adata.var``, after checking it is a stratified function table."""
    if not isinstance(adata, AnnData):
        msg = f"adata must be an AnnData such as mdata[{BY_TAXON_KEY!r}], not {type(adata).__name__}"
        raise TypeError(msg)
    var = cast("pd.DataFrame", adata.var)
    if not all(column in var.columns for column in _STRATIFIED_COLUMNS):
        msg = (
            f"adata needs var columns {list(_STRATIFIED_COLUMNS)}, as the {BY_TAXON_KEY!r} modality has; "
            f"pass mdata[{BY_TAXON_KEY!r}]"
        )
        raise KeyError(msg)
    return var
