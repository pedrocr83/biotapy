"""phyloseq objects saved from R (.RData/.rda/.rds), read natively (decisions/phyloseq-import-route)."""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

from biotapy._core import TreeData, as_csr, infer_x_kind, make_treedata, normalize_ranks, tree_from_phylo

from ._join import _join_to
from ._qiime2 import _TEXT
from ._rdata import load_phyloseq


def read_phyloseq(path: str | Path, *, name: str | None = None) -> TreeData:
    """Read a phyloseq object saved from R, without R.

    Parameters
    ----------
    path
        ``.rds`` file (``saveRDS``) or ``.RData``/``.rda`` file (``save``).
    name
        Object to read from an ``.RData`` file holding several phyloseq objects.

    Returns
    -------
    TreeData
        ``otu_table`` in ``X`` (samples are rows); ``tax_table`` as rank columns in
        ``var``; ``sample_data`` in ``obs``; ``phy_tree`` in ``vart['phylo']``.

    Raises
    ------
    ValueError
        The file holds no phyloseq object, holds several and ``name`` is not given,
        the phyloseq object has no ``otu_table``, or the ``refseq`` slot holds
        sequences: rdata cannot parse a populated refseq yet, so the message gives the
        R fix (``Biostrings::writeXStringSet`` then re-save without the slot).
    KeyError
        ``name`` is not a phyloseq object in the file.

    Notes
    -----
    R equivalent: ``phyloseq::phyloseq``
    Guide: :doc:`/guide/reading_data`

    The OTU table is read densely once (8 bytes x samples x taxa) and stored sparse.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.io.read_phyloseq("ps.rds")  # doctest: +SKIP
    """
    slots = load_phyloseq(Path(path), name=name)
    if slots["refseq"] is not None:
        # A later rdata release might parse Biostrings sequences; until then they must
        # not be dropped silently (R7.4), so a populated refseq is a hard error, not a
        # warned skip (ruling 2026-09-27).
        msg = (
            f"path={str(path)!r}: the refseq slot holds sequences, which rdata cannot parse; "
            'export them first: Biostrings::writeXStringSet(refseq(ps), "refseq.fasta"); '
            'ps@refseq <- NULL; saveRDS(ps, "ps.rds")'
        )
        raise ValueError(msg)
    X, samples, features = _counts(slots["otu_table"], path)
    var = pd.DataFrame(index=features)
    if slots["tax_table"] is not None:
        ranks = normalize_ranks(slots["tax_table"])
        var = _join_to(ranks, pd.Index(features), argument=f"path={str(path)!r} (tax_table)")
    obs = pd.DataFrame(index=samples)
    if slots["sam_data"] is not None:
        sam_data = _text_columns(slots["sam_data"])
        obs = _join_to(sam_data, pd.Index(samples), argument=f"path={str(path)!r} (sam_data)")
    phylo = slots["phy_tree"]
    tree = (
        None
        if phylo is None
        else tree_from_phylo(
            phylo["edge"],
            phylo.get("edge.length"),
            [str(tip) for tip in phylo["tip.label"]],
            argument=f"path={str(path)!r} (phy_tree)",
        )
    )
    return make_treedata(X, obs=obs, var=var, tree=tree, x_kind=infer_x_kind(X), source="io.read_phyloseq")


def _counts(otu: Any, path: str | Path) -> tuple[sp.csr_matrix, list[str], list[str]]:
    if otu is None:
        msg = f"path={str(path)!r}: the phyloseq object has no otu_table"
        raise ValueError(msg)
    array, taxa_are_rows = otu
    rows, columns = ([str(v) for v in array.coords[dim].values] for dim in array.dims)
    values = np.asarray(array.values)
    # Samples are rows (R6.1): transpose once when phyloseq stored taxa as rows.
    return (as_csr(values.T), columns, rows) if taxa_are_rows else (as_csr(values), rows, columns)


def _text_columns(frame: pd.DataFrame) -> pd.DataFrame:
    # dataframe_constructor gives R character columns pd.NA-backed strings; read_qiime2's
    # NaN-backed dtype (._qiime2._TEXT) is what round-trips through h5td (h5py cannot
    # write an all-missing object column, and pd.NA behaves differently there than NaN).
    text = [column for column in frame.columns if pd.api.types.is_string_dtype(frame[column])]
    return frame.astype(dict.fromkeys(text, _TEXT)) if text else frame
