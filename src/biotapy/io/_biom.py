"""BIOM tables (JSON 1.0 and HDF5 2.1) through biom-format."""

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp

from biotapy._core import TreeData, as_csr, make_treedata, split_lineage, tree_from_newick


def read_biom(path: str | Path, *, tree: str | Path | None = None) -> TreeData:
    """Read a BIOM table (JSON 1.0 or HDF5 2.1) with samples as rows.

    Parameters
    ----------
    path
        BIOM file; JSON or HDF5 is detected from the content.
    tree
        Newick file whose tips are the table's observation ids.

    Returns
    -------
    TreeData
        Counts in ``X``; observation ``taxonomy`` metadata as rank columns in
        ``var``; sample metadata in ``obs``; the tree in ``vart['phylo']``.

    Raises
    ------
    biom.exception.TableException
        The file repeats a sample or observation id (raised by biom-format).

    Warns
    -----
    UserWarning
        Tree tips and table features differ; only shared features are kept.

    Notes
    -----
    R equivalent: ``phyloseq::import_biom``
    Guide: :doc:`/guide/reading_data`

    ``X`` is read as counts; BIOM does not record whether it holds counts.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biom
    >>> from biom.util import biom_open
    >>> import biotapy as bt
    >>> toy = bt.datasets.toy()
    >>> path = Path(tempfile.mkdtemp()) / "toy.biom"
    >>> table = biom.Table(toy.X.T, list(toy.var_names), list(toy.obs_names))
    >>> with biom_open(str(path), "w") as handle:
    ...     table.to_hdf5(handle, "example")
    >>> bt.io.read_biom(path).shape
    (6, 8)
    """
    X, obs, var = _biom_parts(biom.load_table(str(path)))
    phylo = None if tree is None else tree_from_newick(Path(tree).read_text())
    return make_treedata(X, obs=obs, var=var, tree=phylo, x_kind="counts", source="io.read_biom")


def _biom_parts(table: biom.Table) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]:
    """Counts as samples x features, plus obs and var frames, from a biom Table."""
    samples = [str(i) for i in table.ids(axis="sample")]
    features = [str(i) for i in table.ids(axis="observation")]
    # BIOM stores features x samples; transpose once so samples are rows (R6.1).
    X = as_csr(table.matrix_data.T)
    obs = _metadata_frame(table.metadata(axis="sample"), samples)
    var = _taxonomy_frame(table.metadata(axis="observation"), features)
    return X, obs, var


def _metadata_frame(metadata: Sequence[Mapping[str, object]] | None, ids: list[str]) -> pd.DataFrame:
    if metadata is None:
        return pd.DataFrame(index=ids)
    return pd.DataFrame([dict(entry) for entry in metadata], index=ids)


def _taxonomy_frame(metadata: Sequence[Mapping[str, object]] | None, ids: list[str]) -> pd.DataFrame:
    if metadata is None:
        return pd.DataFrame(index=ids)
    return split_lineage(pd.Series([_lineage(entry) for entry in metadata], index=ids, dtype=object))


def _lineage(entry: Mapping[str, object]) -> str | float:
    # HDF5 gives a list of ranks; JSON gives whatever the producer wrote (list or string).
    value = entry.get("taxonomy") or entry.get("Taxonomy")
    if isinstance(value, str):
        return value
    if isinstance(value, Iterable):
        return "; ".join(str(part) for part in value)
    return np.nan
