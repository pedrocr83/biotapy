"""BIOM tables (JSON 1.0 and HDF5 2.1) through biom-format."""

from collections.abc import Iterable, Mapping, Sequence
from importlib.metadata import version
from pathlib import Path
from typing import Literal, cast

import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from biom.util import biom_open

from biotapy._core import RANKS, TreeData, as_csr, infer_x_kind, make_treedata, split_lineage, tree_from_newick


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
        The table in ``X``; observation ``taxonomy`` metadata as rank columns in
        ``var``; sample metadata in ``obs``, with empty values as NaN; the
        tree in ``vart['phylo']``.

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

    BIOM does not record what ``X`` holds, so ``uns['biotapy']['x_kind']`` is
    inferred from the values: whole numbers are ``"counts"``, rows that each
    sum to 1 are ``"relative"``, anything else is ``"abundance"``.

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
    return make_treedata(X, obs=obs, var=var, tree=phylo, x_kind=infer_x_kind(X), source="io.read_biom")


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
    # BIOM has no null: write_biom (and others) write a missing value as "".
    return pd.DataFrame([dict(entry) for entry in metadata], index=ids).replace("", np.nan)


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


def write_biom(adata: AnnData, path: str | Path, *, fmt: Literal["hdf5", "json"] = "hdf5") -> None:
    """Write ``X`` with taxonomy and sample metadata as a BIOM table.

    Parameters
    ----------
    adata
        Samples x features; rank columns in ``var`` become observation
        ``taxonomy`` metadata, ``obs`` columns become sample metadata.
    path
        Output file.
    fmt
        ``"hdf5"`` (BIOM 2.1) or ``"json"`` (BIOM 1.0).

    Notes
    -----
    R equivalent: ``biomformat::write_biom``
    Guide: :doc:`/guide/reading_data`

    BIOM has no slot for a tree, layers or embeddings: a TreeData's tree and
    everything outside ``X``, rank columns and ``obs`` are not written.
    Sample metadata is written as text; missing values are written as empty
    strings and read back as NaN. Missing ranks are written as bare prefixes
    (``g__``) so every rank keeps its place.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> path = Path(tempfile.mkdtemp()) / "toy.biom"
    >>> bt.io.write_biom(bt.datasets.toy(), path)
    >>> bt.io.read_biom(path).shape
    (6, 8)
    """
    table = biom.Table(
        # biotapy keeps samples as rows; BIOM stores features x samples.
        as_csr(adata.X).T,
        [str(i) for i in adata.var_names],
        [str(i) for i in adata.obs_names],
        # anndata types .var/.obs as DataFrame | Dataset2D (its lazy/backed variant); biotapy's
        # data model (data-model-slots) guarantees a real DataFrame here.
        observation_metadata=_taxonomy_metadata(cast("pd.DataFrame", adata.var)),
        sample_metadata=_sample_metadata(cast("pd.DataFrame", adata.obs)),
    )
    generated_by = f"biotapy {version('biotapy')}"
    if fmt == "json":
        Path(path).write_text(table.to_json(generated_by))
        return
    with biom_open(str(path), "w") as handle:
        table.to_hdf5(handle, generated_by)


def _taxonomy_metadata(var: pd.DataFrame) -> list[dict[str, list[str]]] | None:
    ranks = [rank for rank in RANKS if rank in var.columns]
    if not ranks:
        return None
    # Prefixed values keep their rank: BIOM HDF5 drops empty list entries on read.
    values = var[ranks].astype("string").fillna("")
    return [
        {"taxonomy": [f"{rank[0]}__{value}" for rank, value in zip(ranks, row, strict=True)]}
        for row in values.itertuples(index=False)
    ]


def _sample_metadata(obs: pd.DataFrame) -> list[dict[str, str]] | None:
    if obs.columns.empty:
        return None
    values = obs.astype("string").fillna("")
    return [dict(zip(map(str, values.columns), row, strict=True)) for row in values.itertuples(index=False)]
