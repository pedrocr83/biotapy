"""MuData files that keep a TreeData modality's trees (decisions/multiomics-as-mudata)."""

import os
from collections.abc import MutableMapping
from typing import cast

import h5py
import mudata
from anndata import AnnData
from mudata import MuData

from biotapy._core import TREE_SLOTS_ATTR, TreeData, read_tree_slots, write_tree_slots


def write_h5mu(mdata: MuData, path: str | os.PathLike[str]) -> None:
    """Write a MuData to an ``.h5mu`` file, keeping each TreeData modality's trees.

    Parameters
    ----------
    mdata
        The data; its modalities may be TreeData or AnnData.
    path
        Output file; an existing file is overwritten.

    Raises
    ------
    TypeError
        ``mdata`` is not a MuData.
    ValueError
        ``mdata`` is backed by a file; load it into memory first.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/multiomics`

    ``mudata.MuData.write_h5mu`` writes a TreeData as an AnnData and drops its
    trees. This function writes the file with mudata, then adds each TreeData
    modality's trees (``obst``, ``vart``, ``label``, ``allow_overlap``,
    ``alignment``) under its group. ``mudata.read_h5mu`` still reads the file and
    ignores the trees. Only ``.h5mu`` is supported, not zarr. ``mdata`` is not
    changed; the write goes through a copy, so it needs memory for one more copy
    of the data.

    Examples
    --------
    >>> import tempfile, pathlib
    >>> import biotapy as bt
    >>> mdata = bt.io.to_mudata({"taxa": bt.datasets.toy()})
    >>> path = pathlib.Path(tempfile.mkdtemp()) / "study.h5mu"
    >>> bt.io.write_h5mu(mdata, path)
    >>> list(bt.io.read_h5mu(path)["taxa"].vart)
    ['phylo']
    """
    if not isinstance(mdata, MuData):
        msg = f"mdata must be a MuData, got {type(mdata).__name__}"
        raise TypeError(msg)
    if mdata.isbacked:
        msg = "mdata is backed by a file; load it into memory first (mudata.read_h5mu(path) without backed=True)"
        raise ValueError(msg)
    # mudata's writer calls strings_to_categoricals on every modality, which would change the caller's obs/var.
    staged = mdata.copy()
    staged.write_h5mu(path)
    trees = {name: modality for name, modality in staged.mod.items() if isinstance(modality, TreeData)}
    if not trees:
        return
    with h5py.File(path, "a") as handle:
        for name, tdata in trees.items():
            write_tree_slots(handle["mod"][name], tdata)


def read_h5mu(path: str | os.PathLike[str]) -> MuData:
    """Read an ``.h5mu`` file, giving back the TreeData modalities that ``write_h5mu`` wrote.

    Parameters
    ----------
    path
        An ``.h5mu`` file, written by ``write_h5mu`` or by mudata alone.

    Returns
    -------
    MuData
        A modality is a TreeData, with its trees, where the file holds them, and an
        AnnData otherwise.

    Raises
    ------
    ValueError
        A modality's trees were written by a newer biotapy (a later layout version).

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/multiomics`

    Runs ``mudata.read_h5mu``, then rebuilds each modality marked by ``write_h5mu``.
    Only ``.h5mu`` is supported, not zarr. String ``obs`` and ``var`` columns come
    back as categoricals, as with ``mudata.read_h5mu``.

    Examples
    --------
    >>> import tempfile, pathlib
    >>> import biotapy as bt
    >>> path = pathlib.Path(tempfile.mkdtemp()) / "study.h5mu"
    >>> bt.io.write_h5mu(bt.io.to_mudata({"taxa": bt.datasets.toy()}), path)
    >>> type(bt.io.read_h5mu(path)["taxa"]).__name__
    'TreeData'
    """
    mdata = mudata.read_h5mu(path)
    modalities = cast("MutableMapping[str, AnnData | MuData]", mdata.mod)
    with h5py.File(path, "r") as handle:
        for name, modality in list(modalities.items()):
            group = handle["mod"][name]
            if TREE_SLOTS_ATTR in group.attrs and isinstance(modality, AnnData) and not isinstance(modality, TreeData):
                modalities[name] = read_tree_slots(group, modality)
    return mdata
