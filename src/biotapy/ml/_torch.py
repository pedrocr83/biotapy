"""A samples x features table as a PyTorch dataset (extra ``torch``, decisions/optional-heavy-dependencies)."""

import operator
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from biotapy._core import as_csr, import_optional

if TYPE_CHECKING:
    from torch import Tensor
    from torch.utils.data import Dataset

# What a dataset item is read from: CSR referenced as is, or a dense array.
Table = sp.csr_matrix | npt.NDArray[Any]


def to_torch(adata: AnnData, *, label_key: str | None = None, layer: str | None = None) -> "Dataset":
    """A PyTorch dataset over the samples, one row of features per item.

    Parameters
    ----------
    adata
        Samples x features.
    label_key
        An ``obs`` column to pair with each row. A category, string or bool
        column becomes int64 codes in category order (sorted values for a
        string column), for a classifier; a numeric column becomes float32,
        for a regression. By default an item is the row alone.
    layer
        Read ``layers[layer]``, such as ``"clr"`` from :func:`biotapy.pp.clr`,
        instead of ``X``.

    Returns
    -------
    Dataset
        A map-style :class:`torch.utils.data.Dataset` of ``adata.n_obs`` items
        in ``obs`` order. Item ``i`` is sample ``i``'s features as a 1-D
        float32 tensor, or with ``label_key`` the pair ``(features, label)``.

    Raises
    ------
    ImportError
        torch is not installed: ``pip install 'biotapy[torch]'``.
    KeyError
        ``label_key`` is not a column of ``obs``, or ``layer`` is not a layer.
    ValueError
        The ``label_key`` column has a missing value.
    TypeError
        The table is neither a NumPy array nor a SciPy sparse matrix.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/machine_learning`

    A row is densified when its item is read, so a sparse table is never dense
    in full: an item costs 4 bytes x features, and a
    :class:`torch.utils.data.DataLoader` stacks items into batches, shuffling
    them if asked. The dataset holds the table it was given instead of
    copying it (a sparse table in another format than CSR is converted to CSR
    once), and reads the labels once, at construction. Do not modify ``adata``
    while using the dataset. Each item is a new tensor, so editing it leaves
    ``adata`` unchanged. Label codes follow
    ``pd.Categorical(adata.obs[label_key]).categories``.

    Examples
    --------
    >>> import biotapy as bt
    >>> from torch.utils.data import DataLoader
    >>> dataset = bt.ml.to_torch(bt.datasets.toy(), label_key="group")
    >>> features, labels = next(iter(DataLoader(dataset, batch_size=4)))
    >>> features.shape, features.dtype, labels.tolist()
    (torch.Size([4, 8]), torch.float32, [0, 0, 0, 1])
    """
    table = _table(adata, layer)
    labels = None if label_key is None else _labels(adata, label_key)
    return _dataset(table, labels)


def _dataset(table: Table, labels: npt.NDArray[np.int64] | npt.NDArray[np.float32] | None) -> "Dataset":
    """The dataset over a table and its labels; module-level so the dataset pickles for DataLoader workers."""
    # torch is the extra `torch`, so the Dataset subclass is defined only once it imports (rules.md R4.6, R3.6);
    # mypy treats torch as Any (pyproject.toml), so the subclassing needs the ignore with or without torch installed.
    torch: Any = import_optional("torch", extra="torch")
    targets = None if labels is None else torch.from_numpy(labels)

    class AnnDataDataset(torch.utils.data.Dataset):  # type: ignore[misc]
        def __len__(self) -> int:
            return int(table.shape[0])

        def __getitem__(self, index: int) -> "Tensor | tuple[Tensor, Tensor]":
            # One row at a time, so the full table is never dense (rules.md R6.2).
            position = operator.index(index)  # a slice would silently give one CSR row but several dense rows
            row = table[position].toarray()[0] if isinstance(table, sp.csr_matrix) else table[position]
            features = torch.from_numpy(np.array(row, dtype=np.float32))
            return features if targets is None else (features, targets[position].clone())

        def __reduce__(self) -> tuple[Any, tuple[Table, Any]]:
            # A class defined in a function cannot be pickled, which spawn and forkserver workers need.
            return (_dataset, (table, labels))

    return AnnDataDataset()


def _table(adata: AnnData, layer: str | None) -> Table:
    """``X`` or ``layers[layer]``, referenced: CSR or dense as is, any other sparse format as CSR."""
    if layer is not None and layer not in adata.layers:
        msg = f"layer={layer!r} is not in adata.layers"
        raise KeyError(msg)
    values = adata.X if layer is None else adata.layers[layer]
    if isinstance(values, np.ndarray):
        return values
    if sp.issparse(values):
        return as_csr(values)
    name = "adata.X" if layer is None else f"layers[{layer!r}]"
    msg = f"{name} is a {type(values).__name__}; to_torch reads a NumPy array or a SciPy sparse matrix"
    raise TypeError(msg)


def _labels(adata: AnnData, label_key: str) -> npt.NDArray[np.int64] | npt.NDArray[np.float32]:
    """``obs[label_key]`` as int64 codes in category order, or float32 for a numeric column."""
    if label_key not in adata.obs.columns:
        msg = f"label_key={label_key!r} is not a column of obs"
        raise KeyError(msg)
    # anndata types obs columns as Series | DataArray (its lazy variant); the data model guarantees a Series.
    values = cast("pd.Series", adata.obs[label_key])
    missing = int(values.isna().sum())
    if missing:
        msg = f"label_key={label_key!r} has {missing} missing value(s); drop those samples or fill them first"
        raise ValueError(msg)
    if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
        return values.to_numpy(dtype=np.float32, copy=True)
    return pd.Categorical(values).codes.astype(np.int64)
