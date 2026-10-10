"""Sample embeddings from models that plugins register in the entry-point group ``biotapy.embeddings``."""

import re
from collections.abc import Iterator
from importlib.metadata import entry_points
from typing import cast

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
from anndata import AnnData

# The entry-point group a package registers an embedding model in (decisions/embedding-plugins).
GROUP = "biotapy.embeddings"


def embed(adata: AnnData, model: str, *, inplace: bool = False) -> npt.NDArray[np.floating] | None:
    """One embedding per sample from a model that a plugin provides.

    Parameters
    ----------
    adata
        Samples x features, in the form the model reads (``"mgm"``: counts or
        relative abundances with ``var["genus"]``).
    model
        The name a plugin registers in the entry-point group
        ``biotapy.embeddings``: letters, digits, ``_``, ``-`` or ``.``, so that
        ``obsm["X_<model>"]`` is a plain key. biotapy registers ``"mgm"``.
    inplace
        Write the embedding to ``obsm[f"X_{model}"]`` and return ``None``.

    Returns
    -------
    numpy.ndarray or None
        Samples x dimensions, a float array in ``obs`` order.

    Raises
    ------
    KeyError
        No installed plugin registers ``model``; the message lists those that do.
    TypeError
        The plugin returns something other than a plain NumPy array.
    ValueError
        ``model`` has other characters, two installed packages register
        ``model``, or the plugin's array is not 2-D, float and finite with one
        row per sample. Errors about the result name the plugin.
    Exception
        What the plugin itself raises passes through with its type and a note
        that names the plugin. ``"mgm"`` raises ``ImportError`` without the
        extra, ``KeyError`` without ``var["genus"]`` and ``ValueError`` for
        negative or non-finite ``X``.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/machine_learning`

    A plugin is a callable ``embed(adata)`` that returns a samples x
    dimensions NumPy array and leaves ``adata`` unchanged. A package registers
    it under a name in its ``pyproject.toml``::

        [project.entry-points."biotapy.embeddings"]
        mymodel = "mypackage.module:embed"

    biotapy finds it when ``embed`` is called, so installing the package is
    enough. It checks the array before returning or storing it.

    ``"mgm"`` is MGM, the Microbial General Model (MIT licence), a GPT-2
    pretrained on genus profiles from MGnify. It needs the extra ``mgm``
    (``pip install 'biotapy[mgm]'``) and, at the first call, downloads the
    pretrained model (33 MB, microformer-mgm 0.5.8's wheel from PyPI, checked
    against its SHA-256) into biotapy's data cache: ``BIOTAPY_DATA_DIR`` if
    set, else pooch's per-user cache. Each feature's ``var["genus"]`` is read
    as MGM reads ``g__<genus>``: the name up to its first character other than
    a letter, digit or underscore (``Escherichia-Shigella`` is
    ``Escherichia``). Features whose genus is missing or outside MGM's 9,665
    genera are left out, with one warning that counts them, and features of
    the same genus are summed, so ``pp.tax_glom(tdata, "genus")`` first
    changes nothing. Then, as MGM's own preprocessing does, each sample becomes
    relative abundances, its genera are sorted by abundance standardised with
    MGM's per-genus mean and standard deviation, and the sentence ``<bos>``,
    genera, ``<eos>`` is cut to 512 tokens. A sample without a single known
    genus is embedded from ``<bos> <eos>``, with a warning naming it. The
    embedding is the mean of the model's last hidden layer over the sample's
    tokens (256 float32 values), the mean pooling MGM's authors use for the
    pretrained model; it matches MGM 0.5.8's own forward pass to 2e-6, on its first call
    in a process as on every other. The
    model runs on the CPU, one sample at a time: batches were slower there and
    needed up to 1.5 GB more memory. Cite Zhang et al. (2026) when you publish
    results that use it.

    References
    ----------
    Zhang H, Zhang Y, Kang Z, Xiong J, Yang R, Ning K (2026) MGM as a
    large-scale pretrained foundation model for microbiome analyses in diverse
    contexts. Adv Sci 13:e13333.

    Examples
    --------
    >>> import biotapy as bt
    >>> genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
    >>> bt.ml.embed(genera, "mgm").shape
    (6, 256)
    """
    if not re.fullmatch(r"[\w.-]+", model):
        msg = f"model={model!r} must be letters, digits, '_', '-' or '.', so that obsm['X_<model>'] is a plain key"
        raise ValueError(msg)
    found = [point for point in entry_points(group=GROUP) if point.name == model]
    if not found:
        installed = sorted({point.name for point in entry_points(group=GROUP)})
        msg = f"model={model!r} is not an installed embedding plugin; installed: {installed}"
        raise KeyError(msg)
    if len(found) > 1:
        packages = sorted({point.dist.name if point.dist else point.value for point in found})
        msg = f"model={model!r} is registered by several packages ({', '.join(packages)}); uninstall all but one"
        raise ValueError(msg)
    # A plugin is not trusted: a wrong row count would pair embeddings with the wrong samples (decisions/embedding-plugins).
    try:
        result = found[0].load()(adata)
    except Exception as err:
        err.add_note(f"while using the embedding plugin {model!r} ({found[0].value})")
        raise
    result = _checked(model, adata, result)
    if not inplace:
        return result
    adata.obsm[f"X_{model}"] = result
    return None


def _checked(model: str, adata: AnnData, result: object) -> npt.NDArray[np.floating]:
    """``result`` if it is a float, finite, samples x dimensions array that shares no memory with ``adata``'s X, layers, obsm, varm, obsp, varp or top-level uns, else an error."""
    if type(result) is not np.ndarray:  # a masked array would hide NaN from the check below, a matrix is always 2-D
        msg = f"plugin {model!r} returned a {type(result).__name__}, not a NumPy array"
        raise TypeError(msg)
    if result.ndim != 2 or result.shape[0] != adata.n_obs or result.shape[1] == 0:
        msg = (
            f"plugin {model!r} returned shape {result.shape}; expected ({adata.n_obs}, dimensions), one row per sample"
        )
        raise ValueError(msg)
    if not np.issubdtype(result.dtype, np.floating):
        msg = f"plugin {model!r} returned {result.dtype} values; expected floats"
        raise ValueError(msg)
    if not np.isfinite(result).all():
        msg = f"plugin {model!r} returned NaN or infinite values"
        raise ValueError(msg)
    if any(np.shares_memory(result, held) for held in _arrays(adata)):
        return result.copy()  # editing the embedding must not edit the caller's table
    return result


def _arrays(adata: AnnData) -> Iterator[npt.NDArray[np.generic]]:
    """The arrays adata holds: X, layers, obsm, varm, obsp, varp and the arrays in uns (a sparse matrix by its data)."""
    for held in (
        adata.X,
        *adata.layers.values(),
        *adata.obsm.values(),
        *adata.varm.values(),
        *adata.obsp.values(),
        *adata.varp.values(),
        *adata.uns.values(),
    ):
        if sp.issparse(held):
            yield cast("sp.csr_matrix", held).data
        elif isinstance(held, np.ndarray):
            yield held
