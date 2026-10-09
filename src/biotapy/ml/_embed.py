"""Sample embeddings from models that plugins register in the entry-point group ``biotapy.embeddings``."""

from importlib.metadata import entry_points

import numpy as np
import numpy.typing as npt
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
        ``biotapy.embeddings``. biotapy registers ``"mgm"``.
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
        The plugin returns something other than a NumPy array.
    ValueError
        Two installed packages register ``model``, or the plugin's array is not
        2-D, float and finite with one row per sample. Errors about the result
        name the plugin.

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
    pretrained model; it matches MGM 0.5.8's own forward pass to 2e-6. The
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
    found = [point for point in entry_points(group=GROUP) if point.name == model]
    if not found:
        installed = sorted({point.name for point in entry_points(group=GROUP)})
        msg = f"model={model!r} is not an installed embedding plugin; installed: {installed}"
        raise KeyError(msg)
    if len(found) > 1:
        packages = sorted(point.dist.name if point.dist else point.value for point in found)
        msg = f"model={model!r} is registered by several packages ({', '.join(packages)}); uninstall all but one"
        raise ValueError(msg)
    # A plugin is not trusted: a wrong row count would pair embeddings with the wrong samples (decisions/embedding-plugins).
    result = found[0].load()(adata)
    if not isinstance(result, np.ndarray):
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
    if not inplace:
        return result
    adata.obsm[f"X_{model}"] = result
    return None
