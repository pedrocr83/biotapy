"""Datasets downloaded once and cached with pooch: phyloseq's examples and the ENZYME files."""

from functools import cache
from typing import cast

import pooch

from biotapy._core import TreeData
from biotapy.io import read_phyloseq

# Pinned to one phyloseq commit so the SHA-256 hashes stay valid.
_BASE_URL = "https://raw.githubusercontent.com/joey711/phyloseq/8a6c2350b985afb909428d396081605e9b3e2f0b/data/"
_REGISTRY = {
    "GlobalPatterns.RData": "sha256:bea90c3c48275ea874e0c9400b133da1647e4cddd11b39f89a3d8ffd78512d2d",
    "enterotype.RData": "sha256:0701dd010023344a917bd31680f78580c076bf039befe829830dc43bfc56b8db",
    "esophagus.RData": "sha256:0b06d9c35f2e694c34461308af149eb54419453fcb98763de20ab61980b87e46",
    # ENZYME keeps no old releases, so no hash can stay valid: the first download is
    # cached for good, and enzyme() records the release it read (datasets/_enzyme.py).
    "enzyme.dat": None,
    "enzclass.txt": None,
}
_URLS = {
    "enzyme.dat": "https://ftp.expasy.org/databases/enzyme/enzyme.dat",
    "enzclass.txt": "https://ftp.expasy.org/databases/enzyme/enzclass.txt",
}


@cache
def _pooch() -> pooch.Pooch:
    # BIOTAPY_DATA_DIR overrides the per-user cache directory. pooch ships no py.typed
    # marker, so mypy --strict infers Any for the untyped `create`; cast it back to Pooch.
    return cast(
        pooch.Pooch,
        pooch.create(
            path=pooch.os_cache("biotapy"), base_url=_BASE_URL, registry=_REGISTRY, urls=_URLS, env="BIOTAPY_DATA_DIR"
        ),
    )


def _fetch(name: str) -> str:
    return str(_pooch().fetch(name))


def global_patterns() -> TreeData:
    """GlobalPatterns: 26 samples from 9 environments, 19,216 OTUs, with taxonomy and a tree.

    Downloaded once (435 kB) from phyloseq's repository and cached.

    Returns
    -------
    TreeData
        Counts in ``X``; sample types in ``obs``; kingdom to species in ``var``;
        the tree in ``vart['phylo']``.

    Notes
    -----
    R equivalent: ``utils::data``
    Guide: :doc:`/guide/datasets`

    In R: ``data(GlobalPatterns, package = "phyloseq")``.

    References
    ----------
    Caporaso JG et al. (2011) Global patterns of 16S rRNA diversity at a depth of millions
    of sequences per sample. PNAS 108:4516-4522.

    Examples
    --------
    >>> import biotapy as bt
    >>> bt.datasets.global_patterns().shape  # doctest: +SKIP
    (26, 19216)
    """
    return read_phyloseq(_fetch("GlobalPatterns.RData"))


def enterotype() -> TreeData:
    """enterotype: 280 gut samples, 553 genera, relative abundances; no tree.

    Returns
    -------
    TreeData
        Relative abundances in ``X`` (``x_kind == 'relative'``); sample metadata in ``obs``.

    Notes
    -----
    R equivalent: ``utils::data``
    Guide: :doc:`/guide/datasets`

    In R: ``data(enterotype, package = "phyloseq")``.

    References
    ----------
    Arumugam M et al. (2011) Enterotypes of the human gut microbiome. Nature 473:174-180.

    Examples
    --------
    >>> import biotapy as bt
    >>> bt.datasets.enterotype().shape  # doctest: +SKIP
    (280, 553)
    """
    return read_phyloseq(_fetch("enterotype.RData"))


def esophagus() -> TreeData:
    """esophagus: 3 esophageal biopsies, 58 OTUs, with a tree; no taxonomy or sample data.

    Downloaded once (2 kB) from phyloseq's repository and cached.

    Returns
    -------
    TreeData
        Counts in ``X`` (samples ``B``, ``C``, ``D``) and the tree in ``vart['phylo']``;
        ``obs`` and ``var`` have no columns.

    Notes
    -----
    R equivalent: ``utils::data``
    Guide: :doc:`/guide/datasets`

    In R: ``data(esophagus, package = "phyloseq")``.

    References
    ----------
    Pei Z et al. (2004) Bacterial biota in the human distal esophagus. PNAS 101:4250-4255.

    Examples
    --------
    >>> import biotapy as bt
    >>> bt.datasets.esophagus().shape  # doctest: +SKIP
    (3, 58)
    """
    return read_phyloseq(_fetch("esophagus.RData"))
