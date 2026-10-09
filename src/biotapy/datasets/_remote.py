"""Datasets downloaded once and cached with pooch: phyloseq's examples, the ENZYME files and the HMP2 tables."""

from functools import cache

import pooch

from biotapy._core import TreeData, make_pooch
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
    # HMP2 (IBDMDB) products: HUMAnN 3 pathways and MetaPhlAn 3 profiles of 2018-05-04, metadata of 2018-08-20.
    "pathabundances_3.tsv.gz": "sha256:dd983871b0e155255844b91ec10d50fb09230d2f4e915464ab680fa3a9c9ddb3",
    "taxonomic_profiles_3.tsv.gz": "sha256:d790ff15e46d61ca0cadc55d9f918de4e3415d7f97c992ac37610aaee02117ed",
    "hmp2_metadata_2018-08-20.csv": "sha256:656b7bd97660ddb875548805e30bede31f2d1208293f7170d2d5755e33862ec9",
}
_IBDMDB = "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/"
_URLS = {
    "enzyme.dat": "https://ftp.expasy.org/databases/enzyme/enzyme.dat",
    "enzclass.txt": "https://ftp.expasy.org/databases/enzyme/enzclass.txt",
    "pathabundances_3.tsv.gz": f"{_IBDMDB}products/HMP2/MGX/2018-05-04/pathabundances_3.tsv.gz",
    "taxonomic_profiles_3.tsv.gz": f"{_IBDMDB}products/HMP2/MGX/2018-05-04/taxonomic_profiles_3.tsv.gz",
    "hmp2_metadata_2018-08-20.csv": f"{_IBDMDB}metadata/hmp2_metadata_2018-08-20.csv",
}


@cache
def _pooch() -> pooch.Pooch:
    return make_pooch(_BASE_URL, _REGISTRY, urls=_URLS)


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
