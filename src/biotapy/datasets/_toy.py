"""Tiny in-memory dataset for docstring examples and tests."""

import numpy as np
import pandas as pd

from biotapy._core import TreeData, make_treedata, tree_from_edges

_COUNTS = np.array(
    [
        [10, 5, 20, 30, 0, 2, 1, 0],
        [8, 7, 25, 22, 3, 0, 0, 1],
        [12, 4, 18, 35, 1, 5, 2, 0],
        [2, 1, 5, 10, 0, 40, 15, 3],
        [0, 2, 3, 12, 2, 38, 20, 5],
        [1, 0, 4, 8, 1, 45, 12, 2],
    ],
    dtype=np.int64,
)
_FIRMICUTES = ("Bacteria", "Firmicutes", "Clostridia")
_BACTEROIDOTA = ("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales")
_PROTEOBACTERIA = ("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae")
_TAXONOMY = [
    (*_FIRMICUTES, "Lachnospirales", "Lachnospiraceae", "Blautia"),
    (*_FIRMICUTES, "Lachnospirales", "Lachnospiraceae", "Roseburia"),
    (*_FIRMICUTES, "Oscillospirales", "Ruminococcaceae", "Faecalibacterium"),
    (*_BACTEROIDOTA, "Bacteroidaceae", "Bacteroides"),
    (*_BACTEROIDOTA, "Bacteroidaceae", "Bacteroides"),
    (*_BACTEROIDOTA, "Prevotellaceae", "Prevotella"),
    (*_PROTEOBACTERIA, "Escherichia"),
    (*_PROTEOBACTERIA, None),
]
_EDGES = [
    ("root", "n1", 0.1),
    ("root", "n2", 0.1),
    ("root", "n3", 0.2),
    ("n1", "n4", 0.05),
    ("n4", "f1", 0.1),
    ("n4", "f2", 0.12),
    ("n1", "f3", 0.2),
    ("n2", "n5", 0.05),
    ("n5", "f4", 0.02),
    ("n5", "f5", 0.03),
    ("n2", "f6", 0.15),
    ("n3", "f7", 0.1),
    ("n3", "f8", 0.12),
]
_RANK_COLUMNS = ["kingdom", "phylum", "class", "order", "family", "genus"]


def toy() -> TreeData:
    """Six gut samples x eight features with taxonomy and a phylogeny.

    Built in memory, so examples and tests never download anything.

    Returns
    -------
    TreeData
        Counts in ``X``; ``obs['group']`` is ``A`` (s1-s3) or ``B`` (s4-s6);
        ``var`` holds kingdom to genus, with ``f8`` unassigned at genus;
        ``vart['phylo']`` carries branch lengths.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/data_model`

    Examples
    --------
    >>> import biotapy as bt
    >>> bt.datasets.toy().shape
    (6, 8)
    """
    obs = pd.DataFrame({"group": pd.Categorical(["A"] * 3 + ["B"] * 3)}, index=[f"s{i}" for i in range(1, 7)])
    var = pd.DataFrame(_TAXONOMY, index=[f"f{i}" for i in range(1, 9)], columns=_RANK_COLUMNS)
    return make_treedata(
        _COUNTS, obs=obs, var=var, tree=tree_from_edges(_EDGES), x_kind="counts", source="datasets.toy"
    )
