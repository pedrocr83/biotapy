"""Synthetic benchmark data: a sparse count table with taxonomy and a balanced tree."""

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Benchmarks build their TreeData with biotapy's own constructors, as the readers do (contracts/tree-access).
from biotapy._core import TreeData, make_treedata, tree_from_edges

N_OBS, N_VARS, DENSITY, SEED = 5_000, 50_000, 0.02, 0
# Features per group at each rank: 1,000 genera, 100 families, 20 orders, 10 classes and 5 phyla at 50,000 features.
_RANK_SIZES = {"phylum": 10_000, "class": 5_000, "order": 2_500, "family": 500, "genus": 50}


def synthetic(n_obs: int = N_OBS, n_vars: int = N_VARS) -> TreeData:
    """``n_obs`` x ``n_vars`` counts from 1 to 99 at 2% density, seed 0, with taxonomy and a tree."""
    rng = np.random.default_rng(SEED)
    X = sp.random(
        n_obs, n_vars, density=DENSITY, format="csr", rng=rng, data_rvs=lambda size: rng.integers(1, 100, size)
    )
    features = [f"f{j}" for j in range(n_vars)]
    ranks = {rank: [f"{rank[0]}{j // size}" for j in range(n_vars)] for rank, size in _RANK_SIZES.items()}
    var = pd.DataFrame({"kingdom": "Bacteria"} | ranks, index=features)
    obs = pd.DataFrame(index=[f"s{i}" for i in range(n_obs)])
    tree = tree_from_edges(_balanced_edges(features, rng))
    return make_treedata(X, obs=obs, var=var, tree=tree, x_kind="counts", source="benchmarks")


def _balanced_edges(tips: list[str], rng: np.random.Generator) -> list[tuple[str, str, float]]:
    """Join neighbours pairwise until one root is left; branch lengths uniform in [0.01, 0.1)."""
    level, edges = list(tips), []
    while len(level) > 1:
        pairs = [level[i : i + 2] for i in range(0, len(level), 2)]
        level = [f"n{len(edges)}_{k}" for k in range(len(pairs))]
        for parent, pair in zip(level, pairs, strict=True):
            edges += [(parent, child, float(rng.uniform(0.01, 0.1))) for child in pair]
    return edges
