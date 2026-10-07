"""Synthetic benchmark data: a count table with taxonomy and a tree, a genus table, a function table, taxon traits."""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from mudata import MuData

# Benchmarks build their TreeData with biotapy's own constructors, as the readers do (contracts/tree-access).
from biotapy._core import TreeData, make_function_mudata, make_treedata, tree_from_edges

N_OBS, N_VARS, DENSITY, SEED = 5_000, 50_000, 0.02, 0
# Features per group at each rank: 1,000 genera, 100 families, 20 orders, 10 classes and 5 phyla at 50,000 features.
_RANK_SIZES = {"phylum": 10_000, "class": 5_000, "order": 2_500, "family": 500, "genus": 50}
# HMP2's pathway table has 478 community and 21,635 stratified rows over 1,638 samples, 53% and 7% non-zero.
N_SAMPLES, N_FUNCTIONS, N_STRATA, N_GROUPS = 1_600, 500, 43, 50
# Tian et al.'s functional redundancy at 2,000 taxa, as measured when fn.functional_redundancy was written.
N_TAXA, N_GENES, N_ABUNDANCE_SAMPLES = 2_000, 2_500, 100
# PhILR densifies X: 1,000 samples keep it at 400 MB on the 50,000-feature table.
N_PHILR_OBS = 1_000
# Differential abundance on a large genus-level study: 2,000 samples in two groups, 10,000 features, 30% non-zero.
N_DA_OBS, N_DA_VARS, DA_DENSITY = 2_000, 10_000, 0.3


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


def synthetic_function() -> MuData:
    """A 1,600-sample pathway table: 500 functions, each with 43 strata (22,000 rows), CPM-like values, seed 0."""
    rng = np.random.default_rng(SEED)
    functions = [f"F{j}" for j in range(N_FUNCTIONS)]
    strata = [f"{function}|g__G{k}.s__G{k}_sp" for function in functions for k in range(N_STRATA)]
    community = sp.random(N_SAMPLES, len(functions), density=0.5, format="csr", rng=rng)
    by_taxon = sp.random(N_SAMPLES, len(strata), density=0.07, format="csr", rng=rng)
    X = sp.hstack([community, by_taxon], format="csr") * 1_000
    obs = pd.DataFrame(index=[f"s{i}" for i in range(N_SAMPLES)])
    return make_function_mudata(X, obs=obs, row_ids=pd.Index(functions + strata), x_kind="cpm", source="benchmarks")


def function_groups() -> pd.DataFrame:
    """Each of the 500 functions in two of 50 groups, so ``func_glom`` sums many-to-many."""
    children = [f"F{j}" for j in range(N_FUNCTIONS)]
    parents = [f"G{j % N_GROUPS}" for j in range(N_FUNCTIONS)] + [f"G{(j + 25) % N_GROUPS}" for j in range(N_FUNCTIONS)]
    return pd.DataFrame({"child": children * 2, "parent": parents, "level": "group"})


def synthetic_traits() -> tuple[AnnData, pd.DataFrame]:
    """100 samples x 2,000 taxa at 5% density, and the taxa's copy numbers of 2,500 genes (70% zeros), seed 0."""
    rng = np.random.default_rng(SEED)
    taxa = [f"t{j}" for j in range(N_TAXA)]
    X = sp.random(N_ABUNDANCE_SAMPLES, N_TAXA, density=0.05, format="csr", rng=rng)
    adata = AnnData(
        X=X, obs=pd.DataFrame(index=[f"s{i}" for i in range(N_ABUNDANCE_SAMPLES)]), var=pd.DataFrame(index=taxa)
    )
    copies = rng.integers(1, 6, size=(N_TAXA, N_GENES)) * (rng.random((N_TAXA, N_GENES)) >= 0.7)
    return adata, pd.DataFrame(copies, index=taxa, columns=[f"g{k}" for k in range(N_GENES)])


def synthetic_genera() -> AnnData:
    """2,000 samples x 10,000 features, 30% non-zero counts from 1 to 99, ``obs["group"]`` a / b, seed 0."""
    rng = np.random.default_rng(SEED)
    X = sp.random(
        N_DA_OBS, N_DA_VARS, density=DA_DENSITY, format="csr", rng=rng, data_rvs=lambda size: rng.integers(1, 100, size)
    )
    group = pd.Categorical(np.repeat(["a", "b"], N_DA_OBS // 2))
    obs = pd.DataFrame({"group": group}, index=[f"s{i}" for i in range(N_DA_OBS)])
    return AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"f{j}" for j in range(N_DA_VARS)]))
