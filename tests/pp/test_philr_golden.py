from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt
from biotapy._core import get_skbio_tree

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


def test_philr_matches_philr_philr():
    # R: philr::philr(x, tree, pseudocount = 0.5) on the 293 GlobalPatterns taxa with more than 3 reads in over
    # half of the samples. R and biotapy name nodes differently, so a balance is matched by its partition:
    # the taxa in its numerator, then in its denominator.
    sbp = pd.read_csv(GOLDEN / "philr_sbp.csv.gz", dtype={"balance": str, "taxon_id": str})
    golden = pd.read_csv(GOLDEN / "philr.csv.gz", dtype={"sample_id": str, "balance": str})
    r_names = {
        (frozenset(group.loc[group["sign"] > 0, "taxon_id"]), frozenset(group.loc[group["sign"] < 0, "taxon_id"])): name
        for name, group in sbp.groupby("balance")
    }
    gp = bt.datasets.global_patterns()
    # TreeData subsetting keeps one-child nodes, which ape::drop.tip removed in R.
    tdata = gp[:, gp.var_names.isin(sbp["taxon_id"])].copy()
    out = bt.pp.philr(tdata).obsm["X_philr"]
    tree = get_skbio_tree(tdata)
    ours = {
        tuple(frozenset(tip.name for tip in child.tips(include_self=True)) for child in tree.find(name).children): name
        for name in out.columns
    }
    # Equal partitions with equal numerators: the children are in R's order, so the signs agree.
    assert ours.keys() == r_names.keys()
    renamed = out.rename(columns={name: r_names[key] for key, name in ours.items()})
    expected = golden.pivot(index="sample_id", columns="balance", values="value")
    # atol: a balance between taxa that are all absent from a sample is about 1e-16 on both sides.
    np.testing.assert_allclose(renamed.loc[expected.index, expected.columns], expected, rtol=1e-7, atol=1e-12)
