from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


def test_relative_matches_phyloseq_transform_sample_counts():
    # OTU ids look numeric; read them as text so they match var_names.
    golden = pd.read_csv(GOLDEN / "relative.csv.gz", dtype={"sample_id": str, "taxon_id": str})
    out = bt.pp.relative(bt.datasets.global_patterns())
    rel = out.layers["relative"].tocoo()
    ours = pd.DataFrame({"sample_id": out.obs_names[rel.row], "taxon_id": out.var_names[rel.col], "value": rel.data})
    merged = golden.merge(ours, on=["sample_id", "taxon_id"], how="outer", suffixes=("_r", "_py"), indicator=True)
    assert (merged["_merge"] == "both").all()
    np.testing.assert_allclose(merged["value_py"], merged["value_r"], rtol=1e-7)


def test_clr_matches_vegan_decostand():
    # R: vegan::decostand(x, "clr", pseudocount = 0.5) on all of GlobalPatterns; every 100th taxon is kept.
    golden = pd.read_csv(GOLDEN / "clr.csv.gz", dtype={"sample_id": str, "taxon_id": str})
    out = bt.pp.clr(bt.datasets.global_patterns())
    rows, cols = out.obs_names.get_indexer(golden["sample_id"]), out.var_names.get_indexer(golden["taxon_id"])
    assert (rows >= 0).all() and (cols >= 0).all()
    np.testing.assert_allclose(out.layers["clr"][rows, cols], golden["value"], rtol=1e-7)
