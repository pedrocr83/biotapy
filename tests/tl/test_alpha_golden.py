from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


def test_alpha_matches_phyloseq_estimate_richness():
    golden = pd.read_csv(GOLDEN / "alpha.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
    out = bt.tl.alpha(bt.datasets.global_patterns())
    assert list(out.index) == list(golden.index)
    for metric, column in [
        ("observed_features", "Observed"),
        ("shannon", "Shannon"),
        ("simpson", "Simpson"),
        ("chao1", "Chao1"),
    ]:
        np.testing.assert_allclose(out[metric], golden[column], rtol=1e-7, err_msg=metric)


def test_faith_pd_matches_picante_pd_with_root():
    golden = pd.read_csv(GOLDEN / "alpha_faith_pd.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
    out = bt.tl.alpha(bt.datasets.global_patterns(), metrics=["faith_pd"])
    assert list(out.index) == list(golden.index)
    np.testing.assert_allclose(out["faith_pd"], golden["pd"], rtol=1e-7)
