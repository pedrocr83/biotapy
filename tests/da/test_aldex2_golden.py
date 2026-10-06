from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
# Marker r only: it needs R, and the r-bridge CI job gives it the network job's pooch cache.
pytestmark = pytest.mark.r


def test_aldex2_matches_aldex2_aldex(benchmark):
    # R: set.seed(1165433077); ALDEx2::aldex(counts, as.character(host), mc.samples = 128, test = "t", effect = TRUE,
    # denom = "all"). 1165433077 is the integer seed=20260927 becomes, and reference="human" keeps R's sorted level
    # order, so the Monte Carlo draws are the same and the numbers match to rounding. The integer is NumPy's
    # default_rng(20260927).integers(2**31 - 1): a NumPy change to that stream fails this test loudly; recompute it and
    # re-export the golden (tests/r/export_golden.R).
    golden = pd.read_csv(GOLDEN / "aldex2.csv.gz", dtype={"taxon_id": str}).set_index("taxon_id")
    out = bt.da.aldex2(benchmark, "host", reference="human", seed=20260927)
    assert out.index.tolist() == golden.index.tolist()
    assert (out["contrast"] == "other vs human").all()
    np.testing.assert_allclose(out["effect"], golden["diff_btw"], rtol=1e-7)
    np.testing.assert_allclose(out["pvalue"], golden["we_ep"], rtol=1e-7)
    np.testing.assert_allclose(out["qvalue"], false_discovery_control(golden["we_ep"]), rtol=1e-7)
