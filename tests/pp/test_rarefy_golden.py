from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


def test_rarefy_drops_the_same_samples_as_phyloseq():
    # R and NumPy generators differ, so only invariants are compared (contracts/r-golden-parity).
    golden = pd.read_csv(GOLDEN / "rarefy.csv.gz", dtype={"sample_id": str})
    tdata = bt.datasets.global_patterns()
    with pytest.warns(UserWarning, match="TRRsed1"):
        out = bt.pp.rarefy(tdata, depth=100_000, seed=0)
    assert list(out.obs_names) == golden["sample_id"].tolist()
    np.testing.assert_array_equal(np.asarray(out.X.sum(axis=1)).ravel(), golden["sample_sum"].to_numpy())
    assert (out.X.toarray() <= tdata[out.obs_names, out.var_names].X.toarray()).all()
