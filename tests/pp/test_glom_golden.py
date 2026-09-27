from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


@pytest.mark.parametrize("rank", ["phylum", "genus"])
def test_tax_glom_matches_phyloseq(rank):
    golden = pd.read_csv(GOLDEN / f"tax_glom_{rank}.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
    out = bt.pp.tax_glom(bt.datasets.global_patterns(), rank)
    assert list(out.obs_names) == list(golden.index)
    assert list(out.var_names) == list(golden.columns)
    np.testing.assert_allclose(out.X.toarray(), golden.to_numpy(), rtol=1e-7)
