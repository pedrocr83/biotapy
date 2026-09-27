from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden"
pytestmark = [pytest.mark.golden, pytest.mark.network]


def _square(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"sample_id": str}).set_index("sample_id")


@pytest.mark.parametrize("metric", ["braycurtis", "jaccard"])
def test_beta_matches_phyloseq_distance(metric):
    # jaccard: phyloseq::distance(GP, "jaccard", binary = TRUE).
    golden = _square(GOLDEN / "global_patterns" / f"beta_{metric}.csv.gz")
    out = bt.tl.beta(bt.datasets.global_patterns(), metric=metric)
    assert list(out.index) == list(golden.index) and list(out.columns) == list(golden.columns)
    np.testing.assert_allclose(out.to_numpy(), golden.to_numpy(), rtol=1e-7)


@pytest.mark.parametrize("dataset", ["global_patterns", "esophagus"])
@pytest.mark.parametrize("weighted", [False, True])
def test_unifrac_matches_phyloseq_unifrac(dataset, weighted):
    name = "weighted" if weighted else "unweighted"
    golden = _square(GOLDEN / dataset / f"unifrac_{name}.csv.gz")
    out = bt.tl.unifrac(getattr(bt.datasets, dataset)(), weighted=weighted)
    assert list(out.index) == list(golden.index)
    np.testing.assert_allclose(out.to_numpy(), golden.to_numpy(), rtol=1e-7)
