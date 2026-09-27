from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


def test_permanova_matches_vegan_adonis2():
    golden = pd.read_csv(GOLDEN / "permanova_sampletype.csv.gz").iloc[0]
    tdata = bt.datasets.global_patterns()
    bt.tl.beta(tdata, inplace=True)
    result = bt.tl.permanova(tdata, "SampleType", permutations=9999, seed=0)
    np.testing.assert_allclose(result["test statistic"], golden["f"], rtol=1e-7)
    assert result["number of groups"] == golden["df"] + 1
    assert abs(result["p-value"] - golden["p"]) <= 0.02
