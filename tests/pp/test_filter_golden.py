from pathlib import Path

import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


@pytest.mark.parametrize(
    ("column", "threshold"), [("prevalence", {"min_prevalence": 0.1}), ("total", {"min_total": 5})]
)
def test_filter_features_matches_phyloseq_filter_taxa(column, threshold):
    # prevalence: filter_taxa(GP, function(x) sum(x > 0) >= 0.1 * length(x)); total: function(x) sum(x) >= 5.
    golden = pd.read_csv(GOLDEN / "filter_features.csv.gz", dtype={"taxon_id": str})
    out = bt.pp.filter_features(bt.datasets.global_patterns(), **threshold)
    assert list(out.var_names) == golden.loc[golden[column], "taxon_id"].tolist()
