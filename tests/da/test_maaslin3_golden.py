from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
# Marker r only: it needs R, and the r-bridge CI job gives it the network job's pooch cache.
pytestmark = pytest.mark.r


@pytest.mark.parametrize(("formula", "covariates"), [("host", ()), ("host + log_depth", ("log_depth",))])
def test_maaslin3_matches_maaslin3_maaslin3(benchmark, formula, covariates):
    # R: set.seed(1165433077); maaslin3(counts, meta, output, formula, min_prevalence = 0, evaluate_only = "abundance",
    # warn_prevalence = FALSE, subtract_median = TRUE), the host rows. The seed is the integer seed=20260927 becomes, so
    # the median test's simulations are the same and the numbers match to rounding.
    golden = pd.read_csv(GOLDEN / "maaslin3.csv.gz", dtype={"taxon_id": str}).query("formula == @formula")
    golden = golden.set_index("taxon_id").reindex(benchmark.var_names)
    out = bt.da.maaslin3(benchmark, "host", covariates=covariates, reference="other", seed=20260927)
    assert (out["contrast"] == "human vs other").all()
    # The 36 genera whose fit reports an error (no read in one host group): maaslin3 leaves them out of its q-values.
    untested = golden["error"]
    assert untested.sum() == 36 and out["effect"].isna().equals(untested.rename("effect"))
    ours, theirs = out[~untested], golden[~untested]
    for column, name in {"effect": "coef", "se": "stderr", "pvalue": "pval"}.items():
        np.testing.assert_allclose(ours[column], theirs[name], rtol=1e-7, err_msg=column)
    # BH over the group's p-values; maaslin3's qval pools them with the covariate's (20 calls against 45 with log_depth).
    np.testing.assert_allclose(ours["qvalue"], false_discovery_control(theirs["pval"]), rtol=1e-7)
