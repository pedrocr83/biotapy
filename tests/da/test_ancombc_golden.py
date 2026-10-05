from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control, spearmanr

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


@pytest.mark.parametrize(("formula", "covariates"), [("host", ()), ("host + log_depth", ("log_depth",))])
def test_ancombc2_matches_ancombc_ancombc2(benchmark, formula, covariates):
    # R: ANCOMBC::ancombc2(counts, meta_data = meta, fix_formula, p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
    # pseudo_sens = FALSE, struc_zero = FALSE)$res, natural logs, on the 636 genera.
    golden = pd.read_csv(GOLDEN / "ancombc2.csv.gz", dtype={"taxon_id": str}).query("formula == @formula")
    golden = golden.set_index("taxon_id")
    out = bt.da.ancombc2(benchmark, "host", covariates=covariates, reference="other")
    assert out.index.tolist() == golden.index.tolist()
    # The 36 genera with no read in one host group: R reports NA with p = 1, biotapy NaN.
    untested = golden["lfc"].isna()
    assert untested.sum() == 36 and (golden.loc[untested, "p"] == 1).all()
    assert out["effect"].isna().equals(untested.rename("effect"))
    ours, theirs = out[~untested], golden[~untested]
    # The bias E-M stops at R's 100-iteration cap before converging on the host-only model, and scikit-bio stops on
    # a slightly different iterate: every effect is shifted by 0.002-0.012 log2 there (1e-8 with log_depth).
    np.testing.assert_allclose(ours["effect"], theirs["lfc"] / np.log(2), atol=0.015)
    np.testing.assert_allclose(ours["se"], theirs["se"] / np.log(2), rtol=2e-3)
    np.testing.assert_allclose(ours["pvalue"], theirs["p"], atol=0.02)
    assert spearmanr(ours["effect"], theirs["lfc"]).statistic > 0.9999
    # R's q counts the untested genera with p = 1; the schema's BH runs over the tested ones only.
    expected = false_discovery_control(theirs["p"])
    np.testing.assert_array_equal(ours["qvalue"] < 0.05, expected < 0.05)
