from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control, spearmanr

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


# Per model: (formula, covariates, (effect atol in log2, se rtol, pvalue atol, Spearman floor)).
# host: the bias E-M stops at R's 100-iteration cap before converging, and scikit-bio stops on a slightly different
# iterate, so every effect is shifted by 0.002-0.012 log2 (measured max 0.0121), se by 9.0e-4 relative and p by
# 0.0122; the shift is near constant, hence the rank floor. With log_depth the E-M converges on both sides and
# agreement is 1.5e-8 (effect), 1.7e-10 (se, relative) and 1.1e-8 (p), checked at 1e-6.
MODELS = [
    pytest.param(("host", (), (0.015, 2e-3, 0.02, 0.9999)), id="host"),
    pytest.param(("host + log_depth", ("log_depth",), (1e-6, 1e-6, 1e-6, 0.999999)), id="host+log_depth"),
]


@pytest.mark.parametrize("case", MODELS)
def test_ancombc2_matches_ancombc_ancombc2(benchmark, case):
    # R: ANCOMBC::ancombc2(counts, meta_data = meta, fix_formula, p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
    # pseudo_sens = FALSE, struc_zero = FALSE)$res, natural logs, on the 636 genera.
    formula, covariates, (effect_atol, se_rtol, p_atol, rank_floor) = case
    golden = pd.read_csv(GOLDEN / "ancombc2.csv.gz", dtype={"taxon_id": str}).query("formula == @formula")
    golden = golden.set_index("taxon_id")
    out = bt.da.ancombc2(benchmark, "host", covariates=covariates, reference="other")
    assert out.index.tolist() == golden.index.tolist()
    # The 36 genera with no read in one host group: R reports NA with p = 1, biotapy NaN.
    untested = golden["lfc"].isna()
    assert untested.sum() == 36 and (golden.loc[untested, "p"] == 1).all()
    assert out["effect"].isna().equals(untested.rename("effect"))
    ours, theirs = out[~untested], golden[~untested]
    np.testing.assert_allclose(ours["effect"], theirs["lfc"] / np.log(2), atol=effect_atol)
    np.testing.assert_allclose(ours["se"], theirs["se"] / np.log(2), rtol=se_rtol)
    np.testing.assert_allclose(ours["pvalue"], theirs["p"], atol=p_atol)
    assert spearmanr(ours["effect"], theirs["lfc"]).statistic > rank_floor
    # R's q counts the untested genera with p = 1; the schema's BH runs over the tested ones only.
    expected = false_discovery_control(theirs["p"])
    np.testing.assert_array_equal(ours["qvalue"] < 0.05, expected < 0.05)
