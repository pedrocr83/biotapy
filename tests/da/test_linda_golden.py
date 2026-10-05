from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
pytestmark = [pytest.mark.golden, pytest.mark.network]


@pytest.mark.parametrize(("formula", "covariates"), [("host", ()), ("host + log_depth", ("log_depth",))])
def test_linda_matches_microbiomestat_linda(benchmark, formula, covariates):
    # R: MicrobiomeStat::linda(counts, meta, "~<formula>", is.winsor = FALSE)$output$hosthuman on the 636 genera.
    golden = pd.read_csv(GOLDEN / "linda.csv.gz", dtype={"taxon_id": str}).query("formula == @formula")
    out = bt.da.linda(benchmark, "host", covariates=covariates, reference="other")
    assert out.index.tolist() == golden["taxon_id"].tolist()
    assert (out["contrast"] == "human vs other").all()
    columns = {"effect": "log2FoldChange", "se": "lfcSE", "pvalue": "pvalue", "qvalue": "padj"}
    for ours, theirs in columns.items():
        np.testing.assert_allclose(out[ours], golden[theirs], rtol=1e-7, err_msg=ours)
