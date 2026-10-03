from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import biotapy as bt

TESTS = Path(__file__).parents[1]
pytestmark = pytest.mark.golden


CASES = {
    "renorm_relab": ("pathabundance.tsv", "relab", True),
    "renorm_cpm": ("pathabundance.tsv", "cpm", True),
    "renorm_relab_nospecial": ("pathabundance.tsv", "relab", False),
    "renorm_cpm_genefamilies": ("genefamilies.tsv", "cpm", True),
}


@pytest.mark.parametrize("golden_name", CASES)
def test_renorm_matches_humann_renorm_table(golden_name):
    source, units, special = CASES[golden_name]
    golden = pd.read_csv(TESTS / "golden" / "humann" / f"{golden_name}.csv.gz", index_col="sample_id")
    mdata = bt.io.read_humann(TESTS / "data" / "humann" / source)
    # The fixture's S3 has no abundance, so HUMAnN and renorm both warn and keep it zero.
    with pytest.warns(UserWarning, match="S3"):
        out = bt.fn.renorm(mdata, units, special=special)
    mods = [out[key] for key in ("function", "function_by_taxon")]
    table = pd.concat([pd.DataFrame(m.X.toarray(), index=m.obs_names, columns=m.var_names) for m in mods], axis=1)
    assert sorted(table.columns) == sorted(golden.columns)
    # Looser than the 1e-7 default: humann_renorm_table prints six significant digits (%.6g).
    np.testing.assert_allclose(table[golden.columns].to_numpy(), golden.to_numpy(), rtol=5e-6)
