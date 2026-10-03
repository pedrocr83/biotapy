from pathlib import Path

import mudata
import numpy as np
import pandas as pd
import pytest

import biotapy as bt

TESTS = Path(__file__).parents[1]
pytestmark = pytest.mark.golden


def _as_table(mdata: mudata.MuData) -> pd.DataFrame:
    """Both modalities side by side, samples x HUMAnN row ids (``ID`` or ``ID|taxon``)."""
    mods = [mdata[key] for key in ("function", "function_by_taxon")]
    return pd.concat([pd.DataFrame(m.X.toarray(), index=m.obs_names, columns=m.var_names) for m in mods], axis=1)


@pytest.mark.parametrize("agg", ["sum", "mean"])
def test_func_glom_matches_humann_regroup_table(agg):
    golden = pd.read_csv(TESTS / "golden" / "humann" / f"regroup_{agg}.csv.gz", index_col="sample_id")
    mdata = bt.io.read_humann(TESTS / "data" / "humann" / "genefamilies.tsv")
    hierarchy = bt.fn.load_hierarchy(TESTS / "data" / "humann" / "regroup_map.tsv", "group")
    out = mudata.MuData(
        {key: bt.fn.func_glom(mod, "group", hierarchy=hierarchy, agg=agg) for key, mod in mdata.mod.items()}
    )
    table = _as_table(out)
    assert sorted(table.columns) == sorted(golden.columns)
    np.testing.assert_allclose(table[golden.columns].to_numpy(), golden.to_numpy(), rtol=1e-7)


def test_golden_files_come_from_humann_3_9():
    assert (TESTS / "golden" / "humann" / "VERSIONS.txt").read_text().startswith("humann 3.9\n")
