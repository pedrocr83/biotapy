import numpy as np
import pandas as pd
import pytest

import biotapy as bt


# Session scope: GlobalPatterns takes seconds to load and glom; tests only read it.
@pytest.fixture(scope="session")
def benchmark():
    """The exit-gate data the golden files use: GlobalPatterns genera in >= 20% of samples, ``host`` and ``log_depth``."""
    tdata = bt.pp.filter_features(bt.pp.tax_glom(bt.datasets.global_patterns(), "genus"), min_prevalence=0.2)
    human = tdata.obs["SampleType"].isin(["Feces", "Skin", "Tongue"])
    tdata.obs["host"] = pd.Categorical(np.where(human, "human", "other"))
    tdata.obs["log_depth"] = np.log(np.asarray(tdata.X.sum(axis=1)).ravel())
    return tdata
