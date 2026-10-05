"""The one result table every da method returns (contracts/data-model-slots, DA results)."""

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import false_discovery_control


def result(
    features: "pd.Index[str]",
    *,
    effect: npt.NDArray[np.float64],
    se: npt.NDArray[np.float64],
    pvalue: npt.NDArray[np.float64],
    method: str,
    contrast: str,
) -> pd.DataFrame:
    """The schema table: one row per feature, BH ``qvalue`` over the finite p-values, ``direction`` = sign(effect)."""
    tested = np.isfinite(pvalue)
    qvalue = np.full(pvalue.shape, np.nan)
    qvalue[tested] = false_discovery_control(pvalue[tested], method="bh")
    return pd.DataFrame(
        {
            "effect": effect,
            "se": se,
            "pvalue": pvalue,
            "qvalue": qvalue,
            "direction": np.sign(np.nan_to_num(effect)).astype(np.int8),
            "method": method,
            "contrast": contrast,
        },
        index=pd.Index(features, name="feature"),
    )
