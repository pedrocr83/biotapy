"""The one result table every da method returns (contracts/data-model-slots, DA results)."""

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import false_discovery_control

COLUMNS = ("effect", "se", "pvalue", "qvalue", "direction", "method", "contrast")


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


def validate_result(table: object, *, arg: str) -> pd.DataFrame:
    """``table`` if it is a result table of one method and contrast; ``arg`` names it in errors."""
    if not isinstance(table, pd.DataFrame):
        msg = f"{arg} must be a result table of a bt.da method, got {type(table).__name__}"
        raise TypeError(msg)
    missing = [column for column in COLUMNS if column not in table.columns]
    if missing:
        msg = f"{arg} lacks the result columns {missing}; pass tables that bt.da methods return"
        raise ValueError(msg)
    floats = table[["effect", "se", "pvalue", "qvalue"]].dtypes
    if not all(pd.api.types.is_float_dtype(dtype) for dtype in floats) or not pd.api.types.is_integer_dtype(
        table["direction"]
    ):
        msg = f"{arg}: effect, se, pvalue and qvalue must be floats and direction an integer"
        raise ValueError(msg)
    if not table.index.is_unique:
        msg = f"{arg} repeats features {table.index[table.index.duplicated()].unique()[:3].tolist()}"
        raise ValueError(msg)
    _check_values(table, arg=arg)
    return table


def _check_values(table: pd.DataFrame, *, arg: str) -> None:
    """Raise unless p-values, q-values and effects are missing together, p and q are probabilities, direction is sign(effect), one method and contrast."""
    probabilities = table[["pvalue", "qvalue"]]
    if not (probabilities.isna() | probabilities.ge(0) & probabilities.le(1)).all().all():
        msg = f"{arg}: pvalue and qvalue must lie between 0 and 1 (NaN for an untested feature)"
        raise ValueError(msg)
    if not table["pvalue"].isna().equals(table["qvalue"].isna()):
        msg = f"{arg}: qvalue must be NaN exactly where pvalue is"
        raise ValueError(msg)
    if not table["effect"].isna().equals(table["pvalue"].isna()):
        msg = f"{arg}: effect must be NaN exactly where pvalue is (an untested feature has no estimate)"
        raise ValueError(msg)
    if not (table["direction"].to_numpy() == np.sign(np.nan_to_num(table["effect"].to_numpy()))).all():
        msg = f"{arg}: direction must be the sign of effect, 0 where effect is NaN"
        raise ValueError(msg)
    for column in ("method", "contrast"):
        if table[column].nunique(dropna=False) != 1:
            msg = f"{arg} must hold one {column}, found {table[column].unique().tolist()}"
            raise ValueError(msg)
