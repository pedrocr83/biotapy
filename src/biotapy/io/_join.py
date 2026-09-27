"""Checked joins of side files (taxonomy, metadata) onto a table's ids."""

import numpy as np
import pandas as pd

from biotapy._core import warn_user


def _join_to(frame: pd.DataFrame, ids: "pd.Index[str]", *, argument: str) -> pd.DataFrame:
    """``frame`` reindexed to the table's ``ids``; rows for other ids are ignored.

    ``argument`` names the input in messages, e.g. ``"taxonomy='taxonomy.qza'"``.
    Repeated ids in ``frame`` or no shared id raise ``ValueError``; table ids
    missing from ``frame`` get NaN and one ``UserWarning`` with their count.
    """
    repeated = frame.index[frame.index.duplicated()].unique().tolist()
    if repeated:
        msg = f"{argument} repeats ids: {repeated[:5]}"
        raise ValueError(msg)
    shared = ids.isin(frame.index)
    if not shared.any():
        msg = (
            f"{argument} shares no ids with the table; "
            f"table ids: {ids[:3].tolist()}, its ids: {frame.index[:3].tolist()}"
        )
        raise ValueError(msg)
    n_missing = int(np.count_nonzero(~shared))
    if n_missing:
        warn_user(f"{argument}: {n_missing} of {len(ids)} table ids have no entry; they get NaN")
    return frame.reindex(ids)
