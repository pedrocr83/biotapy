"""Function hierarchies from the user's own mapping files (decisions/no-bundled-kegg)."""

import gzip
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


def load_hierarchy(
    path: str | Path, level: str, *, layout: Literal["parent_first", "child_first"] = "parent_first"
) -> pd.DataFrame:
    r"""Read a function mapping file into an edge table for ``func_glom``.

    Parameters
    ----------
    path
        A tab-separated mapping file with no header line, gzip (``.gz``)
        allowed. Each line holds one id and then one or more ids it maps
        to or from (see ``layout``); a line may repeat its first id.
    level
        Name for the level the parents form, for example ``"pathway"``;
        ``func_glom`` selects edges by it.
    layout
        ``"parent_first"``: ``parent<TAB>child<TAB>child...``, the format of
        ``humann_regroup_table --custom`` and PICRUSt2 mapping files.
        ``"child_first"``: ``child<TAB>parent<TAB>parent...``, the format of
        ``humann_regroup_table --reversed`` and of two-column ``child, parent``
        tables.

    Returns
    -------
    pandas.DataFrame
        Columns ``child``, ``parent``, ``level`` and ``parent_name`` (NaN:
        mapping files name no parents), one row per distinct pair.
        ``attrs["source"]`` is the file's path.

    Raises
    ------
    ValueError
        ``layout`` is not one of the two names, or a line holds a single id.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/function`

    biotapy ships no KEGG or MetaCyc mapping: their licences forbid
    redistributing them. Point ``path`` at files you are licensed to use.
    For EC numbers, ``bt.datasets.enzyme()`` gives the open ENZYME
    hierarchy.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> path = Path(tempfile.mkdtemp()) / "map.tsv"
    >>> _ = path.write_text("P1\tK1\tK2\nP2\tK2\n")
    >>> bt.fn.load_hierarchy(path, "pathway")[["child", "parent"]].values.tolist()
    [['K1', 'P1'], ['K2', 'P1'], ['K2', 'P2']]
    """
    if layout not in ("parent_first", "child_first"):
        msg = f"layout={layout!r} must be 'parent_first' or 'child_first'"
        raise ValueError(msg)
    path = Path(path)
    opener = gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")
    with opener as handle:
        rows = [(number, line.rstrip("\r\n").split("\t")) for number, line in enumerate(handle, start=1)]
    cells = [(number, [cell for cell in row if cell]) for number, row in rows if any(row)]
    short = [number for number, row in cells if len(row) < 2]
    if short:
        msg = f"path={str(path)!r}: every line needs an id and at least one id it maps to; lines {short[:3]} have one"
        raise ValueError(msg)
    pairs = [(row[0], other) for _, row in cells for other in row[1:]]
    first, other = ("parent", "child") if layout == "parent_first" else ("child", "parent")
    edges = pd.DataFrame(pairs, columns=[first, other], dtype=str)[["child", "parent"]].drop_duplicates()
    edges = edges.reset_index(drop=True).assign(level=level, parent_name=np.nan)
    edges["parent_name"] = edges["parent_name"].astype(pd.StringDtype(na_value=np.nan))
    edges.attrs["source"] = str(path)
    return edges
