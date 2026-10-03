"""Strict reading of the tab-separated tables HUMAnN, MetaPhlAn and PICRUSt2 write."""

import csv
import gzip
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def _leading_lines(path: Path, *, argument: str) -> tuple[list[str], str]:
    """The ``#`` lines that open ``path`` (gzip is read directly) and the first line after them (``""`` if none).

    A file that is not UTF-8 text, or a ``.gz`` that is not gzip, raises ``ValueError`` naming ``argument``.
    """
    comments: list[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith("#"):
                    return comments, line
                comments.append(line)
    except (UnicodeDecodeError, gzip.BadGzipFile) as error:
        msg = f"{argument} is not a UTF-8 text table (or a valid .gz): {error}"
        raise ValueError(msg) from error
    return comments, ""


def _read_table(path: Path, header: str, *, skiprows: int, argument: str, text: int = 1) -> pd.DataFrame:
    """The table whose header line is ``header``, indexed by its first column.

    ``skiprows`` lines precede the header. The first ``text`` columns are read
    as text, so ids such as ``0042`` stay as written. ``argument`` names the
    input in messages, e.g. ``"path='table.tsv'"``. An empty file, repeated
    column names, an empty column name (other than the first), a data row with more cells than the
    header or a row with no (or a blank) id raise ``ValueError``; a short row's missing cells are NaN, for the
    caller to check in the columns it reads (``_numbers``).
    """
    if header.strip() == "":
        msg = f"{argument} is empty (it has no header line)"
        raise ValueError(msg)
    names = header.rstrip("\r\n").split("\t")
    # The first cell names the id column and is never used: pandas and R's write.table(col.names=NA) leave it empty.
    if any(not name.strip() for name in names[1:]):
        msg = f"{argument} has an empty column name in its header (a trailing tab?)"
        raise ValueError(msg)
    repeated = sorted(name for name, count in Counter(names).items() if count > 1)
    if repeated:
        msg = f"{argument} repeats column names {repeated[:3]}"
        raise ValueError(msg)
    try:
        table = pd.read_csv(
            path,
            sep="\t",
            skiprows=skiprows,
            index_col=0,
            dtype=dict.fromkeys(range(text), str),
            quoting=csv.QUOTE_NONE,
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError, gzip.BadGzipFile) as error:
        msg = f"{argument} is not a valid tab-separated table: {error}"
        raise ValueError(msg) from error
    # pandas shifts the header over when the first data row is longer, so compare with the header's own cells.
    if table.shape[1] != len(names) - 1:
        msg = f"{argument} has a data row with more cells than the header"
        raise ValueError(msg)
    if np.any(table.index.isna() | (table.index.str.strip() == "")):
        msg = f"{argument} has a data row with no id (an empty first cell)"
        raise ValueError(msg)
    return table


def _numbers(frame: pd.DataFrame, *, argument: str, nonnegative: bool = False) -> np.ndarray:
    """``frame``'s values as a float64 array, raising ``ValueError`` naming ``argument`` on a non-number, a gap or an infinite value.

    With ``nonnegative``, a negative value raises too (abundances cannot be negative).
    """
    try:
        values = frame.to_numpy(dtype=np.float64)
    except ValueError as error:
        msg = f"{argument} has a value that is not a number: {error}"
        raise ValueError(msg) from error
    if np.isnan(values).any():
        msg = f"{argument} has a missing or NaN value (a data row with fewer cells than the header, or an empty cell)"
        raise ValueError(msg)
    if not np.isfinite(values).all():
        msg = f"{argument} has a value that is not finite (inf)"
        raise ValueError(msg)
    if nonnegative and (values < 0).any():
        msg = f"{argument} has a negative value; abundances cannot be negative"
        raise ValueError(msg)
    return values
