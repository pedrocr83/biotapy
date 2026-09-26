"""Canonical taxonomic ranks and their normalization (contracts/data-model-slots)."""

import re

import numpy as np
import pandas as pd
from anndata import AnnData

RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")

RANK_ALIASES = {"domain": "kingdom"}
# Greengenes/RESCRIPT "k__"-style or SILVA "D_0__"-style rank prefixes.
_PREFIX = re.compile(r"^(?:(?P<letter>[kdpcofgs])__|D_(?P<level>\d+)__)")
_PREFIX_RANK = {
    "k": "kingdom",
    "d": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
_MISSING = ("", "NA")


def split_ranks(adata: AnnData, rank: str) -> tuple[list[str], list[str]]:
    """Split the rank columns of ``var`` into up-to-and-including ``rank`` and below it."""
    present = [r for r in RANKS if r in adata.var.columns]
    if rank not in present:
        msg = f"rank={rank!r} is not a taxonomy column; available ranks: {present}"
        raise KeyError(msg)
    cut = present.index(rank) + 1
    return present[:cut], present[cut:]


def _canonical(column: object) -> object:
    name = RANK_ALIASES.get(str(column).strip().lower(), str(column).strip().lower())
    return name if name in RANKS else column


def normalize_ranks(frame: pd.DataFrame) -> pd.DataFrame:
    """Canonical lowercase rank columns with missing values as NaN (contracts/data-model-slots)."""
    out = frame.rename(columns=_canonical)
    for rank in [column for column in out.columns if column in RANKS]:
        values = out[rank].astype("string").str.strip().str.replace(_PREFIX, "", regex=True)
        out[rank] = values.astype(object).mask(values.isna() | values.isin(_MISSING), np.nan)
    return out


def _rank_and_value(part: str, position: int) -> tuple[str | None, str]:
    match = _PREFIX.match(part)
    if match is None:
        return (RANKS[position] if position < len(RANKS) else None), part
    if match["letter"] is not None:
        return _PREFIX_RANK[match["letter"]], part[match.end() :]
    level = int(match["level"])
    return (RANKS[level] if level < len(RANKS) else None), part[match.end() :]


def _parse_lineage(text: str) -> dict[str, str]:
    ranks: dict[str, str] = {}
    for position, part in enumerate(text.split(";")):
        rank, value = _rank_and_value(part.strip(), position)
        if rank is not None:
            ranks[rank] = value
    return ranks


def split_lineage(lineage: pd.Series) -> pd.DataFrame:
    """Split ``;``-separated lineages (``k__Bacteria; p__Firmicutes``) into rank columns."""
    frame = pd.DataFrame([_parse_lineage(t) if isinstance(t, str) else {} for t in lineage], index=lineage.index)
    depth = max((RANKS.index(column) + 1 for column in frame.columns), default=0)
    return normalize_ranks(frame.reindex(columns=list(RANKS[:depth])))
