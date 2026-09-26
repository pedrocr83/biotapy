"""Canonical taxonomic ranks (contracts/data-model-slots)."""

from anndata import AnnData

RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")


def split_ranks(adata: AnnData, rank: str) -> tuple[list[str], list[str]]:
    """Split the rank columns of ``var`` into up-to-and-including ``rank`` and below it."""
    present = [r for r in RANKS if r in adata.var.columns]
    if rank not in present:
        msg = f"rank={rank!r} is not a taxonomy column; available ranks: {present}"
        raise KeyError(msg)
    cut = present.index(rank) + 1
    return present[:cut], present[cut:]
