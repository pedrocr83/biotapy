"""DADA2 sequence tables and taxonomy, as written by R's write.csv / write.table / saveRDS."""

import re
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from biotapy._core import (
    TreeData,
    infer_x_kind,
    make_treedata,
    normalize_ranks,
    relabel_tips,
    tree_from_newick,
    tree_tips,
)

from ._join import _join_to
from ._rdata import read_matrix_rds

_SEQUENCE = re.compile(r"^[ACGTN]+$")


def read_dada2(seqtab: str | Path, *, taxa: str | Path | None = None, tree: str | Path | None = None) -> TreeData:
    r"""Read a DADA2 sequence table, with optional taxonomy and tree.

    Parameters
    ----------
    seqtab
        CSV, TSV or ``.rds`` of ``seqtab``/``seqtab.nochim``: samples x
        sequences, row names in the first column (``write.csv(seqtab.nochim,
        ...)`` or ``saveRDS(seqtab.nochim, ...)``). Sample names are kept
        verbatim as text (``001`` stays ``001``, ``NA`` is a name). Any
        ``.csv`` suffix means CSV (``seqtab.csv.gz`` works); ``.rds`` means an
        R matrix; otherwise TSV.
    taxa
        CSV, TSV or ``.rds`` of ``assignTaxonomy``/``addSpecies`` output:
        sequences x ranks. Sequences it lists that are not in ``seqtab`` are
        ignored.
    tree
        Newick file whose tips are DNA sequences, as DADA2 workflows produce
        (e.g. a tree built from ``seqtab``'s column names). Tips are relabeled
        to the matching ASV ids; sequence tips not in ``seqtab`` are pruned.

    Returns
    -------
    TreeData
        The table in ``X``; features named ``ASV1..n`` with the sequence in
        ``var['sequence']``; rank columns in ``var``; the tree in ``vart['phylo']``.

    Raises
    ------
    ValueError
        ``seqtab`` has no columns or its column names are not DNA sequences
        (the table is transposed); ``taxa`` shares no sequence with
        ``seqtab`` or repeats one; ``tree`` is not valid Newick, shares no tip
        with the table, or has a tip that is not a DNA sequence (an ASV id,
        say: biotapy numbers ASVs by column order, which need not match yours).

    Warns
    -----
    UserWarning
        Sequences missing from ``taxa`` (they get NaN); tree tips and table
        features differ (only shared features are kept).

    Notes
    -----
    R equivalent: ``phyloseq::phyloseq``
    Guide: :doc:`/guide/reading_data`

    ``uns['biotapy']['x_kind']`` is inferred from the values: non-negative whole numbers
    are ``"counts"`` (as DADA2 writes them), rows that each sum to 1 are
    ``"relative"``, anything else is ``"abundance"``.

    The text table is read densely once and stored sparse: while reading, it
    takes about 8 bytes x samples x ASVs (80 MB for 1,000 samples x 10,000
    ASVs).

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> path = Path(tempfile.mkdtemp()) / "seqtab.csv"
    >>> _ = path.write_text('"","ACGT","TTGA"\n"S1",3,0\n"S2",1,4\n')
    >>> bt.io.read_dada2(path).var_names.tolist()
    ['ASV1', 'ASV2']
    """
    counts = _read_table(Path(seqtab), argument=f"seqtab={str(seqtab)!r}")
    sequences = [str(column) for column in counts.columns]
    if not sequences or not all(_SEQUENCE.match(sequence) for sequence in sequences):
        msg = (
            f"seqtab={str(seqtab)!r} must be samples x sequences like DADA2's seqtab, "
            f"with DNA sequences as column names; found columns {sequences[:3]}"
        )
        raise ValueError(msg)
    names = [f"ASV{i}" for i in range(1, len(sequences) + 1)]
    var = pd.DataFrame({"sequence": sequences}, index=names)
    if taxa is not None:
        taxa_argument = f"taxa={str(taxa)!r}"
        taxa_table = _read_table(Path(taxa), argument=taxa_argument)
        ranks = _join_to(normalize_ranks(taxa_table), pd.Index(sequences), argument=taxa_argument)
        var = ranks.set_axis(names).join(var)
    phylo = None
    if tree is not None:
        phylo = tree_from_newick(Path(tree).read_text(), argument=f"tree={str(tree)!r}")
        _require_sequence_tips(tree_tips(phylo), tree)
        phylo = relabel_tips(phylo, dict(zip(sequences, names, strict=True)))
    obs, X = pd.DataFrame(index=counts.index), counts.to_numpy()
    # Independent of the input format (CSV gives int64; a .rds matrix keeps R's int32).
    if np.issubdtype(X.dtype, np.integer):
        X = X.astype(np.int64)
    return make_treedata(X, obs=obs, var=var, tree=phylo, x_kind=infer_x_kind(X), source="io.read_dada2")


def _read_table(path: Path, *, argument: str) -> pd.DataFrame:
    """A table whose first column holds row names, as R's write.csv writes it, or a plain .rds matrix.

    ``argument`` names the input in error messages, e.g. ``"seqtab='seqtab.rds'"``.
    """
    if path.suffix.lower() == ".rds":
        return read_matrix_rds(path, argument=argument)
    # Row names are ids: read them as text ("001" stays "001", a sample named "NA" is not
    # missing). Rank cells saying "NA" stay text here; normalize_ranks maps them to NaN.
    frame = pd.read_csv(
        path,
        sep="," if ".csv" in [suffix.lower() for suffix in path.suffixes] else "\t",
        index_col=0,
        dtype={0: str},
        keep_default_na=False,
        na_values=[""],
    )
    # write.csv leaves the corner cell empty; pandas would name the index "Unnamed: 0".
    return frame.rename_axis(index=None)


def _require_sequence_tips(tips: Iterable[str], tree: str | Path) -> None:
    others = [tip for tip in tips if not _SEQUENCE.match(tip)]
    if others:
        msg = f"tree={str(tree)!r} tips must be DNA sequences, as DADA2 workflows write them; found {others[:3]}"
        raise ValueError(msg)
