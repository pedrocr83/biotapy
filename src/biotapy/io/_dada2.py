"""DADA2 sequence tables and taxonomy, as written by R's write.csv / write.table."""

import re
from pathlib import Path

import pandas as pd

from biotapy._core import TreeData, make_treedata, normalize_ranks, relabel_tips, tree_from_newick

_SEQUENCE = re.compile(r"^[ACGTN]+$")


def read_dada2(seqtab: str | Path, taxa: str | Path | None = None, *, tree: str | Path | None = None) -> TreeData:
    r"""Read a DADA2 sequence table, with optional taxonomy and tree.

    Parameters
    ----------
    seqtab
        CSV or TSV of ``seqtab``/``seqtab.nochim``: samples x sequences, row
        names in the first column (``write.csv(seqtab.nochim, ...)``).
    taxa
        CSV or TSV of ``assignTaxonomy``/``addSpecies`` output: sequences x ranks.
    tree
        Newick file whose tips are sequences or ASV ids.

    Returns
    -------
    TreeData
        Counts in ``X``; features named ``ASV1..n`` with the sequence in
        ``var['sequence']``; rank columns in ``var``; the tree in ``vart['phylo']``.

    Raises
    ------
    ValueError
        ``seqtab``'s column names are not DNA sequences (the table is transposed).

    Notes
    -----
    R equivalent: ``phyloseq::phyloseq``
    Guide: :doc:`/guide/reading_data`

    The text table is read densely once and stored sparse; DADA2 tables are
    small enough for this. ``.rds`` input is not supported yet.

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
    counts = _read_table(Path(seqtab))
    sequences = [str(column) for column in counts.columns]
    if not all(_SEQUENCE.match(sequence) for sequence in sequences):
        msg = f"seqtab={str(seqtab)!r} must be samples x sequences like DADA2's seqtab; its columns are not sequences"
        raise ValueError(msg)
    names = [f"ASV{i}" for i in range(1, len(sequences) + 1)]
    var = pd.DataFrame({"sequence": sequences}, index=names)
    if taxa is not None:
        ranks = normalize_ranks(_read_table(Path(taxa))).reindex(sequences)
        var = ranks.set_axis(names).join(var)
    phylo = None
    if tree is not None:
        phylo = relabel_tips(tree_from_newick(Path(tree).read_text()), dict(zip(sequences, names, strict=True)))
    obs = pd.DataFrame(index=counts.index)
    return make_treedata(counts.to_numpy(), obs=obs, var=var, tree=phylo, x_kind="counts", source="io.read_dada2")


def _read_table(path: Path) -> pd.DataFrame:
    """A table whose first column holds row names, as R's write.csv writes it."""
    # write.csv leaves the corner cell empty; pandas would name the index "Unnamed: 0".
    return pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", index_col=0).rename_axis(index=None)
