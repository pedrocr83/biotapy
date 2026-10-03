"""MetaPhlAn 3 and 4 taxonomic profiles, per sample or merged by merge_metaphlan_tables.py."""

from pathlib import Path

import numpy as np
import pandas as pd

from biotapy._core import RANKS, RELATIVE_TOLERANCE, TreeData, make_treedata, normalize_ranks, split_lineage

from ._table import _header, _numbers, _read_table

# NCBI taxid columns: a MetaPhlAn 3 merged table keeps one beside its samples.
_TAXID_COLUMNS = ("NCBI_tax_id", "clade_taxid")


def read_metaphlan(path: str | Path) -> TreeData:
    r"""Read a MetaPhlAn profile, or several merged, into a samples x clades table.

    Parameters
    ----------
    path
        A MetaPhlAn 3 or 4 profile (``-t rel_ab``, the default, or
        ``-t rel_ab_w_read_stats``), or a table of several merged by
        ``merge_metaphlan_tables.py`` or with one column per sample. Gzip
        (``.gz``) is read directly.

    Returns
    -------
    TreeData
        Relative abundances in ``X`` (MetaPhlAn's percentages divided by
        100), one feature per leaf clade: the deepest row of each lineage,
        such as MetaPhlAn 4's SGBs (``t__SGB1871``) or MetaPhlAn 3's species.
        ``var_names`` are the leaf's last name without its rank prefix
        (``SGB1871``, ``Bacteroides_ovatus``); ``var`` holds the rank columns
        ``kingdom`` to ``species``. ``UNCLASSIFIED`` (``UNKNOWN`` in older
        tables) stays a feature with every rank NaN, so each sample sums to 1.
        A single profile's sample is named after its file, and ``_profile``
        is removed from sample names, as ``merge_metaphlan_tables.py`` does.
        There is no tree.

    Raises
    ------
    ValueError
        The file is empty, not UTF-8 text, or a ``.gz`` that is not valid
        gzip; the header repeats a column name or has an empty sample column
        name; a data row has more cells than the header, or no (or a blank)
        id; an abundance is missing (a short row or an empty cell), not a
        number, negative or not finite; the leaf clades of a sample do not
        sum to 100 (rows removed, or a table that is not a profile, such as
        one with ``;`` lineages); or leaf names or sample names repeat.
        Messages name ``path``.

    Notes
    -----
    R equivalent: ``mia::importMetaPhlAn``
    Guide: :doc:`/guide/reading_data`

    A profile lists every rank, and a clade's abundance is the sum of its
    children's, so keeping only the leaves keeps all the abundance once.
    ``bt.pp.tax_glom`` gives back the higher ranks; it drops
    ``UNCLASSIFIED`` unless ``dropna=False``. ``uns['biotapy']['x_kind']`` is
    ``"relative"``. NCBI taxids, ``additional_species``, coverage and read
    estimates are not read.

    The reader builds one dense clades x samples ``float64`` array of every
    row the file prints before keeping the leaves as CSR (about 12 MB for
    HMP2's 932-row x 1,638-sample table).

    References
    ----------
    Blanco-Míguez A et al. (2023) Extending and improving metagenomic taxonomic profiling with
    uncharacterized species using MetaPhlAn 4. Nature Biotechnology 41:1633-1644.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> rows = ["k__Bacteria\t2\t90.0\t", "k__Bacteria|g__Bacteroides\t2|816\t90.0\t", "UNCLASSIFIED\t-1\t10.0\t"]
    >>> text = "#mpa_vJan25\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n" + "\n".join(rows)
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = Path(tmp) / "S1_profile.tsv"
    ...     _ = path.write_text(text)
    ...     tdata = bt.io.read_metaphlan(path)
    >>> tdata.obs_names.tolist(), tdata.var_names.tolist(), tdata.X.toarray().tolist()
    (['S1'], ['Bacteroides', 'UNCLASSIFIED'], [[0.9, 0.1]])
    """
    path = Path(path)
    argument = f"path={str(path)!r}"
    header, skiprows, is_profile = _header(path, argument=argument)
    table = _read_table(path, header, skiprows=skiprows, argument=argument)
    values = _sample_columns(table, path, is_profile=is_profile)
    # A leaf is a clade no other clade descends from, through any ancestor: MetaPhlAn can omit an
    # intermediate rank's row (the 4.0.6 fixture has no o__Corynebacteriales), so direct parents are not enough.
    lineages = [clade.split("|") for clade in table.index]
    ancestors = {"|".join(lineage[:depth]) for lineage in lineages for depth in range(1, len(lineage))}
    leaf = np.array([clade not in ancestors for clade in table.index], dtype=bool)
    percent = _numbers(values, argument=argument, nonnegative=True)
    X = percent[leaf].T / 100
    totals = X.sum(axis=1)
    # Only a column that is zero on every row is an empty sample; zero leaves under nonzero internal rows lost reads.
    bad = ~np.isclose(totals, 1, rtol=0, atol=RELATIVE_TOLERANCE) & percent.any(axis=0)
    if bad.any():
        shown = {
            name: round(float(total) * 100, 3)
            for name, total in zip(values.columns[bad][:3], totals[bad][:3], strict=True)
        }
        msg = (
            f"{argument}: the leaf clades of {int(bad.sum())} sample(s) do not sum to 100%: {shown}; "
            "read_metaphlan reads whole MetaPhlAn profiles with '|'-separated lineages"
        )
        raise ValueError(msg)
    clades = table.index[leaf]
    # A clade lineage is "k__A|p__B|..."; UNCLASSIFIED and UNKNOWN have none, so every rank is NaN.
    lineage = pd.Series(clades.str.replace("|", ";"), index=clades).where(clades.str.contains("__", regex=False))
    names = clades.str.split("|").str[-1].str.replace(r"^[a-z]__", "", regex=True).rename(None)
    # split_lineage stops at the deepest rank present; a profile always has all seven, missing ones NaN.
    var = normalize_ranks(split_lineage(lineage).reindex(columns=list(RANKS))).set_axis(names)
    obs = pd.DataFrame(index=values.columns.str.replace("_profile", "", regex=False))
    try:
        return make_treedata(X, obs=obs, var=var, tree=None, x_kind="relative", source="io.read_metaphlan")
    except ValueError as error:  # repeated leaf or sample names
        msg = f"{argument}: {error}"
        raise ValueError(msg) from error


def _sample_columns(table: pd.DataFrame, path: Path, *, is_profile: bool) -> pd.DataFrame:
    """The abundance columns: a profile's one, named after its file, or a merged table's samples.

    A profile's header is a ``#`` line; a merged table's is not, so a merged
    sample may be called ``relative_abundance``.
    """
    for column in ("relative_abundance", "Metaphlan2_Analysis"):
        if is_profile and column in table.columns:
            return table[[column]].set_axis([Path(path.name.removesuffix(".gz")).stem], axis=1)
    return table.drop(columns=[column for column in _TAXID_COLUMNS if column in table.columns])
