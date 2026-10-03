"""HUMAnN 3 and 4 tables: gene families, reactions, pathway abundance, and their regrouped or renormalised forms."""

import re
from pathlib import Path

import pandas as pd
from mudata import MuData

from biotapy._core import XKind, make_function_mudata

from ._table import _leading_lines, _numbers, _read_table

# Sample-column suffixes: HUMAnN's own ("_Abundance-RPKs", "_Abundance"), renorm --update-snames'
# ("-CPM", "-RELAB"), and the file names humann_join_tables uses when every file names its
# sample alike ("<sample>_pathabundance_cpm", as in the HMP2 merged tables).
_SUFFIX = re.compile(r"(?:_Abundance|_genefamilies|_pathabundance)?(?:[-_](?:RPKs|CPM|RELAB|cpm|relab))?$")
_COVERAGE = re.compile(r"_(?:Coverage|pathcoverage)")
# Units a header names, in this order: renorm --update-snames rewrites the sample columns but
# keeps the first cell, so a "-RELAB" column outranks an "Adjusted CPMs" first cell.
_UNITS: tuple[tuple[re.Pattern[str], XKind], ...] = (
    (re.compile(r"[-_](?:RELAB|relab)(?:\t|$)"), "relative"),
    (re.compile(r"[-_](?:CPM|cpm)(?:\t|$)|Adjusted CPMs"), "cpm"),
    (re.compile(r"RPKs(?:\t|$)"), "rpk"),
)


def read_humann(path: str | Path) -> MuData:
    r"""Read one HUMAnN table into community and per-taxon modalities.

    Parameters
    ----------
    path
        A HUMAnN 3 or 4 output table, per sample or merged by
        ``humann_join_tables``: gene families, reactions or pathway
        abundance, as written or after ``humann_regroup_table`` /
        ``humann_renorm_table``. Gzip (``.gz``) is read directly.

    Returns
    -------
    MuData
        Two modalities over the same samples:

        - ``"function"``: community rows (no ``|``), ``var`` columns
          ``name`` and ``special``;
        - ``"function_by_taxon"``: stratified rows (``ID|taxon``), ``var``
          columns ``function``, ``name``, ``taxon``, ``genus``, ``species``
          and ``special``.

        ``var_names`` are the row ids without their ``": name"`` part.
        ``special`` flags ``UNMAPPED``, ``READS_UNMAPPED``, ``UNINTEGRATED``
        and ``UNGROUPED``, which stay as features. Sample names lose the
        suffix HUMAnN adds (``_Abundance-RPKs``, ``_Abundance``, ``-CPM``,
        ``-RELAB``, or a joined file name's ``_pathabundance_cpm``).

    Raises
    ------
    ValueError
        The file is empty, not UTF-8 text, or a ``.gz`` that is not valid
        gzip; it is a pathway coverage table; the header repeats a column
        name or has an empty sample column name; a data row has more cells
        than the header, fewer cells, or no (or a blank) id; a value is
        missing, not a number, negative or not finite; a row id holds more
        than one ``|``; or sample names repeat once their suffix is removed.
        Messages name ``path``.

    Notes
    -----
    R equivalent: ``mia::importHUMAnN``
    Guide: :doc:`/guide/reading_data`

    ``uns['biotapy']['x_kind']`` comes from the header: ``-RELAB`` or
    ``_relab`` is ``"relative"``; ``-CPM``, ``_cpm`` or ``Adjusted CPMs`` is
    ``"cpm"``; ``RPKs`` is ``"rpk"``. A header without a unit (pathway
    abundance) is ``"abundance"``, even when the values are whole numbers.
    Renormalise with ``humann_renorm_table --update-snames`` so the header
    names the new unit.

    Read one table per call: a gene family table and a pathway table both
    hold ``UNMAPPED``, so they cannot share a modality.

    The reader builds one dense rows x samples ``float64`` array before
    converting to CSR (about 290 MB for HMP2's 22,113-row x 1,638-sample
    pathway table).

    The community and stratified rows are kept apart because a pathway's
    community abundance is not the sum of its strata.

    References
    ----------
    Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of
    diverse microbial communities with bioBakery 3. eLife 10:e65088.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> path = Path(tempfile.mkdtemp()) / "genefamilies.tsv"
    >>> _ = path.write_text("# Gene Family\tS1_Abundance-RPKs\nUNMAPPED\t2.0\nK1\t4.0\nK1|g__A.s__A_b\t4.0\n")
    >>> mdata = bt.io.read_humann(path)
    >>> mdata["function"].var_names.tolist(), mdata["function_by_taxon"].var_names.tolist()
    (['UNMAPPED', 'K1'], ['K1|g__A.s__A_b'])
    """
    path = Path(path)
    argument = f"path={str(path)!r}"
    # HUMAnN's rule: the last "#" line is the header; with none, the first line is. Not _table._header, which
    # skips a last "#" line without a tab and would then read a HUMAnN table's first data row as its header.
    comments, first = _leading_lines(path, argument=argument)
    header = comments[-1] if comments else first
    if _COVERAGE.search(header):
        msg = f"{argument} is a pathway coverage table; read_humann reads abundance tables"
        raise ValueError(msg)
    table = _read_table(path, header, skiprows=max(len(comments) - 1, 0), argument=argument)
    X = _numbers(table, argument=argument, nonnegative=True).T
    obs = pd.DataFrame(index=table.columns.str.replace(_SUFFIX, "", regex=True))
    # HUMAnN never writes raw counts: a table whose header names no unit holds pathway abundances.
    unit = next((kind for pattern, kind in _UNITS if pattern.search(header.rstrip("\n"))), None)
    x_kind = unit or "abundance"
    try:
        return make_function_mudata(X, obs=obs, row_ids=table.index, x_kind=x_kind, source="io.read_humann")
    except ValueError as error:  # repeated sample or row ids, or a row id with two "|"
        msg = f"{argument}: {error}"
        raise ValueError(msg) from error
