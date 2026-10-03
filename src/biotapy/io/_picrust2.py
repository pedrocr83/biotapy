"""PICRUSt2 predictions: metagenome and pathway tables, their contributions, and per-ASV trait tables."""

from pathlib import Path

import pandas as pd
import scipy.sparse as sp
from mudata import MuData

from biotapy._core import make_function_mudata

from ._table import _leading_lines, _numbers, _read_table

# PICRUSt2 writes EC numbers as "EC:1.1.1.1"; ENZYME, HUMAnN and bt.datasets.enzyme write "1.1.1.1".
_EC_PREFIX = "EC:"
# The contribution columns read_picrust2 uses; the long table's first column is "sample".
_CONTRIB_COLUMNS = ("function", "taxon", "taxon_function_abun")


def read_picrust2(path: str | Path, *, contrib: str | Path | None = None) -> MuData:
    r"""Read a PICRUSt2 prediction, and optionally its contributions, into community and per-taxon modalities.

    Parameters
    ----------
    path
        An unstratified PICRUSt2 table: ``pred_metagenome_unstrat.tsv.gz``
        (EC, KO or another trait) or ``path_abun_unstrat.tsv.gz``.
    contrib
        The matching long-format contributions, written with
        ``--stratified``: ``pred_metagenome_contrib.tsv.gz`` or
        ``path_abun_contrib.tsv.gz``. Its ``taxon_function_abun`` column
        becomes the per-taxon modality.

    Returns
    -------
    MuData
        Two modalities over ``path``'s samples, as ``bt.io.read_humann``
        returns them:

        - ``"function"``: one feature per function, ``var`` columns ``name``
          and ``special``;
        - ``"function_by_taxon"``: one feature per function and taxon
          (``1.1.1.1|ASV1``), ``var`` columns ``function``, ``name``,
          ``taxon``, ``genus``, ``species`` and ``special``; no features
          without ``contrib``.

        ``taxon`` is the ASV id (or ``RARE``, PICRUSt2's group of rare
        ASVs), so ``genus`` and ``species`` are NaN. ``EC:`` is removed from
        EC numbers. ``uns['biotapy']['x_kind']`` is ``"abundance"``.

    Raises
    ------
    ValueError
        A file is empty or malformed (a value that is not a number, a missing
        value or id, a data row with more cells than the header, repeated
        column names); ``path`` has a ``description`` column; ``contrib``
        lacks a ``sample``, ``function``, ``taxon`` or
        ``taxon_function_abun`` column, repeats a sample, function and taxon,
        or names a sample or function that ``path`` lacks; a value is negative;
        a taxon written ``NA`` (read as missing: "no function or taxon").
        Messages name the file's argument.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/reading_data`

    PICRUSt2's predictions are marker-normalised read counts weighted by
    predicted gene copy numbers, not counts, so ``bt.pp.rarefy`` refuses
    them. A sample's gene-family contributions sum to its unstratified value;
    a pathway's need not, so neither modality is derived from the other.

    ``EC:`` is removed so the ids match ENZYME's and HUMAnN's
    (``bt.datasets.enzyme``); remove it from a PICRUSt2 mapping file too
    before regrouping with ``bt.fn.func_glom``.

    The unstratified table is read into one dense functions x samples
    ``float64`` array, and ``contrib`` into one pandas table of all its rows,
    before the result is stored sparse.

    References
    ----------
    Douglas GM et al. (2020) PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38:685-688.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> folder = Path(tempfile.mkdtemp())
    >>> _ = (folder / "unstrat.tsv").write_text("function\tS1\nEC:1.1.1.1\t6.0\n")
    >>> columns = "sample\tfunction\ttaxon\ttaxon_abun\ttaxon_rel_abun\tgenome_function_count"
    >>> columns += "\ttaxon_function_abun\ttaxon_rel_function_abun\tnorm_taxon_function_contrib\n"
    >>> _ = (folder / "contrib.tsv").write_text(columns + "S1\tEC:1.1.1.1\tASV1\t3.0\t100.0\t2\t6.0\t200.0\t1.0\n")
    >>> mdata = bt.io.read_picrust2(folder / "unstrat.tsv", contrib=folder / "contrib.tsv")
    >>> mdata["function"].var_names.tolist(), mdata["function_by_taxon"].var_names.tolist()
    (['1.1.1.1'], ['1.1.1.1|ASV1'])
    """
    path = Path(path)
    argument = f"path={str(path)!r}"
    table = _read(path, argument=argument)
    if "description" in table.columns:
        msg = f"{argument} has a 'description' column (add_descriptions.py output); pass the table without it"
        raise ValueError(msg)
    X = sp.csr_matrix(_numbers(table, argument=argument, nonnegative=True).T)
    row_ids = table.index.str.removeprefix(_EC_PREFIX)
    if contrib is not None:
        by_taxon, keys = _contributions(Path(contrib), samples=table.columns, functions=row_ids)
        X, row_ids = sp.hstack([X, by_taxon], format="csr"), row_ids.append(keys)
    try:
        return make_function_mudata(
            X, obs=pd.DataFrame(index=table.columns), row_ids=row_ids, x_kind="abundance", source="io.read_picrust2"
        )
    except ValueError as error:  # repeated functions, or a function or taxon id holding "|"
        msg = f"{argument}, contrib={None if contrib is None else str(contrib)!r}: {error}"
        raise ValueError(msg) from error


def _contributions(path: Path, *, samples: pd.Index, functions: pd.Index) -> tuple[sp.csr_matrix, pd.Index]:
    """The long contribution table as a samples x (function, taxon) matrix, and its ``function|taxon`` ids."""
    argument = f"contrib={str(path)!r}"
    table = _read(path, argument=argument, text=3)
    missing = [column for column in _CONTRIB_COLUMNS if column not in table.columns]
    if table.index.name != "sample" or missing:
        msg = f"{argument} needs PICRUSt2's long-format columns sample, {', '.join(_CONTRIB_COLUMNS)}; found {[table.index.name, *table.columns]}"
        raise ValueError(msg)
    values = _numbers(table[["taxon_function_abun"]], argument=argument, nonnegative=True).ravel()
    ids = table[["function", "taxon"]]
    if ids.isna().any().any():
        msg = f"{argument} has a row with no function or taxon"
        raise ValueError(msg)
    function = ids["function"].str.removeprefix(_EC_PREFIX)
    rows = samples.get_indexer(table.index)
    for name, unknown in (("samples", table.index[rows < 0]), ("functions", function[~function.isin(functions)])):
        if len(unknown):
            msg = f"{argument} names {name} that path lacks: {sorted(set(unknown))[:3]}"
            raise ValueError(msg)
    codes, keys = pd.factorize(function + "|" + ids["taxon"])
    if pd.Series(rows * len(keys) + codes).duplicated().any():
        msg = f"{argument} repeats a sample, function and taxon"
        raise ValueError(msg)
    matrix = sp.csr_matrix((values, (rows, codes)), shape=(len(samples), len(keys)))
    matrix.eliminate_zeros()
    return matrix, pd.Index(keys)


def read_picrust2_traits(path: str | Path) -> pd.DataFrame:
    r"""Read PICRUSt2's predicted gene copy numbers per ASV.

    Parameters
    ----------
    path
        A per-sequence trait table, such as ``EC_predicted.tsv.gz`` or
        ``KO_predicted.tsv.gz``: one row per ASV, one column per function.

    Returns
    -------
    pandas.DataFrame
        ASVs x functions, ``float64`` copy numbers, indexed by the ASV ids as
        written (``0042`` stays ``0042``). ``EC:`` is removed from EC numbers,
        as in ``bt.io.read_picrust2``; a ``metadata_NSTI`` column is dropped.

    Raises
    ------
    ValueError
        The file is empty or malformed (a value that is not a number, a
        missing value or id, a negative value, a data row with more cells than
        the header, repeated column names, or function ids that coincide once ``EC:`` is
        removed), or an ASV id repeats. Messages name
        ``path``.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/reading_data`

    The table describes genomes, not samples, so it is a DataFrame rather
    than a modality of the samples' MuData. It is read into one dense
    ASVs x functions ``float64`` array.

    References
    ----------
    Douglas GM et al. (2020) PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38:685-688.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> path = Path(tempfile.mkdtemp()) / "EC_predicted.tsv"
    >>> _ = path.write_text("sequence\tEC:1.1.1.1\tEC:2.7.1.1\tmetadata_NSTI\nASV1\t1\t2\t0.03\n")
    >>> bt.io.read_picrust2_traits(path)
          1.1.1.1  2.7.1.1
    ASV1      1.0      2.0
    """
    path = Path(path)
    argument = f"path={str(path)!r}"
    table = _read(path, argument=argument)
    table = table.drop(columns=[column for column in table.columns if column == "metadata_NSTI"])
    repeated = table.index[table.index.duplicated()].unique().tolist()
    if repeated:
        msg = f"{argument} repeats ASV ids: {repeated[:3]}"
        raise ValueError(msg)
    values = _numbers(table, argument=argument, nonnegative=True)
    functions = table.columns.str.removeprefix(_EC_PREFIX)
    repeated_functions = functions[functions.duplicated()].unique().tolist()
    if repeated_functions:
        msg = f"{argument} repeats function ids after removing 'EC:': {repeated_functions[:3]}"
        raise ValueError(msg)
    return pd.DataFrame(values, index=table.index.rename(None), columns=functions)


def _read(path: Path, *, argument: str, text: int = 1) -> pd.DataFrame:
    """A PICRUSt2 table: its header is the first line that does not start with ``#``."""
    comments, first = _leading_lines(path, argument=argument)
    return _read_table(path, first, skiprows=len(comments), argument=argument, text=text)
