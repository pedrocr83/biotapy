"""QIIME 2 artifacts (.qza) and metadata files, read without a QIIME 2 install."""

import csv
import tempfile
import zipfile
from pathlib import Path

import biom
import numpy as np
import pandas as pd

from biotapy._core import TreeData, make_treedata, split_lineage, tree_from_newick

from ._biom import _biom_parts
from ._join import _join_to

_ID_HEADERS_ANY_CASE = frozenset({"id", "sampleid", "sample id", "sample-id", "featureid", "feature id", "feature-id"})
_ID_HEADERS_EXACT = frozenset({"#SampleID", "#Sample ID", "#OTUID", "#OTU ID", "sample_name"})


def read_qiime2(
    table: str | Path,
    *,
    taxonomy: str | Path | None = None,
    tree: str | Path | None = None,
    metadata: str | Path | None = None,
) -> TreeData:
    """Read QIIME 2 artifacts into one TreeData, without QIIME 2 installed.

    Parameters
    ----------
    table
        ``FeatureTable[Frequency]`` artifact (``.qza``).
    taxonomy
        ``FeatureData[Taxonomy]`` artifact; replaces any taxonomy in the table.
        Features it lists that are not in the table are ignored.
    tree
        ``Phylogeny[Rooted]`` or ``Phylogeny[Unrooted]`` artifact.
    metadata
        QIIME 2 sample metadata file (TSV); samples not in the table are ignored
        (metadata often covers more samples than one table).

    Returns
    -------
    TreeData
        Counts in ``X``; rank columns and ``confidence`` in ``var``; metadata
        in ``obs``; the tree in ``vart['phylo']``.

    Raises
    ------
    ValueError
        An artifact lacks the expected payload, the metadata has no ID header,
        or a metadata column declared ``numeric`` holds a non-numeric value.

    Notes
    -----
    R equivalent: ``qiime2R::qza_to_phyloseq``
    Guide: :doc:`/guide/reading_data`

    Artifacts are recognized by their payload file, not by ``metadata.yaml``.

    Examples
    --------
    >>> import tempfile
    >>> import zipfile
    >>> from pathlib import Path
    >>> import biom
    >>> from biom.util import biom_open
    >>> import biotapy as bt
    >>> toy = bt.datasets.toy()
    >>> folder = Path(tempfile.mkdtemp())
    >>> with biom_open(str(folder / "t.biom"), "w") as handle:
    ...     biom.Table(toy.X.T, list(toy.var_names), list(toy.obs_names)).to_hdf5(handle, "example")
    >>> with zipfile.ZipFile(folder / "table.qza", "w") as archive:
    ...     archive.write(folder / "t.biom", "0000/data/feature-table.biom")
    >>> bt.io.read_qiime2(folder / "table.qza").shape
    (6, 8)
    """
    with tempfile.TemporaryDirectory() as directory:
        workdir = Path(directory)
        X, obs, var = _biom_parts(
            biom.load_table(str(_payload(table, "feature-table.biom", workdir, argument="table")))
        )
        if taxonomy is not None:
            ranks = _taxonomy(_payload(taxonomy, "taxonomy.tsv", workdir, argument="taxonomy"))
            var = _join_to(ranks, var.index, argument=f"taxonomy={str(taxonomy)!r}")
        phylo = (
            None if tree is None else tree_from_newick(_payload(tree, "tree.nwk", workdir, argument="tree").read_text())
        )
    if metadata is not None:
        obs = _join_to(_metadata(Path(metadata)), obs.index, argument=f"metadata={str(metadata)!r}")
    return make_treedata(X, obs=obs, var=var, tree=phylo, x_kind="counts", source="io.read_qiime2")


def _payload(artifact: str | Path, filename: str, directory: Path, *, argument: str) -> Path:
    """Extract ``<uuid>/data/<filename>`` from a .qza into ``directory``."""
    with zipfile.ZipFile(artifact) as archive:
        members = [name for name in archive.namelist() if Path(name).parts[1:] == ("data", filename)]
        if len(members) != 1:
            msg = f"{argument}={str(artifact)!r} is not a QIIME 2 artifact holding data/{filename}"
            raise ValueError(msg)
        return Path(archive.extract(members[0], directory))


def _taxonomy(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
    ranks = split_lineage(frame["Taxon"])
    if "Confidence" in frame.columns:
        ranks["confidence"] = pd.to_numeric(frame["Confidence"], errors="coerce")
    return ranks


def _is_id_header(cell: str) -> bool:
    return cell in _ID_HEADERS_EXACT or cell.lower() in _ID_HEADERS_ANY_CASE


def _metadata_rows(path: Path) -> tuple[list[str], list[list[str]], list[str] | None]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = [[cell.strip() for cell in row] for row in csv.reader(handle, delimiter="\t")]
    rows = [row for row in rows if any(row)]
    while rows and rows[0][0].startswith("#") and not _is_id_header(rows[0][0]):
        rows = rows[1:]
    if not rows or not _is_id_header(rows[0][0]):
        msg = f"metadata={str(path)!r} has no QIIME 2 ID header (e.g. 'sample-id', 'id', '#SampleID')"
        raise ValueError(msg)
    header, body = rows[0], rows[1:]
    types = next((_padded_row(row, header, path) for row in body if row[0] == "#q2:types"), None)
    data = [_padded_row(row, header, path) for row in body if not row[0].startswith("#")]
    return header, data, types


def _padded_row(row: list[str], header: list[str], path: Path) -> list[str]:
    if len(row) > len(header):
        # A spreadsheet export can add empty trailing cells; QIIME 2 drops those and only
        # errors when an extra cell actually holds data.
        if any(cell for cell in row[len(header) :]):
            msg = f"metadata={str(path)!r} row {row[0]!r} has {len(row)} cells but the header has {len(header)}"
            raise ValueError(msg)
        return row[: len(header)]
    return row + [""] * (len(header) - len(row))


def _typed(values: pd.Series, declared: str, path: Path) -> pd.Series:
    values = values.mask(values == "", np.nan)
    if declared == "categorical":
        return values
    numeric = pd.to_numeric(values, errors="coerce")
    bad = values[values.notna() & numeric.isna()].unique().tolist()
    if declared == "numeric" and bad:
        msg = f"metadata={str(path)!r} column {values.name!r} is declared numeric but holds {bad[:5]}"
        raise ValueError(msg)
    # QIIME 2 infers numeric when every present value parses as a number.
    return values if bad else numeric


def _metadata(path: Path) -> pd.DataFrame:
    header, data, types = _metadata_rows(path)
    ids = [row[0] for row in data]
    frame = pd.DataFrame([row[1:] for row in data], index=ids, columns=header[1:], dtype=object)
    declared = [t.lower() for t in types[1:]] if types else [""] * len(frame.columns)
    return pd.DataFrame(
        {c: _typed(frame[c], d, path) for c, d in zip(frame.columns, declared, strict=True)}, index=frame.index
    )
