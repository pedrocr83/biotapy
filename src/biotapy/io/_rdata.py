"""R data files through rdata: phyloseq objects and plain matrices (.RData/.rda/.rds)."""

import warnings
from collections.abc import Callable, Mapping
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pandas as pd
import rdata
import xarray as xr
from rdata.conversion import DEFAULT_CLASS_MAP, convert_attrs, convert_char, dataframe_constructor
from rdata.parser import RData, RObjectType

# An unset "...OrNULL" S4 slot is serialised as this marker string without a class,
# so no constructor sees it; the phyloseq constructor maps it to None.
_R_NULL = "\x01NULL\x01"
PHYLOSEQ_SLOTS = ("otu_table", "tax_table", "sam_data", "phy_tree", "refseq")
# rdata does not export its ConstructorDict alias; constructors take Any because callable
# parameters are contravariant, and narrow with cast.
_Constructors = Mapping[str | bytes, Callable[[Any, Mapping[str, Any]], Any]]


def _otu_table(obj: Any, attrs: Mapping[str, Any]) -> tuple[xr.DataArray, bool]:
    return cast(xr.DataArray, obj), bool(attrs["taxa_are_rows"][0])


def _char_matrix(obj: Any, attrs: Mapping[str, Any]) -> pd.DataFrame:
    # rdata leaves character matrices flat (column-major) and drops dim/dimnames.
    rows, columns = attrs["dimnames"]
    values = np.reshape(np.asarray(obj, dtype=object), attrs["dim"], order="F")
    return pd.DataFrame(values, index=rows, columns=columns)


def _passthrough(obj: Any, attrs: Mapping[str, Any]) -> Any:
    return obj


def _or_none(value: Any) -> Any:
    return None if isinstance(value, str) and value == _R_NULL else value


def _phyloseq(obj: Any, attrs: Mapping[str, Any]) -> dict[str, Any]:
    slots = cast(SimpleNamespace, obj)
    return {name: _or_none(getattr(slots, name)) for name in PHYLOSEQ_SLOTS}


_PHYLOSEQ: _Constructors = {
    **DEFAULT_CLASS_MAP,
    "otu_table": _otu_table,
    "taxonomyTable": _char_matrix,
    "sample_data": dataframe_constructor,
    "phylo": _passthrough,
    "phyloseq": _phyloseq,
}


def _is_phyloseq(value: object) -> bool:
    return isinstance(value, dict) and tuple(value) == PHYLOSEQ_SLOTS


def _refseq_error(path: str | Path, cause: str | None) -> ValueError:
    """Build the one refseq/parse-failure message, shared by both raise sites.

    ``cause`` is rdata's own text when the whole file failed to parse for some reason
    (unknown format, an unsupported version, an unimplemented node such as RAW, an
    unknown ALTREP class - rdata 1.1.0 has no RAW branch, so a populated refseq, a
    Biostrings DNAStringSet, is one such case, but never the only one a bare
    ``NotImplementedError`` could mean). ``None`` means the file parsed fine and it is
    specifically the phyloseq object's own refseq slot that is populated.
    """
    fix = 'Biostrings::writeXStringSet(refseq(ps), "refseq.fasta"); ps@refseq <- NULL; saveRDS(ps, "ps.rds")'
    if cause is None:
        msg = f"path={str(path)!r}: the refseq slot holds sequences, which rdata cannot parse; export them first: {fix}"
    else:
        msg = (
            f"path={str(path)!r}: rdata cannot parse this file ({cause}); if it is a phyloseq "
            f"object with a populated refseq slot, export the sequences first: {fix}"
        )
    return ValueError(msg)


def load_phyloseq(path: Path, *, name: str | None) -> dict[str, Any]:
    """The phyloseq object in ``path``: an ``.rds`` file, or ``name`` (or the only one) in an ``.RData``."""
    try:
        with warnings.catch_warnings():
            # extension only steers rdata's own suffix-consistency UserWarnings, never what
            # gets parsed (rdata.parser._parser.parse_data, confirmed empirically, R2.2):
            # this reader tells RDS from RDATA by content below, so a mismatched or
            # upper-case suffix must never warn (ruling 2026-09-27).
            warnings.simplefilter("ignore", UserWarning)
            parsed = rdata.parser.parse_file(path)
        converted = rdata.conversion.convert(parsed, _PHYLOSEQ)
    except NotImplementedError as error:
        raise _refseq_error(path, str(error)) from error
    # R's save() writes a tagged pairlist (RObjectType.LIST) of name -> object at the top
    # level; saveRDS() never does, even for a named R list (a VECSXP, not a pairlist) -
    # confirmed empirically against toy.rds (S4) vs toy.RData (LIST).
    is_environment = parsed.object.info.type is RObjectType.LIST
    if name is not None and not is_environment:
        msg = f"name={name!r} selects an object in an .RData/.rda file; path={str(path)!r} holds a single object"
        raise ValueError(msg)
    objects = cast(dict[str, Any], converted) if is_environment else {"": converted}
    return _select_phyloseq(path, objects, name=name)


def _select_phyloseq(path: Path, objects: Mapping[str, Any], *, name: str | None) -> dict[str, Any]:
    found = {key: value for key, value in objects.items() if _is_phyloseq(value)}
    if name is not None:
        if name not in found:
            msg = f"name={name!r} is not a phyloseq object in path={str(path)!r}; found {sorted(found)}"
            raise KeyError(msg)
        return cast(dict[str, Any], found[name])
    if not found:
        kinds = sorted({type(value).__name__ for value in objects.values()})
        msg = f"path={str(path)!r} holds no phyloseq object; found {kinds}"
        raise ValueError(msg)
    if len(found) != 1:
        msg = f"path={str(path)!r} holds {len(found)} phyloseq objects; pass name= with one of {sorted(found)}"
        raise ValueError(msg)
    return cast(dict[str, Any], next(iter(found.values())))


def read_matrix_rds(path: Path, *, argument: str) -> pd.DataFrame:
    """A plain R matrix saved with saveRDS, rows and columns named from its dimnames.

    ``argument`` names the input in the ``ValueError`` raised when the file holds
    something other than a matrix, e.g. ``"seqtab='seqtab.rds'"``.
    """
    try:
        # rdata infers the extension from path.suffix by default, case-sensitively, and
        # warns twice for e.g. ".RDS"; passing it lower-cased avoids that false positive.
        parsed = rdata.parser.parse_file(path, extension=path.suffix.lower())
        frame = _matrix_frame(parsed, argument)
    except NotImplementedError as error:
        msg = f"{argument} must be a matrix saved with saveRDS; rdata cannot parse this file ({error})"
        raise ValueError(msg) from error
    # xarray names a plain DataArray's axes dim_0/dim_1; a plain R matrix has no axis names.
    return frame.rename_axis(index=None, columns=None)


def _matrix_frame(parsed: RData, argument: str) -> pd.DataFrame:
    if parsed.object.info.type is RObjectType.STR:
        attrs = convert_attrs(parsed.object, lambda node: rdata.conversion.convert(node))
        if attrs.get("dim") is None or "dimnames" not in attrs or len(attrs["dim"]) != 2:
            msg = f"{argument} must be a matrix saved with saveRDS; the file holds a character vector or a matrix without dimnames"
            raise ValueError(msg)
        flat = [convert_char(cell, default_encoding=None, force_default_encoding=False) for cell in parsed.object.value]
        return _char_matrix(flat, attrs)
    # _PHYLOSEQ's constructors are reused (not duplicated) so a phyloseq-shaped .rds
    # converts quietly to a dict instead of rdata warning about missing constructors
    # before the isinstance check below rejects it as not a matrix.
    obj = rdata.conversion.convert(parsed, _PHYLOSEQ)
    if not isinstance(obj, xr.DataArray) or obj.ndim != 2:
        msg = f"{argument} must be a matrix saved with saveRDS; the file holds {type(obj).__name__}"
        raise ValueError(msg)
    # DataArray.to_pandas() is typed to also return a Series/DataArray for other ndims;
    # ndim == 2 is checked above, so this is a DataFrame at runtime.
    return cast(pd.DataFrame, obj.to_pandas())
