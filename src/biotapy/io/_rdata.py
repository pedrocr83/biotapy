"""R data files through rdata: phyloseq objects and plain matrices (.RData/.rda/.rds)."""

from collections.abc import Callable, Mapping
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pandas as pd
import rdata
import xarray as xr
from rdata.conversion import DEFAULT_CLASS_MAP, convert_attrs, convert_char, dataframe_constructor
from rdata.parser import RObjectType

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


def _refseq_error(path: Path, error: Exception) -> ValueError:
    # rdata 1.1.0's parser has no RAW branch, so a populated refseq (a Biostrings
    # DNAStringSet) fails the whole read; nothing can be skipped (ruling 2026-09-27).
    msg = (
        f"path={str(path)!r}: rdata cannot parse part of this file ({error}); a populated "
        "refseq slot (Biostrings sequences) is the usual cause. In R: "
        'Biostrings::writeXStringSet(refseq(ps), "refseq.fasta"); ps@refseq <- NULL; '
        'saveRDS(ps, "ps.rds")'
    )
    return ValueError(msg)


def load_phyloseq(path: Path, *, name: str | None) -> dict[str, Any]:
    """The phyloseq object in ``path``: an ``.rds`` file, or ``name`` (or the only one) in an ``.RData``."""
    try:
        if path.suffix.lower() == ".rds":
            objects = {"": rdata.read_rds(path, constructor_dict=_PHYLOSEQ)}
        else:
            objects = rdata.read_rda(path, constructor_dict=_PHYLOSEQ)
    except NotImplementedError as error:
        raise _refseq_error(path, error) from error
    found = {key: value for key, value in objects.items() if _is_phyloseq(value)}
    if name is not None:
        if name not in found:
            msg = f"name={name!r} is not a phyloseq object in path={str(path)!r}; found {sorted(found)}"
            raise KeyError(msg)
        return cast(dict[str, Any], found[name])
    if not found:
        msg = f"path={str(path)!r} holds no phyloseq object"
        raise ValueError(msg)
    if len(found) != 1:
        msg = f"path={str(path)!r} holds {len(found)} phyloseq objects; pass name= with one of {sorted(found)}"
        raise ValueError(msg)
    return cast(dict[str, Any], next(iter(found.values())))


def read_matrix_rds(path: Path) -> pd.DataFrame:
    """A plain R matrix saved with saveRDS, rows and columns named from its dimnames."""
    parsed = rdata.parser.parse_file(path)
    if parsed.object.info.type is RObjectType.STR:
        attrs = convert_attrs(parsed.object, lambda node: rdata.conversion.convert(node))
        flat = [convert_char(cell, default_encoding=None, force_default_encoding=False) for cell in parsed.object.value]
        return _char_matrix(flat, attrs)
    # DataArray.to_pandas() is typed to also return a Series/DataArray for other ndims;
    # an R matrix is always 2-D, so this is a DataFrame at runtime.
    return cast(pd.DataFrame, cast(xr.DataArray, rdata.conversion.convert(parsed)).to_pandas())
