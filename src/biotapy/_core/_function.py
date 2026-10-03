"""Function tables: HUMAnN-style row ids and the two-modality MuData (contracts/data-model-slots)."""

import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData

from ._matrix import as_csr
from ._slots import XKind, add_provenance
from ._taxonomy import normalize_ranks
from ._tree import _with_str_ids

FUNCTION_KEY = "function"
BY_TAXON_KEY = "function_by_taxon"
# Rows HUMAnN writes for what it could not map or integrate; humann_renorm_table's --special list.
UNGROUPED = "UNGROUPED"
SPECIAL_FEATURES = ("UNMAPPED", "READS_UNMAPPED", "UNINTEGRATED", UNGROUPED)
# The specials humann_regroup_table passes through as themselves (master; 3.9 lacks READS_UNMAPPED).
PROTECTED_FEATURES = ("UNMAPPED", "READS_UNMAPPED", "UNINTEGRATED")
# HUMAnN strata: g__Genus.s__Species, optionally .t__SGB<id> (HUMAnN 4), or "unclassified".
_GENUS = r"(?:^|\.)g__(?P<genus>[^.]+)"
_SPECIES = r"(?:^|\.)s__(?P<species>.+?)(?:\.t__|$)"


def function_var(row_ids: "pd.Index[str]") -> pd.DataFrame:
    """Split ``ID: name|stratum`` row ids into ``var`` columns, indexed by ``ID`` or ``ID|stratum``.

    Columns: ``function`` (the id), ``name`` (NaN when the row has none),
    ``taxon`` (the stratum, NaN on community rows), ``genus`` and ``species``
    parsed from the stratum, and ``special``. A row id with more than one
    ``|`` raises ``ValueError``.
    """
    bad = row_ids[(pd.Series(row_ids, dtype=str).str.count(r"\|") > 1).to_numpy()]
    if len(bad):
        msg = f"row ids may hold one '|' (function|taxon); found {bad[:3].tolist()}"
        raise ValueError(msg)
    parts = pd.Series(row_ids, dtype=str).str.split("|")
    head = parts.str[0].str.split(": ", n=1)
    # The pandas str dtype, never object: anndata's writer rejects an all-NaN object column (contracts/data-model-slots).
    text = pd.StringDtype(na_value=np.nan)
    var = pd.DataFrame({"function": head.str[0], "name": head.str[1], "taxon": parts.str[1]}).astype(text)
    ranks = normalize_ranks(pd.concat([var["taxon"].str.extract(_GENUS), var["taxon"].str.extract(_SPECIES)], axis=1))
    var = pd.concat([var, ranks], axis=1)
    var["special"] = var["function"].isin(SPECIAL_FEATURES)
    var.index = var["function"].where(var["taxon"].isna(), var["function"] + "|" + var["taxon"]).to_numpy()
    return var


def make_function_mudata(
    X: object, *, obs: pd.DataFrame, row_ids: "pd.Index[str]", x_kind: XKind, source: str
) -> MuData:
    """Split a samples x rows function table into its ``function`` and ``function_by_taxon`` modalities.

    ``row_ids`` are HUMAnN-style (``ID: name|stratum``), one per column of
    ``X``. Ids become unique strings, as in ``make_treedata``.
    """
    obs = _with_str_ids(obs, "obs")
    var = _with_str_ids(function_var(row_ids), "var")
    matrix = as_csr(X)
    stratified = var["taxon"].notna().to_numpy()
    modalities = {
        FUNCTION_KEY: (~stratified, ["name", "special"]),
        BY_TAXON_KEY: (stratified, ["function", "name", "taxon", "genus", "species", "special"]),
    }
    mods = {}
    for key, (mask, columns) in modalities.items():
        index = np.flatnonzero(mask)
        mod = AnnData(X=matrix[:, index], obs=obs.copy(), var=var.iloc[index][columns])
        mod.uns["biotapy"] = {"x_kind": x_kind}
        add_provenance(mod, source)
        mods[key] = mod
    return MuData(mods)
