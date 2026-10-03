---
type: Module
title: fn
description: Function hierarchies, aggregation along them and HUMAnN-style renormalisation over function tables; owns no reader and no download.
resource: /src/biotapy/fn/
paths: ["src/biotapy/fn/**"]
tags: [fn, function, humann]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-03T16:20:52Z }
commit: 020efbb
status: stable
---

# Responsibility

Owns `bt.fn.*`: reading a user's mapping file into an edge table
(`_hierarchy.py:load_hierarchy`), aggregating features to one level of such a
table (`_glom.py:func_glom`) and rescaling a function MuData to community
totals (`_renorm.py:renorm`). Both verbs reproduce HUMAnN 3.9
(`humann_regroup_table`, `humann_renorm_table`); the golden tests compare them.

Does NOT own: reading HUMAnN tables or building the two-modality MuData
(`io.read_humann`, [io](/modules/io.md), over `_core.make_function_mudata`);
any download (the ENZYME hierarchy is `datasets.enzyme`,
[datasets](/modules/datasets.md)).

# Entry points

- `_hierarchy.py:load_hierarchy` - local mapping file (no header, optional
  `.gz`, BOM allowed) to an edge table, in either `parent_first` or
  `child_first` layout.
- `_glom.py:func_glom` - one AnnData (a `"function"` or a `"function_by_taxon"`
  modality, or any AnnData whose `var_names` are function ids) plus an edge
  table to a samples x groups AnnData. Call it once per modality.
- `_glom.py:_pairs` - private; the membership builder: one `(feature, group)`
  row per parent a feature has at the level, else the feature itself if it is
  protected, else `UNGROUPED`. It is where the HUMAnN semantics live.
- `_renorm.py:renorm` - the whole MuData in, a copy out with both function
  modalities divided by each sample's community total.

# Invariants

- **Edge table**: columns `child`, `parent`, `level`, optional `parent_name`,
  one row per distinct `(child, parent)` pair (`_glom.py:_EDGE_COLUMNS`,
  `_glom.py:_edges_at`). A feature's ancestors at every level are rows, not
  just its direct parent, so any level is reachable from any member
  (`datasets/_enzyme.py:_ancestors`). `load_hierarchy` instead gives exactly
  the pairs the file lists, under the one `level` the caller names. `func_glom` raises `TypeError` for a non-DataFrame `hierarchy`.
- **HUMAnN 3.9 semantics** (design notes 1-2 of
  [phase-2-function](/roadmap/phase-2-function.md)): a feature counts in full
  toward every parent at the level; unmapped features sum into `UNGROUPED`
  (per taxon when stratified); `UNMAPPED`, `READS_UNMAPPED`, `UNINTEGRATED`
  pass through (`_core._function.py:PROTECTED_FEATURES`); groups are sorted by
  name; `agg="mean"` divides by members present. `_glom.py:_pairs`,
  `_glom.py:func_glom`.
- **Nothing-maps guard** only on unstratified input: no plain feature a child
  at the level raises `ValueError` with example ids (an id-format mismatch such
  as `EC:1.1.1.1` against `1.1.1.1`). A stratified modality never raises for this: a
  pathway need not have strata, and HUMAnN writes `UNGROUPED|<taxon>` there.
  `_glom.py:_pairs` (`stratified=`), `tests/fn/test_glom.py`.
- **`x_kind` after `func_glom`**: kept for a sum in which every feature has at
  most one parent; set to `"abundance"` after a mean or a many-to-many sum,
  for any input kind, so `pp.rarefy` and `tl.alpha` refuse the result.
  `_glom.py:func_glom`.
- **`renorm` divides stratified rows by the community total** of the
  `"function"` modality (specials included unless `special=False`), not by the
  stratum's own pathway total. It raises `ValueError` unless the two modalities
  hold the same samples in the same order, warns (`_core.warn_user`) on a
  zero-total sample and leaves it zero, and sets `x_kind` to `relative` or
  `cpm`. `_renorm.py:_function_modalities`, `_renorm.py:renorm`,
  `_renorm.py:_rescaled`.
- **Outputs go through `_core`'s slot helpers**: `func_glom` builds its result
  with `_core.replace_features`, `renorm` subsets with
  `_core.feature_subset`; layers, `obsm`, `obsp`, `varm`, `varp` are dropped
  ([data-model-slots](/contracts/data-model-slots.md), Propagation).
  Neither returns a view; neither mutates its input.
- **Provenance**: `func_glom` records `hierarchy.attrs["source"]` as
  `hierarchy` and, when the table has one, `attrs["license"]` as `license`,
  so the ENZYME CC BY notice travels into the user's h5ad
  (`_glom.py:func_glom`).
- **`load_hierarchy` input rules**: every cell is stripped; blank lines and
  lines starting `#` are skipped; a line with a single id, an empty cell
  before the last id, or a file with no edges raises `ValueError` naming
  `path`. `_hierarchy.py:_read_rows`.
- **No network code**: `fn` imports only `biotapy._core` (plus mudata for
  `renorm`'s type), enforced by import-linter's `tl | fn` layer
  ([module-boundaries](/contracts/module-boundaries.md)) and by
  [no-bundled-kegg](/decisions/no-bundled-kegg.md).
- Text columns of `var` use the pandas `str` dtype, never `object`
  (`_glom.py:_group_var`), so h5ad/h5mu can be written.

# Dependencies

- [core](/modules/core.md): `sum_pairs`, `replace_features`,
  `feature_subset`, `add_provenance`, `as_csr`, `warn_user`, and the function
  constants (`FUNCTION_KEY`, `BY_TAXON_KEY`, `SPECIAL_FEATURES`,
  `PROTECTED_FEATURES`, `UNGROUPED`, `RANKS`).
- `mudata` (`_renorm.py`): the `MuData` type and `MuData.update`; see
  [function-tables-as-mudata](/decisions/function-tables-as-mudata.md).

# Verification

`uv run --group test pytest tests/fn -q`. The HUMAnN goldens alone:
`uv run --group test pytest -m golden tests/fn -q`.

# Gotchas

- **Many-to-many inflates totals.** A KO in two pathways counts in both, so
  the output's row sum exceeds the input's. That is HUMAnN's behaviour and the
  reason such an output is labelled `abundance`, not its input kind.
  `_glom.py:func_glom`.
- **`mean` counts members present in the table**, not the group's size in the
  hierarchy (`np.bincount(codes)` over the pairs actually built). Same as
  HUMAnN. `_glom.py:func_glom`.
- **`READS_UNMAPPED` is unknown to HUMAnN 3.9**: `_core._function.py:PROTECTED_FEATURES`
  follows HUMAnN master, where it passes through; 3.9 sums it into
  `UNGROUPED`. The goldens come from 3.9, so no fixture exercises the
  difference.
- **`renorm` replaces `X`**: `counts` are gone afterwards, and a stratified
  `relative` row does not sum to its pathway's community value. HUMAnN's
  `--mode levelwise` is `pp.relative` on each modality, not `renorm`.
  `_renorm.py:_rescaled`.
- **`pd.concat` of hierarchies drops `attrs`** unless every input's are
  identical, so a combined edge table records `hierarchy: null` in provenance
  and no licence. Set `attrs` on the result if provenance should name it.
- **The `str` dtype is needed for h5**: anndata's writer rejects an all-NaN
  `object` column, which `name` is whenever a hierarchy names no parents
  (`load_hierarchy`'s `parent_name` is always NaN). `_hierarchy.py:load_hierarchy`,
  `_glom.py:_group_var`.
- **`special=False` can empty the `"function"` modality** (a table holding
  only `UNMAPPED`/`UNINTEGRATED`); `renorm` then raises a `ValueError` that
  says so rather than "needs community rows". `_renorm.py:renorm`.
