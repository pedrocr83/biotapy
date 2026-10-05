---
type: Contract
title: Data-model slots
description: Which AnnData/TreeData slot holds what, the exact result keys, the x_kind and provenance conventions, and which slots feature-changing operations drop.
tags: [data-model, api]
status: stable
paths: ["src/biotapy/_core/**", "src/biotapy/io/**", "src/biotapy/pp/**", "src/biotapy/tl/**", "src/biotapy/fn/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T14:00:00Z }
commit: 0f0c9d2
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: treedata
    resource: https://github.com/YosefLab/treedata
    title: treedata source (0.3.1)
  - id: phyloseq-glom
    resource: https://github.com/joey711/phyloseq/blob/master/R/transform_filter-methods.R
    title: phyloseq tax_glom source
---

# Statement

## Slots
Extends the spec's data-model table with exact keys.[^spec]

| Slot | Holds | Keys |
|---|---|---|
| `X` | samples x features, `scipy.sparse.csr_matrix` | kind recorded in `uns["biotapy"]["x_kind"]` |
| `layers` | same-shape transforms of `X` | `relative` (sparse CSR); `clr` from `pp.clr` (dense float64 `ndarray`: CLR has no zeros) |
| `obs` | sample metadata; `tl` per-sample results with `inplace=True` | `alpha_<metric>` (e.g. `alpha_shannon`) |
| `var` | taxonomy, one lowercase column per rank; sequences; QIIME 2 assignment confidence | ranks from `kingdom, phylum, class, order, family, genus, species`; `sequence`; `confidence` (float, from a QIIME 2 `FeatureData[Taxonomy]` artifact's `Confidence` column); function tables: see Function tables |
| `vart` | phylogeny as `networkx.DiGraph`, leaves = `var_names`, edge attribute `length` | `phylo` only |
| `obsm` | ordinations and embeddings | `X_pcoa`, `X_nmds`, `X_<plugin>` |
| `obsp` | sample-sample distance matrices | metric name: `braycurtis`, `jaccard`, `unweighted_unifrac`, `weighted_unifrac` |
| `uns["biotapy"]` | biotapy metadata, nothing else | `x_kind`, `provenance`, `pcoa` (`eigenvalues`, `proportion_explained`), `nmds` (`stress`) |

`pl` reads these slots and writes none ([pl](/modules/pl.md)); `pl.contributions`
draws `fn.contributions`' table for the stratified function modality.

## Conventions
1. **Missing taxonomy** is `NaN`. Readers convert `""`, whitespace, `"NA"`, and
   bare prefixes (`"g__"`) to `NaN`, strip `k__`-style prefixes, and map rank
   aliases (`domain` -> `kingdom`) to the canonical lowercase names.
   Rank columns use `pd.StringDtype(na_value=np.nan)` - pandas 3's `str`
   dtype, spelled out because `astype("str")` on pandas 2.3 gives `object`
   with every missing value turned into the text `"<NA>"` - never
   `object`: anndata's h5ad/h5td writer rejects an all-`NaN` `object` column,
   so a reader whose `species` rank is entirely missing could not be saved
   (`_core/_taxonomy.py:normalize_ranks`). This holds for every reader's
   output; it does not (yet) hold after `pp.tax_glom`, whose emptied ranks
   below the target rank are pandas' default float64 NaN, not the `str`
   dtype (`pp/_glom.py:tax_glom`, `out.var[below] = np.nan`) - a known,
   deferred inconsistency, not a second convention.
2. **`x_kind`** is one of `counts`, `relative`, `rpk`, `cpm`, `abundance`.
   Readers always set it; most infer it from the values by `_core.infer_x_kind`
   (`_core/_slots.py`): whole numbers are `counts`; otherwise, if every
   nonzero row sums to 1 within `1e-3`, `relative`; otherwise `abundance`.
   Three readers are exceptions. `io.read_humann` reads it from the table header (`RPKs` ->
   `rpk`; `CPM`, `_cpm` or `Adjusted CPMs` -> `cpm`; `RELAB`, `_relab` ->
   `relative`) and labels a header without a unit `abundance`, never `counts`.
   `io.read_metaphlan` divides MetaPhlAn's percentages by 100 and sets `relative`, after
   checking that every sample's leaf clades sum to 1 within `1e-3` (`_core.RELATIVE_TOLERANCE`);
   a sample whose whole column is zero is exempt and stays all zero.
   `io.read_picrust2` sets `abundance`: PICRUSt2's values are read counts divided by predicted
   marker copies and multiplied by gene copies.
   `fn.renorm` rescales `X` (and may drop the special rows), setting `x_kind` to `relative` or `cpm`;
   its `relative` stratified rows do not sum to 1: they are shares of the community total.
   The other readers infer it because their formats record no unit (BIOM, QIIME 2 `RelativeFrequency`, a DADA2
   text table), and labeling proportions `counts` would misdescribe them to
   every function that reads `x_kind`. Missing key means `counts`. Functions that need raw counts
   (`pp.rarefy`; `tl.alpha` for `observed_features` and `chao1`, which
   phyloseq's `estimate_richness` refuses on non-integers; and
   `tl.unifrac(weighted=True)`) call `_core.require_counts`, which raises
   unless `x_kind` is `counts` *and* every stored value in `X` is a whole
   number (`_slots.py:require_counts`, through `infer_x_kind`'s rule, O(nnz)).
   The value check matters because fractions are otherwise truncated
   silently: `pp.rarefy` casts `X` to int64 for `subsample_counts`, and
   scikit-bio 0.7.4's tree code (Faith PD, UniFrac; `_nodes_by_counts`)
   casts abundances to int64. `tl.alpha` gives `faith_pd` presence/absence,
   so it needs no counts.
3. **Provenance** is `uns["biotapy"]["provenance"]`: a list of JSON strings
   `{"step", "version", "params"}`, appended by `_core.add_provenance`.
   JSON strings, not dicts, because h5ad cannot store a list of dicts.
4. **Trees** are created only by `_core` ([tree-access](/contracts/tree-access.md))
   with `label=None`, so TreeData adds no `tree` column to `var`.
5. **Ids are unique strings.** `_core.make_treedata` casts `obs`/`var` index
   values to `str` and raises `ValueError` on duplicates, naming them (up to
   5), so integer or mixed-type ids from a reader (e.g. unquoted BIOM JSON
   ids) never collide silently.

## Function tables
`io.read_humann` and `io.read_picrust2` return a `MuData` built by
`_core.make_function_mudata` with two modalities over the same samples, each
an `AnnData` with its own copy of `obs`:

| Modality | Features | `var` columns |
|---|---|---|
| `"function"` | community rows, e.g. `PWY-1` | `name`, `special` |
| `"function_by_taxon"` | stratified rows, e.g. `PWY-1\|g__Bacteroides.s__Bacteroides_ovatus` | `function`, `name`, `taxon`, `genus`, `species`, `special` |

`var_names` drop the row's `": name"`. `special` flags `UNMAPPED`,
`READS_UNMAPPED`, `UNINTEGRATED` and `UNGROUPED`, which stay features. Text
columns use the pandas `str` dtype, as rank columns do. Both modalities
always exist; either may have 0 features. The community values are not
the sum of their strata for pathways, which is why there are two.

In a PICRUSt2 table the stratified rows come from the long contribution
table (`taxon_function_abun`), and the modality has 0 features when no
contribution table is read. `taxon` is the ASV id as written, or `RARE`
(PICRUSt2's group of rare ASVs: an ordinary stratum, not `special`), so
`genus` and `species` are NaN. `EC:` is removed from EC numbers, in both
modalities and in `io.read_picrust2_traits`' columns, so ids match ENZYME's
and HUMAnN's.

`datasets.hmp2` returns a function table with a third modality, `"taxa"`
(MetaPhlAn species from `io.read_metaphlan`), over the same samples, and its
sample metadata in the global `obs`, pushed into every modality
(`datasets/_hmp2.py:hmp2`). `fn` verbs take the modality they need;
`fn.renorm` keeps any other modality (`fn/_renorm.py:renorm`).

## Taxonomic profiles (MetaPhlAn)
`io.read_metaphlan` keeps one feature per leaf clade: a row that no other
row descends from through any ancestor (MetaPhlAn can omit an intermediate
rank's row; the real 4.0.6 fixture has no `o__Corynebacteriales`). A clade's
row is the sum of its leaves', so every read counts once and `pp.tax_glom`
rebuilds the higher ranks. `var_names` are the leaf's last name without its
rank prefix (`SGB1871` from `t__SGB1871`); `t__` has no rank column. Rank
columns are always the seven, `kingdom` to `species`, even for a genus-level table or an all-`UNCLASSIFIED`
profile (the missing ranks are NaN, `io/_metaphlan.py:read_metaphlan`). `UNCLASSIFIED` (`UNKNOWN` in older
tables) stays a feature with every rank NaN, so every non-empty sample sums to 1.

## Propagation
| Operation | Keeps | Drops |
|---|---|---|
| Feature-changing (`pp.filter_features`, `pp.tax_glom`, `fn.func_glom`, `fn.renorm`, `pp.rarefy`) | `obs`, `var` rows kept, `vart` (pruned by TreeData), `uns["biotapy"]["x_kind"]` and `["provenance"]` | all `layers`, `obsm`, `obsp`, `varm`, `varp`, `uns["biotapy"]["pcoa"]`, `["nmds"]`, other `uns` keys |
| Sample-only (`pp.filter_samples`) | everything, subset by AnnData indexing; a kept `obsm` ordination and its `pcoa`/`nmds` summary still reflect the dropped samples, so recompute them | nothing |
| Layer-adding (`pp.relative`, `pp.clr`) | everything | nothing; adds one layer |

Feature-changing operations go through `_core.feature_subset`, or `_core.replace_features` when the new
features are groups rather than a subset; both keep only `KEPT_META`, the
single definition of the "Drops" column.

## Aggregation semantics (`tax_glom`)
Matches phyloseq:[^phyloseq-glom] features are grouped by the full lineage up
to the rank (not the rank value alone, so `uncultured` genera in different
families stay separate); the representative ("archetype") is the most abundant
feature, first on ties; ranks below the target become `NaN`. The kept tree is
the archetypes' subtree; TreeData keeps unary nodes, which leaves root-to-tip
path lengths, and therefore Faith PD and UniFrac, unchanged.[^treedata]

## Aggregation semantics (`func_glom`)
Matches `humann_regroup_table` (HUMAnN 3.9, `--ungrouped Y --protected Y`):
a feature counts in full toward every parent it has at the level; features
with none are summed into `UNGROUPED` (per taxon when `var` has `taxon`);
`UNMAPPED`, `READS_UNMAPPED` and `UNINTEGRATED` pass through
(`READS_UNMAPPED` as in HUMAnN master; 3.9 sums it into `UNGROUPED`);
`agg="mean"` divides by the members present. Groups are sorted by name.
`var` holds
`name` (from the hierarchy's `parent_name`) and `special`, plus `function`,
`taxon` and rank columns for a stratified input. The input's `x_kind` is
kept only for a sum in which every feature has at most one parent at the
level; a mean, or a sum with a feature in several parents, sets it to
`abundance`, so `require_counts` refuses it (a read would count once per
parent, and proportions would no longer sum to 1).

# Why
A slot whose meaning depends on which function wrote it cannot be trusted by
the next function. Dropping derived slots on feature changes prevents stale
distances or ordinations from being plotted against new data.

# Enforced by
- `tests/core/test_slots.py` (Phase 1, task 1.2) for `feature_subset` and provenance.
- `tests/datasets/test_toy.py` (task 1.3): h5td round-trip keeps every convention.
- Per-function tests assert the documented keys.

[^spec]: Python Microbiome Toolkit development report, section Data model
[^treedata]: treedata source (0.3.1)
[^phyloseq-glom]: phyloseq tax_glom source
