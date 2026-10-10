---
type: Module
title: io
description: File readers and writers (BIOM, and h5mu files that keep a TreeData modality's trees) for BIOM, QIIME 2 artifacts, DADA2 sequence tables, phyloseq objects and MetaPhlAn profiles (each a TreeData through _core.make_treedata), HUMAnN and PICRUSt2 tables (a MuData through _core.make_function_mudata) and PICRUSt2 per-ASV trait tables (a DataFrame).
resource: /src/biotapy/io/
paths: ["src/biotapy/io/**"]
tags: [io]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-10T23:40:56Z }
commit: bd5f929
status: stable
---

# Responsibility

Owns the `bt.io.*` verbs that move data between files and a TreeData (or, for
a function table, a MuData):
`read_biom`/`write_biom` (`_biom.py`), `read_qiime2` (`_qiime2.py`),
`read_dada2` (`_dada2.py`), `read_phyloseq` (`_phyloseq.py`, `_rdata.py`),
`read_humann` (`_humann.py`, a MuData), `read_metaphlan` (`_metaphlan.py`, a
TreeData with no tree), `read_picrust2` (`_picrust2.py`, a MuData) and
`read_picrust2_traits` (`_picrust2.py`, a DataFrame), plus two private helpers shared across
readers: the checked join (`_join.py`) and the strict TSV reading of the three
function-table readers (`_table.py`). `to_mudata` (`_mudata.py`) reads no file:
it combines tables already in memory into one MuData. Does
NOT own downloaded example datasets (`datasets.global_patterns`/`enterotype`,
[datasets](/modules/datasets.md), Task 1.11): those call `read_phyloseq`.

# Entry points

- `_biom.py:read_biom` / `_biom.py:write_biom` - BIOM 1.0 JSON and BIOM 2.1
  HDF5, both directions, through biom-format.
- `_biom.py:_biom_parts` - private; turns a loaded `biom.Table` into `(X, obs,
  var)`. Used by `read_biom` and reused by `read_qiime2` (see Gotchas).
- `_qiime2.py:read_qiime2` - a `.qza` feature table plus optional taxonomy,
  tree and QIIME 2 metadata TSV, with no QIIME 2 install.
- `_dada2.py:read_dada2` - a DADA2 `seqtab`/`seqtab.nochim` CSV, TSV or
  `.rds` (`saveRDS`), plus optional `assignTaxonomy`/`addSpecies` taxonomy
  (same three formats) and a Newick tree.
- `_phyloseq.py:read_phyloseq` - a phyloseq object saved from R
  (`.rds`/`.RData`), read natively through `rdata` (no R, no rpy2); `name=`
  selects one object from an `.RData` holding several.
- `_rdata.py:load_phyloseq` / `_rdata.py:read_matrix_rds` - private; the
  `rdata` `constructor_dict` for phyloseq's five S4 slots, and a plain R
  matrix reader (numeric or character) shared with `read_dada2`'s `.rds`
  input (Task 1.9b).
- `_humann.py:read_humann` - one HUMAnN 3 or 4 table (gene families,
  reactions, pathway abundance; per sample or joined; `.gz` allowed) to a
  two-modality MuData, `function` and `function_by_taxon`, through
  `_core.make_function_mudata`. Reads one table per call.
- `_metaphlan.py:read_metaphlan` - one MetaPhlAn 3 or 4 profile, or several
  merged, to a samples x leaf-clades TreeData; `_metaphlan.py:_sample_columns`
  picks the abundance columns.
- `_picrust2.py:read_picrust2` - an unstratified PICRUSt2 table, plus
  optionally its long contribution table (`contrib=`,
  `_picrust2.py:_contributions`), to the same two-modality MuData as
  `read_humann`. `_picrust2.py:read_picrust2_traits` - a per-ASV copy-number
  table to an ASVs x functions DataFrame. `_picrust2.py:_read` and
  `_picrust2.py:_check_first_cell` are their private header reading and
  wrong-table check.
- `_mudata.py:to_mudata` - a mapping of modality name -> AnnData to one MuData
  of the samples every modality has, in the first modality's order
  ([multiomics-as-mudata](/decisions/multiomics-as-mudata.md)).
  `_mudata.py:_check_modality` is its private per-entry check.
- `_h5mu.py:write_h5mu` / `_h5mu.py:read_h5mu` - an `.h5mu` file that keeps a
  TreeData modality's trees: mudata writes the file from `mdata.copy()`, then
  `_core.write_tree_slots` adds the trees under each TreeData modality's
  group; `read_h5mu` runs `mudata.read_h5mu` and rebuilds each marked
  modality with `_core.read_tree_slots`.
- `_table.py:_leading_lines` / `_table.py:_header` / `_table.py:_read_table` /
  `_table.py:_numbers` - private; the one strict reading of the tab-separated
  tables, shared by the three readers above (the structure in `_read_table`,
  the values in `_numbers`). Every rejection raises `ValueError` naming the
  argument.
- `_join.py:_join_to` - private; the one checked reindex-join for a side file
  (taxonomy, metadata, taxa) onto a table's ids, used by every reader above
  that joins a side file.

# Invariants

- Every TreeData reader builds its result only through `_core.make_treedata`;
  `read_humann` and `read_picrust2` build their MuData only through
  `_core.make_function_mudata`.
  `_biom.py:read_biom`, `_qiime2.py:read_qiime2` and `_dada2.py:read_dada2`
  set `x_kind` by calling `_core.infer_x_kind` on the freshly parsed matrix,
  never a literal; `_phyloseq.py:read_phyloseq` does the same.
- `read_humann` takes `x_kind` from the table header instead, because HUMAnN
  records the unit there: `RPKs` is `rpk`, `CPM`/`_cpm`/`Adjusted CPMs` is
  `cpm`, `RELAB`/`_relab` is `relative`, and a header with no unit (pathway
  abundance) is `abundance`, never `counts`, even when every value is whole.
  `_humann.py:_UNITS`, `_humann.py:read_humann`. The header is HUMAnN's rule:
  the last `#` line, or the first line when there is none
  (`_humann.py:read_humann`); the first cell can be stale after
  `humann_renorm_table --update-snames`, so a `-RELAB` or `-CPM` sample column
  outranks it (`_humann.py:_UNITS` order).
- `read_metaphlan` and `read_picrust2` also set `x_kind` themselves, never
  inferring it: `relative` after dividing percentages by 100, and
  `abundance` (PICRUSt2's values are not counts)
  (`_metaphlan.py:read_metaphlan`, `_picrust2.py:read_picrust2`).
- `read_metaphlan` keeps one feature per leaf: a clade no other clade
  descends from through *any* ancestor, not just its direct parent
  (`_metaphlan.py:read_metaphlan`). A sample whose leaves do not sum to 100%
  (`_core.RELATIVE_TOLERANCE`) raises, unless its whole column is zero; there
  is no `rank=`. It always gives all seven rank columns, whatever the deepest
  rank in the file (`_metaphlan.py:read_metaphlan`, via `_core.RANKS`).
- `read_picrust2` removes the `EC:` prefix in both modalities and in
  `read_picrust2_traits`; `RARE` is an ordinary taxon, not `special`
  (`_picrust2.py:_EC_PREFIX`, `_core/_function.py:SPECIAL_FEATURES`). It
  raises when `contrib` has no rows for a sample whose total in `path` is
  nonzero, but accepts a function subset (`_picrust2.py:_contributions`).
- Both PICRUSt2 readers check the first header cell and name the other reader
  when it is wrong: `function`, `pathway`, `#OTU ID` or `OTU ID` for
  `read_picrust2`; `sequence` for `read_picrust2_traits`
  (`_picrust2.py:_PREDICTION_FIRST`, `_picrust2.py:_TRAIT_FIRST`,
  `_picrust2.py:_check_first_cell`).
- `_table.py` reads `utf-8-sig`, so a byte-order mark never hides a leading
  `#`. It accepts an empty first (index-name) cell, and raises on: an empty
  file; repeated column names; an empty name after the first; a row with
  more cells than the header; a blank or missing id; a non-number or missing
  value; a non-finite value; a negative value when `_numbers` is called with
  `nonnegative=` (every reader does); and a decode, gzip, zlib or EOF error
  (`_table.py:_read_table`, `_table.py:_numbers`, `_table.py:_leading_lines`).
- Header rules: `read_metaphlan` and the PICRUSt2 readers use
  `_table.py:_header` (the last `#` line when it holds a tab, else the first
  line after the `#` lines). `read_humann` keeps its own, which takes the
  last `#` line whether or not it holds a tab: the shared rule would read a
  HUMAnN table's first data row as the header. They agree on every real file
  (`_humann.py:read_humann`, `_table.py:_header`).
- Sample names lose HUMAnN's suffixes (`_Abundance-RPKs`, `-CPM`, `-RELAB`, a
  joined file's `_pathabundance_cpm`) (`_humann.py:_SUFFIX`). Pathway coverage
  tables are refused (`_humann.py:_COVERAGE`).
- Every reader's error names its file argument, including the ones raised in
  `_core` (repeated sample or row ids, a row id with two `|`), which
  `read_humann`, `read_metaphlan` and `read_picrust2` re-raise with the path
  and the original as `__cause__`. `_humann.py:read_humann`,
  `_metaphlan.py:read_metaphlan`, `_picrust2.py:read_picrust2`.
- biom-format stores observations (features) x samples; `_biom_parts`
  transposes exactly once so biotapy's samples-as-rows convention holds from
  there on. `_biom.py:_biom_parts`.
- `write_biom` writes every rank as a prefixed value (`k__`..`g__`, including
  a bare prefix for a missing rank) because BIOM's HDF5 reader drops empty
  taxonomy list entries on read; a bare prefix stays non-empty so position and
  rank survive the round trip. `_biom.py:_taxonomy_metadata`.
- BIOM has no null: `write_biom` writes a missing `obs` value as `""`
  (`_biom.py:_sample_metadata`), and `read_biom` maps `""` back to NaN on the
  way in (`_biom.py:_metadata_frame`), so a missing value survives a
  `write_biom` -> `read_biom` round trip as NaN, not `""`.
- `read_qiime2` recognizes a `.qza`'s payload by the path inside the zip
  (`<uuid>/data/<filename>`), not by parsing `metadata.yaml` - avoids a YAML
  runtime dependency for a check the payload path already gives for free.
  `_qiime2.py:_payload`.
- QIIME 2 metadata ID-header matching is mixed-case on purpose: modern headers
  (`sample-id`, `id`, ...) match case-insensitively, legacy headers
  (`#SampleID`, ...) must match exactly, mirroring QIIME 2's own metadata
  parser. `_qiime2.py:_is_id_header`, `_qiime2.py:_ID_HEADERS_ANY_CASE`,
  `_qiime2.py:_ID_HEADERS_EXACT`.
- `_join.py:_join_to` is the one checked join for a side file onto a table's
  ids: a duplicate id in the side file or zero shared ids raise `ValueError`
  naming the argument; some-but-not-all table ids covered gives one
  `UserWarning` naming the count, and ids in the side file but not the table
  are silently ignored (metadata commonly covers more samples than one
  table). `_join.py:_join_to`.
- `read_dada2` requires a supplied tree's tips to be DNA sequences (checked
  against `_dada2.py:_SEQUENCE`), never bare ASV ids, before relabeling them
  to `ASV<n>`: an ASV-id tip raises naming `tree=`, because biotapy's ASV
  numbering (table column order) has no way to be verified against a
  user-supplied id. `_dada2.py:_require_sequence_tips`, `_dada2.py:read_dada2`.
- `_dada2.py:_read_table` reads the seqtab's row-name column as text verbatim
  (`dtype={0: str}`, `keep_default_na=False`, `na_values=[""]`), so a
  numeric-looking or literal `"NA"` sample name survives unparsed; rank cells
  spelling `"NA"` are still normalized to NaN downstream by
  `_core.normalize_ranks`, not by `_read_table` itself.

- `to_mudata` keeps the intersection of the modalities' sample ids
  (`Index.intersection(sort=False)`, so the first modality's order), never
  MuData's union. Dropped samples give one `UserWarning` naming how many each
  modality loses; no shared sample, an empty mapping or a repeated sample id
  raises `ValueError`. A non-string or empty name, or a non-AnnData value (a
  MuData function table included, which must go in as `{**table.mod, ...}`),
  raises `TypeError`. Each modality is `mod[shared].copy()`, so the result
  shares nothing with the inputs (rules.md R6.5) and a TreeData stays a
  TreeData. `_mudata.py:to_mudata`, `_mudata.py:_check_modality`.
- `to_mudata` writes no `uns["biotapy"]["provenance"]` entry:
  `_core.add_provenance` creates `uns["biotapy"]` with `x_kind="counts"` when
  it is missing, which would mislabel a metabolite table. Modality names are
  documented, not enforced. R equivalents: `MultiAssayExperiment::MultiAssayExperiment`
  with `MultiAssayExperiment::intersectColumns`, named in the docstring.
  `_mudata.py:to_mudata`; [phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md), design note 2.

# Dependencies

- [core](/modules/core.md): `make_function_mudata`, `XKind`, `make_treedata`, `infer_x_kind`, `split_lineage`,
  `normalize_ranks`, `RELATIVE_TOLERANCE`, `tree_from_newick`, `tree_from_phylo`, `tree_tips`,
  `relabel_tips`, `warn_user`, `TreeData`, `RANKS`, `as_csr`.
- Third-party: `biom-format` (`_biom.py`, and `_qiime2.py` via
  `_biom_parts`). No CPython 3.14 wheels yet (biocore/biom-format#1004): on
  3.14 it builds from source and needs a C compiler, noted in `README.md`.
- Third-party: `rdata` (`>=1.1,<2`) and `xarray` (`_rdata.py`), read via a
  custom `constructor_dict` - see
  [phyloseq-import-route](/decisions/phyloseq-import-route.md).

# Verification

`uv run --group test pytest tests/io -q` (`tests/io/test_mudata.py` for `to_mudata`) plus `uvx prek run --all-files`.

# Gotchas

- `MuData.write_h5mu` drops a TreeData modality's `vart` and reads it back as
  an AnnData; `io.write_h5mu` / `io.read_h5mu` keep it. Both write mudata's
  file first, so a plain mudata reader still opens it. `write_h5mu` works on
  `mdata.copy()` because mudata's writer calls `strings_to_categoricals` on
  its input (R3.3), at the cost of one more copy in memory; a non-MuData
  raises `TypeError` and a backed MuData `ValueError` (naming `mdata`) before
  any copy. String `obs`/`var` columns come back as categoricals (mudata's
  behaviour), so tests compare with `check_dtype=False`. h5mu only, not
  zarr. `read_h5mu` leaves a modality that already comes back as a TreeData
  alone (a future mudata that keeps it). `tests/io/test_mudata.py` still pins
  plain mudata's behaviour and fails the day mudata keeps the tree. Once the
  upstream PRs ship (mudata#211, treedata#103), `mudata.read_h5mu` may itself
  return TreeData modalities; `read_h5mu` already leaves those as they are. `mudata.to_mudata` exists with another meaning (it splits
  one AnnData by a column): biotapy's is always `bt.io.to_mudata`.
  [multiomics-as-mudata](/decisions/multiomics-as-mudata.md).
- `read_humann` reads the whole table into one dense rows x samples
  `float64` array before converting to CSR (about 290 MB for HMP2's 22,113 x
  1,638 pathway table); `_humann.py:read_humann`. `read_metaphlan` and
  `read_picrust2` build a dense array too (`_metaphlan.py:read_metaphlan`,
  `_picrust2.py:read_picrust2`), and `read_picrust2` holds the contribution
  table whole in memory (`_picrust2.py:_contributions`). R6.2's no-densify
  rule is met by each docstring `Notes` stating the cost, not by avoiding the
  array.
- MetaPhlAn may omit an intermediate rank's row (the 4.0.6 fixture has no
  `o__Corynebacteriales`), which is why the leaf rule looks at any ancestor
  (`_metaphlan.py:read_metaphlan`). `pp.tax_glom` drops `UNCLASSIFIED` by
  default (`pp/_glom.py:tax_glom`, `dropna=True`); pass `dropna=False` to keep
  it.
- PICRUSt2 coverage tables look like abundance tables, so `read_picrust2`
  cannot tell them apart; the docs say so, nothing detects it (slice 2B
  decision 4 in [phase-2-function](/roadmap/phase-2-function.md)). The
  Checkpoint B review's silent misreads (a `#OTU ID` table read from its
  first data row, an `EC_predicted` table read as a metagenome) are why
  `_picrust2.py:_read` uses the shared header rule and
  `_picrust2.py:_check_first_cell` exists.
- pandas renames a repeated header cell (`S1` -> `S1.1`) and, under
  `usecols`, silently drops the extra cells of a longer row, so
  `_table.py:_read_table` reads every column and compares the header's own
  cells.
- anndata 0.13.4 turns `str` columns into categoricals when writing h5mu, in
  place; the round-trip tests compare values, not dtypes
  (`tests/io/test_picrust2.py:test_round_trips_through_h5mu`).
- `read_picrust2` has no R golden; `read_metaphlan`'s parity is MetaPhlAn's
  own rows (`pp.tax_glom` equals what the profile prints,
  `tests/io/test_metaphlan.py`), see
  [r-golden-parity](/contracts/r-golden-parity.md), statement 8.
- `read_humann` has no R golden. Its R equivalent (`mia::importHUMAnN`) is named
  for the Coming-from-R table, and its parity comes from the HUMAnN goldens of
  `fn.func_glom` and `fn.renorm` that read through it
  ([r-golden-parity](/contracts/r-golden-parity.md), statement 8).
- `read_humann` reads one table per call: a gene family table and a pathway
  table both hold `UNMAPPED`, so they cannot share a modality. h5mu round-trips
  the result only because the text `var` columns are the `str` dtype
  (`_core/_function.py:function_var`).
- `io/_qiime2.py` imports the private `_biom_parts` straight from its sibling
  `io/_biom.py` (`_qiime2.py`'s `from ._biom import _biom_parts`) rather than duplicating the BIOM-parsing
  logic. This is a deliberate exception to
  [module-boundaries](/contracts/module-boundaries.md) rule 1, now amended to
  allow a private import between topic files of the *same* subpackage; it
  would still be a violation across subpackages.
- A `#q2:types` row is padded through the same `_padded_row` logic as data
  rows and zipped with `strict=True` (`_qiime2.py:_metadata_rows`): a types
  row shorter than the header used to `zip(..., strict=False)` and silently
  drop the trailing declared-but-untyped columns from the result (fixed at
  Checkpoint B1, C2).
- A metadata column declared `numeric` (`#q2:types`) raises `ValueError`
  naming the column and the bad values on a non-numeric entry; only an
  *undeclared* column falls back to text when it fails to parse as numbers
  everywhere - `pd.to_numeric(errors="coerce")` alone would have silently
  turned bad values into NaN instead of rejecting the file the way QIIME 2
  itself does. `_qiime2.py:_typed`.
- biom-format ships no `py.typed`; rather than a blanket
  `disallow_untyped_calls = false`, mypy strict mode is narrowed with
  `untyped_calls_exclude` (`[tool.mypy]` in `pyproject.toml`, which lists
  `"biom"` among a few other untyped libraries), so only calls into
  biom-format itself are exempt, not the rest of `io`.
- A populated phyloseq `refseq` slot (Biostrings sequences) raises
  `ValueError` naming the R fix, rather than being read or skipped: rdata
  1.1.0's parser has no `RAW` branch, so the whole file fails to parse
  (`NotImplementedError`) and nothing can be read around it.
  `_rdata.py:_refseq_error`, `_phyloseq.py:read_phyloseq`. See
  [phyloseq-import-route](/decisions/phyloseq-import-route.md).
- `_rdata.py:_refseq_error` takes an optional `cause`: `load_phyloseq` passes
  rdata's own text for *any* caught `NotImplementedError` (unknown file
  format, an unsupported version, RAW/WEAKREF/DOT, an unknown ALTREP class),
  so a file that fails to parse for a reason unrelated to `refseq` (a CSV
  renamed `.RData`, say) is told so rather than being blamed on a populated
  refseq slot it may not even have (Checkpoint B fix F1).
- `_rdata.py:load_phyloseq` tells an `.rds` file from an `.RData`/`.rda` one
  by content, not by suffix: `save()` writes a tagged pairlist
  (`RObjectType.LIST`) of name -> object at the top level, `saveRDS()` never
  does. This makes the reader immune to a mismatched or wrong-case suffix
  (`ps.RData` holding what `saveRDS()` wrote reads fine, no warning); `name=`
  on a file that content-detects as single-object raises `ValueError` naming
  `name=`, and the "holds no phyloseq object" message lists the Python types
  found (Checkpoint B fix F4).
- An absent phyloseq slot (`tax_table`, `sam_data`, `phy_tree`, `refseq`)
  deserializes as the literal string `"\x01NULL\x01"` with no R class, so it
  never reaches a constructor - the NULL check lives in the `phyloseq`
  constructor itself (`_rdata.py:_or_none`), not in each slot's constructor.
- `_rdata.py:read_matrix_rds` takes a keyword-only `argument` and raises
  `ValueError` naming it when the file's converted object is not a 2-D
  `xr.DataArray` (an `.rds` that holds something other than a plain matrix,
  e.g. a phyloseq object passed as `read_dada2`'s `seqtab=`). It converts
  with `_rdata.py:_PHYLOSEQ` (reused, not duplicated) rather than rdata's
  default `constructor_dict`, so a phyloseq-shaped `.rds` converts quietly to
  a `dict` instead of rdata warning once per missing S4-slot constructor
  before the `ValueError` fires (Task 1.9b fix round 1). Its STR branch (a
  character-typed `.rds`) raises the same argument-naming `ValueError`,
  instead of a raw `KeyError: 'dimnames'`, when the object has no `dim` or no
  `dimnames` attribute - a plain character vector (DADA2's `taxa=` given the
  wrong file) has neither (Checkpoint B fix F3). A `NotImplementedError` rdata
  raises while parsing or converting (e.g. a DNAStringSet `.rds`) is wrapped
  into a `ValueError` naming `argument`, with rdata's own text and no mention
  of refseq (Checkpoint B fix F3).
- `rdata.parser.parse_file` infers a file's format from `path.suffix`
  case-sensitively when no `extension=` is passed, so an upper-case `.RDS`
  path made it warn twice ("Unknown file type", "Wrong extension"). Fixed by
  passing `extension=path.suffix.lower()` in `_rdata.py:read_matrix_rds`,
  which still branches on the file's suffix (it only ever reads a `.rds`, for
  DADA2). `_rdata.py:load_phyloseq` no longer branches on suffix at all
  (previous paragraph): it parses once and hides only rdata's
  suffix-consistency `UserWarning`s, matched by their exact text
  (`_rdata.py:_SUFFIX_WARNING`), since none of them are meaningful once format
  comes from content. A filter on every `UserWarning` would also hide rdata's
  "Tag not implemented ... and ignored", which means an attribute was dropped
  (rules.md R7.4); if rdata rewords the three messages, they surface again
  instead of staying hidden (Task 1.9b fix round 1; superseded for
  `load_phyloseq` by Checkpoint B fix F4).
- `read_dada2` casts `X` to `np.int64` when it holds integers, so counts are
  the same dtype whether they came from a CSV (pandas' default `int64`) or an
  `.rds` matrix (R's 32-bit integer, `int32`) (Task 1.9b fix round 1).
  `_phyloseq.py:_counts` does the same for `read_phyloseq`, whose `otu_table`
  can arrive as R's `int32` when the R object's storage mode is integer
  (Checkpoint B fix F5).
- `_phyloseq.py:_text_columns` casts only `pd.StringDtype`/`object` columns of
  `sam_data` to biotapy's NaN-backed text dtype, not every column
  `pd.api.types.is_string_dtype` accepts - that predicate is also true for a
  Categorical of strings on pandas 3.0.6, which used to silently flatten an R
  factor (order, `ordered`, unused levels all lost) into plain text
  (Checkpoint B fix F2).
