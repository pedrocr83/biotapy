# Knowledge bundle log

## 2026-09-27
* **Update**: Task 1.14 done: `pp.rarefy(adata, *, depth=, seed=)` subsamples
  every sample to `depth` reads without replacement, via
  `skbio.stats.subsample_counts`. `depth` defaults to the smallest non-zero
  sample depth, so all-zero samples are skipped by default; samples with
  strictly fewer than `depth` reads are dropped with one `UserWarning` naming
  up to 5; a sample at exactly `depth` is kept. Features left all-zero after
  subsampling are dropped through `_core.feature_subset`, and `X` becomes
  `int64`. New `src/biotapy/pp/_rarefy.py`, `tests/pp/test_rarefy.py`,
  `tests/pp/test_rarefy_golden.py` (matches which samples `phyloseq::
  rarefy_even_depth` drops on GlobalPatterns; counts themselves are not
  compared, per r-golden-parity's rarefaction row) and
  `docs/guide/filtering.md`. `pyproject.toml`'s `untyped_calls_exclude` gains
  `skbio.stats._subsample`. [pp](modules/pp.md) and [core](modules/core.md)
  updated; ticked in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.13 done: `pp.filter_features(adata, *, min_prevalence=,
  min_total=)` and `pp.filter_samples(adata, min_depth)`, both inclusive
  thresholds, `ValueError` when nothing passes. `filter_features` goes
  through `_core.feature_subset`; `filter_samples` subsets with AnnData
  indexing so every slot (including `obsp` distances) survives. New
  `src/biotapy/pp/_filter.py`, `tests/pp/test_filter.py`,
  `tests/pp/test_filter_golden.py` (matches `phyloseq::filter_taxa` on
  GlobalPatterns exactly) and `docs/guide/filtering.md`. [pp](modules/pp.md)
  and [core](modules/core.md) updated; ticked in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.13 step 1: fixed a pre-existing `_core.feature_subset`
  bug found while prototyping slice 1C - anndata 0.13 lists `X` itself as
  `layers[None]`, so deleting every `layers` key deleted `X` too;
  `feature_subset` now skips the `None` key. New test
  `tests/core/test_slots.py:test_feature_subset_keeps_x`. Gotcha added to
  [core](modules/core.md).
* **Update**: Task 1.15a done: added CRAN `picante` to the `biotapy-golden`
  image (own commit); `export_golden.R` gains the slice 1C block (filtering,
  rarefaction, alpha/Faith PD, Bray-Curtis/Jaccard, UniFrac on GlobalPatterns
  and esophagus, PCoA, NMDS, PERMANOVA); 15 new golden `.csv.gz` files under
  `tests/golden/global_patterns/` and `tests/golden/esophagus/`,
  `tests/golden/VERSIONS.txt` gains `vegan`, `ape` and `picante` lines. Two
  container runs are bit-identical and the 14 pre-existing golden/fixture
  files are unchanged. [r-golden-parity](contracts/r-golden-parity.md)
  statement 1 and [regenerate-golden-files](playbooks/regenerate-golden-files.md)
  (glob and two Common mistakes) updated to match; ticked in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Checkpoint B closed: the user reviewed slice 1B (PR #6 merged at their request) and approved the slice 1C plan as written; last box ticked in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Expanded slice 1C of [phase-1-core](roadmap/phase-1-core.md) into TDD steps (rules.md R1.2a): tasks 1.15a, 1.13, 1.14, 1.15b, 1.15c, 1.15, 1.16, 1.17 and Checkpoint C, prototyped against R goldens; scikit-learn>=1.8 (runtime) and picante (R image) approved. Awaiting user approval.
* **Update**: PR #6 CI green on all 19 checks, including the first GitHub run of the `network` job (6 passed); Checkpoint B CI box ticked in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 6: golden files hold derived numbers only, which may be complete for a dataset; the user accepted `relative.csv.gz` (all GlobalPatterns proportions, AGPL-3 source) for this BSD-3 repo.
* **Update**: Checkpoint B whole-branch review (Opus) of slice 1B stage 2 and its fix pass are done (F1-F8 plus a narrowed rdata warning filter); ticked in [phase-1-core](roadmap/phase-1-core.md). The PR's CI (network/golden job) and the user's review remain.
* **Update** (Checkpoint B fix pass, F1-F8): fixed 8 controller-ruled findings
  from the whole-branch review of stage 2. `_rdata.py:_refseq_error` now
  takes an optional `cause` so a rdata `NotImplementedError` unrelated to
  `refseq` (unknown format, an unsupported version, RAW/WEAKREF/DOT, an
  unknown ALTREP class) is no longer blamed on it (F1).
  `_rdata.py:read_matrix_rds`'s STR branch names `argument` instead of
  raising a raw `KeyError: 'dimnames'` for a character vector or a
  dimnames-less matrix, and wraps rdata's own `NotImplementedError` the same
  way, without mentioning refseq (F3). `_rdata.py:load_phyloseq` tells
  `.rds` from `.RData`/`.rda` by content (`parsed.object.info.type`), not
  suffix, fixing a false "holds no phyloseq object" plus two spurious
  warnings for a mismatched suffix, and now lists the Python types found and
  rejects `name=` on a single-object file naming `name=` (F4).
  `_phyloseq.py:_text_columns` casts only `pd.StringDtype`/`object` columns,
  not everything `pd.api.types.is_string_dtype` accepts (also true for a
  pandas 3.0.6 Categorical), so an R factor's order, `ordered` flag and
  unused levels survive (F2). `_phyloseq.py:_counts` casts an integer
  `otu_table` to `int64`, matching `read_dada2` (F5). The upper-case-`.RDS`
  tests now assert only `UserWarning`, not "no warning at all" (F6); a new
  test checks `read_phyloseq`'s tree against `bt.datasets.toy()`'s by
  pairwise tip-to-tip path length (F7). New fixtures
  `tests/data/phyloseq/ordered_factor.rds` (integer `otu_table`, an ordered
  factor with an unused level - `phyloseq::sample_data()` silently drops
  unused levels, so the factor is set directly on the S4 slot) and
  `tests/data/dada2/char_vector.rds` (a plain character vector, no
  dim/dimnames), generated by the unchanged `biotapy-golden` image
  (bit-identical two-run check; every pre-existing golden file and fixture
  kept its sha256 hash). Updated [io](modules/io.md) (description,
  `make_treedata` invariant, 5 new Gotchas) and [core](modules/core.md)
  (`tree_from_phylo`'s `argument` parameter); re-checked
  [datasets](modules/datasets.md) against `_remote.py`, nothing false found.
  Fixed `docs/guide/reading_data.md` (the "Example datasets" section now
  follows the general `x_kind` paragraph rather than containing it; pooch's
  cache hit does no network access at all, confirmed against the installed
  pooch's `download_action`/`fetch`, not just "checks it is still there").
  Fixed `_remote.py`'s `R equivalent:` lines: ``phyloseq::data(...)`` is not
  valid R; now ``utils::data``, with the real call
  (``data(GlobalPatterns, package = "phyloseq")``) in the docstring body.
  [phase-1-core](roadmap/phase-1-core.md): added amended-2026-09-27 pointers
  at Task 1.10 and stage-2 review focus item 3 (the historical warn-and-skip
  plan is superseded by the `ValueError` in
  [phyloseq-import-route](decisions/phyloseq-import-route.md)'s Amendment,
  without rewriting the plan text); ticked Checkpoint B's Knowledge box only.
  [regenerate-golden-files](playbooks/regenerate-golden-files.md) no longer
  claims the Dockerfile pins phyloseq's version - it pins the Bioconductor
  release (3.22); phyloseq's resulting version (1.54.2) is only recorded in
  `VERSIONS.txt`. `tests/r/export_golden.R`'s usage comment now includes
  `-e HOME=/tmp`, matching the playbook.
* **Creation**: Task 1.12b done: `tests/pp/test_transform_golden.py` and
  `tests/pp/test_glom_golden.py` compare `bt.pp.relative` and
  `bt.pp.tax_glom` (phylum, genus) against the phyloseq golden files from
  Task 1.12a, on `bt.datasets.global_patterns()`. All 3 passed first time
  against the existing golden files, no source change. Ticked Task 1.12b's
  steps in [phase-1-core](roadmap/phase-1-core.md).
* **Creation**: Task 1.11 done: `bt.datasets.global_patterns()` and
  `bt.datasets.enterotype()` (`src/biotapy/datasets/_remote.py`) download
  GlobalPatterns/enterotype from phyloseq's repository through `pooch`,
  pinned to one commit and a SHA-256, and read them with `bt.io.read_phyloseq`.
  Added `pooch` to `[project] dependencies` (approved 2026-09-27) and its
  mypy override (P8: `follow_untyped_imports`, no `py.typed` marker in
  1.9.0). Added the `network` CI job (`pytest -m "network or golden"`,
  cached `BIOTAPY_DATA_DIR`) to [phase-1-core](roadmap/phase-1-core.md)'s
  workflow, extended `tests/test_ci.py`. Updated
  [modules/datasets.md](modules/datasets.md) (entry points, invariants,
  dependencies, the `+SKIP` doctest gotcha),
  [modules/io.md](modules/io.md) (the stale "not yet written" note) and
  [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md)
  (the pooch sentence). Ticked Task 1.11's steps in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.9b fix round 1 (I1, R3.5): `_rdata.py:read_matrix_rds`
  gained a keyword-only `argument` and raises `ValueError` naming it when the
  file's converted object isn't a 2-D `xr.DataArray` - a phyloseq `.rds`
  passed as `read_dada2`'s `seqtab=`/`taxa=` used to crash with a raw
  `AttributeError`. Converting with the existing `_PHYLOSEQ` constructor dict
  instead of rdata's default avoids rdata warning once per missing S4-slot
  constructor before that check runs. Added a Gotcha to
  [modules/io.md](modules/io.md).
* **Update**: Task 1.9b fix round 1 (I2, R12.1): `modules/io.md`'s
  `_dada2.py:read_dada2` entry point said "CSV or TSV"; now names `.rds`
  (`saveRDS`) too, matching Task 1.9b.
* **Update**: Task 1.9b fix round 1 (deferred minors 4-6, addendum): an
  upper-case `.RDS` path no longer makes `read_matrix_rds` emit two spurious
  warnings (`rdata.parser.parse_file` now gets
  `extension=path.suffix.lower()`; `load_phyloseq` was checked and needs no
  change, since `rdata.read_rds` always passes the literal `extension=".rds"`
  internally, confirmed empirically before deciding not to touch it);
  `read_dada2` now casts integer `X` to `np.int64` so counts have the same
  dtype from CSV and `.rds` input. Confirmed by temporarily reverting
  `read_matrix_rds`'s `rename_axis` fix that
  `test_read_dada2_rds_matches_csv`'s `obs.index.name` assertion is the only
  one that observes it - `var.index.name` (reset unconditionally by
  `.set_axis` in `read_dada2`) and a candidate `var.columns.name` check both
  stayed `None` either way for these fixtures, so the existing
  `var.index.name` assertion was kept rather than replaced. Added two more
  Gotchas to [modules/io.md](modules/io.md).
* **Update**: Task 1.9b done: `bt.io.read_dada2`'s `seqtab` and `taxa` accept
  `.rds` files, dispatched in `_dada2.py:_read_table` to `_rdata.py:read_matrix_rds`
  (Task 1.10). Controller ruling: fixed a deferred minor from Task 1.10 in the
  same change — `read_matrix_rds` now drops xarray's anonymous `dim_0`/`dim_1`
  axis names (`.rename_axis(index=None, columns=None)`) so a `.rds`-read
  frame's index/columns match a CSV-read one's `None` names; no other change
  needed (`.rds` counts already arrive as a NumPy integer dtype, and R `NA` in
  a character matrix already arrives as NaN, both confirmed against the
  fixtures from 1.12a). Ticked Task 1.9b step checkboxes in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.10 done: `bt.io.read_phyloseq` (`src/biotapy/io/_phyloseq.py`,
  `src/biotapy/io/_rdata.py`) reads a phyloseq object saved from R
  (`.rds`/`.RData`) natively through `rdata`, with `name=` selecting one
  object from an `.RData` holding several; `_core.tree_from_phylo`
  (`_tree.py`) builds a tree from ape's `phylo` edge matrix, reusing
  `tree_from_newick`'s collision-free internal-node naming. Added runtime
  deps `rdata` (`>=1.1,<2`) and `xarray` (approved 2026-09-26; recorded in
  [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md)).
  Controller ruling (2026-09-27): a populated `refseq` slot raises
  `ValueError` naming the R fix instead of the plan's warned skip, because
  Task 1.12a's `with_refseq.rds` fixture showed rdata 1.1.0 cannot parse the
  file at all in that case; amended
  [phyloseq-import-route](decisions/phyloseq-import-route.md) with a dated
  `# Amendment 2026-09-27` section recording the evidence. Updated
  [tree-access](contracts/tree-access.md) (statement 2 and a `tree_from_phylo`
  gotcha), [modules/core.md](modules/core.md) and [modules/io.md](modules/io.md)
  (entry points, dependencies, the refseq/NULL-sentinel gotchas), and
  `docs/guide/reading_data.md`/`docs/api.md`. Ticked Task 1.10's steps in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.12a done: added the pinned R golden-file image
  (`tests/r/Dockerfile`, `rocker/r-ver:4.5.3` + Bioconductor 3.22 +
  phyloseq/Biostrings), `tests/r/export_golden.R` writing the GlobalPatterns
  golden files (`relative`, `tax_glom_phylum`, `tax_glom_genus`, gzip CSV per
  the pyarrow decline below) and R-only phyloseq/DADA2 fixtures built from
  biotapy's `toy()` numbers, and `tests/test_data_files.py` enforcing the
  1 MB cap (R6.6). The plan's Dockerfile needed a controller-approved fix
  mid-task: `rocker/r-ver:4.5.3` ships no `zlib.h`/curl/xml2/ssl/glpk/gmp
  development headers, so `BiocManager::install()` failed to compile
  `XVector` -> `Biostrings` -> `phyloseq` from source yet still exited 0
  (`docker build` looked green on a broken image); added an `apt-get`
  system-library layer and a `stopifnot(requireNamespace(...))` build-time
  check so a broken install now fails the build. Two docker builds against
  the fixed Dockerfile were needed in total (first attempt failed without
  the fix, second succeeded); the bit-identical rerun check passed on the
  first attempt with the fix. Added the
  [regenerate-golden-files](playbooks/regenerate-golden-files.md) playbook
  and its line in [playbooks/index.md](playbooks/index.md);
  [r-golden-parity](contracts/r-golden-parity.md) now states golden tests
  carry both `golden` and `network` (they run in the network CI job), that
  the image installs only what current golden files need, and that golden
  files are gzip CSV, not parquet (pyarrow declined 2026-09-27). Ticked Task
  1.12a's steps in [phase-1-core](roadmap/phase-1-core.md).
* **Verification**: `human:pedrocr83` approved the slice 1B stage-2 plan, the runtime dependency pooch and a local Docker build/run for the R golden image; declined pyarrow, so golden files are gzip CSV ([phase-1-core](roadmap/phase-1-core.md)).
* **Update**: [phase-1-core](roadmap/phase-1-core.md): Checkpoint B1 closed (PR #5 merged as c64a138, CI green incl. Python 3.14; user go-ahead). Slice 1B stage 2 expanded into TDD steps for approval: 1.12a (R container, golden files, R-written fixtures), 1.10 (read_phyloseq), 1.9b (DADA2 .rds), 1.11 (datasets via pooch + network CI job), 1.12b (golden tests); no third-party data committed.
* **Update** (Checkpoint B1 residual N1): [data-model-slots](contracts/data-model-slots.md)
  convention 1 names the rank dtype as `pd.StringDtype(na_value=np.nan)`:
  `astype("str")` means that dtype only on pandas 3, and on pandas 2.3 it
  turned missing ranks into the text `"<NA>"`.
* **Creation** (Checkpoint B1 knowledge step): [io](modules/io.md) Module
  concept for `src/biotapy/io/`; added its line to
  [modules/index.md](modules/index.md).
* **Update** (Checkpoint B1 knowledge step): [core](modules/core.md) gains
  entries for `tree_from_newick`, `tree_tips`, `relabel_tips`,
  `split_lineage`, `normalize_ranks`, `infer_x_kind` and `warn_user`,
  `make_treedata`'s id/alignment invariants, the rank-column `str` dtype, and
  scikit-bio as a `_tree.py` dependency with its measured import cost
  (~0.5s of `import biotapy`'s ~1.0-1.1s; a lazy import waits for Task
  1.21's asv measurement, R10.1).
* **Update** (Checkpoint B1 knowledge step, reviewer M4):
  [data-model-slots](contracts/data-model-slots.md) adds `confidence`
  (float, from QIIME 2 taxonomy) to the `var` keys, and states honestly that
  `pp.tax_glom`'s emptied lower ranks are float64 NaN, not the `str` dtype
  readers emit - a deferred inconsistency, not a second convention.
* **Update** (Checkpoint B1 knowledge step, reviewer M5):
  [module-boundaries](contracts/module-boundaries.md) rule 1 now allows a
  private import between topic files of the *same* subpackage (`io/_qiime2.py`
  importing `_biom_parts` from `io/_biom.py`), never across subpackages.
* **Update**: refreshed `commit` to 43d6efb on
  [tree-access](contracts/tree-access.md) and
  [function-shape](contracts/function-shape.md) after re-checking both
  against the `io` code with no content change needed.
* **Update**: [phase-1-core](roadmap/phase-1-core.md): Task 1.9's Interfaces
  line notes the Checkpoint B1 signature change (keyword-only options,
  sequence-only tree tips) over its historical plan code block; Task 1.10's
  `refseq` bullet now matches the approved
  [phyloseq-import-route](decisions/phyloseq-import-route.md) decision
  (warned and skipped, not written to `var["sequence"]`); Stage 2's heading
  gains a one-line note that `.rds` input for `read_dada2` is planned there;
  ticked Checkpoint B1 items 1-2 (review/fix pass and knowledge done).
* **Correction**: the `disallow_untyped_calls = false` entry below (reviewer
  M8) is superseded: Task 1.7c's fix round replaced it with
  `untyped_calls_exclude = ["biom"]`, scoping the strict-mode exemption to
  calls into biom-format itself rather than every call in `biotapy.io._biom`.
* **Update** (Checkpoint B1 fix M2): [tree-access](contracts/tree-access.md)
  gains a gotcha: malformed Newick raises `UnrecognizedFormatError` or
  `NewickFormatError` from scikit-bio, and `tree_from_newick` turns exactly
  those two into a `ValueError` naming its `argument`.
* **Update** (Checkpoint B1 fixes I5, I6): [tree-access](contracts/tree-access.md)
  statement 2 adds `tree_tips`; a new gotcha records that `relabel_tips`
  raises instead of merging nodes when a new name already exists, and that
  `read_dada2` tree tips must be sequences. `read_dada2` is now
  `read_dada2(seqtab, *, taxa=None, tree=None)`.
* **Update** (Checkpoint B1 fix I4): [data-model-slots](contracts/data-model-slots.md)
  convention 2 now states that readers infer `x_kind` with
  `_core.infer_x_kind` (counts, relative within `1e-3`, else abundance)
  instead of always writing `counts`.
* **Update** (Checkpoint B1 fix C1): [data-model-slots](contracts/data-model-slots.md)
  convention 1 now states that rank columns use the pandas `str` dtype, so
  reader output with an all-missing rank saves to h5td/h5ad.
* **Update**: Task 1.9 done: added `bt.io.read_dada2(seqtab, taxa=None, *,
  tree=None) -> TreeData`, the last reader in slice 1B stage 1, reading a
  DADA2 sequence table (CSV or TSV, samples x sequences) as R's `write.csv`
  writes it, with optional `assignTaxonomy`/`addSpecies` taxonomy and an
  optional Newick tree. A transposed table (columns not DNA sequences) raises
  `ValueError` naming the argument. Sequences become `ASV1..n` ids with the
  original sequence kept in `var["sequence"]`; taxonomy reuses
  `normalize_ranks`; a tree's sequence-named tips are relabeled to the
  matching ASV ids by the new `_core.relabel_tips(tree, names)` (thin wrapper
  over `nx.relabel_nodes`) before `make_treedata` attaches it. `_core.relabel_tips`
  is exported and added to
  [tree-access](contracts/tree-access.md) statement 2. Updated
  `docs/guide/reading_data.md` (new "DADA2" section) and `docs/api.md`;
  ticked Task 1.9's steps in [phase-1-core](roadmap/phase-1-core.md).

## 2026-09-26
* **Update**: Task 1.8 done: added `bt.io.read_qiime2(table, *, taxonomy=None,
  tree=None, metadata=None) -> TreeData`, reading QIIME 2 `.qza` artifacts
  (feature table, taxonomy, tree) and a QIIME 2 sample-metadata TSV without a
  QIIME 2 install. A `.qza` is read as a plain zip: `_payload` extracts
  `<uuid>/data/<filename>` and raises `ValueError` naming the argument and
  the expected payload if it is missing, so artifacts are recognized by
  content rather than by `metadata.yaml` (no YAML dependency, per stage 1
  design). The feature table reuses `_biom_parts` (Task 1.7c); taxonomy reuses
  `split_lineage` and adds an optional `confidence` column from `Confidence`;
  the tree reuses `tree_from_newick`. The metadata TSV parser recognizes both
  modern (case-insensitive `sample-id`/`id`/...) and legacy (case-sensitive
  `#SampleID`/...) ID headers, skips leading `#` comments and blank rows,
  applies an optional `#q2:types` row, and otherwise infers numeric columns
  the same way QIIME 2 does (every present value parses as a number).
  `mypy --strict` passed with no annotation changes needed. Updated
  `docs/guide/reading_data.md` (new "QIIME 2" section) and `docs/api.md`;
  ticked Task 1.8's steps in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.7d done: added `bt.io.write_biom(adata, path, *,
  fmt="hdf5") -> None`, writing `X`, rank columns and `obs` back out as a BIOM
  2.1 HDF5 or BIOM 1.0 JSON table. Rank columns are written as prefixed
  values (`k__` .. `g__`), including bare prefixes for missing ranks, because
  BIOM's HDF5 reader drops empty taxonomy entries on read and a bare prefix
  stays truthy so `read_biom` can still map it back to its rank. No tree,
  layers or embeddings are written; sample metadata is written as text.
  Confirmed against the installed biom-format source (`biom/table.py`
  `general_formatter`/`vlen_list_of_str_parser`) that HDF5 accepts str sample
  metadata and keeps `"g__"` on a round trip. Updated
  `docs/guide/reading_data.md` and `docs/api.md`; ticked Task 1.7d's steps in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.7c done: added the `io` subpackage and
  `bt.io.read_biom(path, *, tree=None) -> TreeData`, the first file reader,
  reading both BIOM dialects (JSON 1.0, HDF5 2.1) through biom-format.
  `_biom_parts` (private, reused by Task 1.8) transposes the matrix once so
  samples are rows, casts ids to `str`, turns `taxonomy` metadata into rank
  columns via `split_lineage`, and puts sample metadata straight into `obs`.
  biom-format itself rejects duplicate ids; `make_treedata` handles the
  tree/table alignment. Added runtime dependency biom-format (`>=2.1.16`) and
  a mypy override for it, plus `disallow_untyped_calls = false` scoped to
  `biotapy.io._biom` (biom-format's plain functions resolve to concrete
  untyped defs, unlike treedata/skbio, so strict flags every call into them).
  Added `docs/guide/reading_data.md` and an "Input and output" API block.
  Updated [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md);
  ticked Task 1.7c's steps in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.7b done: `_core._taxonomy` gained `normalize_ranks`
  (canonical lowercase rank columns, `domain` -> `kingdom`, `k__`/`D_0__`
  prefixes stripped, `""`/whitespace/`"NA"`/bare-prefix values -> NaN) and
  `split_lineage` (splits `;`-separated Greengenes/RESCRIPT/SILVA/`D_n__`/
  unprefixed lineages into rank columns by prefix or position, columns run to
  the deepest rank seen). Both exported from `_core`. No contract edit:
  [data-model-slots](contracts/data-model-slots.md) convention 1 already
  described this behaviour. Ticked Task 1.7b's steps in
  [phase-1-core](roadmap/phase-1-core.md).
* **Verification**: `human:pedrocr83` approved [phyloseq-import-route](decisions/phyloseq-import-route.md) (native `rdata` route, `refseq` warned and skipped) and the runtime deps `rdata` + `xarray` for stage 2; decision now `stable`, index line synced, spike files deleted.
* **Update**: Task 1.7a done: `_core._tree.tree_from_newick` parses Newick via
  scikit-bio (`convert_underscores=False`, unique internal names, NaN for
  missing lengths); `make_treedata` now casts obs/var ids to `str`, raises on
  duplicates, and aligns a tree with the table (keep the shared features, one
  `UserWarning` naming both counts, `ValueError` if nothing is shared). Added
  runtime dependency scikit-bio (`>=0.7.4,<0.8`) and a mypy override for it.
  Updated [tree-access](contracts/tree-access.md),
  [data-model-slots](contracts/data-model-slots.md) (new convention 5) and
  [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md);
  ticked Task 1.7a's steps in [phase-1-core](roadmap/phase-1-core.md).
* **Creation**: [phyloseq-import-route](decisions/phyloseq-import-route.md)
  (Task 1.6 spike, draft): `rdata.read_rda` + a ~35-line `constructor_dict`
  reads GlobalPatterns/enterotype/esophagus with zero shape mismatches and
  zero residual warnings; recommends the native `rdata` route over an R
  export script, with `refseq`/`XStringSet` deferred (no example file has a
  populated `refseq` to test against). Spike code kept in
  `.superpowers/sdd/phase-1-core/spike/` pending user approval (Step 6).
* **Verification**: `human:pedrocr83` approved the slice 1B stage-1 plan, runtime deps scikit-bio and biom-format, keeping Python 3.14 (biom-format built from source until wheels ship), and a throwaway `rdata` env for the 1.6 spike; recorded in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: [phase-1-core](roadmap/phase-1-core.md): slice 1B stage 1 (tasks 1.6-1.9) expanded into TDD steps for approval (rules.md R1.2a); `io.write_biom` renumbered 1.7b -> 1.7d; `.rds` input for DADA2 deferred to the 1.6 decision; stage 2 (1.10-1.12) expanded after it.
* **Update**: [phase-1-core](roadmap/phase-1-core.md): Checkpoint A closed; slice 1A merged to `master` (79d89bd, CI green); the user gave the go-ahead for slice 1B.
* **Creation** (Checkpoint A): [core](modules/core.md), [pp](modules/pp.md) and
  [datasets](modules/datasets.md) Module concepts for Slice 1A; added
  [modules/index.md](modules/index.md) and linked it from `index.md`.
* **Update**: [phase-1-core](roadmap/phase-1-core.md): ticked Checkpoint A's
  review and module-concepts items; fixed Task 1.1's `as_csr` interface to
  `as_csr(X: object)` (code moved in commit fc26baa, the doc had not).
  Refreshed `commit` to `0fdbd4d` on it and, after re-checking each against
  Slice 1A with no content change needed, on
  [data-model-slots](contracts/data-model-slots.md),
  [module-boundaries](contracts/module-boundaries.md),
  [tree-access](contracts/tree-access.md),
  [function-shape](contracts/function-shape.md),
  [r-golden-parity](contracts/r-golden-parity.md),
  [engine-parity](contracts/engine-parity.md),
  [add-a-function](playbooks/add-a-function.md) and
  [cut-a-release](playbooks/cut-a-release.md).
* **Update**: [phase-1-core](roadmap/phase-1-core.md): ticked Task 1.5 (`pp.tax_glom`) step checkboxes.
* **Update**: [phase-1-core](roadmap/phase-1-core.md): ticked Task 1.4 (`pp.relative`) step checkboxes.
* **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) records treedata, networkx and types-networkx (Phase 1, task 1.3).
* **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) records scipy, pandas and their type stubs (Phase 1, task 1.1).
* **Update**: [phase-1-core](roadmap/phase-1-core.md): Phase 1 runs subagent-driven (user's choice). Approved scipy and pandas (runtime) and pandas-stubs, scipy-stubs, types-networkx (dev) for `mypy --strict`. Added Task 1.7b `io.write_biom` at the user's request.
* **Update**: Phase 0 exit gate passed; [phase-0-foundation](roadmap/phase-0-foundation.md) is `done` and [phase-1-core](roadmap/phase-1-core.md) is `in-progress`. `biotapy 0.0.1` released; the release failure and fix are recorded in the Phase 0 Deviations and in [cut-a-release](playbooks/cut-a-release.md).
* **Update**: Phase 0 final-review follow-ups. Added a Deviations section to [phase-0-foundation](roadmap/phase-0-foundation.md); [cut-a-release](playbooks/cut-a-release.md) now pushes only the tag; [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) records numpy and session-info2; [maintain-knowledge](playbooks/maintain-knowledge.md) narrows its paths and warns against squash merges; [module-boundaries](contracts/module-boundaries.md) names the CI lint job. Refreshed `commit` to b77a226 on every concept with `paths` after re-checking it against the branch.
* **Creation**: Added the [cut-a-release](playbooks/cut-a-release.md) playbook (Phase 0, task 0.8).
* **Verification**: `human:pedrocr83` confirmed [package-name-biotapy](decisions/package-name-biotapy.md) (with Python >= 3.12), [pure-by-default](decisions/pure-by-default.md), [docs-okf-and-sphinx](decisions/docs-okf-and-sphinx.md) and [r-bridge-before-ports](decisions/r-bridge-before-ports.md); all now `stable`. Phase 0 dependencies and the cookiecutter-scverse v0.8.0 template approved. Phase 0 runs inline (Native); execution method to be re-asked at Phase 1.
* **Initialization**: Bundle created from the spec (`plan.md`) by `claude-code/claude-opus-5-5`: roadmap phases 0-5, six contracts, nine decisions, two playbooks. All concepts unverified; decisions marked `draft` await user confirmation.
