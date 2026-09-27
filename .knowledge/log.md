# Knowledge bundle log

## 2026-09-27
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
