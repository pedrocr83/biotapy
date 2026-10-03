---
type: Phase
title: Phase 2 - Function as a first-class hierarchy (0.2)
description: HUMAnN 4, PICRUSt2 and MetaPhlAn readers; KO -> module -> pathway hierarchies with func_glom; stratified taxa-to-function links; functional redundancy.
tags: [roadmap, fn, io]
status: stable
release: "0.2"
phase_state: in-progress
effort: ~4 weeks part-time
depends_on: [/roadmap/phase-1-core.md]
paths: ["src/biotapy/fn/**", "src/biotapy/io/**", "src/biotapy/_core/**"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-03T11:14:15Z }
commit: 2b9fc24
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: humann
    resource: https://github.com/biobakery/humann
    title: HUMAnN
---

# Goal
Differentiator 1: function gets the same glom / filter / plot verbs as taxonomy,
and HUMAnN stratified output keeps its taxa-to-function link.[^spec]

# Entry criteria
- Phase 1 exit gate passed and 0.1 released.

# Design notes (resolve before task 2.1)
- **Stratified data needs two tables.** HUMAnN pathway abundances are not
  additive over taxa, so the unstratified table cannot be derived from the
  stratified one.[^humann] `read_humann` therefore returns a `MuData` with
  modalities `function` (unstratified) and `function_by_taxon` (stratified;
  `var` columns `function`, `taxon`). This pulls `mudata` into core deps in
  Phase 2 instead of Phase 4 - needs user approval.
- **Hierarchies are many-to-many** (one KO in several pathways). `func_glom`
  uses a sparse membership matrix, so a KO's abundance counts toward every
  parent. This matches HUMAnN regrouping and is stated in the docs.
- `X` for these inputs is not raw counts: readers set
  `uns["biotapy"]["x_kind"]` (see [data-model-slots](/contracts/data-model-slots.md)).

# Tasks
Each task follows [add-a-function](/playbooks/add-a-function.md). Expand into
TDD steps (superpowers:writing-plans) when the phase starts.

- [ ] **2.1 Membership aggregation kernel** - `_core/_groupby.py`:
  `membership_matrix(child: pd.Index, edges: pd.DataFrame) -> sp.csr_matrix`;
  `sum_by` from Phase 1 reused with it. Tests: one-to-one equals Phase 1
  `sum_by`; many-to-many double-counts by design; unmapped children dropped and counted in a warning.
- [ ] **2.2 `io.read_metaphlan(path) -> TreeData`** - lineage `k__|p__|...|t__`
  split into rank columns; only leaf rows kept (no double counting of parent
  clades); `x_kind="relative"`. Tests on MetaPhlAn 4 fixture.
- [ ] **2.3 `io.read_humann(genefamilies=None, pathabundance=None, pathcoverage=None) -> MuData`** -
  unstratified + stratified modalities; `UNMAPPED`/`UNINTEGRATED` kept as
  features, flagged in `var["special"]`. Tests on HUMAnN 4 fixture; roundtrip
  sums against the fixture's own totals.
- [ ] **2.4 `io.read_picrust2(dir) -> MuData`** - `pred_metagenome_unstrat`,
  `path_abun_unstrat`, contributions file as `function_by_taxon`.
- [ ] **2.5 `fn.load_hierarchy(source: str | Path, *, levels: Sequence[str]) -> pd.DataFrame`** -
  edge table `child, parent, level`; path or URL via pooch; no bundled KEGG
  ([no-bundled-kegg](/decisions/no-bundled-kegg.md)). Licensing of each open source checked and written in the docstring.
- [ ] **2.6 `fn.func_glom(adata, level, *, hierarchy) -> AnnData`** - KO -> module -> pathway, EC -> MetaCyc.
  Property test: total abundance of mapped features x membership count preserved.
- [ ] **2.7 `fn.contributions(mdata, function, *, top=None) -> pd.DataFrame`** - taxa x samples contribution to one function.
- [ ] **2.8 `fn.functional_redundancy(mdata, *, method="contributional") -> pd.DataFrame`** -
  method from Tian et al. 2020 (contributional diversity); method note in `docs/methods/`.
- [ ] **2.9 `pl.contributions(...)`** - stacked bar of 2.7; reads data, computes nothing new.
- [ ] **2.10 Tutorial** - shotgun + HUMAnN on the anchor cohort (HMP2/IBDMDB proposed), executed notebook in `docs/tutorials/`.
- [ ] **2.11 Knowledge** - `Module` concepts for `fn` and the new readers; decision concept for the MuData-in-Phase-2 change.

# Exit gate
- [ ] Tutorial 2.10 runs in CI.
- [ ] `func_glom` and renormalisation results cross-checked against
  `humann_regroup_table` / `humann_renorm_table` outputs on the fixture (golden files).
- [ ] All Phase 1 gates still green.

# Risks
- HUMAnN 4 output format drift -> fixtures pinned to a HUMAnN version written in the test file.
- Mapping-file licensing -> 2.5 blocks on a per-source licence check.

[^spec]: Python Microbiome Toolkit development report, sections Positioning and Roadmap
[^humann]: HUMAnN
