---
type: Phase
title: Phase 2 - Function as a first-class hierarchy (0.2)
description: "HUMAnN 3/4, PICRUSt2 and MetaPhlAn readers; user-supplied and ENZYME hierarchies with func_glom; HUMAnN-parity renorm; stratified taxa-to-function links; functional redundancy (Tian 2020)."
tags: [roadmap, fn, io]
status: stable
release: "0.2"
phase_state: in-progress
effort: ~4 weeks part-time
depends_on: [/roadmap/phase-1-core.md]
paths: ["src/biotapy/fn/**", "src/biotapy/io/**", "src/biotapy/_core/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-03T17:56:00Z }
commit: 731617d
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: humann
    resource: https://github.com/biobakery/humann
    title: HUMAnN
---

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Slice 2A has full TDD steps;
> slices 2B-2D are outlines, expanded with superpowers:writing-plans when
> reached (rules.md R1.2a).

**Goal:** 0.2 makes function a first-class hierarchy: HUMAnN, MetaPhlAn and
PICRUSt2 outputs read into the data model, function aggregated and
renormalised with HUMAnN's own semantics (golden-tested against HUMAnN), the
taxa-to-function link kept and plotted, and functional redundancy computed.[^spec]

**Architecture:**
- A function table is a `MuData` with two modalities over the same samples:
  `"function"` (community rows) and `"function_by_taxon"` (stratified rows).
  `_core` builds it (`make_function_mudata`), so `io` readers and `datasets`
  produce the identical layout.
- `fn` verbs are pure and take one `AnnData` modality (`func_glom`) or the
  `MuData` when both modalities are needed (`renorm`). Many-to-many
  aggregation is one sparse product in `_core` (`sum_pairs`).
- Hierarchies are plain edge tables (`child, parent, level, parent_name`).
  `fn.load_hierarchy` parses the user's local files and never touches the
  network; the only download, ENZYME (CC BY 4.0), lives in `datasets` with the
  other pooch fetches.
- Golden files come from HUMAnN 3.9 itself, run by a pinned `uv run` script;
  no Docker and no database are needed.

**Tech stack:** anndata 0.13 · mudata 0.4.1 · numpy · scipy.sparse · pandas 3 ·
pooch 1.9 · HUMAnN 3.9 (golden files only, never a dependency) ·
pytest/hypothesis.

**Spec:** [plan.md](../../plan.md). Contracts that bind every task:
[function-shape](/contracts/function-shape.md),
[data-model-slots](/contracts/data-model-slots.md),
[module-boundaries](/contracts/module-boundaries.md),
[r-golden-parity](/contracts/r-golden-parity.md). Decision:
[no-bundled-kegg](/decisions/no-bundled-kegg.md). Research (2026-10-03, in the
session scratchpad, summarised in the ledger
`.superpowers/sdd/phase-2-function/progress.md`): A HUMAnN/MetaPhlAn formats
and tools, B PICRUSt2, mapping licences and functional redundancy, C mudata,
tutorial cohort and the golden route.

**How slice 2A was checked.** Every file in slice 2A was written into a
scratch copy of the repository at 2dd8488 (master, after 0.1.0) and passed:
ruff 0.16.9 (`check`, `format --check`), `mypy --strict`, import-linter, the
full pytest run (691 passed, 22 deselected, doctests included; the property
tests also under four extra Hypothesis seeds), coverage
(99% overall; every new or changed file 100%) and `sphinx-build -W`. mudata came from an
ephemeral `uv run --with 'mudata>=0.4'` overlay (0.4.1). The golden files
were generated twice by the export script and are bit-identical. The real
HMP2 tables were read end to end (`pathabundances_3.tsv.gz`, 1,638 samples x
22,113 rows, in 3.0 s; a per-sample `level4ec` table regrouped to ENZYME
classes and renormalised). Implementers still re-check every API they call
(R2.2).

# Goal
Differentiator 1: function gets the same glom / filter / plot verbs as taxonomy,
and HUMAnN stratified output keeps its taxa-to-function link.[^spec]

# Entry criteria
- Phase 1 exit gate passed and 0.1 released.

# Design notes (resolved)

Each note answers one design question from the planning brief, with the
reason. Notes marked **(user)** change a contract, a rule or a ruling and are
repeated under "Decisions for the user".

1. **Renormalisation needs a new function, `fn.renorm(mdata, units, *,
   special=True) -> MuData` (new task 2.12).** HUMAnN's
   `humann_renorm_table` (research A section 2.2, verified by running 3.9):
   - `--units cpm|relab` (default `cpm`): totals 1,000,000 or 1.
   - The level of a row is its number of `|`-separated parts (1 =
     community, 2 = stratified). In the default `--mode community` **every
     row, stratified ones included, is divided by the sample's level-1
     (community) total**. `--mode levelwise` divides each level by its own
     total.
   - With `--special y` (default) the community total includes the special
     rows (`UNMAPPED`, `UNINTEGRATED`, `UNGROUPED`, and in master
     `READS_UNMAPPED`); `--special n` drops those rows, community and
     stratified, before totalling.
   - A zero-total sample: HUMAnN sets the total to 1 and warns, so the
     sample stays zero.
   - Output is printed with `%.6g` (six significant digits).

   `pp.relative` cannot do community mode: it divides each modality by its
   own row sums, which is exactly `--mode levelwise`. So community mode needs
   both modalities at once, hence a `MuData` verb. `renorm` replaces `X`
   (HUMAnN's tables are what users analyse, and `x_kind` exists to say what
   `X` holds), sets `x_kind` to `"relative"` or `"cpm"`, and goes through
   `_core.feature_subset`, so derived slots are dropped as after any feature
   change. `units` is a required argument because HUMAnN's default (`cpm`)
   and biotapy's habit (`relative`) differ; no default means no surprise.
   Levelwise mode gets no option: it is `bt.pp.relative` per modality, which
   the docs say (R2.3). **(user: new task)**

2. **`func_glom` matches `humann_regroup_table`'s defaults (`-u Y -p Y`).**
   Features with no parent at the level are summed into `UNGROUPED` (per
   taxon in the stratified modality, e.g. `UNGROUPED|unclassified`);
   `UNMAPPED`, `READS_UNMAPPED` and `UNINTEGRATED` map to themselves; an input
   `UNGROUPED` row is unmapped, so it joins the new `UNGROUPED`. Reasons:
   - the exit gate compares against HUMAnN's output, which has these rows;
   - dropping unmapped abundance would silently shrink the community total
     that `renorm` divides by;
   - `UNGROUPED` shows the user how much did not map. Dropping it is one
     indexing line.

   This replaces roadmap 2.1's "unmapped children dropped and counted in a
   warning". No warning on a partial mapping (a real UniRef -> EC regroup
   leaves most abundance in `UNGROUPED`, so a warning would fire on every
   call). When **no** non-special feature of an unstratified input maps,
   `func_glom` raises `ValueError` naming three features and three hierarchy
   children, which catches id-format mismatches (`EC:1.1.1.1` vs `1.1.1.1`).
   A stratified input never raises (Checkpoint A ruling): a pathway need not
   have strata, so its rows go to `UNGROUPED|<taxon>`, as HUMAnN writes them,
   and the community call already catches a mismatch. `READS_UNMAPPED`
   is protected, as in HUMAnN master; HUMAnN 3.9 and 4.0.0a2 do not know it,
   so no golden input contains it. **(user: semantics change)**

3. **`x_kind`: no new value; the reader rules change.** Values in use:

   | Source | `x_kind` | Rule |
   |---|---|---|
   | HUMAnN `RPKs` header (`_Abundance-RPKs`, `Adjusted RPKs`) | `rpk` | header |
   | HUMAnN `-CPM` / `_cpm` / `Adjusted CPMs` | `cpm` | header |
   | HUMAnN `-RELAB` / `_relab` | `relative` | header |
   | HUMAnN header without a unit (pathway abundance) | `abundance` | never inferred: whole-number pathway abundances must not become `counts` |
   | MetaPhlAn percentages (2.2) | `relative` | divided by 100 at read; `UNCLASSIFIED` kept so rows sum to 1 (2B decision) |
   | PICRUSt2 predicted metagenomes (2.4) | `abundance` | copy-number-weighted, marker-normalised read counts, never integers |
   | after `fn.renorm` | `relative` / `cpm` | stratified rows do not sum to 1, by design |

   The data-model contract's convention 2 says readers infer `x_kind` from
   the values; that sentence gains the HUMAnN header rule and the MetaPhlAn
   division. Contract wording change only. **(user)**

4. **MuData layout.**
   - Modalities: `"function"` and `"function_by_taxon"` (constants
     `_core.FUNCTION_KEY`, `_core.BY_TAXON_KEY`). Both always exist; a table
     with no stratified rows gives a `function_by_taxon` with 0 features (it
     round-trips through h5mu, verified).
   - `var_names`: the row id without its `": name"` part, so `PWY-1` and
     `PWY-1|g__Bacteroides.s__Bacteroides_ovatus`. These match the hierarchy
     children and HUMAnN's own regrouped ids, and they are unique across the
     two modalities because only stratified ids hold `|`.
   - `var` of `"function"`: `name`, `special`. `var` of
     `"function_by_taxon"`: `function`, `name`, `taxon`, `genus`, `species`,
     `special`. Text columns are pandas' `str` dtype (h5 writers reject an
     all-NaN `object` column: the `name` column of a gene-family table is all
     NaN). `genus`/`species` come from the HUMAnN stratum
     (`g__X.s__Y[.t__SGB…]`); `unclassified` and bare labels give NaN.
   - `obs`: each modality gets its own copy of the same sample index, so
     `AnnData` verbs work on a modality alone. The global `mdata.obs` starts
     with no columns (mudata 0.4 does not pull); users put metadata there and
     `push_obs()` it.
   - One table per `read_humann` call. A gene-family and a pathway table
     both hold `UNMAPPED`; in one modality they would collide (mudata warns
     `var_names are not unique`, verified).
   - The global obs, other modalities and `uns` survive `fn.renorm`.
   - Why two tables at all: a pathway's community abundance is not the sum
     of its strata.[^humann]
   **(user: contract addition)**

5. **Layering.** `fn` sits on the `tl | fn` layer; `io` and `datasets` are
   below and import only `_core`.
   - **`_core` holds the kernel and the function-table builder.**
     `sum_pairs` (many-to-many group sum) joins `sum_by` in
     `_core/_matrix.py`, not a new `_groupby.py`: same family, one topic.
     `replace_features` (new features, old samples, derived slots dropped)
     joins `feature_subset` in `_slots.py` because it needs the private
     `KEPT_META`, and the contract names `_slots.py` the one place for
     propagation rules. `make_function_mudata`, `function_var` and the
     special-row constants are used by `io` and `datasets` (and in 2B by
     `read_picrust2`), so R4.3 puts them in `_core/_function.py`.
     `sum_pairs` and `replace_features` have one consumer (`fn`) today. R4.3
     would keep them in `fn`, but roadmap 2.1 put the kernel in `_core` and
     `sum_by` set the precedent. **(user: judgement call)**
   - **The ENZYME fetch lives in `datasets`, as `bt.datasets.enzyme()`, not
     in `fn.load_hierarchy`.** `datasets` already owns pooch, the
     `BIOTAPY_DATA_DIR` cache, the offline-test pattern and the licensing
     page. Keeping every download there makes "biotapy never fetches KEGG or
     MetaCyc" a structural fact: `fn` has no network code at all.
     `load_hierarchy` becomes a pure local parser with no `source="enzyme"`
     special case. The ledger's ruling says "load_hierarchy built-in fetch =
     ENZYME"; this keeps the substance (ENZYME is the only built-in fetch)
     and moves the entry point. **(user)**
   - `pl.contributions` (2.9) may import `fn` (higher layer).

6. **Licence notices.**
   - `tests/data/humann/NOTICE.txt`: the three files copied from HUMAnN's
     test data (commit e07b3a3), with HUMAnN's MIT text, and the statement
     that the other three fixtures are synthetic. The sdist ships `tests/`,
     so the notice travels with the copies, as MIT requires. The golden CSVs
     are HUMAnN's output on biotapy's synthetic inputs and carry no
     third-party data.
   - `tests/data/enzyme/NOTICE.txt`: the ENZYME excerpt's source, release,
     CC BY 4.0 and the changes made (CC BY requires stating them).
   - `datasets.enzyme()` carries the attribution three ways: the docstring
     `Notes`, `attrs = {"source": "ENZYME release <date>, SIB,
     https://enzyme.expasy.org/", "license": "CC BY 4.0"}`, and the datasets
     guide's Licensing paragraph. `func_glom` copies `hierarchy.attrs["source"]`
     into its provenance entry, so the attribution travels with results.
   - MetaPhlAn's `demo_metaphlan_bugs_list.tsv` (2B) comes from HUMAnN's
     repository and gets its own `tests/data/metaphlan/NOTICE.txt`.
   - ENZYME keeps only its current release online (checked 2026-10-03:
     `ftp.expasy.org/databases/enzyme/` lists one `enzyme.dat`, 2026-09-02).
     A pinned SHA-256 would break for every new user at the next release, so
     the registry hash is `None` and the release date read from the file goes
     into `attrs["source"]`. That is the first unpinned download in biotapy.
     **(user)**

7. **Roadmap corrections.**
   - **2.8:** the method is Tian et al. 2020, *Nature Communications*
     11:6217 (not *Nature Ecology & Evolution*, not "contributional"): FR =
     Gini-Simpson - Rao's Q over a taxa x genes copy-number matrix with
     weighted-Jaccard taxon dissimilarity. With one method, `method=` is a
     speculative option (R2.3) and is dropped:
     `fn.functional_redundancy(adata, *, traits) -> pd.DataFrame`.
   - **2.5:** "KO -> module -> pathway" is user-supplied only (KEGG is
     restricted, research B 2.1); the built-in hierarchy is ENZYME's EC
     class tree. `load_hierarchy` reads local files only (no URL: an
     arbitrary URL is how KEGG would get fetched).
   - **2.10:** the cohort is an HMP2/IBDMDB subset: HUMAnN 3.0 alpha and
     MetaPhlAn 3 outputs in CPM, fetched by pooch with pinned SHA-256,
     never committed. ibdmdb.org states no licence, so the tutorial carries
     a caveat and the citation (Lloyd-Price et al., Nature 2019).
   - **2.3:** `read_humann(path) -> MuData` reads one table per call;
     `pathcoverage` is dropped (HUMAnN 4 no longer writes it).
   - **Frontmatter**: `description` reads "HUMAnN 3/4, PICRUSt2 and MetaPhlAn
     readers; user-supplied and ENZYME hierarchies with func_glom;
     HUMAnN-parity renorm; stratified taxa-to-function links; functional
     redundancy (Tian 2020)." HUMAnN 4 is still an alpha (4.0.0a2), so
     "HUMAnN 4" alone overstates it.

# Global constraints
- Python >= 3.12; `mudata>=0.4` (runtime, approved 2026-10-03); everything
  from Phase 1 stays.
- A function table is a `MuData` with modalities `"function"` and
  `"function_by_taxon"` (design note 4). `X` is CSR, samples x features, in
  every modality.
- **No KEGG or MetaCyc content anywhere:** not in the wheel, the sdist, the
  tests, the fixtures or the docs, and never fetched. `fn` has no network
  code; downloads live only in `datasets` (pooch). Identifiers alone (`K00001`,
  `PWY-5100`) are fine in synthetic fixtures; prefer EC numbers.
- **PICRUSt2 (GPL-3) is never installed, imported or copied**, code or test
  data. Its fixtures are synthetic, built from the documented formats.
- Third-party fixtures only under MIT (HUMAnN) or CC BY 4.0 (ENZYME), each
  directory with a `NOTICE.txt` (design note 6). Every file under 1 MB (R6.6).
- HUMAnN semantics follow 3.9's `humann_regroup_table` /
  `humann_renorm_table` defaults; `READS_UNMAPPED` is protected and special
  as in HUMAnN master.
- HUMAnN golden files: `uv run --no-project --with humann==3.9 --with
  pandas==3.0.6 python tests/humann/export_golden.py`, never in CI; two runs
  must be bit-identical (`mtime=0` gzip).
- Size limits (R5), with one trap the prototype hit: ruff's
  `max-positional-args = 3` (`PLR0917`) applies to private helpers and to
  tests too. Make the fourth argument keyword-only, or parametrize a test
  with one id and unpack a case dict.
- mypy: mudata ships no `py.typed`. It gets `follow_untyped_imports` like
  pooch, and `"mudata"` joins `untyped_calls_exclude` (only for
  `MuData.update`). `MuData.mod` is typed as a read-only `Mapping`; it is a
  `dict` subclass (`ModDict`), and the one write casts it.
- Docs are `nitpicky`: `MuData` in a signature needs the intersphinx entry
  `"mudata": ("https://mudata.scverse.org/stable/", None)` (the readthedocs
  URL 404s; checked 2026-10-03). Docstrings name later functions only as
  ``literals``, never as cross-references.
- Commits stage explicit paths only. Never stage `.claude/`, `.superpowers/`,
  `.worktrees/` or `notebooks/`. Every task's last commit also stages
  `.knowledge/roadmap/phase-2-function.md` with that task's boxes ticked and
  `.knowledge/log.md` with its dated line (R12.4).

# Dependencies to approve (ask at the start of the task named)
| Task | Group | Package | Reason |
|---|---|---|---|
| 2.1b | runtime | mudata `>=0.4` | two-modality function tables; pure Python, BSD-3; its `scverse-misc[settings]` needs are already installed through anndata - approved 2026-10-03 |
| 2.0 | none (tool run, not installed) | humann `==3.9`, pandas `==3.0.6` via `uv run --no-project --with` | golden files - ruling 2026-10-03 |
| 2.2 | none | - | mia not added to the R image (user, 2026-10-03); parity from MetaPhlAn's own rows (r-golden-parity statement 8) |
| 2.8 | none | - | weighted Jaccard comes from SciPy's `braycurtis` (identity checked) |
| 2.10 | none | - | the cohort loader uses pooch, already a dependency |

# Review focus
The five ways real users are most likely to get a wrong answer from slice 2A
without an error. Each line names the test that pins it.

1. **The table's ids and the hierarchy's ids differ in form** (PICRUSt2
   `EC:1.1.1.1` vs ENZYME `1.1.1.1`; KEGG REST `ko:K00001` vs HUMAnN
   `K00001`). Expected: `func_glom` raises, naming three ids from each side,
   and never returns an all-`UNGROUPED` table. Test: 2.6
   `test_no_feature_in_the_hierarchy_raises_with_examples`.
2. **Units mislabelled.** Cases: a renormalised HUMAnN table, an HMP2 table
   whose columns are joined file names (`<id>_pathabundance_cpm`, no `#`
   header), a pathway table of whole numbers. Expected: `x_kind` is `cpm`,
   `relative` or `abundance` as the header says, and never `counts`, so
   `pp.rarefy` refuses HUMAnN tables. Tests: 2.3
   `test_x_kind_comes_from_the_header`,
   `test_renormalised_column_names_win_over_the_first_cell` and
   `test_without_a_hash_line_the_first_line_is_the_header`.
3. **Saving fails on real data.** A gene-family table has no names, so
   `var["name"]` is all NaN, and the h5 writers reject an all-NaN `object`
   column. Expected: reader, toy and `func_glom` output round-trip through
   h5mu/h5ad. Tests: 2.1b `test_make_function_mudata_round_trips_through_h5mu`,
   2.3 `test_round_trips_through_h5mu`, 2.6 `test_output_round_trips_through_h5ad`.
4. **Strata scaled on a different base than their community rows**
   (`pp.relative` per modality). Expected: `fn.renorm` divides stratified rows
   by the community total, so a pathway's strata keep their share. Tests:
   2.12 `test_strata_are_divided_by_the_community_total` and the
   `renorm_*` golden tests.
5. **HUMAnN 4 alpha output.** Its header reads `# Gene Family HUMAnN v4…
   Adjusted CPMs`, the unmapped row is `READS_UNMAPPED`, strata end in
   `.t__SGB…` and some strata are bare labels. Expected: read, flagged and
   parsed like 3.x. Tests: 2.1b
   `test_function_var_parses_genus_and_species_from_the_stratum`, 2.3
   `test_reads_humann_4_style_tables`.

# Slices
The roadmap's order (kernel and hierarchy first, readers second) is changed.
The exit gate's golden cross-check runs through the reader. `func_glom`'s
stratified behaviour is defined by the reader's `var` layout. `renorm`, which
the exit gate also requires, needs the `MuData`. Building function tables by
hand in 2A's tests would duplicate the reader and test the wrong thing.

| Slice | Delivers | Tasks | Ends with |
|---|---|---|---|
| **2A - HUMAnN path** | read, regroup and renormalise HUMAnN tables, equal to HUMAnN 3.9 | 2.0 goldens and fixtures · 2.1 `sum_pairs`, `replace_features` · 2.1b function tables in `_core` (+ mudata) · 2.3 `io.read_humann` · 2.3b `datasets.toy_humann` · 2.5a `datasets.enzyme` · 2.5 `fn.load_hierarchy` · 2.6 `fn.func_glom` · 2.12 `fn.renorm` | Checkpoint A |
| **2B - Other readers** | MetaPhlAn and PICRUSt2 into the same data model | 2.2a shared table reading · 2.2 `io.read_metaphlan` · 2.4 `io.read_picrust2` · 2.4b `io.read_picrust2_traits` | Checkpoint B |
| **2C - Analysis** | the taxa-to-function link, redundancy, the plot | 2.7 `fn.contributions` · 2.8 `fn.functional_redundancy` · 2.9 `pl.contributions` | Checkpoint C |
| **2D - Cohort, docs, release** | the HMP2 tutorial in CI, knowledge, benchmarks, 0.2 | 2.10 `datasets.hmp2` + tutorial · 2.11 knowledge · 2.13 Coming-from-R check · 2.14 benchmarks · 2.15 release 0.2 | exit gate |

Execution order inside 2B: **2.2a -> 2.2 -> 2.4 -> 2.4b -> Checkpoint B.**

Execution order inside 2A: **2.0 -> 2.1 -> 2.1b -> 2.3 -> 2.3b -> 2.5a -> 2.5
-> 2.6 -> 2.12 -> Checkpoint A.** The golden export (2.0) needs no biotapy
code and every later golden test reads its files. ENZYME (2.5a) comes before
`load_hierarchy` (2.5) so the function guide never names a function that does
not exist yet.

# Tasks (checklist)
- [x] 2.0 HUMAnN fixtures, notice and golden export
- [x] 2.1 `_core.sum_pairs` and `_core.replace_features`
- [x] 2.1b `_core` function tables and the mudata dependency
- [x] 2.3 `io.read_humann(path) -> MuData`
- [x] 2.3b `datasets.toy_humann() -> MuData`
- [x] 2.5a `datasets.enzyme() -> pd.DataFrame`
- [x] 2.5 `fn.load_hierarchy(path, level, *, layout="parent_first") -> pd.DataFrame`
- [x] 2.6 `fn.func_glom(adata, level, *, hierarchy, agg="sum") -> AnnData`
- [x] 2.12 `fn.renorm(mdata, units, *, special=True) -> MuData`
- [ ] Checkpoint A
- [x] 2.2a io/_table.py: strict table reading shared by the readers
- [x] 2.2 io.read_metaphlan(path) -> TreeData
- [ ] 2.4 io.read_picrust2(path, *, contrib=None) -> MuData and 2.4b io.read_picrust2_traits(path) -> pd.DataFrame (2.4 done, 2.4b pending)
- [ ] Checkpoint B
- [ ] 2.7 `fn.contributions(mdata, function, *, top=None) -> pd.DataFrame` (outline)
- [ ] 2.8 `fn.functional_redundancy(adata, *, traits) -> pd.DataFrame` (outline)
- [ ] 2.9 `pl.contributions(...) -> Axes` (outline)
- [ ] Checkpoint C
- [ ] 2.10 `datasets.hmp2()` and the tutorial (outline)
- [ ] 2.11 Knowledge (outline)
- [ ] 2.13-2.15 Coming-from-R check, benchmarks, release 0.2 (outline)

# Exit gate
- [ ] Tutorial 2.10 runs in CI (docs job, pooch cache).
- [ ] `fn.func_glom` and `fn.renorm` equal `humann_regroup_table` /
  `humann_renorm_table` (HUMAnN 3.9) on the fixtures: `tests/fn/*_golden.py`.
- [ ] All Phase 1 gates still green.

# Risks
- HUMAnN 4 is an alpha; its format may drift -> fixtures and goldens pin
  HUMAnN commit e07b3a3 / version 3.9, written in `NOTICE.txt` and
  `VERSIONS.txt`.
- ENZYME has no versioned archive -> unpinned download; the release read is
  recorded in `attrs` and provenance (design note 6).
- PICRUSt2 file layouts are taken from its code and wiki, never from a run
  -> 2B's synthetic fixtures may differ from a real file; the first real
  user file is the check (ruling 2026-10-03).
- HMP2 hosting is a single Globus endpoint with no stated licence -> pooch
  cache in CI; tutorial caveat.

---
## Slice 2A - HUMAnN path

**Goal:** a HUMAnN user reads any HUMAnN 3 or 4 table, regroups it along a
hierarchy and renormalises it, and gets HUMAnN 3.9's numbers.

### Slice 2A design
- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `_core/_matrix.py` | `sum_pairs`: many-to-many column sums (2.1) |
  | `_core/_slots.py` | `replace_features`: new features, same samples, derived slots dropped (2.1) |
  | `_core/_function.py` | `FUNCTION_KEY`, `BY_TAXON_KEY`, `SPECIAL_FEATURES`, `PROTECTED_FEATURES`, `function_var`, `make_function_mudata` (2.1b) |
  | `io/_humann.py` | `read_humann` (2.3) |
  | `datasets/_toy.py` | `toy_humann` beside `toy` (2.3b) |
  | `datasets/_enzyme.py` | `enzyme`; its two files join `_remote.py`'s pooch registry (2.5a) |
  | `fn/_hierarchy.py` | `load_hierarchy` (2.5) |
  | `fn/_glom.py` | `func_glom` (2.6) |
  | `fn/_renorm.py` | `renorm` (2.12) |
  | `tests/humann/export_golden.py` | the HUMAnN 3.9 golden export (2.0) |

- **The hierarchy edge table** is the one interface between hierarchy
  sources and `func_glom`. Columns are `child`, `parent`, `level` (the
  parent's level) and `parent_name` (optional; it becomes `var["name"]`), one
  row per distinct pair. Every ancestor is a row, not only the direct parent:
  `enzyme()` gives `1.1.1.1` three rows (class, subclass, sub-subclass), so a
  table already grouped to sub-subclasses still reaches classes. `attrs["source"]`
  names where it came from and is copied into `func_glom`'s provenance.
- **Facts the tasks rely on** (measured on the prototype; re-check each, R2.2).
  - HUMAnN 3.9 (`uv run --no-project --with humann==3.9`) runs
    `humann_regroup_table -c map.tsv` and `humann_renorm_table` with no
    database (contrary to research A's reading of the source; the ledger's
    controller check confirms it). 4.0.0a2 gave byte-identical output on the
    research tables.
  - On the synthetic gene families, regroup puts `UniRef90_D` and
    `UniRef90_unknown` into `UNGROUPED` (5.0 in S1), with
    `UNGROUPED|unclassified` from `UniRef90_D|unclassified`.
    `UniRef90_unknown` is not special, so it lands there too. A group whose
    members are all absent (`G3`) gets no row. `-f mean` divides by the
    members present (`G2` = (6 + 4) / 2) and applies to `UNGROUPED` too
    (2.5 = 5 / 2).
  - `humann_renorm_table -u cpm` on the gene families divides `UNMAPPED` 10 by
    the community total 33 (all six community rows, specials included):
    303030 at `%.6g`.
  - pandas' default C float parser keeps about 15 significant digits
    (5.1e-14 relative error measured); `float_precision="round_trip"` is exact
    but takes 5.8 s instead of 2.1 s on the 91 MB HMP2 table. The reader keeps
    the default, like the Phase 1 readers.
  - The HMP2 merged table `pathabundances_3.tsv.gz` has **no `#` header**
    (first cell `Feature\Sample`) and sample columns such as
    `CSM5FZ3N_P_pathabundance_cpm`, with whole-number CPM values. Its
    per-sample tables use `# Pathway<TAB><id>_Abundance-CPM`.
  - mudata 0.4.1: `MuData({...})` keeps references to the modalities; global
    `obs` starts with no columns; a 0-feature modality and h5mu round trips
    work; `mdata.copy()` followed by replacing a modality and `update()` keeps
    the global `obs`.
  - anndata's h5 writer raises `TypeError: Can't implicitly convert non-string
    objects to strings` on an all-NaN `object` column (hit by `func_glom`'s
    `var["name"]` in the prototype), and `ValueError: DataFrame.index.name
    ('function') is also used by a column` when `var.index` keeps the name of
    the `function` column (hit by `function_var`). Both are fixed in the code
    below and pinned by tests.
  - `renorm` divides each stored value by its sample's total instead of
    multiplying by `1 / total`. Hypothesis found a subnormal total (5e-324)
    whose reciprocal overflows to `inf`. (`pp.relative` multiplies by the
    reciprocal; it is not changed here, R1.4. Reported for a later fix.)

### Slice 2A global constraints (in addition to the Phase 2 list)
- Golden tests in `tests/fn/*_golden.py` carry `pytest.mark.golden` only. Their
  inputs are committed fixtures, so they need no network and run in every CI
  job.
- Tests reach functions through `bt.<module>.<fn>`. `_core` unit tests and
  test fixtures built with `_core.make_function_mudata` are the exceptions
  (Phase 1 precedent: `from biotapy._core import get_tree` in `tests/pp`).
  The ENZYME tests monkeypatch `datasets._enzyme._fetch`, the same documented
  exception as `tests/datasets/test_remote.py` (R11.4).
- Run commands: `uv run --group test pytest <path> -q`. Gate before every
  commit: `uvx prek run --all-files`; with docs changes also
  `uv run --group doc sphinx-build -W -b html docs docs/_build/html`.

### Slice 2A review focus
The Phase 2 review focus (above) is slice 2A's. In addition:
- **Purity:** `func_glom` and `renorm` never mutate their input. `renorm`
  copies the `MuData`, and `feature_subset` returns real objects. Tests:
  2.6 and 2.12 `test_input_unchanged`, 2.1 `test_replace_features_does_not_touch_its_input`.
- **Empty modalities** (a table with no strata, or a stratified-only split):
  `read_humann` gives a 0-feature modality, `func_glom` returns 0 groups, and
  `renorm` raises naming the missing community rows. Tests: 2.3
  `test_community_rows_with_no_strata_give_an_empty_modality`, 2.6
  `test_empty_modality_gives_no_groups`, 2.12 `test_stratified_only_table_raises`.

---

### Task 2.0: HUMAnN fixtures, notice and golden export

**Files:** create `tests/data/humann/{multi_sample_genefamilies.tsv,
gene_families.tsv, demo_pathabundance_with_names.tsv}` (copied, MIT),
`tests/data/humann/{genefamilies.tsv, pathabundance.tsv, regroup_map.tsv}`
(synthetic), `tests/data/humann/NOTICE.txt`, `tests/humann/export_golden.py`,
`tests/golden/humann/{regroup_sum, regroup_mean, renorm_relab, renorm_cpm,
renorm_relab_nospecial, renorm_cpm_genefamilies}.csv.gz`,
`tests/golden/humann/VERSIONS.txt`; modify `pyproject.toml` (sdist exclude),
`.knowledge/contracts/r-golden-parity.md`,
`.knowledge/playbooks/regenerate-golden-files.md`.
**Interfaces (produces):** the fixtures and golden files that 2.3, 2.5, 2.6
and 2.12 read. Golden CSV layout: index `sample_id` (`S1`, `S2`, `S3`),
columns = HUMAnN row ids without `": name"` (`UNMAPPED`, `G1`,
`G1|g__Bacteroides.s__Bacteroides_ovatus`, ...).

- [x] **Step 1: Copy HUMAnN's test files at a pinned commit**
  ```bash
  B=https://raw.githubusercontent.com/biobakery/humann/e07b3a3/humann/tests/data
  mkdir -p tests/data/humann
  for f in multi_sample_genefamilies.tsv gene_families.tsv demo_pathabundance_with_names.tsv; do
    curl -sSf -o "tests/data/humann/$f" "$B/$f"
  done
  # Upstream ends these two without a newline; pre-commit's end-of-file-fixer would add it anyway.
  printf '\n' >> tests/data/humann/gene_families.tsv
  printf '\n' >> tests/data/humann/demo_pathabundance_with_names.tsv
  sha256sum tests/data/humann/*.tsv
  ```
  Expected:
  ```text
  ccee707fbb0fdb3843b7c4cd5251841aad6396b7ca1e8c116b20681cb0c34b11  tests/data/humann/demo_pathabundance_with_names.tsv
  abafba07bca02b87f1fc9b656e6194114972199369384c4eef5662be49faa1f7  tests/data/humann/gene_families.tsv
  8a53688414541a78f8137a75c9956fca0a2ee4b2f32a8e2f08193350ce301e76  tests/data/humann/multi_sample_genefamilies.tsv
  ```
- [x] **Step 2: Write the synthetic inputs** (tab-separated; every line of a
  table has the same number of fields). They cover what the copied files lack:
  `UNMAPPED`, `UNINTEGRATED` with strata, a named feature, an unmapped
  feature with a stratum, `UniRef90_unknown`, a many-to-many map with a member
  absent from the table, an all-zero sample (`S3`) and a pathway whose strata
  do not sum to its community value (`PWY-1`: 15 + 2 < 20).

  `tests/data/humann/genefamilies.tsv`:
  ```text
  # Gene Family	S1_Abundance-RPKs	S2_Abundance-RPKs	S3_Abundance-RPKs
  UNMAPPED	10.0	20.0	0.0
  UniRef90_A: alpha protein	8.0	2.0	0.0
  UniRef90_A: alpha protein|g__Bacteroides.s__Bacteroides_ovatus	5.0	2.0	0.0
  UniRef90_A: alpha protein|unclassified	3.0	0.0	0.0
  UniRef90_B	6.0	4.0	0.0
  UniRef90_B|g__Bacteroides.s__Bacteroides_ovatus	6.0	1.0	0.0
  UniRef90_B|g__Blautia.s__Blautia_obeum	0.0	3.0	0.0
  UniRef90_C	4.0	0.5	0.0
  UniRef90_C|g__Blautia.s__Blautia_obeum	4.0	0.5	0.0
  UniRef90_D	2.0	6.0	0.0
  UniRef90_D|unclassified	2.0	6.0	0.0
  UniRef90_unknown	3.0	1.0	0.0
  ```
  `tests/data/humann/pathabundance.tsv`:
  ```text
  # Pathway	S1_Abundance	S2_Abundance	S3_Abundance
  UNMAPPED	40.0	10.0	0.0
  UNINTEGRATED	30.0	25.0	0.0
  UNINTEGRATED|g__Bacteroides.s__Bacteroides_ovatus	20.0	5.0	0.0
  UNINTEGRATED|unclassified	5.0	15.0	0.0
  PWY-1: first pathway	20.0	10.0	0.0
  PWY-1: first pathway|g__Bacteroides.s__Bacteroides_ovatus	15.0	4.0	0.0
  PWY-1: first pathway|g__Blautia.s__Blautia_obeum	2.0	5.0	0.0
  PWY-2: second pathway	10.0	5.0	0.0
  PWY-2: second pathway|g__Blautia.s__Blautia_obeum	10.0	5.0	0.0
  ```
  `tests/data/humann/regroup_map.tsv` (`humann_regroup_table --custom` format:
  group, then members; `G2` repeats; `UniRef90_E` is not in the table):
  ```text
  G1	UniRef90_A	UniRef90_B
  G2	UniRef90_B	UniRef90_C
  G2	UniRef90_E
  G3	UniRef90_E
  ```
  Check: `awk -F'\t' '{print FILENAME, NF}' tests/data/humann/genefamilies.tsv tests/data/humann/pathabundance.tsv | sort -u`
  Expected: `tests/data/humann/genefamilies.tsv 4` and `tests/data/humann/pathabundance.tsv 4`.
- [x] **Step 3: Write the notice** `tests/data/humann/NOTICE.txt`:
  ```text
  multi_sample_genefamilies.tsv, gene_families.tsv and
  demo_pathabundance_with_names.tsv are copied from HUMAnN's test data,
  https://github.com/biobakery/humann/tree/e07b3a3/humann/tests/data
  (commit e07b3a3, 2026-07-10), unchanged except for a final newline added to the
  last two, under HUMAnN's MIT licence:

  The HUMAnN software is licensed under the MIT license.

  Copyright (c) 2014 Harvard School of Public Health

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.

  genefamilies.tsv, pathabundance.tsv and regroup_map.tsv are synthetic, written
  for biotapy (BSD-3-Clause); they follow the formats in HUMAnN's documentation.
  ```
- [x] **Step 4: Write the export script** `tests/humann/export_golden.py`:
  ```python
  """Write the HUMAnN golden files (contracts/r-golden-parity).

  Run from the repository root, never in CI:

      uv run --no-project --with humann==3.9 --with pandas python tests/humann/export_golden.py

  HUMAnN's utility scripts need no database for a custom mapping file. Each run
  rewrites tests/golden/humann/*.csv.gz (samples as rows, row ids without names)
  and tests/golden/humann/VERSIONS.txt.
  """

  import subprocess
  import tempfile
  from importlib.metadata import version
  from pathlib import Path

  import pandas as pd

  ROOT = Path(__file__).resolve().parents[2]
  DATA = ROOT / "tests" / "data" / "humann"
  GOLDEN = ROOT / "tests" / "golden" / "humann"
  RUNS = {
      "regroup_sum": ("humann_regroup_table", "genefamilies.tsv", ["-c", str(DATA / "regroup_map.tsv")]),
      "regroup_mean": ("humann_regroup_table", "genefamilies.tsv", ["-c", str(DATA / "regroup_map.tsv"), "-f", "mean"]),
      "renorm_relab": ("humann_renorm_table", "pathabundance.tsv", ["-u", "relab"]),
      "renorm_cpm": ("humann_renorm_table", "pathabundance.tsv", ["-u", "cpm"]),
      "renorm_relab_nospecial": ("humann_renorm_table", "pathabundance.tsv", ["-u", "relab", "-s", "n"]),
      "renorm_cpm_genefamilies": ("humann_renorm_table", "genefamilies.tsv", ["-u", "cpm"]),
  }


  def to_golden(path: Path) -> pd.DataFrame:
      """A HUMAnN table as samples x row ids, names dropped (``ID: name|taxon`` -> ``ID|taxon``)."""
      table = pd.read_csv(path, sep="\t", index_col=0)
      ids = table.index.str.split("|", n=1)
      table.index = [parts[0].split(": ", 1)[0] + ("|" + parts[1] if len(parts) == 2 else "") for parts in ids]
      table.columns = table.columns.str.replace(r"_Abundance(-RPKs)?$", "", regex=True)
      return table.T.rename_axis("sample_id")


  def main() -> None:
      GOLDEN.mkdir(parents=True, exist_ok=True)
      with tempfile.TemporaryDirectory() as tmp:
          for name, (tool, source, options) in RUNS.items():
              out = Path(tmp) / f"{name}.tsv"
              subprocess.run([tool, "-i", str(DATA / source), "-o", str(out), *options], check=True)
              # mtime=0: two runs write bit-identical files (playbooks/regenerate-golden-files).
              to_golden(out).to_csv(GOLDEN / f"{name}.csv.gz", compression={"method": "gzip", "mtime": 0})
      (GOLDEN / "VERSIONS.txt").write_text(f"humann {version('humann')}\npandas {version('pandas')}\n")


  if __name__ == "__main__":
      main()
  ```
- [x] **Step 5: Run it twice; the second run must be bit-identical**
  ```bash
  uv run --no-project --with humann==3.9 --with pandas==3.0.6 python tests/humann/export_golden.py
  sha256sum tests/golden/humann/* > /tmp/humann-golden.sha
  uv run --no-project --with humann==3.9 --with pandas==3.0.6 python tests/humann/export_golden.py
  sha256sum -c /tmp/humann-golden.sha
  ```
  Expected: the first run builds HUMAnN's 107 MB sdist once (about 1.5 min),
  prints HUMAnN's log lines and `WARNING: Column 3 (S3_Abundance) has zero sum
  at level 1`, and the check reports `OK` for all 7 files.
  `VERSIONS.txt` reads `humann 3.9` / `pandas 3.0.6`. The CSVs
  (`zcat tests/golden/humann/<name>.csv.gz`) must read exactly:

  `regroup_sum`:
  ```text
  sample_id,UNMAPPED,UNGROUPED,UNGROUPED|unclassified,G1,G1|g__Bacteroides.s__Bacteroides_ovatus,G1|g__Blautia.s__Blautia_obeum,G1|unclassified,G2,G2|g__Bacteroides.s__Bacteroides_ovatus,G2|g__Blautia.s__Blautia_obeum
  S1,10.0,5.0,2.0,14.0,11.0,0.0,3.0,10.0,6.0,4.0
  S2,20.0,7.0,6.0,6.0,3.0,3.0,0.0,4.5,1.0,3.5
  S3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
  ```
  `regroup_mean`:
  ```text
  sample_id,UNMAPPED,UNGROUPED,UNGROUPED|unclassified,G1,G1|g__Bacteroides.s__Bacteroides_ovatus,G1|g__Blautia.s__Blautia_obeum,G1|unclassified,G2,G2|g__Bacteroides.s__Bacteroides_ovatus,G2|g__Blautia.s__Blautia_obeum
  S1,10.0,2.5,2.0,7.0,5.5,0.0,3.0,5.0,6.0,2.0
  S2,20.0,3.5,6.0,3.0,1.5,3.0,0.0,2.25,1.0,1.75
  S3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
  ```
  `renorm_relab`:
  ```text
  sample_id,UNMAPPED,UNINTEGRATED,UNINTEGRATED|g__Bacteroides.s__Bacteroides_ovatus,UNINTEGRATED|unclassified,PWY-1,PWY-1|g__Bacteroides.s__Bacteroides_ovatus,PWY-1|g__Blautia.s__Blautia_obeum,PWY-2,PWY-2|g__Blautia.s__Blautia_obeum
  S1,0.4,0.3,0.2,0.05,0.2,0.15,0.02,0.1,0.1
  S2,0.2,0.5,0.1,0.3,0.2,0.08,0.1,0.1,0.1
  S3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
  ```
  `renorm_cpm`:
  ```text
  sample_id,UNMAPPED,UNINTEGRATED,UNINTEGRATED|g__Bacteroides.s__Bacteroides_ovatus,UNINTEGRATED|unclassified,PWY-1,PWY-1|g__Bacteroides.s__Bacteroides_ovatus,PWY-1|g__Blautia.s__Blautia_obeum,PWY-2,PWY-2|g__Blautia.s__Blautia_obeum
  S1,400000,300000,200000,50000,200000,150000,20000,100000,100000
  S2,200000,500000,100000,300000,200000,80000,100000,100000,100000
  S3,0,0,0,0,0,0,0,0,0
  ```
  `renorm_relab_nospecial`:
  ```text
  sample_id,PWY-1,PWY-1|g__Bacteroides.s__Bacteroides_ovatus,PWY-1|g__Blautia.s__Blautia_obeum,PWY-2,PWY-2|g__Blautia.s__Blautia_obeum
  S1,0.666667,0.5,0.0666667,0.333333,0.333333
  S2,0.666667,0.266667,0.333333,0.333333,0.333333
  S3,0.0,0.0,0.0,0.0,0.0
  ```
  `renorm_cpm_genefamilies`:
  ```text
  sample_id,UNMAPPED,UniRef90_A,UniRef90_A|g__Bacteroides.s__Bacteroides_ovatus,UniRef90_A|unclassified,UniRef90_B,UniRef90_B|g__Bacteroides.s__Bacteroides_ovatus,UniRef90_B|g__Blautia.s__Blautia_obeum,UniRef90_C,UniRef90_C|g__Blautia.s__Blautia_obeum,UniRef90_D,UniRef90_D|unclassified,UniRef90_unknown
  S1,303030.0,242424.0,151515.0,90909.1,181818.0,181818.0,0.0,121212.0,121212.0,60606.1,60606.1,90909.1
  S2,597015.0,59701.5,59701.5,0.0,119403.0,29850.7,89552.2,14925.4,14925.4,179104.0,179104.0,29850.7
  S3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
  ```
- [x] **Step 6: Keep the script out of the sdist**, like `tests/r`
  (`pyproject.toml`, `[tool.hatch] build.targets.sdist.exclude`): add
  `"/tests/humann",` after `"/tests/r",`. The fixtures and `NOTICE.txt` stay
  in the sdist, which is what MIT asks of copies.
- [x] **Step 7: Knowledge** (contract wording, user-approved; R12.1).
  - `.knowledge/contracts/r-golden-parity.md`:
    - `description`: append "; HUMAnN-parity functions are tested the same
      way against files exported from pinned HUMAnN 3.9".
    - after statement 1, add statement 1b:
      ````markdown
      1b. HUMAnN golden files (`fn.func_glom`, `fn.renorm`) are produced by
         `tests/humann/export_golden.py`, run with
         `uv run --no-project --with humann==3.9 --with pandas==3.0.6`, never
         in CI. HUMAnN's utility scripts are pure Python and need no
         database, so there is no container. Output:
         `tests/golden/humann/<name>.csv.gz` (samples as rows, row ids without
         their `": name"`) and `tests/golden/humann/VERSIONS.txt`.
      ````
    - comparison table, new row: `| HUMAnN parity (func_glom, renorm) |
      elementwise, matched by row id | rtol=1e-7; renorm rtol=5e-6, because
      humann_renorm_table prints %.6g |`.
    - statement 6, last sentence becomes: "Test fixtures under `tests/data/`
      stay synthetic, except small files copied under a permissive licence
      with a `NOTICE.txt` beside them (`tests/data/humann`: HUMAnN's MIT test
      data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0)."
    - "Enforced by", new bullet: "`tests/fn/*_golden.py`, marker `golden`
      only: their inputs are committed, so they run in every CI job."
  - `.knowledge/playbooks/regenerate-golden-files.md`: add
    `"tests/humann/**"` and `"tests/data/humann/**"` to `paths`, and a section
    before `# Verification`:
    ````markdown
    # HUMAnN golden files
    When HUMAnN's pin (3.9) changes, or a HUMAnN-parity golden is added:
    ```bash
    uv run --no-project --with humann==3.9 --with pandas==3.0.6 python tests/humann/export_golden.py
    sha256sum tests/golden/humann/* > /tmp/humann-golden.sha
    uv run --no-project --with humann==3.9 --with pandas==3.0.6 python tests/humann/export_golden.py
    sha256sum -c /tmp/humann-golden.sha
    ```
    The script writes gzip with `mtime=0`, so reruns are bit-identical. The
    synthetic inputs in `tests/data/humann/` are hand-written: change them
    deliberately, then regenerate; never edit a golden CSV.
    ````
- [x] **Step 8: Run** `uv run --group test pytest tests/test_data_files.py tests/test_knowledge_bundle.py -q`
  Expected: all pass (every new file is under 1 MB).
- [x] **Step 9: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add tests/data/humann tests/humann/export_golden.py tests/golden/humann pyproject.toml \
    .knowledge/contracts/r-golden-parity.md .knowledge/playbooks/regenerate-golden-files.md \
    .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "test(fn): add HUMAnN fixtures and HUMAnN 3.9 golden files"
  ```

### Task 2.1: `_core.sum_pairs` and `_core.replace_features`

**Files:** modify `src/biotapy/_core/_matrix.py`, `src/biotapy/_core/_slots.py`,
`src/biotapy/_core/__init__.py`, `tests/core/test_matrix.py`,
`tests/core/test_slots.py`.
**Interfaces (produces):**
- `sum_pairs(X: sp.csr_matrix, features: npt.NDArray[np.intp], groups: npt.NDArray[np.intp], *, n_groups: int) -> sp.csr_matrix`:
  column sums of `X` into `n_groups` groups from membership pairs; a feature
  in several groups counts fully in each; repeated pairs count once; dtype
  kept.
- `replace_features(adata: AnnData, X: sp.csr_matrix, var: pd.DataFrame) -> AnnData`:
  new `AnnData` with `adata`'s `obs` (copied), the given `X` and `var`, and
  `uns = {"biotapy": {x_kind, provenance}}` from `adata`; `layers`, `obsm`,
  `obsp`, `varm`, `varp` and other `uns` keys are not carried.

- [x] **Step 1: Failing tests.** Append to `tests/core/test_matrix.py` and
  change its import to `from biotapy._core import argmax_by, as_csr, sum_by, sum_pairs`:
  ```python
  def test_sum_pairs_counts_a_feature_in_every_group():
      # f1 belongs to groups 0 and 1, so it counts in full toward both.
      out = sum_pairs(X, np.array([0, 1, 1, 2]), np.array([0, 0, 1, 1]), n_groups=2)
      np.testing.assert_array_equal(out.toarray(), [[3, 5], [9, 11]])


  def test_sum_pairs_counts_a_repeated_pair_once():
      out = sum_pairs(X, np.array([0, 0]), np.array([0, 0]), n_groups=1)
      np.testing.assert_array_equal(out.toarray(), [[1], [4]])


  def test_sum_pairs_with_no_pairs_is_empty_per_group():
      out = sum_pairs(X, np.array([], dtype=np.intp), np.array([], dtype=np.intp), n_groups=2)
      assert out.shape == (2, 2) and out.nnz == 0


  def test_sum_pairs_keeps_integer_dtype():
      assert sum_pairs(X, np.array([0]), np.array([0]), n_groups=1).dtype == np.int64


  @given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 50)), st.data())
  def test_sum_pairs_with_one_group_per_feature_equals_sum_by(dense, data):
      codes = np.array(data.draw(st.lists(st.integers(-1, 2), min_size=dense.shape[1], max_size=dense.shape[1])))
      kept = np.flatnonzero(codes >= 0)
      expected = sum_by(sp.csr_matrix(dense), codes, 3)
      np.testing.assert_array_equal(
          sum_pairs(sp.csr_matrix(dense), kept, codes[kept], n_groups=3).toarray(), expected.toarray()
      )
  ```
  Append to `tests/core/test_slots.py`, adding `replace_features` to its
  `from biotapy._core import (...)` list:
  ```python
  def test_replace_features_keeps_samples_and_drops_derived_slots():
      adata = _adata()
      adata.obs["group"] = ["a", "b"]
      add_provenance(adata, "test.step")
      out = replace_features(adata, sp.csr_matrix(np.ones((2, 1))), pd.DataFrame(index=["g1"]))
      assert out.shape == (2, 1) and list(out.var_names) == ["g1"]
      assert out.obs["group"].tolist() == ["a", "b"]
      assert not out.layers.keys() - {None} and not out.obsm and not out.obsp
      assert set(out.uns) == {"biotapy"} and set(out.uns["biotapy"]) == {"x_kind", "provenance"}


  def test_replace_features_does_not_touch_its_input(assert_unchanged):
      adata = _adata()
      before = adata.copy()
      out = replace_features(adata, sp.csr_matrix(np.ones((2, 1))), pd.DataFrame(index=["g1"]))
      add_provenance(out, "test.after")
      assert_unchanged(before, adata)
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_matrix.py tests/core/test_slots.py -q`
  -> `ImportError: cannot import name 'sum_pairs' from 'biotapy._core'`.
- [x] **Step 3: Implement.** Append to `src/biotapy/_core/_matrix.py`:
  ```python
  def sum_pairs(
      X: sp.csr_matrix, features: npt.NDArray[np.intp], groups: npt.NDArray[np.intp], *, n_groups: int
  ) -> sp.csr_matrix:
      """Sum the columns of ``X`` into groups given ``(feature, group)`` membership pairs.

      A feature listed with several groups counts in full toward each of them
      (many-to-many, as ``humann_regroup_table`` does). The pairs are a set: a
      repeated pair counts once. A feature in no pair is left out.
      """
      indicator = sp.csr_matrix(
          (np.ones(features.size, dtype=X.dtype), (features, groups)),
          shape=(X.shape[1], n_groups),
      )
      # The constructor sums repeated (feature, group) entries; membership is a set, so reset them to 1.
      indicator.data[:] = 1
      return sp.csr_matrix(X @ indicator)
  ```
  In `src/biotapy/_core/_slots.py` add `import scipy.sparse as sp` after
  `import pandas as pd`, and append:
  ```python
  def replace_features(adata: AnnData, X: sp.csr_matrix, var: pd.DataFrame) -> AnnData:
      """A new AnnData with ``adata``'s samples and new features ``var`` (an aggregation's groups).

      Keeps ``obs`` and the ``uns['biotapy']`` keys a feature change keeps; every
      slot derived from the old features is left behind, as in ``feature_subset``.
      """
      meta = adata.uns.get("biotapy", {"x_kind": "counts"})
      uns = {"biotapy": {key: meta[key] for key in KEPT_META if key in meta}}
      # anndata types .obs as DataFrame | Dataset2D (its lazy variant); the data model guarantees a DataFrame.
      return AnnData(X=X, obs=cast("pd.DataFrame", adata.obs).copy(), var=var, uns=uns)
  ```
  In `src/biotapy/_core/__init__.py`: import `sum_pairs` from `._matrix` and
  `replace_features` from `._slots` (the `._slots` import becomes a
  parenthesised one-name-per-line list, as `ruff format` writes it), and add
  both to `__all__` in alphabetical order.
- [x] **Step 4: Run, expect pass** - same command -> `36 passed` (15 in
  `test_matrix.py`, 21 in `test_slots.py`).
- [x] **Step 5: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add src/biotapy/_core/_matrix.py src/biotapy/_core/_slots.py src/biotapy/_core/__init__.py \
    tests/core/test_matrix.py tests/core/test_slots.py .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(core): add many-to-many group sums and feature replacement"
  ```

### Task 2.1b: `_core` function tables and the mudata dependency

**Files:** create `src/biotapy/_core/_function.py`, `tests/core/test_function.py`;
modify `src/biotapy/_core/__init__.py`, `pyproject.toml`, `uv.lock`,
`.knowledge/decisions/optional-heavy-dependencies.md`.
**Interfaces (produces):**
- constants `FUNCTION_KEY = "function"`, `BY_TAXON_KEY = "function_by_taxon"`,
  `SPECIAL_FEATURES = ("UNMAPPED", "READS_UNMAPPED", "UNINTEGRATED", "UNGROUPED")`,
  `PROTECTED_FEATURES = ("UNMAPPED", "READS_UNMAPPED", "UNINTEGRATED")`;
- `function_var(row_ids: pd.Index[str]) -> pd.DataFrame`: columns `function,
  name, taxon, genus, species` (pandas `str` dtype) and `special` (bool),
  indexed by `ID` or `ID|taxon`;
- `make_function_mudata(X: object, *, obs: pd.DataFrame, row_ids: pd.Index[str], x_kind: XKind, source: str) -> MuData`:
  the two-modality layout of design note 4, ids checked as in
  `make_treedata` (`ValueError` naming duplicates), one provenance entry
  `source` per modality.

- [x] **Step 1: Add the approved dependency** (approved 2026-10-03, R9.1):
  `uv add 'mudata>=0.4'`. Expected: `pyproject.toml` gains `"mudata>=0.4",`
  between `matplotlib` and `networkx`, and `uv.lock` gains mudata 0.4.1 with
  `scverse-misc[settings]`, `pydantic-settings` and `python-dotenv`. In
  `[tool.mypy] overrides`, before the pooch entries, add:
  ```toml
  # mudata 0.4.1 ships no py.typed marker either; same treatment.
  { module = "mudata", follow_untyped_imports = true, implicit_reexport = true },
  { module = "mudata.*", follow_untyped_imports = true, implicit_reexport = true },
  ```
- [x] **Step 2: Failing tests** - `tests/core/test_function.py`:
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  from mudata import MuData

  from biotapy._core import BY_TAXON_KEY, FUNCTION_KEY, function_var, make_function_mudata

  IDS = pd.Index(
      [
          "UNMAPPED",
          "PWY-1: first pathway",
          "PWY-1: first pathway|g__Bacteroides.s__Bacteroides_ovatus",
          "PWY-1: first pathway|unclassified",
          "K1|g__Blautia.s__Blautia_obeum.t__SGB4810",
          "UNINTEGRATED|bug1",
      ]
  )


  def test_function_var_splits_id_name_and_taxon():
      var = function_var(IDS)
      assert var.index.tolist() == [
          "UNMAPPED",
          "PWY-1",
          "PWY-1|g__Bacteroides.s__Bacteroides_ovatus",
          "PWY-1|unclassified",
          "K1|g__Blautia.s__Blautia_obeum.t__SGB4810",
          "UNINTEGRATED|bug1",
      ]
      assert var["function"].tolist() == ["UNMAPPED", "PWY-1", "PWY-1", "PWY-1", "K1", "UNINTEGRATED"]
      assert var["name"].tolist()[1:4] == ["first pathway"] * 3 and var["name"].isna().tolist()[4:] == [True, True]


  def test_function_var_parses_genus_and_species_from_the_stratum():
      var = function_var(IDS)
      assert var["genus"].tolist()[2] == "Bacteroides" and var["species"].tolist()[2] == "Bacteroides_ovatus"
      # HUMAnN 4 strata end in .t__SGB<id>; it is not a species name.
      assert var["species"].tolist()[4] == "Blautia_obeum"
      # "unclassified" and a bare label are a taxon with no rank.
      assert var[["genus", "species"]].iloc[[0, 1, 3, 5]].isna().all().all()


  def test_function_var_flags_specials():
      assert function_var(IDS)["special"].tolist() == [True, False, False, False, False, True]


  def test_function_var_text_columns_are_the_str_dtype():
      var = function_var(IDS)
      assert all(isinstance(var[column].dtype, pd.StringDtype) for column in ["function", "name", "taxon", "genus"])


  def test_function_var_rejects_two_bars():
      with pytest.raises(ValueError, match="one '|'"):
          function_var(pd.Index(["K1|g__A|extra"]))


  def test_make_function_mudata_splits_community_and_strata():
      X = np.arange(12, dtype=np.float64).reshape(2, 6)
      mdata = make_function_mudata(X, obs=pd.DataFrame(index=["s1", "s2"]), row_ids=IDS, x_kind="cpm", source="test")
      assert isinstance(mdata, MuData)
      assert mdata[FUNCTION_KEY].var_names.tolist() == ["UNMAPPED", "PWY-1"]
      assert list(mdata[FUNCTION_KEY].var.columns) == ["name", "special"]
      assert list(mdata[BY_TAXON_KEY].var.columns) == ["function", "name", "taxon", "genus", "species", "special"]
      np.testing.assert_array_equal(mdata[BY_TAXON_KEY].X.toarray(), X[:, 2:])
      assert mdata[FUNCTION_KEY].uns["biotapy"]["x_kind"] == "cpm"


  def test_make_function_mudata_rejects_repeated_samples():
      with pytest.raises(ValueError, match="duplicate obs ids"):
          make_function_mudata(
              np.ones((2, 1)), obs=pd.DataFrame(index=["s1", "s1"]), row_ids=pd.Index(["K1"]), x_kind="rpk", source="t"
          )


  def test_make_function_mudata_round_trips_through_h5mu(tmp_path):
      import mudata

      mdata = make_function_mudata(
          np.ones((2, 6)), obs=pd.DataFrame(index=["s1", "s2"]), row_ids=IDS, x_kind="rpk", source="test"
      )
      mdata.write_h5mu(tmp_path / "f.h5mu")
      back = mudata.read_h5mu(tmp_path / "f.h5mu")
      assert back[BY_TAXON_KEY].var["taxon"].tolist() == mdata[BY_TAXON_KEY].var["taxon"].tolist()
      assert back[FUNCTION_KEY].var["name"].isna().tolist() == [True, False]
  ```
- [x] **Step 3: Run, expect failure** - `uv run --group test pytest tests/core/test_function.py -q`
  -> `ImportError: cannot import name 'BY_TAXON_KEY' from 'biotapy._core'`.
- [x] **Step 4: Implement** - `src/biotapy/_core/_function.py`:
  ```python
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
  SPECIAL_FEATURES = ("UNMAPPED", "READS_UNMAPPED", "UNINTEGRATED", "UNGROUPED")
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
  ```
  In `src/biotapy/_core/__init__.py` add, first among the imports:
  ```python
  from ._function import (
      BY_TAXON_KEY,
      FUNCTION_KEY,
      PROTECTED_FEATURES,
      SPECIAL_FEATURES,
      function_var,
      make_function_mudata,
  )
  ```
  and add the six names to `__all__` in alphabetical order (`BY_TAXON_KEY`,
  `FUNCTION_KEY` and `PROTECTED_FEATURES` around `PHYLO_KEY`,
  `SPECIAL_FEATURES` after `RANKS`, `function_var` after `feature_subset`,
  `make_function_mudata` before `make_treedata`).
- [x] **Step 5: Run, expect pass** - same command -> `8 passed`. Then
  `uv run --group dev mypy` -> `Success: no issues found`.
- [x] **Step 6: Knowledge.** In `.knowledge/decisions/optional-heavy-dependencies.md`,
  after the Checkpoint D sentence about threadpoolctl, add: "Phase 2 task 2.1b
  added mudata (`>=0.4`, approved 2026-10-03): a HUMAnN or PICRUSt2 table needs
  a community and a per-taxon modality over the same samples, which neither
  AnnData nor TreeData holds. It is pure Python (BSD-3) and brings
  `scverse-misc[settings]` (pydantic-settings, python-dotenv, pydantic)."
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add pyproject.toml uv.lock src/biotapy/_core/_function.py src/biotapy/_core/__init__.py \
    tests/core/test_function.py .knowledge/decisions/optional-heavy-dependencies.md \
    .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(core): build two-modality function tables with mudata"
  ```

### Task 2.3: `io.read_humann`

**Files:** create `src/biotapy/io/_humann.py`, `tests/io/test_humann.py`,
`docs/guide/function.md`; modify `src/biotapy/io/__init__.py`, `docs/api.md`,
`docs/guide/reading_data.md`, `docs/guide/index.md`, `docs/conf.py`,
`.knowledge/contracts/data-model-slots.md`.
**Interfaces:** consumes `make_function_mudata`, `XKind`; produces
`bt.io.read_humann(path: str | Path) -> MuData` (design notes 3 and 4).

- [x] **Step 1: Failing tests** - `tests/io/test_humann.py`:
  ```python
  import gzip
  import tempfile
  from pathlib import Path

  import mudata
  import numpy as np
  import pytest
  from hypothesis import given
  from hypothesis import strategies as st

  import biotapy as bt

  # HUMAnN fixtures: tests/data/humann/NOTICE.txt says which are HUMAnN's (MIT) and which are synthetic.
  DATA = Path(__file__).parents[1] / "data" / "humann"


  def test_reads_community_and_stratified_rows():
      mdata = bt.io.read_humann(DATA / "genefamilies.tsv")
      function, by_taxon = mdata["function"], mdata["function_by_taxon"]
      assert function.var_names.tolist() == [
          "UNMAPPED",
          "UniRef90_A",
          "UniRef90_B",
          "UniRef90_C",
          "UniRef90_D",
          "UniRef90_unknown",
      ]
      assert by_taxon.n_vars == 6 and by_taxon.var["function"].tolist()[:2] == ["UniRef90_A", "UniRef90_A"]
      np.testing.assert_array_equal(function.X.toarray()[:, 1], [8.0, 2.0, 0.0])
      assert function.var.loc["UniRef90_A", "name"] == "alpha protein"


  def test_samples_are_rows_with_the_unit_suffix_removed():
      mdata = bt.io.read_humann(DATA / "genefamilies.tsv")
      assert mdata["function"].obs_names.tolist() == ["S1", "S2", "S3"]
      assert mdata.obs_names.tolist() == ["S1", "S2", "S3"]


  def test_specials_are_features_flagged_in_var():
      function = bt.io.read_humann(DATA / "pathabundance.tsv")["function"]
      assert function.var["special"].tolist() == [True, True, False, False]
      assert bt.io.read_humann(DATA / "pathabundance.tsv")["function_by_taxon"].var["special"].sum() == 2


  @pytest.mark.parametrize(
      ("name", "kind"),
      [("genefamilies.tsv", "rpk"), ("gene_families.tsv", "cpm"), ("pathabundance.tsv", "abundance")],
  )
  def test_x_kind_comes_from_the_header(name, kind):
      mdata = bt.io.read_humann(DATA / name)
      assert mdata["function"].uns["biotapy"]["x_kind"] == kind
      assert mdata["function_by_taxon"].uns["biotapy"]["x_kind"] == kind


  def test_renormalised_column_names_win_over_the_first_cell(tmp_path):
      path = tmp_path / "relab.tsv"
      path.write_text("# Gene Family HUMAnN v4.0.0.alpha.2 Adjusted CPMs\tS1-RELAB\nREADS_UNMAPPED\t0.25\nK1\t0.75\n")
      mdata = bt.io.read_humann(path)
      assert mdata["function"].uns["biotapy"]["x_kind"] == "relative"
      assert mdata.obs_names.tolist() == ["S1"]


  def test_reads_humann_4_style_tables():
      # HUMAnN's own 4.x-style fixture: READS_UNMAPPED and bare strata such as "bug1".
      mdata = bt.io.read_humann(DATA / "gene_families.tsv")
      assert mdata["function"].var.loc["READS_UNMAPPED", "special"]
      assert mdata["function_by_taxon"].var["genus"].isna().all()
      assert mdata.obs_names.tolist() == ["HUMAnN_test"]


  def test_reads_a_merged_humann_3_table():
      mdata = bt.io.read_humann(DATA / "multi_sample_genefamilies.tsv")
      assert mdata["function"].shape == (2, 25) and mdata["function_by_taxon"].shape == (2, 25)
      assert mdata["function_by_taxon"].var["species"].iloc[0] == "Dialister_invisus"


  def test_community_rows_with_no_strata_give_an_empty_modality():
      mdata = bt.io.read_humann(DATA / "demo_pathabundance_with_names.tsv")
      assert mdata["function"].var_names.tolist() == ["UNMAPPED", "UNINTEGRATED"]
      assert mdata["function_by_taxon"].shape == (1, 0)


  def test_reads_gzip(tmp_path):
      path = tmp_path / "genefamilies.tsv.gz"
      path.write_bytes(gzip.compress((DATA / "genefamilies.tsv").read_bytes()))
      assert bt.io.read_humann(path)["function"].n_vars == 6


  def test_single_sample_table(tmp_path):
      path = tmp_path / "one.tsv"
      path.write_text("# Pathway\tS1_Abundance\nUNMAPPED\t1.5\nPWY-1: x\t2.5\nPWY-1: x|unclassified\t2.5\n")
      mdata = bt.io.read_humann(path)
      assert mdata["function"].shape == (1, 2) and mdata["function_by_taxon"].shape == (1, 1)


  def test_all_zero_sample_and_feature_are_kept(tmp_path):
      path = tmp_path / "zeros.tsv"
      path.write_text("# Gene Family\tS1-RPKs\tS2-RPKs\nK1\t0.0\t3.0\nK2\t0.0\t0.0\nK2|unclassified\t0.0\t0.0\n")
      mdata = bt.io.read_humann(path)
      assert mdata["function"].shape == (2, 2) and mdata["function"].X[0].nnz == 0
      assert mdata["function_by_taxon"].var_names.tolist() == ["K2|unclassified"]


  @given(st.lists(st.floats(0, 1e6, allow_nan=False, width=32), min_size=1, max_size=8))
  def test_values_round_trip_through_a_written_table(values):
      rows = "".join(f"K{i}\t{value!r}\nK{i}|unclassified\t{value!r}\n" for i, value in enumerate(values))
      with tempfile.TemporaryDirectory() as tmp:
          path = Path(tmp) / "t.tsv"
          path.write_text("# Gene Family\tS1_Abundance-RPKs\n" + rows)
          mdata = bt.io.read_humann(path)
      # Looser than exact: pandas' default C float parser keeps about 15 significant digits (5.1e-14
      # relative measured); exact parsing costs 2.7x on the 91 MB HMP2 table (2.1 s vs 5.8 s).
      np.testing.assert_allclose(mdata["function"].X.toarray().ravel(), values, rtol=1e-12)
      np.testing.assert_allclose(mdata["function_by_taxon"].X.toarray().ravel(), values, rtol=1e-12)


  def test_round_trips_through_h5mu(tmp_path):
      mdata = bt.io.read_humann(DATA / "genefamilies.tsv")
      mdata.write_h5mu(tmp_path / "g.h5mu")
      back = mudata.read_h5mu(tmp_path / "g.h5mu")
      assert (back["function_by_taxon"].X != mdata["function_by_taxon"].X).nnz == 0


  def test_without_a_hash_line_the_first_line_is_the_header(tmp_path):
      # The HMP2 merged tables: no "#", sample columns named after the joined files.
      path = tmp_path / "pathabundances_3.tsv"
      path.write_text("Feature\\Sample\tA_P_pathabundance_cpm\tB_P_pathabundance_cpm\nUNMAPPED\t334835\t314137\n")
      mdata = bt.io.read_humann(path)
      assert mdata.obs_names.tolist() == ["A_P", "B_P"]
      assert mdata["function"].uns["biotapy"]["x_kind"] == "cpm"


  def test_pathway_coverage_raises(tmp_path):
      path = tmp_path / "cov.tsv"
      path.write_text("# Pathway\tS1_Coverage\nPWY-1\t0.5\n")
      with pytest.raises(ValueError, match="coverage"):
          bt.io.read_humann(path)


  def test_repeated_samples_after_suffix_removal_raise(tmp_path):
      path = tmp_path / "dup.tsv"
      path.write_text("# Gene Family\tS1_Abundance-RPKs\tS1-RPKs\nK1\t1.0\t2.0\n")
      with pytest.raises(ValueError, match="duplicate obs ids"):
          bt.io.read_humann(path)


  def test_header_only_table_has_no_features(tmp_path):
      path = tmp_path / "empty.tsv"
      path.write_text("# Pathway\tS1_Abundance\tS2_Abundance\n")
      mdata = bt.io.read_humann(path)
      assert mdata["function"].shape == (2, 0) and mdata["function_by_taxon"].shape == (2, 0)
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/io/test_humann.py -q`
  -> `AttributeError: module 'biotapy.io' has no attribute 'read_humann'`.
- [x] **Step 3: Implement** - `src/biotapy/io/_humann.py`:
  ```python
  """HUMAnN 3 and 4 tables: gene families, reactions, pathway abundance, and their regrouped or renormalised forms."""

  import csv
  import gzip
  import re
  from pathlib import Path
  from typing import IO

  import numpy as np
  import pandas as pd
  from mudata import MuData

  from biotapy._core import XKind, make_function_mudata

  # Sample-column suffixes: HUMAnN's own ("_Abundance-RPKs", "_Abundance"), renorm --update-snames'
  # ("-CPM", "-RELAB"), and the file names humann_join_tables uses when every file names its
  # sample alike ("<sample>_pathabundance_cpm", as in the HMP2 merged tables).
  _SUFFIX = re.compile(r"(?:_Abundance|_genefamilies|_pathabundance)?(?:[-_](?:RPKs|CPM|RELAB|cpm|relab))?$")
  _COVERAGE = re.compile(r"_(?:Coverage|pathcoverage)")
  # Units a header names, in this order: renorm --update-snames rewrites the sample columns but
  # keeps the first cell, so a "-RELAB" column outranks an "Adjusted CPMs" first cell.
  _UNITS: tuple[tuple[re.Pattern[str], XKind], ...] = (
      (re.compile(r"[-_](?:RELAB|relab)(?:\t|$)"), "relative"),
      (re.compile(r"[-_](?:CPM|cpm)(?:\t|$)|Adjusted CPMs"), "cpm"),
      (re.compile(r"RPKs(?:\t|$)"), "rpk"),
  )


  def read_humann(path: str | Path) -> MuData:
      r"""Read one HUMAnN table into community and per-taxon modalities.

      Parameters
      ----------
      path
          A HUMAnN 3 or 4 output table, per sample or merged by
          ``humann_join_tables``: gene families, reactions or pathway
          abundance, as written or after ``humann_regroup_table`` /
          ``humann_renorm_table``. Gzip (``.gz``) is read directly.

      Returns
      -------
      MuData
          Two modalities over the same samples:

          - ``"function"``: community rows (no ``|``), ``var`` columns
            ``name`` and ``special``;
          - ``"function_by_taxon"``: stratified rows (``ID|taxon``), ``var``
            columns ``function``, ``name``, ``taxon``, ``genus``, ``species``
            and ``special``.

          ``var_names`` are the row ids without their ``": name"`` part.
          ``special`` flags ``UNMAPPED``, ``READS_UNMAPPED``, ``UNINTEGRATED``
          and ``UNGROUPED``, which stay as features. Sample names lose the
          suffix HUMAnN adds (``_Abundance-RPKs``, ``_Abundance``, ``-CPM``,
          ``-RELAB``, or a joined file name's ``_pathabundance_cpm``).

      Raises
      ------
      ValueError
          The file is a pathway coverage table, or a row id holds more than
          one ``|``; sample names repeat once their suffix is removed.

      Notes
      -----
      R equivalent: ``mia::importHUMAnN``
      Guide: :doc:`/guide/reading_data`

      ``uns['biotapy']['x_kind']`` comes from the header: ``-RELAB`` or
      ``_relab`` is ``"relative"``; ``-CPM``, ``_cpm`` or ``Adjusted CPMs`` is
      ``"cpm"``; ``RPKs`` is ``"rpk"``. A header without a unit (pathway
      abundance) is ``"abundance"``, even when the values are whole numbers.
      Renormalise with ``humann_renorm_table --update-snames`` so the header
      names the new unit.

      Read one table per call: a gene family table and a pathway table both
      hold ``UNMAPPED``, so they cannot share a modality.

      The community and stratified rows are kept apart because a pathway's
      community abundance is not the sum of its strata.

      References
      ----------
      Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of
      diverse microbial communities with bioBakery 3. eLife 10:e65088.

      Examples
      --------
      >>> import tempfile
      >>> from pathlib import Path
      >>> import biotapy as bt
      >>> path = Path(tempfile.mkdtemp()) / "genefamilies.tsv"
      >>> _ = path.write_text("# Gene Family\tS1_Abundance-RPKs\nUNMAPPED\t2.0\nK1\t4.0\nK1|g__A.s__A_b\t4.0\n")
      >>> mdata = bt.io.read_humann(path)
      >>> mdata["function"].var_names.tolist(), mdata["function_by_taxon"].var_names.tolist()
      (['UNMAPPED', 'K1'], ['K1|g__A.s__A_b'])
      """
      path = Path(path)
      header, n_comments = _header(path)
      if _COVERAGE.search(header):
          msg = f"path={str(path)!r} is a pathway coverage table; read_humann reads abundance tables"
          raise ValueError(msg)
      table = pd.read_csv(
          path, sep="\t", skiprows=max(n_comments - 1, 0), index_col=0, dtype={0: str}, quoting=csv.QUOTE_NONE
      )
      X = table.to_numpy(dtype=np.float64).T
      obs = pd.DataFrame(index=table.columns.str.replace(_SUFFIX, "", regex=True))
      # HUMAnN never writes raw counts: a table whose header names no unit holds pathway abundances.
      x_kind = _unit(header) or "abundance"
      return make_function_mudata(X, obs=obs, row_ids=table.index, x_kind=x_kind, source="io.read_humann")


  def _open(path: Path) -> IO[str]:
      if path.suffix == ".gz":
          return gzip.open(path, "rt", encoding="utf-8")
      return path.open(encoding="utf-8")


  def _header(path: Path) -> tuple[str, int]:
      """The header line and how many ``#`` lines precede the data.

      HUMAnN's rule: the last ``#`` line is the header; with none, the first line is.
      """
      header, n_comments = "", 0
      with _open(path) as handle:
          for line in handle:
              if not line.startswith("#"):
                  return (header or line), n_comments
              header, n_comments = line, n_comments + 1
      return header, n_comments


  def _unit(header: str) -> XKind | None:
      return next((kind for pattern, kind in _UNITS if pattern.search(header.rstrip("\n"))), None)
  ```
  `src/biotapy/io/__init__.py`: add `from ._humann import read_humann` after
  the `_dada2` import and `"read_humann"` to `__all__` after `"read_dada2"`.
- [x] **Step 4: Run, expect pass** - same command -> `19 passed`; doctest:
  `uv run --group test pytest src/biotapy/io/_humann.py -q` -> `1 passed`.
- [x] **Step 5: Docs.**
  - `docs/conf.py`, `intersphinx_mapping`: add
    `"mudata": ("https://mudata.scverse.org/stable/", None),` after matplotlib.
  - `docs/api.md`: `io.read_humann` after `io.read_dada2`.
  - `docs/guide/index.md`: `function` after `aggregation` in the toctree.
  - Create `docs/guide/function.md` with its introduction and first section:
    ````markdown
    # Function

    `bt.fn` works on functional profiles - gene families, EC numbers, pathways -
    the way `bt.pp.tax_glom` works on taxonomy. Its rules follow HUMAnN's own
    utility scripts, and biotapy's golden tests compare the two on the same files.

    ## Two tables per HUMAnN file

    A HUMAnN table holds community rows (`PWY-5100`) and the same functions split
    by taxon (`PWY-5100|g__Bacteroides.s__Bacteroides_ovatus`). For pathways the
    community value is not the sum of its strata, so neither can be derived from
    the other. `bt.io.read_humann` therefore returns a
    [MuData](https://mudata.scverse.org/) with two modalities over the same
    samples:

    | Modality | Rows | `var` columns |
    |---|---|---|
    | `"function"` | community rows | `name`, `special` |
    | `"function_by_taxon"` | stratified rows | `function`, `name`, `taxon`, `genus`, `species`, `special` |

    ```python
    import biotapy as bt

    mdata = bt.io.read_humann("pathabundance.tsv")
    mdata["function"]  # samples x pathways
    mdata["function_by_taxon"]  # samples x (pathway, taxon) pairs
    ```

    `UNMAPPED`, `READS_UNMAPPED`, `UNINTEGRATED` and `UNGROUPED` stay features,
    flagged in `var["special"]`, so a sample's total keeps what HUMAnN could not
    assign. Read one table per call: a gene family table and a pathway table both
    hold `UNMAPPED`.
    ````
  - Append to `docs/guide/reading_data.md`:
    ````markdown
    ## HUMAnN

    `bt.io.read_humann` reads one HUMAnN 3 or 4 table - gene families, reactions
    or pathway abundance, per sample or merged with `humann_join_tables`, raw or
    after `humann_regroup_table` / `humann_renorm_table` - into a `MuData` with a
    `"function"` (community) and a `"function_by_taxon"` (stratified) modality;
    the [function guide](function.md) explains why there are two.

    ```python
    import biotapy as bt

    mdata = bt.io.read_humann("sample_genefamilies.tsv")
    ```

    - **Header.** The last line starting with `#` is the header, or the first
      line when none does (the HMP2 merged tables).
    - **Ids.** `UniRef90_X: name|g__Genus.s__Species` becomes the feature
      `UniRef90_X|g__Genus.s__Species`, with `name`, `taxon`, `genus` and
      `species` in `var`.
    - **Samples.** HUMAnN's column suffixes (`_Abundance-RPKs`, `_Abundance`,
      `-CPM`, `-RELAB`, a joined file's `_pathabundance_cpm`) are removed.
    - **Units.** `x_kind` comes from the header: `RPKs` is `"rpk"`, `CPM`
      (`Adjusted CPMs` in HUMAnN 4) is `"cpm"`, `RELAB` is `"relative"`. A
      header without a unit - pathway abundance as HUMAnN writes it - is
      `"abundance"`, never `"counts"`. After renormalising with HUMAnN, pass
      `--update-snames` so the header names the new unit.
    - **Not read.** Pathway coverage tables (HUMAnN 3 only) raise a
      `ValueError`: they are not abundances.
    ````
  - Build: `uv run --group doc sphinx-build -W -b html docs docs/_build/html`
    -> `build succeeded.`
- [x] **Step 6: Knowledge** - `.knowledge/contracts/data-model-slots.md`
  (contract addition, user-approved):
  - Convention 2, after "Readers always set it, inferred from the values by
    `_core.infer_x_kind`": insert "- except `io.read_humann`, which reads it
    from the table header (`RPKs` -> `rpk`; `CPM`, `_cpm` or `Adjusted CPMs`
    -> `cpm`; `RELAB`, `_relab` -> `relative`) and labels a header without a
    unit `abundance`, never `counts`".
  - New section after "Conventions":
    ````markdown
    ## Function tables
    `io.read_humann` (and, from Phase 2 slice 2B, `io.read_picrust2`) returns a
    `MuData` built by `_core.make_function_mudata` with two modalities over the
    same samples, each an `AnnData` with its own copy of `obs`:

    | Modality | Features | `var` columns |
    |---|---|---|
    | `"function"` | community rows, e.g. `PWY-1` | `name`, `special` |
    | `"function_by_taxon"` | stratified rows, e.g. `PWY-1|g__Bacteroides.s__Bacteroides_ovatus` | `function`, `name`, `taxon`, `genus`, `species`, `special` |

    `var_names` drop the row's `": name"`. `special` flags `UNMAPPED`,
    `READS_UNMAPPED`, `UNINTEGRATED` and `UNGROUPED`, which stay features. Text
    columns use the pandas `str` dtype, as rank columns do. Both modalities
    always exist; either may have 0 features. The community values are not
    the sum of their strata for pathways, which is why there are two.
    ````
  - Slots table, `var` row: append "; function tables: see Function tables".
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add src/biotapy/io/_humann.py src/biotapy/io/__init__.py tests/io/test_humann.py \
    docs/conf.py docs/api.md docs/guide/index.md docs/guide/function.md docs/guide/reading_data.md \
    .knowledge/contracts/data-model-slots.md .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(io): read HUMAnN 3 and 4 tables into community and per-taxon modalities"
  ```

### Task 2.3b: `datasets.toy_humann`

**Files:** modify `src/biotapy/datasets/_toy.py`, `src/biotapy/datasets/__init__.py`,
`tests/datasets/test_toy.py`, `docs/api.md`, `docs/guide/datasets.md`,
`.knowledge/contracts/function-shape.md`, `rules.md` (R8.2 wording; user-approved).
**Interfaces:** consumes `make_function_mudata`; produces
`bt.datasets.toy_humann() -> MuData` (6 samples `s1`-`s6`, `obs["group"]`
A/B as in `toy()`; `"function"` 6 features: `UNMAPPED`, `UNGROUPED`,
`1.1.1.1`, `2.7.1.1`, `2.7.1.2`, `3.2.1.4`; `"function_by_taxon"` 7;
`x_kind == "rpk"`). Every later `fn` docstring example uses it.

- [x] **Step 1: Failing tests.** Append to `tests/datasets/test_toy.py` (and
  add `import numpy as np` above `import treedata as td`):
  ```python
  def test_toy_humann_has_both_modalities_over_the_toy_samples():
      mdata = bt.datasets.toy_humann()
      assert mdata["function"].shape == (6, 6) and mdata["function_by_taxon"].shape == (6, 7)
      assert mdata["function"].obs_names.tolist() == bt.datasets.toy().obs_names.tolist()
      assert mdata["function"].obs["group"].tolist() == ["A"] * 3 + ["B"] * 3
      assert mdata["function"].uns["biotapy"]["x_kind"] == "rpk"


  def test_toy_humann_community_rows_are_the_sums_of_their_strata():
      mdata = bt.datasets.toy_humann()
      function, by_taxon = mdata["function"], mdata["function_by_taxon"]
      for name in ["UNGROUPED", "1.1.1.1", "2.7.1.1", "2.7.1.2", "3.2.1.4"]:
          strata = by_taxon[:, (by_taxon.var["function"] == name).to_numpy()].X.sum(axis=1)
          np.testing.assert_array_equal(np.asarray(strata).ravel(), function[:, name].X.toarray().ravel())
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/datasets/test_toy.py -q`
  -> `AttributeError: module 'biotapy.datasets' has no attribute 'toy_humann'`.
- [x] **Step 3: Implement.** In `src/biotapy/datasets/_toy.py`: docstring
  becomes `"""Tiny in-memory datasets for docstring examples and tests."""`;
  add `from mudata import MuData` after the pandas import; the `_core` import
  becomes `from biotapy._core import TreeData, make_function_mudata, make_treedata, tree_from_edges`;
  append:
  ```python
  def toy_humann() -> MuData:
      """The six toy samples' gene families, regrouped to EC numbers as HUMAnN writes them.

      Built in memory, so examples and tests never download anything.

      Returns
      -------
      MuData
          ``"function"``: ``UNMAPPED``, ``UNGROUPED`` and four EC numbers;
          ``"function_by_taxon"``: their seven strata over three species and
          ``unclassified``. ``X`` is in RPK (``x_kind == "rpk"``); each
          modality's ``obs['group']`` is ``A`` (s1-s3) or ``B`` (s4-s6), as in
          ``toy()``.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/datasets`

      Examples
      --------
      >>> import biotapy as bt
      >>> mdata = bt.datasets.toy_humann()
      >>> mdata["function"].shape, mdata["function_by_taxon"].shape
      ((6, 6), (6, 7))
      """
      ids, values = zip(*_HUMANN_ROWS, strict=True)
      obs = pd.DataFrame({"group": pd.Categorical(["A"] * 3 + ["B"] * 3)}, index=[f"s{i}" for i in range(1, 7)])
      X = np.array(values, dtype=np.float64).T
      return make_function_mudata(X, obs=obs, row_ids=pd.Index(ids), x_kind="rpk", source="datasets.toy_humann")
  ```
  with, above it:
  ```python
  # HUMAnN gene families regrouped to EC numbers, in RPK; community rows are the sums of their strata.
  _HUMANN_ROWS = (
      ("UNMAPPED", [20, 25, 18, 30, 22, 27]),
      ("UNGROUPED", [4, 3, 5, 4, 2, 3]),
      ("UNGROUPED|unclassified", [4, 3, 5, 4, 2, 3]),
      ("1.1.1.1: alcohol dehydrogenase", [15, 12, 18, 3, 1, 4]),
      ("1.1.1.1: alcohol dehydrogenase|g__Bacteroides.s__Bacteroides_ovatus", [12, 10, 14, 2, 1, 3]),
      ("1.1.1.1: alcohol dehydrogenase|unclassified", [3, 2, 4, 1, 0, 1]),
      ("2.7.1.1: hexokinase", [1, 0, 2, 9, 11, 8]),
      ("2.7.1.1: hexokinase|g__Blautia.s__Blautia_obeum", [1, 0, 2, 9, 11, 8]),
      ("2.7.1.2: glucokinase", [8, 6, 8, 8, 8, 8]),
      ("2.7.1.2: glucokinase|g__Bacteroides.s__Bacteroides_ovatus", [6, 5, 7, 1, 0, 2]),
      ("2.7.1.2: glucokinase|g__Blautia.s__Blautia_obeum", [2, 1, 1, 7, 8, 6]),
      ("3.2.1.4: cellulase", [0, 1, 0, 6, 5, 7]),
      ("3.2.1.4: cellulase|g__Faecalibacterium.s__Faecalibacterium_prausnitzii", [0, 1, 0, 6, 5, 7]),
  )
  ```
  `src/biotapy/datasets/__init__.py`: `from ._toy import toy, toy_humann`;
  `"toy_humann"` at the end of `__all__`.
- [x] **Step 4: Run, expect pass** - same command -> `5 passed`;
  `uv run --group test pytest src/biotapy/datasets -q` -> doctests pass.
- [x] **Step 5: Docs.** `docs/api.md`: `datasets.toy_humann` after
  `datasets.toy`. `docs/guide/datasets.md`: the first sentence becomes
  "`biotapy.datasets` ships five example datasets, each returning the same
  [data model](data_model.md) every biotapy function relies on."; insert
  before `## global_patterns, enterotype and esophagus`:
  ````markdown
  ## `toy_humann`

  `bt.datasets.toy_humann()` is the toy samples' gene families as HUMAnN writes
  them after regrouping to EC numbers, built in memory: a `MuData` whose
  `"function"` modality holds `UNMAPPED`, `UNGROUPED` and four EC numbers in
  RPK, and whose `"function_by_taxon"` modality holds their seven strata.
  `bt.fn` docstrings use it:

  ```python
  import biotapy as bt

  mdata = bt.datasets.toy_humann()
  mdata["function"].shape  # (6, 6)
  ```
  ````
  Build docs -> `build succeeded.`
- [x] **Step 6: Knowledge** (user-approved wording; R12.1).
  - `.knowledge/contracts/function-shape.md`, statement 6, the bullet
    "Examples use `bt.datasets.toy()`" becomes "Examples use
    `bt.datasets.toy()`, or `bt.datasets.toy_humann()` for function tables
    (both built in memory, no download), and run under doctest in CI."
  - `rules.md` R8.2: "a runnable `Examples` section using
    `bt.datasets.toy()`" becomes "a runnable `Examples` section using
    `bt.datasets.toy()` (or `bt.datasets.toy_humann()` for function tables)".
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add src/biotapy/datasets/_toy.py src/biotapy/datasets/__init__.py tests/datasets/test_toy.py \
    docs/api.md docs/guide/datasets.md .knowledge/contracts/function-shape.md rules.md \
    .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(datasets): add toy_humann, the toy samples' EC gene families"
  ```

### Task 2.5a: `datasets.enzyme`

**Files:** create `src/biotapy/datasets/_enzyme.py`, `tests/datasets/test_enzyme.py`,
`tests/data/enzyme/{enzyme.dat, enzclass.txt, NOTICE.txt}`; modify
`src/biotapy/datasets/_remote.py`, `src/biotapy/datasets/__init__.py`,
`docs/api.md`, `docs/guide/datasets.md`.
**Interfaces:** consumes `_remote._fetch(name) -> str` (same subpackage);
produces `bt.datasets.enzyme() -> pd.DataFrame` (columns `child, parent, level,
parent_name`; levels `class`, `subclass`, `subsubclass`; `attrs["source"]`,
`attrs["license"] == "CC BY 4.0"`). 2.6 consumes it as a `hierarchy`.

- [x] **Step 1: Fixtures** - excerpts of the 02-Sep-2026 release (CC BY 4.0;
  the notice states the changes). `tests/data/enzyme/enzyme.dat`:
  ```text
  CC   -----------------------------------------------------------------------
  CC
  CC   ENZYME nomenclature database
  CC
  CC   -----------------------------------------------------------------------
  CC   Release of 02-Sep-2026
  CC   -----------------------------------------------------------------------
  CC
  CC   Alan Bridge and Kristian Axelsen
  CC   SIB Swiss Institute of Bioinformatics
  CC   Centre Medical Universitaire (CMU)
  CC   1, rue Michel Servet
  CC   1211 Geneva 4
  CC   Switzerland
  CC
  CC   Email: enzyme@expasy.org
  CC
  CC   WWW server: https://enzyme.expasy.org/
  CC
  CC   -----------------------------------------------------------------------
  CC   Copyrighted by the SIB Swiss Institute of Bioinformatics and
  CC   distributed under the Creative Commons Attribution (CC BY 4.0) License
  CC   -----------------------------------------------------------------------
  //
  ID   1.1.1.1
  DE   alcohol dehydrogenase.
  //
  ID   1.1.1.74
  DE   Deleted entry.
  //
  ID   2.7.1.1
  DE   hexokinase.
  //
  ID   2.7.1.2
  DE   glucokinase.
  //
  ID   3.2.1.4
  DE   cellulase.
  //
  ```
  `tests/data/enzyme/enzclass.txt`:
  ```text
  ----------------------------------------------------------------------------
          ENZYME nomenclature database
          SIB Swiss Institute of Bioinformatics; Geneva, Switzerland
  ----------------------------------------------------------------------------

  Description: Definition of enzyme classes, subclasses and sub-subclasses
  Name:        enzclass.txt
  Release:     02-Sep-2026

  ----------------------------------------------------------------------------

  1. -. -.-  Oxidoreductases.
  1. 1. -.-   Acting on the CH-OH group of donors.
  1. 1. 1.-    With NAD(+) or NADP(+) as acceptor.
  2. -. -.-  Transferases.
  2. 7. -.-   Transferring phosphorus-containing groups.
  2. 7. 1.-    Phosphotransferases with an alcohol group as acceptor.
  3. -. -.-  Hydrolases.
  3. 2. -.-   Glycosylases.
  3. 2. 1.-    Glycosidases, i.e. enzymes hydrolyzing O- and S-glycosyl compounds.

  ----------------------------------------------------------------------------
  Copyrighted by the SIB Swiss Institute of Bioinformatics and
  distributed under the Creative Commons Attribution (CC BY 4.0) License
  ----------------------------------------------------------------------------
  ```
  `tests/data/enzyme/NOTICE.txt`:
  ```text
  enzyme.dat and enzclass.txt are excerpts of the ENZYME nomenclature database,
  release of 02-Sep-2026 (https://ftp.expasy.org/databases/enzyme/), copyrighted
  by the SIB Swiss Institute of Bioinformatics and distributed under the Creative
  Commons Attribution 4.0 International licence (CC BY 4.0,
  https://creativecommons.org/licenses/by/4.0/).

  Changes: enzyme.dat keeps its header and five entries (1.1.1.1, 1.1.1.74,
  2.7.1.1, 2.7.1.2, 3.2.1.4), each cut to its ID and DE lines; enzclass.txt keeps
  its header, footer and the nine class lines above those entries.
  ```
- [x] **Step 2: Failing tests** - `tests/datasets/test_enzyme.py`:
  ```python
  from pathlib import Path

  import pytest

  import biotapy as bt
  from biotapy.datasets import _enzyme

  # Excerpts of ENZYME (CC BY 4.0); tests/data/enzyme/NOTICE.txt gives the source and changes.
  DATA = Path(__file__).parents[1] / "data" / "enzyme"


  @pytest.fixture
  def offline(monkeypatch):
      # As in test_remote.py: the only offline route to the loader is its private _fetch (R11.4).
      monkeypatch.setattr(_enzyme, "_fetch", lambda name: str(DATA / name))


  def test_lists_every_ancestor_of_an_ec_number(offline):
      edges = bt.datasets.enzyme()
      rows = edges[edges["child"] == "1.1.1.1"]
      assert rows[["parent", "level"]].values.tolist() == [
          ["1.-.-.-", "class"],
          ["1.1.-.-", "subclass"],
          ["1.1.1.-", "subsubclass"],
      ]


  def test_internal_ids_reach_the_levels_above_them(offline):
      edges = bt.datasets.enzyme()
      assert edges[edges["child"] == "2.7.1.-"]["parent"].tolist() == ["2.-.-.-", "2.7.-.-"]
      assert edges[edges["child"] == "2.-.-.-"].empty


  def test_parents_carry_enzclass_names(offline):
      edges = bt.datasets.enzyme()
      names = edges.drop_duplicates("parent").set_index("parent")["parent_name"]
      assert names["1.-.-.-"] == "Oxidoreductases"
      assert names["3.2.1.-"] == "Glycosidases, i.e. enzymes hydrolyzing O- and S-glycosyl compounds"


  def test_deleted_entries_keep_their_place(offline):
      assert "1.1.1.74" in set(bt.datasets.enzyme()["child"])


  def test_attrs_name_the_release_and_licence(offline):
      attrs = bt.datasets.enzyme().attrs
      assert "02-Sep-2026" in attrs["source"] and attrs["license"] == "CC BY 4.0"


  @pytest.mark.network
  def test_enzyme_downloads_and_parses():
      edges = bt.datasets.enzyme()
      assert edges["child"].nunique() > 8000 and set(edges["level"]) == {"class", "subclass", "subsubclass"}
      assert edges.attrs["source"].startswith("ENZYME release ")
  ```
- [x] **Step 3: Run, expect failure** - `uv run --group test pytest tests/datasets/test_enzyme.py -q`
  -> `ImportError: cannot import name '_enzyme' from 'biotapy.datasets'`.
- [x] **Step 4: Implement.** In `src/biotapy/datasets/_remote.py`:
  - module docstring: `"""Datasets downloaded once and cached with pooch: phyloseq's examples and the ENZYME files."""`;
  - end `_REGISTRY` with
    ```python
        # ENZYME keeps no old releases, so no hash can stay valid: the first download is
        # cached for good, and enzyme() records the release it read (datasets/_enzyme.py).
        "enzyme.dat": None,
        "enzclass.txt": None,
    }
    _URLS = {
        "enzyme.dat": "https://ftp.expasy.org/databases/enzyme/enzyme.dat",
        "enzclass.txt": "https://ftp.expasy.org/databases/enzyme/enzclass.txt",
    }
    ```
  - pass `urls=_URLS` to `pooch.create` (pooch 1.9 `create(..., urls=...)`
    overrides `base_url` per file; a `None` hash is never checked,
    `pooch.hashes.hash_matches`).

  `src/biotapy/datasets/_enzyme.py`:
  ```python
  """The ENZYME (EC) hierarchy from the SIB Swiss Institute of Bioinformatics, CC BY 4.0."""

  import re
  from pathlib import Path

  import numpy as np
  import pandas as pd

  from ._remote import _fetch

  LEVELS = ("class", "subclass", "subsubclass")
  # enzclass.txt lines: "1. 1. 1.-    With NAD(+) or NADP(+) as acceptor."
  _CLASS_LINE = re.compile(r"^(\d+\.\s*[\d-]+\.\s*[\d-]+\.-)\s+(.*?)\.?\s*$", re.MULTILINE)
  _ID_LINE = re.compile(r"^ID   (\S+)$", re.MULTILINE)
  _RELEASE = re.compile(r"^CC   Release of (.+)$", re.MULTILINE)
  _URL = "https://enzyme.expasy.org/"


  def enzyme() -> pd.DataFrame:
      """The ENZYME EC hierarchy as an edge table for ``bt.fn.func_glom``.

      Downloaded once (9.6 MB) from ``ftp.expasy.org`` and cached.

      Returns
      -------
      pandas.DataFrame
          One row per EC number and ancestor: columns ``child`` (an EC number
          such as ``1.1.1.1``, or an internal id such as ``1.1.1.-``),
          ``parent`` (its ancestor), ``level`` (the parent's level:
          ``"class"``, ``"subclass"`` or ``"subsubclass"``) and
          ``parent_name`` (from ``enzclass.txt``; NaN when ENZYME names none).
          ``attrs["source"]`` names the ENZYME release read and
          ``attrs["license"]`` is ``"CC BY 4.0"``.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/datasets`

      Every ancestor is listed, not only the direct parent, so a table of EC
      numbers and one already grouped to sub-subclasses both reach any level.
      Deleted and transferred entries keep their number and so their place.

      ENZYME keeps only its current release online, so the first download is
      cached for good and ``attrs["source"]`` records which release it was.
      Delete the cached ``enzyme.dat`` and ``enzclass.txt`` to take a newer one.

      ENZYME is copyrighted by the SIB Swiss Institute of Bioinformatics and
      distributed under the Creative Commons Attribution 4.0 (CC BY 4.0)
      licence; cite it when you publish results that use it.

      References
      ----------
      Bairoch A (2000) The ENZYME database in 2000. Nucleic Acids Res 28:304-305.

      Examples
      --------
      >>> import biotapy as bt
      >>> edges = bt.datasets.enzyme()  # doctest: +SKIP
      >>> edges.query("child == '1.1.1.1'")["parent"].tolist()  # doctest: +SKIP
      ['1.-.-.-', '1.1.-.-', '1.1.1.-']
      """
      entries = Path(_fetch("enzyme.dat")).read_text(encoding="utf-8")
      classes = Path(_fetch("enzclass.txt")).read_text(encoding="utf-8")
      names = {re.sub(r"\s", "", ec): name for ec, name in _CLASS_LINE.findall(classes)}
      leaves = pd.Series(_ID_LINE.findall(entries), dtype=str)
      edges = _ancestors(leaves)
      internal = pd.Series(sorted(set(names) | set(edges["parent"])), dtype=str)
      edges = pd.concat([edges, _ancestors(internal)], ignore_index=True)
      edges["parent_name"] = edges["parent"].map(names).astype(pd.StringDtype(na_value=np.nan))
      release = _RELEASE.search(entries)
      edges.attrs["source"] = f"ENZYME release {release[1] if release else 'unknown'}, SIB, {_URL}"
      edges.attrs["license"] = "CC BY 4.0"
      return edges


  def _ancestors(ids: pd.Series) -> pd.DataFrame:
      """Every (id, ancestor, ancestor's level) row; ``1.1.1.1`` has three, ``1.-.-.-`` none."""
      depth = 4 - ids.str.count(r"\.-")
      parts = ids.str.split(".")
      frames = [
          pd.DataFrame(
              {"child": ids[depth > k], "parent": parts[depth > k].str[:k].str.join(".") + ".-" * (4 - k), "level": level}
          )
          for k, level in enumerate(LEVELS, start=1)
      ]
      return pd.concat(frames, ignore_index=True)
  ```
  `src/biotapy/datasets/__init__.py`: `from ._enzyme import enzyme` first;
  `"enzyme"` after `"enterotype"` in `__all__`.
- [x] **Step 5: Run, expect pass** - same command -> `5 passed, 1 deselected`
  (the `network` test). Once, by hand, with a scratch cache:
  `BIOTAPY_DATA_DIR=$(mktemp -d) uv run --group test pytest -m network tests/datasets/test_enzyme.py -q`
  -> `1 passed` (26,143 rows, 8,875 children, every parent named, for the
  02-Sep-2026 release).
- [x] **Step 6: Docs.** `docs/api.md`: `datasets.enzyme` after
  `datasets.enterotype`. `docs/guide/datasets.md`: first sentence becomes
  "`biotapy.datasets` ships six datasets: five examples, each returning the
  same [data model](data_model.md) every biotapy function relies on, and the
  ENZYME hierarchy for `bt.fn.func_glom`."; insert before `## Licensing`:
  ````markdown
  ## `enzyme`

  `bt.datasets.enzyme()` downloads ENZYME's `enzyme.dat` (9.6 MB) and
  `enzclass.txt` from the SIB Swiss Institute of Bioinformatics and returns the
  EC hierarchy as an edge table: every EC number with its sub-subclass,
  subclass and class, named from `enzclass.txt`. ENZYME keeps only its current
  release online, so no hash can be pinned: the first download is cached for
  good, and `attrs["source"]` records which release it was. Delete the cached
  files to take a newer release.
  ````
  and in `## Licensing` replace "phyloseq's authors under AGPL-3. biotapy
  itself is" with "phyloseq's authors under AGPL-3. `enzyme` downloads ENZYME,
  copyrighted by the SIB Swiss Institute of Bioinformatics and distributed
  under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); cite it when
  you publish results that use it. biotapy itself is". Build docs ->
  `build succeeded.`
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add src/biotapy/datasets/_enzyme.py src/biotapy/datasets/_remote.py src/biotapy/datasets/__init__.py \
    tests/datasets/test_enzyme.py tests/data/enzyme docs/api.md docs/guide/datasets.md \
    .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(datasets): add the ENZYME EC hierarchy (CC BY 4.0) via pooch"
  ```

### Task 2.5: `fn.load_hierarchy`

**Files:** create `src/biotapy/fn/__init__.py`, `src/biotapy/fn/_hierarchy.py`,
`tests/fn/test_hierarchy.py`; modify `src/biotapy/__init__.py`,
`tests/test_docstrings.py`, `docs/api.md`, `docs/guide/function.md`,
`.knowledge/decisions/no-bundled-kegg.md`.
**Interfaces:** produces
`bt.fn.load_hierarchy(path: str | Path, level: str, *, layout: Literal["parent_first", "child_first"] = "parent_first") -> pd.DataFrame`
(columns `child, parent, level, parent_name`; `parent_name` all NaN;
`attrs["source"] == str(path)`). Local files only, gzip allowed, no header
line. This creates the `fn` subpackage (R4.8: in the phase that fills it).

- [x] **Step 1: Failing tests** - `tests/fn/test_hierarchy.py`:
  ```python
  import gzip
  from pathlib import Path

  import pytest

  import biotapy as bt

  DATA = Path(__file__).parents[1] / "data" / "humann"


  def test_parent_first_lines_become_child_parent_rows():
      edges = bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group")
      assert edges[["child", "parent"]].values.tolist() == [
          ["UniRef90_A", "G1"],
          ["UniRef90_B", "G1"],
          ["UniRef90_B", "G2"],
          ["UniRef90_C", "G2"],
          ["UniRef90_E", "G2"],
          ["UniRef90_E", "G3"],
      ]
      assert list(edges.columns) == ["child", "parent", "level", "parent_name"]
      assert set(edges["level"]) == {"group"} and edges["parent_name"].isna().all()


  def test_child_first_reads_two_column_tables(tmp_path):
      path = tmp_path / "pairs.tsv"
      path.write_text("ko:K00001\tpath:map00010\nko:K00001\tpath:map00071\nko:K00002\tpath:map00010\n")
      edges = bt.fn.load_hierarchy(path, "pathway", layout="child_first")
      assert edges[["child", "parent"]].values.tolist() == [
          ["ko:K00001", "path:map00010"],
          ["ko:K00001", "path:map00071"],
          ["ko:K00002", "path:map00010"],
      ]


  def test_repeated_lines_and_blank_lines_give_distinct_pairs(tmp_path):
      path = tmp_path / "map.tsv"
      path.write_text("P1\tK1\n\nP1\tK1\tK2\t\n")
      assert len(bt.fn.load_hierarchy(path, "pathway")) == 2


  def test_reads_gzip(tmp_path):
      path = tmp_path / "map.tsv.gz"
      path.write_bytes(gzip.compress(b"P1\tK1\tK2\n"))
      assert len(bt.fn.load_hierarchy(path, "pathway")) == 2


  def test_attrs_name_the_file():
      assert bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group").attrs["source"].endswith("regroup_map.tsv")


  def test_a_line_with_one_id_is_named(tmp_path):
      path = tmp_path / "map.tsv"
      path.write_text("P1\tK1\nP2\n")
      with pytest.raises(ValueError, match=r"lines \[2\]"):
          bt.fn.load_hierarchy(path, "pathway")


  def test_unknown_layout_raises():
      with pytest.raises(ValueError, match="layout="):
          bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group", layout="wide")
  ```
  In `tests/test_docstrings.py`, `test_every_public_subpackage_is_covered`
  expects `{"datasets", "fn", "io", "pl", "pp", "tl"}`.
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/fn/test_hierarchy.py tests/test_docstrings.py -q`
  -> `AttributeError: module 'biotapy' has no attribute 'fn'`, and the
  subpackage-set assertion fails.
- [x] **Step 3: Implement** - `src/biotapy/fn/_hierarchy.py`:
  ```python
  """Function hierarchies from the user's own mapping files (decisions/no-bundled-kegg)."""

  import gzip
  from pathlib import Path
  from typing import Literal

  import numpy as np
  import pandas as pd


  def load_hierarchy(
      path: str | Path, level: str, *, layout: Literal["parent_first", "child_first"] = "parent_first"
  ) -> pd.DataFrame:
      r"""Read a function mapping file into an edge table for ``bt.fn.func_glom``.

      Parameters
      ----------
      path
          A tab-separated mapping file with no header line, gzip (``.gz``)
          allowed. Each line holds one id and then one or more ids it maps
          to or from (see ``layout``); a line may repeat its first id.
      level
          Name for the level the parents form, for example ``"pathway"``;
          ``func_glom`` selects edges by it.
      layout
          ``"parent_first"``: ``parent<TAB>child<TAB>child...``, the format of
          ``humann_regroup_table --custom`` and PICRUSt2 mapping files.
          ``"child_first"``: ``child<TAB>parent<TAB>parent...``, the format of
          ``humann_regroup_table --reversed`` and of two-column ``child, parent``
          tables.

      Returns
      -------
      pandas.DataFrame
          Columns ``child``, ``parent``, ``level`` and ``parent_name`` (NaN:
          mapping files name no parents), one row per distinct pair.
          ``attrs["source"]`` is the file's path.

      Raises
      ------
      ValueError
          ``layout`` is not one of the two names, or a line holds a single id.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/function`

      biotapy ships no KEGG or MetaCyc mapping: their licences forbid
      redistributing them. Point ``path`` at files you are licensed to use.
      For EC numbers, ``bt.datasets.enzyme()`` gives the open ENZYME
      hierarchy.

      Examples
      --------
      >>> import tempfile
      >>> from pathlib import Path
      >>> import biotapy as bt
      >>> path = Path(tempfile.mkdtemp()) / "map.tsv"
      >>> _ = path.write_text("P1\tK1\tK2\nP2\tK2\n")
      >>> bt.fn.load_hierarchy(path, "pathway")[["child", "parent"]].values.tolist()
      [['K1', 'P1'], ['K2', 'P1'], ['K2', 'P2']]
      """
      if layout not in ("parent_first", "child_first"):
          msg = f"layout={layout!r} must be 'parent_first' or 'child_first'"
          raise ValueError(msg)
      path = Path(path)
      pairs = [(cells[0], other) for cells in _rows(path) for other in cells[1:]]
      first, other = ("parent", "child") if layout == "parent_first" else ("child", "parent")
      edges = pd.DataFrame(pairs, columns=[first, other], dtype=str)[["child", "parent"]].drop_duplicates()
      edges = edges.reset_index(drop=True).assign(level=level, parent_name=np.nan)
      edges["parent_name"] = edges["parent_name"].astype(pd.StringDtype(na_value=np.nan))
      edges.attrs["source"] = str(path)
      return edges


  def _rows(path: Path) -> list[list[str]]:
      """Non-blank lines split on tabs, empty cells dropped; a line with one id raises."""
      opener = gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")
      with opener as handle:
          rows = [(number, line.rstrip("\r\n").split("\t")) for number, line in enumerate(handle, start=1)]
      cells = [(number, [cell for cell in row if cell]) for number, row in rows if any(row)]
      short = [number for number, row in cells if len(row) < 2]
      if short:
          msg = f"path={str(path)!r}: every line needs an id and at least one id it maps to; lines {short[:3]} have one"
          raise ValueError(msg)
      return [row for _, row in cells]
  ```
  `src/biotapy/fn/__init__.py`:
  ```python
  from ._hierarchy import load_hierarchy

  __all__ = ["load_hierarchy"]
  ```
  `src/biotapy/__init__.py`: `from . import datasets, fn, io, pl, pp, tl` and
  `"fn"` after `"datasets"` in `__all__`. import-linter needs no change: `fn`
  is already on the `tl | fn` layer.
- [x] **Step 4: Run, expect pass** - same command -> all pass (7 in
  `test_hierarchy.py`).
- [x] **Step 5: Docs.** `docs/api.md`, new section before `## Tools`:
  ````markdown
  ## Function

  ```{eval-rst}
  .. module:: biotapy.fn
  .. currentmodule:: biotapy

  .. autosummary::
      :toctree: generated

      fn.load_hierarchy
  ```
  ````
  Append to `docs/guide/function.md`:
  ````markdown
  ## Where hierarchies come from

  - **EC numbers:** `bt.datasets.enzyme()` downloads the open ENZYME hierarchy
    (CC BY 4.0) once and caches it; see the [datasets guide](datasets.md).
  - **Your own mapping files:** `bt.fn.load_hierarchy(path, level)` reads
    `humann_regroup_table --custom` files and PICRUSt2 mapping files
    (`parent<TAB>child<TAB>...`, the default `layout="parent_first"`), or
    `child<TAB>parent...` tables with `layout="child_first"`.

  biotapy ships and downloads no KEGG or MetaCyc mapping: their licences do not
  allow it to. If you hold a KEGG or MetaCyc licence, export the mapping you
  need and read it with `load_hierarchy`.
  ````
  Build docs -> `build succeeded.`
- [x] **Step 6: Knowledge** - `.knowledge/decisions/no-bundled-kegg.md`
  (R12.1: it says loaders take "a path or URL" and names MetaCyc as open):
  - `description`: "Function hierarchies come from the user's local files or
    from ENZYME (CC BY 4.0, `bt.datasets.enzyme`); KEGG and MetaCyc are never
    shipped or fetched."
  - Decision bullets become: "`fn.load_hierarchy` reads local files only, no
    URL. The one built-in download is ENZYME's EC hierarchy, in `datasets`
    with the other pooch fetches. No KEGG or MetaCyc mapping is committed,
    packaged, downloaded or used as a fixture: MetaCyc has been
    subscription-only since 2024, and KEGG's REST API is for academic users
    only. Test fixtures are synthetic maps and a CC BY 4.0 ENZYME excerpt."
  - Consequences: "KEGG and MetaCyc users point `load_hierarchy` at files they
    are licensed to use."
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add src/biotapy/fn src/biotapy/__init__.py tests/fn/test_hierarchy.py tests/test_docstrings.py \
    docs/api.md docs/guide/function.md .knowledge/decisions/no-bundled-kegg.md \
    .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(fn): read function hierarchies from local mapping files"
  ```

### Task 2.6: `fn.func_glom`

**Files:** create `src/biotapy/fn/_glom.py`, `tests/fn/test_glom.py`,
`tests/fn/test_glom_golden.py`; modify `src/biotapy/fn/__init__.py`,
`tests/fn/test_hierarchy.py`, `tests/datasets/test_enzyme.py`, `docs/api.md`,
`docs/guide/function.md`, `.knowledge/contracts/data-model-slots.md`.
**Interfaces:** consumes `sum_pairs`, `replace_features`, `add_provenance`,
`as_csr`, `PROTECTED_FEATURES`, `SPECIAL_FEATURES`, `RANKS`, the edge table
(2.5, 2.5a); produces
`bt.fn.func_glom(adata: AnnData, level: str, *, hierarchy: pd.DataFrame, agg: Literal["sum", "mean"] = "sum") -> AnnData`
(design note 2). 2.7 and 2.10 consume its output.

- [x] **Step 1: Failing tests** - `tests/fn/test_glom.py`:
  ```python
  import anndata as ad
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from anndata import AnnData
  from hypothesis import given
  from hypothesis import strategies as st

  import biotapy as bt

  EC = pd.DataFrame(
      {
          "child": ["1.1.1.1", "2.7.1.1", "2.7.1.2", "2.7.1.1", "2.7.1.2"],
          "parent": ["1.-.-.-", "2.-.-.-", "2.-.-.-", "kinase", "kinase"],
          "level": ["class", "class", "class", "role", "role"],
          "parent_name": ["Oxidoreductases", "Transferases", "Transferases", np.nan, np.nan],
      }
  )


  def _function():
      return bt.datasets.toy_humann()["function"]


  def _by_taxon():
      return bt.datasets.toy_humann()["function_by_taxon"]


  def _column(adata, name):
      return adata[:, name].X.toarray().ravel()


  def test_sums_children_into_parents():
      out = bt.fn.func_glom(_function(), "class", hierarchy=EC)
      assert out.var_names.tolist() == ["1.-.-.-", "2.-.-.-", "UNGROUPED", "UNMAPPED"]
      np.testing.assert_array_equal(
          _column(out, "2.-.-.-"), _column(_function(), "2.7.1.1") + _column(_function(), "2.7.1.2")
      )


  def test_many_to_many_counts_a_child_in_every_parent():
      hierarchy = pd.DataFrame({"child": ["2.7.1.1", "2.7.1.1"], "parent": ["P1", "P2"], "level": "pathway"})
      out = bt.fn.func_glom(_function(), "pathway", hierarchy=hierarchy)
      np.testing.assert_array_equal(_column(out, "P1"), _column(_function(), "2.7.1.1"))
      np.testing.assert_array_equal(_column(out, "P2"), _column(_function(), "2.7.1.1"))


  def test_unmapped_features_go_to_ungrouped_and_protected_specials_pass_through():
      out = bt.fn.func_glom(_function(), "class", hierarchy=EC)
      # 3.2.1.4 has no parent here; the input's own UNGROUPED row joins it.
      np.testing.assert_array_equal(
          _column(out, "UNGROUPED"), _column(_function(), "UNGROUPED") + _column(_function(), "3.2.1.4")
      )
      np.testing.assert_array_equal(_column(out, "UNMAPPED"), _column(_function(), "UNMAPPED"))
      assert out.var["special"].tolist() == [False, False, True, True]


  def test_parent_names_become_var_name():
      out = bt.fn.func_glom(_function(), "class", hierarchy=EC)
      assert out.var["name"].tolist()[:2] == ["Oxidoreductases", "Transferases"] and out.var["name"].isna().tolist()[
          2:
      ] == [True, True]


  def test_stratified_input_is_grouped_per_taxon():
      out = bt.fn.func_glom(_by_taxon(), "role", hierarchy=EC)
      assert out.var_names.tolist() == [
          "UNGROUPED|g__Bacteroides.s__Bacteroides_ovatus",
          "UNGROUPED|g__Faecalibacterium.s__Faecalibacterium_prausnitzii",
          "UNGROUPED|unclassified",
          "kinase|g__Bacteroides.s__Bacteroides_ovatus",
          "kinase|g__Blautia.s__Blautia_obeum",
      ]
      assert out.var.loc["kinase|g__Blautia.s__Blautia_obeum", ["function", "taxon", "genus"]].tolist() == [
          "kinase",
          "g__Blautia.s__Blautia_obeum",
          "Blautia",
      ]
      by_taxon = _by_taxon()
      np.testing.assert_array_equal(
          _column(out, "kinase|g__Blautia.s__Blautia_obeum"),
          _column(by_taxon, "2.7.1.1|g__Blautia.s__Blautia_obeum")
          + _column(by_taxon, "2.7.1.2|g__Blautia.s__Blautia_obeum"),
      )


  def test_mean_divides_by_members_present():
      hierarchy = pd.DataFrame({"child": ["2.7.1.1", "2.7.1.2", "9.9.9.9"], "parent": "P", "level": "pathway"})
      out = bt.fn.func_glom(_function(), "pathway", hierarchy=hierarchy, agg="mean")
      np.testing.assert_allclose(
          _column(out, "P"), (_column(_function(), "2.7.1.1") + _column(_function(), "2.7.1.2")) / 2
      )


  def test_keeps_obs_and_x_kind_and_drops_derived_slots():
      function = bt.pp.relative(_function())
      out = bt.fn.func_glom(function, "class", hierarchy=EC)
      assert out.obs["group"].tolist() == function.obs["group"].tolist()
      assert out.uns["biotapy"]["x_kind"] == "rpk" and "relative" not in out.layers
      assert '"step": "fn.func_glom"' in out.uns["biotapy"]["provenance"][-1]


  def test_output_round_trips_through_h5ad(tmp_path):
      out = bt.fn.func_glom(_by_taxon(), "role", hierarchy=EC)
      out.write_h5ad(tmp_path / "glom.h5ad")
      back = ad.read_h5ad(tmp_path / "glom.h5ad")
      assert back.var["taxon"].tolist() == out.var["taxon"].tolist() and back.var["name"].isna().all()


  def test_input_unchanged(assert_unchanged):
      function = _function()
      before = function.copy()
      bt.fn.func_glom(function, "class", hierarchy=EC)
      assert_unchanged(before, function)


  def test_all_zero_sample_stays_zero():
      function = _function()
      dense = function.X.toarray()
      dense[0] = 0
      function.X = sp.csr_matrix(dense)
      assert bt.fn.func_glom(function, "class", hierarchy=EC).X[0].nnz == 0


  def test_all_zero_feature_keeps_its_group():
      function = _function()
      dense = function.X.toarray()
      dense[:, function.var_names.get_loc("1.1.1.1")] = 0
      function.X = sp.csr_matrix(dense)
      out = bt.fn.func_glom(function, "class", hierarchy=EC)
      assert "1.-.-.-" in out.var_names and _column(out, "1.-.-.-").sum() == 0


  def test_single_sample():
      out = bt.fn.func_glom(_function()[:1].copy(), "class", hierarchy=EC)
      assert out.shape == (1, 4)


  def test_empty_modality_gives_no_groups():
      empty = _by_taxon()[:, []].copy()
      assert bt.fn.func_glom(empty, "class", hierarchy=EC).shape == (6, 0)


  def test_no_feature_in_the_hierarchy_raises_with_examples():
      hierarchy = pd.DataFrame({"child": ["EC:2.7.1.1"], "parent": ["P"], "level": "pathway"})
      with pytest.raises(
          ValueError, match=r"features: \['UNGROUPED', '1.1.1.1', '2.7.1.1'\], children: \['EC:2.7.1.1'\]"
      ):
          bt.fn.func_glom(_function(), "pathway", hierarchy=hierarchy)


  def test_unknown_level_names_the_levels():
      with pytest.raises(KeyError, match=r"its levels: \['class', 'role'\]"):
          bt.fn.func_glom(_function(), "pathway", hierarchy=EC)


  def test_missing_hierarchy_column_raises():
      with pytest.raises(KeyError, match="missing"):
          bt.fn.func_glom(_function(), "class", hierarchy=EC.drop(columns="level"))


  def test_unknown_agg_raises():
      with pytest.raises(ValueError, match="agg="):
          bt.fn.func_glom(_function(), "class", hierarchy=EC, agg="median")


  @st.composite
  def _cases(draw):
      n_obs, n_var = draw(st.integers(1, 4)), draw(st.integers(1, 6))
      x = np.array(draw(st.lists(st.integers(0, 20), min_size=n_obs * n_var, max_size=n_obs * n_var))).reshape(
          n_obs, n_var
      )
      pairs = draw(st.sets(st.tuples(st.integers(0, n_var - 1), st.integers(0, 2)), min_size=1))
      return x, sorted(pairs)


  @given(_cases())
  def test_total_equals_each_feature_times_its_parent_count(case):
      # The roadmap's invariant: mapped abundance x membership count is preserved; unmapped goes to UNGROUPED.
      x, pairs = case
      adata = AnnData(
          X=sp.csr_matrix(x),
          obs=pd.DataFrame(index=[f"s{i}" for i in range(x.shape[0])]),
          var=pd.DataFrame(index=[f"f{j}" for j in range(x.shape[1])]),
      )
      hierarchy = pd.DataFrame(
          {"child": [f"f{j}" for j, _ in pairs], "parent": [f"P{g}" for _, g in pairs], "level": "l"}
      )
      out = bt.fn.func_glom(adata, "l", hierarchy=hierarchy)
      parents = np.bincount([j for j, _ in pairs], minlength=x.shape[1])
      weights = np.where(parents > 0, parents, 1)
      np.testing.assert_array_equal(np.asarray(out.X.sum(axis=1)).ravel(), x @ weights)
  ```
  `tests/fn/test_glom_golden.py`:
  ```python
  from pathlib import Path

  import mudata
  import numpy as np
  import pandas as pd
  import pytest

  import biotapy as bt

  TESTS = Path(__file__).parents[1]
  pytestmark = pytest.mark.golden


  def _as_table(mdata: mudata.MuData) -> pd.DataFrame:
      """Both modalities side by side, samples x HUMAnN row ids (``ID`` or ``ID|taxon``)."""
      mods = [mdata[key] for key in ("function", "function_by_taxon")]
      return pd.concat([pd.DataFrame(m.X.toarray(), index=m.obs_names, columns=m.var_names) for m in mods], axis=1)


  @pytest.mark.parametrize("agg", ["sum", "mean"])
  def test_func_glom_matches_humann_regroup_table(agg):
      golden = pd.read_csv(TESTS / "golden" / "humann" / f"regroup_{agg}.csv.gz", index_col="sample_id")
      mdata = bt.io.read_humann(TESTS / "data" / "humann" / "genefamilies.tsv")
      hierarchy = bt.fn.load_hierarchy(TESTS / "data" / "humann" / "regroup_map.tsv", "group")
      out = mudata.MuData(
          {key: bt.fn.func_glom(mod, "group", hierarchy=hierarchy, agg=agg) for key, mod in mdata.mod.items()}
      )
      table = _as_table(out)
      assert sorted(table.columns) == sorted(golden.columns)
      np.testing.assert_allclose(table[golden.columns].to_numpy(), golden.to_numpy(), rtol=1e-7)


  def test_golden_files_come_from_humann_3_9():
      assert (TESTS / "golden" / "humann" / "VERSIONS.txt").read_text().startswith("humann 3.9\n")
  ```
  Append to `tests/fn/test_hierarchy.py`:
  ```python
  def test_round_trips_through_func_glom():
      edges = bt.fn.load_hierarchy(DATA / "regroup_map.tsv", "group")
      out = bt.fn.func_glom(bt.io.read_humann(DATA / "genefamilies.tsv")["function"], "group", hierarchy=edges)
      assert out.var_names.tolist() == ["G1", "G2", "UNGROUPED", "UNMAPPED"]
  ```
  Append to `tests/datasets/test_enzyme.py`:
  ```python
  def test_feeds_func_glom(offline):
      out = bt.fn.func_glom(bt.datasets.toy_humann()["function"], "class", hierarchy=bt.datasets.enzyme())
      assert out.var_names.tolist() == ["1.-.-.-", "2.-.-.-", "3.-.-.-", "UNGROUPED", "UNMAPPED"]
      assert out.var.loc["2.-.-.-", "name"] == "Transferases"
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/fn tests/datasets/test_enzyme.py -q`
  -> `AttributeError: module 'biotapy.fn' has no attribute 'func_glom'`.
- [x] **Step 3: Implement** - `src/biotapy/fn/_glom.py`:
  ```python
  """Aggregation along a function hierarchy, with humann_regroup_table's semantics."""

  from typing import Literal, cast

  import numpy as np
  import pandas as pd
  import scipy.sparse as sp
  from anndata import AnnData

  from biotapy._core import (
      PROTECTED_FEATURES,
      RANKS,
      SPECIAL_FEATURES,
      add_provenance,
      as_csr,
      replace_features,
      sum_pairs,
  )

  UNGROUPED = "UNGROUPED"
  _EDGE_COLUMNS = ("child", "parent", "level")
  # var columns that describe a stratum's taxon; they are the same for every member of a group.
  _TAXON_COLUMNS = ("taxon", *RANKS)


  def func_glom(adata: AnnData, level: str, *, hierarchy: pd.DataFrame, agg: Literal["sum", "mean"] = "sum") -> AnnData:
      """Aggregate function features to one level of a hierarchy.

      Each feature's abundance counts in full toward every parent it has at
      ``level`` (one KO in two pathways counts in both), as
      ``humann_regroup_table`` does.

      Parameters
      ----------
      adata
          Samples x functions: one modality of ``bt.io.read_humann``'s result,
          or any AnnData whose ``var_names`` are function ids. A
          ``"function_by_taxon"`` modality (``var`` columns ``function`` and
          ``taxon``) is grouped per taxon.
      level
          The ``level`` value of the hierarchy rows to use, for example
          ``"pathway"`` or ``"class"``.
      hierarchy
          Edge table with columns ``child``, ``parent`` and ``level`` (and
          optionally ``parent_name``), as ``bt.fn.load_hierarchy`` and
          ``bt.datasets.enzyme`` return.
      agg
          ``"sum"``, or ``"mean"`` over the members present in ``adata``.

      Returns
      -------
      AnnData
          Samples x groups, sorted by name. Features with no parent at
          ``level`` are summed into ``UNGROUPED`` (per taxon when stratified);
          ``UNMAPPED``, ``READS_UNMAPPED`` and ``UNINTEGRATED`` pass through
          unchanged. ``var`` holds ``name`` (from ``parent_name``) and
          ``special``; a stratified input also keeps ``function``, ``taxon``
          and its rank columns. ``obs`` and ``uns['biotapy']['x_kind']`` are
          kept; ``layers``, ``obsm``, ``obsp``, ``varm`` and ``varp`` are
          dropped because they described the old features.

      Raises
      ------
      KeyError
          ``hierarchy`` lacks a required column, or has no row at ``level``.
      ValueError
          ``agg`` is not ``"sum"`` or ``"mean"``; no feature of ``adata``
          is a child at ``level``.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/function`

      Matches ``humann_regroup_table`` (HUMAnN 3.9) with its defaults
      ``--ungrouped Y --protected Y``; the golden tests compare the two. As
      there, a mean divides by the number of members present in the table,
      not by the group's size in ``hierarchy``.

      References
      ----------
      Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of
      diverse microbial communities with bioBakery 3. eLife 10:e65088.

      Examples
      --------
      >>> import pandas as pd
      >>> import biotapy as bt
      >>> edges = pd.DataFrame(
      ...     {"child": ["2.7.1.1", "2.7.1.2"], "parent": ["2.7.1.-", "2.7.1.-"], "level": "subsubclass"}
      ... )
      >>> out = bt.fn.func_glom(bt.datasets.toy_humann()["function"], "subsubclass", hierarchy=edges)
      >>> out.var_names.tolist()
      ['2.7.1.-', 'UNGROUPED', 'UNMAPPED']
      """
      if agg not in ("sum", "mean"):
          msg = f"agg={agg!r} must be 'sum' or 'mean'"
          raise ValueError(msg)
      edges = _edges_at(hierarchy, level)
      var = cast("pd.DataFrame", adata.var)
      function = var["function"] if "function" in var.columns else var.index.to_series()
      pairs = _pairs(function.astype(str).to_numpy(), edges, level=level)
      key = pairs["group"]
      if "taxon" in var.columns:
          key = key + "|" + var["taxon"].astype(str).to_numpy()[pairs["feature"]]
      codes, labels = pd.factorize(key, sort=True)
      X = sum_pairs(as_csr(adata.X), pairs["feature"].to_numpy(dtype=np.intp), codes, n_groups=labels.size)
      if agg == "mean":
          X = sp.csr_matrix(X @ sp.diags(1.0 / np.bincount(codes, minlength=labels.size)))
      out = replace_features(adata, X, _group_var(var, pairs.assign(code=codes), labels, edges=edges))
      add_provenance(out, "fn.func_glom", level=level, agg=agg, hierarchy=hierarchy.attrs.get("source"))
      return out


  def _edges_at(hierarchy: pd.DataFrame, level: str) -> pd.DataFrame:
      missing = [column for column in _EDGE_COLUMNS if column not in hierarchy.columns]
      if missing:
          msg = f"hierarchy needs columns {list(_EDGE_COLUMNS)}; missing {missing}"
          raise KeyError(msg)
      edges = hierarchy[hierarchy["level"] == level]
      if edges.empty:
          msg = f"hierarchy has no row at level={level!r}; its levels: {sorted(hierarchy['level'].unique().tolist())}"
          raise KeyError(msg)
      return edges.drop_duplicates(["child", "parent"])


  def _pairs(function: np.ndarray, edges: pd.DataFrame, *, level: str) -> pd.DataFrame:
      """(feature, group) rows: each parent of a feature, else itself if protected, else UNGROUPED."""
      features = pd.DataFrame({"feature": np.arange(function.size), "child": function})
      protected = features["child"].isin(PROTECTED_FEATURES)
      mapped = features[~protected].merge(edges[["child", "parent"]], on="child")
      if not protected.all() and mapped.empty:
          msg = (
              f"no feature of adata is a child at level={level!r}; "
              f"features: {features['child'][~protected][:3].tolist()}, children: {edges['child'][:3].tolist()}"
          )
          raise ValueError(msg)
      alone = features[protected | ~features["feature"].isin(mapped["feature"])]
      alone = alone.assign(parent=alone["child"].where(alone["child"].isin(PROTECTED_FEATURES), UNGROUPED))
      rows = pd.concat([mapped, alone], ignore_index=True)
      return rows[["feature", "parent"]].rename(columns={"parent": "group"})


  def _group_var(var: pd.DataFrame, pairs: pd.DataFrame, labels: pd.Index, *, edges: pd.DataFrame) -> pd.DataFrame:
      """One var row per group (``pairs["code"]``): function, name, special, and a stratified input's taxon columns."""
      group = pairs.groupby("code")["group"].first()
      first = pairs.groupby("code")["feature"].min().to_numpy()
      names = edges.drop_duplicates("parent").set_index("parent").get("parent_name")
      # The pandas str dtype, never object: anndata's writer rejects an all-NaN object column (contracts/data-model-slots).
      text = pd.StringDtype(na_value=np.nan)
      out = pd.DataFrame(index=labels)
      if "taxon" in var.columns:
          out["function"] = pd.array(group.to_numpy(), dtype=text)
      out["name"] = pd.array(group.map(names).to_numpy() if names is not None else [np.nan] * len(group), dtype=text)
      for column in [column for column in _TAXON_COLUMNS if column in var.columns]:
          out[column] = var[column].array[first]
      out["special"] = group.isin(SPECIAL_FEATURES).to_numpy()
      return out
  ```
  `src/biotapy/fn/__init__.py`: `from ._glom import func_glom` first;
  `__all__ = ["func_glom", "load_hierarchy"]`.
- [x] **Step 4: Run, expect pass** - same command -> all pass (18 in
  `test_glom.py`, 3 golden, 8 in `test_hierarchy.py`, 6 in `test_enzyme.py`;
  the `network` test deselected).
  Doctest: `uv run --group test pytest src/biotapy/fn -q`.
- [x] **Step 5: Docs.** `docs/api.md`: `fn.func_glom` before
  `fn.load_hierarchy`. In `docs/guide/function.md`, insert before
  `## Where hierarchies come from`:
  ````markdown
  ## Aggregating along a hierarchy

  `bt.fn.func_glom` sums functions into their parents at one level of a
  hierarchy, with `humann_regroup_table`'s rules:

  - **Many-to-many.** A function with two parents counts in full toward both,
    so a level's total can exceed the input's.
  - **`UNGROUPED`.** Functions with no parent at that level are summed into
    `UNGROUPED`, per taxon in a stratified table.
  - **Specials pass through.** `UNMAPPED`, `READS_UNMAPPED` and `UNINTEGRATED`
    keep their own rows.
  - **`agg="mean"`** divides by the members present in the table, not by the
    group's size in the hierarchy.

  Apply it to each modality to keep the pair together:

  ```python
  import mudata

  import biotapy as bt

  mdata = bt.datasets.toy_humann()
  edges = bt.datasets.enzyme()
  by_class = mudata.MuData({key: bt.fn.func_glom(mod, "class", hierarchy=edges) for key, mod in mdata.mod.items()})
  ```

  A hierarchy is an edge table with columns `child`, `parent`, `level` and,
  optionally, `parent_name`, which becomes `var["name"]`. Every ancestor is a
  row, not only the direct parent, so a table already grouped one level up
  still reaches the levels above.
  ````
  Build docs -> `build succeeded.`
- [x] **Step 6: Knowledge** - `.knowledge/contracts/data-model-slots.md`
  (user-approved):
  - Propagation table, "Feature-changing" row: add `fn.func_glom` to the
    operations.
  - The sentence "Feature-changing operations go through `_core.feature_subset`,
    the single place that implements the "Drops" column." becomes "...go
    through `_core.feature_subset`, or `_core.replace_features` when the new
    features are groups rather than a subset; both keep only `KEPT_META`, the
    single definition of the "Drops" column."
  - New section after "Aggregation semantics (`tax_glom`)":
    ````markdown
    ## Aggregation semantics (`func_glom`)
    Matches `humann_regroup_table` (HUMAnN 3.9, `--ungrouped Y --protected Y`):
    a feature counts in full toward every parent it has at the level; features
    with none are summed into `UNGROUPED` (per taxon when `var` has `taxon`);
    `UNMAPPED`, `READS_UNMAPPED` and `UNINTEGRATED` pass through; `agg="mean"`
    divides by the members present. Groups are sorted by name. `var` holds
    `name` (from the hierarchy's `parent_name`) and `special`, plus `function`,
    `taxon` and rank columns for a stratified input.
    ````
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add src/biotapy/fn/_glom.py src/biotapy/fn/__init__.py tests/fn/test_glom.py tests/fn/test_glom_golden.py \
    tests/fn/test_hierarchy.py tests/datasets/test_enzyme.py docs/api.md docs/guide/function.md \
    .knowledge/contracts/data-model-slots.md .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(fn): add func_glom with humann_regroup_table semantics"
  ```

### Task 2.12: `fn.renorm`

**Files:** create `src/biotapy/fn/_renorm.py`, `tests/fn/test_renorm.py`,
`tests/fn/test_renorm_golden.py`; modify `src/biotapy/fn/__init__.py`,
`pyproject.toml` (mypy `untyped_calls_exclude`), `docs/api.md`,
`docs/guide/function.md`, `.knowledge/contracts/data-model-slots.md`.
**Interfaces:** consumes `feature_subset`, `add_provenance`, `as_csr`,
`FUNCTION_KEY`, `BY_TAXON_KEY`, `SPECIAL_FEATURES`; produces
`bt.fn.renorm(mdata: MuData, units: Literal["relab", "cpm"], *, special: bool = True) -> MuData`
(design note 1).

- [x] **Step 1: Failing tests** - `tests/fn/test_renorm.py`:
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from hypothesis import given
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays

  import biotapy as bt
  from biotapy._core import make_function_mudata


  def _row_totals(adata):
      return np.asarray(adata.X.sum(axis=1)).ravel()


  def test_community_rows_sum_to_one_per_sample():
      out = bt.fn.renorm(bt.datasets.toy_humann(), "relab")
      np.testing.assert_allclose(_row_totals(out["function"]), 1.0)
      assert out["function"].uns["biotapy"]["x_kind"] == "relative"


  def test_strata_are_divided_by_the_community_total():
      mdata = bt.datasets.toy_humann()
      out = bt.fn.renorm(mdata, "cpm")
      totals = _row_totals(mdata["function"])
      np.testing.assert_allclose(
          out["function_by_taxon"].X.toarray(), mdata["function_by_taxon"].X.toarray() / totals[:, None] * 1e6
      )
      assert out["function_by_taxon"].uns["biotapy"]["x_kind"] == "cpm"


  def test_special_false_drops_specials_before_totalling():
      out = bt.fn.renorm(bt.datasets.toy_humann(), "relab", special=False)
      assert not out["function"].var["special"].any() and not out["function_by_taxon"].var["special"].any()
      np.testing.assert_allclose(_row_totals(out["function"]), 1.0)


  def test_keeps_global_obs_and_other_modalities():
      mdata = bt.datasets.toy_humann()
      mdata.obs["subject"] = ["a", "b", "c", "d", "e", "f"]
      out = bt.fn.renorm(mdata, "relab")
      assert out.obs["subject"].tolist() == ["a", "b", "c", "d", "e", "f"]
      assert '"step": "fn.renorm"' in out["function"].uns["biotapy"]["provenance"][-1]


  def test_input_unchanged(assert_unchanged):
      mdata = bt.datasets.toy_humann()
      before = mdata.copy()
      bt.fn.renorm(mdata, "relab", special=False)
      for key in ("function", "function_by_taxon"):
          assert_unchanged(before[key], mdata[key])


  def test_all_zero_sample_stays_zero():
      mdata = bt.datasets.toy_humann()
      for key in ("function", "function_by_taxon"):
          dense = mdata[key].X.toarray()
          dense[0] = 0
          mdata.mod[key].X = sp.csr_matrix(dense)
      out = bt.fn.renorm(mdata, "relab")
      assert out["function"].X[0].nnz == 0 and np.isfinite(out["function_by_taxon"].X.toarray()).all()


  def test_single_sample():
      out = bt.fn.renorm(bt.datasets.toy_humann()[:1].copy(), "relab")
      np.testing.assert_allclose(_row_totals(out["function"]), [1.0])


  def test_unknown_units_raise():
      with pytest.raises(ValueError, match="units="):
          bt.fn.renorm(bt.datasets.toy_humann(), "tpm")


  def test_missing_modality_names_it():
      mdata = bt.datasets.toy_humann()
      del mdata.mod["function_by_taxon"]
      with pytest.raises(KeyError, match="function_by_taxon"):
          bt.fn.renorm(mdata, "relab")


  def test_stratified_only_table_raises():
      mdata = bt.datasets.toy_humann()
      mdata.mod["function"] = mdata["function"][:, []].copy()
      with pytest.raises(ValueError, match="community"):
          bt.fn.renorm(mdata, "relab")


  @given(arrays(np.float64, st.tuples(st.integers(1, 4), st.integers(1, 5)), elements=st.floats(0, 1e6)))
  def test_community_totals_are_one_or_zero(dense):
      obs = pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])])
      ids = pd.Index([f"K{j}" for j in range(dense.shape[1])] + [f"K{j}|unclassified" for j in range(dense.shape[1])])
      mdata = make_function_mudata(np.hstack([dense, dense]), obs=obs, row_ids=ids, x_kind="rpk", source="test")
      totals = _row_totals(bt.fn.renorm(mdata, "relab")["function"])
      np.testing.assert_allclose(totals, np.where(dense.sum(axis=1) > 0, 1.0, 0.0), rtol=1e-12)
  ```
  `tests/fn/test_renorm_golden.py`:
  ```python
  from pathlib import Path

  import numpy as np
  import pandas as pd
  import pytest

  import biotapy as bt

  TESTS = Path(__file__).parents[1]
  pytestmark = pytest.mark.golden


  CASES = {
      "renorm_relab": ("pathabundance.tsv", "relab", True),
      "renorm_cpm": ("pathabundance.tsv", "cpm", True),
      "renorm_relab_nospecial": ("pathabundance.tsv", "relab", False),
      "renorm_cpm_genefamilies": ("genefamilies.tsv", "cpm", True),
  }


  @pytest.mark.parametrize("golden_name", CASES)
  def test_renorm_matches_humann_renorm_table(golden_name):
      source, units, special = CASES[golden_name]
      golden = pd.read_csv(TESTS / "golden" / "humann" / f"{golden_name}.csv.gz", index_col="sample_id")
      out = bt.fn.renorm(bt.io.read_humann(TESTS / "data" / "humann" / source), units, special=special)
      mods = [out[key] for key in ("function", "function_by_taxon")]
      table = pd.concat([pd.DataFrame(m.X.toarray(), index=m.obs_names, columns=m.var_names) for m in mods], axis=1)
      assert sorted(table.columns) == sorted(golden.columns)
      # Looser than the 1e-7 default: humann_renorm_table prints six significant digits (%.6g).
      np.testing.assert_allclose(table[golden.columns].to_numpy(), golden.to_numpy(), rtol=5e-6)
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/fn/test_renorm.py tests/fn/test_renorm_golden.py -q`
  -> `AttributeError: module 'biotapy.fn' has no attribute 'renorm'`.
- [x] **Step 3: Implement** - `src/biotapy/fn/_renorm.py`:
  ```python
  """Renormalisation of HUMAnN tables, with humann_renorm_table's community semantics."""

  from typing import Literal, cast

  import numpy as np
  from anndata import AnnData
  from mudata import MuData

  from biotapy._core import BY_TAXON_KEY, FUNCTION_KEY, SPECIAL_FEATURES, add_provenance, as_csr, feature_subset

  _SCALE = {"relab": 1.0, "cpm": 1e6}
  _X_KIND: dict[str, Literal["relative", "cpm"]] = {"relab": "relative", "cpm": "cpm"}


  def renorm(mdata: MuData, units: Literal["relab", "cpm"], *, special: bool = True) -> MuData:
      """Rescale both modalities so each sample's community total is 1 (or one million).

      Parameters
      ----------
      mdata
          ``bt.io.read_humann``'s result: modalities ``"function"`` and
          ``"function_by_taxon"``.
      units
          ``"relab"`` (totals 1) or ``"cpm"`` (totals 1,000,000).
      special
          Keep ``UNMAPPED``, ``READS_UNMAPPED``, ``UNINTEGRATED`` and
          ``UNGROUPED`` (community and per-taxon rows) and count them in the
          total; ``False`` drops them first.

      Returns
      -------
      MuData
          A copy whose two function modalities have ``X`` divided by each
          sample's total over the ``"function"`` modality, and
          ``uns['biotapy']['x_kind']`` set to ``"relative"`` or ``"cpm"``.
          Their ``layers``, ``obsm``, ``obsp``, ``varm`` and ``varp`` are
          dropped, as after any change to ``X``'s features. Other modalities
          and the global ``obs`` are copied unchanged.

      Raises
      ------
      KeyError
          ``mdata`` lacks one of the two modalities.
      ValueError
          ``units`` is not ``"relab"`` or ``"cpm"``; the ``"function"``
          modality has no feature.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/function`

      Matches ``humann_renorm_table`` (HUMAnN 3.9) in its default community
      mode: stratified rows are divided by the community total, so a
      pathway's strata need not sum to its community value. For HUMAnN's
      ``--mode levelwise``, use ``bt.pp.relative`` on each modality. A sample
      with a zero total stays zero.

      References
      ----------
      Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of
      diverse microbial communities with bioBakery 3. eLife 10:e65088.

      Examples
      --------
      >>> import biotapy as bt
      >>> out = bt.fn.renorm(bt.datasets.toy_humann(), "relab")
      >>> round(float(out["function"].X[0].sum()), 6)
      1.0
      """
      if units not in _SCALE:
          msg = f"units={units!r} must be 'relab' or 'cpm'"
          raise ValueError(msg)
      community, by_taxon = (_kept(mod, special=special) for mod in _function_modalities(mdata))
      if community.n_vars == 0:
          msg = f"mdata[{FUNCTION_KEY!r}] has no feature to total; renorm needs the community (unstratified) rows"
          raise ValueError(msg)
      totals = np.asarray(as_csr(community.X).sum(axis=1), dtype=np.float64).ravel()
      out = mdata.copy()
      # ModDict is a dict; mudata types the property as a read-only Mapping.
      mods = cast("dict[str, AnnData | MuData]", out.mod)
      for key, mod in ((FUNCTION_KEY, community), (BY_TAXON_KEY, by_taxon)):
          mods[key] = _rescaled(mod, totals, units=units, special=special)
      out.update()
      return out


  def _function_modalities(mdata: MuData) -> tuple[AnnData, AnnData]:
      community, by_taxon = mdata.mod.get(FUNCTION_KEY), mdata.mod.get(BY_TAXON_KEY)
      if not isinstance(community, AnnData) or not isinstance(by_taxon, AnnData):
          msg = f"mdata needs AnnData modalities {[FUNCTION_KEY, BY_TAXON_KEY]}, as bt.io.read_humann returns them"
          raise KeyError(msg)
      return community, by_taxon


  def _kept(adata: AnnData, *, special: bool) -> AnnData:
      """``adata`` without its special rows unless ``special``; a real object, never a view."""
      function = adata.var["function"] if "function" in adata.var.columns else adata.var_names.to_series()
      keep = np.ones(adata.n_vars, dtype=bool) if special else ~function.isin(SPECIAL_FEATURES).to_numpy()
      return feature_subset(adata, np.flatnonzero(keep))


  def _rescaled(adata: AnnData, totals: np.ndarray, *, units: str, special: bool) -> AnnData:
      """Divide each stored value by its sample's total (not by a reciprocal, which overflows for tiny totals)."""
      X = as_csr(adata.X).astype(np.float64)
      row_totals = np.repeat(totals, np.diff(X.indptr))
      X.data = np.divide(X.data, row_totals, out=np.zeros_like(X.data), where=row_totals > 0) * _SCALE[units]
      adata.X = X
      adata.uns["biotapy"]["x_kind"] = _X_KIND[units]
      add_provenance(adata, "fn.renorm", units=units, special=special)
      return adata
  ```
  `src/biotapy/fn/__init__.py`: add `from ._renorm import renorm`;
  `__all__ = ["func_glom", "load_hierarchy", "renorm"]`. `pyproject.toml`:
  `untyped_calls_exclude = [ "biom", "mudata", "skbio.stats._subsample", "sklearn", "threadpoolctl" ]`,
  and the comment above it names "mudata's MuData.update" beside
  threadpoolctl's `threadpool_limits`.
- [x] **Step 4: Run, expect pass** - same command -> `15 passed` (11 + 4
  golden). `uv run --group dev mypy` -> no issues.
- [x] **Step 5: Docs.** `docs/api.md`: `fn.renorm` after `fn.load_hierarchy`.
  Append to `docs/guide/function.md`:
  ````markdown
  ## Renormalising

  `bt.fn.renorm(mdata, "relab")` (or `"cpm"`) divides every row of both
  modalities by the sample's community total, as `humann_renorm_table` does by
  default. Stratified rows are divided by the community total too, so a
  pathway's strata keep their share of the community. `special=False` drops
  the special rows before totalling, as `--special n` does.

  HUMAnN's `--mode levelwise`, where each modality is scaled by its own total,
  is `bt.pp.relative` applied to each modality.
  ````
  Build docs -> `build succeeded.`
- [x] **Step 6: Knowledge** - `.knowledge/contracts/data-model-slots.md`:
  Propagation table, "Feature-changing" row: add `fn.renorm` ("rescales `X`
  and may drop the special rows; sets `x_kind` to `relative` or `cpm`").
  Convention 2: "`fn.renorm`'s `relative` stratified rows do not sum to 1:
  they are shares of the community total."
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  git add src/biotapy/fn/_renorm.py src/biotapy/fn/__init__.py tests/fn/test_renorm.py tests/fn/test_renorm_golden.py \
    pyproject.toml docs/api.md docs/guide/function.md .knowledge/contracts/data-model-slots.md \
    .knowledge/roadmap/phase-2-function.md .knowledge/log.md
  git commit -m "feat(fn): add renorm with humann_renorm_table community semantics"
  ```

### Checkpoint A - review slice 2A
- [ ] Review the whole slice (superpowers:requesting-code-review) against
  every contract, the pure-by-default decision, no-bundled-kegg, and the
  Phase 2 and slice 2A review focus; then a fix pass, one commit per finding,
  each with a test.
- [ ] Run the exit-gate check now: `uv run --group test pytest -m golden tests/fn -q`
  (all HUMAnN goldens pass) and the full `uv run --group test pytest`.
- [x] Knowledge (codebase-map templates; R12.2-R12.4):
  - **Create `.knowledge/modules/fn.md`** (`type: Module`, `paths:
    ["src/biotapy/fn/**"]`).
    - **Responsibility:** function hierarchies, aggregation and
      renormalisation; owns no reader and no download.
    - **Entry points:** `load_hierarchy`, `func_glom`, `renorm` (+
      `_glom.py:_pairs`, the membership builder).
    - **Invariants:**
      - HUMAnN 3.9 semantics (design notes 1-2);
      - edge table columns and transitive rows;
      - `UNGROUPED`/protected rules;
      - outputs through `replace_features` / `feature_subset`;
      - hierarchy `attrs["source"]` copied into provenance;
      - no network code.
    - **Dependencies:** `_core` (`sum_pairs`, `replace_features`,
      `feature_subset`, the function constants), mudata.
    - **Verification:** `uv run --group test pytest tests/fn -q`.
    - **Gotchas:**
      - many-to-many inflates totals;
      - `mean` counts members present;
      - `READS_UNMAPPED` is unknown to HUMAnN 3.9;
      - `renorm` replaces `X`;
      - the `str` dtype needed for h5.
  - **Create `.knowledge/decisions/function-tables-as-mudata.md`**
    (`type: Decision`). Context: pathway strata are not additive. Options
    weighed (research C section 1.4): MuData, `dict` of AnnData, `uns`,
    `obsm`, one mixed AnnData. Decision: MuData with the two modalities, in
    Phase 2 rather than Phase 4. Consequences: h5mu drops a TreeData
    modality's tree (no function modality has one), and mudata joins the
    dependencies.
  - **Update `.knowledge/modules/io.md`** (`read_humann`, `_humann.py`
    helpers, the header rule), **`datasets.md`** (`toy_humann`, `enzyme`,
    the unpinned ENZYME hash) and **`core.md`** (`_function.py`,
    `sum_pairs`, `replace_features`). Add `fn.md` to `modules/index.md` and
    the decision to `decisions/index.md`, with log lines.
- [ ] Push the branch and open the PR, after the user approves that push
  (R13.3). CI green, including docs and the network job.
- [ ] Ask the user to review slice 2A before slice 2B.

---
## Slice 2B - Other readers


**Goal:** a MetaPhlAn user reads a profile, or a merged table, into a
`TreeData` whose rows sum to 1 with nothing counted twice. A PICRUSt2 user
reads a prediction and its contributions into the same two-modality `MuData`
as HUMAnN, and reads the per-ASV copy numbers that 2.8 needs.

**How this slice was checked.** Every file below was written into a scratch
clone of the repository at `2d47f8d`, the `phase-2a` tip. The clone was then
replayed as one commit per task (2.2a, 2.2, 2.4, 2.4b), and each committed
state was gated:

| Task state | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `sphinx-build -W` |
|---|---|---|---|
| 2.2a | passed | 752 passed, 22 deselected | (no docs change) |
| 2.2 | passed | 784 passed, 22 deselected | build succeeded |
| 2.4 | passed | 810 passed, 22 deselected | build succeeded |
| 2.4b | passed | 818 passed, 22 deselected | build succeeded |

- prek covers ruff 0.16.9 check and format, `mypy --strict`, import-linter
  and pyproject-fmt.
- Coverage on the final state: `_table.py`, `_humann.py`, `_metaphlan.py` and
  `_picrust2.py` are each 100%; overall 99%.
- The new property tests also passed under Hypothesis seeds 1, 2 and 3.
- The RED counts in each Step 2 were reproduced by unexporting the new
  readers.
- Real data: HMP2's `taxonomic_profiles_3.tsv.gz` (MetaPhlAn, merged, 933 rows
  x 1,638 samples) reads to 1,638 x 579 in 0.5 s. Its leaf sums run
  0.9999993-1.0000008. All 1,638 sample names equal those `read_humann` gives
  for HMP2's `pathabundances_3.tsv.gz`.
- APIs were checked in the installed versions: pandas 3.0.6, numpy 2.5.3,
  scipy 1.18.1, anndata 0.13.4, mudata 0.4.1, treedata 0.3.1.

### Slice 2B design

- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `io/_table.py` (new, 2.2a) | `_read_table`, `_numbers`: the strict TSV reading every function-table reader shares |
  | `io/_humann.py` (2.2a) | `read_humann`, now on `_table.py`; its own `_read_table` is removed |
  | `io/_metaphlan.py` (2.2) | `read_metaphlan`, `_header` |
  | `io/_picrust2.py` (2.4, 2.4b) | `read_picrust2`, `_contributions`, `read_picrust2_traits`, `_first_line` |
  | `_core/__init__.py` (2.2) | exports `RELATIVE_TOLERANCE` (already defined in `_slots.py`) |

- **Shared parsing (brief question 7).** The three readers need the same
  strict checks: an empty file, a row longer than the header, a gap, a
  non-number, and every message naming the file. Those checks are
  `_humann.py:_read_table` today, so they move to `io/_table.py` rather than
  being copied (R4.3). They stay in `io`, not `_core`, because one
  subpackage uses them: module-boundaries rule 2 sends a helper to `_core`
  only when two subpackages use it. `io/_join.py` is the precedent, a private
  topic file shared by `io` readers. The split into `_read_table` (structure)
  and `_numbers` (values) exists for two readers:
  - a MetaPhlAn profile has text columns it does not read (taxids, `-`
    coverage, `additional_species`);
  - a PICRUSt2 contribution table has three id columns.

  Both parse only the columns they use as numbers.

  `_read_table` gains one check, for repeated column names. pandas renames a
  second `S1` to `S1.1` without a word (measured: `read_humann` read two
  samples, `S1` and `S1.1`).

- **MetaPhlAn leaf rule (brief question 1).** A row becomes a feature when no
  other row descends from it through any ancestor. That is the deepest row of
  each lineage: SGBs (`t__`) in MetaPhlAn 4, species in MetaPhlAn 3 and HMP2,
  or the one rank of a `--tax_lev` table. Reasons:
  - **No double counting, nothing lost.** MetaPhlAn computes a clade's
    coverage as the sum of its children's (`metaphlan.py:compute_coverage`)
    and each row's abundance as coverage over the leaves' total. So every
    internal row is the sum of its leaves, and the leaves alone hold every
    read once.
  - **Parity with MetaPhlAn.** `pp.tax_glom` of the leaves to any rank gives
    back MetaPhlAn's own rows for that rank. The fixture test checks this at
    all seven ranks.
  - **"Any ancestor", not "direct parent".** The real 4.0.6 fixture has no
    `o__Corynebacteriales` row. With a direct-parent rule,
    `c__Actinomycetia` would be a leaf and 52.86% would be counted twice.
  - **"Only `t__` rows" (the roadmap's first idea) fails** on MetaPhlAn 3 and
    HMP2, which have no `t__` rows.
  - **No `rank=` parameter (R2.3).** `pp.tax_glom` already gives any higher
    rank, and equals MetaPhlAn's rows when it does.

  The guard is a check, not a warning. Each sample's leaves (and
  `UNCLASSIFIED`) must sum to 1 within `_core.RELATIVE_TOLERANCE` (1e-3, the
  data model's own definition of `relative`); otherwise the reader raises a
  `ValueError` naming the samples. The check catches:
  - rows removed by hand;
  - lineages the rule cannot read (GTDB-style `;`, which would sum to
    n_ranks x 100%);
  - tables that are not profiles (marker tables).

  It replaces the outline's proposed warning (R7.4: a table whose label would
  be false must not load). All-zero samples are allowed, as `infer_x_kind`
  ignores them.

- **Feature ids and `var["sgb"]` (brief question 2).**
  - `var_names` are the leaf's last name without its prefix: `SGB1871`,
    `Bacteroides_ovatus`, `UNCLASSIFIED`. Repeats raise, naming `path`.
  - Ranks `kingdom`..`species` come from `_core.split_lineage`. It already
    ignores `t__`, which sits at position 7, past `RANKS`.
  - **No `sgb` column (R2.3).** The SGB id is already the feature id.
    MetaPhlAn 2 and 3 used `t__` for strains, not SGBs, so a column named
    `sgb` would mislabel them. Nothing in 0.2 reads it: 2.10 uses MetaPhlAn
    3 species.

- **Formats (brief question 3).**

  | Input | Header | Samples | Handled by |
  |---|---|---|---|
  | Per-sample profile, `-t rel_ab` (3.x, 4.x) | last `#` line, `#clade_name ... relative_abundance ...` | one, named after the file | the `relative_abundance` column |
  | `-t rel_ab_w_read_stats` | the same, with `clade_taxid`, `coverage`, `estimated_number_of_reads_from_the_clade` | one | the same column; `-` coverage never parsed |
  | Merged by `merge_metaphlan_tables.py` (4.x) | a `#mpa_v...` line (no tab), then `clade_name<TAB>samples` | columns | first line after the `#` lines |
  | Merged MetaPhlAn 3 | the same, with `NCBI_tax_id` second | columns | `NCBI_tax_id` / `clade_taxid` dropped |
  | HMP2 `taxonomic_profiles_3` | no `#`; `Feature\Sample` corner; CRLF; `UNKNOWN` row | columns | same rule |

  - **The header rule** is HUMAnN's (the last `#` line), except that a `#`
    line with no tab is never the header. That excludes the database line.
  - **Sample names:** `_profile` is removed, as `merge_metaphlan_tables.py`
    removes it from file names. So a profile read alone and the same profile
    in a merged table get one name. HMP2's taxonomic and functional tables
    then share all 1,638 names, which 2.10 needs.
  - **Empty taxids** are never read. **Short rows** (before 4.2.3 MetaPhlAn
    could drop the empty `additional_species` cell, CHANGELOG 4.2.3) are
    fine; only the columns read must be complete.
  - **Values:** percentages divided by 100, so `x_kind="relative"` (design
    note 3).
  - **`UNCLASSIFIED` and `UNKNOWN`** stay features with every rank NaN.
    `pp.tax_glom` drops them unless `dropna=False`; the guide says so.

- **PICRUSt2's `EC:` prefix (brief question 4). Strip it.** PICRUSt2 writes
  `EC:1.1.1.1`; ENZYME, HUMAnN and `bt.datasets.enzyme` write `1.1.1.1`.
  - Keeping it would make `func_glom` against ENZYME raise its nothing-maps
    error (good), but the fix would fall on the user. They would have to
    rewrite `var_names` and `var["function"]` in two modalities. Stripping a
    hierarchy's `child` column is one line on one table.
  - It is stripped in both modalities and in the trait table's columns. KO
    and pathway ids are untouched.
  - No provenance parameter records the strip (R2.3); the docstring and the
    guide state it.

- **PICRUSt2 inputs and signature (brief question 5).**
  `read_picrust2(path, *, contrib=None) -> MuData`.
  - `path` is the unstratified table (`pred_metagenome_unstrat`,
    `path_abun_unstrat`).
  - `contrib` is its long contribution table (`pred_metagenome_contrib`,
    `path_abun_contrib`).
  - **One call takes both files.** PICRUSt2 writes the community and
    stratified rows to two files, which are one table split in two.
    `read_humann`'s one-table-per-call rule was about not mixing a gene
    family table with a pathway table in one modality, and still holds: a
    `contrib` naming a function or sample that `path` lacks raises.
  - **Accepted:** unstratified tables alone (the per-taxon modality then has 0
    features, like a HUMAnN table without strata), and the long contribution
    format, located by the column names `sample`, `function`, `taxon`,
    `taxon_function_abun`. It therefore takes both the 9-column gene-family
    file and the pathway file, whose column count the sources disagree on (8
    or 9).
  - **Rejected with a `ValueError`:**
    - the legacy wide stratified table (`--wide_table`; deprecated by
      PICRUSt2 as memory-hungry);
    - an `add_descriptions.py` `description` column;
    - the 4-column pathway coverage contributions.
  - **Dropped from the outline: `taxa=` (`seqtab_norm`, R2.3).** 2.8 takes any
    samples x ASVs `AnnData`, such as the user's ASV table through
    `read_biom`, `read_dada2` or `read_qiime2`. Whether Tian's `p` should be
    16S-copy-corrected is a 2.8 question; if it should, 2.8 adds the reader.

- **`x_kind` and additivity (brief question 6).** Both modalities are
  `"abundance"`, so `pp.rarefy` refuses them.
  - **Gene families are additive.** Per research B section 1.4, read from
    PICRUSt2's source (not re-read here: GPL): a gene family's unstratified
    value is sum over ASVs of `seqtab_norm x copies`, and a contribution's
    `taxon_function_abun` is `taxon_abun x genome_function_count` on the same
    rounded `seqtab_norm`. So a sample's contributions sum to its
    unstratified value, up to float rounding.
  - **Pathways are not.** Pathway abundances come from MinPath/HUMAnN-style
    key-reaction means (FAQ), as in HUMAnN.
  - **So the reader checks nothing about additivity.** It reads both tables
    as written, and the test pins that pathway strata are kept unchanged.

- **`RARE`.** PICRUSt2 groups ASVs under `--min_reads`/`--min_samples` as
  taxon `RARE` (seen in its own test data's `taxon` column). It is an
  ordinary stratum, not `special`. Marking it special, as the outline
  proposed, would make `fn.renorm(special=False)` drop real abundance.

- **`read_picrust2_traits(path) -> pd.DataFrame`** (2.4b) reads ASVs x
  functions copy numbers. Ids are read as text (PICRUSt2's own test data uses
  `2593338844` and `2568526487.0`), `EC:` is stripped, and `metadata_NSTI`
  (the last column when `hsp.py -n` ran) is dropped. Repeated ASV ids raise.
  It is a DataFrame because its rows are genomes, not samples, so it cannot
  be a modality.

- **Facts the tasks rely on** (checked 2026-10-03; re-check any API you call,
  R2.2).
  - **MetaPhlAn** (`metaphlan.py` and `merge_metaphlan_tables.py` from
    MetaPhlAn master, MIT, read in research A's copy):
    - `relative_abundances()` divides each clade's coverage by the leaves'
      total, and internal coverage is the children's sum;
    - `UNCLASSIFIED` is `(1 - fraction_mapped) x 100`, first, with taxid
      `-1`;
    - values are rounded to 5-7 decimals;
    - `merge_metaphlan_tables.py` writes the first input's database line,
      then `clade_name<TAB>samples`, keeps only `relative_abundance`, and
      names samples `splitext(basename)[0].replace('_profile', '')`.
  - **The fixture** `demo_metaphlan_bugs_list.tsv`: HUMAnN e07b3a3
    `humann/tests/data/`, 3,679 bytes, SHA-256
    `5b5ea7d9a86e481346c2fd9622a04ff172d0d91e9b8ceeb3f9df1f57b959fbc3`. It is
    `-t rel_ab_w_read_stats` with five `#` lines and has no
    `o__Corynebacteriales` row.
  - **PICRUSt2** (header lines only, read from `picrust/picrust2` master
    `tests/test_data/`; no data copied):

    | File | Header |
    |---|---|
    | `pred_metagenome_unstrat` | `function<TAB>samples` |
    | `path_abun_unstrat` | `pathway<TAB>samples` |
    | contributions | `sample, function, taxon, taxon_abun, taxon_rel_abun, genome_function_count, taxon_function_abun, taxon_rel_function_abun, norm_taxon_function_contrib` |
    | pathway coverage contributions | 4 columns: `sample, function, taxon, genome_function_count` |
    | `hsp.py` output | `sequence<TAB>EC:...<TAB>metadata_NSTI` |
    | `add_descriptions.py` output | `function<TAB>description<TAB>samples` |
    | legacy wide | `function<TAB>sequence<TAB>samples` |

    The `taxon` column holds `RARE` in `expected_metagenome_contrib_rare`.
  - **pandas 3.0.6:**
    - `read_csv` renames a repeated header cell (`S1` -> `S1.1`);
    - with `usecols`, it silently drops the extra cells of a longer row,
      which is why `_read_table` reads every column;
    - its default float parser keeps about 15 significant digits.
  - **anndata 0.13.4** turns `str` columns into categoricals when writing
    h5mu, in place. Round-trip tests compare values, not dtypes.

### Slice 2B global constraints (in addition to the Phase 2 list)
- No new dependency. `mia` stays out of the R image (user, 2026-10-03), so the
  readers get no R golden: parity comes from invariants and the tools' own
  output (r-golden-parity statement 8).
- PICRUSt2 (GPL-3) is never installed, imported or copied. Its test tables
  are synthetic strings in `tests/io/test_picrust2.py`, written from the
  documented headers, and no file is committed for it. The MetaPhlAn fixture
  is the one copied file, with its `NOTICE.txt`.
- Every error a user can hit names the file through its argument
  (`path='...'` or `contrib='...'`). That includes errors raised in `_core`,
  which are re-raised with the argument and the original as `__cause__`.
- Tests call `bt.io.*` only; `_table.py` is tested through the readers
  (R4.9).
- Diff blocks below are for reading, not `git apply`: their blank context
  lines have lost their leading space, which prek's trailing-whitespace hook
  would strip anyway.
- Run commands: `uv run --group test pytest <path> -q`. Gate before every
  commit: `uvx prek run --all-files`. With docs changes, also run
  `BIOTAPY_DATA_DIR=<scratchpad>/pooch uv run --group doc sphinx-build -W -b html docs docs/_build/html`.
  The tutorials download phyloseq datasets: this plan's own check wrote
  `~/.cache/biotapy` once without the variable, and that directory was
  removed.
- Branch: `phase-2b` from `master` once PR #15 has merged, else from
  `phase-2a`'s tip. Push, PR and merge-commit on green are approved for Phase
  2 slice branches (user, 2026-10-03).

### Slice 2B review focus
The five ways real users are most likely to get a wrong answer from these
readers without an error, each pinned by a test in its task:

1. **Abundance counted twice, or lost**, from internal rows (including a rank
   MetaPhlAn skipped) or from a table the leaf rule cannot read. Expected:
   leaves only; a sample whose leaves do not sum to 100% raises. Tests (2.2):
   - `test_reads_the_leaf_clades_of_a_metaphlan_4_profile` (the skipped
     order);
   - `test_glom_to_each_rank_matches_metaphlans_own_rows`;
   - `test_leaves_that_do_not_sum_to_100_raise_naming_the_path`.
2. **Percentages taken as counts, or `UNCLASSIFIED` dropped**, so samples sum
   to less than 1. Expected: `x_kind == "relative"` and `UNCLASSIFIED` kept
   with no rank. Tests (2.2): `test_reads_the_leaf_clades_of_a_metaphlan_4_profile`,
   `test_unclassified_is_kept_with_every_rank_missing`.
3. **PICRUSt2's `EC:` ids never match a hierarchy.** Expected: ids stripped,
   and `func_glom` against an ENZYME-style edge table maps. Test (2.4):
   `test_ec_ids_match_enzyme_style_hierarchies`.
4. **Mismatched PICRUSt2 files** (an EC unstratified table with a KO
   contribution table, a pathway table with gene contributions, or a
   different sample set) silently give a stratified modality that does not
   belong to the community one. Expected: a `ValueError` naming `contrib` and
   up to three unknown ids. Test (2.4):
   `test_malformed_inputs_raise_naming_the_file[unknown-function]` and
   `[unknown-sample]`.
5. **Ids rewritten by parsing**: numeric-looking ASV ids (`0042`,
   `2593338844`), a repeated sample header, or MetaPhlAn's `_profile` suffix
   breaking the join with HUMAnN sample names. Expected: ids as written,
   repeats raise, `_profile` removed. Tests:
   - 2.4 `test_var_matches_the_humann_layout` (`0042`);
   - 2.4b `test_traits_are_asvs_by_functions`;
   - 2.2a `test_repeated_sample_columns_raise_naming_the_path`;
   - 2.2 `test_reads_an_hmp2_style_table`.

Execution order: **2.2a -> 2.2 -> 2.4 -> 2.4b -> Checkpoint B.** 2.2a comes
first because both new readers import `io/_table.py`. 2.4b extends 2.4's
file.

---

### Task 2.2a: shared strict table reading in `io`

**Files:** create `src/biotapy/io/_table.py`; modify
`src/biotapy/io/_humann.py`, `tests/io/test_humann.py`.

**Interfaces:**
- Consumes nothing new.
- Produces:
  - `_read_table(path: Path, header: str, *, skiprows: int, argument: str, text: int = 1) -> pd.DataFrame`
    (indexed by the first column; raises `ValueError` starting with
    `argument`);
  - `_numbers(frame: pd.DataFrame, *, argument: str) -> np.ndarray`
    (float64, same shape as `frame`).

  Both are private to `io`; 2.2, 2.4 and 2.4b import them with
  `from ._table import _numbers, _read_table`.

**Will not touch:** `read_humann`'s header rule, `_SUFFIX`, `_UNITS`, its
messages other than the shared ones, or any other reader (BIOM, QIIME 2,
DADA2 and phyloseq parse differently and are not changed, R1.4).

- [x] **Step 1: Failing tests.** Append to `tests/io/test_humann.py`:

```python
def test_repeated_sample_columns_raise_naming_the_path(tmp_path):
    # pandas would rename the second "S1" to "S1.1" and read two samples.
    path = tmp_path / "twice.tsv"
    path.write_text("# Pathway\tS1\tS1\nPWY-1\t1.0\t2.0\n")
    with pytest.raises(ValueError, match=r"twice\.tsv.*repeats column names \['S1'\]"):
        bt.io.read_humann(path)


def test_row_without_an_id_raises_naming_the_path(tmp_path):
    path = tmp_path / "no_id.tsv"
    path.write_text("# Pathway\tS1\nPWY-1\t1.0\n\t2.0\n")
    with pytest.raises(ValueError, match=r"no_id\.tsv.*no id"):
        bt.io.read_humann(path)
```

- [x] **Step 2: Run, expect failure.** `uv run --group test pytest tests/io/test_humann.py -q`
  -> `2 failed, 29 passed`:
  - `test_repeated_sample_columns_raise_naming_the_path` fails with
    `DID NOT RAISE`;
  - the no-id test fails on its message, because today's text is `var ids
    must not be missing`.
- [x] **Step 3: Implement.** Create `src/biotapy/io/_table.py`:

```python
"""Strict reading of the tab-separated tables HUMAnN, MetaPhlAn and PICRUSt2 write."""

import csv
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


def _read_table(path: Path, header: str, *, skiprows: int, argument: str, text: int = 1) -> pd.DataFrame:
    """The table whose header line is ``header``, indexed by its first column.

    ``skiprows`` lines precede the header. The first ``text`` columns are read
    as text, so ids such as ``0042`` stay as written. ``argument`` names the
    input in messages, e.g. ``"path='table.tsv'"``. An empty file, repeated
    column names, a data row with more cells than the header or a row with no
    id raise ``ValueError``; a short row's missing cells are NaN, for the
    caller to check in the columns it reads (``_numbers``).
    """
    names = header.rstrip("\r\n").split("\t")
    repeated = sorted(name for name, count in Counter(names).items() if count > 1)
    if repeated:
        msg = f"{argument} repeats column names {repeated[:3]}"
        raise ValueError(msg)
    try:
        table = pd.read_csv(
            path,
            sep="\t",
            skiprows=skiprows,
            index_col=0,
            dtype=dict.fromkeys(range(text), str),
            quoting=csv.QUOTE_NONE,
        )
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
        msg = f"{argument} is not a valid tab-separated table: {error}"
        raise ValueError(msg) from error
    # pandas shifts the header over when the first data row is longer, so compare with the header's own cells.
    if table.shape[1] != len(names) - 1:
        msg = f"{argument} has a data row with more cells than the header"
        raise ValueError(msg)
    if np.any(table.index.isna()):
        msg = f"{argument} has a data row with no id (an empty first cell)"
        raise ValueError(msg)
    return table


def _numbers(frame: pd.DataFrame, *, argument: str) -> np.ndarray:
    """``frame``'s values as a float64 array, raising ``ValueError`` naming ``argument`` on a non-number or a gap."""
    try:
        values = frame.to_numpy(dtype=np.float64)
    except ValueError as error:
        msg = f"{argument} has a value that is not a number: {error}"
        raise ValueError(msg) from error
    if np.isnan(values).any():
        msg = f"{argument} has a missing or NaN value (a data row with fewer cells than the header, or an empty cell)"
        raise ValueError(msg)
    return values
```

  Then change `src/biotapy/io/_humann.py` exactly as this diff does:
  - imports;
  - one `argument` string for every message;
  - the two calls;
  - the `Raises` text;
  - the old `_read_table` deleted.

```diff
diff --git a/src/biotapy/io/_humann.py b/src/biotapy/io/_humann.py
index af5b6f4..6cec334 100644
--- a/src/biotapy/io/_humann.py
+++ b/src/biotapy/io/_humann.py
@@ -1,16 +1,16 @@
 """HUMAnN 3 and 4 tables: gene families, reactions, pathway abundance, and their regrouped or renormalised forms."""

-import csv
 import gzip
 import re
 from pathlib import Path

-import numpy as np
 import pandas as pd
 from mudata import MuData

 from biotapy._core import XKind, make_function_mudata

+from ._table import _numbers, _read_table
+
 # Sample-column suffixes: HUMAnN's own ("_Abundance-RPKs", "_Abundance"), renorm --update-snames'
 # ("-CPM", "-RELAB"), and the file names humann_join_tables uses when every file names its
 # sample alike ("<sample>_pathabundance_cpm", as in the HMP2 merged tables).
@@ -56,11 +56,11 @@ def read_humann(path: str | Path) -> MuData:
     Raises
     ------
     ValueError
-        The file is empty or a pathway coverage table; a value is not a
-        number; a data row has more cells than the header, fewer cells, or an
-        empty cell (any missing value); a row id holds more than one ``|``;
-        or sample names repeat once their suffix is removed. Messages name
-        ``path``.
+        The file is empty or a pathway coverage table; the header repeats a
+        column name; a value is not a number; a data row has more cells than
+        the header, fewer cells, an empty cell (any missing value) or no id;
+        a row id holds more than one ``|``; or sample names repeat once their
+        suffix is removed. Messages name ``path``.

     Notes
     -----
@@ -101,6 +101,7 @@ def read_humann(path: str | Path) -> MuData:
     (['UNMAPPED', 'K1'], ['K1|g__A.s__A_b'])
     """
     path = Path(path)
+    argument = f"path={str(path)!r}"
     # HUMAnN's rule: the last "#" line is the header; with none, the first line is.
     header, n_comments = "", 0
     opener = gzip.open if path.suffix == ".gz" else open
@@ -111,9 +112,10 @@ def read_humann(path: str | Path) -> MuData:
                 break
             header, n_comments = line, n_comments + 1
     if _COVERAGE.search(header):
-        msg = f"path={str(path)!r} is a pathway coverage table; read_humann reads abundance tables"
+        msg = f"{argument} is a pathway coverage table; read_humann reads abundance tables"
         raise ValueError(msg)
-    table, X = _read_table(path, header, max(n_comments - 1, 0))
+    table = _read_table(path, header, skiprows=max(n_comments - 1, 0), argument=argument)
+    X = _numbers(table, argument=argument).T
     obs = pd.DataFrame(index=table.columns.str.replace(_SUFFIX, "", regex=True))
     # HUMAnN never writes raw counts: a table whose header names no unit holds pathway abundances.
     unit = next((kind for pattern, kind in _UNITS if pattern.search(header.rstrip("\n"))), None)
@@ -121,27 +123,5 @@ def read_humann(path: str | Path) -> MuData:
     try:
         return make_function_mudata(X, obs=obs, row_ids=table.index, x_kind=x_kind, source="io.read_humann")
     except ValueError as error:  # repeated sample or row ids, or a row id with two "|"
-        msg = f"path={str(path)!r}: {error}"
-        raise ValueError(msg) from error
-
-
-def _read_table(path: Path, header: str, skiprows: int) -> tuple[pd.DataFrame, np.ndarray]:
-    """Read the table and its values, raising ``ValueError`` that names ``path`` on a malformed file."""
-    try:
-        table = pd.read_csv(path, sep="\t", skiprows=skiprows, index_col=0, dtype={0: str}, quoting=csv.QUOTE_NONE)
-    except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
-        msg = f"path={str(path)!r} is not a valid HUMAnN table: {error}"
+        msg = f"{argument}: {error}"
         raise ValueError(msg) from error
-    # pandas shifts the header over when the first data row is longer, so compare with the header's own cells.
-    if table.shape[1] != header.rstrip("\n").count("\t"):
-        msg = f"path={str(path)!r} has a data row with more cells than the header"
-        raise ValueError(msg)
-    try:
-        values = table.to_numpy(dtype=np.float64)
-    except ValueError as error:
-        msg = f"path={str(path)!r} has a value that is not a number: {error}"
-        raise ValueError(msg) from error
-    if np.isnan(values).any():
-        msg = f"path={str(path)!r} has a missing or NaN value (a data row with fewer cells than the header, or an empty cell)"
-        raise ValueError(msg)
-    return table, values.T
```

- [x] **Step 4: Run, expect pass.** The same command gives `31 passed`. The
  doctest `uv run --group test pytest src/biotapy/io/_humann.py -q` gives `1 passed`.
- [x] **Step 5: Bookkeeping.**
  - Tick this task's boxes and its line under "# Tasks (checklist)".
  - Add to `.knowledge/log.md`, under a new
    `## <date> (Phase 2, slice 2B)` heading at the top:
    `- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.2a done: io/_table.py holds the strict TSV reading shared by the HUMAnN, MetaPhlAn and PICRUSt2 readers; read_humann now rejects repeated column names.`
  - No contract changes. `modules/io.md` is refreshed at Checkpoint B.
- [x] **Step 6: Gate and commit.**

```bash
uvx prek run --all-files
git add src/biotapy/io/_table.py src/biotapy/io/_humann.py tests/io/test_humann.py \
  .knowledge/roadmap/phase-2-function.md .knowledge/log.md
git commit -m "refactor(io): share read_humann's strict table reading; reject repeated column names"
uv run --group test pytest -q -W error::UserWarning   # 752 passed, 22 deselected
```

### Task 2.2: `io.read_metaphlan`

**Files:**
- Create:
  - `src/biotapy/io/_metaphlan.py`;
  - `tests/io/test_metaphlan.py`;
  - `tests/data/metaphlan/demo_metaphlan_bugs_list.tsv`;
  - `tests/data/metaphlan/NOTICE.txt`.
- Modify:
  - `src/biotapy/io/__init__.py`, `src/biotapy/_core/__init__.py`;
  - `docs/api.md`, `docs/guide/reading_data.md`;
  - `.knowledge/contracts/data-model-slots.md`,
    `.knowledge/contracts/r-golden-parity.md`.

**Interfaces:**
- Consumes:
  - 2.2a's `_read_table` and `_numbers`;
  - `_core.make_treedata(X, *, obs, var, tree, x_kind, source)`,
    `_core.split_lineage(lineage: pd.Series) -> pd.DataFrame`, `_core.RANKS`;
  - `_core.RELATIVE_TOLERANCE` (exported here).
- Produces: `bt.io.read_metaphlan(path: str | Path) -> TreeData`. 2.10 reads
  HMP2's taxonomic table with it.

**Will not touch:** `split_lineage` and `normalize_ranks` (they already
ignore `t__`), `pp.tax_glom`'s `dropna` behaviour, and `RANKS` (no `strain`
or `sgb` rank).

- [x] **Step 1: Fixture and notice.** Copy the file from HUMAnN's repository
  at the pinned commit and check its hash:

```bash
mkdir -p tests/data/metaphlan
curl -sfL https://raw.githubusercontent.com/biobakery/humann/e07b3a3/humann/tests/data/demo_metaphlan_bugs_list.tsv \
  -o tests/data/metaphlan/demo_metaphlan_bugs_list.tsv
sha256sum tests/data/metaphlan/demo_metaphlan_bugs_list.tsv
# 5b5ea7d9a86e481346c2fd9622a04ff172d0d91e9b8ceeb3f9df1f57b959fbc3
```

  Create `tests/data/metaphlan/NOTICE.txt`:

```text
demo_metaphlan_bugs_list.tsv is copied unchanged from HUMAnN's test data,
https://github.com/biobakery/humann/tree/e07b3a3/humann/tests/data
(commit e07b3a3, 2026-07-10). It is a MetaPhlAn 4.0.6 profile
(-t rel_ab_w_read_stats, database mpa_vOct22_CHOCOPhlAnSGB_202212), used under
HUMAnN's MIT licence:

The HUMAnN software is licensed under the MIT license.

Copyright (c) 2014 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

The other MetaPhlAn tables in tests/io/test_metaphlan.py are synthetic,
written for biotapy (BSD-3-Clause) from the formats in MetaPhlAn's source.
```

- [x] **Step 2: Failing tests.** Create `tests/io/test_metaphlan.py`:

```python
import gzip
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import treedata as td
from hypothesis import given
from hypothesis import strategies as st

import biotapy as bt
from biotapy._core import RANKS

# tests/data/metaphlan/NOTICE.txt: a real MetaPhlAn 4.0.6 profile from HUMAnN's test data (MIT).
DEMO = Path(__file__).parents[1] / "data" / "metaphlan" / "demo_metaphlan_bugs_list.tsv"
LINEAGE = "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae"
# A default MetaPhlAn 4.2 profile (synthetic): four "#" lines, the header, UNCLASSIFIED first.
PROFILE_4_2 = (
    "#mpa_vJan25_CHOCOPhlAnSGB_202503\n#metaphlan s1.fastq -o s1_profile.tsv\n#1000 reads processed\n"
    "#SampleID\tMetaphlan_Analysis\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
    "UNCLASSIFIED\t-1\t20.0\t\n"
    "k__Bacteria\t2\t80.0\t\n"
    f"{LINEAGE}\t2|976|200643|171549|815\t80.0\t\n"
    f"{LINEAGE}|g__Bacteroides\t2|976|200643|171549|815|816\t80.0\t\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus\t2|976|200643|171549|815|816|28116\t50.0\t\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus|t__SGB1871\t2|976|200643|171549|815|816|28116|\t50.0\t\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_SGB1\t2|976|200643|171549|815|816|\t30.0\ts__Bacteroides_x\n"
    f"{LINEAGE}|g__Bacteroides|s__Bacteroides_SGB1|t__SGB1\t2|976|200643|171549|815|816||\t30.0\t\n"
)


def write(tmp_path, text, name="table.tsv"):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_reads_the_leaf_clades_of_a_metaphlan_4_profile():
    # The file skips o__Corynebacteriales; its parent c__Actinomycetia must still not count as a leaf.
    tdata = bt.io.read_metaphlan(DEMO)
    assert tdata.obs_names.tolist() == ["demo_metaphlan_bugs_list"]
    assert tdata.var_names.tolist() == ["SGB2091", "SGB1871", "SGB2301", "SGB1814"]
    np.testing.assert_allclose(tdata.X.toarray(), [[0.5285687, 0.346816, 0.1214543, 0.0031609]], rtol=1e-12)
    assert tdata.uns["biotapy"]["x_kind"] == "relative"


def test_var_holds_the_ranks_of_each_leaf():
    var = bt.io.read_metaphlan(DEMO).var
    assert var.columns.tolist() == list(RANKS)
    assert var.loc["SGB1871"].tolist() == [
        "Bacteria",
        "Bacteroidetes",
        "Bacteroidia",
        "Bacteroidales",
        "Bacteroidaceae",
        "Bacteroides",
        "Bacteroides_ovatus",
    ]
    assert var.loc["SGB2091", "species"] == "Corynebacterium_SGB2091"


@pytest.mark.parametrize("rank", RANKS)
def test_glom_to_each_rank_matches_metaphlans_own_rows(rank):
    # MetaPhlAn's parity check (contracts/r-golden-parity, statement 8): a clade's row is the sum of its leaves.
    rows = pd.read_csv(DEMO, sep="\t", skiprows=5, index_col=0)["relative_abundance"]
    last = rows.index.str.split("|").str[-1]
    at_rank = last.str.startswith(f"{rank[0]}__")
    expected = (rows[at_rank] / 100).set_axis(last[at_rank].str[3:])
    out = bt.pp.tax_glom(bt.io.read_metaphlan(DEMO), rank)
    # MetaPhlAn 4.0.6 printed no o__Corynebacteriales row, so compare the clades it printed.
    got = pd.Series(out.X.toarray()[0], index=out.var[rank].to_numpy())[expected.index]
    # MetaPhlAn prints each percentage rounded to 5 decimals, so a parent and the sum of its leaves differ by up to
    # (number of leaves) x 5e-6 % (1.1e-5 % measured for Bacteroidaceae).
    pd.testing.assert_series_equal(got.sort_index(), expected.sort_index(), check_names=False, rtol=0, atol=1e-6)


def test_unclassified_is_kept_with_every_rank_missing(tmp_path):
    tdata = bt.io.read_metaphlan(write(tmp_path, PROFILE_4_2, "s1_profile.tsv"))
    assert tdata.obs_names.tolist() == ["s1"]
    assert tdata.var_names.tolist() == ["UNCLASSIFIED", "SGB1871", "SGB1"]
    assert tdata.var.loc["UNCLASSIFIED"].isna().all()
    np.testing.assert_allclose(tdata.X.toarray(), [[0.2, 0.5, 0.3]])


def test_reads_a_merged_metaphlan_4_table(tmp_path):
    text = (
        "#mpa_vJan25_CHOCOPhlAnSGB_202503\nclade_name\ts1\ts2_profile\n"
        "UNCLASSIFIED\t20.0\t0.0\nk__Bacteria\t80.0\t100.0\n"
        f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus\t80.0\t100.0\n"
        f"{LINEAGE}|g__Bacteroides|s__Bacteroides_ovatus|t__SGB1871\t80.0\t100.0\n"
    )
    tdata = bt.io.read_metaphlan(write(tmp_path, text))
    assert tdata.obs_names.tolist() == ["s1", "s2"]
    np.testing.assert_allclose(tdata.X.toarray(), [[0.2, 0.8], [0.0, 1.0]])


def test_reads_a_merged_metaphlan_3_table_and_drops_its_taxids(tmp_path):
    # MetaPhlAn 3 leaves are species; a leaf may stop above species, which is then NaN.
    text = (
        "#mpa_v30_CHOCOPhlAn_201901\nclade_name\tNCBI_tax_id\tA\tB\n"
        "k__Bacteria\t2\t100.0\t100.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae\t2|1239|186801|186802|186803\t100.0\t100.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Roseburia\t2|1239|186801|186802|186803|841\t40.0\t0.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia\t2|1239|186801|186802|186803|572511\t60.0\t100.0\n"
        "k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Lachnospiraceae|g__Blautia|s__Blautia_obeum\t2|1239|186801|186802|186803|572511|40520\t60.0\t100.0\n"
    )
    tdata = bt.io.read_metaphlan(write(tmp_path, text))
    assert tdata.obs_names.tolist() == ["A", "B"]
    assert tdata.var_names.tolist() == ["Roseburia", "Blautia_obeum"]
    assert pd.isna(tdata.var.loc["Roseburia", "species"]) and tdata.var.loc["Roseburia", "genus"] == "Roseburia"


def test_reads_an_hmp2_style_table(tmp_path):
    # HMP2's taxonomic_profiles_3.tsv.gz: no "#" line, a "Feature\Sample" corner, UNKNOWN, CRLF line ends.
    text = (
        "Feature\\Sample\tCSM5FZ3N_P_profile\tCSM5FZ4M_profile\r\n"
        "UNKNOWN\t0\t100\r\nk__Bacteria\t100\t0\r\nk__Bacteria|p__Firmicutes\t100\t0\r\n"
    )
    tdata = bt.io.read_metaphlan(write(tmp_path, text, "taxonomic_profiles_3.tsv"))
    assert tdata.obs_names.tolist() == ["CSM5FZ3N_P", "CSM5FZ4M"]
    assert tdata.var_names.tolist() == ["UNKNOWN", "Firmicutes"]
    np.testing.assert_array_equal(tdata.X.toarray(), [[0.0, 1.0], [1.0, 0.0]])


def test_a_single_rank_table_reads_that_rank(tmp_path):
    # metaphlan --tax_lev s writes species rows only.
    text = "#mpa_v\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
    text += f"{LINEAGE}|g__B|s__B_a\t\t70.0\t\n{LINEAGE}|g__B|s__B_b\t\t30.0\t\n"
    assert bt.io.read_metaphlan(write(tmp_path, text)).var_names.tolist() == ["B_a", "B_b"]


def test_reads_read_stats_profiles_and_short_rows(tmp_path):
    # -t rel_ab_w_read_stats has a "-" coverage for UNCLASSIFIED; before 4.2.3 some rows lacked additional_species.
    text = (
        "#mpa_v\n#clade_name\tclade_taxid\trelative_abundance\tcoverage\testimated_number_of_reads_from_the_clade\n"
        "UNCLASSIFIED\t-1\t10.0\t-\t50\nk__Bacteria\t2\t90.0\n"
    )
    np.testing.assert_allclose(bt.io.read_metaphlan(write(tmp_path, text)).X.toarray(), [[0.1, 0.9]])


def test_reads_gzip(tmp_path):
    path = tmp_path / "demo_profile.tsv.gz"
    path.write_bytes(gzip.compress(DEMO.read_bytes()))
    tdata = bt.io.read_metaphlan(path)
    assert tdata.obs_names.tolist() == ["demo"] and tdata.n_vars == 4


def test_all_zero_sample_and_feature_are_kept(tmp_path):
    text = "clade_name\tA\tB\nk__Bacteria\t100.0\t0.0\nk__Bacteria|p__F\t100.0\t0.0\nk__Bacteria|p__G\t0.0\t0.0\n"
    tdata = bt.io.read_metaphlan(write(tmp_path, text))
    assert tdata.shape == (2, 2) and tdata.X[1].nnz == 0
    assert tdata.var_names.tolist() == ["F", "G"]


def test_an_empty_profile_has_no_features(tmp_path):
    text = "#mpa_v\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
    assert bt.io.read_metaphlan(write(tmp_path, text, "empty_profile.tsv")).shape == (1, 0)


@pytest.mark.parametrize(
    ("rows", "total"),
    [
        (f"k__Bacteria\t100.0\n{LINEAGE}|g__B\t60.0\n", 60.0),  # a leaf row removed
        ("d__Bacteria\t100.0\nd__Bacteria;p__Firmicutes\t100.0\n", 200.0),  # GTDB-style ";" lineages
    ],
    ids=["rows-removed", "semicolon-lineages"],
)
def test_leaves_that_do_not_sum_to_100_raise_naming_the_path(tmp_path, rows, total):
    path = write(tmp_path, "clade_name\tS1\n" + rows, "partial.tsv")
    with pytest.raises(ValueError, match=rf"partial\.tsv.*1 sample\(s\) do not sum to 100%: \{{'S1': {total}\}}"):
        bt.io.read_metaphlan(path)


@pytest.mark.parametrize(
    "text",
    [
        "clade_name\tS1\nk__Bacteria\tabc\n",
        "#mpa\n#clade_name\tNCBI_tax_id\trelative_abundance\nk__Bacteria\t2\t\n",
        "clade_name\tS1\nk__Bacteria\t50.0\t50.0\n",
        "clade_name\tS1\tS1\nk__Bacteria\t100.0\t100.0\n",
        "",
    ],
    ids=["not-a-number", "missing-value", "longer-row", "repeated-sample", "empty-file"],
)
def test_malformed_tables_raise_naming_the_path(tmp_path, text):
    with pytest.raises(ValueError, match=r"bad\.tsv"):
        bt.io.read_metaphlan(write(tmp_path, text, "bad.tsv"))


def test_repeated_leaf_names_raise_naming_the_path(tmp_path):
    text = "clade_name\tS1\nk__A|g__X\t50.0\nk__B|g__X\t50.0\n"
    with pytest.raises(ValueError, match=r"twice\.tsv.*duplicate var ids: \['X'\]"):
        bt.io.read_metaphlan(write(tmp_path, text, "twice.tsv"))


def test_round_trips_through_h5td(tmp_path):
    tdata = bt.io.read_metaphlan(DEMO)
    tdata.write_h5td(tmp_path / "m.h5td")
    back = td.read_h5td(tmp_path / "m.h5td")
    pd.testing.assert_frame_equal(back.var, tdata.var)


@given(st.lists(st.floats(0.01, 100, allow_nan=False), min_size=1, max_size=6))
def test_leaves_keep_every_percentage(weights):
    leaves = [round(100 * weight / sum(weights), 5) for weight in weights]
    clades = {f"{LINEAGE}|g__G{i % 2}|s__G{i % 2}_s{i}|t__SGB{i}": value for i, value in enumerate(leaves)}
    rows = dict(clades)
    for clade, value in clades.items():
        parts = clade.split("|")
        for depth in range(1, len(parts)):
            rows["|".join(parts[:depth])] = rows.get("|".join(parts[:depth]), 0) + value
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "t.tsv"
        path.write_text("clade_name\tS1\n" + "".join(f"{clade}\t{value!r}\n" for clade, value in rows.items()))
        tdata = bt.io.read_metaphlan(path)
    assert tdata.var_names.tolist() == [f"SGB{i}" for i in range(len(leaves))]
    np.testing.assert_allclose(tdata.X.toarray().ravel(), np.array(leaves) / 100, rtol=1e-12)
```

- [x] **Step 3: Run, expect failure.** `uv run --group test pytest tests/io/test_metaphlan.py -q`
  -> `28 failed`: `AttributeError: module 'biotapy.io' has no attribute 'read_metaphlan'`.
- [x] **Step 4: Implement.**
  - In `src/biotapy/_core/__init__.py`, import `RELATIVE_TOLERANCE` from
    `._slots`, before `XKind`, and add `"RELATIVE_TOLERANCE"` to `__all__`
    after `"RANKS"`:

```diff
diff --git a/src/biotapy/_core/__init__.py b/src/biotapy/_core/__init__.py
index dec9722..c466b8f 100644
--- a/src/biotapy/_core/__init__.py
+++ b/src/biotapy/_core/__init__.py
@@ -13,6 +13,7 @@ from ._matrix import argmax_by, as_csr, sum_by, sum_pairs
 from ._optional import import_optional
 from ._rng import as_generator
 from ._slots import (
+    RELATIVE_TOLERANCE,
     XKind,
     add_provenance,
     feature_subset,
@@ -43,6 +44,7 @@ __all__ = [
     "PHYLO_KEY",
     "PROTECTED_FEATURES",
     "RANKS",
+    "RELATIVE_TOLERANCE",
     "SPECIAL_FEATURES",
     "UNGROUPED",
     "TreeData",
```

  - Create `src/biotapy/io/_metaphlan.py`. `_header` is a single-use helper
    that R4.4 allows: inlined, `read_metaphlan` has 36 statements, over
    `PLR0915`'s 30. The leaf rule is inlined, because it fits.

```python
"""MetaPhlAn 3 and 4 taxonomic profiles, per sample or merged by merge_metaphlan_tables.py."""

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from biotapy._core import RELATIVE_TOLERANCE, TreeData, make_treedata, split_lineage

from ._table import _numbers, _read_table

# NCBI taxid columns: a MetaPhlAn 3 merged table keeps one beside its samples.
_TAXID_COLUMNS = ("NCBI_tax_id", "clade_taxid")


def read_metaphlan(path: str | Path) -> TreeData:
    r"""Read a MetaPhlAn profile, or several merged, into a samples x clades table.

    Parameters
    ----------
    path
        A MetaPhlAn 3 or 4 profile (``-t rel_ab``, the default, or
        ``-t rel_ab_w_read_stats``), or a table of several merged by
        ``merge_metaphlan_tables.py`` or with one column per sample. Gzip
        (``.gz``) is read directly.

    Returns
    -------
    TreeData
        Relative abundances in ``X`` (MetaPhlAn's percentages divided by
        100), one feature per leaf clade: the deepest row of each lineage,
        such as MetaPhlAn 4's SGBs (``t__SGB1871``) or MetaPhlAn 3's species.
        ``var_names`` are the leaf's last name without its rank prefix
        (``SGB1871``, ``Bacteroides_ovatus``); ``var`` holds the rank columns
        ``kingdom`` to ``species``. ``UNCLASSIFIED`` (``UNKNOWN`` in older
        tables) stays a feature with every rank NaN, so each sample sums to 1.
        A single profile's sample is named after its file, and ``_profile``
        is removed from sample names, as ``merge_metaphlan_tables.py`` does.
        There is no tree.

    Raises
    ------
    ValueError
        The file is empty or malformed (a value that is not a number, a
        missing value, a data row with more cells than the header, repeated
        column names); the leaf clades of a sample do not sum to 100 (rows
        removed, or a table that is not a profile, such as one with ``;``
        lineages); or leaf names or sample names repeat. Messages name
        ``path``.

    Notes
    -----
    R equivalent: ``mia::importMetaPhlAn``
    Guide: :doc:`/guide/reading_data`

    A profile lists every rank, and a clade's abundance is the sum of its
    children's, so keeping only the leaves keeps all the abundance once.
    ``bt.pp.tax_glom`` gives back the higher ranks; it drops
    ``UNCLASSIFIED`` unless ``dropna=False``. ``uns['biotapy']['x_kind']`` is
    ``"relative"``. NCBI taxids, ``additional_species``, coverage and read
    estimates are not read.

    References
    ----------
    Blanco-Míguez A et al. (2023) Extending and improving metagenomic taxonomic profiling with
    uncharacterized species using MetaPhlAn 4. Nature Biotechnology 41:1633-1644.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> path = Path(tempfile.mkdtemp()) / "S1_profile.tsv"
    >>> rows = ["k__Bacteria\t2\t90.0\t", "k__Bacteria|g__Bacteroides\t2|816\t90.0\t", "UNCLASSIFIED\t-1\t10.0\t"]
    >>> _ = path.write_text(
    ...     "#mpa_vJan25\n#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n" + "\n".join(rows)
    ... )
    >>> tdata = bt.io.read_metaphlan(path)
    >>> tdata.obs_names.tolist(), tdata.var_names.tolist(), tdata.X.toarray().tolist()
    (['S1'], ['Bacteroides', 'UNCLASSIFIED'], [[0.9, 0.1]])
    """
    path = Path(path)
    argument = f"path={str(path)!r}"
    header, skiprows = _header(path)
    table = _read_table(path, header, skiprows=skiprows, argument=argument)
    if "relative_abundance" in table.columns:
        values = table[["relative_abundance"]].set_axis([Path(path.name.removesuffix(".gz")).stem], axis=1)
    else:
        values = table.drop(columns=[column for column in _TAXID_COLUMNS if column in table.columns])
    # A leaf is a clade no other clade descends from, through any ancestor: MetaPhlAn can omit an
    # intermediate rank's row (the 4.0.6 fixture has no o__Corynebacteriales), so direct parents are not enough.
    lineages = [clade.split("|") for clade in table.index]
    ancestors = {"|".join(lineage[:depth]) for lineage in lineages for depth in range(1, len(lineage))}
    leaf = np.array([clade not in ancestors for clade in table.index], dtype=bool)
    X = _numbers(values, argument=argument)[leaf].T / 100
    totals = X.sum(axis=1)
    bad = (totals != 0) & (np.abs(totals - 1) > RELATIVE_TOLERANCE)
    if bad.any():
        shown = {
            name: round(float(total) * 100, 3)
            for name, total in zip(values.columns[bad][:3], totals[bad][:3], strict=True)
        }
        msg = (
            f"{argument}: the leaf clades of {int(bad.sum())} sample(s) do not sum to 100%: {shown}; "
            "read_metaphlan reads whole MetaPhlAn profiles with '|'-separated lineages"
        )
        raise ValueError(msg)
    clades = table.index[leaf]
    # A clade lineage is "k__A|p__B|..."; UNCLASSIFIED and UNKNOWN have none, so every rank is NaN.
    lineage = pd.Series(clades.str.replace("|", ";"), index=clades).where(clades.str.contains("__", regex=False))
    var = split_lineage(lineage).set_axis(clades.str.split("|").str[-1].str.replace(r"^[a-z]__", "", regex=True))
    obs = pd.DataFrame(index=values.columns.str.replace("_profile", "", regex=False))
    try:
        return make_treedata(X, obs=obs, var=var, tree=None, x_kind="relative", source="io.read_metaphlan")
    except ValueError as error:  # repeated leaf or sample names
        msg = f"{argument}: {error}"
        raise ValueError(msg) from error


def _header(path: Path) -> tuple[str, int]:
    """The header line and how many lines precede it.

    A profile's header is its last ``#`` line (``#clade_name...``); a merged
    table's is the first line after its ``#`` lines, the first of which names
    the database and holds no tab.
    """
    comments: list[str] = []
    header = ""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                header = line
                break
            comments.append(line)
    if comments and "\t" in comments[-1]:
        return comments[-1], len(comments) - 1
    return header, len(comments)
```

  - `src/biotapy/io/__init__.py` becomes:

```python
from ._biom import read_biom, write_biom
from ._dada2 import read_dada2
from ._humann import read_humann
from ._metaphlan import read_metaphlan
from ._phyloseq import read_phyloseq
from ._qiime2 import read_qiime2

__all__ = ["read_biom", "read_dada2", "read_humann", "read_metaphlan", "read_phyloseq", "read_qiime2", "write_biom"]
```

- [x] **Step 5: Run, expect pass.** The same command gives `28 passed`. The
  doctest `uv run --group test pytest src/biotapy/io/_metaphlan.py -q` gives
  `1 passed`.
- [x] **Step 6: Docs.**
  - `docs/api.md`: add `io.read_metaphlan` after `io.read_humann`.
  - `docs/guide/reading_data.md`: replace the introduction's second
    paragraph, which says every reader infers `x_kind` (false since
    `read_humann`), and append the MetaPhlAn section:

````diff
diff --git a/docs/guide/reading_data.md b/docs/guide/reading_data.md
index 1c08141..fdeab8b 100644
--- a/docs/guide/reading_data.md
+++ b/docs/guide/reading_data.md
@@ -3,11 +3,13 @@
 Readers turn a file format into a `TreeData` that follows the one
 [data model](data_model.md) every biotapy function relies on.

-No format below records whether its table holds counts or proportions, so
-every reader infers `uns["biotapy"]["x_kind"]` from the values: whole numbers
-are `"counts"`; otherwise, if every sample with a nonzero total sums to 1
-(within `1e-3`), `"relative"`; anything else is `"abundance"`. Functions
-that need raw counts check this and refuse proportions.
+BIOM, QIIME 2, DADA2 and phyloseq tables do not record whether they hold
+counts or proportions, so their readers infer `uns["biotapy"]["x_kind"]` from
+the values: whole numbers are `"counts"`; otherwise, if every sample with a
+nonzero total sums to 1 (within `1e-3`), `"relative"`; anything else is
+`"abundance"`. The HUMAnN, MetaPhlAn and PICRUSt2 readers take it from the
+format instead, as their sections say. Functions that need raw counts check
+`x_kind` and refuse anything else.

 ## Example datasets

@@ -282,3 +284,37 @@ mdata = bt.io.read_humann("sample_genefamilies.tsv")
   `--update-snames` so the header names the new unit.
 - **Not read.** Pathway coverage tables (HUMAnN 3 only) raise a
   `ValueError`: they are not abundances.
+
+## MetaPhlAn
+
+`bt.io.read_metaphlan` reads a MetaPhlAn 3 or 4 profile, or a table of
+several merged by `merge_metaphlan_tables.py`, into a `TreeData` of relative
+abundances.
+
+```python
+import biotapy as bt
+
+tdata = bt.io.read_metaphlan("merged_abundance_table.tsv")
+```
+
+- **Leaf clades.** A profile lists every rank, kingdom to SGB, and a clade's
+  abundance is the sum of its children's. The reader keeps the leaves - the
+  rows no other row descends from: SGBs (`t__SGB1871`) in MetaPhlAn 4,
+  species in MetaPhlAn 3 - so every read is counted once. Use
+  `bt.pp.tax_glom` for higher ranks; its sums equal MetaPhlAn's own rows.
+- **Ids and ranks.** A leaf's last name without its prefix becomes the
+  feature id (`SGB1871`, `Bacteroides_ovatus`), and its lineage fills
+  `kingdom` to `species` in `var`.
+- **Units.** Percentages are divided by 100, so `x_kind` is `"relative"`.
+- **`UNCLASSIFIED`.** MetaPhlAn 4.2 estimates the share of reads from
+  unknown organisms (`UNKNOWN` in older tables). It stays a feature with no
+  rank, so each sample sums to 1. `bt.pp.tax_glom` drops it unless you pass
+  `dropna=False`.
+- **Samples.** A single profile's sample is named after its file;
+  `_profile` is removed from sample names, as `merge_metaphlan_tables.py`
+  does when it names columns.
+- **Checked.** Each sample's leaves must sum to 100%. A table with rows
+  removed, or one that is not a MetaPhlAn profile (GTDB-style `;` lineages,
+  marker tables), raises a `ValueError` instead of double counting.
+- **Not read.** NCBI taxids, `additional_species`, coverage and estimated
+  read counts.
````

  - Build: `BIOTAPY_DATA_DIR=<scratchpad>/pooch uv run --group doc sphinx-build -W -b html docs docs/_build/html`
    -> `build succeeded.`
- [x] **Step 7: Contracts** (user-approved as decisions 1, 2 and 7).
  - `data-model-slots.md`: one sentence in convention 2, and a new section
    "Taxonomic profiles (MetaPhlAn)" before "Propagation".
  - `r-golden-parity.md`: statements 6 and 8.

  The exact text (ignore the `docs/api.md` hunk, which is Step 6's):

```diff
diff --git a/.knowledge/contracts/data-model-slots.md b/.knowledge/contracts/data-model-slots.md
index 1e2ab00..3785146 100644
--- a/.knowledge/contracts/data-model-slots.md
+++ b/.knowledge/contracts/data-model-slots.md
@@ -59,6 +59,8 @@ Extends the spec's data-model table with exact keys.[^spec]
    The exception is `io.read_humann`, which reads it from the table header (`RPKs` ->
    `rpk`; `CPM`, `_cpm` or `Adjusted CPMs` -> `cpm`; `RELAB`, `_relab` ->
    `relative`) and labels a header without a unit `abundance`, never `counts`.
+   `io.read_metaphlan` divides MetaPhlAn's percentages by 100 and sets `relative`, after
+   checking that every sample's leaf clades sum to 1 within `1e-3` (`_core.RELATIVE_TOLERANCE`).
    `fn.renorm` rescales `X` (and may drop the special rows), setting `x_kind` to `relative` or `cpm`;
    its `relative` stratified rows do not sum to 1: they are shares of the community total.
    The other readers infer it because their formats record no unit (BIOM, QIIME 2 `RelativeFrequency`, a DADA2
@@ -100,6 +102,16 @@ columns use the pandas `str` dtype, as rank columns do. Both modalities
 always exist; either may have 0 features. The community values are not
 the sum of their strata for pathways, which is why there are two.

+## Taxonomic profiles (MetaPhlAn)
+`io.read_metaphlan` keeps one feature per leaf clade: a row that no other
+row descends from through any ancestor (MetaPhlAn can omit an intermediate
+rank's row; the real 4.0.6 fixture has no `o__Corynebacteriales`). A clade's
+row is the sum of its leaves', so every read counts once and `pp.tax_glom`
+rebuilds the higher ranks. `var_names` are the leaf's last name without its
+rank prefix (`SGB1871` from `t__SGB1871`); `t__` has no rank column. Rank
+columns run `kingdom` to `species`. `UNCLASSIFIED` (`UNKNOWN` in older
+tables) stays a feature with every rank NaN, so samples sum to 1.
+
 ## Propagation
 | Operation | Keeps | Drops |
 |---|---|---|
diff --git a/.knowledge/contracts/r-golden-parity.md b/.knowledge/contracts/r-golden-parity.md
index 324ce0d..05c6d53 100644
--- a/.knowledge/contracts/r-golden-parity.md
+++ b/.knowledge/contracts/r-golden-parity.md
@@ -55,7 +55,8 @@ sources:
    accepted this for this BSD-3 repository on 2026-09-27. Test fixtures under `tests/data/`
    stay synthetic, except small files copied under a permissive licence
    with a `NOTICE.txt` beside them (`tests/data/humann`: HUMAnN's MIT test
-   data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0).
+   data; `tests/data/metaphlan`: a MetaPhlAn 4.0.6 profile from HUMAnN's MIT
+   test data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0).
 7. `pl` functions have an R equivalent but no golden test. They draw numbers
    that `tl` stores, and `tl`'s golden tests check those numbers (controller
    ruling 2026-09-27; rules.md R11.2).
@@ -64,7 +65,11 @@ sources:
    `io.read_humann` (R equivalent `mia::importHUMAnN`) is checked this way:
    the HUMAnN golden tests of `fn.func_glom` and `fn.renorm` read their
    inputs through it and compare with HUMAnN's own output. mia is not added
-   to the R image (user-approved 2026-10-03).
+   to the R image (user-approved 2026-10-03). `io.read_metaphlan` (R equivalent
+   `mia::importMetaPhlAn`) is checked against MetaPhlAn's own output: its
+   leaves, grouped by `pp.tax_glom` to each rank, equal the clade rows the
+   profile prints (`tests/io/test_metaphlan.py`, atol `1e-6` because
+   MetaPhlAn rounds each percentage to 5 decimals).

 # Why
 R and NumPy random generators differ, so stochastic outputs can never match
diff --git a/docs/api.md b/docs/api.md
index 8da4273..4d5008a 100644
--- a/docs/api.md
+++ b/docs/api.md
@@ -14,6 +14,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
     io.read_biom
     io.read_dada2
     io.read_humann
+    io.read_metaphlan
     io.read_phyloseq
     io.read_qiime2
     io.write_biom
```

- [x] **Step 8: Bookkeeping.**
  - Tick this task's boxes and its checklist line.
  - In `data-model-slots.md` and `r-golden-parity.md`, set
    `generated: { by: claude-code/<your model id>, at: <UTC now> }` and
    `commit:` to `git rev-parse --short HEAD`. Their `description`s do not
    change.
  - Add to `.knowledge/log.md`:
    `- **Update**: [data-model-slots](contracts/data-model-slots.md) convention 2 (MetaPhlAn percentages / 100, checked to sum to 1) and new section "Taxonomic profiles (MetaPhlAn)"; [r-golden-parity](contracts/r-golden-parity.md) statement 6 lists tests/data/metaphlan and statement 8 gives read_metaphlan's parity (tax_glom equals MetaPhlAn's own rows); [phase-2-function](roadmap/phase-2-function.md) task 2.2 done.`
- [x] **Step 9: Gate and commit.**

```bash
uvx prek run --all-files
git add src/biotapy/io/_metaphlan.py src/biotapy/io/__init__.py src/biotapy/_core/__init__.py \
  tests/io/test_metaphlan.py tests/data/metaphlan/demo_metaphlan_bugs_list.tsv tests/data/metaphlan/NOTICE.txt \
  docs/api.md docs/guide/reading_data.md .knowledge/contracts/data-model-slots.md \
  .knowledge/contracts/r-golden-parity.md .knowledge/roadmap/phase-2-function.md .knowledge/log.md
git commit -m "feat(io): read MetaPhlAn 3 and 4 profiles and merged tables as leaf clades"
uv run --group test pytest -q -W error::UserWarning   # 784 passed, 22 deselected
```

### Task 2.4: `io.read_picrust2`

**Files:**
- Create: `src/biotapy/io/_picrust2.py`, `tests/io/test_picrust2.py`.
- Modify:
  - `src/biotapy/io/__init__.py`;
  - `docs/api.md`, `docs/guide/reading_data.md`, `docs/guide/function.md`;
  - `.knowledge/contracts/data-model-slots.md`,
    `.knowledge/contracts/r-golden-parity.md`.

**Interfaces:**
- Consumes:
  - 2.2a's `_read_table(..., text=3)` and `_numbers`;
  - `_core.make_function_mudata(X, *, obs, row_ids, x_kind, source) -> MuData`,
    which splits `function|taxon` ids into the two modalities.
- Produces:
  - `bt.io.read_picrust2(path: str | Path, *, contrib: str | Path | None = None) -> MuData`.
    Modalities `"function"` and `"function_by_taxon"`; `var` exactly as
    `read_humann`'s; `taxon` = ASV id or `RARE`.
  - `fn.contributions` (2.7) and `pl.contributions` (2.9) read this layout.

**Will not touch:** `_core/_function.py` (`function_var` already gives NaN
genus and species for ASV ids), `SPECIAL_FEATURES` (`RARE` is not special),
and `fn.func_glom`.

- [x] **Step 1: Failing tests.** Create `tests/io/test_picrust2.py`:

```python
import gzip
import tempfile
from pathlib import Path

import mudata
import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

import biotapy as bt

# Synthetic PICRUSt2 2.6 outputs (PICRUSt2 is GPL-3: nothing is copied from it). Column names follow its
# documented headers. Abundances: S1 holds ASV1 (10) and 0042 (4), S2 holds 0042 (5.5) and RARE (2), S3 is empty.
# Copy numbers: ASV1 has EC:1.1.1.1 x1 and EC:2.7.1.1 x2; 0042 has EC:2.7.1.1 x1 and EC:3.2.1.1 x1; RARE has
# EC:1.1.1.1 x1. Each unstratified value is the sum of its contributions; EC:4.1.1.1 is predicted nowhere.
UNSTRAT = (
    "function\tS1\tS2\tS3\n"
    "EC:1.1.1.1\t10.0\t2.0\t0.0\n"
    "EC:2.7.1.1\t24.0\t5.5\t0.0\n"
    "EC:3.2.1.1\t4.0\t5.5\t0.0\n"
    "EC:4.1.1.1\t0.0\t0.0\t0.0\n"
)
CONTRIB_COLUMNS = (
    "sample\tfunction\ttaxon\ttaxon_abun\ttaxon_rel_abun\tgenome_function_count"
    "\ttaxon_function_abun\ttaxon_rel_function_abun\tnorm_taxon_function_contrib\n"
)
CONTRIB = CONTRIB_COLUMNS + (
    "S1\tEC:1.1.1.1\tASV1\t10.0\t71.43\t1\t10.0\t71.43\t1.0\n"
    "S1\tEC:2.7.1.1\tASV1\t10.0\t71.43\t2\t20.0\t142.86\t0.833\n"
    "S1\tEC:2.7.1.1\t0042\t4.0\t28.57\t1\t4.0\t28.57\t0.167\n"
    "S1\tEC:3.2.1.1\t0042\t4.0\t28.57\t1\t4.0\t28.57\t1.0\n"
    "S2\tEC:1.1.1.1\tRARE\t2.0\t26.67\t1\t2.0\t26.67\t1.0\n"
    "S2\tEC:2.7.1.1\t0042\t5.5\t73.33\t1\t5.5\t73.33\t1.0\n"
    "S2\tEC:3.2.1.1\t0042\t5.5\t73.33\t1\t5.5\t73.33\t1.0\n"
)
# Pathways (8 columns): a pathway's contributions need not sum to its unstratified value.
PATH_UNSTRAT = "pathway\tS1\tS2\nPWY-1\t3.5\t1.0\nPWY-2\t0.5\t2.0\n"
PATH_CONTRIB = (
    "sample\tfunction\ttaxon\ttaxon_abun\ttaxon_rel_abun\tgenome_function_count\ttaxon_function_abun\ttaxon_rel_function_abun\n"
    "S1\tPWY-1\tASV1\t10.0\t71.43\t0.5\t5.0\t35.71\n"
    "S1\tPWY-1\t0042\t4.0\t28.57\t0.25\t1.0\t7.14\n"
    "S2\tPWY-2\t0042\t5.5\t73.33\t0.5\t2.75\t36.67\n"
)


def write(tmp_path, text, name):
    path = tmp_path / name
    path.write_text(text)
    return path


@pytest.fixture
def mdata(tmp_path):
    unstrat = write(tmp_path, UNSTRAT, "pred_metagenome_unstrat.tsv")
    return bt.io.read_picrust2(unstrat, contrib=write(tmp_path, CONTRIB, "pred_metagenome_contrib.tsv"))


def test_reads_community_and_per_taxon_modalities(mdata):
    function, by_taxon = mdata["function"], mdata["function_by_taxon"]
    assert mdata.obs_names.tolist() == ["S1", "S2", "S3"]
    assert function.var_names.tolist() == ["1.1.1.1", "2.7.1.1", "3.2.1.1", "4.1.1.1"]
    assert by_taxon.var_names.tolist() == [
        "1.1.1.1|ASV1",
        "2.7.1.1|ASV1",
        "2.7.1.1|0042",
        "3.2.1.1|0042",
        "1.1.1.1|RARE",
    ]
    np.testing.assert_array_equal(by_taxon.X.toarray()[1], [0.0, 0.0, 5.5, 5.5, 2.0])


def test_var_matches_the_humann_layout(mdata):
    assert mdata["function"].var.columns.tolist() == ["name", "special"]
    var = mdata["function_by_taxon"].var
    assert var.columns.tolist() == ["function", "name", "taxon", "genus", "species", "special"]
    assert var["taxon"].tolist() == ["ASV1", "ASV1", "0042", "0042", "RARE"]
    assert var[["name", "genus", "species"]].isna().all().all() and not var["special"].any()


def test_x_kind_is_abundance_so_rarefy_refuses_it(mdata):
    assert mdata["function"].uns["biotapy"]["x_kind"] == "abundance"
    assert mdata["function_by_taxon"].uns["biotapy"]["x_kind"] == "abundance"
    with pytest.raises(ValueError, match="raw counts"):
        bt.pp.rarefy(mdata["function"], depth=5, seed=0)


def test_gene_family_contributions_sum_to_the_community_values(mdata):
    by_taxon = mdata["function_by_taxon"]
    summed = pd.DataFrame(by_taxon.X.toarray(), columns=by_taxon.var["function"]).T.groupby(level=0).sum().T
    community = pd.DataFrame(mdata["function"].X.toarray(), columns=mdata["function"].var_names)
    pd.testing.assert_frame_equal(summed, community[summed.columns], check_names=False)


def test_pathway_contributions_are_kept_as_written(tmp_path):
    unstrat = write(tmp_path, PATH_UNSTRAT, "path_abun_unstrat.tsv")
    mdata = bt.io.read_picrust2(unstrat, contrib=write(tmp_path, PATH_CONTRIB, "path_abun_contrib.tsv"))
    np.testing.assert_array_equal(mdata["function"].X.toarray(), [[3.5, 0.5], [1.0, 2.0]])
    np.testing.assert_array_equal(mdata["function_by_taxon"].X.toarray(), [[5.0, 1.0, 0.0], [0.0, 0.0, 2.75]])


def test_without_contrib_the_per_taxon_modality_is_empty(tmp_path):
    mdata = bt.io.read_picrust2(write(tmp_path, UNSTRAT, "unstrat.tsv"))
    assert mdata["function"].shape == (3, 4) and mdata["function_by_taxon"].shape == (3, 0)


def test_ec_ids_match_enzyme_style_hierarchies(mdata):
    edges = pd.DataFrame({"child": ["1.1.1.1", "2.7.1.1"], "parent": ["1.-.-.-", "2.-.-.-"], "level": "class"})
    out = bt.fn.func_glom(mdata["function"], "class", hierarchy=edges)
    assert out.var_names.tolist() == ["1.-.-.-", "2.-.-.-", "UNGROUPED"]


def test_ko_ids_are_kept(tmp_path):
    mdata = bt.io.read_picrust2(write(tmp_path, "function\tS1\nK00001\t1.5\n", "ko.tsv"))
    assert mdata["function"].var_names.tolist() == ["K00001"]


def test_all_zero_sample_and_feature_are_kept(mdata):
    assert mdata["function"].X[2].nnz == 0 and mdata["function_by_taxon"].X[2].nnz == 0
    assert mdata["function"].X[:, 3].nnz == 0


def test_single_sample(tmp_path):
    columns = CONTRIB_COLUMNS
    mdata = bt.io.read_picrust2(
        write(tmp_path, "function\tS1\nEC:1.1.1.1\t6.0\n", "u.tsv"),
        contrib=write(tmp_path, columns + "S1\tEC:1.1.1.1\tASV1\t3\t100\t2\t6.0\t200\t1\n", "c.tsv"),
    )
    assert mdata["function"].shape == (1, 1) and mdata["function_by_taxon"].shape == (1, 1)


def test_reads_gzip(tmp_path):
    unstrat, contrib = tmp_path / "u.tsv.gz", tmp_path / "c.tsv.gz"
    unstrat.write_bytes(gzip.compress(UNSTRAT.encode()))
    contrib.write_bytes(gzip.compress(CONTRIB.encode()))
    assert bt.io.read_picrust2(unstrat, contrib=contrib)["function_by_taxon"].n_vars == 5


def test_round_trips_through_h5mu(tmp_path, mdata):
    mdata.write_h5mu(tmp_path / "p.h5mu")
    back = mudata.read_h5mu(tmp_path / "p.h5mu")
    assert (back["function_by_taxon"].X != mdata["function_by_taxon"].X).nnz == 0
    assert back["function_by_taxon"].var["taxon"].tolist() == ["ASV1", "ASV1", "0042", "0042", "RARE"]


ROW = "\t1\t1\t1\t1.0\t1\t1\n"  # the six numeric contribution cells, taxon_function_abun = 1.0
# case -> (unstrat text, contrib text, the message expected)
MALFORMED = {
    "contrib-columns": (UNSTRAT, CONTRIB.replace("taxon_function_abun", "abun"), r"contrib=.*c\.tsv.*long-format"),
    "legacy-wide-contrib": (UNSTRAT, "function\tsequence\tS1\nEC:1.1.1.1\tASV1\t1.0\n", r"c\.tsv.*long-format"),
    "unknown-sample": (UNSTRAT, CONTRIB + "S9\tEC:1.1.1.1\tASV1" + ROW, r"c\.tsv.*samples that path lacks: \['S9'\]"),
    "unknown-function": (UNSTRAT, CONTRIB + "S1\tK00001\tASV1" + ROW, r"functions that path lacks: \['K00001'\]"),
    "repeated-row": (UNSTRAT, CONTRIB + "S1\tEC:1.1.1.1\tASV1" + ROW, r"c\.tsv.*repeats a sample, function and taxon"),
    "missing-taxon": (UNSTRAT, CONTRIB + "S1\tEC:1.1.1.1\t" + ROW, r"c\.tsv.*no function or taxon"),
    "non-number": (UNSTRAT, CONTRIB + "S1\tEC:4.1.1.1\tASV1\t1\t1\t1\tx\t1\t1\n", r"c\.tsv.*not a number"),
    "bar-in-taxon": (UNSTRAT, CONTRIB + "S1\tEC:4.1.1.1\tA|B" + ROW, r"path=.*u\.tsv.*contrib=.*c\.tsv.*one '\|'"),
    "repeated-sample": (UNSTRAT.replace("S3", "S2"), CONTRIB, r"path=.*u\.tsv.*repeats column names \['S2'\]"),
    "description": (
        "function\tdescription\tS1\nEC:1.1.1.1\tAlcohol dehydrogenase\t1.0\n",
        CONTRIB,
        r"u\.tsv.*'description'",
    ),
    "short-row": (UNSTRAT + "EC:1.1.1.1\t1\t1\n", CONTRIB, r"u\.tsv.*missing or NaN"),
}


@pytest.mark.parametrize("case", MALFORMED)
def test_malformed_inputs_raise_naming_the_file(tmp_path, case):
    unstrat, contrib, message = MALFORMED[case]
    with pytest.raises(ValueError, match=message):
        bt.io.read_picrust2(write(tmp_path, unstrat, "u.tsv"), contrib=write(tmp_path, contrib, "c.tsv"))


@given(
    st.lists(
        st.tuples(st.integers(0, 2), st.integers(0, 2), st.integers(0, 3), st.floats(0.01, 1e6, allow_nan=False)),
        min_size=1,
        max_size=12,
        unique_by=lambda row: row[:3],
    )
)
def test_each_contribution_lands_in_its_cell(rows):
    # (sample, function, taxon, value) rows in any order; each must appear once, at (sample, function|taxon).
    contrib = CONTRIB_COLUMNS + "".join(
        f"S{s}\tEC:{f}.1.1.1\tASV{t}\t1\t1\t1\t{value!r}\t1\t1\n" for s, f, t, value in rows
    )
    with tempfile.TemporaryDirectory() as tmp:
        unstrat = Path(tmp) / "u.tsv"
        unstrat.write_text("function\tS0\tS1\tS2\n" + "".join(f"EC:{f}.1.1.1\t0\t0\t0\n" for f in range(3)))
        (Path(tmp) / "c.tsv").write_text(contrib)
        by_taxon = bt.io.read_picrust2(unstrat, contrib=Path(tmp) / "c.tsv")["function_by_taxon"]
    X = pd.DataFrame(by_taxon.X.toarray(), index=by_taxon.obs_names, columns=by_taxon.var_names)
    got = [X.loc[f"S{s}", f"{f}.1.1.1|ASV{t}"] for s, f, t, _ in rows]
    # Looser than exact: pandas' default C float parser keeps about 15 significant digits (as in test_humann.py).
    np.testing.assert_allclose(got, [value for *_, value in rows], rtol=1e-12)
    assert by_taxon.X.nnz == len(rows)
```

- [x] **Step 2: Run, expect failure.** `uv run --group test pytest tests/io/test_picrust2.py -q`
  -> `17 failed, 7 errors`, all `AttributeError: module 'biotapy.io' has no attribute 'read_picrust2'`.
  The 7 errors are the tests that use the `mdata` fixture.
- [x] **Step 3: Implement.** Create `src/biotapy/io/_picrust2.py`.
  `_contributions` is the single-use helper R4.4 allows: inlined,
  `read_picrust2` has 37 statements and complexity 9. `_first_line` serves
  both files here, and the trait reader in 2.4b.

```python
"""PICRUSt2 predictions: metagenome and pathway tables and their contributions."""

import gzip
from pathlib import Path

import pandas as pd
import scipy.sparse as sp
from mudata import MuData

from biotapy._core import make_function_mudata

from ._table import _numbers, _read_table

# PICRUSt2 writes EC numbers as "EC:1.1.1.1"; ENZYME, HUMAnN and bt.datasets.enzyme write "1.1.1.1".
_EC_PREFIX = "EC:"
# The contribution columns read_picrust2 uses; the long table's first column is "sample".
_CONTRIB_COLUMNS = ("function", "taxon", "taxon_function_abun")


def read_picrust2(path: str | Path, *, contrib: str | Path | None = None) -> MuData:
    r"""Read a PICRUSt2 prediction, and optionally its contributions, into community and per-taxon modalities.

    Parameters
    ----------
    path
        An unstratified PICRUSt2 table: ``pred_metagenome_unstrat.tsv.gz``
        (EC, KO or another trait) or ``path_abun_unstrat.tsv.gz``.
    contrib
        The matching long-format contributions, written with
        ``--stratified``: ``pred_metagenome_contrib.tsv.gz`` or
        ``path_abun_contrib.tsv.gz``. Its ``taxon_function_abun`` column
        becomes the per-taxon modality.

    Returns
    -------
    MuData
        Two modalities over ``path``'s samples, as ``bt.io.read_humann``
        returns them:

        - ``"function"``: one feature per function, ``var`` columns ``name``
          and ``special``;
        - ``"function_by_taxon"``: one feature per function and taxon
          (``1.1.1.1|ASV1``), ``var`` columns ``function``, ``name``,
          ``taxon``, ``genus``, ``species`` and ``special``; no features
          without ``contrib``.

        ``taxon`` is the ASV id (or ``RARE``, PICRUSt2's group of rare
        ASVs), so ``genus`` and ``species`` are NaN. ``EC:`` is removed from
        EC numbers. ``uns['biotapy']['x_kind']`` is ``"abundance"``.

    Raises
    ------
    ValueError
        A file is empty or malformed (a value that is not a number, a missing
        value or id, a data row with more cells than the header, repeated
        column names); ``path`` has a ``description`` column; ``contrib``
        lacks a ``sample``, ``function``, ``taxon`` or
        ``taxon_function_abun`` column, repeats a sample, function and taxon,
        or names a sample or function that ``path`` lacks. Messages name the
        file's argument.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/reading_data`

    PICRUSt2's predictions are marker-normalised read counts weighted by
    predicted gene copy numbers, not counts, so ``bt.pp.rarefy`` refuses
    them. A sample's gene-family contributions sum to its unstratified value;
    a pathway's need not, so neither modality is derived from the other.

    ``EC:`` is removed so the ids match ENZYME's and HUMAnN's
    (``bt.datasets.enzyme``); remove it from a PICRUSt2 mapping file too
    before regrouping with ``bt.fn.func_glom``.

    The unstratified table is read into one dense functions x samples
    ``float64`` array, and ``contrib`` into one pandas table of all its rows,
    before the result is stored sparse.

    References
    ----------
    Douglas GM et al. (2020) PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38:685-688.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> folder = Path(tempfile.mkdtemp())
    >>> _ = (folder / "unstrat.tsv").write_text("function\tS1\nEC:1.1.1.1\t6.0\n")
    >>> columns = "sample\tfunction\ttaxon\ttaxon_abun\ttaxon_rel_abun\tgenome_function_count"
    >>> columns += "\ttaxon_function_abun\ttaxon_rel_function_abun\tnorm_taxon_function_contrib\n"
    >>> _ = (folder / "contrib.tsv").write_text(columns + "S1\tEC:1.1.1.1\tASV1\t3.0\t100.0\t2\t6.0\t200.0\t1.0\n")
    >>> mdata = bt.io.read_picrust2(folder / "unstrat.tsv", contrib=folder / "contrib.tsv")
    >>> mdata["function"].var_names.tolist(), mdata["function_by_taxon"].var_names.tolist()
    (['1.1.1.1'], ['1.1.1.1|ASV1'])
    """
    path = Path(path)
    argument = f"path={str(path)!r}"
    table = _read_table(path, _first_line(path), skiprows=0, argument=argument)
    if "description" in table.columns:
        msg = f"{argument} has a 'description' column (add_descriptions.py output); pass the table without it"
        raise ValueError(msg)
    X = sp.csr_matrix(_numbers(table, argument=argument).T)
    row_ids = table.index.str.removeprefix(_EC_PREFIX)
    if contrib is not None:
        by_taxon, keys = _contributions(Path(contrib), samples=table.columns, functions=row_ids)
        X, row_ids = sp.hstack([X, by_taxon], format="csr"), row_ids.append(keys)
    try:
        return make_function_mudata(
            X, obs=pd.DataFrame(index=table.columns), row_ids=row_ids, x_kind="abundance", source="io.read_picrust2"
        )
    except ValueError as error:  # repeated functions, or a function or taxon id holding "|"
        msg = f"{argument}, contrib={None if contrib is None else str(contrib)!r}: {error}"
        raise ValueError(msg) from error


def _contributions(path: Path, *, samples: pd.Index, functions: pd.Index) -> tuple[sp.csr_matrix, pd.Index]:
    """The long contribution table as a samples x (function, taxon) matrix, and its ``function|taxon`` ids."""
    argument = f"contrib={str(path)!r}"
    table = _read_table(path, _first_line(path), skiprows=0, argument=argument, text=3)
    missing = [column for column in _CONTRIB_COLUMNS if column not in table.columns]
    if table.index.name != "sample" or missing:
        msg = f"{argument} needs PICRUSt2's long-format columns sample, {', '.join(_CONTRIB_COLUMNS)}; found {[table.index.name, *table.columns]}"
        raise ValueError(msg)
    values = _numbers(table[["taxon_function_abun"]], argument=argument).ravel()
    ids = table[["function", "taxon"]]
    if ids.isna().any().any():
        msg = f"{argument} has a row with no function or taxon"
        raise ValueError(msg)
    function = ids["function"].str.removeprefix(_EC_PREFIX)
    rows = samples.get_indexer(table.index)
    for name, unknown in (("samples", table.index[rows < 0]), ("functions", function[~function.isin(functions)])):
        if len(unknown):
            msg = f"{argument} names {name} that path lacks: {sorted(set(unknown))[:3]}"
            raise ValueError(msg)
    codes, keys = pd.factorize(function + "|" + ids["taxon"])
    if pd.Series(rows * len(keys) + codes).duplicated().any():
        msg = f"{argument} repeats a sample, function and taxon"
        raise ValueError(msg)
    matrix = sp.csr_matrix((values, (rows, codes)), shape=(len(samples), len(keys)))
    return matrix, pd.Index(keys)


def _first_line(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return handle.readline()
```

  `src/biotapy/io/__init__.py` becomes (ruff format wraps `__all__`):

```python
from ._biom import read_biom, write_biom
from ._dada2 import read_dada2
from ._humann import read_humann
from ._metaphlan import read_metaphlan
from ._phyloseq import read_phyloseq
from ._picrust2 import read_picrust2
from ._qiime2 import read_qiime2

__all__ = [
    "read_biom",
    "read_dada2",
    "read_humann",
    "read_metaphlan",
    "read_phyloseq",
    "read_picrust2",
    "read_qiime2",
    "write_biom",
]
```

- [x] **Step 4: Run, expect pass.** The same command gives `24 passed`. The
  doctest `uv run --group test pytest src/biotapy/io/_picrust2.py -q` gives
  `1 passed`.
- [x] **Step 5: Docs and contracts.**
  - `docs/api.md`: add `io.read_picrust2` after `io.read_phyloseq`.
  - `docs/guide/reading_data.md`: append the PICRUSt2 section.
  - `docs/guide/function.md`: add "PICRUSt2 tables" before "Aggregating
    along a hierarchy".
  - `data-model-slots.md`: convention 2 and the Function tables section.
  - `r-golden-parity.md`: statements 6 and 8.

  The exact text:

````diff
diff --git a/.knowledge/contracts/data-model-slots.md b/.knowledge/contracts/data-model-slots.md
index 3785146..a15e4bf 100644
--- a/.knowledge/contracts/data-model-slots.md
+++ b/.knowledge/contracts/data-model-slots.md
@@ -61,6 +61,8 @@ Extends the spec's data-model table with exact keys.[^spec]
    `relative`) and labels a header without a unit `abundance`, never `counts`.
    `io.read_metaphlan` divides MetaPhlAn's percentages by 100 and sets `relative`, after
    checking that every sample's leaf clades sum to 1 within `1e-3` (`_core.RELATIVE_TOLERANCE`).
+   `io.read_picrust2` sets `abundance`: PICRUSt2's values are read counts divided by predicted
+   marker copies and multiplied by gene copies.
    `fn.renorm` rescales `X` (and may drop the special rows), setting `x_kind` to `relative` or `cpm`;
    its `relative` stratified rows do not sum to 1: they are shares of the community total.
    The other readers infer it because their formats record no unit (BIOM, QIIME 2 `RelativeFrequency`, a DADA2
@@ -87,9 +89,9 @@ Extends the spec's data-model table with exact keys.[^spec]
    ids) never collide silently.

 ## Function tables
-`io.read_humann` (and, from Phase 2 slice 2B, `io.read_picrust2`) returns a
-`MuData` built by `_core.make_function_mudata` with two modalities over the
-same samples, each an `AnnData` with its own copy of `obs`:
+`io.read_humann` and `io.read_picrust2` return a `MuData` built by
+`_core.make_function_mudata` with two modalities over the same samples, each
+an `AnnData` with its own copy of `obs`:

 | Modality | Features | `var` columns |
 |---|---|---|
@@ -102,6 +104,13 @@ columns use the pandas `str` dtype, as rank columns do. Both modalities
 always exist; either may have 0 features. The community values are not
 the sum of their strata for pathways, which is why there are two.

+In a PICRUSt2 table the stratified rows come from the long contribution
+table (`taxon_function_abun`), and the modality has 0 features when no
+contribution table is read. `taxon` is the ASV id as written, or `RARE`
+(PICRUSt2's group of rare ASVs: an ordinary stratum, not `special`), so
+`genus` and `species` are NaN. `EC:` is removed from EC numbers in both
+modalities, so ids match ENZYME's and HUMAnN's.
+
 ## Taxonomic profiles (MetaPhlAn)
 `io.read_metaphlan` keeps one feature per leaf clade: a row that no other
 row descends from through any ancestor (MetaPhlAn can omit an intermediate
diff --git a/.knowledge/contracts/r-golden-parity.md b/.knowledge/contracts/r-golden-parity.md
index 05c6d53..ba1e738 100644
--- a/.knowledge/contracts/r-golden-parity.md
+++ b/.knowledge/contracts/r-golden-parity.md
@@ -56,7 +56,9 @@ sources:
    stay synthetic, except small files copied under a permissive licence
    with a `NOTICE.txt` beside them (`tests/data/humann`: HUMAnN's MIT test
    data; `tests/data/metaphlan`: a MetaPhlAn 4.0.6 profile from HUMAnN's MIT
-   test data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0).
+   test data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0). PICRUSt2
+   (GPL-3) fixtures are always synthetic, written from its documented
+   column headers.
 7. `pl` functions have an R equivalent but no golden test. They draw numbers
    that `tl` stores, and `tl`'s golden tests check those numbers (controller
    ruling 2026-09-27; rules.md R11.2).
@@ -69,7 +71,8 @@ sources:
    `mia::importMetaPhlAn`) is checked against MetaPhlAn's own output: its
    leaves, grouped by `pp.tax_glom` to each rank, equal the clade rows the
    profile prints (`tests/io/test_metaphlan.py`, atol `1e-6` because
-   MetaPhlAn rounds each percentage to 5 decimals).
+   MetaPhlAn rounds each percentage to 5 decimals). `io.read_picrust2` has no
+   R equivalent; invariants on synthetic files check it.

 # Why
 R and NumPy random generators differ, so stochastic outputs can never match
diff --git a/docs/api.md b/docs/api.md
index 4d5008a..ec1787b 100644
--- a/docs/api.md
+++ b/docs/api.md
@@ -16,6 +16,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
     io.read_humann
     io.read_metaphlan
     io.read_phyloseq
+    io.read_picrust2
     io.read_qiime2
     io.write_biom
 ```
diff --git a/docs/guide/function.md b/docs/guide/function.md
index 8df302f..9703b28 100644
--- a/docs/guide/function.md
+++ b/docs/guide/function.md
@@ -31,6 +31,26 @@ flagged in `var["special"]`, so a sample's total keeps what HUMAnN could not
 assign. Read one table per call: a gene family table and a pathway table both
 hold `UNMAPPED`.

+## PICRUSt2 tables
+
+`bt.io.read_picrust2` gives PICRUSt2 predictions the same two modalities. Its
+`"function_by_taxon"` features are functions per ASV (`2.7.1.1|ASV1`), read
+from the long contribution table. For gene families a sample's
+contributions sum to its community value; for pathways they need not, as in
+HUMAnN. EC numbers lose PICRUSt2's `EC:` prefix, so they match
+`bt.datasets.enzyme()`:
+
+```python
+import biotapy as bt
+
+mdata = bt.io.read_picrust2("pred_metagenome_unstrat.tsv.gz", contrib="pred_metagenome_contrib.tsv.gz")
+by_class = bt.fn.func_glom(mdata["function"], "class", hierarchy=bt.datasets.enzyme())
+```
+
+PICRUSt2's own mapping files write `EC:1.1.1.1`; remove the prefix from the
+edge table's `child` column before regrouping with them, or `func_glom`
+raises because nothing maps.
+
 ## Aggregating along a hierarchy

 `bt.fn.func_glom` sums functions into their parents at one level of a
diff --git a/docs/guide/reading_data.md b/docs/guide/reading_data.md
index fdeab8b..b91f760 100644
--- a/docs/guide/reading_data.md
+++ b/docs/guide/reading_data.md
@@ -318,3 +318,32 @@ tdata = bt.io.read_metaphlan("merged_abundance_table.tsv")
   marker tables), raises a `ValueError` instead of double counting.
 - **Not read.** NCBI taxids, `additional_species`, coverage and estimated
   read counts.
+
+## PICRUSt2
+
+`bt.io.read_picrust2` reads a PICRUSt2 prediction into the same two-modality
+`MuData` as `bt.io.read_humann`: the unstratified table becomes
+`"function"`, and the long-format contributions written with `--stratified`
+become `"function_by_taxon"`, one feature per function and ASV.
+
+```python
+import biotapy as bt
+
+mdata = bt.io.read_picrust2(
+    "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
+    contrib="EC_metagenome_out/pred_metagenome_contrib.tsv.gz",
+)
+pathways = bt.io.read_picrust2("pathways_out/path_abun_unstrat.tsv.gz")
+```
+
+- **Ids.** `EC:1.1.1.1` becomes `1.1.1.1`, the form ENZYME, HUMAnN and
+  `bt.datasets.enzyme` use. KO and pathway ids are kept.
+- **Taxa.** A stratified feature's `taxon` is the ASV id as written (`0042`
+  stays `0042`), or `RARE` for the rare ASVs PICRUSt2 groups together.
+- **Units.** `x_kind` is `"abundance"`: the values are read counts divided by
+  predicted marker copies and multiplied by gene copies, not counts.
+- **Not read.** The deprecated wide stratified table (`--wide_table`) and
+  tables with an `add_descriptions.py` `description` column raise a
+  `ValueError`; pass the table without descriptions. PICRUSt2's pathway
+  coverage tables have the same layout as its abundance tables and cannot be
+  told apart, so pass only abundance tables.
````

  Build the docs as in 2.2 Step 6 -> `build succeeded.`
- [x] **Step 6: Bookkeeping.**
  - Tick this task's boxes and the 2.4 part of its checklist line.
  - Bump `generated` and `commit` in `data-model-slots.md` and
    `r-golden-parity.md`.
  - Add to `.knowledge/log.md`:
    `- **Update**: [data-model-slots](contracts/data-model-slots.md) convention 2 (PICRUSt2 is abundance) and Function tables (read_picrust2: stratified rows from the contribution table, taxon = ASV id or RARE, EC: removed); [r-golden-parity](contracts/r-golden-parity.md) statements 6 and 8 (synthetic PICRUSt2 fixtures, no R equivalent); [phase-2-function](roadmap/phase-2-function.md) task 2.4 done.`
- [x] **Step 7: Gate and commit.**

```bash
uvx prek run --all-files
git add src/biotapy/io/_picrust2.py src/biotapy/io/__init__.py tests/io/test_picrust2.py \
  docs/api.md docs/guide/reading_data.md docs/guide/function.md .knowledge/contracts/data-model-slots.md \
  .knowledge/contracts/r-golden-parity.md .knowledge/roadmap/phase-2-function.md .knowledge/log.md
git commit -m "feat(io): read PICRUSt2 predictions and contributions into community and per-taxon modalities"
uv run --group test pytest -q -W error::UserWarning   # 810 passed, 22 deselected
```

### Task 2.4b: `io.read_picrust2_traits`

**Why a separate function (user-approved new public function, Phase 2
decision 4):** the per-ASV copy-number table has taxa, not samples, as rows.
It cannot be a modality of the samples' MuData, because mudata's `axis=0`
would add the ASVs to the sample index.

**Files:** modify `src/biotapy/io/_picrust2.py`, `src/biotapy/io/__init__.py`,
`tests/io/test_picrust2.py`, `docs/api.md`, `docs/guide/reading_data.md`,
`.knowledge/contracts/data-model-slots.md`, `.knowledge/contracts/r-golden-parity.md`.

**Interfaces:**
- Consumes `_read_table`, `_numbers` and 2.4's `_first_line`, `_EC_PREFIX`.
- Produces `bt.io.read_picrust2_traits(path: str | Path) -> pd.DataFrame`:
  - ASVs x functions, `float64`;
  - index = ASV ids as text, unnamed; columns = function ids with `EC:`
    removed;
  - no `metadata_NSTI`.

  2.8 calls `fn.functional_redundancy(adata, *, traits=read_picrust2_traits(...))`.

**Will not touch:** `read_picrust2`.

- [ ] **Step 1: Failing tests.** Append to `tests/io/test_picrust2.py`:

```python
TRAITS = (
    "sequence\tEC:1.1.1.1\tEC:2.7.1.1\tEC:3.2.1.1\tmetadata_NSTI\n"
    "ASV1\t1\t2\t0\t0.03\n"
    "0042\t0\t1\t1\t0.12\n"
    "ASV9\t0\t0\t0\t1.5\n"
)


def test_traits_are_asvs_by_functions(tmp_path):
    traits = bt.io.read_picrust2_traits(write(tmp_path, TRAITS, "EC_predicted.tsv"))
    assert traits.index.tolist() == ["ASV1", "0042", "ASV9"]
    assert traits.columns.tolist() == ["1.1.1.1", "2.7.1.1", "3.2.1.1"]
    assert traits.dtypes.eq(np.float64).all() and traits.loc["ASV9"].eq(0).all()


def test_traits_read_gzip_and_keep_ko_ids(tmp_path):
    path = tmp_path / "KO_predicted.tsv.gz"
    path.write_bytes(gzip.compress(b"sequence\tK00001\tK00002\nASV1\t1\t0\n"))
    assert bt.io.read_picrust2_traits(path).columns.tolist() == ["K00001", "K00002"]


@pytest.mark.parametrize(
    ("text", "message"),
    [
        (TRAITS + "ASV1\t1\t1\t1\t0.1\n", r"repeats ASV ids: \['ASV1'\]"),
        (TRAITS + "ASV2\t1\tx\t1\t0.1\n", "not a number"),
        (TRAITS + "ASV2\t1\t1\n", "missing or NaN"),
        ("", "not a valid tab-separated table"),
    ],
    ids=["repeated-asv", "non-number", "short-row", "empty-file"],
)
def test_malformed_traits_raise_naming_the_path(tmp_path, text, message):
    with pytest.raises(ValueError, match=rf"traits\.tsv.*{message}"):
        bt.io.read_picrust2_traits(write(tmp_path, text, "traits.tsv"))
```

- [ ] **Step 2: Run, expect failure.** `uv run --group test pytest tests/io/test_picrust2.py -q -k traits`
  -> `6 failed`: `AttributeError: module 'biotapy.io' has no attribute 'read_picrust2_traits'`.
- [ ] **Step 3: Implement.**
  - In `src/biotapy/io/_picrust2.py`, change the module docstring to
    `"""PICRUSt2 predictions: metagenome and pathway tables, their contributions, and per-ASV trait tables."""`.
  - Insert this function before `def _first_line`:

```python
def read_picrust2_traits(path: str | Path) -> pd.DataFrame:
    r"""Read PICRUSt2's predicted gene copy numbers per ASV.

    Parameters
    ----------
    path
        A per-sequence trait table, such as ``EC_predicted.tsv.gz`` or
        ``KO_predicted.tsv.gz``: one row per ASV, one column per function.

    Returns
    -------
    pandas.DataFrame
        ASVs x functions, ``float64`` copy numbers, indexed by the ASV ids as
        written (``0042`` stays ``0042``). ``EC:`` is removed from EC numbers,
        as in ``bt.io.read_picrust2``; a ``metadata_NSTI`` column is dropped.

    Raises
    ------
    ValueError
        The file is empty or malformed (a value that is not a number, a
        missing value or id, a data row with more cells than the header,
        repeated column names), or an ASV id repeats. Messages name ``path``.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/reading_data`

    The table describes genomes, not samples, so it is a DataFrame rather
    than a modality of the samples' MuData. It is read into one dense
    ASVs x functions ``float64`` array.

    References
    ----------
    Douglas GM et al. (2020) PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38:685-688.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import biotapy as bt
    >>> path = Path(tempfile.mkdtemp()) / "EC_predicted.tsv"
    >>> _ = path.write_text("sequence\tEC:1.1.1.1\tEC:2.7.1.1\tmetadata_NSTI\nASV1\t1\t2\t0.03\n")
    >>> bt.io.read_picrust2_traits(path)
          1.1.1.1  2.7.1.1
    ASV1      1.0      2.0
    """
    path = Path(path)
    argument = f"path={str(path)!r}"
    table = _read_table(path, _first_line(path), skiprows=0, argument=argument)
    table = table.drop(columns=[column for column in table.columns if column == "metadata_NSTI"])
    repeated = table.index[table.index.duplicated()].unique().tolist()
    if repeated:
        msg = f"{argument} repeats ASV ids: {repeated[:3]}"
        raise ValueError(msg)
    values = _numbers(table, argument=argument)
    return pd.DataFrame(values, index=table.index.rename(None), columns=table.columns.str.removeprefix(_EC_PREFIX))
```

  - `src/biotapy/io/__init__.py` becomes:

```python
from ._biom import read_biom, write_biom
from ._dada2 import read_dada2
from ._humann import read_humann
from ._metaphlan import read_metaphlan
from ._phyloseq import read_phyloseq
from ._picrust2 import read_picrust2, read_picrust2_traits
from ._qiime2 import read_qiime2

__all__ = [
    "read_biom",
    "read_dada2",
    "read_humann",
    "read_metaphlan",
    "read_phyloseq",
    "read_picrust2",
    "read_picrust2_traits",
    "read_qiime2",
    "write_biom",
]
```

- [ ] **Step 4: Run, expect pass.** `uv run --group test pytest tests/io/test_picrust2.py -q`
  gives `30 passed`. The doctests
  `uv run --group test pytest src/biotapy/io/_picrust2.py -q` give `2 passed`.
- [ ] **Step 5: Docs and contracts.**
  - `docs/api.md`: add `io.read_picrust2_traits` after `io.read_picrust2`.
  - `docs/guide/reading_data.md`: append the paragraph below to the PICRUSt2
    section.
  - Make the two contract edits shown.

````diff
diff --git a/.knowledge/contracts/data-model-slots.md b/.knowledge/contracts/data-model-slots.md
index a15e4bf..d3aac34 100644
--- a/.knowledge/contracts/data-model-slots.md
+++ b/.knowledge/contracts/data-model-slots.md
@@ -108,8 +108,9 @@ In a PICRUSt2 table the stratified rows come from the long contribution
 table (`taxon_function_abun`), and the modality has 0 features when no
 contribution table is read. `taxon` is the ASV id as written, or `RARE`
 (PICRUSt2's group of rare ASVs: an ordinary stratum, not `special`), so
-`genus` and `species` are NaN. `EC:` is removed from EC numbers in both
-modalities, so ids match ENZYME's and HUMAnN's.
+`genus` and `species` are NaN. `EC:` is removed from EC numbers, in both
+modalities and in `io.read_picrust2_traits`' columns, so ids match ENZYME's
+and HUMAnN's.

 ## Taxonomic profiles (MetaPhlAn)
 `io.read_metaphlan` keeps one feature per leaf clade: a row that no other
diff --git a/.knowledge/contracts/r-golden-parity.md b/.knowledge/contracts/r-golden-parity.md
index ba1e738..2422543 100644
--- a/.knowledge/contracts/r-golden-parity.md
+++ b/.knowledge/contracts/r-golden-parity.md
@@ -71,8 +71,9 @@ sources:
    `mia::importMetaPhlAn`) is checked against MetaPhlAn's own output: its
    leaves, grouped by `pp.tax_glom` to each rank, equal the clade rows the
    profile prints (`tests/io/test_metaphlan.py`, atol `1e-6` because
-   MetaPhlAn rounds each percentage to 5 decimals). `io.read_picrust2` has no
-   R equivalent; invariants on synthetic files check it.
+   MetaPhlAn rounds each percentage to 5 decimals). `io.read_picrust2` and
+   `io.read_picrust2_traits` have no R equivalent; invariants on synthetic
+   files check them.

 # Why
 R and NumPy random generators differ, so stochastic outputs can never match
diff --git a/docs/api.md b/docs/api.md
index ec1787b..41a2c83 100644
--- a/docs/api.md
+++ b/docs/api.md
@@ -17,6 +17,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
     io.read_metaphlan
     io.read_phyloseq
     io.read_picrust2
+    io.read_picrust2_traits
     io.read_qiime2
     io.write_biom
 ```
diff --git a/docs/guide/reading_data.md b/docs/guide/reading_data.md
index b91f760..2f41546 100644
--- a/docs/guide/reading_data.md
+++ b/docs/guide/reading_data.md
@@ -347,3 +347,7 @@ pathways = bt.io.read_picrust2("pathways_out/path_abun_unstrat.tsv.gz")
   `ValueError`; pass the table without descriptions. PICRUSt2's pathway
   coverage tables have the same layout as its abundance tables and cannot be
   told apart, so pass only abundance tables.
+
+`bt.io.read_picrust2_traits` reads the predicted gene copy numbers per ASV
+(`EC_predicted.tsv.gz`, `KO_predicted.tsv.gz`) into an ASVs x functions
+`DataFrame`. It describes genomes, not samples, so it is not a modality.
````

  Build the docs -> `build succeeded.`
- [ ] **Step 6: Bookkeeping.**
  - Tick this task's boxes and the rest of its checklist line.
  - Bump `generated` and `commit` in both contracts.
  - Add to `.knowledge/log.md`:
    `- **Update**: [data-model-slots](contracts/data-model-slots.md) Function tables and [r-golden-parity](contracts/r-golden-parity.md) statement 8 name io.read_picrust2_traits (EC: removed from its columns; no R equivalent); [phase-2-function](roadmap/phase-2-function.md) task 2.4b done.`
- [ ] **Step 7: Gate and commit.**

```bash
uvx prek run --all-files
git add src/biotapy/io/_picrust2.py src/biotapy/io/__init__.py tests/io/test_picrust2.py docs/api.md \
  docs/guide/reading_data.md .knowledge/contracts/data-model-slots.md .knowledge/contracts/r-golden-parity.md \
  .knowledge/roadmap/phase-2-function.md .knowledge/log.md
git commit -m "feat(io): read PICRUSt2 per-ASV copy numbers"
uv run --group test pytest -q -W error::UserWarning   # 818 passed, 22 deselected
```

### Checkpoint B - review slice 2B
- [ ] **Review the whole slice** with superpowers:requesting-code-review,
  against:
  - data-model-slots, function-shape, module-boundaries and r-golden-parity;
  - no-bundled-kegg and function-tables-as-mudata;
  - the slice 2B review focus above.

  Then a fix pass, one commit per finding, each with a test. Reviewers may
  run the readers on real files they hold. Never commit those files, and
  never copy a PICRUSt2 output into `tests/` (GPL).
- [ ] **Run the gates** on the committed tree:
  - `uvx prek run --all-files`;
  - `uv run --group test pytest -q -W error::UserWarning`;
  - `uv run --group test pytest -m golden tests/fn -q`;
  - `BIOTAPY_DATA_DIR=<scratchpad>/pooch uv run --group doc sphinx-build -W -b html docs docs/_build/html`.

  Confirm `~/.cache/biotapy` does not exist.
- [ ] **Knowledge** (codebase-map templates; R12.2-R12.4).
  - **Update `.knowledge/modules/io.md`.**
    - Responsibility: add `read_metaphlan` (TreeData), `read_picrust2`
      (MuData) and `read_picrust2_traits` (DataFrame).
    - Entry points:
      - `_metaphlan.py:read_metaphlan`, `_metaphlan.py:_header`;
      - `_picrust2.py:read_picrust2`, `_picrust2.py:_contributions`,
        `_picrust2.py:read_picrust2_traits`, `_picrust2.py:_first_line`;
      - `_table.py:_read_table`, `_table.py:_numbers`, replacing the
        `_humann.py:_read_table` entry.
    - Invariants:
      - the leaf rule (any ancestor) and the 100% check;
      - `x_kind` set by the reader, never inferred, for MetaPhlAn
        (`relative`) and PICRUSt2 (`abundance`);
      - `EC:` removed;
      - `RARE` not special;
      - every message names its argument;
      - `_read_table` rejects repeated column names.
    - Gotchas:
      - MetaPhlAn may omit an intermediate rank's row (4.0.6 fixture);
      - pandas renames repeated header cells and drops extra cells under
        `usecols`;
      - PICRUSt2 coverage tables look like abundance tables;
      - the contribution table is held whole in memory;
      - `pp.tax_glom` drops `UNCLASSIFIED` by default;
      - the h5mu writer turns `str` columns into categoricals in place.
    - Update its `description`, and the copy in `modules/index.md`.
  - **Update `.knowledge/modules/core.md`:** `RELATIVE_TOLERANCE` is
    exported, with `io/_metaphlan.py` as its second consumer.
  - **Update `.knowledge/decisions/function-tables-as-mudata.md`:** add
    `src/biotapy/io/_picrust2.py` to `paths`. In Context, say that PICRUSt2's
    two files fill the same two modalities.
  - Add log lines; bump `generated` and `commit` on each changed concept;
    tick this box.
- [ ] **Push** `phase-2b`, open the PR, and merge-commit on green (approved
  2026-10-03). CI must be green, including docs and the network job.
- [ ] **Ask the user to review slice 2B** before slice 2C.

### Slice 2B decisions for the user
Each changes a contract, the public surface or an earlier slice, or is a
judgement call. The recommended answer comes first.

1. **MetaPhlAn leaf rule: structural, any ancestor, with a hard 100% check.**
   - No `rank=` parameter: `pp.tax_glom` gives higher ranks and matches
     MetaPhlAn's own rows.
   - A sample whose leaves do not sum to 100% (1e-3) raises instead of
     warning, as the outline proposed.
   - Contract addition: data-model-slots "Taxonomic profiles (MetaPhlAn)".
2. **No `var["sgb"]` column.** The SGB id is the feature id (`SGB1871`, prefix
   removed), and `t__` gets no rank column. The roadmap had proposed a
   contract addition for `sgb`.
3. **Strip PICRUSt2's `EC:` prefix at read**, in both modalities and in the
   trait table, so ids match ENZYME, HUMAnN and `bt.datasets.enzyme`. Users
   regrouping with PICRUSt2's own `EC:`-prefixed mapping files strip their
   hierarchy's `child` column instead, and `func_glom`'s nothing-maps error
   tells them to. Not recorded in provenance (R2.3).
4. **`read_picrust2(path, *, contrib=None)` signature.**
   - One call reads an unstratified table and its contribution table.
   - The outline's `taxa=` (`seqtab_norm`) is dropped (R2.3): 2.8 takes any
     samples x ASVs AnnData.
   - Rejected with a clear `ValueError`:
     - the legacy wide stratified table;
     - `description` columns from `add_descriptions.py` (the alternative is
       reading them into `var["name"]`);
     - coverage contributions.
   - Pathway coverage unstratified tables cannot be detected; the docs say
     so.
5. **`RARE` is an ordinary taxon, not `special`.** The outline proposed a
   special flag, which `fn.renorm(special=False)` would drop as if it were
   unmapped abundance.
6. **New task 2.2a changes a slice 2A reader.** `read_humann`'s strict
   parsing moves to `io/_table.py` (R4.3, three readers), and `read_humann`
   gains two behaviours:
   - a repeated sample column now raises; today pandas silently reads `S1`
     and `S1.1`;
   - a row with no id raises "no id" instead of `_core`'s "var ids must not
     be missing" (both name `path`).
7. **`_core` exports `RELATIVE_TOLERANCE`.** MetaPhlAn's 100% check uses the
   data model's own definition of `relative`.
8. **MetaPhlAn sample names.** A single profile is named after its file, and
   `_profile` is removed from every sample name, as
   `merge_metaphlan_tables.py` does. The alternative is the profile's
   `#SampleID` line, which defaults to the useless `Metaphlan_Analysis`. With
   this rule, HMP2's MetaPhlAn and HUMAnN tables share all 1,638 names.
9. **Parity without mia.** `read_metaphlan`'s check is MetaPhlAn's own output
   (`tax_glom` equals the printed rows), and r-golden-parity statement 8
   grows to say so. The PICRUSt2 readers have no R equivalent and rest on
   invariants and synthetic tables.
10. **Forward note for 2.8.** With `taxa=` dropped, 2.8 must decide whether
    Tian's `p` needs 16S-copy-corrected abundances (`seqtab_norm`). If so, it
    adds a reader for that file (a new public function to approve then).
11. **Reported, not fixed (R1.4): a slice 2A flaky test.** Under `--cov`,
    `tests/fn/test_renorm.py::test_strata_are_raw_values_over_the_community_total`
    took 205.9 ms once against Hypothesis's 200 ms deadline (`DeadlineExceeded`).
    It passes without coverage. If CI's coverage job runs Hypothesis tests,
    this can fail at random. Suggested separate fix: `@settings(deadline=None)`
    on that test.

### Slice 2B self-review
Run against the brief, this concept's outline and the writing-plans
checklist.

1. **Spec coverage.**

   | Brief or outline item | Where it lands |
   |---|---|
   | MetaPhlAn leaf rule, `rank=`, additive totals (q1) | design; 2.2 tests "leaf clades", "glom ... own rows", "do not sum to 100" |
   | `var["sgb"]`, `t__` ranks (q2) | design; decision 2; 2.2 `test_var_holds_the_ranks_of_each_leaf` |
   | Single vs merged, read stats, empty taxids, `#` lines (q3) | 2.2 tests for 4.2 profile, merged 4.x, merged 3 with taxids, HMP2, read stats with short rows, the 4.0.6 fixture (empty taxids) |
   | `EC:` prefix (q4) | design; 2.4 `test_ec_ids_match_enzyme_style_hierarchies`; 2.4b column test |
   | PICRUSt2 inputs and signature (q5) | design; 2.4 tests for unstrat only, contrib, pathways, legacy wide and description rejected |
   | `x_kind`, additivity (q6) | 2.4 `test_x_kind_is_abundance_so_rarefy_refuses_it`, `test_gene_family_contributions_sum_to_the_community_values`, `test_pathway_contributions_are_kept_as_written` |
   | Module placement, shared parsing (q7) | design; Task 2.2a |
   | Binding: no mia, PICRUSt2 never copied, the MetaPhlAn fixture and notice, the same MuData, traits DataFrame, / 100 with `UNCLASSIFIED`, `path` in every error, malformed input rejected | global constraints; 2.2 Step 1; 2.4 and 2.4b; the malformed-input tests in 2.2a, 2.2, 2.4 and 2.4b |
   | R11.2 edge cases (all-zero sample and feature, missing rank, single sample) and properties | 2.2 `test_all_zero_sample_and_feature_are_kept`, the MetaPhlAn 3 genus leaf (species NaN), the single profiles, `test_leaves_keep_every_percentage`; 2.4 `test_all_zero_sample_and_feature_are_kept`, `test_single_sample`, `test_each_contribution_lands_in_its_cell`; 2.4b the all-zero `ASV9` |

   Purity tests do not apply: readers take paths. No golden test exists,
   because mia was declined; statement 8's MetaPhlAn parity test takes its
   place.

   Gaps found and fixed while writing:
   - the skipped `o__` row in the real fixture (the leaf rule became "any
     ancestor");
   - pandas' silent rename of a repeated header cell (2.2a);
   - numeric-looking ASV ids (`text=3`);
   - the reading guide's false "every reader infers `x_kind`" paragraph.
2. **Placeholder scan.** Every step shows the full file, the exact diff or
   the exact text, rendered from the scratch commits that passed the gates.
   The only values left to fill at run time are the implementer's model id,
   the UTC time, the log heading's date, `HEAD`'s short hash in bookkeeping
   steps, and `<scratchpad>` (the session's scratchpad path).
3. **Type consistency.**
   - The shared helpers are used with the same signature in 2.2a, 2.2, 2.4
     and 2.4b: `_read_table(path, header, *, skiprows, argument, text=1)` and
     `_numbers(frame, *, argument)`.
   - `read_picrust2` returns `make_function_mudata`'s layout, so `var`
     columns match `read_humann`, as 2.7 and 2.9 expect.
   - `read_picrust2_traits`' index uses the same text ids as
     `var["taxon"]`.
4. **Review focus.** Each of the five items names tests that exist in the
   code above; checked by searching this section for each name.
5. **Known residual risks.**
   - PICRUSt2 formats come from its headers and research B, never from a
     run. The first real user file is the check (Phase 2 risk list).
   - `path_abun_contrib`'s column count (8 per research B, 9 in PICRUSt2's
     per-sequence test data) does not matter, because the reader keys on
     column names.
   - The additivity formula was not re-read here (GPL source).
   - `mia::importMetaPhlAn`'s exact behaviour was never seen (its source
     could not be fetched), so nothing claims equality with it beyond the
     Coming-from-R mapping.

---

## Slice 2C - Analysis (outline)

### Task 2.7: `fn.contributions`
**Interface:** `bt.fn.contributions(mdata: MuData, function: str, *, top: int | None = None) -> pd.DataFrame`.
Returns samples x taxa (samples as rows, like `tl.alpha`): the
`"function_by_taxon"` rows whose `var["function"] == function`, columns =
`taxon`. With `top`, the `top` taxa by total are kept and the rest summed
into `"other"`. R equivalent: none.

**Design questions:**
1. **Raw values or shares of the community value?** For pathways the strata
   do not sum to the community value. Proposed: raw values, the strata as
   HUMAnN wrote them. A share of the community value is a display choice; it
   goes into `pl.contributions` only if the tutorial needs it (R2.3).
2. **Accepting `func_glom` output** (stratified groups): yes, since it has
   the same `var` columns. Test it.
3. **A function absent from the table:** `KeyError` naming up to three ids
   that contain the string.

**Tests:** row sums equal the sum of that function's strata (Hypothesis, with
and without `top`); an all-zero sample; one sample; `unclassified` kept;
purity.

### Task 2.8: `fn.functional_redundancy`
**Interface:** `bt.fn.functional_redundancy(adata: AnnData, *, traits: pd.DataFrame) -> pd.DataFrame`.
- `adata`: samples x taxa abundances (PICRUSt2's `"taxa"` modality, or any
  taxa table whose `var_names` match `traits.index`).
- `traits`: taxa x genes copy numbers (2.4b).
- Returns samples x `["taxonomic_diversity", "functional_diversity",
  "redundancy", "normalized_redundancy"]`. Tian 2020 terms:
  - TD = Gini-Simpson, 1 - Σp²;
  - FD = Rao's Q, pᵀDp;
  - FR = TD - FD;
  - nFR = FR / TD.
- R equivalent: none. References: Tian L et al. (2020) Deciphering
  functional redundancy in the human microbiome. *Nat Commun* 11:6217.

**Method facts (checked):**
- Tian's taxon dissimilarity is weighted Jaccard,
  d = 1 - Σmin(Gᵢ,Gⱼ)/Σmax(Gᵢ,Gⱼ). For non-negative vectors it equals
  2·BC/(1+BC), where BC is SciPy's `pdist(G, "braycurtis")`. Verified
  numerically 2026-10-03 with scipy 1.18.1. So the kernel is one SciPy
  call (R2.1), not a custom loop.
- The reference implementation (MATLAB, liangtian85/FR) has no licence. Do
  not read it into the code; implement from the paper's equations (CC BY
  4.0 article).

**Design questions:**
1. **Taxa without traits, or with all-zero gene vectors.** BC of zero
   against zero is 0 in SciPy; against non-zero it is 1. Proposed: drop taxa
   absent from `traits` with a `warn_user` count. Keep all-zero taxa
   (d = 1 to others).
2. **Memory.** D is N x N dense (10,000 ASVs = 800 MB). Proposed: restrict to
   taxa present in at least one sample, compute D once, and state the cost in
   `Notes` (R6.2). Measure before any optimisation (R10.1).
3. **Normalising p:** per sample over the taxa kept; an all-zero sample
   gives NaN for all four (never raises), as `tl.alpha` does.
4. **A method note in docs:** a "Functional redundancy" section in
   `docs/guide/function.md` with the equations (R8.3: method math lives in
   `docs/`), instead of the roadmap's non-existent `docs/methods/`.

**Tests (no golden: no R equivalent):**
- hand-computed toy cases from the paper's equations;
- Hypothesis properties: 0 <= FD <= TD; 0 <= FR <= TD; identical genomes
  give FR = TD (D = 0); disjoint gene sets give FR = 0 (D = 1 off the
  diagonal);
- permutation of taxa leaves results unchanged;
- purity.

### Task 2.9: `pl.contributions`
**Interface:** `bt.pl.contributions(mdata: MuData, function: str, *, top: int = 8, ax: Axes | None = None, plot_kwargs: dict[str, object] | None = None) -> Axes`.
A stacked bar per sample of `fn.contributions(mdata, function, top=top)`.
It computes nothing new: `pl` may import `fn`, a lower layer. R equivalent:
none.

**Design questions:**
1. Reuse `pl/_abundance.py`'s stacked-bar drawing (same subpackage) rather
   than a second implementation (R4.3). Check its private helper's
   signature first.
2. Optional grouping/faceting by an `obs` column, as `pl.bar` does: only if
   the tutorial needs it (R2.3).

**Tests:** returns `Axes`; one bar per sample, one segment per taxon;
`other` present with `top`; purity; no golden (r-golden-parity statement 7).

### Checkpoint C
Review; `modules/fn.md` and `modules/pl.md` updates; ask the user to review
before 2D.

---

## Slice 2D - Cohort, docs and release (outline)

### Task 2.10: `datasets.hmp2()` and the tutorial
**Interface (proposed):** `bt.datasets.hmp2() -> MuData` with modalities
`"function"`, `"function_by_taxon"` (pathway abundance, CPM) and `"taxa"`
(MetaPhlAn 3 profiles via `read_metaphlan`), and global `obs` from the HMP2
metadata (diagnosis CD / UC / nonIBD, participant, week). It is a fixed,
documented subset of samples.

**Sources (pinned SHA-256, research C):**
- `.../HMP2/MGX/2018-05-04/pathabundances_3.tsv.gz` (13.3 MB):
  `dd983871b0e155255844b91ec10d50fb09230d2f4e915464ab680fa3a9c9ddb3`;
- `.../taxonomic_profiles_3.tsv.gz` (0.54 MB):
  `d790ff15e46d61ca0cadc55d9f918de4e3415d7f97c992ac37610aaee02117ed`;
- `metadata/hmp2_metadata_2018-08-20.csv` (9.1 MB): hash to compute;
- base `https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/`.

Never committed (R6.6; no licence stated).

**Design questions:**
1. **Subset rule** (deterministic): e.g. the first stool sample of each of
   40 participants, balanced by diagnosis. The loader reads all 1,638 columns
   (3 s, measured) and subsets, or reads `usecols`.
2. **Sample ids.** 2A's reader already strips `_pathabundance_cpm` and sets
   `x_kind="cpm"`. The metadata's `External ID` joins on `CSM5FZ3N_P`; check
   the MetaPhlAn table's column names.
3. **A `func_glom` demo.** HMP2 has no merged gene-family table. Options:
   - fetch about 10 per-sample `*_humann3.tar.bz2` (1.4 MB each) for their
     `level4ec` tables, then group them to ENZYME classes;
   - demonstrate on pathways with a user-style map. MetaCyc pathway classes
     are restricted, so the map would be synthetic.

   Proposed: the per-sample EC route (open data end to end).
4. **CI.** The docs job sets `BIOTAPY_DATA_DIR` and caches it; the notebook
   runs under `nb_execution_mode = "cache"` with
   `nb_execution_raise_on_error = True`. Budget: under 2 minutes.
5. **Tutorial content** (`docs/tutorials/function.md`):
   - read, renorm and group;
   - contributions plot for a butyrate pathway;
   - functional redundancy needs PICRUSt2 traits, which HMP2 lacks. Show it
     on the toy, or leave it out of the tutorial. Decide here.
   - Licence caveat and citation: Lloyd-Price J et al. (2019) Multi-omics of
     the gut microbial ecosystem in inflammatory bowel diseases. *Nature*
     569:655-662.

### Task 2.11: Knowledge
`modules/fn.md` and `modules/io.md` final pass (MetaPhlAn, PICRUSt2); the
data-model contract's MetaPhlAn/PICRUSt2 rules; `roadmap/index.md` when the
phase closes; log lines. (The MuData decision and the first `fn` concept
were written at Checkpoint A.)

### Task 2.13: Coming-from-R check
The table is generated (R8.4). New rows come from docstrings:
`mia::importHUMAnN`, `mia::importMetaPhlAn`. Check the rendered page. Decide
whether `docs/_data/r_idioms.toml`'s "not in 0.1" entries become "not in
0.2"; `tests/test_coming_from_r.py` asserts the literal "not in 0.1".

### Task 2.14: Benchmarks
asv `benchmarks/benchmarks/fn.py`, measuring only (R10.1):
- `func_glom` on a synthetic 1,600 x 22,000 stratified table with a
  many-to-many map;
- `read_humann` on the same written to TSV;
- `functional_redundancy` at N = 2,000 taxa (O(N²) memory).

No optimisation without a profile.

### Task 2.15: Release 0.2
Follow [cut-a-release](/playbooks/cut-a-release.md):
- version 0.2.0 and the CHANGELOG;
- the exit gate ticked;
- the tag pushed only with the user's approval (R13.3).

---

# Decisions for the user
Each changes a contract or rule, adds a public function or dependency, or is
a judgement call. Recommended answer first.

1. **ENZYME's entry point.** It is `bt.datasets.enzyme()`, not
   `fn.load_hierarchy("enzyme")`, so every download stays in `datasets` and
   `fn` has no network code (design note 5). The ledger's ruling named
   `load_hierarchy`.
2. **ENZYME's hash is unpinned** (`None`). ENZYME keeps no old releases. The
   release read goes into `attrs["source"]` and `func_glom`'s provenance.
   This is biotapy's first unpinned download.
3. **New task 2.12 `fn.renorm(mdata, units, *, special=True)`.** The exit
   gate requires renorm parity, but the roadmap has no task for it. It
   replaces `X` and sets `x_kind`; `units` has no default.
4. **New public functions in 2A:**
   - `datasets.toy_humann()`: in-memory doctest data for `fn`;
   - `datasets.enzyme()`.

   In 2B: `io.read_picrust2_traits()` (taxa x genes cannot be a modality of
   a samples MuData).
5. **`func_glom` semantics.** Match HUMAnN: `UNGROUPED` and protected
   specials, no warning on partial mapping, and a `ValueError` when nothing
   maps. This replaces roadmap 2.1's "dropped and counted in a warning".
6. **Signature changes from the roadmap:**
   - `read_humann(path)`: one table per call, no `pathcoverage`;
   - `load_hierarchy(path, level, *, layout)`: local files only, no URL;
   - `func_glom(adata, level, *, hierarchy, agg="sum")`;
   - `functional_redundancy(adata, *, traits)` with no `method`.
7. **Contract and rule changes:**
   - `data-model-slots`: the Function tables section, the HUMAnN `x_kind`
     header rule, `replace_features`, the `func_glom` aggregation semantics,
     and the `renorm` propagation;
   - `r-golden-parity`: HUMAnN goldens via `uv run`; third-party fixtures
     with a NOTICE;
   - `function-shape` and **rules.md R8.2**: examples may use `toy_humann()`;
   - `no-bundled-kegg`: local paths only, and ENZYME as the one download.
8. **`x_kind` rules:**
   - HUMAnN from the header, never `counts`;
   - MetaPhlAn divided by 100 into `relative`, keeping `UNCLASSIFIED` (2B);
   - PICRUSt2 `abundance`.
9. **Slice order.** `read_humann` moves into 2A, so the golden cross-check
   runs through the reader. MetaPhlAn and PICRUSt2 readers form 2B.
10. **`_core` placement** of `sum_pairs` and `replace_features`. Each has one
    consumer today, which is in tension with R4.3. They follow roadmap 2.1
    and the `sum_by` precedent.
11. Slice 2B decisions: see the slice's own list.
12. **Frontmatter `description`** of phase-2-function.md, as proposed in
    design note 7.

# Self-review
Run against the brief, the roadmap concept and the writing-plans checklist.

1. **Spec coverage.**

   | Roadmap item | Where it lands |
   |---|---|
   | 2.1 | Tasks 2.1 / 2.1b |
   | 2.2 | 2B outline |
   | 2.3 | Task 2.3 |
   | 2.4 | 2B outline, with 2.4b |
   | 2.5 | Tasks 2.5a / 2.5 |
   | 2.6 | Task 2.6 |
   | 2.7-2.9 | 2C |
   | 2.10, 2.11 | 2D |
   | "renormalisation cross-checked" (exit gate) | Task 2.12 plus golden tests |
   | "Licensing of each open source checked and written in the docstring" | `enzyme()` docstring and NOTICE files |

   Brief design questions 1-7 are answered in design notes 1-7. Gaps found
   and fixed while writing:
   - renorm had no task (added as 2.12);
   - doctests for `fn` had no function data (added `toy_humann`);
   - the HMP2 table has no `#` header (the reader now handles it, with a
     test).
2. **Placeholder scan.** Every slice 2A step carries the full file or the
   exact lines to add, rendered from the prototype that passed every gate. No
   step says "add tests" without code. Outline slices are allowed open
   questions by the brief; each question gives a proposed answer.
3. **Type consistency.** The same names and types are used throughout:
   - `sum_pairs(..., *, n_groups)` (keyword-only, as ruff `PLR0917`
     demanded);
   - `replace_features(adata, X, var)`;
   - `make_function_mudata(X, *, obs, row_ids, x_kind, source)`;
   - the edge table `child, parent, level, parent_name`;
   - `FUNCTION_KEY` / `BY_TAXON_KEY` = `"function"` / `"function_by_taxon"`.

   2C's interfaces use the same modality names and `var["function"]`.
4. **Review focus.** Each of the five Phase 2 items and the two slice 2A
   items names a test that exists in the rendered code. Checked by
   searching for each test name in this file.
5. **Known residual risks.**
   - The intermediate states between tasks were not each re-run in the
     prototype. The final state was, and each task's tests depend only on
     earlier tasks. Two tests that need `func_glom` move to Task 2.6 for
     that reason (`test_round_trips_through_func_glom`,
     `test_feeds_func_glom`).
   - `uv add` in 2.1b was not run here: auto-mode denied installing into a
     venv. The prototype used an ephemeral `--with` overlay instead, so the
     exact `uv.lock` diff is unverified.

[^spec]: Python Microbiome Toolkit development report, sections Positioning and Roadmap
[^humann]: HUMAnN
