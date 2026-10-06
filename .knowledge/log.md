# Knowledge bundle log

## 2026-10-06 (Phase 3, slice 3C)
- **Update**: [data-model-slots](contracts/data-model-slots.md) DA result rows: `se` is NaN for every `da.aldex2` feature, `effect` units per method (MaAsLin 3 per SD and minus its median, ALDEx2 median `diff.btw`), and only `effect`, `pvalue`, `qvalue` are NaN together for an untested feature (Checkpoint C fix, I1).
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.11 fix round 1: the `r-bridge` job has a 30-minute timeout and logs the R and package versions; Checkpoint C must read `CFFI_MODE.API` and `19 passed` in its log.
- **Update**: [r-bridge-before-ports](decisions/r-bridge-before-ports.md) consequence: the `r-bridge` CI job runs the `r` tests. [phase-3-stats](roadmap/phase-3-stats.md) task 3.11 done.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.7 fix round 1: `da/_design._column` returns unordered categoricals (an ordered group shrank MaAsLin 3's effect by 1/sqrt(2)), and the `da.maaslin3` docs state the median and the `reference` swap exactly.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md): `da.maaslin3` is compared with `maaslin3::maaslin3` (abundance model, median subtracted) elementwise, seed-matched. [phase-3-stats](roadmap/phase-3-stats.md) task 3.7 done.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `maaslin3` for the `da.maaslin3` golden file.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.6 fix round 2: `r_function` catches R warning conditions in R (`withCallingHandlers`) instead of rpy2's console output, so `message()` stays a message.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md): a seeded Monte Carlo bridge matches R exactly when biotapy passes R its seed integer (the generic DA row and the Why section said otherwise), and the golden's integer is recorded as NumPy-stream dependent. [phase-3-stats](roadmap/phase-3-stats.md) task 3.6 fix round 1: duplicate names raise, R warnings and errors surface through `da/_r.py`, `seed` is checked before R loads.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md): `da.aldex2` is compared with `ALDEx2::aldex` elementwise (the golden's seed is the one biotapy derives, so the Monte Carlo draws match); R bridge golden tests carry the marker `r` only. [phase-3-stats](roadmap/phase-3-stats.md) task 3.6 done.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `ALDEx2` for the `da.aldex2` golden file.

## 2026-10-06 (Phase 3, slice 3C plan)
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) expands slice 3C (tasks 3.6 `da.aldex2` with the extra `r` and `da/_r.py`, 3.7 `da.maaslin3`, 3.11 CI job `r-bridge`, Checkpoint C) into full TDD steps, prototyped and gated per commit; ticks Checkpoint B (PR #22 merged, slice 3B approved).

## 2026-10-05 (Phase 3, slice 3C)
- **Update**: [da](modules/da.md) gotcha: LinDA's single densify is allowed by rules.md R6.2, whose wording the user approved after Checkpoint B (a native method whose algorithm needs the full table may densify once).

## 2026-10-05 (Phase 3, Checkpoint B)
- **Create**: [da](modules/da.md): the module concept for `linda`, `ancombc2` and `consensus` (one result schema, NaN means not tested, strict `q < alpha`, no filtering, no formulas) with the gotchas the review found (ANCOM-BC2 not antisymmetric in `reference`, no residual degrees of freedom, numeric-group units, prefixed scikit-bio errors, LinDA's densify against R6.2 left unresolved); description copied into [modules/index.md](modules/index.md), and [index](index.md) lists `da`.
- **Update**: [pl](modules/pl.md) gains `consensus` (Responsibility, Entry points, and the invariant that it reads `da.consensus`'s columns, imports nothing from `da` and keeps its legend one row above the axes); description copied into [modules/index.md](modules/index.md).
- **Update**: [core](modules/core.md) lists `da._design.dense_counts` among `require_counts`' callers; [tl](modules/tl.md) says non-negative whole-number counts; [module-boundaries](contracts/module-boundaries.md) says `pl` draws what `da` gives.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) ticks Checkpoint B's review (0 Critical / 2 Important / 5 Minor; fix pass `9f45c28..fd9e94a` plus `927e5ae`; re-review 7/7 addressed), exit-gate check (4 golden/network passed; 1174 passed in the full run) and knowledge boxes; task 3.14 now updates `modules/da.md` for the R bridges; push and merge stay open.
- **Verification**: re-checked against `6ade269..927e5ae` and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [phase-1-core](roadmap/phase-1-core.md), [phase-2-function](roadmap/phase-2-function.md), [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md), [io](modules/io.md), [data-model-slots](contracts/data-model-slots.md), [engine-parity](contracts/engine-parity.md), [function-shape](contracts/function-shape.md), [r-golden-parity](contracts/r-golden-parity.md), [tree-access](contracts/tree-access.md), [add-a-function](playbooks/add-a-function.md), [cut-a-release](playbooks/cut-a-release.md), [regenerate-golden-files](playbooks/regenerate-golden-files.md).

## 2026-10-05 (Phase 3, slice 3B)
- **Update**: Checkpoint B fixes. [da-consensus-agreement](decisions/da-consensus-agreement.md): consensus counts depend on `reference` through ANCOM-BC2 (104 against 112 on GlobalPatterns `host`). [data-model-slots](contracts/data-model-slots.md): `require_counts` wants non-negative whole numbers, as `infer_x_kind` does. [phase-3-stats](roadmap/phase-3-stats.md) slice 3B design prose states what the code does after the fix rounds (schema checks, `model`'s errors, ANCOM-BC2's untested rule).
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.9 code and tests follow the second review fix: `pl.consensus`'s legend is one row above the axes (a default save clips it beside them), and its Notes say the effect ranking mixes units for a numeric `group`.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.9 code and tests follow the review fix: `pl.consensus` validates `table` and `top` (TypeError/ValueError naming the argument), draws its legend outside the axes and documents row order and `top` above 30.
- **Update**: [function-shape](contracts/function-shape.md) (and rules.md R3.2): `pl.consensus` takes the table `da.consensus` returns. [phase-3-stats](roadmap/phase-3-stats.md) task 3.9 done.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.8 code and tests follow the second fix: `da.consensus` accepts NumPy scalars for `alpha` and `min_methods`.
- **Update**: [da-consensus-agreement](decisions/da-consensus-agreement.md): a call whose `effect` is exactly 0 has no direction (counts in `n_significant`, no consensus, no conflict); [phase-3-stats](roadmap/phase-3-stats.md) task 3.8 code and tests follow the review fixes (`validate_result` counts NaN method/contrast and requires NaN `effect` where `pvalue` is NaN; `consensus` type-checks `alpha` and `min_methods` and refuses an empty `results`).
- **Create**: [da-consensus-agreement](decisions/da-consensus-agreement.md): what "methods agree" means in `da.consensus` (strict `q < alpha`, one BH, untested is not "not significant", same sign, conflict), the options rejected, and that the user picks the methods; listed in [decisions/index.md](decisions/index.md).
- **Update**: [function-shape](contracts/function-shape.md) (and rules.md R3.2): `da.consensus` takes `da` result tables instead of an AnnData. [phase-3-stats](roadmap/phase-3-stats.md) task 3.8 done.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md): `da.ancombc2`'s tolerances are per model (`host` loose, `host + log_depth` 1e-6); [phase-3-stats](roadmap/phase-3-stats.md) task 3.4 code, tests and counts follow the fix round (a feature with no residual degrees of freedom is untested, scikit-bio errors are re-raised with the function name).
- **Update**: [r-golden-parity](contracts/r-golden-parity.md): `da.ancombc2` is compared with `ANCOMBC::ancombc2` at the measured tolerances. [phase-3-stats](roadmap/phase-3-stats.md) task 3.4 done.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `ANCOMBC` for the `da.ancombc2` golden file, with CVXR 1.0-15 and `libgsl27`. [regenerate-golden-files](playbooks/regenerate-golden-files.md) Common mistakes: why CVXR is pinned.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) Task 3.5 code and tests synced with the review fixes: `da._design` also raises for fewer than two features, a constant `group` or covariate column, non-list `covariates` and a non-string `reference`, and densifies `X` without a second copy.
- **Update**: [data-model-slots](contracts/data-model-slots.md) gains "DA results": the schema `da` methods return, written by no slot; `paths` gains `src/biotapy/da/**`. [r-golden-parity](contracts/r-golden-parity.md): `da.linda` is compared elementwise with `MicrobiomeStat::linda`. [phase-3-stats](roadmap/phase-3-stats.md) tasks 3.3 and 3.5 done.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs CRAN `MicrobiomeStat` for the `da.linda` golden file.
- **Update**: [data-model-slots](contracts/data-model-slots.md) and [core](modules/core.md) state that `x_kind` `"counts"` means non-negative whole numbers (`infer_x_kind` and `require_counts` reject a table holding a negative value; Task 3.B0); [phase-3-stats](roadmap/phase-3-stats.md) ticks 3.B0.

## 2026-10-05 (Phase 3, slice 3B plan)
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) expands slice 3B (tasks 3.5 `da.linda` with the result schema, 3.4 `da.ancombc2`, 3.8 `da.consensus`, 3.9 `pl.consensus`, Checkpoint B) into full TDD steps, prototyped and gated per commit; adds Task 3.B0 (`fix(core)`: counts are non-negative whole numbers) approved with the plan; ticks Checkpoint A (PR #21 merged, slice 3A approved).

## 2026-10-05 (Phase 3, Checkpoint A)
- **Update**: [pp](modules/pp.md) documents `clr`, `philr`, `pseudocounted` and `_binary_tree` (pseudocount rule, dense `layers["clr"]`, `obsm["X_philr"]` layout and sign, one-child and wide nodes) and the gotchas: `tree_basis`/`TreeNode.prune` conventions, `toy()`'s three-child root, plain slicing keeping stale derived slots, balance names, peak memory, `pseudocounted` moving to `_core` if 3B reuses it; description copied into [modules/index.md](modules/index.md).
- **Update**: [core](modules/core.md) `get_skbio_tree`'s `split_root` and its one caller, `pp.philr`.
- **Update**: [pl](modules/pl.md) invariant: `bar` and `heatmap` raise `ValueError` on a table with a negative value (`_common.py:table`); NaN is not checked.
- **Update**: [regenerate-golden-files](playbooks/regenerate-golden-files.md) Common mistakes: `philr` needs `libuv1` in the R image.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) ticks Checkpoint A's review (0 Critical / 1 Important / 5 Minor; fix pass `c920c7b..0de6cca` plus `6ade269`; re-review 14/14), exit-gate check (8 golden/network passed; 1039 passed in the full run) and knowledge boxes; push and the user's review stay open.
- **Verification**: re-checked against `c015d3d..6ade269` and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [phase-1-core](roadmap/phase-1-core.md), [phase-2-function](roadmap/phase-2-function.md), [data-model-slots](contracts/data-model-slots.md), [engine-parity](contracts/engine-parity.md), [function-shape](contracts/function-shape.md), [module-boundaries](contracts/module-boundaries.md), [r-golden-parity](contracts/r-golden-parity.md), [tree-access](contracts/tree-access.md), [add-a-function](playbooks/add-a-function.md), [cut-a-release](playbooks/cut-a-release.md).

## 2026-10-05 (Phase 3, slice 3A)
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) Task 3.1's `pp.clr` Notes block states the measured peak memory (3.1x to 4.8x one dense array), matching the code and the transforms guide.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.1 and 3.2 code, test and docs blocks match the Checkpoint A fix pass (`pseudocounted(columns=)` reorders while sparse, philr's wide-node message counts every wide node, narrowed mypy exclude, 57 passed). [data-model-slots](contracts/data-model-slots.md): `X_philr` has one column per internal node with two children, named as `vart["phylo"]` names it. [r-golden-parity](contracts/r-golden-parity.md): PhILR atol reason is about 1e-16 on both sides.
- **Update**: [data-model-slots](contracts/data-model-slots.md): `obsm["X_philr"]` from `pp.philr` and its embedding-adding propagation row. [r-golden-parity](contracts/r-golden-parity.md): PhILR balances are matched by partition. [tree-access](contracts/tree-access.md): `get_skbio_tree(split_root=False)` and the child-order guarantee PhILR's signs rely on. [phase-3-stats](roadmap/phase-3-stats.md) task 3.2 done.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.1 code blocks follow the review fixes: `pseudocount` type check, all-zero row and peak-memory notes.
- **Update**: [data-model-slots](contracts/data-model-slots.md): `layers["clr"]` from `pp.clr` is a dense float64 array; `pp.clr` joins the layer-adding row. [r-golden-parity](contracts/r-golden-parity.md): CLR is compared elementwise. [phase-3-stats](roadmap/phase-3-stats.md) task 3.1 done.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `philr` (and `libuv1`) for the `pp.philr` golden file.
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.0 done: CLR and PhILR golden files from vegan 2.7.3 and philr 1.36.0.

## 2026-10-05 (Phase 3 plan)
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) carries the user-approved Phase 3 plan: resolved design notes, global constraints, dependencies, review focus, slices 3A-3D, slice 3A in full TDD steps (3.0 goldens, 3.1 pp.clr, 3.2 pp.philr, Checkpoint A) and later slices as outlines, decisions and self-review; new description, paths and sources, copied into the [roadmap index](roadmap/index.md).

## 2026-10-05 (release 0.2.0)
- **Update**: Phase 2 closed after biotapy 0.2.0 reached PyPI (tag v0.2.0, release workflow run 37316580850). [phase-2-function](roadmap/phase-2-function.md) is `phase_state: done` with every Task 2.15 step and exit-gate item ticked; [phase-3-stats](roadmap/phase-3-stats.md) is `phase_state: in-progress`; the [roadmap index](roadmap/index.md) lists Phase 3 as active.
- **Update**: [cut-a-release](playbooks/cut-a-release.md) step 2 says how to write `## [Unreleased]` from the git log when it is empty, and new step 2c moves the "not in X.Y" labels.
- **Verification**: re-checked against the 0.2.0 version bump and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [module-boundaries](contracts/module-boundaries.md), [tree-access](contracts/tree-access.md).
- **Update**: [phase-2-function](roadmap/phase-2-function.md) ticks Checkpoint D (push, merge as PR #18 `d1b89b6`, Read the Docs, the user's go), the exit gate (docs job run 37305490515, golden and Phase 1 gates on PR #18) and Task 2.15 Steps 1-8.

## 2026-10-05 (Phase 2, Checkpoint D knowledge)
- **Update**: [datasets](modules/datasets.md) documents `hmp2` (entry point, pinned IBDMDB files fetched before any is parsed, first-metagenome ordering by `week_num`, `visit_num`, `External ID`, pushed metadata, the 290 MB dense read and 1 GB peak, 23 MB cold download, mudata's `convert_dtypes` on the global `obs`, MuData not surviving `pickle`, synthetic offline fixtures); description copied into [modules/index.md](modules/index.md).
- **Update**: [data-model-slots](contracts/data-model-slots.md) Function tables: `datasets.hmp2` adds a `"taxa"` modality and pushes its metadata into every modality; `fn.renorm` keeps other modalities.
- **Update**: [function-tables-as-mudata](decisions/function-tables-as-mudata.md) forward note: `datasets.hmp2` already holds `function`, `function_by_taxon` and `taxa` side by side; `paths` gains `datasets/_hmp2.py`.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) ticks the Checkpoint A-C boxes left open after slices 2A-2C merged (PRs #15-#17), each with its record, Checkpoint D's review (0 Critical / 0 Important / 7 Minor; fix pass `9f7b271..ef2fe8b`; re-review 6/6) and Task 2.11.
- **Verification**: re-checked against `407cc19..ef2fe8b` and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [phase-1-core](roadmap/phase-1-core.md), [engine-parity](contracts/engine-parity.md), [function-shape](contracts/function-shape.md), [module-boundaries](contracts/module-boundaries.md), [r-golden-parity](contracts/r-golden-parity.md), [add-a-function](playbooks/add-a-function.md).

## 2026-10-05 (Phase 2, slice 2D)
- **Update**: [phase-2-function](roadmap/phase-2-function.md) slice 2D blocks match the Checkpoint D fix pass: `hmp2` orders a participant's metagenomes by `week_num`, `visit_num`, `External ID` (C3007 keeps `CSM5MCVB_P`; PWY-5676 in 113 samples, 25 without strata), fetches all three files before parsing, notes its 1 GB peak; the synthetic fixture pins values per sample and uses neutral ids; the tutorial says what the per-species means are.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) Task 2.10b's tutorial block matches the review fix (empty bars for the 24 samples without per-species PWY-5676 rows; links to the guide's sections).
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.10 done: `datasets.hmp2` returns each HMP2 participant's first metagenome as pathway and species modalities with the metadata.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.10b done: the HMP2 function tutorial runs on every docs build; notebook cells may take 300 s.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.13 done: phyloseq calls without an equivalent read "not in 0.2"; the mia importer rows are pinned by a test.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.14 done: asv baselines for `func_glom`, `read_humann` and `functional_redundancy` in docs/performance.md.

## 2026-10-05 (slice 2D plan)
* **Update**: [phase-2-function](roadmap/phase-2-function.md) carries the user-approved slice 2D plan in full TDD steps (2.10 datasets.hmp2, 2.10b tutorial, 2.13 Coming-from-R, 2.14 benchmarks, Checkpoint D with 2.11 knowledge, 2.15 release 0.2.0) and its decisions; header note, checklist, slices table and decision 11 updated.

## 2026-10-05 (Phase 2, Checkpoint C)
- **Update**: [fn](modules/fn.md) documents `contributions` and `functional_redundancy` (entry points, invariants incl. the boundary checks, the SciPy zero-vector NaN, 16S correction with PICRUSt2's NSTI cutoff applied by the recipe, memory and time cost, `divide_rows` rule); description copied into [modules/index.md](modules/index.md).
- **Update**: [pl](modules/pl.md) documents `contributions`, `_stack`, the `fn` import and `_colors` never giving a group the NA/"other" grey (`pl.bar` colours changed for 8+ groups plus NA); description copied into [modules/index.md](modules/index.md).
- **Update**: [core](modules/core.md) lists `divide_rows` with its three consumers and the divide-don't-multiply gotcha.
- **Update**: [pp](modules/pp.md): `relative` divides each float64-summed value through `divide_rows` instead of multiplying by a reciprocal.
- **Update**: [phase-2-function](roadmap/phase-2-function.md): slice 2C design, decisions and risks carry the corrected NSTI facts and the guide's current recipe; Task 2.8's `_redundancy.py`, test block and guide text equal the files at `f5236e8`; decisions 8 and 11 note the colour fix and `divide_rows`; Checkpoint C's review, gates and knowledge boxes ticked with the review record, and its draft diff replaced by a summary.
- **Update**: [r-golden-parity](contracts/r-golden-parity.md), [add-a-function](playbooks/add-a-function.md), [module-boundaries](contracts/module-boundaries.md), [data-model-slots](contracts/data-model-slots.md): `pl` now draws `tl`, `pp` and `fn` results, and `pl.contributions` has no R equivalent.
- **Verification**: re-checked against `c49aeb4..f5236e8` and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [phase-1-core](roadmap/phase-1-core.md) (its old `pp.relative` block is historical, not rewritten), [phase-3-stats](roadmap/phase-3-stats.md), [function-tables-as-mudata](decisions/function-tables-as-mudata.md), [engine-parity](contracts/engine-parity.md), [function-shape](contracts/function-shape.md).

## 2026-10-05 (Checkpoint C fix pass)
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.8 `_redundancy.py` code block and task 2.9 `_abundance.py` diff now match the fix pass (float64 sample totals, `adata` type check, PICRUSt2 wording, `toy()` example, Notes on the samples x taxa products, explicit taxon count in `pl.contributions`).

## 2026-10-03 (Phase 2, slice 2C)
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.7 done: `fn.contributions` returns one function's strata as a samples x taxa table.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.8 done: `fn.functional_redundancy` computes Tian et al. 2020's TD, FD, FR and nFR per sample.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.8 code and tests now match the fix round (in-place condensed distances with a measured peak note, sample-total overflow and duplicate taxon id errors, `divide_rows`, nullable-NA test).
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.9 done: `pl.contributions` draws `fn.contributions`' table as stacked bars, sharing `bar`'s drawing (`_abundance.py:_stack`).

## 2026-10-03 (slice 2C plan)
* **Update**: [phase-2-function](roadmap/phase-2-function.md) carries the user-approved slice 2C plan in full TDD steps (2.7 fn.contributions, 2.8 fn.functional_redundancy, 2.9 pl.contributions, Checkpoint C) and its decisions, with the pp.relative subnormal fix approved as a separate commit; checklist signatures, header note, decision 11 and the 2.14 outline updated.

## 2026-10-03 (Phase 2, Checkpoint B)
- **Update**: [io](modules/io.md) documents `read_metaphlan`, `read_picrust2` and `read_picrust2_traits` and the shared `_table.py` (header rule, strict checks, `utf-8-sig`, first-cell check, contribution-sample check, seven rank columns), replacing the removed `_humann.py:_read_table`; its description is copied into [modules/index.md](modules/index.md).
- **Update**: [core](modules/core.md) lists `RELATIVE_TOLERANCE` as exported (second consumer `io.read_metaphlan`), `make_function_mudata` also used by `io.read_picrust2`, and `normalize_ranks` used by `io.read_phyloseq` and `io.read_metaphlan`.
- **Update**: [function-tables-as-mudata](decisions/function-tables-as-mudata.md) adds `_picrust2.py` to `paths` and says PICRUSt2's two files fill the same two modalities.
- **Update**: [data-model-slots](contracts/data-model-slots.md): three readers set `x_kind` themselves; an all-zero sample is exempt from the MetaPhlAn 100% check; rank columns are always the seven.
- **Update**: [phase-2-function](roadmap/phase-2-function.md): slice 2B design, decision 6 and the Checkpoint B box describe the final code (no `_first_line`; `_header` in `_table.py`; `read_humann`'s own header rule), the stale `repeats column names` quotes gain their colon, `_numbers` gains `nonnegative=`, Task 2.8 gets the `RARE` forward note, and the Checkpoint B review, gates and Knowledge boxes are ticked.
- **Verification**: re-checked against HEAD with no content change needed, and bumped `commit`: [engine-parity](contracts/engine-parity.md), [function-shape](contracts/function-shape.md), [module-boundaries](contracts/module-boundaries.md), [r-golden-parity](contracts/r-golden-parity.md), [add-a-function](playbooks/add-a-function.md), [phase-0-foundation](roadmap/phase-0-foundation.md), [phase-1-core](roadmap/phase-1-core.md) and [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md).

## 2026-10-03 (Phase 2, slice 2B)
- **Update**: [data-model-slots](contracts/data-model-slots.md) Function tables and [r-golden-parity](contracts/r-golden-parity.md) statement 8 name io.read_picrust2_traits (EC: removed from its columns; no R equivalent); [phase-2-function](roadmap/phase-2-function.md) task 2.4b done.
- **Update**: [data-model-slots](contracts/data-model-slots.md) convention 2 (PICRUSt2 is abundance) and Function tables (read_picrust2: stratified rows from the contribution table, taxon = ASV id or RARE, EC: removed); [r-golden-parity](contracts/r-golden-parity.md) statements 6 and 8 (synthetic PICRUSt2 fixtures, no R equivalent); [phase-2-function](roadmap/phase-2-function.md) task 2.4 done.
- **Update**: [data-model-slots](contracts/data-model-slots.md) convention 2 (MetaPhlAn percentages / 100, checked to sum to 1) and new section "Taxonomic profiles (MetaPhlAn)"; [r-golden-parity](contracts/r-golden-parity.md) statement 6 lists tests/data/metaphlan and statement 8 gives read_metaphlan's parity (tax_glom equals MetaPhlAn's own rows); [phase-2-function](roadmap/phase-2-function.md) task 2.2 done.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) task 2.2a done: io/_table.py holds the strict TSV reading shared by the HUMAnN, MetaPhlAn and PICRUSt2 readers; read_humann now rejects repeated column names.

## 2026-10-03 (slice 2B plan)
* **Update**: [phase-2-function](roadmap/phase-2-function.md) carries the user-approved slice 2B plan in full TDD steps (2.2a shared strict table reading, 2.2 read_metaphlan, 2.4 read_picrust2, 2.4b read_picrust2_traits, Checkpoint B) and its decisions; checklist, slices table, dependency row 2.2 and decision 11 updated.

## 2026-10-03 (slice 2B start)
* **Refresh**: phase-0-foundation, phase-1-core, r-golden-parity, add-a-function, commit and `generated` only, against 4ef5314; the renorm warn-once test fix (2d47f8d) and the Hypothesis deadline (4ef5314) contradict nothing they state.

## 2026-10-03 (Phase 2, Checkpoint A knowledge)
- **Create**: [fn](modules/fn.md) (`load_hierarchy`, `func_glom`, `renorm`; HUMAnN 3.9 semantics, the final `x_kind`, guard and alignment rules) and [function-tables-as-mudata](decisions/function-tables-as-mudata.md) (options weighed, h5mu tree loss, mudata dependency, the Phase 4 task 4.1 forward note); both added to their `index.md`.
- **Update**: [core](modules/core.md): `_function.py`, `sum_pairs`, `replace_features` as a second Propagation implementer, `warn_user`'s third caller, and mudata among the third-party imports (review: three statements were false).
- **Update**: [io](modules/io.md): `read_humann`, its header unit rule, `path` in every error, no R golden; the "every reader goes through `make_treedata`" invariant and the Responsibility line now cover the MuData reader. `modules/index.md` entry reworded.
- **Update**: [datasets](modules/datasets.md): `toy_humann`, `enzyme` and its unpinned hash; description and `modules/index.md` entry no longer say "TreeData objects" only.
- **Update**: [pure-by-default](decisions/pure-by-default.md) table: `datasets.enzyme` returns a `pd.DataFrame`; the root `index.md` Modules line names `fn`.
- **Update**: [phase-2-function](roadmap/phase-2-function.md) design note 7 (frontmatter now edited) and the Checkpoint A Knowledge box; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.1 notes `function` + `function_by_taxon`.
- **Recheck**: phase-0, phase-1, engine-parity, module-boundaries, tree-access, add-a-function and cut-a-release flagged stale by the slice's `src/` and `pyproject.toml` changes; nothing they state was false, so only `commit` and `generated` moved. phase-0's pytest marker snippet now reads "R or HUMAnN golden files", as `pyproject.toml` does.

## 2026-10-03 (Phase 2, Checkpoint A fix pass)
- **Update**: [data-model-slots](contracts/data-model-slots.md) `func_glom` section: a mean, or a sum with a feature in several parents, sets `x_kind` to `abundance` (review I2).
- **Update**: [phase-2-function](roadmap/phase-2-function.md) design note 2: the nothing-maps guard applies to unstratified input only (review Minor 2).
- **Update**: [function-shape](contracts/function-shape.md) Examples bullet: "a function that reads a file may write a small temp file in its example" replaces "a reader's example", matching rules.md R8.2 (user-approved; review Minor 13).
- **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 8: a reader's R parity may come from a tool's golden files downstream (`io.read_humann` through the HUMAnN goldens; mia not added to the R image; user-approved, review F6). Its `description` now carries the HUMAnN clause its `contracts/index.md` entry already had (review Minor 14).
- **Update**: [data-model-slots](contracts/data-model-slots.md) `func_glom` section notes that `READS_UNMAPPED` passes through as in HUMAnN master while 3.9 sums it into `UNGROUPED` (review Minor 15); `paths` gains `src/biotapy/fn/**` (review Minor 14).
- **Update**: [regenerate-golden-files](playbooks/regenerate-golden-files.md) writes its checksum files to the git-ignored `build/`, not `/tmp` (review Minor 12).
- **Update**: [phase-2-function](roadmap/phase-2-function.md) ticks task 2.5's steps 2-7 (done in 03bb529 and its fix round); the roadmap `index.md` Phase 2 entry now carries the concept's description (review Minor 14).

## 2026-10-03 (Phase 2, task 2.12)
- **Update**: [data-model-slots](contracts/data-model-slots.md) lists `fn.renorm` as feature-changing and notes in convention 2 that renormalised stratified rows are shares of the community total; roadmap `phase-2-function` ticks task 2.12.

## 2026-10-03 (Phase 2, task 2.6)
- **Update**: [data-model-slots](contracts/data-model-slots.md) lists `fn.func_glom` as feature-changing, names `_core.replace_features` beside `feature_subset`, and gains the "Aggregation semantics (`func_glom`)" section; roadmap `phase-2-function` ticks task 2.6.

## 2026-10-03 (Phase 2, task 2.5 fix round 1)
* **Update**: no concept changed; `fn.load_hierarchy` now rejects empty files and lines with an empty cell, skips `#` lines, reads a UTF-8 BOM.

## 2026-10-03 (Phase 2, task 2.5)
* **Update**: ticked task 2.5 in [phase-2-function](roadmap/phase-2-function.md); `fn.load_hierarchy` reads local mapping files. Changed [no-bundled-kegg](decisions/no-bundled-kegg.md) (and its `decisions/index.md` entry): loaders read local files only, ENZYME is the one built-in download, KEGG and MetaCyc are never shipped or fetched.

## 2026-10-03 (Phase 2, task 2.5a fix round 1)
* **Update**: no concept changed; `datasets.enzyme` now checks that `enzyme.dat` and `enzclass.txt` are the same release and that `enzclass.txt` has class lines.

## 2026-10-03 (Phase 2, task 2.5a)
* **Update**: ticked task 2.5a in [phase-2-function](roadmap/phase-2-function.md); `datasets.enzyme` downloads ENZYME unpinned (`known_hash=None`, user decision 2026-10-03) and raises if `enzyme.dat` has no release line.

## 2026-10-03 (Phase 2, task 2.3b fix round 1)
* **Update**: [function-shape](contracts/function-shape.md) Examples bullet now also allows a reader's example to write a small temp file, matching rules.md R8.2.

## 2026-10-03 (Phase 2, task 2.3b)
* **Update**: [function-shape](contracts/function-shape.md) lets examples use `bt.datasets.toy_humann()` for function tables (rules.md R8.2 widened to match, user-approved); ticked task 2.3b in [phase-2-function](roadmap/phase-2-function.md).

## 2026-10-03 (Phase 2, task 2.3)
* **Update**: [data-model-slots](contracts/data-model-slots.md) adds the Function tables section and the `io.read_humann` exception to the `x_kind` convention; ticked task 2.3 in [phase-2-function](roadmap/phase-2-function.md).

## 2026-10-03 (Phase 2, task 2.1b)
* **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) records the approved mudata dependency; ticked task 2.1b in [phase-2-function](roadmap/phase-2-function.md).

## 2026-10-03 (Phase 2, task 2.1)
* **Update**: [phase-2-function](roadmap/phase-2-function.md) ticks task 2.1 (`_core.sum_pairs`, `_core.replace_features`).

## 2026-10-03 (Phase 2, task 2.0)
* **Update**: [r-golden-parity](contracts/r-golden-parity.md) (statement 1b, HUMAnN parity row and fixtures exception, enforcement bullet; description reworded) and [regenerate-golden-files](playbooks/regenerate-golden-files.md) (HUMAnN section, paths). Added the HUMAnN fixtures with their MIT notice and the HUMAnN 3.9 golden export; ticked task 2.0 in [phase-2-function](roadmap/phase-2-function.md).

## 2026-10-03 (Phase 2 plan)
* **Update**: [phase-2-function](roadmap/phase-2-function.md) carries the
  user-approved Phase 2 plan: resolved design notes (renorm, func_glom
  semantics, x_kind rules, MuData layout, layering, licence notices, roadmap
  corrections), slices 2A-2D with slice 2A in full TDD steps, review focus
  and exit gate; description reworded.

## 2026-10-03 (after 0.1.0)
* **Refresh**: [phase-0-foundation](roadmap/phase-0-foundation.md),
  [phase-1-core](roadmap/phase-1-core.md),
  [add-a-function](playbooks/add-a-function.md) and
  [cut-a-release](playbooks/cut-a-release.md), commit and `generated` only,
  against 45d0946; the setup-uv 10.2.0 bump (Dependabot PR #11) and the pixi
  install lines contradict nothing they state.

## 2026-10-03 (release 0.1.0)
* **Update**: Phase 1 closed after biotapy 0.1.0 reached PyPI (tag v0.1.0,
  release workflow run 37118911248). [phase-1-core](roadmap/phase-1-core.md)
  is `phase_state: done` with every Task 1.23 step and exit-gate item ticked;
  [phase-2-function](roadmap/phase-2-function.md) is `phase_state:
  in-progress`; the [roadmap index](roadmap/index.md) lists Phase 2 as active.
* **Update**: [phase-1-core](roadmap/phase-1-core.md) records the approved
  build floor `hatchling>=1.27` in its dependency table. Refreshed
  [phase-0-foundation](roadmap/phase-0-foundation.md),
  [module-boundaries](contracts/module-boundaries.md),
  [tree-access](contracts/tree-access.md),
  [add-a-function](playbooks/add-a-function.md) and
  [cut-a-release](playbooks/cut-a-release.md), commit and `generated` only,
  against 659f5a0; the hatchling floor and the performance page's suite time
  contradict nothing they state.
* **Refresh**: [phase-0-foundation](roadmap/phase-0-foundation.md),
  [phase-1-core](roadmap/phase-1-core.md),
  [engine-parity](contracts/engine-parity.md),
  [module-boundaries](contracts/module-boundaries.md),
  [tree-access](contracts/tree-access.md) and
  [add-a-function](playbooks/add-a-function.md), commit and `generated` only,
  against 695507e; the changes under their paths (version bump, license and
  classifier metadata, sdist excludes, benchmark docstring figure) contradict
  nothing they state.
* **Update**: [cut-a-release](playbooks/cut-a-release.md) gains Step 2b (update
  the README) and an sdist `pytest` check; phase-1-core ticks Task 1.23
  Steps 3-8 and its first five exit-gate items.

## 2026-10-03
* **Refresh**: [phase-1-core](roadmap/phase-1-core.md) and
  [engine-parity](contracts/engine-parity.md), commit and `generated` only,
  against 809c508; the only change under their paths was a `tl.py` docstring
  timing, which neither concept states.
* **Fix**: [io](modules/io.md) no longer quotes `untyped_calls_exclude = ["biom"]`
  (the real list has four entries); [r-golden-parity](contracts/r-golden-parity.md)
  description now says "computation", since `pl` functions are exempt. Checked
  against 2d0cab6.
* **Update**: refreshed the 14 concepts the Checkpoint D fixes made stale,
  against 2b9fc24, after checking each against the fixes. Content edits:
  [tl](modules/tl.md) (`permanova` runs OpenMP on one thread through
  threadpoolctl), [pl](modules/pl.md) (a missing layer names
  `adata = bt.pp.relative(adata)`; `heatmap`'s `ValueError` names `adata` or
  `layer=`; a richness group with no point gets no legend entry),
  [phase-1-core](roadmap/phase-1-core.md) (Checkpoint D review and knowledge
  boxes ticked). Commit and `generated` only: [core](modules/core.md),
  [io](modules/io.md), [phase-0-foundation](roadmap/phase-0-foundation.md),
  [phase-2-function](roadmap/phase-2-function.md),
  [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md),
  [engine-parity](contracts/engine-parity.md),
  [module-boundaries](contracts/module-boundaries.md),
  [r-golden-parity](contracts/r-golden-parity.md),
  [tree-access](contracts/tree-access.md),
  [add-a-function](playbooks/add-a-function.md),
  [cut-a-release](playbooks/cut-a-release.md).
* **Update**: three concept errors from the Checkpoint D review (M2):
  [function-shape](contracts/function-shape.md) says each public `pp`, `tl`
  and `pl` function has a purity test, not every test;
  [index](index.md)'s Modules line names all six module concepts;
  [data-model-slots](contracts/data-model-slots.md) marks `pp.clr` and
  `layers["clr"]` as Phase 3, not yet written.
* **Update**: [core](modules/core.md) lists `_slots.py:require_categorical`,
  the one validator for categorical groupings; [pl](modules/pl.md) and
  [tl](modules/tl.md) name it among their `_core` dependencies (Checkpoint D,
  I2).
* **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md)
  and [phase-1-core](roadmap/phase-1-core.md)'s dependency table record
  threadpoolctl (`>=3.5`, approved 2026-10-03), which limits `tl.permanova`'s
  OpenMP F-statistic to one thread (Checkpoint D, P1).
* **Create**: [pl](modules/pl.md), the Module concept for the plots (Task 1.22).
  Listed in [modules](modules/index.md), whose [io](modules/io.md) line now
  names phyloseq objects too.
* **Update**: refreshed every stale concept against 1ad037b (Task 1.22), after
  reading the code behind each statement. Content edits:
  [core](modules/core.md) (`add_provenance` stores numpy scalars through
  `.item()`; no asv import-time benchmark exists),
  [io](modules/io.md) (the rdata warning filter hides only the suffix messages;
  stale `_qiime2.py:15` line reference),
  [pp](modules/pp.md) (`tl` exists; every function is pure; `filter_samples`
  keeps stale ordinations; the prevalence boundary note),
  [tl](modules/tl.md) (the `faith_pd` presence copy per chunk),
  [data-model-slots](contracts/data-model-slots.md) (`pl` reads and writes no
  slot; `filter_samples` keeps stale ordinations),
  [function-shape](contracts/function-shape.md) (purity tests cover `pl`),
  [add-a-function](playbooks/add-a-function.md) (`pl` purity; the `R equivalent:`
  line feeds the Coming-from-R page),
  [phase-0-foundation](roadmap/phase-0-foundation.md) (mypy now covers
  `docs/extensions`, a root `conftest.py`, the `network` and `docs` CI jobs),
  [phase-1-core](roadmap/phase-1-core.md) (every task has TDD steps; Task 1.22
  ticked). Commit and `generated` only, no statement stale:
  [engine-parity](contracts/engine-parity.md),
  [module-boundaries](contracts/module-boundaries.md),
  [r-golden-parity](contracts/r-golden-parity.md),
  [tree-access](contracts/tree-access.md),
  [cut-a-release](playbooks/cut-a-release.md),
  [maintain-knowledge](playbooks/maintain-knowledge.md),
  [regenerate-golden-files](playbooks/regenerate-golden-files.md),
  [phase-2-function](roadmap/phase-2-function.md),
  [phase-3-stats](roadmap/phase-3-stats.md),
  [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md).
* **Update**: [tl](modules/tl.md) gained a Gotcha: `faith_pd` at 5,000 x 50,000
  takes about 48 s, nearly all of it scikit-bio re-indexing the tree per chunk
  (asv baseline, Task 1.21; `docs/performance.md`). Task 1.21 is ticked in
  [phase-1-core](roadmap/phase-1-core.md), with its exit-gate item "asv
  baselines recorded".
* **Update**: [engine-parity](contracts/engine-parity.md) statement 5 now names
  the benchmark dataset (`benchmarks/benchmarks/_data.py:synthetic`) and where
  the baselines are recorded (`docs/performance.md`). Task 1.21 added the asv
  suite.
* **Update**: [phase-1-core](roadmap/phase-1-core.md) Task 1.23 Step 5 gained a
  bullet: the release also switches `docs/tutorials/getting_started.md` from the
  GitHub install line to `pip install biotapy`. PyPI only has the 0.0.1
  placeholder, so the page says the GitHub line until then.

## 2026-10-02
* **Update**: [phase-1-core](roadmap/phase-1-core.md) Task 1.23 Step 5 now
  lists only the README changes that need the release itself: the user asked
  for the README to be brought up to date ahead of it, so its Status bullets,
  datasets, data model, pure-by-default and Coming from R text already
  describe slices 1C and 1D.

## 2026-09-28
* **Update**: Task 1.20 done: two MyST text notebooks,
  `docs/tutorials/getting_started.md` (GlobalPatterns, one ordination and one
  richness plot) and `docs/tutorials/phyloseq_analysis.md` (the phyloseq
  analysis vignette redone on GlobalPatterns, enterotype and esophagus
  through `bt.datasets.*`/pooch - no golden CSV read), added to
  `docs/tutorials/index.md`'s toctree. `docs/conf.py` now sets
  `nb_execution_mode = "cache"` and `nb_execution_raise_on_error = True`
  globally (myst-nb 1.4.0/jupyter-cache 1.0.1, already in the doc group);
  `quick_tour.md` dropped its page-level `execution_mode: force` override.
  A new CI job `docs` builds the docs with `uvx hatch run docs:build`
  (Read the Docs' command) against the `network` job's pooch cache and is
  now required by `check`; `tests/test_ci.py` gained
  `test_docs_job_builds_the_docs_with_the_pooch_cache` and
  `test_docs_job_blocks_merges` (RED: `KeyError: 'docs'`; GREEN: 8 passed).
  A clean `uvx hatch run docs:build` executed all three notebooks (7.0 s,
  16.2 s, 2.6 s) and `sphinx-build -W` succeeded with 8 `<img>` in the
  vignette page, checked by eye against the R vignette's figures. A probe
  cell (`bt.pl.scree(bt.datasets.esophagus())`) confirmed the gate bites:
  the build exited 2 with `WARNING: Executing notebook failed:
  CellExecutionError [mystnb.exec]` and a `KeyError` naming `bt.tl.pcoa`;
  the cell was removed and the build re-verified green. Ticked Task 1.20 in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.19b done: `[tool.coverage].report.fail_under = 90` in
  `pyproject.toml`, so CI's `test` job's `cov-report` step (hatch-test's
  `coverage report`) fails below 90% total line coverage - the proxy for
  R11.6, since coverage.py has no per-function view. `tests/test_ci.py`
  gained `test_coverage_below_90_percent_fails_the_test_job`, asserting the
  config key and that a `test` job step names `cov-report`. Measured: 560
  tests, `TOTAL` 99%; every `pl` file 100%; only
  `src/biotapy/datasets/_remote.py` (89%, downloads run only in the network
  job) is under 90%. Ticked Task 1.19b in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.19 done: a local Sphinx extension,
  `docs/extensions/coming_from_r.py`, parses each public function's
  docstring `Notes` section - `r_equivalents(doc)` raises `ValueError`
  unless there is exactly one `R equivalent:` line naming ```` ``pkg::fn`` ````
  items or `none` - and writes the Coming-from-R table
  (`docs/generated/coming_from_r_table.md`, git-ignored) from `rows()`/
  `render()` at the `builder-inited` hook, `write_table`. `rows()` also reads
  `docs/_data/r_idioms.toml` (stdlib `tomllib`, no new dependency) for the 19
  phyloseq accessors that are plain AnnData/TreeData code plus the 5 "not in
  0.1" rows, and raises if a call is mapped by both a docstring and the
  idioms file. `tests/test_docstrings.py` uses the same parser to check
  every public function across `datasets`, `io`, `pl`, `pp`, `tl` (25
  functions: 18 with one R item, 5 with two, 2 with `none`) has a parseable
  `R equivalent:` line, a `Guide:` link on the next line, and an `Examples`
  section; `tests/test_coming_from_r.py` checks the table covers all 31
  phyloseq functions from the phase's exit-gate list and marks the 5
  uncovered ones. The committed page `docs/coming_from_r.md` `{include}`s
  the generated fragment; `docs/conf.py`'s `exclude_patterns` excludes the
  fragment so `-W` does not fail on an unincluded document, and
  `docs/index.md`'s "User guide" toctree gains it. Built table: 48 rows (49
  `<tr>` with the header); the `{func}` roles resolved under `nitpicky`.
  Updated [function-shape](contracts/function-shape.md)'s "Enforced by" to
  describe the parser instead of only naming the test file. Ticked Task 1.19
  in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.18 done (slice 1D's first task): new top-layer package
  `bt.pl`, computing nothing and reading only slots `tl`/`pp` already write.
  `bt.pl.bar(adata, fill, *, x=None, layer=None, ax=None) -> Axes` sums the
  features of each `fill` group before drawing stacked bars (`phyloseq::plot_bar`
  without the per-feature outlines); `bt.pl.heatmap(adata, *, layer=None,
  ax=None) -> Axes` draws the table in `obs`/`var` order on phyloseq's
  `#000033`-`#66CCFF` log colour scale, zeros black
  (`phyloseq::plot_heatmap`); `bt.pl.richness(adata, metric, *, x=None,
  color=None, ax=None) -> Axes` scatters a stored `obs['alpha_<metric>']`
  (`phyloseq::plot_richness`); `bt.pl.ordination(adata, *, basis="pcoa",
  components=(1, 2), color=None, ax=None) -> Axes` scatters a stored PCoA or
  NMDS with axis-label percentages or a stress note
  (`phyloseq::plot_ordination`); `bt.pl.scree(adata, *, ax=None) -> Axes`
  bars the stored `proportion_explained` (`phyloseq::plot_scree`). A missing
  slot raises `KeyError` naming the `tl`/`pp` call that writes it; a numeric
  grouping column raises `TypeError` with the `.astype("category")` hint, as
  `tl.permanova` does. matplotlib (`>=3.8`, resolved 3.11.2) is now a runtime
  dependency, approved 2026-09-27; `pl` imports it only inside its functions
  (pyplot only when it must make a figure), so `import biotapy` still does not
  load it (`tests/pl/test_init.py` pins this; import time unchanged at
  1.2-1.3 s). rules.md R11.2 now excepts `pl` from the golden-test
  requirement (controller ruling 2026-09-27), recorded in
  [r-golden-parity](contracts/r-golden-parity.md) (Statement 7) and
  [add-a-function](playbooks/add-a-function.md) (Step 4). Corrected
  [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md):
  matplotlib was listed as already added in "Phase 0-1" but was not a
  dependency until this task; the Consequences section now names
  `import-without-extras` instead of an all-extras CI job, which does not
  exist yet. Refreshed `commit` to `806bede` on those two contracts/decisions
  and on [add-a-function](playbooks/add-a-function.md). Ticked Task 1.18 in
  [phase-1-core](roadmap/phase-1-core.md).

## 2026-09-27
* **Update**: Checkpoint C closed in [phase-1-core](roadmap/phase-1-core.md):
  the user approved pushing `phase-1c`, PR #8 merged on green CI (19/19,
  including the network job and Python 3.14) as `5f57d27`, and the user
  chose to continue to slice 1D, which is now being planned (R1.2a).
* **Update**: [pure-by-default](decisions/pure-by-default.md) re-verified by
  the user as amended in Task 1.17 (no `key_added`; `tl.permanova` has no
  `inplace`); `verified.at` updated. The user also chose to keep
  `obs["alpha_*"]` through feature changes, as documented (no contract change).
* **Creation** (Checkpoint C): [tl](modules/tl.md) Module concept for Slice 1C
  (the six `bt.tl` functions and `_beta.py:stored_distances`), linked from
  [modules/index.md](modules/index.md). Ticked Checkpoint C's Knowledge item
  in [phase-1-core](roadmap/phase-1-core.md); its review, push and
  user-review items stay open.
* **Update** (Checkpoint C fix C1): corrects the Task 1.15 entry's claim that
  `faith_pd` "also runs on relative abundances": it returned 0 for every
  sample on proportions, and weighted UniFrac put every pair 0 apart, because
  scikit-bio 0.7.4's tree code (`_nodes_by_counts`, shared by Faith PD and
  both UniFrac engines) casts abundances to int64. `_core.require_counts` now
  also raises when `X` holds non-integer values (`infer_x_kind`'s rule on
  `X.data`), which also covers `pp.rarefy` and `tl.alpha`'s
  `observed_features`/`chao1`; `tl.unifrac(weighted=True)` now calls it;
  `tl.alpha` gives `faith_pd` presence/absence, so Faith PD does run on
  relative abundances now. Updated
  [data-model-slots](contracts/data-model-slots.md) (the `x_kind`
  convention) and [core](modules/core.md) (`require_counts`).
* **Update**: Task 1.17 done (slice 1C's last task): `bt.tl.pcoa(adata, *,
  distance="braycurtis", n_components=10, inplace=False) -> tuple[pd.DataFrame,
  pd.DataFrame] | None` wraps `skbio.stats.ordination.pcoa`, asking for at most
  `n_obs - 1` axes so `proportion_explained` divides by the trace like
  `ape::pcoa`'s `Relative_eig`; `bt.tl.nmds(adata, *, distance="braycurtis",
  n_components=2, seed=None, inplace=False) -> tuple[pd.DataFrame, float] |
  None` wraps `sklearn.manifold.MDS` (non-metric SMACOF, `n_init=20`, as
  `vegan::metaMDS`'s default `try = 20`); `bt.tl.permanova(adata, grouping, *,
  distance="braycurtis", permutations=999, seed=None) -> pd.Series` wraps
  `skbio.stats.distance.permanova`, with no `inplace` (it returns a test
  result, not per-sample/per-pair values). New private
  `tl._beta.stored_distances(adata, key) -> DistanceMatrix`, shared by all
  three, raises `KeyError` naming the `bt.tl...` call that writes a missing
  `obsp` key and `ValueError` on NaN distances. New `src/biotapy/tl/_ordination.py`
  and `src/biotapy/tl/_permanova.py`. `_core.feature_subset` now keeps only
  `x_kind` and `provenance` in `uns["biotapy"]` (`_slots.py:KEPT_META`),
  dropping `pcoa`/`nmds` on a feature change. Added runtime dependency
  `scikit-learn>=1.8` (approved 2026-09-27): every `MDS` argument is explicit
  (`metric_mds=False`, `metric="precomputed"`, `n_init=20`, `init="random"`,
  `normalized_stress="auto"`, `random_state=<int drawn from as_generator>`),
  since scikit-learn 1.8 renamed `dissimilarity` to `metric` and changes
  `n_init`'s and `init`'s defaults in 1.9/1.10. Confirmed against the
  installed scikit-bio 0.7.4 and scikit-learn 1.9.1 (R2.2): `pcoa`'s
  `OrdinationResults.samples`/`eigvals`/`proportion_explained`;
  `permanova`'s `seed` accepts a `np.random.Generator` via
  `skbio.util.get_rng`, and its result `Series` is indexed by `"test
  statistic"`, `"p-value"`, `"sample size"`, `"number of groups"`, `"number
  of permutations"`; `MDS.fit_transform`/`stress_` raise no
  FutureWarning/DeprecationWarning with the pinned kwargs. `uv.lock` is
  gitignored in this repo (`/uv.lock`, "resolve fresh in CI and for users; no
  committed lockfile") so it is not staged, unlike the task brief's
  instruction; `uv sync` added only scikit-learn, joblib, threadpoolctl and
  cloudpickle, no other dependency moved. Amended the verified
  [pure-by-default](decisions/pure-by-default.md) decision under the
  controller's approved ruling: dropped the "key overridable by `key_added`"
  clause (no `tl` function has one) and scoped "every `tl` function supports
  both modes" to functions that return per-sample or per-pair values, since
  `tl.permanova` does not; added a Consequences bullet recording that
  exception. Updated [data-model-slots](contracts/data-model-slots.md) (new
  `uns["biotapy"]` keys `pcoa`/`nmds`, and the Propagation table's Feature-changing
  row) and [core](modules/core.md) (`feature_subset`'s invariant). Updated
  [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md):
  scikit-learn is no longer "pending approval". Import time unchanged at
  about 1.3 s (R10.1). Ticked Task 1.17's steps in
  [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.16 done: `bt.tl.beta(adata, *, metric="braycurtis",
  inplace=False) -> pd.DataFrame | None` wraps
  `skbio.diversity.beta_diversity` for `braycurtis` and `jaccard` (on
  presence/absence, matching `phyloseq::distance(physeq, "jaccard", binary =
  TRUE)`); `bt.tl.unifrac(tdata, *, weighted=False, normalized=True,
  inplace=False) -> pd.DataFrame | None` wraps the same function for
  `unweighted_unifrac`/`weighted_unifrac` via `_core.get_skbio_tree`. Both
  densify `X` once (R6.2) and, with `inplace=True`, write
  `obsp["braycurtis" | "jaccard" | "unweighted_unifrac" |
  "weighted_unifrac"]` and return `None`
  ([pure-by-default](decisions/pure-by-default.md)). New
  `src/biotapy/tl/_beta.py`. Confirmed against the installed scikit-bio 0.7.4
  (R2.2): `beta_diversity(metric, counts, ids=..., taxa=..., tree=...,
  **kwargs) -> DistanceMatrix`; `DistanceMatrix.to_data_frame()` has the ids
  on both index and columns; `"jaccard"` is in `_qualitative_metrics` and is
  auto-qualified to presence/absence; `weighted_unifrac`'s own default is
  `normalized=False`, so `unifrac` always passes it explicitly. Two all-zero
  samples are `NaN` apart under Bray-Curtis but `0.0` under both UniFracs and
  under Jaccard, with no `RuntimeWarning` in any case. New tests in
  `tests/tl/test_beta.py` (unit, Hypothesis, purity, the multifurcating-root
  and post-filtering path-length cases) and `tests/tl/test_beta_golden.py`
  (against `beta_*.csv.gz` and `unifrac_*.csv.gz` from Task 1.15a, on
  GlobalPatterns and esophagus); both pass. Appended to
  `docs/guide/diversity.md` and `docs/api.md`. No concept needed a content
  change: the `obsp` keys were already in
  [data-model-slots](contracts/data-model-slots.md); ticked Task 1.16's
  steps in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.15 done: `bt.tl.alpha(adata, *, metrics=(...), inplace=False)
  -> pd.DataFrame | None` wraps `skbio.diversity.alpha_diversity` for
  `observed_features`, `shannon` (natural log), `simpson` (Gini-Simpson),
  `chao1` (bias-corrected) and `faith_pd` (via `_core.get_skbio_tree`).
  `observed_features` and `chao1` require raw counts through
  `_core.require_counts`; the others also run on relative abundances.
  scikit-bio needs dense rows, so `X` is densified in chunks of at most
  `2**20` values (R6.2). With `inplace=True` writes
  `obs["alpha_<metric>"]` and returns `None`
  ([pure-by-default](decisions/pure-by-default.md)). New package `bt.tl`,
  with `src/biotapy/tl/__init__.py` (imports only, R4.1) and
  `tl/_alpha.py`. Confirmed against the installed scikit-bio 0.7.4
  (R2.2): `alpha_diversity(metric, counts, ids=..., **kwargs) -> pd.Series`;
  `shannon(base=None)` defaults to natural log since 0.6.1, equal to
  `base=math.e`; `chao1(bias_corrected=True)` is already the default;
  `simpson` is `1 - sum(p**2)`; `faith_pd(counts, taxa, tree)` takes
  `taxa=`/`tree=` as keyword args; an all-zero row gives 0 for
  `observed_features`/`chao1`/`faith_pd` and NaN for `shannon`/`simpson`,
  with no warning. New tests in `tests/tl/test_alpha.py` (unit,
  Hypothesis, purity, memory-chunking) and
  `tests/tl/test_alpha_golden.py` (against `alpha.csv.gz` and
  `alpha_faith_pd.csv.gz` from Task 1.15a); both pass. Added
  `docs/guide/diversity.md` and a Tools section to `docs/api.md`. Updated
  [data-model-slots](contracts/data-model-slots.md) (`require_counts`
  convention now names `tl.alpha`) and [pp](modules/pp.md) (`tl` ownership
  note points at Slice 1C instead of "later phases"); ticked Task 1.15's
  steps in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.15c done: `_core.get_skbio_tree(adata: AnnData) ->
  skbio.TreeNode` converts the phylogeny in `vart["phylo"]` to a scikit-bio
  `TreeNode` via `nx.bfs_edges`, rooted where the networkx tree is drawn. A
  root with more than two children (the toy tree's has three) keeps its first
  child and moves the rest under one new zero-length node, which changes no
  root-to-tip path length; NaN branch lengths pass through unchanged. A plain
  `AnnData` raises `TypeError` naming what it needs; a `TreeData` without
  `vart["phylo"]` raises `KeyError` from `get_tree`. Exported from `_core`;
  used by the upcoming `tl.alpha` (faith_pd) and `tl.unifrac`, kept in
  `_tree.py` because only that module may import networkx (R4.5). New tests
  in `tests/core/test_tree.py`. Confirmed against the installed scikit-bio
  0.7.4 (`TreeNode.__init__`, `append`, `extend`, `tips`, `find`, `distance`)
  that `append`/`extend` reparent nodes rather than copying them. Updated
  [tree-access](contracts/tree-access.md) (statement 2, a Gotcha) and
  [core](modules/core.md) (entry point, dependencies); ticked Task 1.15c's
  steps in [phase-1-core](roadmap/phase-1-core.md).
* **Update**: Task 1.15b fix round 1: named `esophagus()` alongside
  `global_patterns()`/`enterotype()` in two places the esophagus change had
  left stale - [datasets](modules/datasets.md)'s Gotchas doctest-`+SKIP`
  bullet and `docs/contributing.md`'s network-test paragraph. No behavior
  change.
* **Update**: Task 1.15b done: `bt.datasets.esophagus() -> TreeData` - 3
  esophageal biopsies (samples `B`, `C`, `D`) x 58 OTUs, with a tree in
  `vart["phylo"]` and no taxonomy or sample data, read through
  `bt.io.read_phyloseq`. `_remote.py:_REGISTRY` gains `esophagus.RData`
  (sha256 `0b06d9c3...`, 1,840 B) at the same pinned phyloseq commit;
  `datasets/__init__.py` exports it. `docs/api.md`, `docs/guide/datasets.md`
  and `docs/guide/reading_data.md` updated for the third loader. 1.16's
  UniFrac golden test will use it. [datasets](modules/datasets.md) updated;
  ticked in [phase-1-core](roadmap/phase-1-core.md).
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
* **Update** (Checkpoint C): the whole-slice review (opus) of slice 1C and its
  fix pass (4a5adaf..0dd46d5, plus this follow-up) are done: `tl.permanova`'s
  stray f-string prefix, `_core.require_counts`'s message now naming the NaN
  case, and `tl.alpha`'s docstring/`docs/guide/diversity.md` now counting
  `faith_pd`'s int64 presence copy. Ticked Checkpoint C's review item in
  [phase-1-core](roadmap/phase-1-core.md); its push and user-review items
  stay open.

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
