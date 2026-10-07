---
type: Phase
title: Phase 3 - Differential abundance consensus (0.3)
description: CLR and PhILR transforms; ANCOM-BC2 and LinDA natively and ALDEx2 and MaAsLin 3 through an optional R bridge, all behind one result schema; a consensus table and plot of where methods agree.
tags: [roadmap, da, pp]
status: stable
release: "0.3"
phase_state: in-progress
effort: ~4 weeks part-time
depends_on: [/roadmap/phase-2-function.md]
paths: ["src/biotapy/da/**", "src/biotapy/pp/**", "src/biotapy/pl/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-07T15:30:00Z }
commit: d88e3f6
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: maaslin3
    resource: https://www.nature.com/articles/s41592-025-02923-9
    title: MaAsLin 3, Nature Methods
  - id: philr
    resource: https://doi.org/10.7554/eLife.21887
    title: Silverman et al. 2017, A phylogenetic transform enhances analysis of compositional microbiota data, eLife
  - id: linda
    resource: https://doi.org/10.1186/s13059-022-02655-5
    title: Zhou et al. 2022, LinDA, Genome Biology
  - id: nearing
    resource: https://www.nature.com/articles/s41467-022-28034-z
    title: Nearing et al. 2022, Microbiome differential abundance methods produce different results across 38 datasets, Nature Communications
  - id: pelto
    resource: https://academic.oup.com/bib/article/26/2/bbaf130/8093585
    title: Pelto et al. 2025, Elementary methods provide more replicable results in microbial differential abundance analysis, Briefings in Bioinformatics
---

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Every slice has full TDD steps.

**Goal:** 0.3 adds compositional statistics: CLR and PhILR transforms equal to
vegan and philr, two native differential abundance (DA) methods (ANCOM-BC2,
LinDA) and two R-bridged ones (ALDEx2, MaAsLin 3) that return one result
schema, and a consensus table and plot that show where the methods agree.[^spec]

**Architecture:**
- Transforms are `pp` verbs that add to a copy: `pp.clr` writes
  `layers["clr"]`, `pp.philr` writes `obsm["X_philr"]` (a samples x
  balances table) and keeps the tree beside it.
- `da` is a new top-layer subpackage (beside `pl` and `ml`). Every method
  takes `(adata, group, *, covariates=(), reference=None, ...)` and returns a
  `pd.DataFrame` indexed by feature with the same columns (`effect` in log2,
  `se`, `pvalue`, `qvalue` by Benjamini-Hochberg, `direction`, `method`,
  `contrast`). No method filters features or sees a formula string.
- `da.consensus(results, *, alpha, min_methods)` combines result tables the
  user computed; it never runs a method. `pl.consensus` draws its table.
- R-only methods (ALDEx2, MaAsLin 3) run through rpy2, an optional extra `r`,
  imported inside the function. Golden files come from the pinned R image
  (`tests/r/Dockerfile`), extended one package per commit.

**Tech stack:** anndata 0.13 · treedata 0.3.1 · scikit-bio 0.7.4 (`clr`,
`tree_basis`, `ancombc2`) · numpy · scipy 1.18 (`stats.t`,
`false_discovery_control`) · pandas 3 · matplotlib · rpy2 (extra `r`, slice 3C)
· R 4.5.3 / Bioconductor 3.22 golden image (vegan, philr, ANCOMBC,
MicrobiomeStat, ALDEx2, maaslin3) · pytest/hypothesis.

**Spec:** [plan.md](../../plan.md). Contracts that bind every task:
[function-shape](/contracts/function-shape.md),
[data-model-slots](/contracts/data-model-slots.md),
[module-boundaries](/contracts/module-boundaries.md),
[tree-access](/contracts/tree-access.md),
[r-golden-parity](/contracts/r-golden-parity.md). Decisions:
[r-bridge-before-ports](/decisions/r-bridge-before-ports.md),
[pure-by-default](/decisions/pure-by-default.md),
[optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md).
Research (2026-10-05, in the session scratchpad `phase3-research/`): A methods
and libraries, B R bridge and CI, C consensus, schema and benchmark data. Facts
marked [V] there were run or read at source; this plan relies only on those,
and names the rest under "Risks" and in the slice outlines.

**How slice 3A was checked.** Every file in slice 3A was written into a
scratch clone of the repository at `c015d3d` (master, after 0.2.0) and
replayed as one commit per task (3.0 as its two commits). Each committed state
was gated:

| Task state | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` |
|---|---|---|---|---|
| 3.0 (`test(golden)` commit) | passed | 987 passed, 23 deselected | 30 passed | build succeeded |
| 3.1 | passed | 1006 passed, 24 deselected | 31 passed | build succeeded |
| 3.2 | passed | 1026 passed, 25 deselected | 32 passed | build succeeded |

- prek covers ruff 0.16.9 check and format, `mypy --strict`, import-linter
  and pyproject-fmt. Every run exported `BIOTAPY_DATA_DIR` to a scratch pooch
  cache; `~/.cache/biotapy` was checked absent afterwards.
- Coverage on the final state: `pp/_transform.py`, `pp/_philr.py` and
  `_core/_tree.py` are each 100%.
- The two property tests also passed under Hypothesis seeds 1, 2 and 3.
- The RED counts in each "expect failure" step were reproduced by running the
  new tests against the previous task's commit.
- The golden image was rebuilt from the edited `tests/r/Dockerfile` (tagged
  `biotapy-golden-p3a` locally, never pushed) and `export_golden.R` run twice:
  bit-identical, and every existing golden file unchanged. The first build
  failed at the `requireNamespace("philr")` guard (the `fs` binary needs
  `libuv1`); the Dockerfile below carries the fix.
- APIs checked in the installed versions: scikit-bio 0.7.4 (`clr`,
  `tree_basis`, `TreeNode(name, children=...)`, `postorder`, `preorder`,
  `tips(include_self=True)`, `TreeNode.prune`, which moves a merged child to
  the end of its parent's list), anndata 0.13.4 (a `DataFrame` in `obsm`
  round-trips through h5ad and h5td), treedata 0.3.1, R philr 1.36.0
  (`pseudocount=` exists), vegan 2.7.3 (`decostand(x, "clr", pseudocount=)`).

# Goal
Differentiator 2: several DA methods behind one schema, and a report of where
they agree.[^spec]

# Entry criteria
- Phase 2 exit gate passed and 0.2 released. (Done: v0.2.0, 2026-10-05.)
- [r-bridge-before-ports](/decisions/r-bridge-before-ports.md) confirmed by
  the user (verified 2026-09-26).

# Design notes (resolved)

Each note answers one design question from the planning brief, with the
reason. Notes marked **(user)** change a contract, a rule, a dependency, CI or
a roadmap signature and are repeated under "Decisions for the user".

1. **Slices and order: 3A transforms, 3B native DA and consensus, 3C R
   bridges and CI, 3D docs and release. (user)**
   - 3A needs no `da` code and no new dependency. Its goldens need only the R
     image, so it starts at once.
   - 3B gives a complete feature without R: the schema, ANCOM-BC2, LinDA,
     `da.consensus` and `pl.consensus`. Consensus comes before the bridges
     because it consumes result tables and needs no method of its own; two
     native methods are enough to test it.
   - 3C adds the extra `r`, ALDEx2, MaAsLin 3 and the CI job that runs them;
     it plugs into the schema 3B fixed.
   - 3D writes method pages, the DA guide and the exit-gate notebook, the
     Coming-from-R check, benchmarks, knowledge and the release.

2. **`pp.clr(adata, *, pseudocount=0.5) -> AnnData` writes `layers["clr"]`.**
   - scikit-bio 0.7.4's `clr(mat)` takes a dense, strictly positive matrix
     (sparse input raises `TypeError`, an AnnData raises `ValueError`; research
     A section 2). The wrapper adds the pseudocount to a dense float64 copy of
     `X` and calls it: a thin wrapper (R2.1), densified once (R6.2), with the
     memory cost in `Notes`. CLR has no zeros, so the layer is a dense
     `ndarray`; the contract says so.
   - Pseudocount convention: added to **every** value, zeros or not. This is
     what `mia::transformAssay(method = "clr", pseudocount = x)` does (its
     `.apply_pseudocount`) and what it delegates to,
     `vegan::decostand(x, "clr", pseudocount = x)`. The golden target is
     therefore vegan, already in the image: no image change for CLR. mia's
     default `pseudocount = FALSE` fails on zeros; the roadmap's 0.5 is
     LinDA's default and suits counts.
   - The roadmap's `layer=None` is dropped (R2.3): no test or use case needs
     to read a layer, and CLR is scale invariant, so `layers["relative"]`
     would only differ by the pseudocount's scale. **(user: signature)**
   - A pseudocount larger than the smallest non-zero value of `X` warns
     (`UserWarning`). With counts the smallest value is 1 and 0.5 is silent;
     with relative abundances or CPM the default 0.5 would swamp every rare
     feature without an error. Negative or non-finite `X` raises; zeros with
     `pseudocount=0` raise naming the fix.

3. **`pp.philr(tdata, *, pseudocount=0.5) -> TreeData` writes
   `obsm["X_philr"]`. (user)** PhILR is Silverman et al.'s phylogenetic
   isometric log-ratio transform.[^philr]
   - **Slot.** The balances are a coordinate system for samples (D-1 columns
     for D features), which is what `obsm` holds: "ordinations and
     embeddings", keys `X_<name>`. A new AnnData of balances would break the
     data-model contract twice: `X` is a sparse samples x features table with
     an `x_kind`, and balances are dense, signed and have no `x_kind`. In
     `obsm` the tree stays beside the balances, so every column name is a
     node of `vart["phylo"]`, and a feature change drops them like every
     other derived slot. The value is a `pd.DataFrame` (samples x nodes), so
     the node names travel with it; it round-trips through h5ad and h5td
     (checked).
   - **Basis.** scikit-bio's `tree_basis(tree)` builds the sparse orthonormal
     basis (R2.1). It needs a strictly bifurcating tree, puts a node's
     **first child in the denominator** (philr puts it in the numerator) and
     orders columns as `tree.tips()`. The wrapper negates the product, so a
     balance is positive when the first child is more abundant, as in
     `philr::philr`; the golden test checks that every balance's numerator
     and denominator taxa equal R's, which pins the sign.
   - **Tree rules.** Subsetting a TreeData keeps one-child nodes
     ([tree-access](/contracts/tree-access.md) gotcha), so any filtered table
     has them. They define no balance and are skipped, as `ape::drop.tip`
     removes them in R. scikit-bio's `TreeNode.prune()` would do it but moves
     the merged child to the end of its parent's list, which flips signs; the
     wrapper rebuilds the tree in one postorder pass instead. A node with
     three or more children has no single balance: `pp.philr` raises naming
     it, **including a three-child root** (an unrooted tree). That needs
     `_core.get_skbio_tree(adata, *, split_root=False)`: by default the
     helper resolves a wide root for Faith PD and UniFrac, which changes no
     path length there but here would invent a root balance. R's philr
     requires a rooted binary tree too (mia warns and proceeds; philr's
     `phylo2sbp` then silently ignores a third child). The alternative,
     resolving the root as UniFrac does, is listed under the decisions.
   - **Weights.** Uniform part and ILR weights only, the code defaults of
     `philr::philr`. The paper's `enorm.x.gm.counts` / `blw.sqrt` options are
     not offered (R2.3: no use case yet; each is about 20 lines over the same
     basis when one appears).
   - **Golden.** The R image gains Bioconductor philr 1.36.0 (GPL-3, run only
     in the image). The export takes the 293 GlobalPatterns taxa with more
     than 3 reads in over half the samples (`prune_taxa`, so R's tree is
     binary and rooted), runs `philr(x, tree, pseudocount = 0.5)` and writes
     the balances and each balance's sign partition. The Python test subsets
     GlobalPatterns to the same taxa (keeping one-child nodes, which tests the
     skip), matches balances by partition and compares values elementwise.
   - `toy()`'s root has three children, so the docstring example uses
     `toy()[:, :6]` (f1-f6 sit under a binary root) and says why.

4. **No formula strings in 0.3: `group`, `covariates`, `reference`. (user:
   signatures)** The roadmap's `da.<method>(adata, formula, ...)` leaves
   undefined which coefficient a one-row-per-feature table reports, and
   patsy (which scikit-bio's ANCOM-BC2 uses) silently drops rows with a
   missing value and picks the alphabetically first level as reference (the
   research run got the sign of CD vs nonIBD flipped this way). Every method
   therefore takes:
   - `group: str`, the `obs` column whose effect is reported. A
     categorical, string or bool column must have exactly two levels (the
     effect is the other level against `reference`); a numeric column gives a
     slope. More levels raise `ValueError` naming them and saying to subset.
   - `covariates: Sequence[str] = ()`, `obs` columns adjusted for (numeric as
     is, categorical as indicator columns against their first level).
   - `reference: str | None = None`, the reference level of a categorical
     `group`. `None` takes the first category of `pd.Categorical(obs[group])`
     (declared order for a categorical, sorted for strings), the rule of R
     factors and patsy, and the result's `contrast` column (note 5) says which
     way round the effect is. Passing it for a numeric `group` raises.
   - Missing values in `group` or a covariate raise, naming the column and
     the count (no silent row drops). `obs` column names that are not
     identifiers are quoted where a method builds a formula (patsy `Q()`, R
     backticks), so users never quote.
   - Random effects and interactions are out of 0.3: ALDEx2's t-test and
     scikit-bio's `ancombc2` cannot take them, so a consensus over them would
     compare different models.
   - Consequence: patsy, formulaic and statsmodels are never imported by
     biotapy (LinDA builds its design matrix with pandas; note 7), so none is
     declared (note 11).

5. **The result schema (`da/_schema.py`).** Every `da` method returns a
   `pd.DataFrame` indexed by `var_names` (index name `feature`), one row per
   feature, in `var_names` order:

   | Column | dtype | Meaning |
   |---|---|---|
   | `effect` | float64 | log2 fold change of `group` (non-reference level vs reference, or per unit / per SD of a numeric `group`; note 7 for LinDA) |
   | `se` | float64 | standard error of `effect` (log2); NaN for methods without one (ALDEx2) |
   | `pvalue` | float64 | the method's p-value |
   | `qvalue` | float64 | Benjamini-Hochberg over the features the method tested |
   | `direction` | int8 | sign of `effect`: -1, 0 or 1; 0 when `effect` is NaN |
   | `method` | str | `"ancombc2"`, `"linda"`, `"aldex2"` or `"maaslin3"` |
   | `contrast` | str | `"<level> vs <reference>"`, or the column name for a numeric `group` |

   - **Units:** log2 everywhere. ANCOM-BC2 reports natural logs (scikit-bio
     0.7.4 renamed its column `Log(FC)`); its `effect` and `se` are divided
     by `ln 2`. LinDA, ALDEx2 (`diff.btw`) and MaAsLin 3 (abundance `coef`)
     are already log2.
   - **One correction:** `qvalue` is BH everywhere. scikit-bio's and R
     ANCOMBC's default is Holm, so the wrapper passes `p_adjust="bh"`; LinDA
     uses `scipy.stats.false_discovery_control(p, method="bh")`; ALDEx2's
     `we.eBH` and MaAsLin 3's `qval_individual` (abundance model only) are
     BH. [Superseded by slice 3C decisions 3 and 4: neither is BH of the
     group's p-values; both bridges carry the p-value and `_schema.result`
     computes the BH `qvalue`.] Mixed Holm/BH q-values would make "significant" mean different
     things per column.
   - **NaN:** a feature a method could not test keeps its row with NaN
     `effect`, `se`, `pvalue`, `qvalue` and `direction` 0. Methods never drop
     features and never filter by prevalence: filtering is the user's one
     `pp.filter_features` call before any method, so every method tests the
     same features and BH families match. Each wrapper turns off its
     library's own filter (R ANCOMBC `prv_cut = 0`, MaAsLin 3
     `min_prevalence = 0`; LinDA's prevalence and library-size filters are
     not ported).
   - Methods' own significance flags (`Signif`, `reject`, `passed_ss`) are
     not carried: consensus recomputes calls from `qvalue`.
   - `contrast` is added to the roadmap's columns because `da.consensus` must
     refuse tables that compare different things, and the plot labels the
     direction with it.
   - `validate_result(frame) -> pd.DataFrame` checks columns, dtypes, a
     unique index, `qvalue`/`pvalue` in [0, 1] or NaN and `direction ==
     sign(effect)`. Every method calls it on its output and `da.consensus`
     on its inputs, which is how the schema is tested through the public API
     (R4.9).

6. **"Agree" and consensus: a two-step API. (user: signature)**
   - **Definition** (written as the Decision concept
     `decisions/da-consensus-agreement.md` in task 3.8): per feature `f` and
     method `m`, `called(f, m) = qvalue < alpha` (strict, so `q == alpha` is
     not called; LinDA's own `reject` uses `<=`, which is why calls are
     recomputed). `n_tested(f)` counts methods with a finite `pvalue`;
     `n_significant(f)` the methods that call it. `consensus(f)` is true when
     `n_significant >= min_methods` **and** every calling method has the same
     non-zero `direction`. `conflict(f)` is true when calling methods
     disagree in sign; such a feature is never a consensus hit. NaN is "not
     tested", never "not significant".
   - **API:** `da.consensus(results, *, alpha=0.05, min_methods=2) ->
     pd.DataFrame`, where `results` is a sequence of schema tables the user
     computed:
     ```python
     results = [bt.da.ancombc2(t, "host"), bt.da.linda(t, "host"), bt.da.aldex2(t, "host", seed=0)]
     table = bt.da.consensus(results)
     bt.pl.consensus(table)
     ```
     The roadmap's one-call `da.consensus(adata, formula, *, methods,
     alpha, min_methods, n_jobs)` needs, with `group`, `covariates`,
     `reference` and a `seed` for ALDEx2 and MaAsLin 3, eight arguments
     (R5 caps them at six, `PLR0913`). It would also need a policy for one
     method failing (R7.4 forbids dropping it silently), a way to pass
     per-method options, and it would hide which methods were chosen. The
     two-step form keeps the method list in the user's code, written before
     the results are seen, which is what the literature asks of a consensus
     report (Nearing et al. recommend a consensus;[^nearing] Pelto et al. and
     the OMA book warn that trying methods until one agrees is selective
     reporting[^pelto]). The guide says so.
   - `n_jobs` and joblib are dropped (R2.3, R10.1): the native methods take
     under 2 s on the benchmark (research C section 6); rpy2 in joblib workers
     must start R in each process; the R methods have their own `cores`.
   - Inputs are validated: at least `min_methods` tables, `1 <= min_methods
     <= len(results)`, unique `method` values, one `contrast` across tables
     (else `ValueError` naming them), `0 < alpha < 1`. Features are the union
     of the indexes (a table missing a feature counts as not tested).
   - Output, indexed by feature: `effect_<method>` and `qvalue_<method>` per
     input in input order, then `n_tested`, `n_significant`, `direction`
     (the shared sign of the calls, 0 when none or in conflict), `consensus`,
     `conflict`; `attrs` are not used.

7. **Methods in 0.3.**
   - **ANCOM-BC2: `da.ancombc2(adata, group, *, covariates=(),
     reference=None) -> pd.DataFrame`. (user: rename)** It wraps scikit-bio
     0.7.4's `ancombc2` (Lin & Peddada 2024), not `ancombc` (v1): v2 accepts
     zeros (treated as missing with `pseudocount=0`, as R's `ancombc2`
     default `pseudo = 0`), while v1 requires strictly positive input. The
     name follows the method and the R function the golden comes from
     (`ANCOMBC::ancombc2`), so `da.ancombc` becomes `da.ancombc2`. The
     roadmap's `alpha` is dropped from every method: it only sets the
     `Signif` flag the schema does not carry. scikit-bio rejects sparse
     input (`np.asarray` of a sparse matrix), so `X` is densified once
     (R6.2); `require_counts` applies. Golden: R `ancombc2(..., p_adj_method
     = "BH", prv_cut = 0, lib_cut = 0, pseudo_sens = FALSE, struc_zero =
     FALSE, iter_control = list(tol = 1e-5, max_iter = 100))`, the settings
     mapped to scikit-bio's. Whether scikit-bio equals R elementwise is
     unknown (R's `conservative = TRUE` variance has no scikit-bio argument);
     task 3.4 measures it and uses elementwise tolerances if it holds, or the
     contract's DA rule (sign agreement and Spearman rank correlation of
     effects) with the measured numbers if not.
   - **LinDA: `da.linda(adata, group, *, covariates=(), reference=None) ->
     pd.DataFrame`, native.**[^linda] Fixed effects only (a formula with `(1|id)`
     needs lme4's Satterthwaite df, which no installed Python library
     gives). The algorithm, from MicrobiomeStat 1.4 `linda` (read at source;
     reimplemented, no GPL code copied): numeric covariates and a numeric
     `group` are scaled (`scale()`, so a numeric effect is per SD);
     adaptive zero handling (if any coefficient of `log(library size) ~
     design` has p <= 0.1, zeros are imputed by library size, else 0.5 is
     added); `W = log2(Y)` centred per sample; one OLS fit per feature
     (`numpy.linalg.lstsq`, df = n - p); the bias of each coefficient is the
     mode of `sqrt(n) * beta` (modeest's `mlv(method = "meanshift", kernel =
     "gaussian")`, ported: Gaussian mean shift from the shorth start with
     `bw.nrd0` bandwidth); `effect = beta - bias`, `p = 2 * t.sf(|effect /
     se|, df)`, BH. Winsorisation: MicrobiomeStat defaults to `is.winsor =
     TRUE, outlier.pct = 0.03`, the original LinDA package to none; 0.3 does
     not winsorise (R2.3), and the golden calls `MicrobiomeStat::linda(...,
     is.winsor = FALSE)`, so it is exact. **(user)** numpy and scipy suffice
     (research A's prototype matched `statsmodels.OLS` standard errors), so
     statsmodels, which roadmap 3.5 asked about, is not needed.
   - **ALDEx2: `da.aldex2(adata, group, *, mc_samples=128, reference=None,
     seed=None) -> pd.DataFrame`, rpy2 bridge (extra `r`).** No Python port
     exists. scikit-bio's `dirmult_ttest` resembles ALDEx2 but is not it:
     different outputs, no parity tested, and it would put another method
     under ALDEx2's name. **(user)** A two-level `group` only (ALDEx2's
     Welch t-test; its `glm` route for covariates is out of 0.3, so no
     `covariates` argument). `effect = diff.btw` (median CLR difference,
     log2; ALDEx2's standardised `effect` is not log2), `se` NaN, `pvalue =
     we.ep`, `qvalue = we.eBH` [superseded by slice 3C decision 3: `qvalue`
     is BH of `we.ep`]. Seeding: `as_generator(seed)` draws one
     integer, passed to R's `set.seed` before the call (R's RNG is global to
     the embedded R; documented in `Notes`) [superseded by the Checkpoint C
     fix: the call restores R's `.Random.seed`, or its absence, on exit];
     `useMC = FALSE`.
   - **MaAsLin 3: `da.maaslin3(adata, group, *, covariates=(),
     reference=None, seed=None) -> pd.DataFrame`, rpy2 bridge.**[^maaslin3] Abundance
     model only (`evaluate_only = "abundance"`): the prevalence model reports
     log-odds, which cannot share an `effect` column with log2 fold changes.
     `min_prevalence = 0`, plots off, output to a temporary directory read
     back from `all_results.tsv`. `effect = coef`, `se = stderr`,
     `pvalue = pval_individual`, `qvalue = qval_individual`. [Superseded by
     slice 3C decisions 4-6: the results come from the R return value, not
     `all_results.tsv`; `effect` is the coefficient minus the median
     coefficient; `qvalue` is BH over the group's p-values.]
   - The `alpha` and `formula` arguments of roadmap 3.4-3.7 are replaced as
     above. **(user: signatures)**

8. **R side. (user: dependency, CI)**
   - **Extra `r = ["rpy2>=3.6.8"]`** (slice 3C). rpy2 3.6.8 is GPL-2.0+;
     ALDEx2 and MicrobiomeStat are GPL-3, ANCOMBC Artistic-2.0, maaslin3
     MIT, R itself GPL. biotapy ships none of it: rpy2 is imported inside the
     two bridge functions through `import_optional("rpy2.robjects",
     extra="r")`, and the user installs R and the R packages. A BSD-3 package
     that imports an optional GPL package on the user's machine is the
     common pattern; bundling them (a wheel or image that includes rpy2 or
     R) would not be. Not legal advice; the user decides. The docs say the
     extra pulls GPL software.
   - rpy2-rinterface 3.6.7 ships no Linux wheels (the sdist builds against
     the installed R, needs R headers and a compiler) and rpy2's classifiers
     stop at Python 3.13: the bridges are documented Linux/macOS-first, and
     the CI job runs Python 3.13.
   - A missing R may surface as an error other than `ImportError` when
     rpy2 initialises ([UNVERIFIED], no R on the planning host).
     `da/_r.py` (one module uses it, so it stays in `da`, R4.3) imports rpy2
     through `import_optional`, then loads each R package with
     `rpy2.robjects.packages.importr` and turns a failure into an
     `ImportError` naming the package and `BiocManager::install(...)`. No
     silent fallback (R7.4).
   - Conversion: dense only (rpy2 has no scipy.sparse -> `Matrix` converter;
     `rpy2-Matrix` is GitHub-only 0.0.3). `X` is densified once per call
     (R6.2), features as rows for ALDEx2, through a local
     `(default_converter + pandas2ri.converter).context()`, never global
     activation.
   - **Golden image additions**, each its own Dockerfile commit with its
     `requireNamespace` guard: philr (3A), ANCOMBC and MicrobiomeStat (3B),
     ALDEx2 and maaslin3 (3C), all from Bioconductor 3.22 / the image's P3M
     snapshot, versions recorded in `tests/golden/VERSIONS.txt`. Goldens are
     generated offline in the image as today; CI never builds it.
   - **Tests.** Native goldens (`da.linda`, `da.ancombc2`) carry `golden` +
     `network` and run in the existing network job. Bridge tests carry `r`
     (deselected by default) and compare the bridge's output with the
     golden from the image: a check that the glue (conversion, column
     mapping, sign, seeding) is right. The glue's pure-Python parts (schema
     mapping from a canned R table, argument checks, the missing-R error)
     are tested without R by monkeypatching the private R call, the
     documented exception already used for `datasets._enzyme._fetch`, so
     coverage holds without R. [Superseded by slice 3C decision 13: the
     tests replace the rpy2 modules (`fake_rpy2` in `tests/da/conftest.py`),
     so the bridges' own code runs unchanged.]
   - **CI job `r-bridge` (task 3.11), recommended.** The decision
     r-bridge-before-ports already says "CI needs a job with R available for
     bridge tests", and without it the bridge code is never run in CI.
     Linux only, Python 3.13, `r-lib/actions/setup-r@v2` with `r-version:
     4.5.3` and P3M binaries, `setup-r-dependencies@v2` with
     `extra-packages: bioc::ALDEx2, bioc::maaslin3` and its cache, the apt
     headers the Dockerfile lists, then `uv run --group test --extra r pytest
     -m r`. Added to `check.needs`, with a `tests/test_ci.py` assertion in
     the existing style. Cost: a cold run compiles lme4-family packages
     (estimated 10-30 min, [UNVERIFIED]); cached runs a few minutes. Actions
     are pinned by SHA as the rest of `test.yaml`, looked up when the task is
     expanded. The alternative (no job; bridge tests run only locally) is
     cheaper and leaves the bridges untested in CI.

9. **`pl.consensus(table, *, top=30, ax=None) -> Axes`: a dot matrix.** Rows
   are features called by at least one method, sorted by `n_significant`
   then mean |effect|, at most `top`; columns are methods. A filled dot
   coloured by direction marks a call, a hollow dot a tested feature that was
   not called, nothing an untested one; consensus rows get a bold label.
   One Axes, as every `pl` function returns (R3.3). An UpSet plot needs two
   panels and `upsetplot` (not installed; would be a new dependency).
   `pl` and `da` are on the same layer and may not import each other, so
   `pl.consensus` reads the table's columns and imports nothing from `da`.

10. **Exit-gate notebook (`docs/tutorials/differential_abundance.md`).**
    GlobalPatterns glommed to genus, filtered once (`min_prevalence=0.2`),
    human samples (Feces, Skin, Tongue: 9) vs the 17 others (14 environmental, 3 mock communities): research C
    ran three scikit-bio methods on exactly this in 1.3 s with 5-43 hits each
    and visible disagreement. It is already a pooch dataset (no new
    download, hash or licence) and holds counts, which every method needs.
    The notebook says the contrast is a demonstration, not a biological
    claim. It runs ANCOM-BC2 and LinDA live: the docs build (CI docs job,
    Read the Docs) has no R, so the ALDEx2 and MaAsLin 3 calls are shown as
    a non-executed code block with their extra and the numbers the
    `r-bridge` job checks, pinned by a test, and the CI `r-bridge` job
    runs the same four-method consensus as a test. **(user)** An HMP2
    contrast was rejected: research C found 0-1 hits per method there, a
    poor demonstration.

11. **Dependencies.** Python: only the extra `r` (rpy2) in 3C. statsmodels,
    patsy, formulaic and joblib are all installed transitively through
    scikit-bio or scikit-learn, but biotapy imports none of them (notes 4, 6,
    7), so none is declared (R9.1 is about what biotapy declares; declaring
    an unused package would be R2.3's speculation). R packages live only in
    the golden image and the CI job.

12. **Roadmap corrections.**
    - 3.1: `pp.clr(adata, *, pseudocount=0.5)`; no `layer`.
    - 3.2: `pp.philr(tdata, *, pseudocount=0.5) -> TreeData`, `obsm["X_philr"]`.
    - 3.3 merges into 3.5 (the schema ships with its first consumer, so it
      is tested through a public function, R4.9).
    - 3.4: `da.ancombc2(adata, group, *, covariates=(), reference=None)`,
      golden vs `ANCOMBC::ancombc2`.
    - 3.5: `da.linda(adata, group, *, covariates=(), reference=None)`; the
      golden target is `MicrobiomeStat::linda` (no separate LinDA package is
      on CRAN or Bioconductor); no statsmodels.
    - 3.6: `da.aldex2(adata, group, *, mc_samples=128, reference=None, seed=None)`.
    - 3.7: `da.maaslin3(adata, group, *, covariates=(), reference=None, seed=None)`.
    - 3.8: `da.consensus(results, *, alpha=0.05, min_methods=2)`.
    - 3.9: `pl.consensus(table, *, top=30, ax=None)`.
    - New: 3.0 (CLR and PhILR goldens), 3.10b (exit-gate notebook), 3.12
      (Coming-from-R check), 3.13 (benchmarks), 3.14 (knowledge, at
      Checkpoint D), 3.15 (release 0.3.0).
    - Frontmatter `description` as above; `paths` gains `src/biotapy/pl/**`.

# Global constraints
- Python >= 3.12; no new runtime dependency in Phase 3; the extra `r`
  (rpy2) only after the user approves it at the start of task 3.6.
- `X` stays CSR. Densify once, only inside a wrapper whose library needs
  dense input (scikit-bio `clr`, `tree_basis` products, `ancombc2`; rpy2),
  with a comment and the memory cost in `Notes` (R6.2).
- Results go only to `layers["clr"]` and `obsm["X_philr"]`; `da` functions
  return tables and write nothing (R3.3, R6.3).
- Every `da` method returns `validate_result`'s schema: log2 `effect`, BH
  `qvalue`, one row per feature in `var_names` order, never a filtered
  subset.
- Calls are `qvalue < alpha`, strict.
- rpy2 and every R package are imported inside the function through
  `import_optional` (R4.6); `import biotapy` and the `import-without-extras`
  CI job stay green without them.
- R's random state is set only through `set.seed(<int>)` derived from
  `as_generator(seed)` right before an R call; never NumPy's global state.
- Golden files: `tests/r/export_golden.R`, run twice in the image built from
  `tests/r/Dockerfile`, bit-identical; each new R package is its own
  Dockerfile commit; files < 1 MB; never in CI.
- Size limits (R5) include ruff's `max-positional-args = 3` (`PLR0917`) and
  `max-args = 6` (`PLR0913`), which shaped `da.consensus`.
- mypy: scikit-bio's `tree_basis` is unannotated, so
  `skbio.stats.composition._base` joins `untyped_calls_exclude` (as
  `skbio.stats._subsample` did); rpy2 ships no `py.typed`
  ([UNVERIFIED]) and gets the same `follow_untyped_imports` override when 3C
  adds it. [Superseded by slice 3C decision 11: 3C adds no override.]
- Commits stage explicit paths only. Never stage `.claude/`, `.superpowers/`,
  `.worktrees/`, `notebooks/` or `build/`. Every task's last commit also
  stages `.knowledge/roadmap/phase-3-stats.md` with that task's boxes ticked
  and `.knowledge/log.md` with its dated line (R12.4).
- Every pytest, sphinx or Python run exports `BIOTAPY_DATA_DIR` to a
  scratch pooch cache; `~/.cache/biotapy` must not appear.

# Dependencies to approve (ask at the start of the task named)
| Task | Group | Package | Reason |
|---|---|---|---|
| 3.0 | R golden image only | Bioconductor `philr` (GPL-3) + apt `libuv1` | the `pp.philr` golden; run only in the image |
| 3.4 | R golden image only | Bioconductor `ANCOMBC` (Artistic-2.0) | the `da.ancombc2` golden |
| 3.5 | R golden image only | CRAN `MicrobiomeStat` (GPL-3) | the `da.linda` golden |
| 3.6 | extra `r` | `rpy2>=3.6.8` (GPL-2.0+) | the ALDEx2 and MaAsLin 3 bridges; the roadmap's planned extra |
| 3.6, 3.7 | R golden image and CI job | Bioconductor `ALDEx2` (GPL-3), `maaslin3` (MIT) | bridge goldens and the `r` tests |
| 3.11 | CI | `r-lib/actions/setup-r`, `setup-r-dependencies` (v2, pinned by SHA) | the `r-bridge` job |
| - | none | statsmodels, patsy, formulaic, joblib | not imported (design notes 4, 6, 7): nothing to declare |

# Review focus
The five ways real users are most likely to get a wrong answer from Phase 3
without an error. Each line names the test that pins it.

1. **A pseudocount on the wrong scale.** The default 0.5 on relative
   abundances or CPM swamps every rare feature. Expected: a `UserWarning`
   naming the smallest value, from `pp.clr` and `pp.philr`. Test: 3.1
   `test_clr_pseudocount_above_the_smallest_value_warns` (`pp.philr` shares
   the helper).
2. **An unrooted or multifurcating tree in PhILR.** Expected: `ValueError`
   naming the node, never an invented balance. Tests: 3.2
   `test_three_child_root_raises`, `test_internal_polytomy_raises`.
3. **A filtered TreeData.** Subsetting keeps one-child nodes. Expected: they
   are skipped with every child order kept, so signs equal philr's. Tests:
   3.2 `test_one_child_nodes_are_skipped_keeping_child_order` and
   `test_philr_matches_philr_philr` (its taxa are a subset).
4. **The reference level flips the sign.** patsy and R pick the
   alphabetically first level. Expected: `reference=` sets it, `contrast`
   says it, and a missing value raises instead of dropping rows. Tests (3B):
   3.5 `test_reference_sets_the_sign`, `test_missing_group_value_raises`;
   3.4 the same two.
5. **Methods that call "significant" differently.** Holm vs BH, `<=` vs `<`,
   per-method prevalence filters. Expected: BH everywhere, strict `<`, no
   filtering inside methods, and consensus refuses tables comparing
   different contrasts. Tests (3B): 3.8 `test_q_equal_to_alpha_is_not_called`,
   `test_results_with_different_contrasts_raise`; 3.4
   `test_qvalue_is_benjamini_hochberg`; 3.5 `test_no_feature_is_dropped`.

# Slices
| Slice | Delivers | Tasks | Ends with |
|---|---|---|---|
| **3A - Transforms** | `pp.clr` and `pp.philr`, equal to vegan and philr | 3.0 CLR and PhILR goldens · 3.1 `pp.clr` · 3.2 `pp.philr` | Checkpoint A |
| **3B - Native DA and consensus** | the schema, ANCOM-BC2, LinDA, the consensus table and plot, all without R | 3.5 `da.linda` with the schema (3.3) · 3.4 `da.ancombc2` · 3.8 `da.consensus` and the agreement decision · 3.9 `pl.consensus` | Checkpoint B |
| **3C - R bridges** | the extra `r`, ALDEx2, MaAsLin 3, R in CI | 3.6 `da.aldex2` (+ `r` extra, `da/_r.py`) · 3.7 `da.maaslin3` · 3.11 CI `r-bridge` job | Checkpoint C |
| **3D - Docs and release** | method pages, the DA guide, the exit-gate notebook, 0.3 | 3.10 docs · 3.10b tutorial · 3.12 Coming-from-R check · 3.13 benchmarks · Checkpoint D with 3.14 knowledge · 3.15 release 0.3.0 | exit gate |

Execution order inside 3A: **3.0 -> 3.1 -> 3.2 -> Checkpoint A.** The
goldens need no biotapy code and both golden tests read them; `pp.philr`
reuses `pp.clr`'s pseudocount helper.

Execution order inside 3B: **3.B0 -> 3.5 -> 3.4 -> 3.8 -> 3.9 -> Checkpoint B.**
Execution order inside 3C: **3.6 -> 3.7 -> 3.11 -> Checkpoint C.**
LinDA first: it is native end to end, so the schema's every column (with a
real `se`) is fixed by code biotapy owns before a wrapper maps a library's
columns onto it.

# Tasks (checklist)
- [x] 3.0 CLR and PhILR golden files (philr in the R image)
- [x] 3.1 `pp.clr(adata, *, pseudocount=0.5) -> AnnData`
- [x] 3.2 `pp.philr(tdata, *, pseudocount=0.5) -> TreeData`
- [x] Checkpoint A (PR #21 merged; the user approved slice 3A on 2026-10-05)
- [x] 3.B0 `fix(core)`: count tables hold non-negative whole numbers (approved with slice 3B)
- [x] 3.3 Result schema `da/_schema.py` (delivered inside 3.5)
- [x] 3.5 `da.linda(adata, group, *, covariates=(), reference=None) -> pd.DataFrame`
- [x] 3.4 `da.ancombc2(adata, group, *, covariates=(), reference=None) -> pd.DataFrame`
- [x] 3.8 `da.consensus(results, *, alpha=0.05, min_methods=2) -> pd.DataFrame` and the agreement decision
- [x] 3.9 `pl.consensus(table, *, top=30, ax=None) -> Axes`
- [x] Checkpoint B (PR #22 merged; the user approved slice 3B on 2026-10-05)
- [x] 3.6 `da.aldex2(adata, group, *, mc_samples=128, reference=None, seed=None) -> pd.DataFrame` and the extra `r`
- [x] 3.7 `da.maaslin3(adata, group, *, covariates=(), reference=None, seed=None) -> pd.DataFrame`
- [x] 3.11 CI job `r-bridge` for `-m r` tests
- [x] Checkpoint C (PR #23 merged; the user approved slice 3C on 2026-10-07)
- [x] 3.10 Method pages in `docs/methods/` and the DA guide
- [x] 3.10b `docs/tutorials/differential_abundance.md`, the exit-gate notebook
- [x] 3.12 Coming-from-R check
- [x] 3.13 asv benchmarks for `pp.philr`, `da.linda`, `da.ancombc2`
- [ ] Checkpoint D
- [ ] 3.14 Knowledge (the DA docs and benchmarks; da.md was updated for the bridges at Checkpoint C)
- [ ] 3.15 Release 0.3.0

# Exit gate
- [ ] Consensus report on one benchmark dataset (GlobalPatterns genus,
  human-associated vs the other samples: environmental and mock communities),
  executed notebook in docs (CI docs job).
- [x] Per-method agreement with the R reference, per
  [r-golden-parity](/contracts/r-golden-parity.md): `pp.clr`, `pp.philr`,
  `da.linda`, `da.ancombc2` in the network job; `da.aldex2`, `da.maaslin3` in
  the `r-bridge` job. Proven on PR #23 (run `37479113322`): `network`
  job `36 passed`, `r-bridge` job `24 passed`.
- [ ] All Phase 1 and 2 gates still green.

# Risks
- **rpy2 on Linux builds from source and lists Python <= 3.13** -> the CI job
  pins 3.13; the bridges are documented Linux/macOS-first; Windows is not
  tested.
- **R package drift.** The image is Bioconductor 3.22 (release is 3.23);
  maaslin3 1.4.0 changed its imports and may change numbers -> the CI job
  pins R 4.5.3 so Bioconductor resolves to 3.22 ([UNVERIFIED] that pak picks
  3.22 for R 4.5), and bridge goldens compare sign and rank, not digits.
- **scikit-bio's ANCOM-BC2 may not equal R's** (`conservative` variance,
  iteration control) -> measured in 3.4 before the tolerance is fixed; if
  only ranks agree, the docstring says so.
- **LinDA's mode estimator** (modeest `shorth` start, `bw.nrd0`) is ported
  from source never run here -> the `MicrobiomeStat::linda` golden is exact,
  so a mismatch fails loudly.
- **Consensus misread as a licence to shop for methods** -> the guide and
  `Notes` say to fix `methods` before looking; the two-step API keeps the
  list in the user's code.
- **GPL in the `r` extra** -> never bundled; the user rules on the stance
  (decision 14).

---
## Slice 3A - Transforms

**Goal:** a user turns counts into CLR values and PhILR balances that equal
`vegan::decostand(x, "clr")` / `mia::transformAssay` and `philr::philr`, and
is stopped, with the fix named, when the pseudocount or the tree cannot give
a meaningful answer.

### Slice 3A design
- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `tests/r/Dockerfile` | philr and `libuv1` in the golden image (3.0) |
  | `tests/r/export_golden.R` | the CLR and PhILR golden sections (3.0) |
  | `tests/golden/global_patterns/{clr,philr,philr_sbp}.csv.gz` | the goldens (3.0) |
  | `pp/_transform.py` | `clr` beside `relative`, and `pseudocounted`, the shared pseudocount check (3.1) |
  | `pp/_philr.py` | `philr` and `_binary_tree` (3.2) |
  | `_core/_tree.py` | `get_skbio_tree(adata, *, split_root=True)` (3.2) |

  `pseudocounted` is used by two topic files of `pp`; module-boundaries rule 1
  lets a topic file import a private helper from another topic file of the
  same subpackage (precedent: `tl/_permanova.py` imports `stored_distances`
  from `tl/_beta.py`), so it stays in `pp` (R4.3: one module uses it).
- **Facts the tasks rely on** (measured on the prototype; re-check each,
  R2.2).
  - `skbio.stats.composition.clr(mat, axis=-1, validate=True)` returns a dense
    `ndarray`; `tree_basis(tree)` returns `(coo_array (D-1, D), node names)`
    in level order, raises `ValueError: Not a strictly bifurcating tree.`, and
    an `IndexError` on a single tip (hence the two-feature check). On a
    two-tip tree `(a,b)r` the basis is `[[-0.7071, 0.7071]]`: the first child
    is negative.
  - `ndarray @ coo_array.T` returns an `ndarray`.
  - `TreeNode.prune()` collapses one-child nodes but appends the child at the
    end of the parent's children (its source has a TODO saying so): on
    `toy()` without `f2`, `n1`'s children became `(f3, f1)`.
  - `bt.datasets.toy()[:, :6]` is a TreeData whose tree is `root -> (n1, n2)`
    (n3 pruned with f7, f8); `toy()[:, ["f1", "f3", ...]]` keeps `n4` with one
    child.
  - GlobalPatterns' stored tree has a two-child root and no polytomy; its
    edge order follows ape's edge matrix, so `get_skbio_tree`'s children are
    in R's order (the golden's partitions match with equal numerators).
  - R: philr 1.36.0 (Bioconductor 3.22) has `pseudocount=`; `phylo2sbp` gives
    the sign partition; `makeNodeLabel` names nodes `n1...` (biotapy's are
    `n0...` and keep their pre-filter names, which is why the golden matches
    by partition, not name). philr's ggtree -> treeio -> fs chain installs fs
    as a P3M binary that needs `libuv.so.1` at load time.
  - vegan 2.7.3 `decostand(x, "clr", pseudocount = 0.5)` adds 0.5 to every
    entry.
  - Smallest non-zero relative abundance in `toy()` is 0.013 (1/77): the
    warning test pins that number.

### Slice 3A global constraints (in addition to the Phase 3 list)
- Golden tests carry `golden` + `network` (their input, GlobalPatterns, is
  downloaded by pooch) and run in the network CI job.
- Tests reach functions through `bt.pp.<fn>`. `_core` unit tests and test
  fixtures built with `_core` (`TreeData`, `tree_from_edges`,
  `get_skbio_tree`) are the exceptions, as in Phase 1 and 2.
- Run commands: `uv run --group test pytest <path> -q`; golden tests with
  `-m "golden or network"`. Gate before every commit: `uvx prek run
  --all-files` (new files must be staged first: prek's `--all-files` only
  sees tracked files); with docs changes also `uv run --group doc
  sphinx-build -W -b html docs docs/_build/html`.

### Slice 3A review focus
The Phase 3 review focus items 1-3 are slice 3A's. In addition:
- **Purity:** both functions copy; `X`, `layers`, `vart` of the input are
  untouched. Tests: 3.1 `test_clr_keeps_x_and_input`, 3.2
  `test_keeps_x_tree_and_input`.
- **Saving:** the `obsm` DataFrame round-trips through h5td; a later feature
  filter drops it. Tests: 3.2 `test_round_trips_through_h5td`,
  `test_feature_filter_drops_the_balances`.

---

### Task 3.0: CLR and PhILR golden files

**Files:** modify `tests/r/Dockerfile`, `tests/r/export_golden.R`,
`tests/golden/VERSIONS.txt`, `.knowledge/contracts/r-golden-parity.md`;
create `tests/golden/global_patterns/clr.csv.gz`,
`tests/golden/global_patterns/philr.csv.gz`,
`tests/golden/global_patterns/philr_sbp.csv.gz`.
**Interfaces (produces):**
- `clr.csv.gz`: `sample_id, taxon_id, value` for every 100th GlobalPatterns
  taxon (193 taxa x 26 samples), CLR of all 19,216 taxa with pseudocount 0.5.
- `philr.csv.gz`: `sample_id, balance, value` (292 balances x 26 samples).
- `philr_sbp.csv.gz`: `balance, taxon_id, sign` (+1 numerator, -1
  denominator; 4,256 rows), R's balance names `n1...`.

- [x] **Step 1: Ask** the user to approve philr (GPL-3) and `libuv1` in the
  golden image (dependency table; they live only in the image).
- [x] **Step 2: Add philr to the image.** In `tests/r/Dockerfile`:
  ```diff
  @@ -14,6 +14,10 @@ RUN Rscript -e 'install.packages("BiocManager")' \
   RUN Rscript -e 'BiocManager::install("phyloseq", version = "3.22", ask = FALSE, update = FALSE)'
   # picante (CRAN, same P3M snapshot): the Faith PD golden file. vegan and ape come with phyloseq.
   RUN Rscript -e 'install.packages("picante")'
  -RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE))'
  +# philr (Bioconductor, GPL-3): the pp.philr golden file; biotapy never calls it. Its ggtree -> treeio -> fs
  +# chain installs fs as a P3M binary that links libuv at load time.
  +RUN apt-get update && apt-get install -y --no-install-recommends libuv1 && rm -rf /var/lib/apt/lists/*
  +RUN Rscript -e 'BiocManager::install("philr", version = "3.22", ask = FALSE, update = FALSE)'
  +RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE))'
   WORKDIR /work
   CMD ["Rscript", "tests/r/export_golden.R"]
  ```
  In `.knowledge/contracts/r-golden-parity.md` statement 1, replace the
  three lines from "`Biostrings`, `vegan` and `ape`, plus CRAN `picante`"
  to "(rules.md R2.3)." with:
  ```text
   `Biostrings`, `vegan` and `ape`, plus CRAN `picante` for Faith PD and
   Bioconductor `philr` for `pp.philr`, with the `libuv1` runtime library its
   `fs` binary loads); a new golden function that needs another package adds
   it in its own commit (rules.md R2.3).
  ```
  Add under a new heading `## <date of the commit> (Phase 3, slice 3A)` at
  the top of `.knowledge/log.md`:
  ```text
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `philr` (and `libuv1`) for the `pp.philr` golden file.
  ```
  Build: `docker build -t biotapy-golden tests/r`. Expected: the build ends
  with the `stopifnot(requireNamespace(...))` layer passing;
  `docker run --rm biotapy-golden Rscript -e 'packageVersion("philr")'` prints
  `1.36.0`, and ape 5.8.1, vegan 2.7.3, phyloseq 1.54.2 are unchanged
  (`update = FALSE`). Without the `libuv1` line the guard fails with
  `requireNamespace("philr", quietly = TRUE) is not TRUE` after
  `libuv.so.1: cannot open shared object file`.
- [x] **Step 3: Commit the image change on its own** (contract statement 1):
  ```bash
  git add tests/r/Dockerfile .knowledge/contracts/r-golden-parity.md .knowledge/log.md
  git commit -m "build(r): add philr to the golden image"
  ```
- [x] **Step 4: Export.** In `tests/r/export_golden.R`, insert before
  `## Synthetic phyloseq fixtures: biotapy's toy() numbers, no third-party data`:
  ```r
  ## Slice 3A golden files: CLR and PhILR (GlobalPatterns)
  # mia::transformAssay(method = "clr", pseudocount = 0.5) delegates to this call, which adds 0.5 to every count.
  clr <- vegan::decostand(samples_as_rows(GlobalPatterns), "clr", pseudocount = 0.5)
  # Every 100th taxon keeps the file small; each value still depends on all 19,216 taxa.
  every_100th <- seq(1, ncol(clr), by = 100)
  write_golden(
    data.frame(
      sample_id = rep(rownames(clr), times = length(every_100th)),
      taxon_id = rep(colnames(clr)[every_100th], each = nrow(clr)),
      value = as.vector(clr[, every_100th])
    ),
    file.path(gp, "clr.csv.gz")
  )
  # The 293 taxa with more than 3 reads in over half of the samples. prune_taxa (ape::drop.tip) leaves a rooted
  # binary tree; makeNodeLabel names the internal nodes, which philr uses as balance names.
  gp_philr <- filter_taxa(GlobalPatterns, function(x) sum(x > 3) > 0.5 * length(x), TRUE)
  philr_tree <- ape::makeNodeLabel(phy_tree(gp_philr), method = "number", prefix = "n")
  stopifnot(ape::is.rooted(philr_tree), ape::is.binary(philr_tree))
  balances <- suppressMessages(philr::philr(samples_as_rows(gp_philr), philr_tree, pseudocount = 0.5))
  write_golden(
    data.frame(
      sample_id = rep(rownames(balances), times = ncol(balances)),
      balance = rep(colnames(balances), each = nrow(balances)),
      value = as.vector(balances)
    ),
    file.path(gp, "philr.csv.gz")
  )
  # Each balance's sequential binary partition: +1 for the taxa in its numerator, -1 for its denominator.
  sbp <- philr::phylo2sbp(philr_tree)
  signs <- which(sbp != 0, arr.ind = TRUE)
  write_golden(
    data.frame(balance = colnames(sbp)[signs[, "col"]], taxon_id = rownames(sbp)[signs[, "row"]], sign = sbp[signs]),
    file.path(gp, "philr_sbp.csv.gz")
  )
  ```
  and add `philr` to the `VERSIONS.txt` lines at the end:
  ```r
  writeLines(c(
    R.version.string,
    paste("Bioconductor", as.character(BiocManager::version())),
    paste0("phyloseq ", packageVersion("phyloseq")),
    paste0("vegan ", packageVersion("vegan")),
    paste0("ape ", packageVersion("ape")),
    paste0("picante ", packageVersion("picante")),
    paste0("philr ", packageVersion("philr"))
  ), "tests/golden/VERSIONS.txt")
  ```
- [x] **Step 5: Run twice, check bit identity** (playbook
  regenerate-golden-files):
  ```bash
  mkdir -p build
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
  sha256sum tests/golden/*/*.csv.gz tests/data/phyloseq/* tests/data/dada2/* > build/golden-run1.sha
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
  sha256sum -c build/golden-run1.sha
  git status --short
  ```
  Expected: every line `OK`; `git status` shows only `M tests/golden/VERSIONS.txt`
  (one new line, `philr 1.36.0`) and the three new files, about 28 KB
  (`clr`), 85 KB (`philr`) and 15 KB (`philr_sbp`). Every existing golden
  file is byte-identical. The run prints `Found more than one class "phylo"
  in cache` (phyloseq and tidytree both define it); harmless.
- [x] **Step 6: Gate and commit.** `uv run --group test pytest
  tests/test_data_files.py -q` -> passes (each file < 1 MB, VERSIONS names
  R and Bioconductor). Then `uvx prek run --all-files` and:
  ```bash
  git add tests/r/export_golden.R tests/golden/VERSIONS.txt tests/golden/global_patterns/clr.csv.gz \
    tests/golden/global_patterns/philr.csv.gz tests/golden/global_patterns/philr_sbp.csv.gz \
    .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  git commit -m "test(golden): export CLR and PhILR golden files"
  ```
  (The log gets this task's line only if a concept changed; the roadmap tick
  is the concept change here: add `- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.0 done: CLR and PhILR golden files from vegan 2.7.3 and philr 1.36.0.`)

### Task 3.1: `pp.clr`

**Files:** modify `src/biotapy/pp/_transform.py`, `src/biotapy/pp/__init__.py`,
`tests/pp/test_transform.py`, `tests/pp/test_transform_golden.py`,
`docs/guide/transforms.md`, `docs/api.md`, `docs/contributing.md`,
`.knowledge/contracts/data-model-slots.md`,
`.knowledge/contracts/r-golden-parity.md`.
**Interfaces:**
- Consumes: `tests/golden/global_patterns/clr.csv.gz` (3.0); `_core`
  `add_provenance`, `as_csr`, `warn_user`.
- Produces: `bt.pp.clr(adata: AnnData, *, pseudocount: float = 0.5) -> AnnData`
  (a copy with `layers["clr"]`, dense float64); private
  `pp._transform.pseudocounted(adata: AnnData, pseudocount: float, *, func: str,
  columns: npt.NDArray[np.intp] | None = None) -> npt.NDArray[np.float64]` (dense
  `X + pseudocount`, its `columns` in that order if given; raises on a bad
  pseudocount, a negative or non-finite `X`, or a remaining zero, naming
  `func`; warns when `pseudocount` exceeds the smallest non-zero value), used
  by 3.2.

- [x] **Step 1: Failing tests.** In `tests/pp/test_transform.py` add
  `import pytest` after `import pandas as pd`, and append:
  ```python
  def _clr_by_hand(dense: np.ndarray, pseudocount: float) -> np.ndarray:
      logs = np.log(dense + pseudocount)
      return logs - logs.mean(axis=1, keepdims=True)


  def test_clr_is_the_log_minus_the_mean_log(make_adata):
      dense = np.array([[1.0, 3.0, 0.0], [0.0, 2.0, 6.0]])
      np.testing.assert_allclose(bt.pp.clr(make_adata(dense)).layers["clr"], _clr_by_hand(dense, 0.5), rtol=1e-12)


  def test_clr_is_dense_float64_and_rows_sum_to_zero():
      out = bt.pp.clr(bt.datasets.toy())
      assert isinstance(out.layers["clr"], np.ndarray) and out.layers["clr"].dtype == np.float64
      np.testing.assert_allclose(out.layers["clr"].sum(axis=1), 0.0, atol=1e-12)


  def test_clr_keeps_x_and_input(assert_unchanged):
      tdata = bt.datasets.toy()
      before = tdata.copy()
      out = bt.pp.clr(tdata)
      assert_unchanged(before, tdata)
      assert (out.X != tdata.X).nnz == 0 and type(out) is type(tdata)


  def test_clr_all_zero_sample_is_all_zero(make_adata):
      out = bt.pp.clr(make_adata(np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])))
      np.testing.assert_array_equal(out.layers["clr"][0], 0.0)


  def test_clr_all_zero_feature_is_finite(make_adata):
      out = bt.pp.clr(make_adata(np.array([[0.0, 4.0, 1.0], [0.0, 2.0, 3.0]])))
      assert np.isfinite(out.layers["clr"]).all()
      np.testing.assert_allclose(out.layers["clr"].sum(axis=1), 0.0, atol=1e-12)


  def test_clr_single_sample():
      out = bt.pp.clr(bt.datasets.toy()[:1].copy())
      assert out.layers["clr"].shape == (1, 8)


  def test_clr_without_pseudocount_needs_no_zeros(make_adata):
      dense = np.array([[1.0, 2.0], [3.0, 5.0]])
      np.testing.assert_allclose(
          bt.pp.clr(make_adata(dense), pseudocount=0).layers["clr"], _clr_by_hand(dense, 0.0), rtol=1e-12
      )


  def test_clr_zero_without_pseudocount_raises(make_adata):
      with pytest.raises(ValueError, match="pass pseudocount > 0 to pp.clr"):
          bt.pp.clr(make_adata(np.array([[0.0, 2.0], [3.0, 5.0]])), pseudocount=0)


  @pytest.mark.parametrize("pseudocount", [-0.5, np.nan, np.inf])
  def test_clr_bad_pseudocount_raises(pseudocount):
      with pytest.raises(ValueError, match="pseudocount must be a finite number >= 0"):
          bt.pp.clr(bt.datasets.toy(), pseudocount=pseudocount)


  @pytest.mark.parametrize("bad", [-1.0, np.nan])
  def test_clr_negative_or_missing_value_raises(make_adata, bad):
      with pytest.raises(ValueError, match="pp.clr needs finite, non-negative values in X"):
          bt.pp.clr(make_adata(np.array([[bad, 2.0], [3.0, 5.0]])))


  def test_clr_pseudocount_above_the_smallest_value_warns():
      relative = bt.pp.relative(bt.datasets.toy())
      relative.X = relative.layers["relative"]
      with pytest.warns(UserWarning, match=r"pseudocount=0.5 is larger than the smallest non-zero value in X \(0.013\)"):
          bt.pp.clr(relative)


  def test_clr_on_counts_does_not_warn(recwarn):
      bt.pp.clr(bt.datasets.toy())
      assert not [w for w in recwarn if issubclass(w.category, UserWarning)]


  def test_clr_records_provenance():
      entries = bt.pp.clr(bt.datasets.toy(), pseudocount=1).uns["biotapy"]["provenance"]
      assert json.loads(entries[-1]) == {"step": "pp.clr", "version": bt.__version__, "params": {"pseudocount": 1}}


  @given(
      arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 1000)),
      # At most 1, the smallest non-zero count, so no call warns.
      st.floats(0.01, 1),
      st.floats(0.01, 100),
  )
  def test_clr_rows_sum_to_zero_and_ignore_scale(dense, pseudocount, scale):
      adata = ad.AnnData(
          X=sp.csr_matrix(dense),
          obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
          var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
      )
      out = bt.pp.clr(adata, pseudocount=pseudocount).layers["clr"]
      np.testing.assert_allclose(out.sum(axis=1), 0.0, atol=1e-9)
      # CLR is scale invariant: scaling X and the pseudocount together changes nothing.
      scaled = adata.copy()
      scaled.X = sp.csr_matrix(dense * scale)
      np.testing.assert_allclose(
          bt.pp.clr(scaled, pseudocount=pseudocount * scale).layers["clr"], out, rtol=1e-9, atol=1e-9
      )


  @pytest.mark.parametrize("pseudocount", [True, False, "0.5", None])
  def test_clr_non_numeric_pseudocount_raises(pseudocount):
      with pytest.raises(TypeError, match="pseudocount must be a real number"):
          bt.pp.clr(bt.datasets.toy(), pseudocount=pseudocount)
  ```
  Append to `tests/pp/test_transform_golden.py`:
  ```python
  def test_clr_matches_vegan_decostand():
      # R: vegan::decostand(x, "clr", pseudocount = 0.5) on all of GlobalPatterns; every 100th taxon is kept.
      golden = pd.read_csv(GOLDEN / "clr.csv.gz", dtype={"sample_id": str, "taxon_id": str})
      out = bt.pp.clr(bt.datasets.global_patterns())
      rows, cols = out.obs_names.get_indexer(golden["sample_id"]), out.var_names.get_indexer(golden["taxon_id"])
      assert (rows >= 0).all() and (cols >= 0).all()
      np.testing.assert_allclose(out.layers["clr"][rows, cols], golden["value"], rtol=1e-7)
  ```
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/pp/test_transform.py -q` -> `17 failed, 9 passed`
  (`AttributeError: module 'biotapy.pp' has no attribute 'clr'`);
  `uv run --group test pytest tests/pp/test_transform_golden.py -q -m "golden or network"`
  -> `1 failed, 1 passed`.
- [x] **Step 3: Implement.** `src/biotapy/pp/_transform.py` becomes
  (`relative` unchanged):
  ```python
  """Per-sample transforms: add one layer, keep everything else."""

  import numpy as np
  import numpy.typing as npt
  from anndata import AnnData
  from skbio.stats.composition import clr as skbio_clr

  from biotapy._core import add_provenance, as_csr, divide_rows, warn_user


  def relative(adata: AnnData) -> AnnData:
      """Add per-sample relative abundance as ``layers['relative']``.

      Parameters
      ----------
      adata
          Samples x features; ``X`` holds counts or another non-negative abundance.

      Returns
      -------
      AnnData
          A copy of ``adata`` (a TreeData stays a TreeData) with
          ``layers['relative']``; ``X`` is unchanged.

      Notes
      -----
      R equivalent: ``phyloseq::transform_sample_counts``, ``mia::transformAssay``
      Guide: :doc:`/guide/transforms`

      All-zero samples stay all-zero, where phyloseq returns ``NaN``.

      Examples
      --------
      >>> import biotapy as bt
      >>> out = bt.pp.relative(bt.datasets.toy())
      >>> round(float(out.layers["relative"][0].sum()), 6)
      1.0
      """
      X = as_csr(adata.X)
      sums = np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
      out = adata.copy()
      out.layers["relative"] = divide_rows(X, sums)
      add_provenance(out, "pp.relative")
      return out


  def clr(adata: AnnData, *, pseudocount: float = 0.5) -> AnnData:
      """Add the centred log-ratio transform of each sample as ``layers['clr']``.

      Parameters
      ----------
      adata
          Samples x features; ``X`` holds counts or another non-negative abundance.
      pseudocount
          Added to every value of ``X`` before the logarithm, so zeros have one.

      Returns
      -------
      AnnData
          A copy of ``adata`` (a TreeData stays a TreeData) with ``layers['clr']``,
          a dense float64 array in which every sample sums to 0; ``X`` is unchanged.

      Raises
      ------
      TypeError
          ``pseudocount`` is a bool or not a real number.
      ValueError
          ``pseudocount`` is negative or not finite; ``X`` holds a negative or
          non-finite value; or ``X`` plus ``pseudocount`` still holds a zero.

      Warns
      -----
      UserWarning
          ``pseudocount`` is larger than the smallest non-zero value in ``X``, as when
          the default 0.5 meets relative abundances.

      Notes
      -----
      R equivalent: ``mia::transformAssay``, ``vegan::decostand``
      Guide: :doc:`/guide/transforms`

      Equals ``mia::transformAssay(tse, method = "clr", pseudocount = 0.5)`` and
      ``vegan::decostand(x, "clr", pseudocount = 0.5)``: the pseudocount is added to
      every value, zeros or not, and the logarithm is natural. mia's default
      ``pseudocount = FALSE`` fails on zeros; biotapy defaults to 0.5, the value
      LinDA uses. The transform is scale invariant, so counts and their relative
      abundances give the same result only when the pseudocount is scaled with them.
      An all-zero sample gives an all-zero CLR row.

      CLR has no zeros, so ``X`` is densified once and ``layers['clr']`` is dense:
      8 bytes x samples x features. Peak memory is three to five such arrays (the
      dense copy of ``X``, the transform's temporaries and its output; 3.1x to 4.8x
      measured on 400 x 512 and 2,000 x 2,000 tables, the most when ``X`` is mostly
      non-zero), so budget for that on large tables.

      References
      ----------
      Aitchison J (1986) The Statistical Analysis of Compositional Data. Chapman & Hall.

      Examples
      --------
      >>> import biotapy as bt
      >>> out = bt.pp.clr(bt.datasets.toy())
      >>> round(float(out.layers["clr"][0, 0]), 3)
      1.048
      """
      values = pseudocounted(adata, pseudocount, func="pp.clr")
      out = adata.copy()
      out.layers["clr"] = skbio_clr(values)
      add_provenance(out, "pp.clr", pseudocount=pseudocount)
      return out


  def pseudocounted(
      adata: AnnData, pseudocount: float, *, func: str, columns: npt.NDArray[np.intp] | None = None
  ) -> npt.NDArray[np.float64]:
      """``X`` (its ``columns``, in that order, if given) as a dense float64 array plus ``pseudocount``, checked > 0."""
      if isinstance(pseudocount, bool) or not isinstance(pseudocount, int | float | np.integer | np.floating):
          msg = f"pseudocount must be a real number, got {pseudocount!r}"
          raise TypeError(msg)
      if not np.isfinite(pseudocount) or pseudocount < 0:
          msg = f"pseudocount must be a finite number >= 0, got {pseudocount!r}"
          raise ValueError(msg)
      X = as_csr(adata.X).astype(np.float64)
      if columns is not None:
          # Reordered while sparse, so the dense copy below is the only one.
          X = X[:, columns]
      if not np.all(np.isfinite(X.data)) or np.any(X.data < 0):
          msg = f"{func} needs finite, non-negative values in X"
          raise ValueError(msg)
      positive = X.data[X.data > 0]
      if positive.size and pseudocount > positive.min():
          warn_user(
              f"pseudocount={pseudocount} is larger than the smallest non-zero value in X ({positive.min():.3g}), "
              "so it swamps the rarest features; for relative abundances pass a pseudocount on their scale"
          )
      # scikit-bio's log-ratio functions need dense input (rules.md R6.2): one dense copy of X.
      values = X.toarray()
      values += pseudocount
      if np.any(values <= 0):
          msg = f"X holds zeros, whose logarithm is undefined; pass pseudocount > 0 to {func}"
          raise ValueError(msg)
      return values
  ```
  `src/biotapy/pp/__init__.py`:
  ```python
  from ._filter import filter_features, filter_samples
  from ._glom import tax_glom
  from ._rarefy import rarefy
  from ._transform import clr, relative

  __all__ = ["clr", "filter_features", "filter_samples", "rarefy", "relative", "tax_glom"]
  ```
- [x] **Step 4: Run, expect pass** - the same two commands -> `26 passed`;
  `2 passed`. The property test also under `--hypothesis-seed=1`, `2`, `3`.
- [x] **Step 5: Docs.** Append to `docs/guide/transforms.md`:
  ```markdown
  ## Centred log-ratio (CLR)

  Microbiome counts are compositional: a sequencing run fixes the total, so only ratios between
  features carry information. `bt.pp.clr` takes the logarithm of each value relative to the
  geometric mean of its sample, and stores it in `layers["clr"]`:

  ```python
  out = bt.pp.clr(tdata)  # pseudocount=0.5
  out.layers["clr"][0].sum()  # 0.0, up to rounding
  ```

  The logarithm of zero is undefined, so a pseudocount is added to every value first, zeros or
  not, as `mia::transformAssay(method = "clr", pseudocount = 0.5)` and
  `vegan::decostand(x, "clr", pseudocount = 0.5)` do. The default 0.5 suits counts. On relative
  abundances or CPM it would swamp the rarest features, so biotapy warns when the pseudocount is
  larger than the smallest non-zero value; pass one on the data's scale, for example half that value.

  An all-zero sample gives an all-zero CLR row.

  `layers["clr"]` is dense: CLR has no zeros, so it takes 8 bytes per sample and feature, and the call peaks at three to five such arrays, the most when `X` is mostly non-zero.
  ```
  In `docs/api.md` add `pp.clr` first in the Preprocessing autosummary. In
  `docs/contributing.md` the network-test sentence reads "compare biotapy with
  R on that data: `pp.relative`, `pp.clr`, `pp.tax_glom`, ...".
- [x] **Step 6: Contracts.** In `.knowledge/contracts/data-model-slots.md`:
  - the `layers` row's keys become ``relative` (sparse CSR); `clr` from `pp.clr` (dense float64 `ndarray`: CLR has no zeros)``;
  - the Propagation row `Layer-adding (`pp.relative`; `pp.clr` in Phase 3)` becomes `Layer-adding (`pp.relative`, `pp.clr`)`.

  In `.knowledge/contracts/r-golden-parity.md` statement 4, the first row
  reads "Deterministic numeric (glom sums, relative, CLR, alpha, Bray-Curtis,
  UniFrac)". Log line:
  ```text
  - **Update**: [data-model-slots](contracts/data-model-slots.md): `layers["clr"]` from `pp.clr` is a dense float64 array; `pp.clr` joins the layer-adding row. [r-golden-parity](contracts/r-golden-parity.md): CLR is compared elementwise.
  ```
- [x] **Step 7: Gate and commit**
  ```bash
  uvx prek run --all-files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git add src/biotapy/pp/_transform.py src/biotapy/pp/__init__.py tests/pp/test_transform.py \
    tests/pp/test_transform_golden.py docs/guide/transforms.md docs/api.md docs/contributing.md \
    .knowledge/contracts/data-model-slots.md .knowledge/contracts/r-golden-parity.md \
    .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  git commit -m "feat(pp): add the centred log-ratio transform"
  ```

### Task 3.2: `pp.philr`

**Files:** create `src/biotapy/pp/_philr.py`, `tests/pp/test_philr.py`,
`tests/pp/test_philr_golden.py`; modify `src/biotapy/_core/_tree.py`,
`src/biotapy/pp/__init__.py`, `tests/core/test_tree.py`, `pyproject.toml`,
`docs/guide/transforms.md`, `docs/api.md`, `docs/contributing.md`,
`.knowledge/contracts/data-model-slots.md`,
`.knowledge/contracts/r-golden-parity.md`,
`.knowledge/contracts/tree-access.md`.
**Interfaces:**
- Consumes: `pseudocounted` (3.1); `philr.csv.gz`, `philr_sbp.csv.gz` (3.0).
- Produces: `bt.pp.philr(tdata: TreeData, *, pseudocount: float = 0.5) -> TreeData`
  (a copy with `obsm["X_philr"]`: `pd.DataFrame`, index `obs_names`, one
  column per internal node with two children in preorder); `_core.get_skbio_tree(adata:
  AnnData, *, split_root: bool = True) -> TreeNode` (default unchanged).

- [x] **Step 1: Failing tests.** `tests/pp/test_philr.py`:
  ```python
  import json

  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from hypothesis import given, settings
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays

  import biotapy as bt
  from biotapy._core import TreeData, get_tree, tree_from_edges


  def _toy6() -> TreeData:
      # toy()'s root has three children; its first six features sit under a binary root.
      return bt.datasets.toy()[:, :6].copy()


  def _gmean(values: np.ndarray) -> float:
      return float(np.exp(np.log(values).mean()))


  def test_balances_are_scaled_log_ratios_of_geometric_means():
      tdata = _toy6()
      x = tdata.X.toarray()[0] + 0.5
      balances = bt.pp.philr(tdata).obsm["X_philr"].loc["s1"]
      # n4 = (f1, f2); n1 = (n4, f3); root = (n1, n2): the first child is the numerator, as in philr::philr.
      assert balances["n4"] == pytest.approx(np.sqrt(1 / 2) * np.log(x[0] / x[1]))
      assert balances["n1"] == pytest.approx(np.sqrt(2 * 1 / 3) * np.log(_gmean(x[:2]) / x[2]))
      assert balances["root"] == pytest.approx(np.sqrt(3 * 3 / 6) * np.log(_gmean(x[:3]) / _gmean(x[3:])))


  def test_one_column_per_internal_node_in_preorder():
      out = bt.pp.philr(_toy6()).obsm["X_philr"]
      assert isinstance(out, pd.DataFrame) and out.index.tolist() == [f"s{i}" for i in range(1, 7)]
      assert out.columns.tolist() == ["root", "n1", "n4", "n2", "n5"]


  def test_keeps_x_tree_and_input(assert_unchanged):
      tdata = _toy6()
      before = tdata.copy()
      edges = [(u, v, dict(data)) for u, v, data in get_tree(tdata).edges(data=True)]
      out = bt.pp.philr(tdata)
      assert_unchanged(before, tdata)
      assert list(get_tree(tdata).edges(data=True)) == edges
      assert isinstance(out, TreeData) and "phylo" in out.vart and (out.X != tdata.X).nnz == 0


  def test_one_child_nodes_are_skipped_keeping_child_order():
      # Without f2, n4 has one child: n1 becomes (f1, f3) with f1 still first, so still the numerator.
      tdata = _toy6()[:, ["f1", "f3", "f4", "f5", "f6"]].copy()
      x = tdata.X.toarray()[0] + 0.5
      out = bt.pp.philr(tdata).obsm["X_philr"]
      assert out.columns.tolist() == ["root", "n1", "n2", "n5"]
      assert out.loc["s1", "n1"] == pytest.approx(np.sqrt(1 / 2) * np.log(x[0] / x[1]))


  def test_three_child_root_raises():
      with pytest.raises(ValueError, match=r"the node\(s\) \['root'\] have more than two children"):
          bt.pp.philr(bt.datasets.toy())


  def test_internal_polytomy_raises():
      edges = [("r", "a", 1.0), ("r", "p", 1.0), ("p", "b", 1.0), ("p", "c", 1.0), ("p", "d", 1.0)]
      tdata = TreeData(
          X=sp.csr_matrix(np.ones((2, 4))),
          obs=pd.DataFrame(index=["s1", "s2"]),
          var=pd.DataFrame(index=list("abcd")),
          vart={"phylo": tree_from_edges(edges)},
          label=None,
      )
      with pytest.raises(ValueError, match=r"\['p'\] have more than two children"):
          bt.pp.philr(tdata)


  def test_feature_outside_the_tree_raises():
      tdata = TreeData(
          X=sp.csr_matrix(np.ones((2, 3))),
          obs=pd.DataFrame(index=["s1", "s2"]),
          var=pd.DataFrame(index=["a", "b", "extra"]),
          vart={"phylo": tree_from_edges([("r", "a", 1.0), ("r", "b", 1.0)])},
          label=None,
      )
      with pytest.raises(ValueError, match=r"every feature to be a tip of the tree; 1 feature\(s\) are not"):
          bt.pp.philr(tdata)


  def test_fewer_than_two_features_raise():
      with pytest.raises(ValueError, match="at least two features, got 1"):
          bt.pp.philr(_toy6()[:, :1].copy())


  def test_plain_anndata_raises():
      with pytest.raises(TypeError, match="needs a TreeData"):
          bt.pp.philr(_toy6().to_adata())


  def test_all_zero_sample_has_zero_balances():
      tdata = _toy6()
      dense = tdata.X.toarray()
      dense[0] = 0
      tdata.X = sp.csr_matrix(dense)
      np.testing.assert_array_equal(bt.pp.philr(tdata).obsm["X_philr"].loc["s1"], 0.0)


  def test_all_zero_feature_is_finite():
      tdata = _toy6()
      dense = tdata.X.toarray()
      dense[:, 0] = 0
      tdata.X = sp.csr_matrix(dense)
      assert np.isfinite(bt.pp.philr(tdata).obsm["X_philr"].to_numpy()).all()


  def test_single_sample():
      assert bt.pp.philr(_toy6()[:1].copy()).obsm["X_philr"].shape == (1, 5)


  def test_philr_pseudocount_above_the_smallest_value_warns():
      relative = _toy6()
      relative.X = bt.pp.relative(relative).layers["relative"]
      with pytest.warns(UserWarning, match=r"pseudocount=0.5 is larger than the smallest non-zero value in X \(0.0"):
          bt.pp.philr(relative)


  def test_zero_without_pseudocount_raises():
      with pytest.raises(ValueError, match="pass pseudocount > 0 to pp.philr"):
          bt.pp.philr(_toy6(), pseudocount=0)


  @pytest.mark.parametrize("pseudocount", [True, "0.5", None])
  def test_non_numeric_pseudocount_raises(pseudocount):
      with pytest.raises(TypeError, match="pseudocount must be a real number"):
          bt.pp.philr(_toy6(), pseudocount=pseudocount)


  def test_records_provenance():
      entries = bt.pp.philr(_toy6()).uns["biotapy"]["provenance"]
      assert json.loads(entries[-1])["step"] == "pp.philr" and json.loads(entries[-1])["params"] == {"pseudocount": 0.5}


  def test_round_trips_through_h5td(tmp_path):
      import treedata

      out = bt.pp.philr(_toy6())
      out.write_h5td(tmp_path / "philr.h5td")
      back = treedata.read_h5td(tmp_path / "philr.h5td")
      pd.testing.assert_frame_equal(back.obsm["X_philr"], out.obsm["X_philr"])


  def test_feature_filter_drops_the_balances():
      out = bt.pp.filter_features(bt.pp.philr(_toy6()), min_prevalence=0.5)
      assert "X_philr" not in out.obsm


  @settings(deadline=None)
  @given(arrays(np.int64, st.tuples(st.integers(1, 5), st.just(6)), elements=st.integers(0, 1000)), st.floats(0.01, 1))
  def test_balances_keep_the_clr_distance(dense, pseudocount):
      # An isometric log-ratio basis: each sample's balances have the length of its CLR vector.
      tdata = _toy6()[: dense.shape[0]].copy()
      tdata.X = sp.csr_matrix(dense)
      balances = bt.pp.philr(tdata, pseudocount=pseudocount).obsm["X_philr"].to_numpy()
      clr = bt.pp.clr(tdata, pseudocount=pseudocount).layers["clr"]
      np.testing.assert_allclose(np.linalg.norm(balances, axis=1), np.linalg.norm(clr, axis=1), rtol=1e-9, atol=1e-9)


  def test_wide_node_error_counts_every_wide_node():
      # Five three-child nodes (q, p1-p4) under a binary root: the message lists the first three, in postorder, and counts all five.
      edges = [("r", "q", 1.0), ("r", "p4", 1.0), ("q", "p1", 1.0), ("q", "p2", 1.0), ("q", "p3", 1.0)]
      edges += [(f"p{i + 1}", name, 1.0) for i in range(4) for name in (f"x{i}", f"y{i}", f"z{i}")]
      tips = [f"{letter}{i}" for i in range(4) for letter in "xyz"]
      tdata = TreeData(
          X=sp.csr_matrix(np.ones((2, len(tips)))),
          obs=pd.DataFrame(index=["s1", "s2"]),
          var=pd.DataFrame(index=tips),
          vart={"phylo": tree_from_edges(edges)},
          label=None,
      )
      with pytest.raises(ValueError, match=r"\['p1', 'p2', 'p3'\] have more than two children \(5 in total\)"):
          bt.pp.philr(tdata)


  def test_feature_order_does_not_change_the_balances():
      tdata = _toy6()
      shuffled = tdata[:, ["f4", "f1", "f6", "f3", "f5", "f2"]].copy()
      pd.testing.assert_frame_equal(bt.pp.philr(shuffled).obsm["X_philr"], bt.pp.philr(tdata).obsm["X_philr"])
  ```
  `tests/pp/test_philr_golden.py`:
  ```python
  from pathlib import Path

  import numpy as np
  import pandas as pd
  import pytest

  import biotapy as bt
  from biotapy._core import get_skbio_tree

  GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
  pytestmark = [pytest.mark.golden, pytest.mark.network]


  def test_philr_matches_philr_philr():
      # R: philr::philr(x, tree, pseudocount = 0.5) on the 293 GlobalPatterns taxa with more than 3 reads in over
      # half of the samples. R and biotapy name nodes differently, so a balance is matched by its partition:
      # the taxa in its numerator, then in its denominator.
      sbp = pd.read_csv(GOLDEN / "philr_sbp.csv.gz", dtype={"balance": str, "taxon_id": str})
      golden = pd.read_csv(GOLDEN / "philr.csv.gz", dtype={"sample_id": str, "balance": str})
      r_names = {
          (frozenset(group.loc[group["sign"] > 0, "taxon_id"]), frozenset(group.loc[group["sign"] < 0, "taxon_id"])): name
          for name, group in sbp.groupby("balance")
      }
      gp = bt.datasets.global_patterns()
      # TreeData subsetting keeps one-child nodes, which ape::drop.tip removed in R.
      tdata = gp[:, gp.var_names.isin(sbp["taxon_id"])].copy()
      out = bt.pp.philr(tdata).obsm["X_philr"]
      tree = get_skbio_tree(tdata)
      ours = {
          tuple(frozenset(tip.name for tip in child.tips(include_self=True)) for child in tree.find(name).children): name
          for name in out.columns
      }
      # Equal partitions with equal numerators: the children are in R's order, so the signs agree.
      assert ours.keys() == r_names.keys()
      renamed = out.rename(columns={name: r_names[key] for key, name in ours.items()})
      expected = golden.pivot(index="sample_id", columns="balance", values="value")
      # atol: a balance between taxa that are all absent from a sample is about 1e-16 on both sides.
      np.testing.assert_allclose(renamed.loc[expected.index, expected.columns], expected, rtol=1e-7, atol=1e-12)
  ```
  Append to `tests/core/test_tree.py`:
  ```python
  def test_get_skbio_tree_can_keep_a_wide_root():
      tree = get_skbio_tree(bt.datasets.toy(), split_root=False)
      assert [child.name for child in tree.children] == ["n1", "n2", "n3"]
  ```
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/pp/test_philr.py tests/core/test_tree.py -q`
  -> `18 failed, 33 passed` (`AttributeError: module 'biotapy.pp' has no
  attribute 'philr'`; `TypeError: get_skbio_tree() got an unexpected keyword
  argument 'split_root'`); the golden test with `-m "golden or network"` ->
  `1 failed`.
- [x] **Step 3: Implement.** In `src/biotapy/_core/_tree.py`:
  ```diff
  @@ -65,15 +65,16 @@ def get_tree(tdata: TreeData) -> nx.DiGraph[str]:
       return cast("nx.DiGraph[str]", tdata.vart[PHYLO_KEY])


  -def get_skbio_tree(adata: AnnData) -> TreeNode:
  +def get_skbio_tree(adata: AnnData, *, split_root: bool = True) -> TreeNode:
       """The phylogeny in ``vart['phylo']`` as a scikit-bio ``TreeNode``, rooted where it is drawn.

       scikit-bio's Faith PD and UniFrac accept a root with at most two children. A
       root with more (an unrooted Newick tree, or the toy tree) keeps its first child
       and gets the others under one new zero-length node, which changes no
  -    root-to-tip path length. Missing (NaN) branch lengths stay NaN; scikit-bio
  -    counts them as zero. A plain AnnData raises ``TypeError``; a TreeData without
  -    ``vart['phylo']`` raises ``KeyError``.
  +    root-to-tip path length; ``split_root=False`` leaves it as stored (PhILR must
  +    refuse it). Children keep the order of the stored edges. Missing (NaN) branch
  +    lengths stay NaN; scikit-bio counts them as zero. A plain AnnData raises
  +    ``TypeError``; a TreeData without ``vart['phylo']`` raises ``KeyError``.
       """
       if not isinstance(adata, TreeData):
           msg = f"needs a TreeData with a tree in vart[{PHYLO_KEY!r}], got {type(adata).__name__}"
  @@ -85,7 +86,7 @@ def get_skbio_tree(adata: AnnData) -> TreeNode:
           nodes[child] = TreeNode(name=child, length=tree.edges[parent, child]["length"])
           nodes[parent].append(nodes[child])
       top = nodes[root]
  -    if len(top.children) > 2:
  +    if split_root and len(top.children) > 2:
           split = TreeNode(length=0.0)
           split.extend(top.children[1:])
           top.append(split)
  ```
  `src/biotapy/pp/_philr.py`:
  ```python
  """PhILR: isometric log-ratio balances over the phylogeny."""

  import pandas as pd
  from skbio import TreeNode
  from skbio.stats.composition import clr, tree_basis

  from biotapy._core import TreeData, add_provenance, get_skbio_tree

  from ._transform import pseudocounted


  def philr(tdata: TreeData, *, pseudocount: float = 0.5) -> TreeData:
      """Add the phylogenetic isometric log-ratio transform (PhILR) as ``obsm['X_philr']``.

      Parameters
      ----------
      tdata
          Samples x features with a rooted binary phylogeny in ``vart['phylo']``;
          ``X`` holds counts or another non-negative abundance.
      pseudocount
          Added to every value of ``X`` before the logarithm, so zeros have one.

      Returns
      -------
      TreeData
          A copy of ``tdata`` with ``obsm['X_philr']``: a samples x balances
          ``pandas.DataFrame`` with one column per internal node with two children,
          named after the node in ``vart['phylo']`` (not as R's ``makeNodeLabel``
          names it), in preorder. A balance is positive when its node's first
          child is more abundant than its second.

      Raises
      ------
      TypeError
          ``tdata`` is an AnnData that is not a TreeData.
      KeyError
          ``tdata`` has no ``vart['phylo']``.
      TypeError
          ``pseudocount`` is a bool or not a real number.
      ValueError
          A node of the tree, the root included, has more than two children; a
          feature is not a tip of the tree; there are fewer than two features; or
          ``pseudocount`` or ``X`` is invalid, as in :func:`biotapy.pp.clr`.

      Warns
      -----
      UserWarning
          ``pseudocount`` is larger than the smallest non-zero value in ``X``.

      Notes
      -----
      R equivalent: ``philr::philr``, ``mia::transformAssay``
      Guide: :doc:`/guide/transforms`

      Equals ``philr::philr(x, tree, pseudocount = 0.5)`` with its default uniform part and
      ILR weights, which ``mia::transformAssay(method = "philr")`` calls; the weights are
      not offered. Each balance of node ``i`` is
      ``sqrt(r s / (r + s)) * log(g(first child) / g(second child))``, where ``r`` and
      ``s`` count the taxa under each child and ``g`` is their geometric mean.
      Subsetting a TreeData keeps one-child nodes; they define no balance and are
      skipped, as ``ape::drop.tip`` removes them in R. A wider node has no single
      balance: an unrooted tree's three-child root must be rooted, and polytomies
      resolved, before PhILR (for example with ``ape::multi2di``), as ``philr``
      requires. Balances are Euclidean coordinates: distances between samples equal
      Aitchison distances, whatever the tree.

      ``X`` is densified once (8 bytes x samples x features); the balances take
      8 bytes x samples x (features - 1). Peak memory is about five such arrays
      (4.8x measured on a 400 x 512 table) plus scikit-bio's sparse basis, which
      holds one value per tip under each node (113 MB on GlobalPatterns' 26 x 19,216,
      28 arrays), so budget for that on large tables.

      References
      ----------
      Silverman JD, Washburne AD, Mukherjee S, David LA (2017) A phylogenetic transform enhances
      analysis of compositional microbiota data. eLife 6:e21887.

      Examples
      --------
      >>> import biotapy as bt
      >>> tdata = bt.datasets.toy()[:, :6].copy()  # toy()'s root has three children; f1-f6 sit under two
      >>> out = bt.pp.philr(tdata)
      >>> out.obsm["X_philr"].columns.tolist()
      ['root', 'n1', 'n4', 'n2', 'n5']
      """
      if tdata.n_vars < 2:
          msg = f"pp.philr needs at least two features, got {tdata.n_vars}"
          raise ValueError(msg)
      tree = _binary_tree(get_skbio_tree(tdata, split_root=False))
      tips = pd.Index([tip.name for tip in tree.tips()])
      if len(tips) != tdata.n_vars:
          msg = f"pp.philr needs every feature to be a tip of the tree; {tdata.n_vars - len(tips)} feature(s) are not"
          raise ValueError(msg)
      values = pseudocounted(tdata, pseudocount, func="pp.philr", columns=tdata.var_names.get_indexer(tips))
      basis, nodes = tree_basis(tree)
      # tree_basis puts a node's first child in the denominator; philr::philr puts it in the numerator.
      balances = pd.DataFrame(-(clr(values) @ basis.T), index=tdata.obs_names, columns=nodes)
      out = tdata.copy()
      out.obsm["X_philr"] = balances[[node.name for node in tree.preorder() if not node.is_tip()]]
      add_provenance(out, "pp.philr", pseudocount=pseudocount)
      return out


  def _binary_tree(tree: TreeNode) -> TreeNode:
      """A copy of ``tree`` without one-child nodes, children in order; raise if a node has more than two."""
      kept: dict[int, TreeNode] = {}
      wide: list[str] = []
      for node in tree.postorder(include_self=True):
          children = [kept[id(child)] for child in node.children]
          if len(children) > 2:
              wide.append(str(node.name))
          # TreeNode.prune would also drop one-child nodes, but it moves the child to the end of its parent's list.
          kept[id(node)] = children[0] if len(children) == 1 else TreeNode(node.name, children=children)
      if wide:
          msg = (
              f"pp.philr needs a rooted binary tree, but the node(s) {wide[:3]} have more than two children "
              f"({len(wide)} in total); "
              "root the tree and resolve its polytomies first (for example with ape::multi2di in R)"
          )
          raise ValueError(msg)
      return kept[id(tree)]
  ```
  `src/biotapy/pp/__init__.py`:
  ```python
  from ._filter import filter_features, filter_samples
  from ._glom import tax_glom
  from ._philr import philr
  from ._rarefy import rarefy
  from ._transform import clr, relative

  __all__ = ["clr", "filter_features", "filter_samples", "philr", "rarefy", "relative", "tax_glom"]
  ```
  In `pyproject.toml` (`[tool.mypy]`; `tree_basis` is unannotated, so mypy
  strict reports `Call to untyped function "tree_basis" in typed context`):
  ```diff
  @@ -176,9 +176,16 @@ lint.pylint.max-statements = 30
   files = [ "docs/extensions", "src/biotapy" ]
   python_version = "3.12"
   # biom-format and scikit-learn ship no type annotations at all, nor do scikit-bio's
  -# subsample_counts, threadpoolctl's threadpool_limits and mudata's MuData.update; exempt only calls into them (not our
  -# own code) from strict's disallow_untyped_calls (mypy/checkexpr.py matches by callee fullname).
  -untyped_calls_exclude = [ "biom", "mudata", "skbio.stats._subsample", "sklearn", "threadpoolctl" ]
  +# subsample_counts and tree_basis, threadpoolctl's threadpool_limits and mudata's MuData.update; exempt only calls into
  +# them (not our own code) from strict's disallow_untyped_calls (mypy/checkexpr.py matches by callee fullname).
  +untyped_calls_exclude = [
  +  "biom",
  +  "mudata",
  +  "skbio.stats._subsample",
  +  "skbio.stats.composition._base.tree_basis",
  +  "sklearn",
  +  "threadpoolctl",
  +]
   # rules.md R7.2: strict on the package; tests stay unannotated
   strict = true
   overrides = [
  ```
- [x] **Step 4: Run, expect pass** - the same commands -> `57 passed`; golden
  `1 passed`. Property test also under three extra Hypothesis seeds.
- [x] **Step 5: Docs.** Append to `docs/guide/transforms.md`:
  ```markdown
  ## PhILR

  `bt.pp.philr` turns a sample into balances along its phylogeny: one per internal node with two
  children, the log-ratio of the geometric means of the taxa under the node's two children, scaled so
  that the balances are orthonormal coordinates. They go to `obsm["X_philr"]`, a table with one
  column per balance:

  ```python
  tdata = bt.datasets.toy()[:, :6].copy()
  out = bt.pp.philr(tdata)  # pseudocount=0.5, as in pp.clr
  out.obsm["X_philr"]  # columns root, n1, n4, n2, n5
  ```

  A balance is positive when the node's first child is more abundant than its second, as in
  `philr::philr` with its default uniform weights. Distances between samples' balances equal
  their Aitchison distances, the Euclidean distances between their CLR vectors.

  PhILR needs a rooted binary tree. Filtering features leaves nodes with one child; those define
  no balance and are skipped. A node with three or more children raises: that includes the root
  of an unrooted tree (`bt.datasets.toy()` has one, which is why the example keeps f1-f6). Root
  the tree and resolve its polytomies before reading it, for example with `ape::multi2di` in R, or
  in Python with scikit-bio, which resolves them arbitrarily, as `multi2di` does:

  ```python
  from skbio import TreeNode

  tree = TreeNode.read("tree.nwk")  # root it first (outgroup or midpoint) if it is unrooted: your choice
  tree.bifurcate()  # in place; every node ends with two children
  tree.write("binary.nwk")
  # then read your table again with tree="binary.nwk"
  ```

  Balance names are the node names in `vart["phylo"]`, which differ from the ones `philr` in R
  makes with `makeNodeLabel`; to compare balances across tools, match them by their numerator and
  denominator taxa, not by name.

  Filtering features with `bt.pp.filter_features` drops `layers["clr"]` and `obsm["X_philr"]`;
  plain `adata[:, ...]` slicing keeps the old values, so transform after subsetting.
  ```python
  tdata = bt.datasets.toy()[:, :6].copy()
  out = bt.pp.philr(tdata)  # pseudocount=0.5, as in pp.clr
  out.obsm["X_philr"]  # columns root, n1, n4, n2, n5
  ```

  A balance is positive when the node's first child is more abundant than its second, as in
  `philr::philr` with its default uniform weights. Distances between samples' balances equal
  their Aitchison distances, the Euclidean distances between their CLR vectors.

  PhILR needs a rooted binary tree. Filtering features leaves nodes with one child; those define
  no balance and are skipped. A node with three or more children raises: that includes the root
  of an unrooted tree (`bt.datasets.toy()` has one, which is why the example keeps f1-f6). Root
  the tree and resolve its polytomies before reading it, for example with `ape::multi2di` in R.
  Filtering features afterwards drops `obsm["X_philr"]`, since the balances described the old
  features.
  ```
  In `docs/api.md` add `pp.philr` after `pp.filter_samples`. In
  `docs/contributing.md` the sentence ends "...against phyloseq, vegan, ape
  and picante, and `pp.philr` against philr."
- [x] **Step 6: Contracts.**
  - `data-model-slots.md`: the `obsm` row's keys gain "; `X_philr` from
    `pp.philr` (a samples x balances `DataFrame`, one column per internal tree
    node with two children, named after the node in `vart["phylo"]`, in preorder; R's
    names differ, so match balances across tools by their taxa)"; after the layer-adding Propagation
    row add `| Embedding-adding (`pp.philr`) | everything | nothing; adds `obsm["X_philr"]`, which a later feature change drops |`.
  - `r-golden-parity.md` statement 4, after the PCoA row: `| PhILR balances | matched by partition (the taxa in each numerator and denominator, so signs must agree too), then elementwise | `rtol=1e-7`; `atol=1e-12`, because a balance between absent taxa is about 1e-16 on both sides |`.
  - `tree-access.md`, the last gotcha gains: "`pp.philr` passes
    `split_root=False`: a split would invent a balance, so it raises instead.
    Children keep the order of the stored edges (ape's edge order for a
    phyloseq tree), which is what makes PhILR's signs equal `philr::philr`'s."

  Log line:
  ```text
  - **Update**: [data-model-slots](contracts/data-model-slots.md): `obsm["X_philr"]` from `pp.philr` and its embedding-adding propagation row. [r-golden-parity](contracts/r-golden-parity.md): PhILR balances are matched by partition. [tree-access](contracts/tree-access.md): `get_skbio_tree(split_root=False)` and the child-order guarantee PhILR's signs rely on.
  ```
- [x] **Step 7: Gate and commit**
  ```bash
  git add src/biotapy/_core/_tree.py src/biotapy/pp/_philr.py src/biotapy/pp/__init__.py tests/core/test_tree.py \
    tests/pp/test_philr.py tests/pp/test_philr_golden.py pyproject.toml docs/guide/transforms.md docs/api.md \
    docs/contributing.md .knowledge/contracts/data-model-slots.md .knowledge/contracts/r-golden-parity.md \
    .knowledge/contracts/tree-access.md .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  uvx prek run --all-files   # after staging: --all-files skips untracked files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git commit -m "feat(pp): add the PhILR transform"
  ```

### Checkpoint A - review slice 3A
- [x] Review the whole slice (superpowers:requesting-code-review) against
  every contract, pure-by-default, the Phase 3 and slice 3A review focus;
  then a fix pass, one commit per finding, each with a test.
  Record: review 0 Critical / 1 Important / 5 Minor (I1: `pl.bar` and
  `pl.heatmap` accepted signed layers such as `layers["clr"]`; fixed in `pl`
  through `pl/_common.py:table` and reported to the user). Fix pass
  `c920c7b..0de6cca` plus controller commit `6ade269`; re-review 14/14
  addressed.
- [x] Run the exit-gate check for 3A: `uv run --group test pytest -m "golden or network" tests/pp -q`
  and the full `uv run --group test pytest`. Record at `6ade269`: 8 passed,
  107 deselected for the golden/network run; 1039 passed, 25 deselected for
  `uv run --group test pytest -q -W error::UserWarning`.
- [x] Knowledge (codebase-map templates; R12.2-R12.4):
  - **Update `.knowledge/modules/pp.md`**: Responsibility gains CLR and
    PhILR; Entry points `_transform.py:clr`, `_transform.py:pseudocounted`,
    `_philr.py:philr`, `_philr.py:_binary_tree`; Invariants: pseudocount
    added to every value, the warning rule, dense `layers["clr"]`,
    `obsm["X_philr"]` as a DataFrame in preorder, philr's sign, one-child
    nodes skipped in order, wide nodes raise; Gotchas: `TreeNode.prune`
    reorders children, `tree_basis`'s first child is the denominator,
    `toy()`'s root has three children; Dependencies add `get_skbio_tree`,
    `warn_user`. Copy the description into `modules/index.md`.
  - **Update `.knowledge/modules/core.md`**: `get_skbio_tree`'s
    `split_root` and its one caller.
  - Log lines for both.
- [ ] Push the branch and open the PR only after the user approves that push
  (R13.3; the standing approval covered Phase 2 slice branches only). CI
  green, including docs and the network job.
- [ ] Ask the user to review slice 3A before slice 3B is expanded.

---
## Slice 3B - Native DA and consensus

**Goal:** without R, a user compares two groups with LinDA and ANCOM-BC2, gets
two tables with the same columns and units, and sees in one table and one plot
where they agree; LinDA equals `MicrobiomeStat::linda` to 1e-12 and ANCOM-BC2
equals `ANCOMBC::ancombc2` within measured tolerances, on the exit-gate
contrast.

**How slice 3B was checked.** Every file below was written into a scratch
clone of the repository at `1e64bd7` (master, after slice 3A) and replayed as
one commit per step that says "commit" (eight commits). Each committed state
was gated:

| Commit | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` |
|---|---|---|---|---|
| 3.5 `build(r)` (MicrobiomeStat) | passed | 1039 passed, 25 deselected | 32 passed | build succeeded |
| 3.5 `test(golden)` (LinDA) | passed | 1040 passed, 25 deselected | 32 passed | build succeeded |
| 3.5 `feat(da)` (LinDA, schema) | passed | 1068 passed, 27 deselected | 34 passed | build succeeded |
| 3.4 `build(r)` (ANCOMBC) | passed | 1068 passed, 27 deselected | 34 passed | build succeeded |
| 3.4 `test(golden)` (ANCOM-BC2) | passed | 1069 passed, 27 deselected | 34 passed | build succeeded |
| 3.4 `feat(da)` (ANCOM-BC2) | passed | 1094 passed, 29 deselected | 36 passed | build succeeded |
| 3.8 `feat(da)` (consensus) | passed | 1124 passed, 29 deselected | 36 passed | build succeeded |
| 3.9 `feat(pl)` (plot) | passed | 1133 passed, 29 deselected | 36 passed | build succeeded |

- prek covers ruff 0.16.9 check and format, `mypy --strict`, import-linter and
  pyproject-fmt. Every run exported `BIOTAPY_DATA_DIR` to a scratch pooch
  cache; `~/.cache/biotapy` was checked absent afterwards.
- Coverage on the final state: `da/_schema.py`, `da/_design.py`,
  `da/_linda.py`, `da/_ancombc.py`, `da/_consensus.py` and `pl/_consensus.py`
  are each 100%.
- The property tests also passed under Hypothesis seeds 1, 2 and 3 (and the
  two DA properties under seeds 4-40).
- One full run at the 3.8 commit, made while other test suites loaded the
  machine, reported `1 failed`; the failing test was not captured, and two
  quiet reruns at the same commit passed (`1124 passed`). The likeliest cause
  is a Hypothesis deadline (200 ms) under load; the consensus property's
  examples take 6-30 ms when the machine is idle. Recorded under the
  self-review's residual risks.
- The RED results in each "expect failure" step were reproduced by running the
  task's tests against the previous commit.
- The golden image was rebuilt from each edited `tests/r/Dockerfile` (tagged
  `biotapy-golden-p3b-ms` after the MicrobiomeStat commit and
  `biotapy-golden-p3b` after the ANCOMBC commit, locally, never pushed).
  `export_golden.R` ran twice on the final image: bit-identical, and every file
  from before slice 3B (38 goldens and fixtures) unchanged; `linda.csv.gz`
  from the MicrobiomeStat-only image equals the final image's byte for byte.
  The first ANCOMBC build failed twice: ANCOMBC 2.12.0's lazy load stops at
  `object 'solve' is not exported by 'namespace:CVXR'` (the CRAN snapshot has
  CVXR 1.8.2), then at `libgsl.so.27: cannot open shared object file`; the
  Dockerfile below carries both fixes.
- **Parity measured** (GlobalPatterns genera in >= 20% of samples, 636
  features, 9 human samples against 17 others; models `host` and
  `host + log_depth`):
  - `da.linda` vs `MicrobiomeStat::linda(is.winsor = FALSE)`: largest
    relative difference 5.3e-13 (`effect`), 5.0e-15 (`se`), 7.9e-13
    (`pvalue`), 7.9e-13 (`qvalue`); largest absolute difference 1.2e-14.
    Exact for practical purposes; the test keeps the contract's `rtol=1e-7`.
  - `da.ancombc2` vs `ANCOMBC::ancombc2`: the same 36 genera untested (NA in
    R with p = 1, NaN in biotapy). With `log_depth`: `effect` within 1.5e-8
    log2, `se` within 1.7e-10 relative, `pvalue` within 1.1e-8. Without it:
    every `effect` is shifted by a near-constant 0.002 to 0.012 log2 (mean
    0.008, the bias estimate), `se` within 9.0e-4 relative, `pvalue` within
    0.012; Spearman correlation of effects 0.999999; signs agree on 599 of
    600 (the other genus's effect is 0.005 log2 here and -0.004 in R); the
    same 208 genera called at `q < 0.05` when R's p-values are corrected over
    the tested genera. The shift comes from the bias E-M, which has not
    converged at its 100-iteration cap on this model in R or scikit-bio:
    scikit-bio at 99, 100, 101, 200 and 400 iterations gives mean shifts of
    0.012, 0.008, 0.004, -0.22 and -0.27 log2 against R's capped result, and R
    itself shifts its effects by -0.28 log2 on average with `em_control =
    list(max_iter = 1000)` (its own q-values then call 220 genera, not 208).
- APIs checked in the installed versions: scikit-bio 0.7.4
  (`ancombc2(table, metadata, formula, ..., p_adjust=None)`, its `result`
  MultiIndex `(FeatureID, Covariate)`, natural-log `Log(FC)`, p = 1 and NaN
  `Log(FC)` for a feature with one observed group level, patsy ordering
  numeric terms after categorical ones), scipy 1.18.1
  (`false_discovery_control(p, method="bh")`, empty input allowed,
  `stats.t.sf`), numpy 2.5.3 (`linalg.lstsq`, `percentile` type 7 = R's
  `IQR`), pandas 3.0.6 (`str` dtype for string columns,
  `Categorical.reorder_categories`, `assert_frame_equal(rtol=)`). R side, read
  at source in the image: MicrobiomeStat 1.4 `linda`, modeest 2.4.0
  `mlv.default`, `meanshift`, `venter(type = "shorth")` and `.deal.ties`,
  statip 0.2.3 `.kernel.gaussian` (`dnorm`), stats `bw.nrd0`, ANCOMBC 2.12.0
  `ancombc2`, `.ancombc2_core`, `.iter_mle`, `.bias_em`, `.data_core`.

### Slice 3B design

- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `tests/r/Dockerfile` | MicrobiomeStat (3.5); ANCOMBC with CVXR 1.0-15 and `libgsl27` (3.4) |
  | `tests/r/export_golden.R` | the benchmark data and the LinDA (3.5) and ANCOM-BC2 (3.4) sections |
  | `tests/golden/global_patterns/{linda,ancombc2}.csv.gz` | the goldens |
  | `da/__init__.py` | imports and `__all__` only (R4.1) |
  | `da/_schema.py` | `result()` builds the table (3.5); `COLUMNS` and `validate_result()` check one (3.8) |
  | `da/_design.py` | `model()`, `design_matrix()`, `dense_counts()`: the shared argument checks and inputs (3.5) |
  | `da/_linda.py` | `linda`, `_fit`, `_mode`, `_shorth` (3.5) |
  | `da/_ancombc.py` | `ancombc2` (3.4) |
  | `da/_consensus.py` | `consensus`, `_check` (3.8) |
  | `pl/_consensus.py` | `consensus` (3.9) |
  | `decisions/da-consensus-agreement.md` | what "agree" means (3.8) |

  `da/_design.py` and `da/_schema.py` serve the `da` topic files only, so they
  stay in `da` (R4.3, module-boundaries rule 2); `pl.consensus` reads the
  consensus table's columns and imports nothing from `da` (same layer).
  `pp.pseudocounted` stays in `pp`: LinDA's zero rule is different (0.5 added
  only when `X` holds a zero, a constant with nothing to validate and no
  scale warning), so the two share no logic and nothing moves to `_core`.
- **The schema module.** `result(features, *, effect, se, pvalue, method,
  contrast)` computes `qvalue` as BH over the finite p-values
  (`scipy.stats.false_discovery_control`) and `direction` as the int8 sign of
  `effect` (0 for NaN), index `feature`, in the given order. Methods call it;
  they do not re-validate their own output (R3.5: validation at the public
  boundary). `validate_result(table, *, arg)`, added with its only consumer in
  3.8, raises `TypeError` for a non-DataFrame and `ValueError` when a schema
  column is missing, a float column or `direction` has another dtype, the
  index repeats a feature, `pvalue`/`qvalue` leave [0, 1], `qvalue` or
  `effect` is not NaN exactly where `pvalue` is, `direction != sign(effect)`
  (0 where `effect` is NaN), or the table holds more than one `method` or
  `contrast` (a NaN `method` or `contrast` counts as another). Extra columns are allowed and ignored.
- **The shared model (`da/_design.py`).** `model(adata, group, *, covariates,
  reference, func)` returns `obs[[group, *covariates]]` with numeric columns
  as float64 and others as a `Categorical` of the levels present (the group's
  `reference` first), and the `contrast` text. It raises when `covariates` is
  not a list of column names or `reference` is not a string (`TypeError`), a
  column is absent (`KeyError`), a column repeats or is constant, a value is
  missing, the group's levels are not two, `reference` is not a level or is
  given for a numeric group, there are fewer than two features, there are no
  more samples than model terms, or the design is rank deficient (collinear covariates: `lstsq` would
  silently return a minimum-norm answer). `design_matrix(frame, *, scale)`
  builds intercept + indicators against the first category + numeric values
  (centred and divided by the SD with `scale=True`, as R's `scale()`), with
  the group always in column 1. `dense_counts(adata, *, func)` runs
  `require_counts`, refuses an empty sample with a message naming it and the
  fix (R's `linda` stops with `missing value where TRUE/FALSE needed`;
  scikit-bio's `ancombc2` with `Input matrix cannot have compositions without
  observed components`) and densifies `X` once.
- **LinDA, exactly as MicrobiomeStat 1.4 computes it with fixed effects and
  `is.winsor = FALSE`.**
  1. Numeric `obs` columns are scaled (`scale()`); categorical columns become
     treatment indicators.
  2. If `X` holds a zero, 0.5 is added to every count. MicrobiomeStat's
     adaptive switch regresses log library size on the design and prints
     "Imputation approach is used." when a p-value is <= 0.1, but then tests
     `zero.handling == 'imputation'` against the value `"Imputation"`, so the
     pseudocount path always runs; the standalone LinDA 0.2.0 does impute.
     Measured in the image: on `host + log_depth` (which prints the
     imputation message) MicrobiomeStat's output equals the pseudocount path
     exactly and differs from the imputation path by up to 1.48 log2.
     biotapy computes what MicrobiomeStat returns (slice 3B decision 1).
  3. `W = log2(Y)` centred per sample; one least-squares fit per feature
     (`numpy.linalg.lstsq`), `df = n - p`, `se = sqrt(sigma2 * (X'X)^-1[1,1])`.
  4. `bias = mlv(sqrt(n) * beta, "meanshift", "gaussian") / sqrt(n)`, ported
     from modeest 2.4.0: bandwidth `bw.nrd0` (with its zero-spread
     fallbacks); start at `shorth` (mean of the shortest window of
     `ceil(n/2)` sorted values; tied windows take the mean of their 1-based
     starts, which R's indexing truncates); Gaussian mean shift until the
     relative step is below `sqrt(eps)`, returning the *previous* value, at
     most 1000 steps.
  5. `effect = beta - bias`, `p = 2 * t.sf(|effect / se|, df)`, BH.
- **ANCOM-BC2.** `skbio.stats.composition.ancombc2(counts, metadata,
  formula, p_adjust=None)` with `counts` a DataFrame of the dense table
  (scikit-bio rejects sparse input) and `metadata` the model frame renamed to
  `x0..xk`, so patsy needs no quoting and categorical columns keep their
  levels (patsy's reference is the first category). The formula is
  `"x0 + x1 + ..."`; the group's row is the covariate named `x0` or
  `x0[...]` (by name: patsy puts numeric terms after categorical ones).
  `effect` and `se` are divided by ln 2; a feature scikit-bio could not fit
  (NaN `Log(FC)`, p = 1, as R does) or with no residual degrees of freedom
  (a finite effect and NaN p) gets NaN in all four floats, so BH runs over the
  tested features only. Settings equal R's `ancombc2(..., p_adj_method = "BH",
  prv_cut = 0, lib_cut = 0, pseudo_sens = FALSE, struc_zero = FALSE)` with its
  other defaults (`pseudo = 0`, `s0_perc = 0.05`, `iter_control` 1e-2/20,
  `em_control` 1e-5/100): scikit-bio's `tol`/`max_iter` are R's `em_control`,
  and its sparse ML step uses R's `iter_control` defaults internally (design
  note 7's `iter_control = list(tol = 1e-5, max_iter = 100)` named the wrong
  list). `alpha` only sets the `Signif` flag the schema drops, and `neg_lb`
  only tunes the structural-zero test, which is off; scikit-bio's standalone
  `struc_zero` is not called.
- **Consensus table.** Columns, per input in order: `effect_<method>`,
  `qvalue_<method>`, `significant_<method>` (bool, `qvalue < alpha`); then
  `n_tested` (methods with a finite `pvalue`), `n_significant`, `direction`
  (int8: the sign every call shares; 0 for no call or a conflict),
  `consensus`, `conflict`. Index: the union of the tables' features in
  first-seen order, named `feature`. `significant_<method>` is added to design
  note 6's list so the plot needs no `alpha` (slice 3B decision 6).
- **`pl.consensus(table, *, top=30, ax=None)`.** Rows: features with
  `n_significant > 0`, sorted by `n_significant` then mean |effect| (both
  descending; stable), at most `top`, the first at the top. Columns: methods
  from the `significant_<method>` columns. A filled dot per call, tab10 red
  for `effect > 0` and blue for `effect < 0`; a hollow grey dot
  (`_common.MISSING_COLOR`) where the method tested but did not call; nothing
  where it did not test. Bold y labels for consensus rows. One scatter per
  kind with a legend entry; `_common.new_axes` and `_common.label_ticks`
  reused. A table with no call raises `ValueError`; a table without the
  consensus columns raises `KeyError`.
- **Goldens (`export_golden.R`).** `tax_glom(GlobalPatterns, "Genus")`,
  `filter_taxa(..., sum(x > 0) >= 0.2 * length(x))` (636 genera: what
  `bt.pp.filter_features(min_prevalence=0.2)` keeps after `bt.pp.tax_glom`,
  checked by the tests' index comparison), `host` = `"human"` for Feces, Skin
  and Tongue (9 samples) and `"other"` for the 17 others, reference
  `"other"`; `log_depth = log(colSums(counts))` as the numeric covariate of
  the second model. LinDA and ANCOM-BC2 are deterministic; the ANCOM-BC2
  section still sets `set.seed(20260927)` because its E-M loop runs under
  `%dorng%`. Files: `linda.csv.gz` (`formula, taxon_id, log2FoldChange,
  lfcSE, pvalue, padj`, 1,272 rows, 47.5 KB) and `ancombc2.csv.gz`
  (`formula, taxon_id, lfc, se, p, q`, natural logs, 45.8 KB).
- **Docs.** One guide page, `docs/guide/differential_abundance.md`, grown by
  each task (the schema and LinDA, ANCOM-BC2, agreement, the plot), listed in
  the guide toctree; `api.md` gains a "Differential abundance" section and
  `pl.consensus`. Slice 3D adds `docs/methods/` (one page per method with its
  model and R defaults), the GPL note for the `r` extra and the exit-gate
  notebook; 3B keeps the method math in the docstrings' Notes and the guide.
- **Tests.** Behaviours every method shares (reference, missing values,
  levels, collinearity, empty samples, raw counts) live in
  `tests/da/test_design.py` and the schema's invariants in
  `tests/da/test_schema.py`, both parametrized over `METHODS`; 3.4 adds
  `bt.da.ancombc2` to `METHODS`, which is how the Phase 3 review focus's
  "3.4 the same two" tests exist. R11.2's "missing taxonomy rank" edge case
  does not apply: `da` reads no rank.

### Slice 3B global constraints (in addition to the Phase 3 list)
- Golden tests carry `golden` + `network` (GlobalPatterns is downloaded by
  pooch) and use the session fixture `benchmark` in `tests/da/conftest.py`,
  which they must not mutate.
- Tests reach functions through `bt.da.<fn>` and `bt.pl.consensus`; no test
  imports `biotapy.da._*`.
- `da` imports nothing from `pl`, and `pl/_consensus.py` nothing from `da`
  (import-linter's layers contract enforces it).
- Run commands: `uv run --group test pytest <path> -q`; golden tests with
  `-m "golden or network"`. Gate before every commit: `uvx prek run
  --all-files` after staging new files; with docs changes also
  `uv run --group doc sphinx-build -W -b html docs docs/_build/html`.
- Commits stage explicit paths only (never `.claude/`, `.superpowers/`,
  `.worktrees/`, `notebooks/`, `build/`). Each task's last commit stages
  `.knowledge/roadmap/phase-3-stats.md` with its boxes ticked and
  `.knowledge/log.md` with its line, under one heading
  `## <date of the commit> (Phase 3, slice 3B)` at the top of the log.

### Slice 3B review focus
The five ways a user is most likely to get a wrong DA answer from slice 3B
without an error, each pinned by a named test:

1. **The sign is the other way round.** R and patsy take the alphabetically
   first level as reference. Expected: `reference=` sets it for every method,
   `contrast` says it, and a string group defaults to its first sorted level.
   Tests: `tests/da/test_design.py::test_reference_sets_the_sign`,
   `test_string_group_takes_the_first_sorted_level_as_reference` (both
   methods).
2. **Samples silently dropped or a model silently degenerate.** Expected: a
   missing value, a third level, a collinear covariate, too few samples or an
   empty sample raise. Tests: `test_design.py::test_missing_group_value_raises`,
   `test_missing_covariate_value_raises`, `test_three_levels_raise`,
   `test_collinear_covariate_raises`, `test_single_sample_raises`,
   `test_all_zero_sample_raises` (both methods).
3. **Units and corrections that differ per column.** ANCOM-BC2 reports
   natural logs and Holm by default. Expected: log2 effects and BH q-values
   from every method. Tests: `tests/da/test_ancombc.py::test_ancombc2_is_scikit_bio_in_log2`,
   `tests/da/test_schema.py::test_qvalue_is_benjamini_hochberg` (both methods).
4. **An untested feature counted as tested.** ANCOM-BC2 cannot fit a feature
   absent from one group; R counts it with p = 1. Expected: NaN, left out of
   BH, and counted as "not tested" by consensus; LinDA drops nothing. Tests:
   `test_ancombc.py::test_feature_absent_from_a_group_is_not_tested`,
   `tests/da/test_linda.py::test_no_feature_is_dropped`,
   `tests/da/test_consensus.py::test_untested_is_not_counted_as_not_significant`.
5. **A consensus that is not one.** Expected: strict `q < alpha`, opposite
   signs are a conflict and never a consensus, and tables of different
   contrasts or the same method raise. Tests: `test_consensus.py::test_q_equal_to_alpha_is_not_called`,
   `test_conflict_is_flagged_and_never_consensus`,
   `test_results_with_different_contrasts_raise`, `test_repeated_method_raises`.

Plus purity (R3.3): `test_linda_keeps_input`, `test_ancombc2_keeps_input`,
`test_consensus_keeps_its_inputs`, `tests/pl/test_consensus.py::test_consensus_keeps_the_table`.

Execution order: **3.B0 -> 3.5 -> 3.4 -> 3.8 -> 3.9 -> Checkpoint B.** 3.B0
is a `_core` fix the user approved with this plan (2026-10-05); it adds two
tests, so every later `Expected: N passed` count in this slice is two higher
than the prototype's table above. LinDA first: it
is native end to end, so the schema and the shared model are fixed by code
biotapy owns before a wrapper maps a library's columns onto them; consensus
needs two methods; the plot draws the consensus table.

---

### Task 3.B0: counts are non-negative (`fix(core)`)

Found while planning 3B (not a 3B task; approved by the user 2026-10-05 as its
own commit, like Phase 2's `pp.relative` fix). `_core.infer_x_kind` calls any
table of whole numbers `"counts"`, negative ones included, and
`_core.require_counts` reuses that rule ("one definition of counts"), so
`pp.rarefy`, `tl.alpha`, `tl.unifrac(weighted=True)` and, from 3.5,
`da.linda` and `da.ancombc2` accept a table holding `-3`. Root fix: counts
are non-negative whole numbers, in the one rule both functions share.

**Files:** modify `src/biotapy/_core/_slots.py`, `tests/core/test_slots.py`,
`.knowledge/contracts/data-model-slots.md` (line ~58), `.knowledge/modules/core.md`
(lines ~56 and ~64), `docs/guide/reading_data.md` (line ~8),
`docs/guide/diversity.md` (line ~48), `docs/guide/filtering.md` (line ~52),
`.knowledge/log.md`.
**Not touched:** every caller of `require_counts`/`infer_x_kind` (they inherit
the fix), `RELATIVE_TOLERANCE`, the relative/abundance branches.
**Interfaces:**
- Consumes: nothing new.
- Produces: `infer_x_kind(X)` returns `"counts"` only when every stored value
  is a whole number `>= 0`; `require_counts(adata, *, func)` raises
  `ValueError` with `"{func} needs raw counts in X, but X holds non-integer,
  negative or missing (NaN) values"` for those tables. Later tasks rely on
  this; 3.5's `dense_counts` calls `require_counts` unchanged.

- [x] **Step 1: Write the failing tests.** Append to `tests/core/test_slots.py`
  (next to the other `infer_x_kind` and `require_counts` tests):
  ```python
  def test_infer_x_kind_negative_whole_numbers_are_not_counts():
      assert infer_x_kind(sp.csr_matrix(np.array([[1.0, -2.0], [2.0, 3.0]]))) == "abundance"


  def test_require_counts_rejects_negative_values():
      adata = _adata()
      adata.X = sp.csr_matrix(np.array([[0.0, -2.0, 1.0], [3.0, 4.0, 5.0]]))
      with pytest.raises(ValueError, match="pp.rarefy needs raw counts in X, but X holds non-integer, negative"):
          require_counts(adata, func="pp.rarefy")
  ```
- [x] **Step 2: Run them, expect failure.**
  Run: `uv run --group test pytest tests/core/test_slots.py -q -k "negative"`
  Expected: 2 failed (`'counts' == 'abundance'`; `DID NOT RAISE`).
- [x] **Step 3: Fix the rule.** In `src/biotapy/_core/_slots.py`:
  - `infer_x_kind`: the counts test becomes
    `if np.all(matrix.data == np.round(matrix.data)) and not np.any(matrix.data < 0):`
    and the docstring's first rule reads "Non-negative whole numbers are
    ``"counts"``".
  - `require_counts`: docstring "every stored value a non-negative whole
    number"; message `f"{func} needs raw counts in X, but X holds non-integer,
    negative or missing (NaN) values"`; the comment keeps naming
    `infer_x_kind`'s rule as the one definition.
  - In `tests/core/test_slots.py`, the existing
    `test_require_counts_rejects_non_integer_values` match string becomes
    `"pp.rarefy needs raw counts in X, but X holds non-integer, negative or missing"`
    (the message changed; the assertion is not weakened, R11.5).
- [x] **Step 4: Run them, expect success.**
  Run: `uv run --group test pytest tests/core/test_slots.py -q`
  Expected: all passed. Then grep the repo for other assertions on the old
  message: `grep -rn "non-integer or missing" src tests docs .knowledge` must
  print only roadmap history (`.knowledge/roadmap/phase-*.md`); update any
  live test or doc it finds.
- [x] **Step 5: Docs and knowledge (R12.1).** "whole numbers are `counts`" ->
  "non-negative whole numbers are `counts`" in
  `.knowledge/contracts/data-model-slots.md`, `docs/guide/reading_data.md`;
  "whole numbers" -> "non-negative whole numbers" in
  `.knowledge/modules/core.md` (both lines), `docs/guide/diversity.md`,
  `docs/guide/filtering.md`. Refresh `generated.at` on the two concepts and
  add a dated line to `.knowledge/log.md` (newest heading). Do not edit
  roadmap history.
- [x] **Step 6: Commit and gate.**
  ```bash
  git add src/biotapy/_core/_slots.py tests/core/test_slots.py .knowledge/contracts/data-model-slots.md \
    .knowledge/modules/core.md .knowledge/log.md docs/guide/reading_data.md docs/guide/diversity.md docs/guide/filtering.md
  git commit -m "fix(core): count tables hold non-negative whole numbers"
  ```
  `git status --short` empty, then the four gates. Expected: prek passed;
  pytest `1041 passed` (1039 + 2); golden/network `32 passed`; sphinx build
  succeeded; `bash scripts/knowledge_stale.sh --against HEAD` reports no
  concept made stale by this commit.

---

### Task 3.5: `da.linda` and the result schema (3.3)

**Files:** modify `tests/r/Dockerfile`, `tests/r/export_golden.R`,
`tests/golden/VERSIONS.txt`, `src/biotapy/__init__.py`,
`tests/test_docstrings.py`, `docs/guide/index.md`, `docs/api.md`,
`docs/contributing.md`, `.knowledge/contracts/r-golden-parity.md`,
`.knowledge/contracts/data-model-slots.md`; create
`tests/golden/global_patterns/linda.csv.gz`, `src/biotapy/da/__init__.py`,
`src/biotapy/da/_schema.py`, `src/biotapy/da/_design.py`,
`src/biotapy/da/_linda.py`, `tests/da/conftest.py`, `tests/da/test_linda.py`,
`tests/da/test_linda_golden.py`, `tests/da/test_design.py`,
`tests/da/test_schema.py`, `docs/guide/differential_abundance.md`.
**Not touched:** `pp/_transform.py` (`pseudocounted` stays), `_core`, `pl`.
**Interfaces:**
- Consumes: `_core.require_counts`, `_core.as_csr`; GlobalPatterns (pooch).
- Produces: `bt.da.linda(adata: AnnData, group: str, *, covariates:
  Sequence[str] = (), reference: str | None = None) -> pd.DataFrame`; private
  `da._schema.result(features: pd.Index[str], *, effect, se, pvalue:
  NDArray[float64], method: str, contrast: str) -> pd.DataFrame`;
  `da._design.model(adata, group, *, covariates, reference, func: str) ->
  tuple[pd.DataFrame, str]`, `design_matrix(frame, *, scale: bool) ->
  NDArray[float64]`, `dense_counts(adata, *, func: str) -> NDArray[float64]`;
  the `benchmark` test fixture; `linda.csv.gz`.

- [x] **Step 1: Approval on record.** MicrobiomeStat (CRAN, GPL-3) in the
  golden image was approved with the Phase 3 plan (Phase 3 decision 6,
  2026-10-05);
  it lives only in the image. Nothing new to ask.
- [x] **Step 2: Add MicrobiomeStat to the image.** `tests/r/Dockerfile`:
  ```diff
  @@ -18,6 +18,8 @@ RUN Rscript -e 'install.packages("picante")'
   # chain installs fs as a P3M binary that links libuv at load time.
   RUN apt-get update && apt-get install -y --no-install-recommends libuv1 && rm -rf /var/lib/apt/lists/*
   RUN Rscript -e 'BiocManager::install("philr", version = "3.22", ask = FALSE, update = FALSE)'
  -RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE))'
  +# MicrobiomeStat (CRAN, GPL-3): the da.linda golden file; biotapy reimplements linda and never calls it.
  +RUN Rscript -e 'install.packages("MicrobiomeStat")'
  +RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE))'
   WORKDIR /work
   CMD ["Rscript", "tests/r/export_golden.R"]
  ```
  In `.knowledge/contracts/r-golden-parity.md` statement 1:
  ```diff
  @@ -20,8 +20,8 @@ sources:
      installs only what the current golden files need (`phyloseq`, which brings
      `Biostrings`, `vegan` and `ape`, plus CRAN `picante` for Faith PD and
      Bioconductor `philr` for `pp.philr`, with the `libuv1` runtime library its
  -   `fs` binary loads); a new golden function that needs another package adds
  -   it in its own commit (rules.md R2.3).
  +   `fs` binary loads, and CRAN `MicrobiomeStat` for `da.linda`); a new golden
  +   function that needs another package adds it in its own commit (rules.md R2.3).
   1b. HUMAnN golden files (`fn.func_glom`, `fn.renorm`) are produced by
      `tests/humann/export_golden.py`, run with
      `uv run --no-project --with humann==3.9 --with pandas==3.0.6`, never
  ```
  Add at the top of `.knowledge/log.md`:
  ```text
  ## <date of the commit> (Phase 3, slice 3B)
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs CRAN `MicrobiomeStat` for the `da.linda` golden file.
  ```
  Build: `docker build -t biotapy-golden tests/r`. Expected: the
  `stopifnot(requireNamespace(...))` layer passes;
  `docker run --rm biotapy-golden Rscript -e 'packageVersion("MicrobiomeStat")'`
  prints `1.4` (modeest 2.4.0 comes with it, as P3M binaries); phyloseq
  1.54.2, vegan 2.7.3, ape 5.8.1, picante 1.8.2, philr 1.36.0 and Matrix 1.7.4
  are unchanged.
- [x] **Step 3: Commit the image change on its own** (contract statement 1):
  ```bash
  git add tests/r/Dockerfile .knowledge/contracts/r-golden-parity.md .knowledge/log.md
  git commit -m "build(r): add MicrobiomeStat to the golden image"
  ```
- [x] **Step 4: Export.** In `tests/r/export_golden.R`, insert the benchmark
  data and the LinDA section before `## Synthetic phyloseq fixtures`, and add
  MicrobiomeStat and modeest to `VERSIONS.txt`:
  ```diff
  @@ -136,6 +136,29 @@ write_golden(
     file.path(gp, "philr_sbp.csv.gz")
   )

  +## Slice 3B golden files: differential abundance (GlobalPatterns genera, human hosts vs the rest)
  +# The genera in at least 20% of samples, as bt.pp.filter_features(min_prevalence=0.2) keeps them; feces, skin and
  +# tongue samples against the other 17; log_depth (log library size) is the numeric covariate of the second model.
  +gp_genus <- filter_taxa(tax_glom(GlobalPatterns, "Genus"), function(x) sum(x > 0) >= 0.2 * length(x), TRUE)
  +da_counts <- t(samples_as_rows(gp_genus))
  +human <- sample_data(gp_genus)$SampleType %in% c("Feces", "Skin", "Tongue")
  +da_meta <- data.frame(
  +  host = factor(ifelse(human, "human", "other"), levels = c("other", "human")),
  +  log_depth = log(colSums(da_counts)),
  +  row.names = colnames(da_counts)
  +)
  +da_formulas <- c("host", "host + log_depth")
  +# is.winsor = FALSE: biotapy does not winsorise. MicrobiomeStat prints "Imputation approach is used." for the second
  +# model but adds the 0.5 pseudocount in both (its switch tests "Imputation" == "imputation").
  +linda_rows <- lapply(da_formulas, function(f) {
  +  out <- MicrobiomeStat::linda(
  +    da_counts, da_meta, paste0("~", f), feature.dat.type = "count", is.winsor = FALSE, verbose = FALSE
  +  )$output$hosthuman
  +  data.frame(formula = f, taxon_id = rownames(out), log2FoldChange = out$log2FoldChange, lfcSE = out$lfcSE,
  +             pvalue = out$pvalue, padj = out$padj)
  +})
  +write_golden(do.call(rbind, linda_rows), file.path(gp, "linda.csv.gz"))
  +
   ## Synthetic phyloseq fixtures: biotapy's toy() numbers, no third-party data
   counts <- rbind(
     c(10, 5, 20, 30, 0, 2, 1, 0), c(8, 7, 25, 22, 3, 0, 0, 1), c(12, 4, 18, 35, 1, 5, 2, 0),
  @@ -197,5 +220,7 @@ writeLines(c(
     paste0("vegan ", packageVersion("vegan")),
     paste0("ape ", packageVersion("ape")),
     paste0("picante ", packageVersion("picante")),
  -  paste0("philr ", packageVersion("philr"))
  +  paste0("philr ", packageVersion("philr")),
  +  paste0("MicrobiomeStat ", packageVersion("MicrobiomeStat")),
  +  paste0("modeest ", packageVersion("modeest"))
   ), "tests/golden/VERSIONS.txt")
  ```
  Run twice and check bit identity (playbook regenerate-golden-files):
  ```bash
  mkdir -p build
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
  sha256sum tests/golden/*/*.csv.gz tests/data/phyloseq/* tests/data/dada2/* > build/golden-run1.sha
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
  sha256sum -c build/golden-run1.sha
  git status --short
  ```
  Expected: every line `OK` (each run takes about two minutes); `git status`
  shows only `M tests/golden/VERSIONS.txt` (two new lines, `MicrobiomeStat
  1.4` and `modeest 2.4.0`), `M tests/r/export_golden.R` and the new
  `tests/golden/global_patterns/linda.csv.gz` (47.5 KB). The run prints, besides
  the known `Found more than one class "phylo"` lines, `Warning message: In
  summary.lm(tmp) : essentially perfect fit` (LinDA's zero-handling check
  regresses log library size on `log_depth`); harmless.
- [x] **Step 5: Gate and commit.** `uv run --group test pytest
  tests/test_data_files.py -q` -> passes. Then `uvx prek run --all-files` and:
  ```bash
  git add tests/r/export_golden.R tests/golden/VERSIONS.txt tests/golden/global_patterns/linda.csv.gz
  git commit -m "test(golden): export the LinDA golden file"
  ```
- [x] **Step 6: Failing tests.** `tests/da/conftest.py`:
  ```python
  import numpy as np
  import pandas as pd
  import pytest

  import biotapy as bt


  # Session scope: GlobalPatterns takes seconds to load and glom; tests only read it.
  @pytest.fixture(scope="session")
  def benchmark():
      """The exit-gate data the golden files use: GlobalPatterns genera in >= 20% of samples, ``host`` and ``log_depth``."""
      tdata = bt.pp.filter_features(bt.pp.tax_glom(bt.datasets.global_patterns(), "genus"), min_prevalence=0.2)
      human = tdata.obs["SampleType"].isin(["Feces", "Skin", "Tongue"])
      tdata.obs["host"] = pd.Categorical(np.where(human, "human", "other"))
      tdata.obs["log_depth"] = np.log(np.asarray(tdata.X.sum(axis=1)).ravel())
      return tdata
  ```
  `tests/da/test_linda_golden.py`:
  ```python
  from pathlib import Path

  import numpy as np
  import pandas as pd
  import pytest

  import biotapy as bt

  GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
  pytestmark = [pytest.mark.golden, pytest.mark.network]


  @pytest.mark.parametrize(("formula", "covariates"), [("host", ()), ("host + log_depth", ("log_depth",))])
  def test_linda_matches_microbiomestat_linda(benchmark, formula, covariates):
      # R: MicrobiomeStat::linda(counts, meta, "~<formula>", is.winsor = FALSE)$output$hosthuman on the 636 genera.
      golden = pd.read_csv(GOLDEN / "linda.csv.gz", dtype={"taxon_id": str}).query("formula == @formula")
      out = bt.da.linda(benchmark, "host", covariates=covariates, reference="other")
      assert out.index.tolist() == golden["taxon_id"].tolist()
      assert (out["contrast"] == "human vs other").all()
      columns = {"effect": "log2FoldChange", "se": "lfcSE", "pvalue": "pvalue", "qvalue": "padj"}
      for ours, theirs in columns.items():
          np.testing.assert_allclose(out[ours], golden[theirs], rtol=1e-7, err_msg=ours)
  ```
  `tests/da/test_linda.py`:
  ```python
  import numpy as np
  import pandas as pd
  import scipy.sparse as sp
  from anndata import AnnData
  from scipy.stats import t as t_dist

  import biotapy as bt


  def _counts_adata(dense, **obs):
      return AnnData(
          X=sp.csr_matrix(dense),
          obs=pd.DataFrame(obs, index=[f"s{i}" for i in range(dense.shape[0])]),
          var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
      )


  def _ols(counts, design):
      """Group coefficient and its standard error for each feature, as R's lm on log2 CLR values."""
      values = counts + 0.5 if (counts == 0).any() else counts
      ratios = np.log2(values) - np.log2(values).mean(axis=1, keepdims=True)
      coef, *_ = np.linalg.lstsq(design, ratios, rcond=None)
      dof = design.shape[0] - design.shape[1]
      sigma2 = ((ratios - design @ coef) ** 2).sum(axis=0) / dof
      return coef[1], np.sqrt(sigma2 * np.linalg.inv(design.T @ design)[1, 1]), dof


  def test_linda_is_least_squares_on_log2_ratios_shifted_by_one_bias():
      tdata = bt.datasets.toy()
      out = bt.da.linda(tdata, "group")
      group = (tdata.obs["group"] == "B").to_numpy(np.float64)
      beta, se, dof = _ols(tdata.X.toarray().astype(float), np.c_[np.ones(6), group])
      np.testing.assert_allclose(out["se"], se, rtol=1e-12)
      # Every feature's effect is its coefficient minus the same bias, the mode of all coefficients.
      np.testing.assert_allclose(np.ptp(beta - out["effect"]), 0, atol=1e-12)
      np.testing.assert_allclose(out["pvalue"], 2 * t_dist.sf(np.abs(out["effect"] / out["se"]), dof), rtol=1e-12)


  def test_bias_correction_recovers_a_fold_change():
      # B repeats A with f5-f7 eight times as abundant: the CLR moves every other feature too, the bias removes that.
      rng = np.random.default_rng(0)
      a = rng.integers(20, 200, size=(4, 8)).astype(float)
      b = a * np.array([1, 1, 1, 1, 1, 8, 8, 8])
      out = bt.da.linda(_counts_adata(np.r_[a, b], g=["a"] * 4 + ["b"] * 4), "g")
      np.testing.assert_allclose(out["effect"].iloc[:5], 0, atol=0.05)
      np.testing.assert_allclose(out["effect"].iloc[5:], 3, atol=0.05)


  def test_pseudocount_is_added_only_when_x_has_a_zero():
      counts = np.array([[5.0, 9, 2], [4, 8, 3], [7, 2, 6], [6, 3, 9]])
      out = bt.da.linda(_counts_adata(counts, g=["a", "a", "b", "b"]), "g")
      _, se, _ = _ols(counts, np.c_[np.ones(4), [0.0, 0, 1, 1]])
      np.testing.assert_allclose(out["se"], se, rtol=1e-12)


  def test_covariates_enter_the_model():
      tdata = bt.datasets.toy()
      tdata.obs["batch"] = pd.Categorical(["x", "y", "y", "x", "y", "x"])
      tdata.obs["age"] = [30.0, 41, 52, 38, 45, 60]
      out = bt.da.linda(tdata, "group", covariates=["batch", "age"])
      group = (tdata.obs["group"] == "B").to_numpy(np.float64)
      batch = (tdata.obs["batch"] == "y").to_numpy(np.float64)
      age = tdata.obs["age"].to_numpy()
      design = np.c_[np.ones(6), group, batch, (age - age.mean()) / age.std(ddof=1)]
      _, se, _ = _ols(tdata.X.toarray().astype(float), design)
      np.testing.assert_allclose(out["se"], se, rtol=1e-12)


  def test_numeric_group_is_per_standard_deviation():
      tdata = bt.datasets.toy()
      tdata.obs["ph"] = [5.1, 5.6, 6.0, 6.8, 7.1, 7.4]
      out = bt.da.linda(tdata, "ph")
      tdata.obs["ph"] *= 10
      pd.testing.assert_frame_equal(bt.da.linda(tdata, "ph"), out, rtol=1e-9)
      assert (out["contrast"] == "ph").all()


  def test_no_feature_is_dropped():
      tdata = bt.datasets.toy()
      dense = tdata.X.toarray()
      dense[:, 4] = 0
      tdata.X = sp.csr_matrix(dense)
      out = bt.da.linda(tdata, "group")
      assert out.index.tolist() == tdata.var_names.tolist()
      # As in MicrobiomeStat: the all-zero feature is 0.5 everywhere, so it is tested like the others.
      assert np.isfinite(out[["effect", "se", "pvalue", "qvalue"]].to_numpy()).all()


  def test_linda_keeps_input(assert_unchanged):
      tdata = bt.datasets.toy()
      before = tdata.copy()
      bt.da.linda(tdata, "group")
      assert_unchanged(before, tdata)
  ```
  `tests/da/test_design.py` (3.4 adds `bt.da.ancombc2` to `METHODS`):
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp

  import biotapy as bt

  METHODS = [bt.da.linda]


  def _toy_with(**columns):
      tdata = bt.datasets.toy()
      for name, values in columns.items():
          tdata.obs[name] = values
      return tdata


  @pytest.mark.parametrize("method", METHODS)
  def test_reference_sets_the_sign(method):
      tdata = bt.datasets.toy()
      a_first, b_first = method(tdata, "group"), method(tdata, "group", reference="B")
      assert (a_first["contrast"] == "B vs A").all() and (b_first["contrast"] == "A vs B").all()
      np.testing.assert_allclose(b_first["effect"], -a_first["effect"], rtol=1e-9, atol=1e-12)
      assert (b_first["direction"] == -a_first["direction"]).all()
      np.testing.assert_allclose(b_first["pvalue"], a_first["pvalue"], rtol=1e-9)


  @pytest.mark.parametrize("method", METHODS)
  def test_string_group_takes_the_first_sorted_level_as_reference(method):
      tdata = _toy_with(site=["gut", "gut", "gut", "air", "air", "air"])
      assert (method(tdata, "site")["contrast"] == "gut vs air").all()


  @pytest.mark.parametrize("method", METHODS)
  def test_bool_group_compares_true_with_false(method):
      tdata = _toy_with(treated=[False] * 3 + [True] * 3)
      np.testing.assert_allclose(method(tdata, "treated")["effect"], method(tdata, "group")["effect"], rtol=1e-12)
      assert (method(tdata, "treated")["contrast"] == "True vs False").all()


  @pytest.mark.parametrize("method", METHODS)
  def test_missing_group_value_raises(method):
      tdata = _toy_with(group=pd.Categorical(["A", None, "A", "B", "B", "B"]))
      with pytest.raises(ValueError, match=r"obs\['group'\] is missing for 1 sample\(s\); drop them first"):
          method(tdata, "group")


  @pytest.mark.parametrize("method", METHODS)
  def test_missing_covariate_value_raises(method):
      tdata = _toy_with(age=[30.0, np.nan, 52, 38, 45, 60])
      with pytest.raises(ValueError, match=r"obs\['age'\] is missing for 1 sample"):
          method(tdata, "group", covariates=["age"])


  @pytest.mark.parametrize("method", METHODS)
  def test_three_levels_raise(method):
      tdata = _toy_with(site=["a", "b", "c", "a", "b", "c"])
      with pytest.raises(ValueError, match=r"compares two levels, but obs\['site'\] has 3: \['a', 'b', 'c'\]"):
          method(tdata, "site")


  @pytest.mark.parametrize("method", METHODS)
  def test_unused_categories_are_not_levels(method):
      tdata = _toy_with(group=pd.Categorical(["A"] * 3 + ["B"] * 3, categories=["A", "B", "C"]))
      assert (method(tdata, "group")["contrast"] == "B vs A").all()


  @pytest.mark.parametrize("method", METHODS)
  def test_unknown_reference_raises(method):
      with pytest.raises(ValueError, match=r"reference='C' is not a level of obs\['group'\], which has \['A', 'B'\]"):
          method(bt.datasets.toy(), "group", reference="C")


  @pytest.mark.parametrize("method", METHODS)
  def test_reference_for_a_numeric_group_raises(method):
      tdata = _toy_with(ph=[5.1, 5.6, 6.0, 6.8, 7.1, 7.4])
      with pytest.raises(ValueError, match=r"reference='A' needs a categorical group, but obs\['ph'\] is numeric"):
          method(tdata, "ph", reference="A")


  @pytest.mark.parametrize("method", METHODS)
  def test_absent_column_raises(method):
      with pytest.raises(KeyError, match=r"\['batch'\] not in obs"):
          method(bt.datasets.toy(), "group", covariates=["batch"])


  @pytest.mark.parametrize("method", METHODS)
  def test_covariates_as_a_string_raises(method):
      with pytest.raises(TypeError, match=r"covariates must be a list of obs columns, such as \['group'\]"):
          method(bt.datasets.toy(), "group", covariates="group")


  @pytest.mark.parametrize("method", METHODS)
  def test_group_among_covariates_raises(method):
      with pytest.raises(ValueError, match="repeat a column"):
          method(bt.datasets.toy(), "group", covariates=["group"])


  @pytest.mark.parametrize("method", METHODS)
  def test_collinear_covariate_raises(method):
      tdata = _toy_with(arm=["x"] * 3 + ["y"] * 3)
      with pytest.raises(ValueError, match="are collinear; drop the covariate that repeats another"):
          method(tdata, "group", covariates=["arm"])


  @pytest.mark.parametrize("method", METHODS)
  def test_single_sample_raises(method):
      tdata = _toy_with(ph=[5.1, 5.6, 6.0, 6.8, 7.1, 7.4])[:1].copy()
      with pytest.raises(ValueError, match=r"needs more samples \(1\) than model terms \(2, with the intercept\)"):
          method(tdata, "ph")


  @pytest.mark.parametrize("method", METHODS)
  def test_all_zero_sample_raises(method):
      tdata = bt.datasets.toy()
      dense = tdata.X.toarray()
      dense[1] = 0
      tdata.X = sp.csr_matrix(dense)
      with pytest.raises(ValueError, match=r"1 sample\(s\) have none \(\['s2'\]\).*bt.pp.filter_samples"):
          method(tdata, "group")


  @pytest.mark.parametrize("method", METHODS)
  def test_relative_abundances_raise(method):
      relative = bt.pp.relative(bt.datasets.toy())
      relative.X = relative.layers["relative"]
      relative.uns["biotapy"]["x_kind"] = "relative"
      with pytest.raises(ValueError, match="needs raw counts in X"):
          method(relative, "group")


  @pytest.mark.parametrize("method", METHODS)
  def test_single_feature_raises(method):
      with pytest.raises(ValueError, match=r"needs at least two features"):
          method(bt.datasets.toy()[:, :1].copy(), "group")


  @pytest.mark.parametrize("method", METHODS)
  def test_constant_numeric_group_raises(method):
      tdata = _toy_with(ph=[5.0] * 6)
      with pytest.raises(ValueError, match=r"group column 'ph' is constant across samples"):
          method(tdata, "ph")


  @pytest.mark.parametrize("method", METHODS)
  @pytest.mark.parametrize("values", [[3.0] * 6, ["x"] * 6])
  def test_constant_covariate_raises(method, values):
      tdata = _toy_with(c=values)
      with pytest.raises(ValueError, match=r"covariates column 'c' is constant across samples; drop it"):
          method(tdata, "group", covariates=["c"])


  @pytest.mark.parametrize("method", METHODS)
  @pytest.mark.parametrize("covariates", [None, [1], ("age", 2)])
  def test_covariates_must_be_a_sequence_of_names(method, covariates):
      with pytest.raises(TypeError, match=r"covariates must be a list of obs columns"):
          method(bt.datasets.toy(), "group", covariates=covariates)


  @pytest.mark.parametrize("method", METHODS)
  def test_reference_must_be_a_string(method):
      tdata = _toy_with(treated=[False] * 3 + [True] * 3)
      with pytest.raises(TypeError, match=r"reference must be the level's name as a string, such as 'False'"):
          method(tdata, "treated", reference=True)
  ```
  `tests/da/test_schema.py` (3.8 appends the `validate_result` tests):
  ```python
  import anndata as ad
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from hypothesis import given, settings
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays
  from scipy.stats import false_discovery_control

  import biotapy as bt

  METHODS = [bt.da.linda]
  COLUMNS = ["effect", "se", "pvalue", "qvalue", "direction", "method", "contrast"]


  @pytest.mark.parametrize("method", METHODS)
  def test_result_has_the_schema_columns_and_dtypes(method):
      out = method(bt.datasets.toy(), "group")
      assert out.columns.tolist() == COLUMNS and out.index.name == "feature"
      assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
      assert (out.dtypes.iloc[:4] == np.float64).all() and out["direction"].dtype == np.int8
      assert pd.api.types.is_string_dtype(out["method"]) and pd.api.types.is_string_dtype(out["contrast"])
      assert (out["method"] == method.__name__).all()


  @pytest.mark.parametrize("method", METHODS)
  def test_qvalue_is_benjamini_hochberg(method):
      out = method(bt.datasets.toy(), "group")
      tested = out["pvalue"].notna()
      np.testing.assert_allclose(
          out.loc[tested, "qvalue"], false_discovery_control(out.loc[tested, "pvalue"]), rtol=1e-12
      )


  @pytest.mark.parametrize("method", METHODS)
  @settings(max_examples=25, deadline=None)
  @given(
      counts=arrays(np.int64, (6, 7), elements=st.integers(0, 50)),
      split=st.integers(2, 4),
  )
  def test_direction_is_the_sign_and_qvalue_bounds_pvalue(method, counts, split):
      counts[:, 0] += 1  # no empty sample
      adata = ad.AnnData(
          X=sp.csr_matrix(counts),
          obs=pd.DataFrame({"g": ["a"] * split + ["b"] * (6 - split)}, index=[f"s{i}" for i in range(6)]),
          var=pd.DataFrame(index=[f"f{i}" for i in range(7)]),
      )
      out = method(adata, "g")
      tested = out["pvalue"].notna()
      assert out.index.tolist() == adata.var_names.tolist()
      assert (out["direction"] == np.sign(out["effect"].fillna(0))).all()
      assert (out.loc[tested, "qvalue"] >= out.loc[tested, "pvalue"] - 1e-15).all()
      assert out.loc[tested, ["pvalue", "qvalue"]].stack().between(0, 1).all()
      assert out.loc[~tested, ["effect", "qvalue"]].isna().all().all()
  ```
  `tests/test_docstrings.py`, the covered subpackages:
  ```diff
  @@ -15,7 +15,7 @@ PUBLIC = coming_from_r.public_functions()


   def test_every_public_subpackage_is_covered():
  -    assert {name.split(".")[1] for name, _ in PUBLIC} == {"datasets", "fn", "io", "pl", "pp", "tl"}
  +    assert {name.split(".")[1] for name, _ in PUBLIC} == {"da", "datasets", "fn", "io", "pl", "pp", "tl"}


   @pytest.mark.parametrize(("name", "function"), PUBLIC, ids=[name for name, _ in PUBLIC])
  ```
- [x] **Step 7: Run, expect failure** -
  `uv run --group test pytest tests/da tests/test_docstrings.py -q --continue-on-collection-errors`
  -> `8 failed, 45 passed, 2 deselected, 2 errors`: `test_design.py` and
  `test_schema.py` fail to collect and the seven `test_linda.py` tests fail with
  `AttributeError: module 'biotapy' has no attribute 'da'`, and
  `test_every_public_subpackage_is_covered` fails on the missing `"da"`.
  `uv run --group test pytest tests/da/test_linda_golden.py -q -m "golden or network"`
  -> `2 failed`.
- [x] **Step 8: Implement.** `src/biotapy/da/__init__.py`:
  ```python
  """Differential abundance: methods that share one result table (contracts/data-model-slots, DA results)."""

  from ._linda import linda

  __all__ = ["linda"]
  ```
  `src/biotapy/da/_schema.py`:
  ```python
  """The one result table every da method returns (contracts/data-model-slots, DA results)."""

  import numpy as np
  import numpy.typing as npt
  import pandas as pd
  from scipy.stats import false_discovery_control


  def result(
      features: "pd.Index[str]",
      *,
      effect: npt.NDArray[np.float64],
      se: npt.NDArray[np.float64],
      pvalue: npt.NDArray[np.float64],
      method: str,
      contrast: str,
  ) -> pd.DataFrame:
      """The schema table: one row per feature, BH ``qvalue`` over the finite p-values, ``direction`` = sign(effect)."""
      tested = np.isfinite(pvalue)
      qvalue = np.full(pvalue.shape, np.nan)
      qvalue[tested] = false_discovery_control(pvalue[tested], method="bh")
      return pd.DataFrame(
          {
              "effect": effect,
              "se": se,
              "pvalue": pvalue,
              "qvalue": qvalue,
              "direction": np.sign(np.nan_to_num(effect)).astype(np.int8),
              "method": method,
              "contrast": contrast,
          },
          index=pd.Index(features, name="feature"),
      )
  ```
  `src/biotapy/da/_design.py`:
  ```python
  """The model every da method fits: ``group``, ``covariates`` and ``reference`` checked once, without formulas."""

  from collections.abc import Sequence
  from typing import cast

  import numpy as np
  import numpy.typing as npt
  import pandas as pd
  from anndata import AnnData

  from biotapy._core import as_csr, require_counts


  def model(
      adata: AnnData, group: str, *, covariates: Sequence[str], reference: str | None, func: str
  ) -> tuple[pd.DataFrame, str]:
      """``obs[[group, *covariates]]`` ready to fit, and the ``contrast`` text.

      Numeric columns become float64; any other column a ``Categorical`` without
      unused levels, the group's ``reference`` first. Raises when ``covariates`` is
      not a list of column names, a column is missing, has missing values or is
      constant, the group does not have two levels, or the design has no residual
      degrees of freedom or collinear columns.
      """
      if (
          not isinstance(covariates, Sequence)
          or isinstance(covariates, str)
          or not all(isinstance(c, str) for c in covariates)
      ):
          example = covariates if isinstance(covariates, str) else "age"
          msg = f"{func}: covariates must be a list of obs columns, such as [{example!r}]"
          raise TypeError(msg)
      names = [group, *covariates]
      absent = [name for name in names if name not in adata.obs.columns]
      if absent:
          msg = f"{func}: {absent} not in obs; group and covariates name obs columns"
          raise KeyError(msg)
      if len(set(names)) < len(names):
          msg = f"{func}: group={group!r} and covariates={list(covariates)} repeat a column"
          raise ValueError(msg)
      # anndata types obs as DataFrame | Dataset2D (its lazy variant); the data model guarantees a DataFrame.
      obs = cast("pd.DataFrame", adata.obs)
      frame = pd.DataFrame({name: _column(obs[name], func=func) for name in names}, index=obs.index)
      frame[group], contrast = _group(frame[group], reference, func=func)
      design = design_matrix(frame, scale=False)
      if design.shape[0] <= design.shape[1]:
          msg = f"{func} needs more samples ({design.shape[0]}) than model terms ({design.shape[1]}, with the intercept)"
          raise ValueError(msg)
      _check_varies(frame, group, func=func)
      if np.linalg.matrix_rank(design) < design.shape[1]:
          msg = f"{func}: group={group!r} and covariates={list(covariates)} are collinear; drop the covariate that repeats another"
          raise ValueError(msg)
      return frame, contrast


  def design_matrix(frame: pd.DataFrame, *, scale: bool) -> npt.NDArray[np.float64]:
      """Intercept, then each column of ``frame``: indicators against the first category, or the values.

      With ``scale=True`` numeric columns are centred and divided by their standard
      deviation, as R's ``scale()``. The group's column is always column 1.
      """
      blocks = [np.ones((len(frame), 1))]
      for name in frame.columns:
          values = frame[name]
          if isinstance(values.dtype, pd.CategoricalDtype):
              codes = values.cat.codes.to_numpy()
              blocks.append((codes[:, None] == np.arange(1, len(values.cat.categories))).astype(np.float64))
          elif scale:
              blocks.append(((values - values.mean()) / values.std()).to_numpy(np.float64)[:, None])
          else:
              blocks.append(values.to_numpy(np.float64)[:, None])
      return np.hstack(blocks)


  def dense_counts(adata: AnnData, *, func: str) -> npt.NDArray[np.float64]:
      """``X`` as a dense float64 array of raw counts with no empty sample."""
      require_counts(adata, func=func)
      if adata.n_vars < 2:
          msg = f"{func} needs at least two features: log-ratio methods compare each feature with the others"
          raise ValueError(msg)
      X = as_csr(adata.X)
      empty = adata.obs_names[np.asarray(X.sum(axis=1)).ravel() == 0]
      if len(empty):
          msg = (
              f"{func} needs reads in every sample; {len(empty)} sample(s) have none ({empty[:3].tolist()}): "
              "drop them first, for example with bt.pp.filter_samples(adata, min_depth=1)"
          )
          raise ValueError(msg)
      # LinDA's log-ratios and scikit-bio's ancombc2 need a dense table (rules.md R6.2): the one dense copy of X.
      return X.toarray().astype(np.float64, copy=False)


  def _check_varies(frame: pd.DataFrame, group: str, *, func: str) -> None:
      """Raise for a constant column; a one-level categorical group is left to ``_group``'s two-level message."""
      for name in frame.columns:
          if frame[name].nunique() > 1 or (name == group and isinstance(frame[name].dtype, pd.CategoricalDtype)):
              continue
          if name == group:
              msg = f"{func}: group column {name!r} is constant across samples; there is nothing to compare"
          else:
              msg = f"{func}: covariates column {name!r} is constant across samples; drop it"
          raise ValueError(msg)


  def _column(values: pd.Series, *, func: str) -> pd.Series:
      """Float64 for a numeric column, else a ``Categorical`` of the levels present; missing values raise."""
      if values.isna().any():
          msg = f"{func}: obs[{values.name!r}] is missing for {int(values.isna().sum())} sample(s); drop them first"
          raise ValueError(msg)
      if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
          return values.astype(np.float64)
      return pd.Series(pd.Categorical(values).remove_unused_categories(), index=values.index, name=values.name)


  def _group(values: pd.Series, reference: str | None, *, func: str) -> tuple[pd.Series, str]:
      """The group column with ``reference`` as its first level, and the contrast text."""
      if not isinstance(values.dtype, pd.CategoricalDtype):
          if reference is not None:
              msg = f"{func}: reference={reference!r} needs a categorical group, but obs[{values.name!r}] is numeric"
              raise ValueError(msg)
          return values, str(values.name)
      levels = [str(level) for level in values.cat.categories]
      if len(levels) != 2:
          msg = f"{func} compares two levels, but obs[{values.name!r}] has {len(levels)}: {levels}; subset the samples first"
          raise ValueError(msg)
      if reference is None:
          reference = levels[0]
      if not isinstance(reference, str):
          msg = f"{func}: reference must be the level's name as a string, such as {levels[0]!r}, not {reference!r}"
          raise TypeError(msg)
      if reference not in levels:
          msg = f"{func}: reference={reference!r} is not a level of obs[{values.name!r}], which has {levels}"
          raise ValueError(msg)
      first = values.cat.categories[levels.index(reference)]
      ordered = values.cat.reorder_categories([first, *[level for level in values.cat.categories if level != first]])
      return ordered, f"{levels[1 - levels.index(reference)]} vs {reference}"
  ```
  `src/biotapy/da/_linda.py`:
  ```python
  """LinDA: linear models on log2 centred log-ratios, corrected for compositional bias."""

  from collections.abc import Sequence

  import numpy as np
  import numpy.typing as npt
  import pandas as pd
  from anndata import AnnData
  from scipy.stats import t as t_dist

  from ._design import dense_counts, design_matrix, model
  from ._schema import result


  def linda(adata: AnnData, group: str, *, covariates: Sequence[str] = (), reference: str | None = None) -> pd.DataFrame:
      """Differential abundance of each feature between two groups by LinDA.

      Parameters
      ----------
      adata
          Samples x features; ``X`` holds raw counts.
      group
          The ``obs`` column whose effect is reported: a categorical, string or bool
          column with two levels, or a numeric column.
      covariates
          ``obs`` columns to adjust for: numeric ones scaled to unit variance, others as
          one indicator per level against their first level.
      reference
          The level of a categorical ``group`` that the other level is compared
          with; by default its first category (sorted values for a string column).

      Returns
      -------
      pandas.DataFrame
          One row per feature, in ``var_names`` order, indexed by ``feature``:
          ``effect`` (log2 fold change, bias-corrected), ``se``, ``pvalue``,
          ``qvalue`` (Benjamini-Hochberg), ``direction`` (sign of ``effect``),
          ``method`` (``"linda"``) and ``contrast`` (``"<level> vs <reference>"``, or
          the column name for a numeric ``group``).

      Raises
      ------
      KeyError
          ``group`` or a covariate is not an ``obs`` column.
      TypeError
          ``covariates`` is not a list of column names, or ``reference`` is not a string.
      ValueError
          ``X`` does not hold raw counts, has an empty sample or fewer than two
          features; a used ``obs`` column has missing values or is constant;
          ``group`` has other than two levels; ``reference`` is not one of them or is
          given for a numeric ``group``; the model has at least as many terms as
          samples, or collinear columns.

      Notes
      -----
      R equivalent: ``MicrobiomeStat::linda``
      Guide: :doc:`/guide/differential_abundance`

      Equals ``MicrobiomeStat::linda(t(counts), meta, "~group + covariates",
      feature.dat.type = "count", is.winsor = FALSE)`` with fixed effects. When ``X``
      holds a zero, 0.5 is added to every count; the log2 counts are centred per
      sample; each feature gets one least-squares fit; the effect's bias is the mode
      of all features' effects, found by a Gaussian mean shift started at the
      shortest half of the data with ``bw.nrd0``'s bandwidth (``modeest::mlv``), and
      is subtracted; p-values are two-sided t-tests on ``n - p`` degrees of freedom.
      Numeric columns are scaled to unit variance first, so a numeric ``group``'s
      effect is per standard deviation.

      MicrobiomeStat 1.4 announces an imputation approach for zeros when library
      size depends on the model, but adds the pseudocount in every case (its switch
      compares ``"Imputation"`` with ``"imputation"``); biotapy computes what it
      returns. Not ported: winsorisation (MicrobiomeStat's default ``is.winsor =
      TRUE``), its prevalence and abundance filters (filter once with
      :func:`biotapy.pp.filter_features` before any method) and random effects.

      The log-ratios have no zeros, so ``X`` is densified once: 8 bytes x samples x
      features. Peak memory is about five such arrays (5.0x measured on 400 x 1,000
      and 200 x 2,000 tables), so budget for that on large tables.

      References
      ----------
      Zhou H, He K, Chen J, Zhang X (2022) LinDA: linear models for differential abundance
      analysis of microbiome compositional data. Genome Biology 23:95.

      Examples
      --------
      >>> import biotapy as bt
      >>> table = bt.da.linda(bt.datasets.toy(), "group")
      >>> round(float(table.loc["f6", "effect"]), 2), table.loc["f6", "contrast"]
      (6.3, 'B vs A')
      """
      frame, contrast = model(adata, group, covariates=covariates, reference=reference, func="da.linda")
      design = design_matrix(frame, scale=True)
      counts = dense_counts(adata, func="da.linda")
      if (counts == 0).any():
          counts += 0.5
      logs = np.log2(counts)
      ratios = logs - logs.mean(axis=1, keepdims=True)
      beta, se, dof = _fit(design, ratios)
      root_n = np.sqrt(adata.n_obs)
      effect = beta - _mode(root_n * beta) / root_n
      pvalue = 2 * t_dist.sf(np.abs(effect / se), dof)
      return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="linda", contrast=contrast)


  def _fit(
      design: npt.NDArray[np.float64], ratios: npt.NDArray[np.float64]
  ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], int]:
      """Least squares of every feature on ``design``: the group coefficient, its standard error, the residual df."""
      coef, *_ = np.linalg.lstsq(design, ratios, rcond=None)
      dof = design.shape[0] - design.shape[1]
      sigma2 = ((ratios - design @ coef) ** 2).sum(axis=0) / dof
      se = np.sqrt(sigma2 * np.linalg.inv(design.T @ design)[1, 1])
      return coef[1], se, dof


  def _mode(values: npt.NDArray[np.float64]) -> float:
      """``modeest::mlv(values, method = "meanshift", kernel = "gaussian")``, as modeest 2.4.0 computes it."""
      sd = values.std(ddof=1)
      iqr = np.subtract(*np.percentile(values, [75, 25]))
      # stats::bw.nrd0, including its fallbacks for a zero spread.
      spread = min(sd, iqr / 1.34) or sd or abs(values[0]) or 1.0
      bandwidth = 0.9 * spread * values.size**-0.2
      mode = _shorth(values)
      for _ in range(1000):
          weights = np.exp(-0.5 * ((values - mode) / bandwidth) ** 2)
          shifted = float(values @ weights / weights.sum())
          # modeest stops on a relative step below sqrt(eps) and returns the previous value, not the shifted one.
          # (It divides by a mode of exactly 0 and fails; an unchanged mode stops it either way.)
          if shifted == mode or (mode != 0 and abs(shifted / mode - 1) < np.sqrt(np.finfo(np.float64).eps)):
              break
          mode = shifted
      return mode


  def _shorth(values: npt.NDArray[np.float64]) -> float:
      """``modeest::shorth``: the mean of the shortest window holding half of the sorted values."""
      ordered = np.sort(values)
      k = int(np.ceil(ordered.size / 2)) - 1
      widths = ordered[k:] - ordered[: ordered.size - k]
      ties = np.flatnonzero(widths == widths.min()) + 1
      # Tied windows: modeest takes the mean of their 1-based starts, and R's indexing truncates the fraction.
      start = int(ties.mean()) - 1
      return float(ordered[start : start + k + 1].mean())
  ```
  `src/biotapy/__init__.py`:
  ```diff
  @@ -1,7 +1,7 @@
   from importlib.metadata import version

  -from . import datasets, fn, io, pl, pp, tl
  +from . import da, datasets, fn, io, pl, pp, tl

  -__all__ = ["__version__", "datasets", "fn", "io", "pl", "pp", "tl"]
  +__all__ = ["__version__", "da", "datasets", "fn", "io", "pl", "pp", "tl"]

   __version__ = version("biotapy")
  ```
- [x] **Step 9: Run, expect pass** - the same two commands -> `81 passed, 2 deselected`;
  `2 passed`. The property test also under `--hypothesis-seed=1`, `2`, `3`.
- [x] **Step 10: Docs.** Create `docs/guide/differential_abundance.md`:
  ````markdown
  # Differential abundance

  Differential abundance (DA) asks which features are more abundant in one group of samples than in
  another. Every `bt.da` method takes the same arguments and returns the same table, so their answers
  can be put side by side.

  ## One question, one table

  A method takes the `obs` column whose effect you want (`group`), the `covariates` to adjust for,
  and the `reference` level that the other level is compared with:

  ```python
  import biotapy as bt

  tdata = bt.datasets.toy()
  table = bt.da.linda(tdata, "group")  # "B vs A": A, the first category, is the reference
  table = bt.da.linda(tdata, "group", reference="B")  # "A vs B": every effect changes sign
  ```

  `group` is a categorical, string or bool column with two levels, or a numeric column, whose effect
  is then a slope. Covariates are numeric, or categorical with one indicator per level against
  their first. There are no formula strings: R and patsy take the alphabetically first level as the
  reference and silently drop samples with a missing value. biotapy takes the reference you give, or
  the first category, says which way round the effect is in the `contrast` column, and raises on a
  missing value, on a group with more than two levels and on collinear covariates.

  The table has one row per feature, in `var_names` order:

  | Column | Meaning |
  |---|---|
  | `effect` | log2 fold change of the other level over the reference, or the slope of a numeric group |
  | `se` | its standard error |
  | `pvalue` | the method's p-value |
  | `qvalue` | Benjamini-Hochberg adjusted p-value, over the features the method tested |
  | `direction` | sign of `effect`: -1, 0 or 1 |
  | `method` | the method's name |
  | `contrast` | `"<level> vs <reference>"`, or the name of a numeric group |

  A feature a method cannot test keeps its row, with NaN values and `direction` 0. Methods never
  filter features: filter once, with `bt.pp.filter_features`, before running any method, so every
  method tests the same features. Methods need raw counts in `X` and a read in every sample
  (`bt.pp.filter_samples(tdata, min_depth=1)` drops empty ones).

  ## LinDA

  `bt.da.linda` fits one linear model per feature to the log2 centred log-ratios of the counts (with
  0.5 added to every count when the table has a zero). Compositionality biases every coefficient by
  the same amount; LinDA estimates that bias as the mode of all features' coefficients and subtracts
  it, so most features end up near zero and the ones that changed stand out:

  ```python
  tdata.obs["age"] = [30, 41, 52, 38, 45, 60]
  table = bt.da.linda(tdata, "group", covariates=["age"])
  table[table["qvalue"] < 0.05]
  ```

  It equals `MicrobiomeStat::linda(..., is.winsor = FALSE)` in R with fixed effects. Numeric columns
  are scaled to unit variance first, as LinDA does, so a numeric group's effect is per standard
  deviation. Not available: winsorisation (MicrobiomeStat's default), random effects such as
  `(1 | subject)`, and LinDA's own prevalence filters.
  ````
  and list it in the guide, the API and the contributing page:
  ````diff
  @@ -14,5 +14,6 @@ function
   filtering
   diversity
   ordination
  +differential_abundance
   plotting
   ```
  ````
  ````diff
  @@ -91,6 +91,18 @@ Public functions are listed here as they ship, from Phase 1 onward.
       tl.unifrac
   ```

  +## Differential abundance
  +
  +```{eval-rst}
  +.. module:: biotapy.da
  +.. currentmodule:: biotapy
  +
  +.. autosummary::
  +    :toctree: generated
  +
  +    da.linda
  +```
  +
   ## Plots

   ```{eval-rst}
  ````
  ```diff
  @@ -59,7 +59,7 @@ These tests download `bt.datasets.global_patterns()`, `bt.datasets.enterotype()`
   and `bt.datasets.esophagus()` through [pooch](https://www.fatiando.org/pooch/)
   and compare biotapy with R on that data: `pp.relative`, `pp.clr`, `pp.tax_glom`, filtering, rarefaction
   (its invariants), alpha and beta diversity, UniFrac, PCoA, NMDS and PERMANOVA against phyloseq,
  -vegan, ape and picante, and `pp.philr` against philr. Set
  +vegan, ape and picante, `pp.philr` against philr, and `da.linda` against MicrobiomeStat. Set
   `BIOTAPY_DATA_DIR` to point the pooch cache somewhere other than the default
   per-user cache directory - CI caches it across runs the same way:

  ```
- [x] **Step 11: Contracts.** `.knowledge/contracts/data-model-slots.md`
  (frontmatter `paths` and a new "DA results" section before "Propagation"):
  ```diff
  @@ -4,7 +4,7 @@ title: Data-model slots
   description: Which AnnData/TreeData slot holds what, the exact result keys, the x_kind and provenance conventions, and which slots feature-changing operations drop.
   tags: [data-model, api]
   status: stable
  -paths: ["src/biotapy/_core/**", "src/biotapy/io/**", "src/biotapy/pp/**", "src/biotapy/tl/**", "src/biotapy/fn/**"]
  +paths: ["src/biotapy/_core/**", "src/biotapy/io/**", "src/biotapy/pp/**", "src/biotapy/tl/**", "src/biotapy/fn/**", "src/biotapy/da/**"]
   generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T15:43:25Z }
   commit: 6ade269
   sources:
  @@ -131,6 +131,26 @@ columns are always the seven, `kingdom` to `species`, even for a genus-level tab
   profile (the missing ranks are NaN, `io/_metaphlan.py:read_metaphlan`). `UNCLASSIFIED` (`UNKNOWN` in older
   tables) stays a feature with every rank NaN, so every non-empty sample sums to 1.

  +## DA results
  +`da` methods write no slot: each returns one `pd.DataFrame` built by
  +`da/_schema.py:result`, indexed by `var_names` (index name `feature`), one row
  +per feature in `var_names` order, never a filtered subset:
  +
  +| Column | dtype | Meaning |
  +|---|---|---|
  +| `effect` | float64 | log2 fold change of `group`'s other level over `reference`, or the slope of a numeric `group` (per standard deviation in `da.linda`, which scales numeric columns) |
  +| `se` | float64 | standard error of `effect` |
  +| `pvalue` | float64 | the method's p-value |
  +| `qvalue` | float64 | Benjamini-Hochberg over the finite p-values (`scipy.stats.false_discovery_control`) |
  +| `direction` | int8 | sign of `effect`; 0 when `effect` is NaN |
  +| `method` | str | the function's name, e.g. `"linda"` |
  +| `contrast` | str | `"<level> vs <reference>"`, or the column name of a numeric `group` |
  +
  +A feature the method cannot test keeps its row with NaN `effect`, `se`,
  +`pvalue` and `qvalue`. `group`, `covariates` and `reference` are checked once
  +for every method by `da/_design.py:model`: no formula strings, missing values
  +raise, a categorical `group` has exactly two levels.
  +
   ## Propagation
   | Operation | Keeps | Drops |
   |---|---|---|
  ```
  `.knowledge/contracts/r-golden-parity.md` statement 4, a row under "DA
  methods":
  ```diff
  @@ -48,6 +48,7 @@ sources:
      | Rarefaction | invariants only: row sums == depth, dropped samples identical, no count exceeds original | exact |
      | HUMAnN parity (func_glom, renorm) | elementwise, matched by row id | rtol=1e-7; renorm rtol=5e-6, because humann_renorm_table prints %.6g |
      | DA methods | sign agreement and rank correlation of effect sizes; exact match only where the R method is deterministic | per method |
  +   | `da.linda` vs `MicrobiomeStat::linda(is.winsor = FALSE)` (deterministic) | `effect`, `se`, `pvalue`, `qvalue` elementwise, matched by taxon | `rtol=1e-7` |

   5. Any looser tolerance is written in the test with a one-line comment giving the reason.
   6. Golden files hold numbers derived from third-party example data, never the
  ```
  Log line (under the slice 3B heading):
  ```text
  - **Update**: [data-model-slots](contracts/data-model-slots.md) gains "DA results": the schema `da` methods return, written by no slot; `paths` gains `src/biotapy/da/**`. [r-golden-parity](contracts/r-golden-parity.md): `da.linda` is compared elementwise with `MicrobiomeStat::linda`. [phase-3-stats](roadmap/phase-3-stats.md) tasks 3.3 and 3.5 done.
  ```
- [x] **Step 12: Gate and commit**
  ```bash
  git add src/biotapy/da/__init__.py src/biotapy/da/_schema.py src/biotapy/da/_design.py src/biotapy/da/_linda.py \
    src/biotapy/__init__.py tests/da/conftest.py tests/da/test_linda.py tests/da/test_linda_golden.py \
    tests/da/test_design.py tests/da/test_schema.py tests/test_docstrings.py docs/guide/differential_abundance.md \
    docs/guide/index.md docs/api.md docs/contributing.md .knowledge/contracts/data-model-slots.md \
    .knowledge/contracts/r-golden-parity.md .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  uvx prek run --all-files   # after staging: --all-files skips untracked files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git commit -m "feat(da): add LinDA and the differential abundance result schema"
  ```

### Task 3.4: `da.ancombc2`

**Files:** modify `tests/r/Dockerfile`, `tests/r/export_golden.R`,
`tests/golden/VERSIONS.txt`, `src/biotapy/da/__init__.py`,
`tests/da/test_design.py`, `tests/da/test_schema.py`, `pyproject.toml`,
`docs/guide/differential_abundance.md`, `docs/api.md`, `docs/contributing.md`,
`.knowledge/contracts/r-golden-parity.md`,
`.knowledge/playbooks/regenerate-golden-files.md`; create
`tests/golden/global_patterns/ancombc2.csv.gz`, `src/biotapy/da/_ancombc.py`,
`tests/da/test_ancombc.py`, `tests/da/test_ancombc_golden.py`.
**Not touched:** `da/_design.py`, `da/_schema.py`, `da/_linda.py`.
**Interfaces:**
- Consumes: `model`, `dense_counts`, `result` (3.5); the `benchmark` fixture;
  `skbio.stats.composition.ancombc2`.
- Produces: `bt.da.ancombc2(adata: AnnData, group: str, *, covariates:
  Sequence[str] = (), reference: str | None = None) -> pd.DataFrame`;
  `ancombc2.csv.gz`.

- [x] **Step 1: Approval on record.** ANCOMBC (Bioconductor, Artistic-2.0) in
  the image was approved with the plan (Phase 3 decision 6). New here and to confirm:
  its install needs CVXR pinned to 1.0-15 (CRAN archive, Apache-2.0, built from
  source in the image) and apt `libgsl27` (slice 3B decision 5); both live only
  in the image.
- [x] **Step 2: Add ANCOMBC to the image.** `tests/r/Dockerfile`:
  ```diff
  @@ -20,6 +20,14 @@ RUN apt-get update && apt-get install -y --no-install-recommends libuv1 && rm -r
   RUN Rscript -e 'BiocManager::install("philr", version = "3.22", ask = FALSE, update = FALSE)'
   # MicrobiomeStat (CRAN, GPL-3): the da.linda golden file; biotapy reimplements linda and never calls it.
   RUN Rscript -e 'install.packages("MicrobiomeStat")'
  -RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE))'
  +# ANCOMBC (Bioconductor, Artistic-2.0): the da.ancombc2 golden file; biotapy calls scikit-bio's ancombc2 instead.
  +# ANCOMBC 2.12.0 imports CVXR::solve, which CVXR 1.8 (the CRAN snapshot's version) removed: pin CVXR 1.0-15, the last
  +# release before it, built from the CRAN archive after its dependencies. The gsl binary among ANCOMBC's dependencies
  +# links libgsl at load time.
  +RUN apt-get update && apt-get install -y --no-install-recommends libgsl27 && rm -rf /var/lib/apt/lists/*
  +RUN Rscript -e 'install.packages(c("Rmpfr", "gmp", "ECOSolveR", "scs", "osqp", "bit64", "cli", "Rcpp", "RcppEigen"))' \
  +    && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/CVXR/CVXR_1.0-15.tar.gz", repos = NULL, type = "source")'
  +RUN Rscript -e 'BiocManager::install("ANCOMBC", version = "3.22", ask = FALSE, update = FALSE)'
  +RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE), requireNamespace("ANCOMBC", quietly = TRUE))'
   WORKDIR /work
   CMD ["Rscript", "tests/r/export_golden.R"]
  ```
  `.knowledge/contracts/r-golden-parity.md` statement 1:
  ```diff
  @@ -20,8 +20,10 @@ sources:
      installs only what the current golden files need (`phyloseq`, which brings
      `Biostrings`, `vegan` and `ape`, plus CRAN `picante` for Faith PD and
      Bioconductor `philr` for `pp.philr`, with the `libuv1` runtime library its
  -   `fs` binary loads, and CRAN `MicrobiomeStat` for `da.linda`); a new golden
  -   function that needs another package adds it in its own commit (rules.md R2.3).
  +   `fs` binary loads, CRAN `MicrobiomeStat` for `da.linda`, and Bioconductor
  +   `ANCOMBC` for `da.ancombc2`, with CRAN's archived CVXR 1.0-15 and the
  +   `libgsl27` runtime library it needs); a new golden function that needs another
  +   package adds it in its own commit (rules.md R2.3).
   1b. HUMAnN golden files (`fn.func_glom`, `fn.renorm`) are produced by
      `tests/humann/export_golden.py`, run with
      `uv run --no-project --with humann==3.9 --with pandas==3.0.6`, never
  ```
  `.knowledge/playbooks/regenerate-golden-files.md`, Common mistakes:
  ```diff
  @@ -86,6 +86,12 @@ uv run --group test pytest tests/test_data_files.py -q
     installs `fs` as a binary that links libuv at load time, so the image needs
     the system package (`tests/r/Dockerfile`). `philr` is installed only for
     `pp.philr`'s golden file; biotapy never calls it.
  +- Installing `ANCOMBC` without pinning CVXR: ANCOMBC 2.12.0 (Bioconductor
  +  3.22) imports `CVXR::solve`, which CVXR 1.8 (the image's CRAN snapshot)
  +  removed, so ANCOMBC's lazy load fails. The Dockerfile builds CVXR 1.0-15 from
  +  the CRAN archive first, and installs `libgsl27`, which the `gsl` binary in
  +  ANCOMBC's dependencies links. `ANCOMBC` is installed only for `da.ancombc2`'s
  +  golden file; biotapy calls scikit-bio's `ancombc2`.
   - Calling `distance()` unqualified: Biostrings, attached after phyloseq, masks
     it with IRanges' generic, so write `phyloseq::distance`.
   - Sharing one `set.seed(20260927)` across sections: every random call needs
  ```
  Log line:
  ```text
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `ANCOMBC` for the `da.ancombc2` golden file, with CVXR 1.0-15 and `libgsl27`. [regenerate-golden-files](playbooks/regenerate-golden-files.md) Common mistakes: why CVXR is pinned.
  ```
  Build: `docker build -t biotapy-golden tests/r`. Expected: the guard layer
  passes; `ANCOMBC 2.12.0`, `CVXR 1.0.15`, `lme4 2.0.1`; every earlier
  package unchanged. Without the CVXR lines the guard fails after
  `Error: object 'solve' is not exported by 'namespace:CVXR'`; without
  `libgsl27`, after `libgsl.so.27: cannot open shared object file`.
- [x] **Step 3: Commit the image change on its own:**
  ```bash
  git add tests/r/Dockerfile .knowledge/contracts/r-golden-parity.md .knowledge/playbooks/regenerate-golden-files.md \
    .knowledge/log.md
  git commit -m "build(r): add ANCOMBC to the golden image"
  ```
- [x] **Step 4: Export.** In `tests/r/export_golden.R`, after the LinDA
  section, and `VERSIONS.txt` gains ANCOMBC:
  ```diff
  @@ -159,6 +159,20 @@ linda_rows <- lapply(da_formulas, function(f) {
   })
   write_golden(do.call(rbind, linda_rows), file.path(gp, "linda.csv.gz"))

  +# The settings scikit-bio's ancombc2 mirrors: BH, no prevalence or library-size filter, no pseudocount sensitivity
  +# analysis, no structural-zero test; pseudo = 0, s0_perc, iter_control and em_control keep R's defaults. The bias E-M
  +# runs in a %dorng% loop, which takes its seeds from R's generator: hence the seed, though no step draws a number.
  +set.seed(20260927)
  +ancombc_rows <- lapply(da_formulas, function(f) {
  +  out <- suppressMessages(ANCOMBC::ancombc2(
  +    data = da_counts, taxa_are_rows = TRUE, meta_data = da_meta, fix_formula = f, p_adj_method = "BH",
  +    prv_cut = 0, lib_cut = 0, pseudo_sens = FALSE, struc_zero = FALSE, verbose = FALSE
  +  ))$res
  +  data.frame(formula = f, taxon_id = out$taxon, lfc = out$lfc_hosthuman, se = out$se_hosthuman,
  +             p = out$p_hosthuman, q = out$q_hosthuman)
  +})
  +write_golden(do.call(rbind, ancombc_rows), file.path(gp, "ancombc2.csv.gz"))
  +
   ## Synthetic phyloseq fixtures: biotapy's toy() numbers, no third-party data
   counts <- rbind(
     c(10, 5, 20, 30, 0, 2, 1, 0), c(8, 7, 25, 22, 3, 0, 0, 1), c(12, 4, 18, 35, 1, 5, 2, 0),
  @@ -222,5 +236,6 @@ writeLines(c(
     paste0("picante ", packageVersion("picante")),
     paste0("philr ", packageVersion("philr")),
     paste0("MicrobiomeStat ", packageVersion("MicrobiomeStat")),
  -  paste0("modeest ", packageVersion("modeest"))
  +  paste0("modeest ", packageVersion("modeest")),
  +  paste0("ANCOMBC ", packageVersion("ANCOMBC")),
  +  paste0("CVXR ", packageVersion("CVXR"), " (CRAN archive, pinned in tests/r/Dockerfile)")
   ), "tests/golden/VERSIONS.txt")
  ```
  Run twice as in 3.5 Step 4. Expected: every line `OK`; `git status` shows
  `M tests/golden/VERSIONS.txt` (two new lines, `ANCOMBC 2.12.0` and `CVXR 1.0.15 (CRAN archive...)`),
  `M tests/r/export_golden.R` and the new `ancombc2.csv.gz` (45.8 KB);
  `linda.csv.gz` and every older file byte-identical. `meta_data` needs both
  columns: ANCOMBC 2.12.0's `data_sanity_check` breaks a one-column
  `meta_data` (`invalid factor level, NA generated`).
- [x] **Step 5: Gate and commit.** `uv run --group test pytest
  tests/test_data_files.py -q` -> passes; `uvx prek run --all-files`;
  ```bash
  git add tests/r/export_golden.R tests/golden/VERSIONS.txt tests/golden/global_patterns/ancombc2.csv.gz
  git commit -m "test(golden): export the ANCOM-BC2 golden file"
  ```
- [x] **Step 6: Failing tests.** `tests/da/test_ancombc.py`:
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from scipy.stats import false_discovery_control
  from skbio.stats.composition import ancombc2

  import biotapy as bt


  def _toy_without(feature, samples):
      tdata = bt.datasets.toy()
      dense = tdata.X.toarray()
      dense[samples, tdata.var_names.get_loc(feature)] = 0
      tdata.X = sp.csr_matrix(dense)
      return tdata


  def test_ancombc2_is_scikit_bio_in_log2():
      tdata = bt.datasets.toy()
      out = bt.da.ancombc2(tdata, "group")
      counts = pd.DataFrame(tdata.X.toarray(), index=tdata.obs_names, columns=tdata.var_names)
      expected = ancombc2(counts, tdata.obs, "group").result.xs("group[T.B]", level="Covariate")
      np.testing.assert_allclose(out["effect"], expected["Log(FC)"] / np.log(2), rtol=1e-12)
      np.testing.assert_allclose(out["se"], expected["SE"] / np.log(2), rtol=1e-12)
      np.testing.assert_allclose(out["pvalue"], expected["pvalue"], rtol=1e-12)


  def test_feature_absent_from_a_group_is_not_tested():
      out = bt.da.ancombc2(_toy_without("f5", [0, 1, 2]), "group")  # no f5 read in group A
      assert out.loc["f5", ["effect", "se", "pvalue", "qvalue"]].isna().all() and out.loc["f5", "direction"] == 0
      tested = out.drop(index="f5")
      assert np.isfinite(tested[["effect", "se", "pvalue", "qvalue"]].to_numpy()).all()
      # scikit-bio and R count it with p = 1; biotapy leaves it out of the correction.
      np.testing.assert_allclose(tested["qvalue"], false_discovery_control(tested["pvalue"]), rtol=1e-12)


  def test_all_zero_feature_is_not_tested():
      out = bt.da.ancombc2(_toy_without("f7", slice(None)), "group")
      assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
      assert out.loc["f7", ["effect", "pvalue"]].isna().all() and out["effect"].notna().sum() == 7


  def test_ancombc2_keeps_input(assert_unchanged):
      tdata = bt.datasets.toy()
      before = tdata.copy()
      bt.da.ancombc2(tdata, "group")
      assert_unchanged(before, tdata)


  def test_feature_with_no_residual_degrees_of_freedom_is_not_tested():
      tdata = bt.datasets.toy()
      dense = tdata.X.toarray()
      dense[:, tdata.var_names.get_loc("f8")] = [0, 4, 0, 0, 0, 3]  # two reads in all: as many as model terms
      tdata.X = sp.csr_matrix(dense)
      out = bt.da.ancombc2(tdata, "group")
      assert out.loc["f8", ["effect", "se", "pvalue", "qvalue"]].isna().all() and out.loc["f8", "direction"] == 0
      assert out.drop(index="f8")["effect"].notna().all()


  def test_scikit_bio_failure_names_the_function():
      tdata = bt.datasets.toy()[:, ["f1", "f2"]].copy()
      tdata.X = sp.csr_matrix(
          np.array([[5, 0], [6, 0], [7, 0], [0, 5], [0, 6], [0, 7]])
      )  # each feature in one group only
      with pytest.raises(ValueError, match=r"da\.ancombc2: scikit-bio could not fit the model: .*estimable"):
          bt.da.ancombc2(tdata, "group")
  ```
  `tests/da/test_ancombc_golden.py`:
  ```python
  from pathlib import Path

  import numpy as np
  import pandas as pd
  import pytest
  from scipy.stats import false_discovery_control, spearmanr

  import biotapy as bt

  GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
  pytestmark = [pytest.mark.golden, pytest.mark.network]


  # Per model: (formula, covariates, (effect atol in log2, se rtol, pvalue atol, Spearman floor)).
  # host: the bias E-M stops at R's 100-iteration cap before converging, and scikit-bio stops on a slightly different
  # iterate, so every effect is shifted by 0.002-0.012 log2 (measured max 0.0121), se by 9.0e-4 relative and p by
  # 0.0122; the shift is near constant, hence the rank floor. With log_depth the E-M converges on both sides and
  # agreement is 1.5e-8 (effect), 1.7e-10 (se, relative) and 1.1e-8 (p), checked at 1e-6.
  MODELS = [
      pytest.param(("host", (), (0.015, 2e-3, 0.02, 0.9999)), id="host"),
      pytest.param(("host + log_depth", ("log_depth",), (1e-6, 1e-6, 1e-6, 0.999999)), id="host+log_depth"),
  ]


  @pytest.mark.parametrize("case", MODELS)
  def test_ancombc2_matches_ancombc_ancombc2(benchmark, case):
      # R: ANCOMBC::ancombc2(counts, meta_data = meta, fix_formula, p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
      # pseudo_sens = FALSE, struc_zero = FALSE)$res, natural logs, on the 636 genera.
      formula, covariates, (effect_atol, se_rtol, p_atol, rank_floor) = case
      golden = pd.read_csv(GOLDEN / "ancombc2.csv.gz", dtype={"taxon_id": str}).query("formula == @formula")
      golden = golden.set_index("taxon_id")
      out = bt.da.ancombc2(benchmark, "host", covariates=covariates, reference="other")
      assert out.index.tolist() == golden.index.tolist()
      # The 36 genera with no read in one host group: R reports NA with p = 1, biotapy NaN.
      untested = golden["lfc"].isna()
      assert untested.sum() == 36 and (golden.loc[untested, "p"] == 1).all()
      assert out["effect"].isna().equals(untested.rename("effect"))
      ours, theirs = out[~untested], golden[~untested]
      np.testing.assert_allclose(ours["effect"], theirs["lfc"] / np.log(2), atol=effect_atol)
      np.testing.assert_allclose(ours["se"], theirs["se"] / np.log(2), rtol=se_rtol)
      np.testing.assert_allclose(ours["pvalue"], theirs["p"], atol=p_atol)
      assert spearmanr(ours["effect"], theirs["lfc"]).statistic > rank_floor
      # R's q counts the untested genera with p = 1; the schema's BH runs over the tested ones only.
      expected = false_discovery_control(theirs["p"])
      np.testing.assert_array_equal(ours["qvalue"] < 0.05, expected < 0.05)
  ```
  `tests/da/test_design.py` runs every shared test on both methods; ANCOM-BC2's
  bias E-M is not exactly antisymmetric, so the reference swap gets a
  per-method tolerance:
  ```diff
  @@ -5,7 +5,9 @@ import scipy.sparse as sp

   import biotapy as bt

  -METHODS = [bt.da.linda]
  +METHODS = [bt.da.ancombc2, bt.da.linda]
  +# ANCOM-BC2's bias E-M starts from asymmetric quantiles, so swapping the reference flips its effects only to ~1e-3.
  +SWAP_ATOL = {"ancombc2": 2e-3, "linda": 1e-12}


   def _toy_with(**columns):
  @@ -20,9 +22,11 @@ def test_reference_sets_the_sign(method):
       tdata = bt.datasets.toy()
       a_first, b_first = method(tdata, "group"), method(tdata, "group", reference="B")
       assert (a_first["contrast"] == "B vs A").all() and (b_first["contrast"] == "A vs B").all()
  -    np.testing.assert_allclose(b_first["effect"], -a_first["effect"], rtol=1e-9, atol=1e-12)
  -    assert (b_first["direction"] == -a_first["direction"]).all()
  -    np.testing.assert_allclose(b_first["pvalue"], a_first["pvalue"], rtol=1e-9)
  +    atol = SWAP_ATOL[method.__name__]
  +    np.testing.assert_allclose(b_first["effect"], -a_first["effect"], rtol=1e-9, atol=atol)
  +    clear = a_first["effect"].abs() > atol
  +    assert (b_first.loc[clear, "direction"] == -a_first.loc[clear, "direction"]).all() and clear.sum() >= 7
  +    np.testing.assert_allclose(b_first["pvalue"], a_first["pvalue"], rtol=1e-9, atol=atol)


   @pytest.mark.parametrize("method", METHODS)
  ```
  `tests/da/test_schema.py`: both methods, and distinct counts in the
  property test (degenerate tables make scikit-bio's E-M fail with
  `ValueError: Quantiles must be in the range [0, 1]`, as R stops with "Zero
  variances have been detected"):
  ```diff
  @@ -10,7 +10,7 @@ from scipy.stats import false_discovery_control

   import biotapy as bt

  -METHODS = [bt.da.linda]
  +METHODS = [bt.da.ancombc2, bt.da.linda]
   COLUMNS = ["effect", "se", "pvalue", "qvalue", "direction", "method", "contrast"]


  @@ -36,11 +36,11 @@ def test_qvalue_is_benjamini_hochberg(method):
   @pytest.mark.parametrize("method", METHODS)
   @settings(max_examples=25, deadline=None)
   @given(
  -    counts=arrays(np.int64, (6, 7), elements=st.integers(0, 50)),
  +    # Distinct counts: no empty sample and no feature fitted exactly, where ANCOM-BC2's variances are 0 and it fails.
  +    counts=arrays(np.int64, (6, 7), elements=st.integers(0, 200), unique=True),
       split=st.integers(2, 4),
   )
   def test_direction_is_the_sign_and_qvalue_bounds_pvalue(method, counts, split):
  -    counts[:, 0] += 1  # no empty sample
       adata = ad.AnnData(
           X=sp.csr_matrix(counts),
           obs=pd.DataFrame({"g": ["a"] * split + ["b"] * (6 - split)}, index=[f"s{i}" for i in range(6)]),
  ```
- [x] **Step 7: Run, expect failure** -
  `uv run --group test pytest tests/da -q --continue-on-collection-errors` ->
  `4 failed, 7 passed, 4 deselected, 2 errors` (`AttributeError: module
  'biotapy.da' has no attribute 'ancombc2'`; `test_design.py` and
  `test_schema.py` fail to collect); `uv run --group test pytest
  tests/da/test_ancombc_golden.py -q -m "golden or network"` -> `2 failed`.
- [x] **Step 8: Implement.** `src/biotapy/da/_ancombc.py`:
  ```python
  """ANCOM-BC2 through scikit-bio: bias-corrected log abundances, one linear model per feature."""

  from collections.abc import Sequence

  import numpy as np
  import pandas as pd
  from anndata import AnnData
  from skbio.stats.composition import ancombc2 as skbio_ancombc2

  from ._design import dense_counts, model
  from ._schema import result


  def ancombc2(
      adata: AnnData, group: str, *, covariates: Sequence[str] = (), reference: str | None = None
  ) -> pd.DataFrame:
      """Differential abundance of each feature between two groups by ANCOM-BC2.

      Parameters
      ----------
      adata
          Samples x features; ``X`` holds raw counts.
      group
          The ``obs`` column whose effect is reported: a categorical, string or bool
          column with two levels, or a numeric column, whose effect is per unit
          (``da.linda`` reports it per standard deviation, as MicrobiomeStat does).
      covariates
          ``obs`` columns to adjust for: numeric ones as they are, others as one
          indicator per level against their first level.
      reference
          The level of a categorical ``group`` that the other level is compared
          with; by default its first category (sorted values for a string column).

      Returns
      -------
      pandas.DataFrame
          One row per feature, in ``var_names`` order, indexed by ``feature``:
          ``effect`` (log2 fold change, bias-corrected), ``se``, ``pvalue``,
          ``qvalue`` (Benjamini-Hochberg over the tested features), ``direction``,
          ``method`` (``"ancombc2"``) and ``contrast``. A feature the model cannot
          test has NaN ``effect``, ``se``, ``pvalue`` and ``qvalue`` and
          ``direction`` 0.

      Raises
      ------
      KeyError
          ``group`` or a covariate is not an ``obs`` column.
      TypeError
          ``covariates`` is a string rather than a list of column names.
      ValueError
          ``X`` does not hold raw counts or has an empty sample; a used ``obs`` column
          has missing values; ``group`` has other than two levels; ``reference`` is
          not one of them or is given for a numeric ``group``; the model has as many
          terms as samples, or collinear columns.

      Notes
      -----
      R equivalent: ``ANCOMBC::ancombc2``
      Guide: :doc:`/guide/differential_abundance`

      Runs scikit-bio's :func:`~skbio.stats.composition.ancombc2` (Lin and Peddada
      2024), which matches ``ANCOMBC::ancombc2(..., fix_formula, p_adj_method =
      "BH", prv_cut = 0, pseudo_sens = FALSE)`` with its other defaults: zeros are
      treated as missing (``pseudo = 0``), ``s0_perc = 0.05``, no structural-zero
      test. ANCOM-BC2 reports natural logs; ``effect`` and ``se`` are divided by
      ln 2. A feature whose zeros leave one level of ``group`` without an observed
      value, or that has no more observed samples than model terms, cannot be
      tested: R and scikit-bio report it with p = 1 or NaN and R counts it in the
      correction, biotapy reports NaN for ``effect``, ``se``, ``pvalue`` and
      ``qvalue`` and leaves it out, so its q-values are BH over the features
      actually tested. scikit-bio rebases a categorical covariate with three or more
      levels when a feature has no read at its first level; only a covariate left
      with one observed level drops the feature.
      Not run: the pseudocount sensitivity analysis (R's default ``pseudo_sens =
      TRUE``), whose ``passed_ss`` flag the result table does not carry, and R's
      prevalence filter (``prv_cut = 0.10``; filter once with
      :func:`biotapy.pp.filter_features` before any method).

      The bias is estimated by an E-M algorithm capped at 100 iterations, R's
      ``em_control`` default, which biotapy keeps. On some data it has not
      converged by then, and the estimate depends on the cap: on GlobalPatterns'
      genera, human hosts against the rest, R moves every effect by about -0.28
      log2 with 1,000 iterations, and calls 220 genera instead of 208.

      scikit-bio needs a dense table, so ``X`` is densified once: 8 bytes x
      samples x features. Peak memory is about seven such arrays (7.2x to 7.5x
      measured on 400 x 1,000 and 200 x 2,000 tables).

      References
      ----------
      Lin H, Peddada SD (2024) Multigroup analysis of compositions of microbiomes with
      covariate adjustments and repeated measures. Nature Methods 21:83-91.

      Examples
      --------
      >>> import biotapy as bt
      >>> table = bt.da.ancombc2(bt.datasets.toy(), "group")
      >>> round(float(table.loc["f6", "effect"]), 2), table.loc["f6", "contrast"]
      (3.76, 'B vs A')
      """
      frame, contrast = model(adata, group, covariates=covariates, reference=reference, func="da.ancombc2")
      counts = pd.DataFrame(
          dense_counts(adata, func="da.ancombc2"), index=adata.obs_names, columns=adata.var_names, copy=False
      )
      # Plain names, so patsy needs no quoting; categorical columns keep their levels, reference first.
      names = [f"x{i}" for i in range(frame.shape[1])]
      try:
          fit = skbio_ancombc2(counts, frame.set_axis(names, axis=1), " + ".join(names), p_adjust=None).result
      except ValueError as err:
          msg = f"da.ancombc2: scikit-bio could not fit the model: {err}"
          raise ValueError(msg) from err
      covariate = next(name for name in fit.index.unique("Covariate") if name == "x0" or name.startswith("x0["))
      table = fit.xs(covariate, level="Covariate").reindex(adata.var_names)
      effect = table["Log(FC)"].to_numpy(np.float64) / np.log(2)
      se = table["SE"].to_numpy(np.float64) / np.log(2)
      pvalue = table["pvalue"].to_numpy(np.float64)
      # scikit-bio reports an unfitted feature with p = 1, and one with no residual degrees of freedom with a finite
      # effect and NaN p; the schema says "not tested" with NaN throughout.
      untested = ~(np.isfinite(effect) & np.isfinite(pvalue))
      effect, se, pvalue = (np.where(untested, np.nan, values) for values in (effect, se, pvalue))
      return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="ancombc2", contrast=contrast)
  ```
  `src/biotapy/da/__init__.py`:
  ```diff
  @@ -1,5 +1,6 @@
   """Differential abundance: methods that share one result table (contracts/data-model-slots, DA results)."""

  +from ._ancombc import ancombc2
   from ._linda import linda

  -__all__ = ["linda"]
  +__all__ = ["ancombc2", "linda"]
  ```
  `pyproject.toml` (`ancombc2` is unannotated; mypy strict reports `Call to
  untyped function "ancombc2" in typed context`):
  ```diff
  @@ -176,12 +176,13 @@ lint.pylint.max-statements = 30
   files = [ "docs/extensions", "src/biotapy" ]
   python_version = "3.12"
   # biom-format and scikit-learn ship no type annotations at all, nor do scikit-bio's
  -# subsample_counts and tree_basis, threadpoolctl's threadpool_limits and mudata's MuData.update; exempt only calls into
  -# them (not our own code) from strict's disallow_untyped_calls (mypy/checkexpr.py matches by callee fullname).
  +# subsample_counts, tree_basis and ancombc2, threadpoolctl's threadpool_limits and mudata's MuData.update; exempt only
  +# calls into them (not our own code) from strict's disallow_untyped_calls (mypy/checkexpr.py matches by callee fullname).
   untyped_calls_exclude = [
     "biom",
     "mudata",
     "skbio.stats._subsample",
  +  "skbio.stats.composition._ancombc.ancombc2",
     "skbio.stats.composition._base.tree_basis",
     "sklearn",
     "threadpoolctl",
  ```
- [x] **Step 9: Run, expect pass** - the same two commands -> `67 passed, 4 deselected` (16 more than the
  prototype: the eight shared model checks of 3.5's fix round run on both methods);
  `2 passed`. The property test also under three extra seeds.
- [x] **Step 10: Docs.** Append to `docs/guide/differential_abundance.md`:
  ````diff
  @@ -58,3 +58,23 @@ It equals `MicrobiomeStat::linda(..., is.winsor = FALSE)` in R with fixed effect
   are scaled to unit variance first, as LinDA does, so a numeric group's effect is per standard
   deviation. Not available: winsorisation (MicrobiomeStat's default), random effects such as
   `(1 | subject)`, and LinDA's own prevalence filters.
  +
  +## ANCOM-BC2
  +
  +`bt.da.ancombc2` runs scikit-bio's ANCOM-BC2 (Lin and Peddada 2024): it estimates each sample's
  +sampling fraction, corrects the log counts for it and the coefficients for their shared bias, and
  +fits one linear model per feature. Zeros are treated as missing rather than given a pseudocount,
  +so a feature with no read in one of the groups cannot be fitted: its row is NaN and it is left out
  +of the Benjamini-Hochberg correction (R's `ANCOMBC::ancombc2` reports it with p = 1 and counts it).
  +
  +```python
  +table = bt.da.ancombc2(tdata, "group")
  +```
  +
  +ANCOM-BC2 reports natural logs; biotapy divides `effect` and `se` by ln 2, so they are log2 like
  +every other method's; a numeric `group`'s effect is per unit, where `da.linda`'s is per standard
  +deviation. The settings are R's `ancombc2(..., p_adj_method = "BH", prv_cut = 0,
  +pseudo_sens = FALSE)`: biotapy does not run R's pseudocount sensitivity analysis (its
  +`passed_ss` flag) or its 10% prevalence filter. On the GlobalPatterns genera that the golden tests
  +use, biotapy's effects are within 0.012 log2 of R's and the significant genera are the same; the
  +small difference comes from the bias estimate, whose iterations stop at R's cap of 100 before they
  +have converged on that data, in R as in scikit-bio.
  ````
  ````diff
  @@ -100,6 +100,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
   .. autosummary::
       :toctree: generated

  +    da.ancombc2
       da.linda
   ```

  ````
  ```diff
  @@ -59,7 +59,7 @@ These tests download `bt.datasets.global_patterns()`, `bt.datasets.enterotype()`
   and `bt.datasets.esophagus()` through [pooch](https://www.fatiando.org/pooch/)
   and compare biotapy with R on that data: `pp.relative`, `pp.clr`, `pp.tax_glom`, filtering, rarefaction
   (its invariants), alpha and beta diversity, UniFrac, PCoA, NMDS and PERMANOVA against phyloseq,
  -vegan, ape and picante, `pp.philr` against philr, and `da.linda` against MicrobiomeStat. Set
  +vegan, ape and picante, `pp.philr` against philr, `da.linda` against MicrobiomeStat and `da.ancombc2` against ANCOMBC. Set
   `BIOTAPY_DATA_DIR` to point the pooch cache somewhere other than the default
   per-user cache directory - CI caches it across runs the same way:

  ```
- [x] **Step 11: Contracts.** `.knowledge/contracts/r-golden-parity.md`
  statement 4:
  ```diff
  @@ -51,6 +51,7 @@ sources:
      | HUMAnN parity (func_glom, renorm) | elementwise, matched by row id | rtol=1e-7; renorm rtol=5e-6, because humann_renorm_table prints %.6g |
      | DA methods | sign agreement and rank correlation of effect sizes; exact match only where the R method is deterministic | per method |
      | `da.linda` vs `MicrobiomeStat::linda(is.winsor = FALSE)` (deterministic) | `effect`, `se`, `pvalue`, `qvalue` elementwise, matched by taxon | `rtol=1e-7` |
  +   | `da.ancombc2` vs `ANCOMBC::ancombc2` (deterministic, but its bias E-M can stop at 100 iterations before converging, on a slightly different iterate in scikit-bio) | the same untested features; `effect`, `se`, `pvalue` elementwise; Spearman correlation of effects; the same calls at `q < 0.05`, with R's p-values corrected over the tested features | `host` model: `effect` atol 0.015 (log2), `se` rtol 2e-3, `pvalue` atol 0.02, Spearman > 0.9999; `host + log_depth`: all three at 1e-6, Spearman > 0.999999 |

   5. Any looser tolerance is written in the test with a one-line comment giving the reason.
   6. Golden files hold numbers derived from third-party example data, never the
  ```
  Log line:
  ```text
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md): `da.ancombc2` is compared with `ANCOMBC::ancombc2` at the measured tolerances. [phase-3-stats](roadmap/phase-3-stats.md) task 3.4 done.
  ```
- [x] **Step 12: Gate and commit**
  ```bash
  git add src/biotapy/da/_ancombc.py src/biotapy/da/__init__.py tests/da/test_ancombc.py tests/da/test_ancombc_golden.py \
    tests/da/test_design.py tests/da/test_schema.py pyproject.toml docs/guide/differential_abundance.md docs/api.md \
    docs/contributing.md .knowledge/contracts/r-golden-parity.md .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  uvx prek run --all-files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git commit -m "feat(da): add ANCOM-BC2 through scikit-bio"
  ```

### Task 3.8: `da.consensus` and the agreement decision

**Files:** modify `src/biotapy/da/_schema.py`, `src/biotapy/da/__init__.py`,
`tests/da/test_schema.py`, `docs/guide/differential_abundance.md`,
`docs/api.md`, `.knowledge/contracts/function-shape.md`, `rules.md`,
`.knowledge/decisions/index.md`; create `src/biotapy/da/_consensus.py`,
`tests/da/test_consensus.py`, `.knowledge/decisions/da-consensus-agreement.md`.
**Not touched:** the methods, `pl`.
**Interfaces:**
- Consumes: result tables from `bt.da.linda` and `bt.da.ancombc2`.
- Produces: `bt.da.consensus(results: Sequence[pd.DataFrame], *, alpha: float
  = 0.05, min_methods: int = 2) -> pd.DataFrame`; private `da._schema.COLUMNS`
  and `validate_result(table: object, *, arg: str) -> pd.DataFrame`; the
  Decision concept.

- [x] **Step 1: Failing tests.** `tests/da/test_consensus.py`:
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  from hypothesis import given, settings
  from hypothesis import strategies as st

  import biotapy as bt


  def _table(method, effect, qvalue, *, features=None, contrast="B vs A"):
      """A result table as a bt.da method writes it; NaN qvalue marks an untested feature."""
      effect, qvalue = np.asarray(effect, dtype=float), np.asarray(qvalue, dtype=float)
      features = features or [f"f{i}" for i in range(len(effect))]
      effect = np.where(np.isnan(qvalue), np.nan, effect)
      return pd.DataFrame(
          {
              "effect": effect,
              "se": np.where(np.isnan(qvalue), np.nan, 0.5),
              "pvalue": qvalue / 2,
              "qvalue": qvalue,
              "direction": np.sign(np.nan_to_num(effect)).astype(np.int8),
              "method": method,
              "contrast": contrast,
          },
          index=pd.Index(features, name="feature"),
      )


  def test_consensus_counts_the_methods_that_call_each_feature():
      a = _table("a", [2.0, -1.0, 0.5, 3.0], [0.01, 0.01, 0.30, 0.02])
      b = _table("b", [1.5, -2.0, 0.4, 2.0], [0.02, 0.20, 0.01, 0.03])
      out = bt.da.consensus([a, b])
      assert out["n_tested"].tolist() == [2, 2, 2, 2]
      assert out["n_significant"].tolist() == [2, 1, 1, 2]
      assert out["direction"].tolist() == [1, -1, 1, 1]
      assert out["consensus"].tolist() == [True, False, False, True]
      assert out["significant_a"].tolist() == [True, True, False, True]


  def test_columns_follow_the_results_order():
      out = bt.da.consensus([_table("b", [1.0], [0.01]), _table("a", [1.0], [0.01])])
      assert out.columns.tolist() == [
          *["effect_b", "qvalue_b", "significant_b", "effect_a", "qvalue_a", "significant_a"],
          *["n_tested", "n_significant", "direction", "consensus", "conflict"],
      ]
      assert out.index.name == "feature" and out["direction"].dtype == np.int8


  def test_q_equal_to_alpha_is_not_called():
      out = bt.da.consensus([_table("a", [1.0], [0.05]), _table("b", [1.0], [0.0499])])
      assert out["significant_a"].tolist() == [False] and out["n_significant"].tolist() == [1]


  def test_min_methods_sets_how_many_must_agree():
      results = [_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01]), _table("c", [1.0], [0.4])]
      assert bt.da.consensus(results, min_methods=2)["consensus"].tolist() == [True]
      assert not bt.da.consensus(results, min_methods=3)["consensus"].any()


  def test_conflict_is_flagged_and_never_consensus():
      out = bt.da.consensus([_table("a", [2.0], [0.01]), _table("b", [-2.0], [0.01]), _table("c", [2.0], [0.01])])
      assert out["conflict"].tolist() == [True] and out["direction"].tolist() == [0]
      assert out["consensus"].tolist() == [False]


  def test_untested_is_not_counted_as_not_significant():
      out = bt.da.consensus([_table("a", [1.0, 1.0], [0.01, np.nan]), _table("b", [1.0, 1.0], [0.01, 0.01])])
      assert out["n_tested"].tolist() == [2, 1] and out["n_significant"].tolist() == [2, 1]
      assert not out["significant_a"].iloc[1] and np.isnan(out["qvalue_a"].iloc[1])


  def test_features_are_the_union_of_the_tables():
      a = _table("a", [1.0, 1.0], [0.01, 0.01], features=["f1", "f2"])
      b = _table("b", [1.0, 1.0], [0.01, 0.01], features=["f2", "f3"])
      out = bt.da.consensus([a, b])
      assert out.index.tolist() == ["f1", "f2", "f3"]
      assert out["n_tested"].tolist() == [1, 2, 1] and out["consensus"].tolist() == [False, True, False]


  def test_results_with_different_contrasts_raise():
      with pytest.raises(ValueError, match=r"different contrasts \['A vs B', 'B vs A'\]"):
          bt.da.consensus([_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01], contrast="A vs B")])


  def test_repeated_method_raises():
      with pytest.raises(ValueError, match=r"repeat the method\(s\) \['a'\]"):
          bt.da.consensus([_table("a", [1.0], [0.01]), _table("a", [1.0], [0.02])])


  @pytest.mark.parametrize("min_methods", [0, 3])
  def test_min_methods_out_of_range_raises(min_methods):
      with pytest.raises(ValueError, match=r"min_methods must be between 1 and the number of results \(2\)"):
          bt.da.consensus([_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01])], min_methods=min_methods)


  @pytest.mark.parametrize("alpha", [0, 1, 1.5])
  def test_alpha_out_of_range_raises(alpha):
      with pytest.raises(ValueError, match="alpha must be between 0 and 1"):
          bt.da.consensus([_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01])], alpha=alpha)


  def test_one_table_instead_of_a_list_raises():
      with pytest.raises(TypeError, match=r"a list of da result tables, such as \[table_a, table_b\]"):
          bt.da.consensus(_table("a", [1.0], [0.01]))


  def test_consensus_of_the_native_methods():
      tdata = bt.datasets.toy()
      out = bt.da.consensus([bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")])
      assert out.index[out["consensus"]].tolist() == ["f6", "f7"] and (out["direction"][out["consensus"]] == 1).all()


  def test_consensus_keeps_its_inputs():
      results = [_table("a", [1.0, -2.0], [0.01, 0.2]), _table("b", [1.0, -1.0], [0.01, np.nan])]
      before = [table.copy() for table in results]
      bt.da.consensus(results)
      for table, copy in zip(results, before, strict=True):
          pd.testing.assert_frame_equal(table, copy)


  @settings(deadline=None)
  @given(
      st.lists(
          st.lists(
              st.tuples(st.floats(-3, 3), st.sampled_from([0.001, 0.04, 0.05, 0.3, np.nan])), min_size=4, max_size=4
          ),
          min_size=1,
          max_size=4,
      ),
      st.integers(1, 4),
  )
  def test_consensus_is_enough_calls_with_one_shared_sign(tables, min_methods):
      """The Decision's definition, both ways: n_significant >= min_methods and every call has the same non-zero sign."""
      min_methods = min(min_methods, len(tables))
      results = [_table(f"m{i}", [effect for effect, _ in rows], [q for _, q in rows]) for i, rows in enumerate(tables)]
      out = bt.da.consensus(results, min_methods=min_methods)
      for row, feature in enumerate(out.index):
          calls = [np.sign(rows[row][0]) for rows in tables if rows[row][1] < 0.05]
          shared = len(set(calls)) == 1 and calls[0] != 0
          assert out["n_significant"].iloc[row] == len(calls)
          assert out["consensus"].iloc[row] == (len(calls) >= min_methods and shared), feature
          assert out["conflict"].iloc[row] == (1 in calls and -1 in calls)
          assert out["direction"].iloc[row] == (calls[0] if shared else 0)


  def test_a_call_with_an_effect_of_exactly_zero_has_no_direction():
      out = bt.da.consensus([_table("a", [0.0], [0.01]), _table("b", [1.0], [0.01])])
      assert out["n_significant"].tolist() == [2] and out["direction"].tolist() == [0]
      assert out["consensus"].tolist() == [False] and out["conflict"].tolist() == [False]


  def test_empty_results_raise():
      with pytest.raises(ValueError, match=r"results is empty; pass at least one method's table"):
          bt.da.consensus([])


  @pytest.mark.parametrize(
      ("option", "value"), [("min_methods", True), ("min_methods", 1.5), ("alpha", "0.05"), ("alpha", True)]
  )
  def test_option_of_the_wrong_type_raises(option, value):
      with pytest.raises(TypeError, match=option):
          bt.da.consensus([_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01])], **{option: value})


  def test_numpy_scalars_are_valid_options():
      results = [_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01])]
      out = bt.da.consensus(results, alpha=np.float32(0.05), min_methods=np.int64(1))
      assert out["consensus"].tolist() == [True]
  ```
  Append to `tests/da/test_schema.py` (the schema checks, through the public
  API):
  ```diff
  @@ -53,3 +53,44 @@ def test_direction_is_the_sign_and_qvalue_bounds_pvalue(method, counts, split):
       assert (out.loc[tested, "qvalue"] >= out.loc[tested, "pvalue"] - 1e-15).all()
       assert out.loc[tested, ["pvalue", "qvalue"]].stack().between(0, 1).all()
       assert out.loc[~tested, ["effect", "qvalue"]].isna().all().all()
  +
  +
  +def _broken(change):
  +    table = bt.da.linda(bt.datasets.toy(), "group")
  +    change(table)
  +    return table
  +
  +
  +def _drop_effect(table, feature):
  +    table.loc[feature, ["effect", "direction"]] = [np.nan, 0]
  +
  +
  +@pytest.mark.parametrize(
  +    ("change", "message"),
  +    [
  +        (lambda t: t.drop(columns="se", inplace=True), r"lacks the result columns \['se'\]"),
  +        (lambda t: t.__setitem__("direction", t["direction"].astype(float)), "direction an integer"),
  +        (lambda t: t.__setitem__("pvalue", t["pvalue"] * 10), "must lie between 0 and 1"),
  +        (lambda t: t.__setitem__("qvalue", t["qvalue"].where(t.index != "f1")), "NaN exactly where pvalue is"),
  +        (lambda t: t.__setitem__("direction", -t["direction"]), "direction must be the sign of effect"),
  +        (lambda t: t.__setitem__("contrast", ["B vs A"] * 7 + ["A vs B"]), "must hold one contrast"),
  +        (lambda t: t.__setitem__("contrast", [np.nan] + ["B vs A"] * 7), "must hold one contrast"),
  +        (lambda t: _drop_effect(t, "f1"), "effect must be NaN exactly where pvalue is"),
  +    ],
  +    ids=["column", "dtype", "range", "missing q", "direction", "contrast", "nan contrast", "missing effect"],
  +)
  +def test_consensus_refuses_tables_that_break_the_schema(change, message):
  +    table = _broken(change)
  +    with pytest.raises(ValueError, match=message):
  +        bt.da.consensus([table, bt.da.ancombc2(bt.datasets.toy(), "group")])
  +
  +
  +def test_consensus_refuses_repeated_features():
  +    table = bt.da.linda(bt.datasets.toy(), "group")
  +    with pytest.raises(ValueError, match=r"results\[0\] repeats features \['f1'\]"):
  +        bt.da.consensus([table.set_axis(["f1"] * 8), bt.da.ancombc2(bt.datasets.toy(), "group")])
  +
  +
  +def test_consensus_refuses_what_is_not_a_table():
  +    with pytest.raises(TypeError, match=r"results\[1\] must be a result table of a bt.da method, got str"):
  +        bt.da.consensus([bt.da.linda(bt.datasets.toy(), "group"), "ancombc2"])
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/da -q`
  -> `26 failed, 67 passed, 4 deselected` (`AttributeError: module
  'biotapy.da' has no attribute 'consensus'`).
- [x] **Step 3: Implement.** `src/biotapy/da/_consensus.py`:
  ```python
  """Where differential abundance methods agree: one row per feature over the result tables the user computed."""

  from collections.abc import Sequence

  import numpy as np
  import pandas as pd

  from ._schema import validate_result


  def consensus(results: Sequence[pd.DataFrame], *, alpha: float = 0.05, min_methods: int = 2) -> pd.DataFrame:
      """Count, per feature, the methods that call it significant and whether they agree on its direction.

      Parameters
      ----------
      results
          Result tables of different ``bt.da`` methods run on the same data, ``group``
          and ``reference``, such as ``[bt.da.ancombc2(t, "host"), bt.da.linda(t, "host")]``.
      alpha
          A method calls a feature significant when its ``qvalue`` is below ``alpha``.
      min_methods
          How many methods must call a feature, all in the same direction, for consensus.

      Returns
      -------
      pandas.DataFrame
          One row per feature in any table (in first-seen order), indexed by ``feature``:
          ``effect_<method>``, ``qvalue_<method>`` and ``significant_<method>`` for each
          table in ``results`` order, then ``n_tested`` (methods with a p-value),
          ``n_significant``, ``direction`` (the sign the calling methods share; 0 when
          none calls the feature or they disagree), ``consensus`` and ``conflict``
          (methods call it in opposite directions).

      Raises
      ------
      TypeError
          ``results`` is a single table, or holds something that is not a table;
          ``alpha`` is not a real number or ``min_methods`` is not an ``int``.
      ValueError
          A table is not a ``bt.da`` result; two tables come from the same method or
          compare different contrasts; ``alpha`` is not between 0 and 1;
          ``min_methods`` is not between 1 and the number of tables; ``results`` is empty.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/differential_abundance`

      Calls are strict: ``qvalue == alpha`` is not significant. A feature a method
      did not test (NaN, or absent from its table) counts as not tested, never as
      not significant. A consensus feature has ``n_significant >= min_methods`` and
      every calling method has the same non-zero direction; a conflict is never a
      consensus. A call whose ``effect`` is exactly 0 has no direction: it counts in
      ``n_significant``, but ``direction`` is then 0, with no conflict and no
      consensus. The rule is recorded in the decision ``da-consensus-agreement``.

      Agreement between methods is a robustness report, not a way to choose a
      method: decide which methods to run before looking at their results.
      Methods that share a model (ANCOM-BC and ANCOM-BC2, say) agree more often for
      that reason alone.

      Examples
      --------
      >>> import biotapy as bt
      >>> tdata = bt.datasets.toy()
      >>> table = bt.da.consensus([bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")])
      >>> table.index[table["consensus"]].tolist()
      ['f6', 'f7']
      """
      if isinstance(results, pd.DataFrame):
          msg = "results must be a list of da result tables, such as [table_a, table_b]; got one table"
          raise TypeError(msg)
      tables = [validate_result(table, arg=f"results[{i}]") for i, table in enumerate(results)]
      methods = _check(tables, alpha=alpha, min_methods=min_methods)
      features = tables[0].index.append([table.index for table in tables[1:]]).unique()
      parts = [table.reindex(features) for table in tables]
      called = np.column_stack([(part["qvalue"] < alpha).to_numpy() for part in parts])
      signs = np.column_stack([np.sign(part["effect"].fillna(0)).to_numpy() for part in parts])
      up, down = (called & (signs > 0)).sum(axis=1), (called & (signs < 0)).sum(axis=1)
      n_significant = called.sum(axis=1)
      direction = np.where((up == n_significant) & (up > 0), 1, np.where((down == n_significant) & (down > 0), -1, 0))
      columns: dict[str, object] = {}
      for method, part, significant in zip(methods, parts, called.T, strict=True):
          columns |= {f"effect_{method}": part["effect"], f"qvalue_{method}": part["qvalue"]}
          columns[f"significant_{method}"] = significant
      out = pd.DataFrame(columns, index=features.rename("feature"))
      out["n_tested"] = np.column_stack([part["pvalue"].notna().to_numpy() for part in parts]).sum(axis=1)
      out["n_significant"] = n_significant
      out["direction"] = direction.astype(np.int8)
      out["consensus"] = (n_significant >= min_methods) & (direction != 0)
      out["conflict"] = (up > 0) & (down > 0)
      return out


  def _check_options(n_tables: int, *, alpha: float, min_methods: int) -> None:
      """Raise unless ``alpha`` is a real number in (0, 1) and ``min_methods`` an int in [1, n_tables]."""
      if n_tables == 0:
          msg = "results is empty; pass at least one method's table"
          raise ValueError(msg)
      if isinstance(alpha, bool) or not isinstance(alpha, int | float | np.integer | np.floating):
          msg = f"alpha must be a real number, got {type(alpha).__name__}"
          raise TypeError(msg)
      if isinstance(min_methods, bool) or not isinstance(min_methods, int | np.integer):
          msg = f"min_methods must be an int, got {type(min_methods).__name__}"
          raise TypeError(msg)
      if not 0 < alpha < 1:
          msg = f"alpha must be between 0 and 1, got {alpha}"
          raise ValueError(msg)
      if not 1 <= min_methods <= n_tables:
          msg = f"min_methods must be between 1 and the number of results ({n_tables}), got {min_methods}"
          raise ValueError(msg)


  def _check(tables: list[pd.DataFrame], *, alpha: float, min_methods: int) -> list[str]:
      """The tables' method names, after checking they can be compared under ``alpha`` and ``min_methods``."""
      _check_options(len(tables), alpha=alpha, min_methods=min_methods)
      methods = [str(table["method"].iloc[0]) for table in tables]
      repeated = sorted({method for method in methods if methods.count(method) > 1})
      if repeated:
          msg = f"results repeat the method(s) {repeated}; pass one table per method"
          raise ValueError(msg)
      contrasts = sorted({str(table["contrast"].iloc[0]) for table in tables})
      if len(contrasts) > 1:
          msg = f"results compare different contrasts {contrasts}; run every method with the same group and reference"
          raise ValueError(msg)
      return methods
  ```
  `src/biotapy/da/_schema.py` gains `COLUMNS` and the validator:
  ```diff
  @@ -5,6 +5,8 @@ import numpy.typing as npt
   import pandas as pd
   from scipy.stats import false_discovery_control

  +COLUMNS = ("effect", "se", "pvalue", "qvalue", "direction", "method", "contrast")
  +

   def result(
       features: "pd.Index[str]",
  @@ -31,3 +33,46 @@ def result(
           },
           index=pd.Index(features, name="feature"),
       )
  +
  +
  +def validate_result(table: object, *, arg: str) -> pd.DataFrame:
  +    """``table`` if it is a result table of one method and contrast; ``arg`` names it in errors."""
  +    if not isinstance(table, pd.DataFrame):
  +        msg = f"{arg} must be a result table of a bt.da method, got {type(table).__name__}"
  +        raise TypeError(msg)
  +    missing = [column for column in COLUMNS if column not in table.columns]
  +    if missing:
  +        msg = f"{arg} lacks the result columns {missing}; pass tables that bt.da methods return"
  +        raise ValueError(msg)
  +    floats = table[["effect", "se", "pvalue", "qvalue"]].dtypes
  +    if not all(pd.api.types.is_float_dtype(dtype) for dtype in floats) or not pd.api.types.is_integer_dtype(
  +        table["direction"]
  +    ):
  +        msg = f"{arg}: effect, se, pvalue and qvalue must be floats and direction an integer"
  +        raise ValueError(msg)
  +    if not table.index.is_unique:
  +        msg = f"{arg} repeats features {table.index[table.index.duplicated()].unique()[:3].tolist()}"
  +        raise ValueError(msg)
  +    _check_values(table, arg=arg)
  +    return table
  +
  +
  +def _check_values(table: pd.DataFrame, *, arg: str) -> None:
  +    """Raise unless p-values, q-values and effects are missing together, p and q are probabilities, direction is sign(effect), one method and contrast."""
  +    probabilities = table[["pvalue", "qvalue"]]
  +    if not (probabilities.isna() | probabilities.ge(0) & probabilities.le(1)).all().all():
  +        msg = f"{arg}: pvalue and qvalue must lie between 0 and 1 (NaN for an untested feature)"
  +        raise ValueError(msg)
  +    if not table["pvalue"].isna().equals(table["qvalue"].isna()):
  +        msg = f"{arg}: qvalue must be NaN exactly where pvalue is"
  +        raise ValueError(msg)
  +    if not table["effect"].isna().equals(table["pvalue"].isna()):
  +        msg = f"{arg}: effect must be NaN exactly where pvalue is (an untested feature has no estimate)"
  +        raise ValueError(msg)
  +    if not (table["direction"].to_numpy() == np.sign(np.nan_to_num(table["effect"].to_numpy()))).all():
  +        msg = f"{arg}: direction must be the sign of effect, 0 where effect is NaN"
  +        raise ValueError(msg)
  +    for column in ("method", "contrast"):
  +        if table[column].nunique(dropna=False) != 1:
  +            msg = f"{arg} must hold one {column}, found {table[column].unique().tolist()}"
  +            raise ValueError(msg)
  ```
  `src/biotapy/da/__init__.py`:
  ```diff
  @@ -1,6 +1,7 @@
   """Differential abundance: methods that share one result table (contracts/data-model-slots, DA results)."""

   from ._ancombc import ancombc2
  +from ._consensus import consensus
   from ._linda import linda

  -__all__ = ["ancombc2", "linda"]
  +__all__ = ["ancombc2", "consensus", "linda"]
  ```
- [x] **Step 4: Run, expect pass** - the same command -> `102 passed, 4 deselected`. The two
  property tests also under three extra seeds.
- [x] **Step 5: Decision concept.** Create
  `.knowledge/decisions/da-consensus-agreement.md` (`generated.at` is the
  commit's UTC time and `commit` the previous commit):
  ```markdown
  ---
  type: Decision
  title: What "methods agree" means in da.consensus
  description: A feature is a consensus hit when at least min_methods methods call it at q < alpha and every calling method gives it the same sign; untested is not "not significant", opposite calls are a conflict, and the user picks the methods.
  tags: [da, statistics, api]
  status: stable
  paths: ["src/biotapy/da/_consensus.py", "src/biotapy/da/_schema.py"]
  generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T20:05:08Z }
  commit: 95a35a1
  sources:
    - id: nearing
      resource: https://www.nature.com/articles/s41467-022-28034-z
      title: Nearing et al. 2022, Microbiome differential abundance methods produce different results across 38 datasets, Nature Communications
    - id: pelto
      resource: https://academic.oup.com/bib/article/26/2/bbaf130/8093585
      title: Pelto et al. 2025, Elementary methods provide more replicable results in microbial differential abundance analysis, Briefings in Bioinformatics
    - id: oma
      resource: https://microbiome.github.io/OMA/docs/devel/pages/differential_abundance.html
      title: Orchestrating Microbiome Analysis, Differential abundance chapter
  ---

  # Context
  Differential abundance methods disagree: on 38 datasets, 14 methods called very
  different sets of features.[^nearing] Phase 3 runs several methods behind one
  result schema (`contracts/data-model-slots`, DA results) and needs one rule for
  "where they agree" that every reader of `da.consensus`'s table can check by hand.
  The literature is split on the remedy: Nearing et al. recommend a consensus of
  several methods; Pelto et al. and the OMA book recommend one elementary method
  and warn that trying methods until one agrees is selective reporting.[^pelto][^oma]

  # Decision
  Per feature `f` and method `m`, over result tables the user computed:

  - `called(f, m)` is `qvalue < alpha`, strictly: `q == alpha` is not called.
    Calls are recomputed from `qvalue`, never taken from a method's own flag
    (LinDA's `reject` uses `<=`; ANCOM-BC2's `Signif` sits on Holm by default).
  - Every `qvalue` is Benjamini-Hochberg over the features that method tested,
    so "called" means the same false discovery rate in every column.
  - `n_tested(f)` counts methods with a finite `pvalue`; a feature a method did
    not test (NaN, or absent from its table) is "not tested", never "not
    significant". `n_significant(f)` counts the calls.
  - `consensus(f)` is true when `n_significant >= min_methods` and every calling
    method has the same non-zero `direction`. `conflict(f)` is true when calling
    methods disagree in sign; a conflict is never a consensus. `direction(f)` is
    the shared sign, 0 when there is no call, a conflict, or a call whose
    `effect` is exactly 0 (a call with no sign: it counts in `n_significant`,
    blocks consensus on that feature, and is not a conflict).
  - `da.consensus(results, *, alpha=0.05, min_methods=2)` only combines tables:
    it never runs a method, so the list of methods is written in the user's code
    before any result is seen. Tables must come from different methods and
    compare the same `contrast`.

  # Rejected
  - **Intersection of all methods** (`min_methods = len(results)` as the only
    rule): one method that cannot test a feature (ANCOM-BC2 on a feature absent
    from one group) would veto it. Still available as `min_methods=len(results)`.
  - **k of n without direction**: two methods calling a feature in opposite
    directions would count as agreement.
  - **Agreement among all methods with an estimate** (every method's sign, called
    or not): a method that tested a feature but found nothing would block a
    consensus on a sign it never claimed.
  - **Rank aggregation** (combining p-values or ranks): needs p-values that are
    comparable across methods, which they are not, and options no use case asks
    for (rules.md R2.3).
  - **A one-call `da.consensus(adata, formula, *, methods, ...)`** that runs the
    methods: eight arguments (rules.md R5 allows six), a policy for a method that
    fails (R7.4 forbids dropping it silently), per-method options, and it hides
    which methods were chosen.

  # Consequences
  - The table is a robustness report, not a way to choose a method; the guide and
    `da.consensus`'s Notes say to fix the methods before looking.
  - Methods that share a model (ANCOM-BC and ANCOM-BC2) agree more often for that
    reason alone; consensus counts methods, not independent evidence.
  - Each method's call is stored as `significant_<method>`, so a reader of the
    table, or a plot of it, needs no `alpha`.

  [^nearing]: Nearing et al. 2022, Microbiome differential abundance methods produce different results across 38 datasets, Nature Communications
  [^pelto]: Pelto et al. 2025, Elementary methods provide more replicable results in microbial differential abundance analysis, Briefings in Bioinformatics
  [^oma]: Orchestrating Microbiome Analysis, Differential abundance chapter
  ```
  and list it in `.knowledge/decisions/index.md`:
  ```diff
  @@ -7,6 +7,7 @@
   * [Python first, compiled code last](python-first-compiled-last.md) - Pure Python/NumPy by default; delegate to compiled libraries; Numba then Rust only for a benchmarked hotspot; never new C/C++.
   * [Heavy dependencies are optional extras](optional-heavy-dependencies.md) - torch, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily.
   * [Knowledge in OKF, user docs in Sphinx](docs-okf-and-sphinx.md) - Contributor and agent knowledge is an OKF v0.2 bundle in .knowledge/; user docs are a Sphinx site in docs/.
  +* [What "methods agree" means in da.consensus](da-consensus-agreement.md) - A feature is a consensus hit when at least min_methods methods call it at q < alpha and every calling method gives it the same sign; untested is not "not significant", opposite calls are a conflict, and the user picks the methods.
   * [R bridge before native ports](r-bridge-before-ports.md) - R-only DA methods ship first through an optional rpy2 bridge; native ports only for the most used.
   * [No bundled KEGG mapping files](no-bundled-kegg.md) - Function hierarchies come from the user's local files or from ENZYME (CC BY 4.0, `bt.datasets.enzyme`); KEGG and MetaCyc are never shipped or fetched.
   * [Package name biotapy](package-name-biotapy.md) - Distribution and import name is biotapy, hosted at github.com/pedrocr83/biotapy.
  ```
- [x] **Step 6: Contract and rule (approved Phase 3 decision 9).**
  `.knowledge/contracts/function-shape.md` statement 1:
  ```diff
  @@ -21,6 +21,8 @@ Every public function in `io`, `datasets`, `pp`, `tl`, `fn`, `da`, `ml`, `pl`:
   1. **Signature**: `verb(data, <required args>, *, <options>) -> <result>`.
      - `data` is annotated with the widest type that works: `AnnData` when no tree
        is needed, `TreeData` when `vart` is read, `MuData` for multi-modal.
  +     `da.consensus` takes a sequence of `da` result tables instead: it combines
  +     results, not data.
      - Everything after the required arguments is keyword-only (`*`).
      - No `**kwargs` pass-through, except a documented `plot_kwargs` in `pl`.
   2. **Return and mutation**: per [pure-by-default](/decisions/pure-by-default.md).
  ```
  `rules.md` R3.2:
  ```diff
  @@ -57,6 +57,7 @@ Full contract: [.knowledge/contracts/function-shape.md](.knowledge/contracts/fun
     `verb(data, <required>, *, <options>) -> result`. Options are keyword-only.
   - **R3.2** Annotate `data` with the widest type that works: `AnnData` unless the
     tree is needed (`TreeData`) or several modalities are (`MuData`).
  +  `da.consensus` takes `da` result tables instead.
   - **R3.3** Purity ([pure-by-default](.knowledge/decisions/pure-by-default.md)):
     `io`/`pp`/`fn`/`da` return new objects and never mutate input; `tl` returns
     its result, and writes to the documented slot only with `inplace=True`;
  ```
  Log lines:
  ```text
  - **Create**: [da-consensus-agreement](decisions/da-consensus-agreement.md): what "methods agree" means in `da.consensus` (strict `q < alpha`, one BH, untested is not "not significant", same sign, conflict), the options rejected, and that the user picks the methods; listed in [decisions/index.md](decisions/index.md).
  - **Update**: [function-shape](contracts/function-shape.md) (and rules.md R3.2): `da.consensus` takes `da` result tables instead of an AnnData. [phase-3-stats](roadmap/phase-3-stats.md) task 3.8 done.
  ```
- [x] **Step 7: Docs.** Append to `docs/guide/differential_abundance.md` and
  list the function:
  ````diff
  @@ -78,3 +78,27 @@ pseudo_sens = FALSE)`: biotapy does not run R's pseudocount sensitivity analysis
   use, biotapy's effects are within 0.012 log2 of R's and the significant genera are the same; the
   small difference comes from the bias estimate, whose iterations stop at R's cap of 100 before they
   have converged on that data, in R as in scikit-bio.
  +
  +## Where methods agree
  +
  +`bt.da.consensus` puts the tables of several methods side by side and counts, for each feature,
  +the methods that call it significant (`qvalue < alpha`, strictly) and whether they agree on its
  +direction:
  +
  +```python
  +results = [bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")]
  +table = bt.da.consensus(results)  # alpha=0.05, min_methods=2
  +table[table["consensus"]]
  +```
  +
  +A feature is a consensus hit when at least `min_methods` methods call it and all of them give it
  +the same sign. Methods that call it in opposite directions mark a `conflict`, which is never a
  +consensus. A method that could not test a feature does not count against it: `n_tested` says how
  +many methods tested each feature. The tables must come from different methods and compare the same
  +`contrast`, so run every method with the same `group` and `reference`.
  +
  +Methods disagree a lot on real data, and the literature is split on what to do about it: Nearing
  +et al. (2022) recommend a consensus of several methods, Pelto et al. (2025) one simple method. Either
  +way, choose the methods before you look at their results; trying methods until one finds what you
  +hoped for is selective reporting, and a consensus table does not protect against it. Methods that
  +share a model, such as ANCOM-BC and ANCOM-BC2, also agree more often for that reason alone.
  ````
  ````diff
  @@ -101,6 +101,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
       :toctree: generated

       da.ancombc2
  +    da.consensus
       da.linda
   ```

  ````
- [x] **Step 8: Gate and commit**
  ```bash
  git add src/biotapy/da/_consensus.py src/biotapy/da/_schema.py src/biotapy/da/__init__.py tests/da/test_consensus.py \
    tests/da/test_schema.py docs/guide/differential_abundance.md docs/api.md rules.md \
    .knowledge/contracts/function-shape.md .knowledge/decisions/da-consensus-agreement.md \
    .knowledge/decisions/index.md .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  uvx prek run --all-files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git commit -m "feat(da): add the consensus table and the agreement decision"
  ```

### Task 3.9: `pl.consensus`

**Files:** modify `src/biotapy/pl/__init__.py`,
`docs/guide/differential_abundance.md`, `docs/api.md`,
`.knowledge/contracts/function-shape.md`, `rules.md`; create
`src/biotapy/pl/_consensus.py`, `tests/pl/test_consensus.py`.
**Not touched:** `pl/_common.py` (reused as is), `da`.
**Interfaces:**
- Consumes: the table `bt.da.consensus` returns (its column names only);
  `pl/_common.py`'s `new_axes`, `label_ticks`, `MISSING_COLOR`.
- Produces: `bt.pl.consensus(table: pd.DataFrame, *, top: int = 30, ax: Axes |
  None = None) -> Axes`.

- [x] **Step 1: Failing tests.** `tests/pl/test_consensus.py` (the `ax`
  fixture is `tests/pl/conftest.py`'s untracked figure):
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  from matplotlib.colors import to_rgba

  import biotapy as bt


  def _consensus_table():
      """da.consensus's layout for three features and two methods, hand-built."""
      return pd.DataFrame(
          {
              "effect_a": [2.0, -1.0, 0.1, 0.5],
              "qvalue_a": [0.01, 0.02, 0.5, np.nan],
              "significant_a": [True, True, False, False],
              "effect_b": [1.0, -3.0, 4.0, 0.2],
              "qvalue_b": [0.01, 0.30, 0.01, 0.9],
              "significant_b": [True, False, True, False],
              "n_tested": [2, 2, 2, 1],
              "n_significant": [2, 1, 1, 0],
              "direction": np.array([1, -1, 1, 0], dtype=np.int8),
              "consensus": [True, False, False, False],
              "conflict": [False] * 4,
          },
          index=pd.Index(["x", "y", "z", "w"], name="feature"),
      )


  def _dots(ax):
      return {collection.get_label(): len(collection.get_offsets()) for collection in ax.collections}


  def test_consensus_draws_one_dot_per_call_and_per_tested_feature(ax):
      assert bt.pl.consensus(_consensus_table(), ax=ax) is ax
      # x is called by both, y by a only (b tested it), z by b only (a tested it); w is called by none and not drawn.
      assert _dots(ax) == {"effect > 0": 3, "effect < 0": 1, "not significant": 2}
      assert [label.get_text() for label in ax.get_xticklabels()] == ["a", "b"]


  def test_rows_are_sorted_by_calls_then_mean_absolute_effect(ax):
      bt.pl.consensus(_consensus_table(), ax=ax)
      # x has two calls; z (mean |effect| 2.05) comes before y (2.0).
      assert [label.get_text() for label in ax.get_yticklabels()] == ["x", "z", "y"]
      assert [label.get_fontweight() for label in ax.get_yticklabels()] == ["bold", "normal", "normal"]
      assert ax.get_ylim() == (2.5, -0.5)


  def test_top_keeps_the_first_rows(ax):
      bt.pl.consensus(_consensus_table(), top=1, ax=ax)
      assert [label.get_text() for label in ax.get_yticklabels()] == ["x"]
      assert _dots(ax) == {"effect > 0": 2}


  def test_nothing_called_raises(ax):
      table = _consensus_table()
      table["n_significant"] = 0
      with pytest.raises(ValueError, match="no method calls any feature significant: nothing to draw"):
          bt.pl.consensus(table, ax=ax)


  def test_top_below_one_raises(ax):
      with pytest.raises(ValueError, match="top must be at least 1, got 0"):
          bt.pl.consensus(_consensus_table(), top=0, ax=ax)


  def test_a_result_table_instead_of_a_consensus_table_raises(ax):
      with pytest.raises(KeyError, match="pl.consensus draws the table bt.da.consensus returns"):
          bt.pl.consensus(bt.da.linda(bt.datasets.toy(), "group"), ax=ax)


  def test_consensus_keeps_the_table(ax):
      table = _consensus_table()
      before = table.copy()
      bt.pl.consensus(table, ax=ax)
      pd.testing.assert_frame_equal(table, before)


  def test_dots_sit_at_their_row_and_method_with_the_colour_of_their_sign(ax):
      bt.pl.consensus(_consensus_table(), ax=ax)
      by_label = {collection.get_label(): collection for collection in ax.collections}
      # Rows are x, z, y (top first); columns are a, b. x is called up by both, y down by a only.
      up, down = by_label["effect > 0"], by_label["effect < 0"]
      assert sorted(map(tuple, up.get_offsets().tolist())) == [(0, 0), (1, 0), (1, 1)]
      assert down.get_offsets().tolist() == [[0, 2]]
      assert up.get_facecolor().tolist() == [list(to_rgba("#d62728"))]
      assert down.get_facecolor().tolist() == [list(to_rgba("#1f77b4"))]


  def test_a_non_table_raises(ax):
      with pytest.raises(TypeError, match="table must be a pandas DataFrame, got list"):
          bt.pl.consensus([1, 2], ax=ax)


  @pytest.mark.parametrize("top", [2.5, True, "3"])
  def test_top_that_is_not_an_integer_raises(ax, top):
      with pytest.raises(TypeError, match="top must be an integer"):
          bt.pl.consensus(_consensus_table(), top=top, ax=ax)


  def test_numpy_integer_top_is_accepted(ax):
      bt.pl.consensus(_consensus_table(), top=np.int64(1), ax=ax)
      assert [label.get_text() for label in ax.get_yticklabels()] == ["x"]


  def test_an_empty_table_raises(ax):
      with pytest.raises(ValueError, match="no method calls any feature significant: nothing to draw"):
          bt.pl.consensus(_consensus_table().iloc[0:0], ax=ax)


  def test_the_legend_is_above_the_axes_and_inside_the_figure(ax):
      bt.pl.consensus(_consensus_table(), ax=ax)
      ax.figure.canvas.draw()
      legend = ax.get_legend().get_window_extent()
      assert legend.y0 >= ax.get_window_extent().y1
      assert ax.figure.bbox.contains(legend.x0, legend.y0)
      assert ax.figure.bbox.contains(legend.x1, legend.y1)
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest
  tests/pl/test_consensus.py -q` -> `7 failed` (`AttributeError: module
  'biotapy.pl' has no attribute 'consensus'`).
- [x] **Step 3: Implement.** `src/biotapy/pl/_consensus.py`:
  ```python
  """Where differential abundance methods agree: a dot per feature and method, from da.consensus's table."""

  from typing import TYPE_CHECKING

  import numpy as np
  import pandas as pd

  from ._common import MISSING_COLOR, label_ticks, new_axes

  if TYPE_CHECKING:
      from matplotlib.axes import Axes

  # tab10's red and blue: a call's sign. A hollow grey dot is a feature the method tested but did not call.
  _KINDS = {
      "effect > 0": {"color": "#d62728"},
      "effect < 0": {"color": "#1f77b4"},
      "not significant": {"facecolors": "none", "edgecolors": MISSING_COLOR},
  }


  def _check(table: pd.DataFrame, *, top: int) -> list[str]:
      """The method names in ``table``, after checking it is a consensus table and ``top`` a positive int."""
      if not isinstance(table, pd.DataFrame):
          msg = f"table must be a pandas DataFrame, got {type(table).__name__}"
          raise TypeError(msg)
      methods = [column.removeprefix("significant_") for column in table.columns if column.startswith("significant_")]
      needed = ["n_significant", "consensus", *[f"{kind}_{m}" for m in methods for kind in ("effect", "qvalue")]]
      absent = [column for column in needed if column not in table.columns]
      if not methods or absent:
          msg = f"table lacks {absent or ['significant_<method>']}; pl.consensus draws the table bt.da.consensus returns"
          raise KeyError(msg)
      if isinstance(top, bool) or not isinstance(top, int | np.integer):
          msg = f"top must be an integer, got {top!r}"
          raise TypeError(msg)
      if top < 1:
          msg = f"top must be at least 1, got {top}"
          raise ValueError(msg)
      return methods


  def consensus(table: pd.DataFrame, *, top: int = 30, ax: "Axes | None" = None) -> "Axes":
      """A dot matrix of which methods call which features, from :func:`biotapy.da.consensus`.

      Parameters
      ----------
      table
          The table :func:`biotapy.da.consensus` returns.
      top
          Draw at most this many features: those called by the most methods, then with
          the largest mean absolute effect, then in table order. A ``top`` above 30
          needs a taller figure, passed as ``ax``.
      ax
          Axes to draw on; by default a new figure's.

      Returns
      -------
      matplotlib.axes.Axes
          One row per feature (top row first) and one column per method: a filled dot
          coloured by the sign of the effect where the method calls the feature, a
          hollow dot where it tested the feature without calling it, nothing where it did
          not test it. Consensus features have bold labels.

      Raises
      ------
      TypeError
          ``table`` is not a DataFrame, or ``top`` is not an integer.
      KeyError
          ``table`` lacks the columns :func:`biotapy.da.consensus` writes.
      ValueError
          ``top`` is below 1, or no method calls any feature.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/differential_abundance`

      One axes, as every ``pl`` function returns; an UpSet plot of the same calls
      needs two panels. For a two-level ``group`` the effects of all methods are log2
      fold changes, which makes their mean comparable when ranking features. For a
      numeric ``group`` they are not (``da.ancombc2`` gives the change per unit,
      ``da.linda`` per standard deviation), so the ranking mixes units: rank by
      ``n_significant`` and read each method's effect on its own.

      Examples
      --------
      >>> import biotapy as bt
      >>> tdata = bt.datasets.toy()
      >>> table = bt.da.consensus([bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")])
      >>> ax = bt.pl.consensus(table)
      >>> [label.get_text() for label in ax.get_yticklabels()]
      ['f6', 'f7', 'f8']
      """
      methods = _check(table, top=top)
      called = table[table["n_significant"] > 0]
      if called.empty:
          msg = "no method calls any feature significant: nothing to draw"
          raise ValueError(msg)
      strength = called[[f"effect_{m}" for m in methods]].abs().mean(axis=1).to_numpy()
      rows = called.iloc[np.lexsort((-strength, -called["n_significant"].to_numpy()))[:top]]
      significant = rows[[f"significant_{m}" for m in methods]].to_numpy(bool)
      effect = rows[[f"effect_{m}" for m in methods]].to_numpy(np.float64)
      tested = rows[[f"qvalue_{m}" for m in methods]].notna().to_numpy()
      y, x = np.indices(significant.shape)
      masks = [significant & (effect > 0), significant & (effect < 0), tested & ~significant]
      ax = new_axes(ax)
      for (label, style), mask in zip(_KINDS.items(), masks, strict=True):
          if mask.any():
              ax.scatter(x[mask], y[mask], label=label, **style)  # type: ignore[arg-type]
      ax.set_xticks(np.arange(len(methods)), labels=methods)
      label_ticks(ax, rows.index.astype(str).tolist(), axis="y")
      for tick, bold in zip(ax.get_yticklabels(), rows["consensus"].to_numpy(bool), strict=False):
          tick.set_fontweight("bold" if bold else "normal")
      ax.set_xlim(-0.5, len(methods) - 0.5)
      ax.set_ylim(len(rows) - 0.5, -0.5)
      # One row above the axes: inside they cover dots, beside them a default save clips them.
      kinds = len(ax.get_legend_handles_labels()[0])
      ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=kinds, frameon=False, borderaxespad=0)
      return ax
  ```
  `src/biotapy/pl/__init__.py`:
  ```diff
  @@ -1,7 +1,8 @@
   """Plots of what tl, pp and fn give; pl computes nothing itself (contracts/module-boundaries)."""

   from ._abundance import bar, contributions, heatmap
  +from ._consensus import consensus
   from ._ordination import ordination, scree
   from ._richness import richness

  -__all__ = ["bar", "contributions", "heatmap", "ordination", "richness", "scree"]
  +__all__ = ["bar", "consensus", "contributions", "heatmap", "ordination", "richness", "scree"]
  ```
- [x] **Step 4: Run, expect pass** - the same command -> `15 passed`; the
  doctest (`uv run --group test pytest src/biotapy/pl/_consensus.py -q`) ->
  `1 passed`.
- [x] **Step 5: Docs, contract and rule.**
  ````diff
  @@ -102,3 +102,18 @@ et al. (2022) recommend a consensus of several methods, Pelto et al. (2025) one
   way, choose the methods before you look at their results; trying methods until one finds what you
   hoped for is selective reporting, and a consensus table does not protect against it. Methods that
   share a model, such as ANCOM-BC and ANCOM-BC2, also agree more often for that reason alone.
  +
  +## Plotting the consensus
  +
  +`bt.pl.consensus` draws the consensus table as a dot matrix: one row per feature that at least one
  +method calls, one column per method. A filled dot is a call, red for a positive effect and blue for
  +a negative one; a hollow dot is a feature the method tested without calling it; no dot, a feature
  +the method could not test. Consensus features have bold labels. Rows are sorted by the number of
  +methods calling them, then by their mean absolute effect, and `top` keeps the first 30:
  +
  +```python
  +ax = bt.pl.consensus(table, top=20)
  +```
  +
  +Feature ids are often accession numbers; to label rows with a rank, rename the table's index first,
  +for example `table.rename(index=tdata.var["genus"])`.
  ````
  ```diff
  @@ -115,6 +115,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
       :toctree: generated

       pl.bar
  +    pl.consensus
       pl.contributions
       pl.heatmap
       pl.ordination
  ```
  ```diff
  @@ -21,8 +21,8 @@ Every public function in `io`, `datasets`, `pp`, `tl`, `fn`, `da`, `ml`, `pl`:
   1. **Signature**: `verb(data, <required args>, *, <options>) -> <result>`.
      - `data` is annotated with the widest type that works: `AnnData` when no tree
        is needed, `TreeData` when `vart` is read, `MuData` for multi-modal.
  -     `da.consensus` takes a sequence of `da` result tables instead: it combines
  -     results, not data.
  +     `da.consensus` takes a sequence of `da` result tables and `pl.consensus` the
  +     table `da.consensus` returns instead: they combine and draw results, not data.
      - Everything after the required arguments is keyword-only (`*`).
      - No `**kwargs` pass-through, except a documented `plot_kwargs` in `pl`.
   2. **Return and mutation**: per [pure-by-default](/decisions/pure-by-default.md).
  ```
  ```diff
  @@ -57,7 +57,7 @@ Full contract: [.knowledge/contracts/function-shape.md](.knowledge/contracts/fun
     `verb(data, <required>, *, <options>) -> result`. Options are keyword-only.
   - **R3.2** Annotate `data` with the widest type that works: `AnnData` unless the
     tree is needed (`TreeData`) or several modalities are (`MuData`).
  -  `da.consensus` takes `da` result tables instead.
  +  `da.consensus` takes `da` result tables and `pl.consensus` its table instead.
   - **R3.3** Purity ([pure-by-default](.knowledge/decisions/pure-by-default.md)):
     `io`/`pp`/`fn`/`da` return new objects and never mutate input; `tl` returns
     its result, and writes to the documented slot only with `inplace=True`;
  ```
  Log line:
  ```text
  - **Update**: [function-shape](contracts/function-shape.md) (and rules.md R3.2): `pl.consensus` takes the table `da.consensus` returns. [phase-3-stats](roadmap/phase-3-stats.md) task 3.9 done.
  ```
- [x] **Step 6: Gate and commit**
  ```bash
  git add src/biotapy/pl/_consensus.py src/biotapy/pl/__init__.py tests/pl/test_consensus.py \
    docs/guide/differential_abundance.md docs/api.md rules.md .knowledge/contracts/function-shape.md \
    .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  uvx prek run --all-files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git commit -m "feat(pl): add the consensus dot matrix"
  ```

### Checkpoint B - review slice 3B
- [x] Review the whole slice (superpowers:requesting-code-review, opus: the
  LinDA port and the ANCOM-BC2 mapping are numerics) against every contract,
  pure-by-default, the Phase 3 and slice 3B review focus, and the measured
  parity; then a fix pass, one commit per finding, each with a test. Record
  the counts (Critical / Important / Minor) and the fix range here.
  Recorded: 0 Critical / 2 Important / 5 Minor; fix pass `9f45c28..fd9e94a`
  plus `927e5ae`; re-review 7/7 addressed.
- [x] Run the exit-gate check for 3B: `uv run --group test pytest -m "golden
  or network" tests/da -q` (expected `4 passed`) and the full
  `uv run --group test pytest -q -W error::UserWarning`; record both counts.
  Recorded at `927e5ae`: `4 passed, 104 deselected` for the golden/network
  run; `1174 passed, 29 deselected` for the full run (1172 before `modules/da.md` added two bundle tests).
- [x] Knowledge (codebase-map templates; R12.2-R12.4):
  - **Create `.knowledge/modules/da.md`** (`type: Module`, moved forward from
    3.14, slice 3B decision 13): Responsibility (two native methods behind one schema,
    the consensus table; `da` writes no slot); Entry points
    `_linda.py:linda`, `_mode`, `_shorth`, `_ancombc.py:ancombc2`,
    `_consensus.py:consensus`, `_design.py:model`, `design_matrix`,
    `dense_counts`, `_schema.py:result`, `validate_result`; Invariants: log2
    `effect`, BH over finite p-values, `var_names` order, NaN = not tested, no
    filtering, no formula strings, missing values and collinearity raise,
    numeric columns scaled in LinDA only, strict `q < alpha`; Gotchas:
    MicrobiomeStat's unreachable imputation branch, modeest returns the
    previous mean-shift value and truncates tied shorth starts, ANCOM-BC2's
    unconverged E-M at 100 iterations, scikit-bio's p = 1 for unfitted
    features, patsy's numeric-after-categorical term order, ANCOM-BC2 not
    exactly antisymmetric in the reference; Verification `uv run --group test
    pytest tests/da -q` and `-m "golden or network" tests/da`. Index entry in
    `modules/index.md`.
  - **Update `.knowledge/modules/pl.md`**: `_consensus.py:consensus` in
    Responsibility and Entry points; invariant: it reads `da.consensus`'s
    columns and imports nothing from `da`.
  - **Update `.knowledge/roadmap/phase-3-stats.md`**: tick Checkpoint B's
    boxes; task 3.14's line says it updates `modules/da.md` for the bridges.
  - **Verification bump** for the concepts `bash scripts/knowledge_stale.sh
    --touched` lists whose statements still hold (on the replay:
    phase-0-foundation, phase-1-core, engine-parity, module-boundaries; pl.md
    is updated above), as Checkpoint A did.
  - Log lines for each.
- [ ] Push `phase-3b`, open the PR and merge (merge commit) when CI is green,
  including docs, the network job and the knowledge report (push, PR and
  merge on green approved for Phase 3 slice branches, 2026-10-05).
- [ ] Ask the user to review slice 3B before slice 3C is expanded (R1.2a).

### Slice 3B decisions for the user
Recommended answer first. Items 1-5 change what the approved plan said or add
an image pin; 6-15 are judgement calls inside the approved design.

1. **LinDA's zeros follow what MicrobiomeStat 1.4 computes: 0.5 is added
   whenever `X` holds a zero.** MicrobiomeStat's adaptive imputation branch
   is unreachable (`"Imputation" == "imputation"` is false), so its output is
   always the pseudocount path, and the golden matches to 8e-13. Design note
   7 said "adaptive zero handling (... zeros are imputed by library size ...)".
   Alternative: port the paper's adaptive rule (standalone LinDA 0.2.0 does
   it); then biotapy differs from `MicrobiomeStat::linda` by up to 1.48 log2
   whenever library size depends on the model (any depth covariate), and the
   golden would need the GitHub-only LinDA package.
2. **ANCOM-BC2 keeps R's iteration defaults** (bias E-M capped at 100,
   scikit-bio's `max_iter=100`), so biotapy reproduces R's default output; on
   the benchmark that output is unconverged (R moves every effect by -0.28
   log2, 208 calls -> 220, at 1,000 iterations), which the docstring states.
   Alternative: raise the cap (a better estimate that departs from R and
   needs a parameter or a constant).
3. **ANCOM-BC2's unfitted features are NaN and left out of BH** (the schema's
   rule), where R and scikit-bio report p = 1 and count them: 36 of 636
   genera on the benchmark; with `log_depth` biotapy calls 65 genera where
   R's own q-values call 63. Alternative: keep p = 1 in the family (R's
   q-values exactly; "NaN = not tested" would no longer hold for ANCOM-BC2).
4. **ANCOM-BC2's golden tolerances are measured, not elementwise-exact:**
   `effect` atol 0.015 log2, `se` rtol 2e-3, `pvalue` atol 0.02, Spearman >
   0.9999, identical calls (recorded in r-golden-parity). Alternative: signs
   and ranks only (weaker). Exactness is not reachable through scikit-bio's
   API (its E-M stops on a different iterate).
5. **Image pins:** CVXR 1.0-15 from the CRAN archive (built in the image) and
   apt `libgsl27`, both only so ANCOMBC 2.12.0 installs on Bioconductor
   3.22's snapshot; `VERSIONS.txt` also records modeest 2.4.0, whose
   `mlv` LinDA's numbers depend on. Alternative: move the image to
   Bioconductor 3.23 / R 4.6 (ANCOMBC 2.14 drops CVXR), a pin bump that
   regenerates every golden.
6. **The consensus table gains `significant_<method>`** beside design note 6's
   `effect_<method>` and `qvalue_<method>`, so `pl.consensus` needs no `alpha`.
7. **`pl.consensus`'s legend says "effect > 0" / "effect < 0"**, not the
   contrast: the consensus table does not carry it (R2.3). Design note 5 said
   the plot labels the direction with the contrast. Alternative: a `contrast`
   column in the consensus table, used as the plot title.
8. **Methods do not re-validate their output.** They build it with
   `_schema.result`; `da.consensus` validates its inputs at the public
   boundary (R3.5), and `validate_result` arrives with it in 3.8 (R2.3).
   Design note 5 said every method calls it. The validator also requires
   `qvalue` NaN exactly where `pvalue` is, and one `method` and `contrast` per
   table.
9. **`pp.pseudocounted` stays in `pp`** (the ledger's "pseudocounted -> _core
   if LinDA reuses it" resolves to "not reused"): LinDA's rule is
   conditional on a zero, a constant, with no argument to validate and no
   scale warning, so the two share no logic (R4.3).
10. **Shared checks beyond design note 4:** `model()` also raises on a rank
    deficient design (collinear covariates, which `lstsq` would answer with an
    arbitrary minimum-norm fit), on no more samples than model terms, on an
    empty sample (`dense_counts`) and on `covariates="age"` given as a
    string (`TypeError`).
11. **No `Q()` quoting:** ANCOM-BC2 renames the model's columns to `x0..xk`
    before patsy (design note 4 said `Q()` and backticks) and finds the
    group's row by name, because patsy orders numeric terms after
    categorical ones.
12. **Benchmark contrast numbers:** 9 human samples (Feces 4, Skin 3, Tongue
    2) against 17 others, which include the 3 Mock communities; design note
    10 says 8 against 18 and calls the others "environmental". The goldens
    name the levels `human` and `other` (reference `other`). For the 3D
    notebook: keep "human vs other" (recommended), or drop Mock (9 vs 14).
13. **`modules/da.md` is created at Checkpoint B**, not in 3.14 (3D), so the
    module is mapped while 3C builds on it; 3.14 then adds the bridges.
14. **Shared behaviours are tested once per method in `tests/da/test_design.py`
    and `tests/da/test_schema.py`** (parametrized over `METHODS`): the Phase 3
    review focus's 3.4/3.5 test names live there, not in `test_linda.py` and
    `test_ancombc.py`.
15. **mypy config:** `skbio.stats.composition._ancombc.ancombc2` joins
    `untyped_calls_exclude` (configuration only, no dependency).

Edits outside this section that the plan commit should make: in "Review
focus" item 4, the tests are `tests/da/test_design.py::test_reference_sets_the_sign`
and `test_missing_group_value_raises` (both methods); in item 5, 3.4's
`test_qvalue_is_benjamini_hochberg` is `tests/da/test_schema.py`'s (both
methods); design notes 7 and 10 carry slice 3B decisions 1, 2 and 12 if approved.

### Slice 3B self-review
1. **Brief coverage.** Task order and boundaries (design, "Execution
   order"); the schema module and validator contract; the "agree" Decision
   text (3.8 Step 5); LinDA as MicrobiomeStat computes it, mode estimator and
   bandwidth included, verified against R (parity above); `pseudocounted`
   (slice 3B decision 9); ANCOM-BC2's argument mapping, defaults, BH, ln -> log2,
   untested rows, densify note and golden with measured tolerance (slice 3B
   decisions 2-4); `group`/`covariates`/`reference` (`da/_design.py`, no patsy import);
   consensus columns, validation, NaN semantics and bounds; `pl.consensus`'s
   dot matrix; the image additions with versions and determinism; docs (one
   guide page; 3D's method pages); Checkpoint B with knowledge, push/PR/merge
   and the question before 3C.
2. **Placeholder scan.** Every step carries the full file or the exact diff,
   rendered from the scratch clone's commits that passed every gate; no
   "TBD". The Decision concept's `generated.at` and `commit` are the only
   values the implementer fills in.
3. **Type consistency.** `model(adata, group, *, covariates, reference,
   func) -> (frame, contrast)`, `design_matrix(frame, *, scale)`,
   `dense_counts(adata, *, func)`, `result(features, *, effect, se, pvalue,
   method, contrast)`, `validate_result(table, *, arg)`; public signatures as
   Phase 3 decisions 7, 9 and 10 give them; consensus columns as in the
   design.
4. **Review focus.** Every test named in the slice 3B review focus exists in
   the rendered code (checked by searching this section for each name); the
   Phase 3 review focus's 3B names exist with the locations given above.
5. **Known residual risks and what was not verified.**
   - ANCOM-BC2 parity rests on one dataset and two models; the near-exact
     `log_depth` model and the shifted `host` model show the E-M sensitivity
     is data-dependent. Another dataset could need wider tolerances.
   - scikit-bio fails on degenerate tables (zero coefficient variances) with
     an unhelpful `ValueError: Quantiles must be in the range [0, 1]`; R stops
     with "Zero variances have been detected". biotapy passes the error
     through (not silent); a clearer message would be a small wrapper.
   - `_core.require_counts` accepted negative whole numbers; LinDA would then
     take the log of a negative value (NumPy warns, results are NaN). Fixed at
     the root by Task 3.B0 (approved by the user 2026-10-05 as its own
     `fix(core)` commit), before Task 3.5.
   - LinDA's mean shift guards a mode of exactly 0 (where modeest divides by
     zero and R errors); only Hypothesis reaches it, with identical
     coefficients.
   - CVXR 1.0-15 is built from the CRAN archive at image build time; if CRAN
     removed the archive URL the build would fail loudly (guard layer).
   - One unexplained `1 failed` in a loaded full run at the 3.8 commit (see
     "How slice 3B was checked"). Ruled at approval (2026-10-05): the
     consensus property test carries `@settings(deadline=None)`, as the
     schema property does.
   - When gating a commit after a later one in the same worktree, delete
     `docs/generated/` first: stale autosummary stubs of later functions make
     `sphinx-build -W` fail (seen once on the replay, not a code problem).

## Slice 3C - R bridges

**Goal:** with R, the R package and `pip install 'biotapy[r]'`, a user runs ALDEx2
and MaAsLin 3 from Python and gets the same seven-column table as `da.linda` and
`da.ancombc2`; with the same seed, the numbers equal R's own call digit for digit;
without R, the call stops with an `ImportError` that says what to install; CI runs
the bridges on Linux, including the four-method consensus on the exit-gate data.

**How slice 3C was checked.** Every file below was written into a scratch clone
of the repository at `e9ec091` (`phase-3c`, after slice 3B and the R6.2 wording
commit) and replayed as one commit per step that says "commit" (seven commits).
Each committed state was gated; the `r` column is `pytest -m r` inside a local
image built `FROM` the golden image of that commit's `tests/r/Dockerfile`, plus
Python 3.13.16 and uv (tagged `biotapy-p3c-py`, never pushed), with the clone
mounted:

| Commit | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` | `pytest -q -m r` (R image) |
|---|---|---|---|---|---|
| 3.6 `build(r)` (ALDEx2) `f10f61f` | passed | 1174 passed, 29 deselected | 36 passed | build succeeded | none selected (no `r` test yet) |
| 3.6 `test(golden)` (ALDEx2) `b0519fe` | passed | 1175 passed, 29 deselected | 36 passed | build succeeded | none selected |
| 3.6 `feat(da)` (ALDEx2, extra `r`) `83f9422` | passed | 1190 passed, 35 deselected | 36 passed | build succeeded | 6 passed |
| 3.7 `build(r)` (maaslin3) `bda991c` | passed | 1190 passed, 35 deselected | 36 passed | build succeeded | 6 passed |
| 3.7 `test(golden)` (MaAsLin 3) `1da6805` | passed | 1191 passed, 35 deselected | 36 passed | build succeeded | 6 passed |
| 3.7 `feat(da)` (MaAsLin 3) `e8f483f` | passed | 1203 passed, 44 deselected | 36 passed | build succeeded | 15 passed |
| 3.11 `ci` (`r-bridge`) `da686e3` | passed | 1206 passed, 45 deselected | 36 passed | build succeeded | 16 passed |

The master count before the slice is 1174 passed, 29 deselected (Checkpoint B);
the hashes are the scratch clone's, not the repository's.

- prek covers ruff 0.16.9 check and format, `mypy --strict`, import-linter,
  pyproject-fmt and zizmor (which checked the new job). Every run exported
  `BIOTAPY_DATA_DIR` to a scratch pooch cache; `~/.cache/biotapy` was checked
  absent afterwards. `docs/generated/` was deleted before each docs build.
- Coverage on the final state, without R (`coverage run -m pytest`, the default
  marker set): `da/_r.py`, `da/_aldex2.py` and `da/_maaslin3.py` are each 100%, as is every other `da` file; the package is at 99% (1206 passed, 45 deselected under coverage). The no-R tests reach every line because they replace
  only the rpy2 modules (design, "Tests").
- The RED results in each "expect failure" step were reproduced by running the
  task's new tests against the previous commit.
- The golden image was rebuilt from each edited `tests/r/Dockerfile` (tagged
  `biotapy-golden-p3c-aldex2` after the ALDEx2 commit and `biotapy-golden-p3c`
  after the maaslin3 commit). A rebuild from the committed maaslin3 Dockerfile
  hit the cache with the same image id; a rebuild from the committed ALDEx2
  Dockerfile recompiled the ANCOMBC and ALDEx2 layers and regenerated all 41
  files byte for byte. Layer times: ALDEx2 122 s (MatrixGenerics, GenomicRanges,
  S4Arrays, SparseArray, DelayedArray, BiocParallel, SummarizedExperiment and
  ALDEx2 compile from source; the CRAN dependencies are P3M binaries), maaslin3
  74 s (SingleCellExperiment, TreeSummarizedExperiment, maaslin3; lme4 was already
  there for ANCOMBC). No apt package was needed. `export_golden.R` ran twice on
  each image: bit-identical (41 files, then 42), and every older file unchanged.
  The ALDEx2 section prints `Warning: stack imbalance in '::', 10 then 12` (ALDEx2's
  namespace loading inside a call; harmless, and the reason the bridge loads
  packages with `importr`, design).
- **Parity measured** (the 636 GlobalPatterns genera of slice 3B's goldens, 9
  human samples against 17 others):
  - `da.aldex2(benchmark, "host", reference="human", seed=20260927)` against
    `set.seed(1165433077); ALDEx2::aldex(counts, host, mc.samples = 128, test = "t",
    effect = TRUE, denom = "all")`: `effect` vs `diff.btw` largest absolute
    difference 5.3e-15 (relative 4.2e-15), `pvalue` vs `we.ep` 5.0e-16 (relative
    8.4e-13). Exact: 1165433077 is the integer biotapy derives from seed 20260927,
    and `reference="human"` keeps R's sorted level order, so the Monte Carlo draws
    are the same [corrected 2026-10-07: the Dirichlet draws do not depend on the labels; only `aldex.effect`'s
    resampling of `diff.btw` follows the label order]. With other seeds (1, 2, 3) the Monte Carlo spread is: Spearman
    correlation of effects 0.993, signs agree on 612-617 of 636 genera, largest
    effect difference 0.64-0.83 log2, p-values within 0.12, 11-12 calls at
    q < 0.05 (the same calls on 635-636 genera). With `reference="other"` (`diff.btw` resampled in
    the other order [corrected 2026-10-07: was "draws in the other order"]) Spearman 0.995 against the negated golden. Calls: BH of
    `we.ep` 11 genera, ALDEx2's own `we.eBH` 19 (slice 3C decision 3).
  - `da.maaslin3(benchmark, "host", covariates, reference="other",
    seed=20260927)` against `set.seed(1165433077); maaslin3::maaslin3(...,
    evaluate_only = "abundance", warn_prevalence = FALSE, subtract_median = TRUE)`,
    models `host` and `host + log_depth`: the same 36 genera untested (fit errors,
    R's p-value NA); `effect` vs `coef` 5.3e-15 (relative 2.1e-14), `se` vs
    `stderr` 5.3e-15, `pvalue` vs `pval_individual` 5.0e-16 (relative 7.2e-13).
    Exact, for the same reason (the median test's 10,000 normal draws). Other
    seeds change only the p-values (largest change 0.0015; 52 calls either way).
    Calls at q < 0.05: 52 (`host`) and 45 (`host + log_depth`) with BH over the
    group's p-values; MaAsLin 3's `qval_individual` gives 52 and 20, because it
    pools the covariate's p-values (slice 3C decision 4).
  - `subtract_median = TRUE` (slice 3C decision 5): the median coefficient is
    -1.70 log2; the effect's sign differs from the raw coefficient's for 109 of
    600 tested genera and for none of the 52 called; signs agree with LinDA's on
    484 of 600 tested genera (437 with the raw coefficient) and on 50 of the 50
    genera both call.
  - **Four methods** (`reference="other"`, `seed=0` for the bridges; 15 s in
    all): calls 208 (ANCOM-BC2), 118 (LinDA), 13 (ALDEx2), 52 (MaAsLin 3);
    `n_tested` 4 for 600 genera and 2 for 36; `n_significant` 0/1/2/3/4 for
    413/114/62/35/12 genera; consensus 223, 109, 47 and 12 genera at
    `min_methods` 1, 2, 3 and 4; no conflict; the 12 four-way calls are all up in
    human hosts.
- **rpy2, measured.** `uv lock` resolves rpy2 3.6.8 with rpy2-rinterface 3.6.7
  and rpy2-robjects 3.6.5 (new transitive packages: cffi, jinja2, tzlocal;
  `uv.lock` is git-ignored, so nothing else is committed). rpy2-rinterface is an
  sdist on Linux and links `R CMD config --ldflags`, which on R 4.5.3 lists
  `-lpcre2-8 -ldeflate -lzstd -llzma -lbz2 -licuuc ...`. Without those `-dev`
  packages its API-mode build fails at the linker and it silently installs the
  ABI mode instead, whose import then fails with `libR.so: cannot open shared
  object file` unless `LD_LIBRARY_PATH` holds R's `lib`. With the packages it
  builds in API mode and imports cleanly; with `RPY2_CFFI_MODE=API` and without
  them the install fails loudly (`Failed to build rpy2-rinterface==3.6.7 ...
  cannot find -lpcre2-8`). Without R, `pip install rpy2==3.6.8` fails ("rpy2 in
  API mode cannot be built without R in the PATH or R_HOME defined"); with R
  removed after the build, `import rpy2.robjects` raises `ModuleNotFoundError`,
  an `ImportError`, so `import_optional` turns it into the "pip install
  'biotapy[r]'" message with the cause chained (resolves design note 8's
  [UNVERIFIED]). `py.typed` ships only in `rpy2/rinterface`,
  `rpy2/rinterface_lib`, `rpy2/rlike` and `rpy2/situation`, not in
  `rpy2/robjects` or `rpy2/` (resolves the [UNVERIFIED] in the global
  constraints; slice 3C decision 11).
- **CI install, simulated.** In a fresh `rocker/r-ver:4.5.3` container on the
  image's dated CRAN snapshot, pak (stable) resolves Bioconductor 3.22 for R 4.5.3
  (`pak::repo_status()`), and `pak::pkg_install(c("bioc::ALDEx2",
  "bioc::maaslin3"), dependencies = NA)` installed 98 packages (116.8 MB) in
  2 min 40 s, 20 of them Bioconductor packages built from source, plus the system
  packages `libpng-dev cmake libuv1-dev libicu-dev`; the versions equal the
  image's. Expected `r-bridge` job time: about 5 minutes cold (R, pak, packages,
  rpy2's build, the tests: 30 s for the `r` tests in the container) and about
  2 minutes with the package cache; [UNVERIFIED] on GitHub's runners.
- APIs checked in the installed versions, R side read at source in the image:
  ALDEx2 1.42.0 (`aldex(reads, conditions, mc.samples = 128, test = "t", effect
  = TRUE, include.sample.summary = FALSE, verbose = FALSE, paired.test = FALSE,
  denom = "all", iterate = FALSE, gamma = NULL)` runs `aldex.clr`, `aldex.ttest`,
  `aldex.effect` in that order; `aldex.clr` removes features with no read, adds
  0.5, fails on a factor `conds` (`'round' not meaningful for factors`) and warns
  below 128 draws; `aldex.ttest`/`aldex.effect` take `as.factor(conditions)`,
  so the levels are sorted, and `diff.btw` is the second level minus the first;
  `aldex.effect` stops when a level has fewer than two samples; `we.ep` is the
  per-draw two-sided Welch p averaged [corrected 2026-10-07: per draw the one-sided Welch p is doubled in each direction and capped at 1, each
  direction is averaged over the draws, and the smaller average is kept], `we.eBH` the per-draw BH values averaged;
  output columns `rab.all, rab.win.<l1>, rab.win.<l2>, diff.btw, diff.win,
  effect, overlap, we.ep, we.eBH, wi.ep, wi.eBH`). maaslin3 1.2.0
  (`maaslin3(input_data, input_metadata, output, formula, ..., min_abundance =
  0, min_prevalence = 0, normalization = "TSS", transform = "LOG", standardize =
  TRUE, median_comparison_abundance = TRUE, subtract_median = FALSE,
  warn_prevalence = TRUE, evaluate_only = NULL, plot_summary_plot = TRUE,
  plot_associations = TRUE, cores = 1, verbosity = "FINEST")`; `evaluate_only`
  with `warn_prevalence = TRUE` stops; `LOG` is `log2` of the non-zero values; a
  factor column not named in `reference` keeps its levels, the first being the
  reference; `reference` is a `"var,level;..."` string split on `,` and `;`;
  numeric columns are z-scored; the median comparison draws `rnorm` 10,000
  times per term; `subtract_median = TRUE` reports `coef - median` with the same
  p-values; `add_qvals` applies BH over every row of both models, setting rows
  with an `error` to NA; with `evaluate_only = "abundance"`, `pval_individual =
  pval`; the return value's `fit_data_abundance$results` holds `feature,
  metadata, value, name, coef, null_hypothesis, stderr, pval_individual,
  qval_individual, pval_joint, qval_joint, error, model, N, N_not_zero`; the
  output folder is written even with plots off; fits go through `pbapply`, whose
  progress bar ignores `verbosity`). rpy2 3.6.8 (`rpy2.robjects.r(code)`,
  `packages.importr`, `packages.PackageNotInstalledError`, an `ImportError`
  subclass; `default_converter + pandas2ri.converter` and `.context()`;
  `pandas2ri` alone does not convert NumPy arrays (`NotImplementedError`); a
  `pd.Series` of strings becomes a character vector; a `DataFrame` keeps its row
  and column names unmangled; a `Categorical` becomes a factor only when its
  categories are strings, otherwise rpy2 warns and sends strings; R `NA`
  character comes back as rpy2's `NA_character_` object, not NaN; R's
  `message()` output goes to the `rpy2.rinterface_lib.callbacks` logger at
  WARNING). r-lib/actions v2 is v2.14.0 at
  `f9a764fea8d5c63df6ef9a5c7795bf7deb5d7e05` (`gh api
  repos/r-lib/actions/git/ref/tags/v2`); `setup-r` takes `r-version`,
  `use-public-rspm`, `cran`; `setup-r-dependencies` takes `packages` (default
  `deps::., any::sessioninfo`, which needs a DESCRIPTION file biotapy does not
  have), `extra-packages`, `dependencies` (an R expression, default `"all"`) and
  caches the library keyed on its lockfile.

### Slice 3C design

- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `tests/r/Dockerfile` | ALDEx2 (3.6), maaslin3 (3.7), one `build(r)` commit each |
  | `tests/r/export_golden.R` | the ALDEx2 (3.6) and MaAsLin 3 (3.7) sections, seeded with the integer biotapy derives |
  | `tests/golden/global_patterns/{aldex2,maaslin3}.csv.gz` | the goldens (19.4 KB, 43.4 KB) |
  | `pyproject.toml` | `optional-dependencies.r = [ "rpy2>=3.6.8" ]` (3.6); nothing for mypy |
  | `da/_r.py` | `r_seed`, `r_function`, `call_r`: the only code that touches rpy2 (3.6) |
  | `da/_aldex2.py` | `aldex2`, `_conditions` (3.6) |
  | `da/_maaslin3.py` | `maaslin3`, `_string_levels` (3.7) |
  | `tests/da/conftest.py` | the `fake_rpy2` fixture (3.6) |
  | `tests/da/test_aldex2.py`, `test_maaslin3.py` | no-R tests and `r` tests per bridge |
  | `tests/da/test_aldex2_golden.py`, `test_maaslin3_golden.py` | the goldens (`r`) |
  | `tests/da/test_consensus.py` | the four-method consensus (`r`, 3.11) |
  | `.github/workflows/test.yaml`, `tests/test_ci.py` | the `r-bridge` job (3.11) |

  `_r.py` serves the two `da` bridges only, so it stays in `da` (R4.3). The
  shared helpers `_design.model`, `_design.dense_counts` and `_schema.result` are
  reused. [Corrected by the fix rounds: they were not unchanged. `dense_counts`
  gained the repeated `var_names`/`obs_names` check, and `_design._column` now
  returns unordered categoricals, which fixed an ordered group's effect for all
  four methods.]
- **The bridge module (`da/_r.py`).** Three functions, each a few lines:
  - `r_seed(seed) -> int`: `int(as_generator(seed).integers(2**31 - 1))`, one
    draw, a non-negative 32-bit integer for `set.seed` (R3.4).
  - `r_function(code, *, package, func)`: imports `rpy2.robjects` and
    `rpy2.robjects.packages` through `import_optional(..., extra="r")`, loads the
    package with `importr` (a `PackageNotInstalledError` becomes `ImportError:
    da.aldex2 needs the R package ALDEx2. Install it in R with:
    BiocManager::install("ALDEx2")`), and returns the R function `code` defines.
    `importr` first matters: a namespace that `::` loads in the middle of a call
    made from rpy2 makes R print "stack imbalance" warnings (measured).
  - `call_r(function, *args)`: calls it inside a local `(default_converter +
    pandas2ri.converter).context()`, never rpy2's global activation, so pandas
    arguments become R vectors and data frames and the R data frame it returns a
    pandas one.
  Each bridge holds its R code as a short closure string that calls `set.seed`
  first and returns a data frame, so all mapping to the schema is Python.
- **ALDEx2.** `aldex2(adata, group, *, mc_samples=128, reference=None,
  seed=None)` (roadmap 3.6's signature). `model(..., covariates=())` checks the
  group; `_conditions` refuses a numeric group (ALDEx2's t-test has no slope)
  and a level in fewer than two samples (ALDEx2 stops there), then sends the
  levels as `"0"` (reference) and `"1"`: ALDEx2 breaks on a factor and sorts
  character labels, which depends on R's locale for level names. `mc_samples`
  is an int >= 1. `X` is densified once (features as rows, as ALDEx2 wants).
  The R closure runs `set.seed(seed)` and `suppressMessages(ALDEx2::aldex(reads,
  conditions, mc.samples, test = "t", effect = TRUE, denom = "all"))`.
  Mapping: `effect = diff.btw` (log2), `se` NaN, `pvalue = we.ep`, rows reindexed
  to `var_names` (ALDEx2 drops a feature with no read: NaN, untested), then
  `_schema.result` computes BH over the finite p-values.
- **MaAsLin 3.** `maaslin3(adata, group, *, covariates=(), reference=None,
  seed=None)` (roadmap 3.7's signature). The model frame goes to R with its
  columns renamed `x0, x1, ...` and the formula `"~ x0 + x1"` (no quoting for
  names such as `"body site"`, as `da.ancombc2` does for patsy); categorical
  columns go as factors with string levels in `model()`'s order, so the group's
  first level is the reference and no `reference` string is built (its
  `"var,level"` syntax breaks on a comma). Bool levels are renamed to strings
  first: rpy2 sends non-string categories as plain strings with a UserWarning,
  maaslin3 then sorts them, and `reference="True"` would silently flip the sign
  (found and fixed in the prototype; review focus 2). The R closure calls
  `maaslin3(counts, metadata, tempfile(), formula, min_abundance = 0,
  min_prevalence = 0, evaluate_only = "abundance", warn_prevalence = FALSE,
  subtract_median = TRUE, plot_summary_plot = FALSE, plot_associations = FALSE,
  cores = 1, verbosity = "ERROR")` with pbapply's progress bar off, deletes the
  folder on exit, and returns `feature, metadata, coef, stderr, pval, failed =
  !is.na(error)` from `fit_data_abundance$results`. Mapping: the rows of
  `metadata == "x0"`, reindexed to `var_names`; a failed fit or a missing row is
  untested (NaN); `effect = coef` (already minus the median), `se = stderr`,
  `pvalue = pval_individual`; BH by `_schema.result` over the group's rows.
- **Seeds.** Both bridges call `set.seed(r_seed(seed))` inside the R closure,
  right before the method. ALDEx2 is random throughout; MaAsLin 3 only in its
  median test (10,000 normal draws per term), so its `effect` and `se` do not
  depend on the seed and its p-values move by about 0.001. The goldens use
  `set.seed(1165433077)`, the integer `seed=20260927` becomes, so the bridges
  are compared elementwise.
- **Tests.** Without R (every CI leg): the `fake_rpy2` fixture puts stand-in
  modules for `rpy2`, `rpy2.robjects`, `rpy2.robjects.packages` and
  `rpy2.robjects.pandas2ri` into `sys.modules` (monkeypatch), whose R functions
  record their arguments and return a canned data frame. biotapy's own code runs
  unchanged: argument checks, labels, orientation, seed derivation, the schema
  mapping and both `ImportError` messages (a missing rpy2 is simulated by `None`
  entries in `sys.modules`, so the test also holds where rpy2 is installed). The
  canned tables are R's real output on `toy()`; an `r` test checks each against
  R. With R (marker `r`): toy tables, seeds, reference swaps, an all-zero
  feature, schema validity through `bt.da.consensus([table], min_methods=1)`, the
  goldens and the four-method consensus. Bridge tests carry `r` only, never
  `golden` or `network` (slice 3C decision 1).
- **CI job `r-bridge`.** Ubuntu 24.04 (the CRAN snapshot's binaries are for
  noble), `r-lib/actions/setup-r` with R 4.5.3 and CRAN at the image's P3M
  snapshot, `setup-r-dependencies` with `packages: bioc::ALDEx2, bioc::maaslin3`
  and hard dependencies only, the six `-dev` packages rpy2's link needs, uv with
  Python 3.13, the pooch cache step of the network job, then `uv run --group test
  --extra r pytest -m r` with `RPY2_CFFI_MODE=API`. Added to `check.needs`.
- **Docs.** The DA guide gains "Methods that run in R" (install, licences, the
  `ImportError`, densifying), with an ALDEx2 (3.6) and a MaAsLin 3 (3.7)
  subsection; `api.md` gains both; `contributing.md` gains "R bridge tests".
  Slice 3D's method pages and tutorial are unchanged by 3C.

### Slice 3C global constraints (in addition to the Phase 3 list)
- rpy2 is reached only through `import_optional` inside `da/_r.py`'s functions;
  no module imports rpy2 at load time, so `import biotapy` and the
  `import-without-extras` job need no R.
- Bridge tests carry `r` (and only `r`); the default run and `-m "golden or
  network"` never need R. Tests reach the bridges through `bt.da.aldex2` and
  `bt.da.maaslin3`; the only stand-in is `fake_rpy2`, at the rpy2 boundary.
- Every bridge golden is seed-matched and compared at the contract's default
  `rtol=1e-7`.
- Run commands: host as in slice 3B; the `r` tests with `uv run --group test
  --extra r pytest -m r` where R and the packages are installed (the prototype
  used the image above with `BIOTAPY_DATA_DIR` mounted).
- Commits stage explicit paths only (never `.claude/`, `.superpowers/`,
  `.worktrees/`, `notebooks/`, `build/`). Each task's last commit stages
  `.knowledge/roadmap/phase-3-stats.md` with its box ticked and
  `.knowledge/log.md` with its line, under the heading `## <date of the commit>
  (Phase 3, slice 3C)` at the top of the log (create it if the date differs).

### Slice 3C review focus
The five ways a user is most likely to get a wrong answer from slice 3C without
an error, each pinned by a named test:

1. **The sign is the other way round in R.** ALDEx2 sorts its labels; maaslin3
   sorts character levels. Expected: `reference` decides the sign in both
   bridges, in any R locale. Tests: `tests/da/test_aldex2.py::test_aldex2_sends_features_as_rows_the_reference_as_0_and_one_seed`,
   `test_reference_sets_the_sign_in_r`;
   `tests/da/test_maaslin3.py::test_maaslin3_gets_counts_plain_column_names_a_formula_and_one_seed`,
   `test_reference_sets_the_sign_in_r`.
2. **A bool group reaches R as unordered strings.** Expected: levels are sent as
   strings in the reference's order. Tests:
   `test_maaslin3.py::test_bool_levels_reach_r_as_strings_in_reference_order`,
   `test_bool_reference_sets_the_sign_in_r`.
3. **R's own q-values mixed with biotapy's.** ALDEx2's `we.eBH` averages per-draw
   corrections; MaAsLin 3's `qval_individual` pools covariates. Expected: BH over
   the group's tested features, as for the native methods. Tests:
   `test_aldex2.py::test_aldex2_maps_aldex_onto_the_schema`,
   `test_maaslin3.py::test_maaslin3_maps_the_group_rows_onto_the_schema`, and
   both goldens.
4. **A feature R could not fit counted as tested.** Expected: ALDEx2's dropped
   rows and MaAsLin 3's fit errors are NaN and left out of BH. Tests:
   `test_aldex2.py::test_feature_aldex2_drops_is_not_tested`,
   `test_all_zero_feature_is_not_tested_in_r`;
   `test_maaslin3.py::test_failed_fit_is_not_tested`,
   `test_feature_without_a_row_is_not_tested`,
   `test_all_zero_feature_is_not_tested_in_r`,
   `tests/da/test_maaslin3_golden.py` (the 36 genera).
5. **An unreproducible or silently different R run.** Expected: one `seed`
   gives one table, the goldens match R digit for digit, and a missing rpy2 or R
   package raises naming the fix. Tests: `test_seed_reproduces_the_table_in_r`
   (both), `test_seed_accepts_a_generator`,
   `test_missing_rpy2_names_the_extra` and
   `test_missing_r_package_names_the_install_line` (both),
   `tests/da/test_aldex2_golden.py::test_aldex2_matches_aldex2_aldex`,
   `test_maaslin3_golden.py::test_maaslin3_matches_maaslin3_maaslin3`.

Plus purity (R3.3): `test_aldex2_keeps_input`, `test_maaslin3_keeps_input`.

Execution order: **3.6 -> 3.7 -> 3.11 -> Checkpoint C.** 3.6 brings the extra,
`_r.py` and the fixture that 3.7 reuses; the CI job runs both bridges and the
consensus over all four methods.

---
### Task 3.6: `da.aldex2` and the extra `r`

**Files:** modify `tests/r/Dockerfile`, `tests/r/export_golden.R`,
`tests/golden/VERSIONS.txt`, `pyproject.toml`, `src/biotapy/da/__init__.py`,
`tests/da/conftest.py`, `docs/guide/differential_abundance.md`, `docs/api.md`,
`docs/contributing.md`, `.knowledge/contracts/r-golden-parity.md`; create
`tests/golden/global_patterns/aldex2.csv.gz`, `src/biotapy/da/_r.py`,
`src/biotapy/da/_aldex2.py`, `tests/da/test_aldex2.py`,
`tests/da/test_aldex2_golden.py`.
**Not touched:** `da/_design.py`, `da/_schema.py`, `da/_linda.py`,
`da/_ancombc.py`, `da/_consensus.py`, `_core` (`import_optional` and
`as_generator` are used as they are), `[tool.mypy]`, `[tool.pytest]` (the `r`
marker is already registered), `.github/`.
**Interfaces:**
- Consumes: `_core.import_optional`, `_core.as_generator`; `model`,
  `dense_counts` (3.5); `result` (3.5); the `benchmark` fixture.
- Produces: `bt.da.aldex2(adata: AnnData, group: str, *, mc_samples: int = 128,
  reference: str | None = None, seed: int | np.random.Generator | None = None)
  -> pd.DataFrame`; private `da._r.r_seed(seed) -> int`,
  `r_function(code: str, *, package: str, func: str) -> Callable[..., Any]`,
  `call_r(function, /, *args: object) -> pd.DataFrame`; the `fake_rpy2` test
  fixture (`.output`, `.installed`, `.calls`); `aldex2.csv.gz`; the extra `r`.

- [x] **Step 1: Approval on record.** The extra `r = ["rpy2>=3.6.8"]` (Phase 3
  decision 14) and ALDEx2 in the image (Phase 3 decision 6) were approved with
  the Phase 3 plan. rpy2 3.6.8 brings rpy2-rinterface 3.6.7, rpy2-robjects
  3.6.5, cffi, jinja2 and tzlocal (all installed only with the extra). Nothing
  new to ask.
- [x] **Step 2: Add ALDEx2 to the image.** `tests/r/Dockerfile`:
  ```diff
  @@ -28,6 +28,8 @@ RUN apt-get update && apt-get install -y --no-install-recommends libgsl27 && rm
   RUN Rscript -e 'install.packages(c("Rmpfr", "gmp", "ECOSolveR", "scs", "osqp", "bit64", "cli", "Rcpp", "RcppEigen"))' \
       && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/CVXR/CVXR_1.0-15.tar.gz", repos = NULL, type = "source")'
   RUN Rscript -e 'BiocManager::install("ANCOMBC", version = "3.22", ask = FALSE, update = FALSE)'
  -RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE), requireNamespace("ANCOMBC", quietly = TRUE))'
  +# ALDEx2 (Bioconductor, GPL-3): the da.aldex2 golden file and the bridge's r tests; biotapy calls it through rpy2.
  +RUN Rscript -e 'BiocManager::install("ALDEx2", version = "3.22", ask = FALSE, update = FALSE)'
  +RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE), requireNamespace("ANCOMBC", quietly = TRUE), requireNamespace("ALDEx2", quietly = TRUE))'
   WORKDIR /work
   CMD ["Rscript", "tests/r/export_golden.R"]
  ```
  `.knowledge/contracts/r-golden-parity.md` statement 1:
  ```diff
  @@ -22,8 +22,9 @@ sources:
      Bioconductor `philr` for `pp.philr`, with the `libuv1` runtime library its
      `fs` binary loads, CRAN `MicrobiomeStat` for `da.linda`, and Bioconductor
      `ANCOMBC` for `da.ancombc2`, with CRAN's archived CVXR 1.0-15 and the
  -   `libgsl27` runtime library it needs); a new golden function that needs another
  -   package adds it in its own commit (rules.md R2.3).
  +   `libgsl27` runtime library it needs, and Bioconductor `ALDEx2` for
  +   `da.aldex2`); a new golden function that needs another package adds it in
  +   its own commit (rules.md R2.3).
   1b. HUMAnN golden files (`fn.func_glom`, `fn.renorm`) are produced by
      `tests/humann/export_golden.py`, run with
      `uv run --no-project --with humann==3.9 --with pandas==3.0.6`, never
  ```
  Log line (top of the slice 3C heading):
  ```text
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `ALDEx2` for the `da.aldex2` golden file.
  ```
  Build: `docker build -t biotapy-golden tests/r`. Expected: the guard layer
  passes; `docker run --rm biotapy-golden Rscript -e 'packageVersion("ALDEx2")'`
  prints `1.42.0` (zCompositions 1.6.1, Rfast 2.1.5.2 as P3M binaries); every
  earlier package unchanged. The ALDEx2 layer takes about two minutes (eight
  Bioconductor packages compile).
- [x] **Step 3: Commit the image change on its own:**
  ```bash
  git add tests/r/Dockerfile .knowledge/contracts/r-golden-parity.md .knowledge/log.md
  git commit -m "build(r): add ALDEx2 to the golden image

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
- [x] **Step 4: Export.** In `tests/r/export_golden.R`, after the ANCOM-BC2
  section, and `VERSIONS.txt` gains ALDEx2. The seed is the integer
  `bt.da.aldex2(..., seed=20260927)` passes to R; check it before running:
  `uv run python -c "import numpy as np; print(np.random.default_rng(20260927).integers(2**31 - 1))"`
  prints `1165433077` (NumPy 2.5.3).
  ```diff
  @@ -173,6 +173,22 @@ ancombc_rows <- lapply(da_formulas, function(f) {
   })
   write_golden(do.call(rbind, ancombc_rows), file.path(gp, "ancombc2.csv.gz"))

  +## Slice 3C golden files: the R bridges, on the same data
  +# ALDEx2 is Monte Carlo. 1165433077 is the integer bt.da.aldex2(..., seed=20260927) passes to set.seed
  +# (np.random.default_rng(20260927).integers(2**31 - 1)), so the bridge test can compare digits, not only ranks.
  +# The integer comes from NumPy's Generator stream: if a NumPy release changes it, test_aldex2_golden.py fails loudly;
  +# recompute the integer with that expression, put it here and in the test comment, and re-export.
  +# The conditions sort "human" before "other", so diff.btw is other - human; the test sets reference="human".
  +set.seed(1165433077)
  +aldex_out <- suppressMessages(ALDEx2::aldex(
  +  da_counts, as.character(da_meta$host), mc.samples = 128, test = "t", effect = TRUE, denom = "all"
  +))
  +write_golden(
  +  data.frame(taxon_id = rownames(aldex_out), diff_btw = aldex_out$diff.btw, we_ep = aldex_out$we.ep,
  +             we_eBH = aldex_out$we.eBH),
  +  file.path(gp, "aldex2.csv.gz")
  +)
  +
   ## Synthetic phyloseq fixtures: biotapy's toy() numbers, no third-party data
   counts <- rbind(
     c(10, 5, 20, 30, 0, 2, 1, 0), c(8, 7, 25, 22, 3, 0, 0, 1), c(12, 4, 18, 35, 1, 5, 2, 0),
  @@ -238,5 +252,6 @@ writeLines(c(
     paste0("MicrobiomeStat ", packageVersion("MicrobiomeStat")),
     paste0("modeest ", packageVersion("modeest")),
     paste0("ANCOMBC ", packageVersion("ANCOMBC")),
  -  paste0("CVXR ", packageVersion("CVXR"), " (CRAN archive, pinned in tests/r/Dockerfile)")
  +  paste0("CVXR ", packageVersion("CVXR"), " (CRAN archive, pinned in tests/r/Dockerfile)"),
  +  paste0("ALDEx2 ", packageVersion("ALDEx2"))
   ), "tests/golden/VERSIONS.txt")
  ```
  Run twice and check bit identity as in 3.5 Step 4. Expected: every line `OK`
  (41 files; each run about three minutes); `git status` shows `M
  tests/golden/VERSIONS.txt` (one new line, `ALDEx2 1.42.0`), `M
  tests/r/export_golden.R` and the new `aldex2.csv.gz` (19.4 KB, 636 rows,
  columns `taxon_id, diff_btw, we_ep, we_eBH`). Besides the known lines the run
  prints `Warning: stack imbalance in '::', 10 then 12` and `... in '<-', 2 then
  4` while ALDEx2's namespace loads; harmless.
- [x] **Step 5: Gate and commit.** `uv run --group test pytest
  tests/test_data_files.py -q` -> passes; `uvx prek run --all-files`;
  ```bash
  git add tests/r/export_golden.R tests/golden/VERSIONS.txt tests/golden/global_patterns/aldex2.csv.gz
  git commit -m "test(golden): export the ALDEx2 golden file

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
- [x] **Step 6: Failing tests.** `tests/da/conftest.py` gains the fixture:
  ```diff
  @@ -1,3 +1,7 @@
  +import contextlib
  +import sys
  +from types import ModuleType, SimpleNamespace
  +
   import numpy as np
   import pandas as pd
   import pytest
  @@ -14,3 +18,56 @@ def benchmark():
       tdata.obs["host"] = pd.Categorical(np.where(human, "human", "other"))
       tdata.obs["log_depth"] = np.log(np.asarray(tdata.X.sum(axis=1)).ravel())
       return tdata
  +
  +
  +class _Converter:
  +    """rpy2's converter as the bridge uses it: ``a + b``, then ``.context()``."""
  +
  +    def __add__(self, other):
  +        return self
  +
  +    def context(self):
  +        return contextlib.nullcontext()
  +
  +
  +@pytest.fixture
  +def fake_rpy2(monkeypatch):
  +    """The rpy2 modules the R bridges import, without R: every R function returns ``fake.output`` and logs its call.
  +
  +    The bridges' own code (seeds, labels, orientation, the schema mapping) runs unchanged; only rpy2 is replaced.
  +    """
  +    fake = SimpleNamespace(output=None, installed={"ALDEx2", "maaslin3"}, calls=[], warnings=[])
  +
  +    def r(code):
  +        def function(*args):
  +            fake.calls.append((code, args))
  +            if isinstance(fake.output, Exception):
  +                raise fake.output
  +            # What r_function's R wrapper returns after rpy2's conversion: the value and the warning messages.
  +            return SimpleNamespace(values=lambda: (fake.output, fake.warnings))
  +
  +        return function
  +
  +    embedded = ModuleType("rpy2.rinterface_lib.embedded")
  +    # rpy2's class for an R error; the bridge turns it into a RuntimeError.
  +    embedded.RRuntimeError = type("RRuntimeError", (RuntimeError,), {})
  +    robjects = ModuleType("rpy2.robjects")
  +    robjects.r = r
  +    robjects.default_converter = _Converter()
  +    packages = ModuleType("rpy2.robjects.packages")
  +    # rpy2's own class is an ImportError subclass too (rpy2.robjects.packages.LibraryError).
  +    packages.PackageNotInstalledError = type("PackageNotInstalledError", (ImportError,), {})
  +
  +    def importr(name):
  +        if name not in fake.installed:
  +            raise packages.PackageNotInstalledError(f'The R package "{name}" is not installed.')
  +
  +    packages.importr = importr
  +    pandas2ri = ModuleType("rpy2.robjects.pandas2ri")
  +    pandas2ri.converter = _Converter()
  +    modules = {"rpy2": ModuleType("rpy2"), "rpy2.robjects": robjects}
  +    modules |= {"rpy2.robjects.packages": packages, "rpy2.robjects.pandas2ri": pandas2ri}
  +    modules["rpy2.rinterface_lib.embedded"] = embedded
  +    for name, module in modules.items():
  +        monkeypatch.setitem(sys.modules, name, module)
  +    return fake
  ```
  `tests/da/test_aldex2.py` (the canned table is R's output on `toy()` with the
  seed `seed=0` becomes, `1826701614`):
  ```python
  import sys
  import warnings

  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from scipy.stats import false_discovery_control

  import biotapy as bt

  # ALDEx2::aldex(t(toy counts), rep(c("0", "1"), each = 3), mc.samples = 128, ...) after set.seed(1826701614), the
  # integer bt.da.aldex2(toy, "group", seed=0) passes to R (6 significant digits; test_aldex2_on_toy_is_the_canned_table checks it in R).
  CANNED = pd.DataFrame(
      {
          "diff.btw": [-3.70251, -2.8971, -2.66143, -1.72769, -0.198703, 4.06776, 3.66243, 2.69778],
          "effect": [-1.69089, -1.24091, -2.09973, -1.91244, -0.0639279, 2.14194, 1.99029, 1.03629],
          "we.ep": [0.08828, 0.125246, 0.0447655, 0.0520091, 0.736825, 0.074296, 0.0783649, 0.172146],
          "we.eBH": [0.219275, 0.269506, 0.154564, 0.161972, 0.880257, 0.28763, 0.296191, 0.46741],
      },
      index=[f"f{i}" for i in range(1, 9)],
  )


  def test_aldex2_maps_aldex_onto_the_schema(fake_rpy2):
      fake_rpy2.output = CANNED.iloc[::-1]  # R's row order does not matter
      out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
      assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
      np.testing.assert_array_equal(out["effect"], CANNED["diff.btw"])
      np.testing.assert_array_equal(out["pvalue"], CANNED["we.ep"])
      assert out["se"].isna().all()
      # BH of the expected p-values, as for every method; ALDEx2's we.eBH averages per-draw corrections instead.
      np.testing.assert_allclose(out["qvalue"], false_discovery_control(CANNED["we.ep"]), rtol=1e-12)
      assert (out["method"] == "aldex2").all() and (out["contrast"] == "B vs A").all()


  def test_aldex2_sends_features_as_rows_the_reference_as_0_and_one_seed(fake_rpy2):
      fake_rpy2.output = CANNED
      bt.da.aldex2(bt.datasets.toy(), "group", mc_samples=64, reference="B", seed=7)
      ((code, (reads, conditions, mc_samples, seed)),) = fake_rpy2.calls
      assert "set.seed(seed)" in code and "ALDEx2::aldex(" in code
      assert reads.shape == (8, 6) and reads.index.tolist() == bt.datasets.toy().var_names.tolist()
      np.testing.assert_array_equal(reads.to_numpy(), bt.datasets.toy().X.toarray().T)
      assert conditions.tolist() == ["1"] * 3 + ["0"] * 3  # reference B is "0", which ALDEx2 sorts first
      assert mc_samples == 64 and seed == int(np.random.default_rng(7).integers(2**31 - 1))


  def test_seed_accepts_a_generator(fake_rpy2):
      fake_rpy2.output = CANNED
      bt.da.aldex2(bt.datasets.toy(), "group", seed=np.random.default_rng(3))
      bt.da.aldex2(bt.datasets.toy(), "group", seed=3)
      assert fake_rpy2.calls[0][1][3] == fake_rpy2.calls[1][1][3]
      generator = np.random.default_rng(3)
      bt.da.aldex2(bt.datasets.toy(), "group", seed=generator)
      bt.da.aldex2(bt.datasets.toy(), "group", seed=generator)
      assert fake_rpy2.calls[2][1][3] != fake_rpy2.calls[3][1][3]  # the caller's generator advances


  def test_feature_aldex2_drops_is_not_tested(fake_rpy2):
      fake_rpy2.output = CANNED.drop(index="f7")  # ALDEx2 removes a feature with no read before its draws
      out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
      assert out.loc["f7", ["effect", "pvalue", "qvalue"]].isna().all() and out.loc["f7", "direction"] == 0
      np.testing.assert_allclose(out["qvalue"].drop("f7"), false_discovery_control(CANNED["we.ep"].drop("f7")))


  def test_missing_rpy2_names_the_extra(monkeypatch):
      for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
          monkeypatch.setitem(sys.modules, name, None)
      with pytest.raises(ImportError, match=r"pip install 'biotapy\[r\]'"):
          bt.da.aldex2(bt.datasets.toy(), "group")


  def test_missing_r_package_names_the_install_line(fake_rpy2):
      fake_rpy2.installed = set()
      with pytest.raises(
          ImportError, match=r'da\.aldex2 needs the R package ALDEx2\. .*BiocManager::install\("ALDEx2"\)'
      ):
          bt.da.aldex2(bt.datasets.toy(), "group")


  @pytest.mark.parametrize(("axis", "fix"), [("var", "var_names_make_unique"), ("obs", "obs_names_make_unique")])
  def test_repeated_names_raise(fake_rpy2, axis, fix):
      tdata = bt.datasets.toy()
      names = getattr(tdata, f"{axis}_names").tolist()
      setattr(tdata, f"{axis}_names", ["x", "x", *names[2:]])
      with pytest.raises(ValueError, match=rf"unique {axis} names.*\['x'\].*adata\.{fix}\(\)"):
          bt.da.aldex2(tdata, "group")
      assert fake_rpy2.calls == []


  def test_seed_is_checked_before_rpy2_is_imported(monkeypatch):
      for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
          monkeypatch.setitem(sys.modules, name, None)
      with pytest.raises(TypeError, match="seed must be"):
          bt.da.aldex2(bt.datasets.toy(), "group", seed="abc")


  def test_each_r_warning_is_re_emitted_at_the_callers_line(fake_rpy2):
      fake_rpy2.output = CANNED
      fake_rpy2.warnings = ["values are unreliable", "second warning"]
      with pytest.warns(UserWarning) as record:
          bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
      assert [str(w.message) for w in record] == [
          "da.aldex2: R warned: values are unreliable",
          "da.aldex2: R warned: second warning",
      ]
      assert record[0].filename == __file__


  def test_r_errors_name_the_function(fake_rpy2):
      fake_rpy2.output = sys.modules["rpy2.rinterface_lib.embedded"].RRuntimeError("Error in f() : boom\n")
      with pytest.raises(RuntimeError, match=r"da\.aldex2: R stopped: Error in f\(\) : boom") as caught:
          bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
      assert isinstance(caught.value.__cause__, sys.modules["rpy2.rinterface_lib.embedded"].RRuntimeError)


  def test_numeric_group_raises(fake_rpy2):
      tdata = bt.datasets.toy()
      tdata.obs["ph"] = [5.1, 5.6, 6.0, 6.8, 7.1, 7.4]
      with pytest.raises(ValueError, match=r"compares two levels, but obs\['ph'\] is numeric"):
          bt.da.aldex2(tdata, "ph")


  def test_level_in_one_sample_raises(fake_rpy2):
      tdata = bt.datasets.toy()
      tdata.obs["arm"] = ["x"] * 5 + ["y"]
      with pytest.raises(ValueError, match=r"needs two samples in each level of obs\['arm'\]"):
          bt.da.aldex2(tdata, "arm")


  @pytest.mark.parametrize(("mc_samples", "error"), [(1.5, TypeError), (True, TypeError), (0, ValueError)])
  def test_mc_samples_must_be_a_positive_int(fake_rpy2, mc_samples, error):
      with pytest.raises(error, match="mc_samples"):
          bt.da.aldex2(bt.datasets.toy(), "group", mc_samples=mc_samples)


  def test_shared_model_checks_apply(fake_rpy2):
      tdata = bt.datasets.toy()
      tdata.obs["site"] = ["a", "b", "c", "a", "b", "c"]
      with pytest.raises(ValueError, match=r"compares two levels, but obs\['site'\] has 3"):
          bt.da.aldex2(tdata, "site")
      dense = tdata.X.toarray()
      dense[1] = 0
      tdata.X = sp.csr_matrix(dense)
      with pytest.raises(ValueError, match=r"1 sample\(s\) have none"):
          bt.da.aldex2(tdata, "group")
      assert fake_rpy2.calls == []


  def test_aldex2_keeps_input(fake_rpy2, assert_unchanged):
      fake_rpy2.output = CANNED
      tdata = bt.datasets.toy()
      before = tdata.copy()
      bt.da.aldex2(tdata, "group", seed=0)
      assert_unchanged(before, tdata)


  @pytest.mark.r
  def test_aldex2_on_toy_is_the_canned_table():
      out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
      # CANNED holds 6 significant digits of R's output, hence rtol=1e-5.
      np.testing.assert_allclose(out["effect"], CANNED["diff.btw"], rtol=1e-5)
      np.testing.assert_allclose(out["pvalue"], CANNED["we.ep"], rtol=1e-5)


  @pytest.mark.r
  def test_seed_reproduces_the_table_in_r():
      tdata = bt.datasets.toy()
      first = bt.da.aldex2(tdata, "group", seed=0)
      pd.testing.assert_frame_equal(bt.da.aldex2(tdata, "group", seed=0), first)
      assert not bt.da.aldex2(tdata, "group", seed=1)["effect"].equals(first["effect"])


  @pytest.mark.r
  def test_reference_sets_the_sign_in_r():
      tdata = bt.datasets.toy()
      a_first, b_first = bt.da.aldex2(tdata, "group", seed=0), bt.da.aldex2(tdata, "group", reference="B", seed=0)
      assert (b_first["contrast"] == "A vs B").all()
      # Swapping the reference reorders ALDEx2's draws, so the effects flip sign but are not exactly opposite.
      assert (a_first["direction"] == -b_first["direction"]).all()
      np.testing.assert_allclose(b_first["effect"], -a_first["effect"], atol=0.7)


  @pytest.mark.r
  def test_all_zero_feature_is_not_tested_in_r():
      tdata = bt.datasets.toy()
      dense = tdata.X.toarray()
      dense[:, tdata.var_names.get_loc("f7")] = 0
      tdata.X = sp.csr_matrix(dense)
      out = bt.da.aldex2(tdata, "group", seed=0)
      assert out.loc["f7", ["effect", "pvalue", "qvalue"]].isna().all()
      assert np.isfinite(out.drop(index="f7")[["effect", "pvalue", "qvalue"]].to_numpy()).all()


  @pytest.mark.r
  def test_few_mc_samples_warn_in_python():
      with pytest.warns(UserWarning, match=r"da\.aldex2: R warned: .*unreliable"):
          bt.da.aldex2(bt.datasets.toy(), "group", mc_samples=16, seed=0)


  @pytest.mark.r
  def test_r_messages_are_not_warnings():
      # A direct test of the bridge: R's message() is progress text, not a warning condition.
      from biotapy.da._r import call_r, r_function

      function = r_function('function() { message("hi"); data.frame(a = 1) }', package="base", func="da.x")
      with warnings.catch_warnings():
          warnings.simplefilter("error")
          assert call_r(function)["a"].tolist() == [1.0]


  @pytest.mark.r
  def test_aldex2_table_passes_the_consensus_checks():
      out = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
      assert bt.da.consensus([out], min_methods=1)["n_tested"].tolist() == [1] * 8
  ```
  `tests/da/test_aldex2_golden.py`:
  ```python
  from pathlib import Path

  import numpy as np
  import pandas as pd
  import pytest
  from scipy.stats import false_discovery_control

  import biotapy as bt

  GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
  # Marker r only: it needs R, and the r-bridge CI job gives it the network job's pooch cache.
  pytestmark = pytest.mark.r


  def test_aldex2_matches_aldex2_aldex(benchmark):
      # R: set.seed(1165433077); ALDEx2::aldex(counts, as.character(host), mc.samples = 128, test = "t", effect = TRUE,
      # denom = "all"). 1165433077 is the integer seed=20260927 becomes, and reference="human" keeps R's sorted level
      # order, so the Monte Carlo draws are the same and the numbers match to rounding. The integer is NumPy's
      # default_rng(20260927).integers(2**31 - 1): a NumPy change to that stream fails this test loudly; recompute it and
      # re-export the golden (tests/r/export_golden.R).
      golden = pd.read_csv(GOLDEN / "aldex2.csv.gz", dtype={"taxon_id": str}).set_index("taxon_id")
      out = bt.da.aldex2(benchmark, "host", reference="human", seed=20260927)
      assert out.index.tolist() == golden.index.tolist()
      assert (out["contrast"] == "other vs human").all()
      np.testing.assert_allclose(out["effect"], golden["diff_btw"], rtol=1e-7)
      np.testing.assert_allclose(out["pvalue"], golden["we_ep"], rtol=1e-7)
      np.testing.assert_allclose(out["qvalue"], false_discovery_control(golden["we_ep"]), rtol=1e-7)
  ```
- [x] **Step 7: Run, expect failure** -
  `uv run --group test pytest tests/da/test_aldex2.py -q` -> `13 failed, 5
  deselected` (`AttributeError: module 'biotapy.da' has no attribute
  'aldex2'`). Where R and rpy2 are installed, the eight `r` tests (seven in
  `test_aldex2.py`, the golden) fail with the same `AttributeError`.
- [x] **Step 8: Implement.** `pyproject.toml` (pyproject-fmt keeps the extra
  after `dependencies`):
  ```diff
  @@ -44,6 +44,8 @@ dependencies = [
     "treedata>=0.3.1,<0.4",
     "xarray",
   ]
  +# The ALDEx2 and MaAsLin 3 bridges (decisions/optional-heavy-dependencies); R and the R packages are the user's install.
  +optional-dependencies.r = [ "rpy2>=3.6.8" ]
   # https://docs.pypi.org/project_metadata/#project-urls
   urls.Documentation = "https://biotapy.readthedocs.io/"
   urls.Homepage = "https://github.com/pedrocr83/biotapy"
  ```
  `src/biotapy/da/_r.py`:
  ```python
  """The rpy2 bridge for the da methods only R implements (decisions/r-bridge-before-ports); extra ``r``."""

  import warnings
  from collections.abc import Callable
  from typing import Any, cast

  import numpy as np
  import pandas as pd

  from biotapy._core import as_generator, import_optional

  # R's set.seed takes a 32-bit signed integer; drawing below 2**31 - 1 keeps it non-negative and never R's NA.
  _SEED_BOUND = 2**31 - 1


  def r_seed(seed: int | np.random.Generator | None) -> int:
      """The integer for R's ``set.seed``, drawn once from ``as_generator(seed)`` (rules.md R3.4)."""
      return int(as_generator(seed).integers(_SEED_BOUND))


  # Runs the function and returns its value with the messages of the R warnings it raised, muffled so that R does not
  # print them again; message() output is not a warning condition and keeps rpy2's default handling.
  _CATCH_WARNINGS = """function(f) function(...) {
    caught <- character()
    value <- withCallingHandlers(f(...), warning = function(cond) {
      caught <<- c(caught, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = caught)
  }"""


  def r_function(code: str, *, package: str, func: str) -> Callable[..., Any]:
      """The R function ``code`` defines; ImportError naming the extra without rpy2, or the install line without ``package``.

      Calling it returns the R function's value, re-emits each R warning as a ``UserWarning`` and turns an R error into a
      ``RuntimeError``, both naming ``func``.
      """
      robjects = import_optional("rpy2.robjects", extra="r")
      rpackages = import_optional("rpy2.robjects.packages", extra="r")
      embedded = import_optional("rpy2.rinterface_lib.embedded", extra="r")
      try:
          # Loaded before the call: a package that `::` loads during it makes R print "stack imbalance" warnings.
          rpackages.importr(package)
      except rpackages.PackageNotInstalledError as err:
          msg = f'{func} needs the R package {package}. Install it in R with: BiocManager::install("{package}")'
          raise ImportError(msg) from err
      function = robjects.r(f"({_CATCH_WARNINGS})({code})")

      def run(*args: object) -> object:
          try:
              value, caught = function(*args).values()
          except embedded.RRuntimeError as err:
              msg = f"{func}: R stopped: {str(err).strip()}"
              raise RuntimeError(msg) from err
          for message in caught:
              warnings.warn(f"{func}: R warned: {message}", stacklevel=4)
          return value

      return run


  def call_r(function: Callable[..., Any], /, *args: object) -> pd.DataFrame:
      """``function(*args)``: pandas arguments become R vectors and data frames, the R data frame it returns pandas."""
      robjects = import_optional("rpy2.robjects", extra="r")
      pandas2ri = import_optional("rpy2.robjects.pandas2ri", extra="r")
      # A local converter, never rpy2's global activation, so the user's own rpy2 session keeps its conversion rules.
      with (robjects.default_converter + pandas2ri.converter).context():
          # ModuleType attributes are Any, so the result needs a cast under mypy --strict (warn_return_any).
          return cast("pd.DataFrame", function(*args))
  ```
  `src/biotapy/da/_aldex2.py`:
  ```python
  """ALDEx2 through the R bridge: Monte Carlo CLR values of two groups compared by Welch's t-test."""

  import numpy as np
  import pandas as pd
  from anndata import AnnData

  from ._design import dense_counts, model
  from ._r import call_r, r_function, r_seed
  from ._schema import result

  # aldex() is ALDEx2's documented entry point and the golden file's call; its progress notes are R messages, which
  # rpy2 would log as warnings.
  _ALDEX = """function(reads, conditions, mc_samples, seed) {
    set.seed(seed)
    suppressMessages(ALDEx2::aldex(reads, conditions, mc.samples = mc_samples, test = "t", effect = TRUE, denom = "all"))
  }"""


  def aldex2(
      adata: AnnData,
      group: str,
      *,
      mc_samples: int = 128,
      reference: str | None = None,
      seed: int | np.random.Generator | None = None,
  ) -> pd.DataFrame:
      """Differential abundance of each feature between two groups by ALDEx2, run in R.

      Parameters
      ----------
      adata
          Samples x features; ``X`` holds raw counts.
      group
          The ``obs`` column whose two levels are compared: a categorical, string or
          bool column with two levels, each in at least two samples.
      mc_samples
          Monte Carlo draws from each sample's Dirichlet posterior (``mc.samples``).
      reference
          The level of ``group`` that the other level is compared with; by default
          its first category (sorted values for a string column). It changes the
          Monte Carlo draws as well as the sign (Notes).
      seed
          Seeds R's random number generator through ``set.seed``; the same seed
          gives the same table.

      Returns
      -------
      pandas.DataFrame
          One row per feature, in ``var_names`` order, indexed by ``feature``:
          ``effect`` (ALDEx2's ``diff.btw``, the median log2 difference between the
          groups), ``se`` (NaN: ALDEx2 has none), ``pvalue`` (``we.ep``),
          ``qvalue`` (Benjamini-Hochberg over the tested features), ``direction``,
          ``method`` (``"aldex2"``) and ``contrast``. A feature with no read in any
          sample is not tested: NaN values and ``direction`` 0.

      Raises
      ------
      ImportError
          rpy2 (the extra ``r``) or the R package ALDEx2 is not installed.
      KeyError
          ``group`` is not an ``obs`` column.
      TypeError
          ``reference`` is not a string, or ``mc_samples`` not an int.
      RuntimeError
          R stops with an error.
      ValueError
          ``X`` does not hold raw counts, has an empty sample, repeated
          ``var_names`` or ``obs_names``, or fewer than two features; ``group`` is numeric, has missing values or other than two
          levels, or a level in fewer than two samples; ``reference`` is not one of
          them; ``mc_samples`` is below 1.

      Notes
      -----
      R equivalent: ``ALDEx2::aldex``
      Guide: :doc:`/guide/differential_abundance`

      Calls ``ALDEx2::aldex(reads, conditions, mc.samples, test = "t", effect =
      TRUE, denom = "all")`` (Fernandes et al. 2014) through rpy2: Dirichlet Monte
      Carlo draws of each sample's proportions (0.5 added to every count), their
      log2 centred log-ratios, and per draw a Welch t-test, whose two-sided
      p-values are averaged into ``we.ep`` [corrected 2026-10-07: one-sided p doubled per direction, averaged over the draws, smaller kept; see the corrected docstring]. ``effect`` is ``diff.btw``, already
      log2; ALDEx2's own ``effect`` column is a standardised size, not a fold
      change, and is not carried. ``qvalue`` is the Benjamini-Hochberg correction
      of ``we.ep``, as for every method; ALDEx2's ``we.eBH`` averages the
      corrected values of the draws instead and is not carried. ALDEx2 compares
      two groups without covariates (its ``glm`` test is not wrapped).

      Swapping ``reference`` does more than flip the sign: ALDEx2 takes its Monte
      Carlo draws in label order [corrected 2026-10-07: only `diff.btw`'s resampling does; the draws and `we.ep` do not], so the same ``seed`` gives different effects
      (by up to 0.25 log2 on ``toy()``, whose effects are about 4 log2 wide), and
      the same p-values on ``toy()``. Each R warning raised during the call, such as
      the one for fewer than 128 ``mc_samples``, is re-emitted as a
      ``UserWarning``, and an R error is raised as a ``RuntimeError``.

      Needs R with ALDEx2 (``BiocManager::install("ALDEx2")``) and ``pip install
      'biotapy[r]'``, which builds rpy2 (GPL-2.0-or-later) against that R.
      ``seed`` becomes one integer for R's ``set.seed``, which sets the random
      state of the R session rpy2 embeds in Python.

      rpy2 converts dense tables only, so ``X`` is densified once: 8 bytes x
      samples x features, plus R's copy and its Monte Carlo draws (about
      ``mc_samples`` times the table).

      References
      ----------
      Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB (2014)
      Unifying the analysis of high-throughput sequencing datasets: characterizing
      RNA-seq, 16S rRNA gene sequencing and selective growth experiments by
      compositional data analysis. Microbiome 2:15.

      Fernandes AD, Macklaim JM, Linn TG, Reid G, Gloor GB (2013) ANOVA-like
      differential gene expression analysis of single-organism and meta-RNA-seq.
      PLoS ONE 8:e67019.

      Examples
      --------
      >>> import biotapy as bt
      >>> table = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)  # doctest: +SKIP
      >>> table.loc["f6", "contrast"]  # doctest: +SKIP
      'B vs A'
      """
      frame, contrast = model(adata, group, covariates=(), reference=reference, func="da.aldex2")
      conditions = _conditions(frame[group], func="da.aldex2")
      if isinstance(mc_samples, bool) or not isinstance(mc_samples, int | np.integer):
          msg = f"da.aldex2: mc_samples must be an int, got {type(mc_samples).__name__}"
          raise TypeError(msg)
      if mc_samples < 1:
          msg = f"da.aldex2: mc_samples must be at least 1, got {mc_samples}"
          raise ValueError(msg)
      seed_for_r = r_seed(seed)
      # rpy2 has no sparse converter (rules.md R6.2): the one dense copy of X, features as rows as ALDEx2 wants them.
      reads = pd.DataFrame(dense_counts(adata, func="da.aldex2").T, index=adata.var_names, columns=adata.obs_names)
      aldex = r_function(_ALDEX, package="ALDEx2", func="da.aldex2")
      table = call_r(aldex, reads, conditions, int(mc_samples), seed_for_r).reindex(adata.var_names)
      effect = table["diff.btw"].to_numpy(np.float64)
      pvalue = table["we.ep"].to_numpy(np.float64)
      se = np.full(adata.n_vars, np.nan)
      return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="aldex2", contrast=contrast)


  def _conditions(values: pd.Series, *, func: str) -> pd.Series:
      """The group as "0" (reference) and "1" labels; raises for a numeric group or a level in fewer than two samples."""
      if not isinstance(values.dtype, pd.CategoricalDtype):
          msg = f"{func} compares two levels, but obs[{values.name!r}] is numeric; ALDEx2's t-test has no slope"
          raise ValueError(msg)
      sizes = values.value_counts()
      if (sizes < 2).any():
          msg = f"{func} needs two samples in each level of obs[{values.name!r}], which has {sizes.to_dict()}"
          raise ValueError(msg)
      # ALDEx2 sorts the labels and reports the second minus the first; "0" < "1" in every R locale, unlike level names.
      return values.cat.codes.astype(str)
  ```
  `src/biotapy/da/__init__.py`:
  ```diff
  @@ -1,7 +1,8 @@
   """Differential abundance: methods that share one result table (contracts/data-model-slots, DA results)."""

  +from ._aldex2 import aldex2
   from ._ancombc import ancombc2
   from ._consensus import consensus
   from ._linda import linda

  -__all__ = ["ancombc2", "consensus", "linda"]
  +__all__ = ["aldex2", "ancombc2", "consensus", "linda"]
  ```
  No mypy override: biotapy reaches rpy2 only through `import_optional`, whose
  `ModuleType` attributes are `Any` (hence the two `cast`s), and the lint
  environment has no rpy2 (it cannot build without R). Slice 3C decision 11.
- [x] **Step 9: Run, expect pass** - `uv run --group test pytest tests/da -q`
  -> `126 passed, 12 deselected`; in the R environment (`uv sync --group test
  --extra r` with R 4.5.3 and ALDEx2 installed) `uv run --group test --extra r
  pytest -m r tests/da/test_aldex2.py tests/da/test_aldex2_golden.py -q` ->
  `8 passed` (the golden about 10 s, GlobalPatterns from the pooch cache).
  If rpy2 imports with "Error importing in API mode", the R link headers are
  missing: install `libpcre2-dev libdeflate-dev libzstd-dev` and rebuild it with
  `RPY2_CFFI_MODE=API` (design, CI job).
- [x] **Step 10: Docs.**
  ````diff
  @@ -86,6 +86,38 @@ here. On the GlobalPatterns genera (`host`) the calls at q < 0.05 go from 208 to
   plus its swap is about -0.38 log2, not 0; `da.consensus` with LinDA goes from 104 to 112 genera.
   Choose `reference` on the biology (the control or baseline level), not to change the results.

  +## Methods that run in R
  +
  +ALDEx2 and MaAsLin 3 exist only in R, so biotapy calls them there through
  +[rpy2](https://rpy2.github.io/). They need R, the R package, and biotapy's `r` extra, which builds
  +rpy2 against that R (Linux and macOS; it fails to install when no R is found):
  +
  +```bash
  +Rscript -e 'install.packages("BiocManager"); BiocManager::install("ALDEx2")'
  +pip install 'biotapy[r]'
  +```
  +
  +rpy2 is GPL-2.0-or-later and the R packages have their own licences; biotapy itself does not ship
  +any of them. Without rpy2 or the R package, the call raises an `ImportError` that names what to
  +install. Each call converts `X` to a dense table once, because rpy2 has no sparse converter.
  +
  +### ALDEx2
  +
  +`bt.da.aldex2` runs `ALDEx2::aldex` (Fernandes et al. 2014): Monte Carlo draws from each sample's
  +Dirichlet posterior, their log2 centred log-ratios, and a Welch t-test per draw, whose p-values are
  +averaged. `effect` is ALDEx2's `diff.btw`, the median log2 difference between the two groups:
  +
  +```python
  +table = bt.da.aldex2(tdata, "group", seed=0)
  +```
  +
  +ALDEx2 compares two groups without covariates, and each group needs two samples. It is random:
  +`seed` sets R's random state, so the same seed gives the same table, and on the GlobalPatterns
  +genera biotapy's numbers equal R's `set.seed(...); aldex(...)` to a relative 8.4e-13 (when `reference` is R's
  +first sorted level, and R is seeded with the integer biotapy derives from `seed`). `qvalue` is the
  +Benjamini-Hochberg correction of ALDEx2's expected p-value `we.ep`, as for every method; ALDEx2's own
  +`we.eBH` averages the corrections of the draws instead and calls more features (19 against 11 on
  +those genera).
  +
  +Swapping `reference` does more than flip the sign: ALDEx2 takes its Monte Carlo draws in label order [corrected 2026-10-07: only `diff.btw`'s resampling does],
  +so the same `seed` gives different effects (up to 0.25 log2 apart on `toy()`, where they are about
  +4 log2 wide), with the same p-values there. Anything R prints during the call, such as the warning
  +for fewer than 128 `mc_samples`, is re-emitted as a Python `UserWarning`; an R error is raised as a
  +`RuntimeError`. Repeated `var_names` or `obs_names` raise: call `adata.var_names_make_unique()` first.
  +
   ## Where methods agree

   `bt.da.consensus` puts the tables of several methods side by side and counts, for each feature,
  ````
  ```diff
  @@ -100,6 +100,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
   .. autosummary::
       :toctree: generated

  +    da.aldex2
       da.ancombc2
       da.consensus
       da.linda
  ```
  ````diff
  @@ -67,6 +67,16 @@ per-user cache directory - CI caches it across runs the same way:
   BIOTAPY_DATA_DIR=.pooch uv run --group test pytest -m "network or golden"
   ```

  +### R bridge tests
  +
  +Tests that call R through rpy2 (`bt.da.aldex2`) carry the marker `r` and are excluded from the
  +runs above. They need R with the packages `tests/r/Dockerfile` installs and the `r` extra; those
  +that read GlobalPatterns also need the pooch cache:
  +
  +```bash
  +BIOTAPY_DATA_DIR=.pooch uv run --group test --extra r pytest -m r
  +```
  +
   ### Regenerating the R golden files

   The golden CSVs under `tests/golden/` and the R-written fixtures under
  ````
- [x] **Step 11: Contract.** `.knowledge/contracts/r-golden-parity.md`
  statement 4 and "Enforced by":
  ```diff
  @@ -52,6 +52,7 @@ sources:
      | HUMAnN parity (func_glom, renorm) | elementwise, matched by row id | rtol=1e-7; renorm rtol=5e-6, because humann_renorm_table prints %.6g |
      | DA methods | sign agreement and rank correlation of effect sizes; exact match where the R method is deterministic, or Monte Carlo and seeded with the integer biotapy derives from its `seed` | per method |
      | `da.linda` vs `MicrobiomeStat::linda(is.winsor = FALSE)` (deterministic) | `effect`, `se`, `pvalue`, `qvalue` elementwise, matched by taxon | `rtol=1e-7` |
  +   | `da.aldex2` vs `ALDEx2::aldex` (Monte Carlo; the golden's `set.seed` is the integer biotapy derives from `seed=20260927`, and `reference="human"` keeps R's level order, so the draws are the same) | `effect` vs `diff.btw` and `pvalue` vs `we.ep` elementwise, matched by taxon; `qvalue` vs BH of `we.ep` | `rtol=1e-7` (measured 4e-15 and 8e-13); 1165433077 is `np.random.default_rng(20260927).integers(2**31 - 1)`, so a NumPy change to that stream fails the test loudly: recompute the integer and re-export |
      | `da.ancombc2` vs `ANCOMBC::ancombc2` (deterministic, but its bias E-M can stop at 100 iterations before converging, on a slightly different iterate in scikit-bio) | the same untested features; `effect`, `se`, `pvalue` elementwise; Spearman correlation of effects; the same calls at `q < 0.05`, with R's p-values corrected over the tested features | `host` model: `effect` atol 0.015 (log2), `se` rtol 2e-3, `pvalue` atol 0.02, Spearman > 0.9999; `host + log_depth`: all three at 1e-6, Spearman > 0.999999 |

   5. Any looser tolerance is written in the test with a one-line comment giving the reason.
  @@ -94,6 +95,10 @@ invariants keeps them honest.
     They also carry `network`, because the golden inputs (e.g. GlobalPatterns)
     are downloaded by pooch; they run in the network CI job
     (`pytest -m "network or golden"`).
  +- `tests/da/test_*_golden.py` of the R bridges, marker `r` only: they need R
  +  and rpy2 and run in the `r-bridge` CI job, which gives them the network job's
  +  pooch cache; marked `golden` or `network`, they would run in the network job,
  +  which has no R.
   - `tests/fn/*_golden.py`, marker `golden` only: their inputs are committed, so
     they run in every CI job.
   - Missing golden file = test error, not skip.
  ```
  Log line:
  ```text
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md): `da.aldex2` is compared with `ALDEx2::aldex` elementwise (the golden's seed is the one biotapy derives, so the Monte Carlo draws match); R bridge golden tests carry the marker `r` only. [phase-3-stats](roadmap/phase-3-stats.md) task 3.6 done.
  ```
- [x] **Step 12: Gate and commit** (tick 3.6 in the checklist first)
  ```bash
  git add pyproject.toml src/biotapy/da/__init__.py src/biotapy/da/_r.py src/biotapy/da/_aldex2.py \
    tests/da/conftest.py tests/da/test_aldex2.py tests/da/test_aldex2_golden.py docs/guide/differential_abundance.md \
    docs/api.md docs/contributing.md .knowledge/contracts/r-golden-parity.md .knowledge/roadmap/phase-3-stats.md \
    .knowledge/log.md
  uvx prek run --all-files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git commit -m "feat(da): add the ALDEx2 bridge and the r extra

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```

- [x] **Step 13: Fix round 1** (review of 0485b86; `fix(da)`). Repeated `var_names` or
  `obs_names` raise in `da._design.dense_counts` (all methods; `test_design.py` and
  `test_aldex2.py::test_repeated_names_raise`); `r_function` wraps the R function in `withCallingHandlers` and returns a closure that
  re-emits each R warning condition as its own `UserWarning` at the caller's line and turns an R error into a
  `RuntimeError`, both naming the function (`test_each_r_warning_is_re_emitted_at_the_callers_line`,
  `test_r_errors_name_the_function`, and in R `test_few_mc_samples_warn_in_python`,
  `test_r_messages_are_not_warnings`); `seed` is checked before rpy2 loads
  (`test_seed_is_checked_before_rpy2_is_imported`); the docstring, guide, golden test, export script and
  contract say that `reference` changes the Monte Carlo draws and that 1165433077 depends on NumPy's
  stream. Expected: `pytest tests/da -q` -> `126 passed, 12 deselected`; `-m r` -> `8 passed`.

---

### Task 3.7: `da.maaslin3`

**Files:** modify `tests/r/Dockerfile`, `tests/r/export_golden.R`,
`tests/golden/VERSIONS.txt`, `src/biotapy/da/__init__.py`,
`docs/guide/differential_abundance.md`, `docs/api.md`, `docs/contributing.md`,
`.knowledge/contracts/r-golden-parity.md`; create
`tests/golden/global_patterns/maaslin3.csv.gz`, `src/biotapy/da/_maaslin3.py`,
`tests/da/test_maaslin3.py`, `tests/da/test_maaslin3_golden.py`.
**Not touched:** `da/_r.py`, `da/_aldex2.py`, `tests/da/conftest.py` (the
fixture already lists `maaslin3` as installed), the shared helpers, `.github/`.
**Interfaces:**
- Consumes: `r_function`, `call_r`, `r_seed`, `fake_rpy2` (3.6); `model`,
  `dense_counts`, `result`; the `benchmark` fixture.
- Produces: `bt.da.maaslin3(adata: AnnData, group: str, *, covariates:
  Sequence[str] = (), reference: str | None = None, seed: int |
  np.random.Generator | None = None) -> pd.DataFrame`; `maaslin3.csv.gz`.

- [x] **Step 1: Approval on record.** maaslin3 (Bioconductor, MIT) in the image
  was approved with the Phase 3 plan. New here: `VERSIONS.txt` also records
  multcomp, whose `glht` computes the median test's p-values (a record, not a
  pin). The behavioural choices (`subtract_median = TRUE`, BH over the group's
  rows) are slice 3C decisions 4 and 5.
- [x] **Step 2: Add maaslin3 to the image.** `tests/r/Dockerfile`:
  ```diff
  @@ -30,6 +30,8 @@ RUN Rscript -e 'install.packages(c("Rmpfr", "gmp", "ECOSolveR", "scs", "osqp", "
   RUN Rscript -e 'BiocManager::install("ANCOMBC", version = "3.22", ask = FALSE, update = FALSE)'
   # ALDEx2 (Bioconductor, GPL-3): the da.aldex2 golden file and the bridge's r tests; biotapy calls it through rpy2.
   RUN Rscript -e 'BiocManager::install("ALDEx2", version = "3.22", ask = FALSE, update = FALSE)'
  -RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE), requireNamespace("ANCOMBC", quietly = TRUE), requireNamespace("ALDEx2", quietly = TRUE))'
  +# maaslin3 (Bioconductor, MIT): the da.maaslin3 golden file and the bridge's r tests; biotapy calls it through rpy2.
  +RUN Rscript -e 'BiocManager::install("maaslin3", version = "3.22", ask = FALSE, update = FALSE)'
  +RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE), requireNamespace("ANCOMBC", quietly = TRUE), requireNamespace("ALDEx2", quietly = TRUE), requireNamespace("maaslin3", quietly = TRUE))'
   WORKDIR /work
   CMD ["Rscript", "tests/r/export_golden.R"]
  ```
  ```diff
  @@ -22,9 +22,9 @@ sources:
      Bioconductor `philr` for `pp.philr`, with the `libuv1` runtime library its
      `fs` binary loads, CRAN `MicrobiomeStat` for `da.linda`, and Bioconductor
      `ANCOMBC` for `da.ancombc2`, with CRAN's archived CVXR 1.0-15 and the
  -   `libgsl27` runtime library it needs, and Bioconductor `ALDEx2` for
  -   `da.aldex2`); a new golden function that needs another package adds it in
  -   its own commit (rules.md R2.3).
  +   `libgsl27` runtime library it needs, Bioconductor `ALDEx2` for `da.aldex2`
  +   and Bioconductor `maaslin3` for `da.maaslin3`); a new golden function that
  +   needs another package adds it in its own commit (rules.md R2.3).
   1b. HUMAnN golden files (`fn.func_glom`, `fn.renorm`) are produced by
      `tests/humann/export_golden.py`, run with
      `uv run --no-project --with humann==3.9 --with pandas==3.0.6`, never
  ```
  Log line:
  ```text
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md) statement 1: the golden image also installs Bioconductor `maaslin3` for the `da.maaslin3` golden file.
  ```
  Build. Expected: the guard passes; `maaslin3 1.2.0`, `multcomp 1.4.30`,
  `lme4 2.0.1` (from ANCOMBC, unchanged); the layer takes about 75 s.
- [x] **Step 3: Commit the image change on its own:**
  ```bash
  git add tests/r/Dockerfile .knowledge/contracts/r-golden-parity.md .knowledge/log.md
  git commit -m "build(r): add maaslin3 to the golden image

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
- [x] **Step 4: Export.** After the ALDEx2 section:
  ```diff
  @@ -187,6 +187,24 @@ write_golden(
     file.path(gp, "aldex2.csv.gz")
   )

  +# MaAsLin 3 tests each coefficient against the median over features by simulation (rnorm): the same seed again. Abundance
  +# model only (warn_prevalence must then be FALSE), no prevalence filter, each coefficient minus that median, no plots.
  +maaslin_rows <- lapply(da_formulas, function(f) {
  +  output <- tempfile()
  +  set.seed(1165433077)
  +  fit <- maaslin3::maaslin3(
  +    as.data.frame(t(da_counts)), da_meta, output, formula = paste0("~", f), min_abundance = 0, min_prevalence = 0,
  +    evaluate_only = "abundance", warn_prevalence = FALSE, subtract_median = TRUE, plot_summary_plot = FALSE,
  +    plot_associations = FALSE, cores = 1, verbosity = "ERROR"
  +  )
  +  unlink(output, recursive = TRUE)
  +  out <- fit$fit_data_abundance$results
  +  out <- out[out$metadata == "host", ]
  +  data.frame(formula = f, taxon_id = out$feature, coef = out$coef, stderr = out$stderr, pval = out$pval_individual,
  +             qval = out$qval_individual, error = !is.na(out$error))
  +})
  +write_golden(do.call(rbind, maaslin_rows), file.path(gp, "maaslin3.csv.gz"))
  +
   ## Synthetic phyloseq fixtures: biotapy's toy() numbers, no third-party data
   counts <- rbind(
     c(10, 5, 20, 30, 0, 2, 1, 0), c(8, 7, 25, 22, 3, 0, 0, 1), c(12, 4, 18, 35, 1, 5, 2, 0),
  @@ -253,5 +271,7 @@ writeLines(c(
     paste0("modeest ", packageVersion("modeest")),
     paste0("ANCOMBC ", packageVersion("ANCOMBC")),
     paste0("CVXR ", packageVersion("CVXR"), " (CRAN archive, pinned in tests/r/Dockerfile)"),
  -  paste0("ALDEx2 ", packageVersion("ALDEx2"))
  +  paste0("ALDEx2 ", packageVersion("ALDEx2")),
  +  paste0("maaslin3 ", packageVersion("maaslin3")),
  +  paste0("multcomp ", packageVersion("multcomp"))
   ), "tests/golden/VERSIONS.txt")
  ```
  Run twice. Measured: the script prints no `OK` lines (only R's `phylo` class notes and the known `stack imbalance` warnings; no maaslin3 warning or message), and the second run's output is byte-identical to the first (about three minutes a run, 37 files under `tests/golden/*/`);
  `M tests/golden/VERSIONS.txt` (`maaslin3 1.2.0`, `multcomp 1.4.30`), `M
  tests/r/export_golden.R`, the new `maaslin3.csv.gz` (43.4 KB, 636 rows per
  model, 36 with `error` TRUE in each); every older file byte-identical.
- [x] **Step 5: Gate and commit.** `tests/test_data_files.py` passes; prek;
  ```bash
  git add tests/r/export_golden.R tests/golden/VERSIONS.txt tests/golden/global_patterns/maaslin3.csv.gz
  git commit -m "test(golden): export the MaAsLin 3 golden file

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
- [x] **Step 6: Failing tests.** `tests/da/test_maaslin3.py` (the canned rows
  are what the R closure returns for `toy()` with `seed=0`):
  ```python
  import sys

  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from scipy.stats import false_discovery_control

  import biotapy as bt

  # What biotapy's R call returns for bt.da.maaslin3(toy, "group", seed=0): maaslin3's abundance results (sorted by its
  # q-value) for the group term x0, 6 significant digits; test_maaslin3_on_toy_is_the_canned_table checks it in R.
  CANNED = pd.DataFrame(
      {
          "feature": ["f6", "f7", "f8", "f1", "f3", "f2", "f4", "f5"],
          "metadata": "x0",
          "coef": [4.64115, 4.38614, 2.45111, -1.85211, -1.50564, -1.0101, -0.625606, 0.625606],
          "stderr": [0.452717, 0.372327, 0.667219, 0.378469, 0.320919, 0.507362, 0.187502, 0.994888],
          "pval": [0.0337446, 0.036925, 0.206379, 0.226817, 0.27898, 0.474965, 0.619162, 0.657504],
          "failed": False,
      }
  )
  EXPECTED = CANNED.set_index("feature").reindex(bt.datasets.toy().var_names)


  def test_maaslin3_maps_the_group_rows_onto_the_schema(fake_rpy2):
      covariate = CANNED.assign(metadata="x1", coef=9.0, pval=1e-9)  # rows of a covariate are not the result
      fake_rpy2.output = pd.concat([covariate, CANNED], ignore_index=True)
      out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
      assert out.index.tolist() == bt.datasets.toy().var_names.tolist()
      np.testing.assert_array_equal(out["effect"], EXPECTED["coef"])
      np.testing.assert_array_equal(out["se"], EXPECTED["stderr"])
      np.testing.assert_array_equal(out["pvalue"], EXPECTED["pval"])
      # BH over the group's p-values only; maaslin3's qval_individual would pool them with the covariate's.
      np.testing.assert_allclose(out["qvalue"], false_discovery_control(EXPECTED["pval"]), rtol=1e-12)
      assert (out["method"] == "maaslin3").all() and (out["contrast"] == "B vs A").all()


  def test_maaslin3_gets_counts_plain_column_names_a_formula_and_one_seed(fake_rpy2):
      fake_rpy2.output = CANNED
      tdata = bt.datasets.toy()
      tdata.obs["body site"] = pd.Categorical(["x", "y", "x", "y", "x", "y"])
      bt.da.maaslin3(tdata, "group", covariates=["body site"], reference="B", seed=5)
      ((code, (counts, metadata, formula, seed)),) = fake_rpy2.calls
      assert 'evaluate_only = "abundance"' in code and "subtract_median = TRUE" in code and "set.seed(seed)" in code
      np.testing.assert_array_equal(counts.to_numpy(), tdata.X.toarray())
      assert counts.index.tolist() == tdata.obs_names.tolist() and counts.columns.tolist() == tdata.var_names.tolist()
      assert formula == "~ x0 + x1" and metadata.columns.tolist() == ["x0", "x1"]
      assert metadata["x0"].cat.categories.tolist() == ["B", "A"]  # maaslin3 keeps a factor's first level as reference
      assert seed == int(np.random.default_rng(5).integers(2**31 - 1))


  def test_bool_levels_reach_r_as_strings_in_reference_order(fake_rpy2):
      fake_rpy2.output = CANNED
      tdata = bt.datasets.toy()
      tdata.obs["treated"] = [False] * 3 + [True] * 3
      out = bt.da.maaslin3(tdata, "treated", reference="True", seed=0)
      assert fake_rpy2.calls[0][1][1]["x0"].cat.categories.tolist() == ["True", "False"]
      assert (out["contrast"] == "False vs True").all()


  def test_failed_fit_is_not_tested(fake_rpy2):
      failed = CANNED["feature"] == "f3"
      fake_rpy2.output = CANNED.assign(failed=failed)
      out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
      assert out.loc["f3", ["effect", "se", "pvalue", "qvalue"]].isna().all() and out.loc["f3", "direction"] == 0
      np.testing.assert_allclose(out["qvalue"].drop("f3"), false_discovery_control(EXPECTED["pval"].drop("f3")))


  def test_feature_without_a_row_is_not_tested(fake_rpy2):
      fake_rpy2.output = CANNED[CANNED["feature"] != "f7"]
      out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
      assert out.loc["f7", ["effect", "pvalue"]].isna().all() and out["effect"].notna().sum() == 7


  def test_numeric_group_reports_its_slope(fake_rpy2):
      fake_rpy2.output = CANNED
      tdata = bt.datasets.toy()
      tdata.obs["ph"] = [5.1, 5.6, 6.0, 6.8, 7.1, 7.4]
      out = bt.da.maaslin3(tdata, "ph", seed=0)
      assert (out["contrast"] == "ph").all()
      assert fake_rpy2.calls[0][1][1]["x0"].dtype == np.float64


  def test_missing_rpy2_names_the_extra(monkeypatch):
      for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
          monkeypatch.setitem(sys.modules, name, None)
      with pytest.raises(ImportError, match=r"pip install 'biotapy\[r\]'"):
          bt.da.maaslin3(bt.datasets.toy(), "group")


  def test_missing_r_package_names_the_install_line(fake_rpy2):
      fake_rpy2.installed = {"ALDEx2"}
      with pytest.raises(ImportError, match=r'BiocManager::install\("maaslin3"\)'):
          bt.da.maaslin3(bt.datasets.toy(), "group")


  def test_shared_model_checks_apply(fake_rpy2):
      tdata = bt.datasets.toy()
      tdata.obs["arm"] = ["x"] * 3 + ["y"] * 3
      with pytest.raises(ValueError, match="are collinear"):
          bt.da.maaslin3(tdata, "group", covariates=["arm"])
      with pytest.raises(ValueError, match=r"obs\['age'\] is missing for 1 sample"):
          tdata.obs["age"] = [30.0, np.nan, 52, 38, 45, 60]
          bt.da.maaslin3(tdata, "group", covariates=["age"])
      assert fake_rpy2.calls == []


  def test_maaslin3_keeps_input(fake_rpy2, assert_unchanged):
      fake_rpy2.output = CANNED
      tdata = bt.datasets.toy()
      before = tdata.copy()
      bt.da.maaslin3(tdata, "group", seed=0)
      assert_unchanged(before, tdata)


  def test_seed_accepts_a_generator(fake_rpy2):
      fake_rpy2.output = CANNED
      bt.da.maaslin3(bt.datasets.toy(), "group", seed=np.random.default_rng(3))
      bt.da.maaslin3(bt.datasets.toy(), "group", seed=3)
      assert fake_rpy2.calls[0][1][3] == fake_rpy2.calls[1][1][3]
      generator = np.random.default_rng(3)
      bt.da.maaslin3(bt.datasets.toy(), "group", seed=generator)
      bt.da.maaslin3(bt.datasets.toy(), "group", seed=generator)
      assert fake_rpy2.calls[2][1][3] != fake_rpy2.calls[3][1][3]  # the caller's generator advances


  @pytest.mark.parametrize(("axis", "fix"), [("var", "var_names_make_unique"), ("obs", "obs_names_make_unique")])
  def test_repeated_names_raise(fake_rpy2, axis, fix):
      tdata = bt.datasets.toy()
      names = getattr(tdata, f"{axis}_names").tolist()
      setattr(tdata, f"{axis}_names", ["x", "x", *names[2:]])
      with pytest.raises(ValueError, match=rf"unique {axis} names.*\['x'\].*adata\.{fix}\(\)"):
          bt.da.maaslin3(tdata, "group")
      assert fake_rpy2.calls == []


  def test_seed_is_checked_before_rpy2_is_imported(monkeypatch):
      for name in ("rpy2", "rpy2.robjects", "rpy2.robjects.packages", "rpy2.robjects.pandas2ri"):
          monkeypatch.setitem(sys.modules, name, None)
      with pytest.raises(TypeError, match="seed must be"):
          bt.da.maaslin3(bt.datasets.toy(), "group", seed="abc")


  def test_each_r_warning_is_re_emitted_at_the_callers_line(fake_rpy2):
      fake_rpy2.output = CANNED
      fake_rpy2.warnings = ["values are unreliable"]
      with pytest.warns(UserWarning) as record:
          bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
      assert [str(w.message) for w in record] == ["da.maaslin3: R warned: values are unreliable"]
      assert record[0].filename == __file__


  def test_r_errors_name_the_function(fake_rpy2):
      fake_rpy2.output = sys.modules["rpy2.rinterface_lib.embedded"].RRuntimeError("Error in f() : boom\n")
      with pytest.raises(RuntimeError, match=r"da\.maaslin3: R stopped: Error in f\(\) : boom"):
          bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)


  def test_ordered_categoricals_reach_r_unordered(fake_rpy2):
      # R codes an ordered factor by polynomial contrasts, which would shrink a two-level effect by 1/sqrt(2).
      fake_rpy2.output = CANNED
      tdata = bt.datasets.toy()
      tdata.obs["group"] = pd.Categorical(tdata.obs["group"], categories=["A", "B"], ordered=True)
      tdata.obs["site"] = pd.Categorical(["x", "y", "z", "x", "y", "z"], categories=["x", "y", "z"], ordered=True)
      bt.da.maaslin3(tdata, "group", covariates=["site"], seed=0)
      metadata = fake_rpy2.calls[0][1][1]
      assert not metadata["x0"].cat.ordered and not metadata["x1"].cat.ordered


  @pytest.mark.r
  def test_maaslin3_on_toy_is_the_canned_table():
      out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
      # CANNED holds 6 significant digits of R's output, hence rtol=1e-5.
      for ours, theirs in {"effect": "coef", "se": "stderr", "pvalue": "pval"}.items():
          np.testing.assert_allclose(out[ours], EXPECTED[theirs], rtol=1e-5, err_msg=ours)


  @pytest.mark.r
  def test_seed_reproduces_the_table_in_r():
      tdata = bt.datasets.toy()
      first = bt.da.maaslin3(tdata, "group", seed=0)
      pd.testing.assert_frame_equal(bt.da.maaslin3(tdata, "group", seed=0), first)
      assert not bt.da.maaslin3(tdata, "group", seed=1)["pvalue"].equals(first["pvalue"])


  @pytest.mark.r
  def test_reference_sets_the_sign_in_r():
      tdata = bt.datasets.toy()
      a_first, b_first = bt.da.maaslin3(tdata, "group", seed=0), bt.da.maaslin3(tdata, "group", reference="B", seed=0)
      assert (b_first["contrast"] == "A vs B").all()
      np.testing.assert_allclose(b_first["effect"], -a_first["effect"], rtol=1e-9, atol=1e-12)
      np.testing.assert_allclose(b_first["se"], a_first["se"], rtol=1e-9)


  @pytest.mark.r
  def test_covariates_and_odd_column_names_reach_r():
      tdata = bt.datasets.toy()
      tdata.obs["body site"] = pd.Categorical(["x", "y", "x", "y", "x", "y"])
      tdata.obs["treated"] = [False] * 3 + [True] * 3
      adjusted = bt.da.maaslin3(tdata, "treated", covariates=["body site"], seed=0)
      assert (adjusted["contrast"] == "True vs False").all() and adjusted["effect"].notna().all()
      assert not np.allclose(adjusted["se"], bt.da.maaslin3(tdata, "treated", seed=0)["se"])


  @pytest.mark.r
  def test_bool_reference_sets_the_sign_in_r():
      tdata = bt.datasets.toy()
      tdata.obs["treated"] = [False] * 3 + [True] * 3
      default, swapped = (
          bt.da.maaslin3(tdata, "treated", seed=0),
          bt.da.maaslin3(tdata, "treated", reference="True", seed=0),
      )
      np.testing.assert_allclose(swapped["effect"], -default["effect"], rtol=1e-9, atol=1e-12)


  @pytest.mark.r
  def test_all_zero_feature_is_not_tested_in_r():
      tdata = bt.datasets.toy()
      dense = tdata.X.toarray()
      dense[:, tdata.var_names.get_loc("f7")] = 0
      tdata.X = sp.csr_matrix(dense)
      out = bt.da.maaslin3(tdata, "group", seed=0)
      assert out.loc["f7", ["effect", "pvalue", "qvalue"]].isna().all()
      assert np.isfinite(out.drop(index="f7")[["effect", "pvalue", "qvalue"]].to_numpy()).all()


  @pytest.mark.r
  def test_maaslin3_table_passes_the_consensus_checks():
      out = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)
      assert bt.da.consensus([out], min_methods=1)["n_tested"].tolist() == [1] * 8


  @pytest.mark.r
  def test_ordered_group_gives_the_unordered_table_in_r():
      tdata = bt.datasets.toy()
      plain = bt.da.maaslin3(tdata, "group", seed=0)
      tdata.obs["group"] = pd.Categorical(tdata.obs["group"], categories=["A", "B"], ordered=True)
      pd.testing.assert_frame_equal(bt.da.maaslin3(tdata, "group", seed=0), plain)
  ```
  `tests/da/test_maaslin3_golden.py`:
  ```python
  from pathlib import Path

  import numpy as np
  import pandas as pd
  import pytest
  from scipy.stats import false_discovery_control

  import biotapy as bt

  GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
  # Marker r only: it needs R, and the r-bridge CI job gives it the network job's pooch cache.
  pytestmark = pytest.mark.r


  @pytest.mark.parametrize(("formula", "covariates"), [("host", ()), ("host + log_depth", ("log_depth",))])
  def test_maaslin3_matches_maaslin3_maaslin3(benchmark, formula, covariates):
      # R: set.seed(1165433077); maaslin3(counts, meta, output, formula, min_prevalence = 0, evaluate_only = "abundance",
      # warn_prevalence = FALSE, subtract_median = TRUE), the host rows. The seed is the integer seed=20260927 becomes, so
      # the median test's simulations are the same and the numbers match to rounding.
      golden = pd.read_csv(GOLDEN / "maaslin3.csv.gz", dtype={"taxon_id": str}).query("formula == @formula")
      golden = golden.set_index("taxon_id").reindex(benchmark.var_names)
      out = bt.da.maaslin3(benchmark, "host", covariates=covariates, reference="other", seed=20260927)
      assert (out["contrast"] == "human vs other").all()
      # The 36 genera whose fit reports an error (no read in one host group): maaslin3 leaves them out of its q-values.
      untested = golden["error"]
      assert untested.sum() == 36 and out["effect"].isna().equals(untested.rename("effect"))
      ours, theirs = out[~untested], golden[~untested]
      for column, name in {"effect": "coef", "se": "stderr", "pvalue": "pval"}.items():
          np.testing.assert_allclose(ours[column], theirs[name], rtol=1e-7, err_msg=column)
      # BH over the group's p-values; maaslin3's qval pools them with the covariate's (20 calls against 45 with log_depth).
      np.testing.assert_allclose(ours["qvalue"], false_discovery_control(theirs["pval"]), rtol=1e-7)
  ```
- [x] **Step 7: Run, expect failure** -
  `uv run --group test pytest tests/da/test_maaslin3.py -q` -> `10 failed, 7
  deselected` (`AttributeError: module 'biotapy.da' has no attribute
  'maaslin3'`); in the R environment the seven `r` tests of that file (eight after fix round 1) fail the same way.
  Fix-round tests (a generator seed, repeated names, the seed checked before rpy2, R warnings
  and errors, as for `da.aldex2`) are in the same file and pass with the implementation.
- [x] **Step 8: Implement.** `src/biotapy/da/_maaslin3.py`:
  ```python
  """MaAsLin 3 through the R bridge: its abundance model, linear models of log2 relative abundance where present."""

  from collections.abc import Sequence

  import numpy as np
  import pandas as pd
  from anndata import AnnData

  from ._design import dense_counts, model
  from ._r import call_r, r_function, r_seed
  from ._schema import result

  # evaluate_only = "abundance" needs warn_prevalence = FALSE (maaslin3 stops otherwise). subtract_median = TRUE reports
  # each coefficient minus the median its p-value is tested against. The output folder maaslin3 insists on is deleted,
  # and pbapply's progress bar, which verbosity does not reach, is off for the call.
  _MAASLIN = """function(counts, metadata, formula, seed) {
    output <- tempfile("biotapy-maaslin3-")
    progress <- pbapply::pboptions(type = "none")
    on.exit({unlink(output, recursive = TRUE); pbapply::pboptions(progress)})
    set.seed(seed)
    fit <- maaslin3::maaslin3(
      counts, metadata, output, formula = formula, min_abundance = 0, min_prevalence = 0, evaluate_only = "abundance",
      warn_prevalence = FALSE, subtract_median = TRUE, plot_summary_plot = FALSE, plot_associations = FALSE,
      cores = 1, verbosity = "ERROR"
    )
    out <- fit$fit_data_abundance$results
    data.frame(feature = out$feature, metadata = out$metadata, coef = out$coef, stderr = out$stderr,
               pval = out$pval_individual, failed = !is.na(out$error))
  }"""


  def maaslin3(
      adata: AnnData,
      group: str,
      *,
      covariates: Sequence[str] = (),
      reference: str | None = None,
      seed: int | np.random.Generator | None = None,
  ) -> pd.DataFrame:
      """Differential abundance of each feature between two groups by MaAsLin 3's abundance model, run in R.

      Parameters
      ----------
      adata
          Samples x features; ``X`` holds raw counts.
      group
          The ``obs`` column whose effect is reported: a categorical, string or bool
          column with two levels, or a numeric column, whose effect is per standard
          deviation (MaAsLin 3 standardises numeric columns).
      covariates
          ``obs`` columns to adjust for: numeric ones standardised, others as one
          indicator per level against their first level.
      reference
          The level of a categorical ``group`` that the other level is compared
          with; by default its first category (sorted values for a string column).
      seed
          Seeds R's random number generator through ``set.seed``; MaAsLin 3's test
          against the median simulates, so the same seed gives the same table.

      Returns
      -------
      pandas.DataFrame
          One row per feature, in ``var_names`` order, indexed by ``feature``:
          ``effect`` (log2 fold change of the relative abundance where the feature
          is present, minus the median described in Notes), ``se``, ``pvalue``,
          ``qvalue`` (Benjamini-Hochberg over the tested features), ``direction``,
          ``method`` (``"maaslin3"``) and ``contrast``. A feature MaAsLin 3 could
          not fit has NaN values and ``direction`` 0.

      Raises
      ------
      ImportError
          rpy2 (the extra ``r``) or the R package maaslin3 is not installed.
      KeyError
          ``group`` or a covariate is not an ``obs`` column.
      TypeError
          ``covariates`` is not a list of column names, or ``reference`` is not a string.
      RuntimeError
          R stops with an error.
      ValueError
          ``X`` does not hold raw counts, has an empty sample, repeated
          ``var_names`` or ``obs_names``, or fewer than two
          features; a used ``obs`` column has missing values, is constant or is
          repeated; ``group`` has other than two levels; ``reference`` is not one of
          them or is given for a numeric ``group``; the model has at least as many
          terms as samples, or collinear columns.

      Notes
      -----
      R equivalent: ``maaslin3::maaslin3``
      Guide: :doc:`/guide/differential_abundance`

      Calls ``maaslin3::maaslin3(counts, metadata, output, formula, min_prevalence
      = 0, evaluate_only = "abundance", warn_prevalence = FALSE, subtract_median =
      TRUE)`` (Nickols et al.) through rpy2, with its other defaults: counts
      are divided by each sample's total (TSS), zeros are left out and the rest
      log2-transformed, one linear model per feature is fitted on the samples where
      it is present, and each coefficient is tested against the median coefficient
      (``median_comparison_abundance = TRUE``), MaAsLin 3's correction for
      compositionality, which draws 10,000 normal samples. That median is taken, as
      MaAsLin 3 computes it, over the features without a fit error whose own
      p-value is below 0.95. ``effect`` is that coefficient minus the median, so
      its sign is the side of the median the test is about; MaAsLin 3's default
      output reports the coefficient itself. Only the abundance model runs: the
      prevalence model's log-odds cannot share a column with fold changes. ``qvalue`` is the
      Benjamini-Hochberg correction of the group's p-values; MaAsLin 3's
      ``qval_individual`` corrects them together with every covariate's. A feature
      whose fit reports an error is not tested, as MaAsLin 3 leaves it out of its
      own correction. Swapping ``reference`` negates ``effect`` exactly, but the
      simulation draws around the coefficients, not their negatives, so the same
      ``seed`` moves ``pvalue`` by up to 1e-3 and ``qvalue`` by up to 2e-3 on the
      GlobalPatterns genera (the calls are the same there). Each R warning raised
      during the call is re-emitted as a ``UserWarning``, and an R error is raised
      as a ``RuntimeError``.

      Needs R with maaslin3 (``BiocManager::install("maaslin3")``) and ``pip
      install 'biotapy[r]'``, which builds rpy2 (GPL-2.0-or-later) against that R.
      ``seed`` becomes one integer for R's ``set.seed``, which sets the random
      state of the R session rpy2 embeds in Python. Plots are off and MaAsLin 3's
      output folder is a temporary one, deleted afterwards.

      rpy2 converts dense tables only, so ``X`` is densified once: 8 bytes x
      samples x features, plus R's copies of the table.

      References
      ----------
      Nickols WA, et al. MaAsLin 3: refining and extending generalized multivariable
      linear models for meta-omic association discovery. Nature Methods,
      doi:10.1038/s41592-025-02923-9.

      Examples
      --------
      >>> import biotapy as bt
      >>> table = bt.da.maaslin3(bt.datasets.toy(), "group", seed=0)  # doctest: +SKIP
      >>> table.loc["f6", "contrast"]  # doctest: +SKIP
      'B vs A'
      """
      frame, contrast = model(adata, group, covariates=covariates, reference=reference, func="da.maaslin3")
      seed_for_r = r_seed(seed)
      # rpy2 has no sparse converter (rules.md R6.2): the one dense copy of X.
      counts = pd.DataFrame(dense_counts(adata, func="da.maaslin3"), index=adata.obs_names, columns=adata.var_names)
      # Plain names in the formula, so an obs column such as "body site" needs no quoting. rpy2 makes an R factor, whose
      # first level maaslin3 keeps as reference, only from string categories; bool ones would arrive as sorted strings.
      metadata = pd.DataFrame(
          {f"x{i}": _string_levels(frame[name]) for i, name in enumerate(frame.columns)}, index=frame.index
      )
      fit = r_function(_MAASLIN, package="maaslin3", func="da.maaslin3")
      table = call_r(fit, counts, metadata, "~ " + " + ".join(metadata.columns), seed_for_r)
      rows = table[table["metadata"] == "x0"].set_index("feature").reindex(adata.var_names)
      untested = rows["failed"].fillna(True).to_numpy(bool) | ~np.isfinite(rows["pval"].to_numpy(np.float64))
      effect, se, pvalue = (np.where(untested, np.nan, rows[c].to_numpy(np.float64)) for c in ("coef", "stderr", "pval"))
      return result(adata.var_names, effect=effect, se=se, pvalue=pvalue, method="maaslin3", contrast=contrast)


  def _string_levels(values: pd.Series) -> pd.Series:
      """A categorical column with its levels as strings, in the same order; a numeric one as it is."""
      if isinstance(values.dtype, pd.CategoricalDtype):
          return values.cat.rename_categories([str(level) for level in values.cat.categories])
      return values
  ```
  `src/biotapy/da/__init__.py`:
  ```diff
  @@ -4,5 +4,6 @@ from ._aldex2 import aldex2
   from ._ancombc import ancombc2
   from ._consensus import consensus
   from ._linda import linda
  +from ._maaslin3 import maaslin3

  -__all__ = ["aldex2", "ancombc2", "consensus", "linda"]
  +__all__ = ["aldex2", "ancombc2", "consensus", "linda", "maaslin3"]
  ```
- [x] **Step 9: Run, expect pass** - `uv run --group test pytest tests/da -q`
  -> `145 passed, 22 deselected`; in the R environment `uv run --group test
  --extra r pytest -m r tests/da -q` -> `18 passed`. Also with `-W
  error::UserWarning`: rpy2's "categories are strings" warning must not appear.
- [x] **Step 10: Docs.**
  ````diff
  @@ -93,7 +93,7 @@ ALDEx2 and MaAsLin 3 exist only in R, so biotapy calls them there through
   rpy2 against that R (Linux and macOS; it fails to install when no R is found):

   ```bash
  -Rscript -e 'install.packages("BiocManager"); BiocManager::install("ALDEx2")'
  +Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("ALDEx2", "maaslin3"))'
   pip install 'biotapy[r]'
   ```

  @@ -118,6 +118,26 @@ Benjamini-Hochberg correction of ALDEx2's expected p-value `we.ep`, as for every
   `we.eBH` averages the corrections of the draws instead and calls more features (19 against 11 on
   those genera).

  +### MaAsLin 3
  +
  +`bt.da.maaslin3` runs MaAsLin 3's abundance model (`maaslin3::maaslin3` with `evaluate_only =
  +"abundance"`): counts become relative abundances, zeros are left out, and one linear model of the
  +log2 abundance per feature is fitted on the samples where the feature is present. Each coefficient
  +is tested against the median coefficient, MaAsLin 3's correction for
  +compositionality (the median over the features without a fit error whose own p-value is below
  +0.95, as MaAsLin 3 computes it); `effect` is the coefficient minus that median, so its sign says on which side
  +of the median the feature moved (MaAsLin 3 reports the coefficient itself unless asked to subtract):
  +
  +```python
  +table = bt.da.maaslin3(tdata, "group", covariates=["age"], seed=0)
  +```
  +
  +Numeric columns are standardised, so a numeric group's effect is per standard deviation, as in
  +`da.linda`. The test against the median simulates, so `seed` makes the p-values reproducible; on the
  +GlobalPatterns genera biotapy's numbers equal R's to 1e-12. The prevalence model, whose effects
  +are log-odds rather than fold changes, is not run. `qvalue` corrects the group's p-values only;
  +MaAsLin 3's `qval_individual` corrects them together with every covariate's, which with one numeric
  +covariate called 20 genera where biotapy calls 45. R warnings and errors surface as for ALDEx2.
  +
  +Swapping `reference` negates `effect` exactly, but the median test's simulation draws around the
  +coefficients rather than their negatives, so the same `seed` moves `pvalue` by up to 1e-3 and
  +`qvalue` by up to 2e-3 on those genera (the calls are the same); `seed` fixes the draws.
  +
   ## Where methods agree

   `bt.da.consensus` puts the tables of several methods side by side and counts, for each feature,
  ````
  ````diff
  @@ -104,6 +104,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
       da.ancombc2
       da.consensus
       da.linda
  +    da.maaslin3
   ```

   ## Plots
  ````
  ```diff
  @@ -69,7 +69,7 @@ BIOTAPY_DATA_DIR=.pooch uv run --group test pytest -m "network or golden"

   ### R bridge tests

  -Tests that call R through rpy2 (`bt.da.aldex2`) carry the marker `r` and are excluded from the
  +Tests that call R through rpy2 (`bt.da.aldex2`, `bt.da.maaslin3`) carry the marker `r` and are excluded from the
   runs above. They need R with the packages `tests/r/Dockerfile` installs and the `r` extra; those
   that read GlobalPatterns also need the pooch cache:

  ```
- [x] **Step 11: Contract.**
  ```diff
  @@ -53,6 +53,7 @@ sources:
      | DA methods | sign agreement and rank correlation of effect sizes; exact match only where the R method is deterministic | per method |
      | `da.linda` vs `MicrobiomeStat::linda(is.winsor = FALSE)` (deterministic) | `effect`, `se`, `pvalue`, `qvalue` elementwise, matched by taxon | `rtol=1e-7` |
      | `da.aldex2` vs `ALDEx2::aldex` (Monte Carlo; the golden's `set.seed` is the integer biotapy derives from `seed=20260927`, and `reference="human"` keeps R's level order, so the draws are the same) | `effect` vs `diff.btw` and `pvalue` vs `we.ep` elementwise, matched by taxon; `qvalue` vs BH of `we.ep` | `rtol=1e-7` (measured 4e-15 and 8e-13) |
  +   | `da.maaslin3` vs `maaslin3::maaslin3(evaluate_only = "abundance", subtract_median = TRUE)` (its median test simulates; the same seed as `da.aldex2`'s golden) | the same untested features (fit errors); `effect` vs `coef`, `se` vs `stderr`, `pvalue` vs `pval_individual` elementwise; `qvalue` vs BH of R's p-values over the group's rows | `rtol=1e-7` (measured 2e-14, 5e-15 and 7e-13), both models |
      | `da.ancombc2` vs `ANCOMBC::ancombc2` (deterministic, but its bias E-M can stop at 100 iterations before converging, on a slightly different iterate in scikit-bio) | the same untested features; `effect`, `se`, `pvalue` elementwise; Spearman correlation of effects; the same calls at `q < 0.05`, with R's p-values corrected over the tested features | `host` model: `effect` atol 0.015 (log2), `se` rtol 2e-3, `pvalue` atol 0.02, Spearman > 0.9999; `host + log_depth`: all three at 1e-6, Spearman > 0.999999 |

   5. Any looser tolerance is written in the test with a one-line comment giving the reason.
  ```
  The "Why" section's seed sentence now reads "(`da.aldex2` and `da.maaslin3` do: the goldens'
  `set.seed(1165433077)` is ... so a change in NumPy's stream fails those goldens loudly".
  Log line:
  ```text
  - **Update**: [r-golden-parity](contracts/r-golden-parity.md): `da.maaslin3` is compared with `maaslin3::maaslin3` (abundance model, median subtracted) elementwise, seed-matched. [phase-3-stats](roadmap/phase-3-stats.md) task 3.7 done.
  ```
- [x] **Step 12: Gate and commit** (tick 3.7)
  ```bash
  git add src/biotapy/da/_maaslin3.py src/biotapy/da/__init__.py tests/da/test_maaslin3.py \
    tests/da/test_maaslin3_golden.py docs/guide/differential_abundance.md docs/api.md docs/contributing.md \
    .knowledge/contracts/r-golden-parity.md .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  uvx prek run --all-files
  uv run --group doc sphinx-build -W -b html docs docs/_build/html
  git commit -m "feat(da): add the MaAsLin 3 bridge

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```

---
- [x] **Step 13: Fix round 1** (review of 4f2b39f; `fix(da)`). `_design._column` returns
  unordered categoricals: R codes an ordered factor by polynomial contrasts, which shrank
  MaAsLin 3's `effect` by 1/sqrt(2) for an ordered two-level group (measured ratio 0.7071 on
  4f2b39f; ANCOM-BC2 and LinDA, which code indicators, were not affected). Tests:
  `test_design.py::test_ordered_categorical_gives_the_unordered_table` (both native methods; they
  passed before), `test_maaslin3.py::test_ordered_categoricals_reach_r_unordered` (RED) and
  `test_ordered_group_gives_the_unordered_table_in_r` (RED). The docstring and guide state the
  median as maaslin3 computes it (error-free features with p < 0.95) and that swapping
  `reference` moves p by up to 1e-3 and q by up to 2e-3. Counts: 1221 passed, 47 deselected;
  `tests/da` 145 passed, 22 deselected; `-m r` 18 passed.

---

### Task 3.11: CI job `r-bridge`

**Files:** modify `.github/workflows/test.yaml`, `tests/test_ci.py`,
`tests/da/test_consensus.py`, `docs/contributing.md`,
`.knowledge/decisions/r-bridge-before-ports.md`.
**Not touched:** the other jobs (the network job's `-m "network or golden"`
selects no bridge test), `pyproject.toml`, `src/`.
**Interfaces:**
- Consumes: the four methods; the `benchmark` fixture; the network job's pooch
  cache key.
- Produces: the `r-bridge` job in `check.needs`; the exit gate's four-method
  consensus as a test.

- [x] **Step 1: Approval on record.** The job and the two actions pinned by SHA
  were approved (Phase 3 decision 15). New and to confirm (slice 3C decision 15):
  the job installs six apt `-dev` packages on the runner (CI only), sets
  `RPY2_CFFI_MODE=API`, runs on `ubuntu-24.04` instead of `ubuntu-latest`, and
  points CRAN at the image's dated P3M snapshot instead of the latest P3M.
- [x] **Step 2: Look up the SHAs** (never guessed):
  `gh api repos/r-lib/actions/git/ref/tags/v2 --jq .object.sha` ->
  `f9a764fea8d5c63df6ef9a5c7795bf7deb5d7e05`, which `gh api
  repos/r-lib/actions/tags` lists as `v2.14.0`. Both actions live in that
  repository, so one SHA pins both. If the tag has moved by execution time, use
  the new SHA and its version in the comment.
- [x] **Step 3: Failing tests.** `tests/test_ci.py`:
  ```diff
  @@ -35,6 +35,34 @@ def test_network_job_blocks_merges():
       assert "network" in WORKFLOW["jobs"]["check"]["needs"]


  +def test_r_bridge_job_runs_the_r_marker():
  +    steps = WORKFLOW["jobs"]["r-bridge"]["steps"]
  +    run = [step for step in steps if step.get("run", "").strip() == "uv run --group test --extra r pytest -m r"]
  +    assert run and "BIOTAPY_DATA_DIR" in run[0]["env"] and run[0]["env"]["RPY2_CFFI_MODE"] == "API"
  +
  +
  +def test_r_bridge_job_installs_the_image_r_and_packages():
  +    steps = {step.get("uses", "").split("@")[0]: step for step in WORKFLOW["jobs"]["r-bridge"]["steps"]}
  +    assert steps["r-lib/actions/setup-r"]["with"]["r-version"] == "4.5.3"
  +    dockerfile = (ROOT / "tests" / "r" / "Dockerfile").read_text(encoding="utf-8")
  +    assert "FROM rocker/r-ver:4.5.3" in dockerfile and 'version = "3.22"' in dockerfile
  +    packages = steps["r-lib/actions/setup-r-dependencies"]["with"]["packages"]
  +    assert {name.strip() for name in packages.split(",")} == {"bioc::ALDEx2", "bioc::maaslin3"}
  +
  +
  +def test_r_bridge_job_uses_the_image_cran_snapshot_and_has_a_timeout():
  +    job = WORKFLOW["jobs"]["r-bridge"]
  +    setup = next(step for step in job["steps"] if step.get("uses", "").startswith("r-lib/actions/setup-r@"))["with"]
  +    # The snapshot comes from the rocker/r-ver:4.5.3 base image (tests/r/Dockerfile names no URL): `R -e 'getOption("repos")'`.
  +    assert setup["cran"] == "https://p3m.dev/cran/__linux__/noble/2026-04-23"
  +    assert setup["use-public-rspm"] is False
  +    assert job["timeout-minutes"] == 30
  +
  +
  +def test_r_bridge_job_blocks_merges():
  +    assert "r-bridge" in WORKFLOW["jobs"]["check"]["needs"]
  +
  +
   def test_coverage_below_90_percent_fails_the_test_job():
       coverage = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["coverage"]
       assert coverage["report"]["fail_under"] == 90
  ```
  `tests/da/test_consensus.py` (an acceptance test: the code exists, so it
  passes where R is):
  ```diff
  @@ -166,3 +166,22 @@ def test_numpy_scalars_are_valid_options():
       results = [_table("a", [1.0], [0.01]), _table("b", [1.0], [0.01])]
       out = bt.da.consensus(results, alpha=np.float32(0.05), min_methods=np.int64(1))
       assert out["consensus"].tolist() == [True]
  +
  +
  +@pytest.mark.r
  +def test_four_methods_on_the_exit_gate_data(benchmark):
  +    # The Phase 3 exit-gate consensus: GlobalPatterns genera, human hosts against the rest, the two native methods and the
  +    # two R bridges. Seeded, so the counts are fixed: calls 208 (ANCOM-BC2), 118 (LinDA), 13 (ALDEx2), 52 (MaAsLin 3).
  +    results = [
  +        bt.da.ancombc2(benchmark, "host", reference="other"),
  +        bt.da.linda(benchmark, "host", reference="other"),
  +        bt.da.aldex2(benchmark, "host", reference="other", seed=0),
  +        bt.da.maaslin3(benchmark, "host", reference="other", seed=0),
  +    ]
  +    assert [int((table["qvalue"] < 0.05).sum()) for table in results] == [208, 118, 13, 52]
  +    table = bt.da.consensus(results)
  +    # ANCOM-BC2 and MaAsLin 3 cannot fit the 36 genera absent from one group; LinDA and ALDEx2 test all 636.
  +    assert table["n_tested"].value_counts().to_dict() == {4: 600, 2: 36}
  +    assert table["n_significant"].value_counts().sort_index().to_dict() == {0: 413, 1: 114, 2: 62, 3: 35, 4: 12}
  +    assert int(table["consensus"].sum()) == 109 and not table["conflict"].any()
  +    assert (table.loc[table["n_significant"] == 4, "direction"] == 1).all()
  ```
- [x] **Step 4: Run, expect failure** - `uv run --group test pytest
  tests/test_ci.py -q` -> `3 failed, 12 passed` (two `KeyError: 'r-bridge'`,
  and `'r-bridge' in ...check.needs` is false; fix round 1 adds a fourth, failing on `KeyError: 'timeout-minutes'`).
- [x] **Step 5: The job.** `.github/workflows/test.yaml`:
  ```diff
  @@ -161,6 +161,52 @@ jobs:
             BIOTAPY_DATA_DIR: ${{ github.workspace }}/.pooch
           run: uv run --group test pytest -m "network or golden"

  +  # Runs the R bridge tests (marker r, decisions/r-bridge-before-ports): R 4.5.3 with Bioconductor 3.22's ALDEx2 and
  +  # maaslin3, as tests/r/Dockerfile pins them, on CRAN's snapshot of that image, so the seed-matched goldens hold. The
  +  # four-method consensus reads GlobalPatterns through the network job's pooch cache.
  +  r-bridge:
  +    runs-on: ubuntu-24.04
  +    timeout-minutes: 30
  +    steps:
  +      - uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1 # v7.0.1
  +        with:
  +          filter: blob:none
  +          fetch-depth: 0
  +          persist-credentials: false
  +      - name: Install R
  +        uses: r-lib/actions/setup-r@f9a764fea8d5c63df6ef9a5c7795bf7deb5d7e05 # v2.14.0
  +        with:
  +          r-version: "4.5.3"
  +          use-public-rspm: false
  +          cran: https://p3m.dev/cran/__linux__/noble/2026-04-23
  +      - name: Install ALDEx2 and maaslin3
  +        uses: r-lib/actions/setup-r-dependencies@f9a764fea8d5c63df6ef9a5c7795bf7deb5d7e05 # v2.14.0
  +        with:
  +          packages: bioc::ALDEx2, bioc::maaslin3
  +          dependencies: '"hard"'
  +      # Not asserted: Bioconductor 3.22 fixes the release line and the seed-matched goldens catch result drift; the log shows
  +      # which versions ran.
  +      - name: Log the R and package versions
  +        run: Rscript -e 'cat(R.version.string, "\n"); for (pkg in c("ALDEx2", "maaslin3")) cat(pkg, as.character(packageVersion(pkg)), "\n")'
  +      # rpy2 builds against R's link flags; without these headers it silently falls back to a mode that cannot load R.
  +      - name: Install the libraries rpy2 links R with
  +        run: sudo apt-get update && sudo apt-get install -y --no-install-recommends libpcre2-dev libdeflate-dev libzstd-dev liblzma-dev libbz2-dev libicu-dev
  +      - name: Install uv
  +        uses: astral-sh/setup-uv@c18668ad3cf93ea998bef934396af7bb5c839dc7 # v10.2.0
  +        with:
  +          python-version: "3.13"
  +      - name: Cache pooch downloads
  +        uses: actions/cache@55cc8345863c7cc4c66a329aec7e433d2d1c52a9 # v6.1.0
  +        with:
  +          path: ${{ github.workspace }}/.pooch
  +          key: pooch-${{ hashFiles('src/biotapy/datasets/_remote.py') }}
  +      # Syncs the environment, so this is where rpy2 builds; the tests step reuses it.
  +      - name: Log the rpy2 version
  +        env:
  +          RPY2_CFFI_MODE: API
  +        run: uv run --group test --extra r python -c "import importlib.metadata as m; print('rpy2', m.version('rpy2'))"
  +      - name: Run the R bridge tests
  +        env:
  +          BIOTAPY_DATA_DIR: ${{ github.workspace }}/.pooch
  +          # Fail the rpy2 build instead of falling back to its ABI mode.
  +          RPY2_CFFI_MODE: API
  +        run: uv run --group test --extra r pytest -m r
  +
     # Builds the docs as Read the Docs does, executing every notebook (Phase 1 exit gate): the
     # phyloseq vignette downloads GlobalPatterns, enterotype and esophagus through the network
     # job's pooch cache.
  @@ -217,6 +258,7 @@ jobs:
         - lint
         - import-without-extras
         - network
  +      - r-bridge
         - docs
       runs-on: ubuntu-latest
       steps:
  ```
- [x] **Step 6: Run, expect pass** - the same command -> `16 passed`; `uvx prek
  run --all-files` passes (zizmor reads the new job); in the R environment
  `uv run --group test --extra r pytest -m r tests/da/test_consensus.py -q` ->
  `1 passed` (about 25 s).
- [x] **Step 7: Docs and decision.**
  ````diff
  @@ -77,6 +77,9 @@ that read GlobalPatterns also need the pooch cache:
   BIOTAPY_DATA_DIR=.pooch uv run --group test --extra r pytest -m r
   ```

  +CI runs them in the `r-bridge` job, with R 4.5.3 and the Bioconductor 3.22 packages the golden image
  +pins.
  +
   ### Regenerating the R golden files

   The golden CSVs under `tests/golden/` and the R-written fixtures under
  ````
  ```diff
  @@ -31,6 +31,8 @@ as open.[^spec]

   # Consequences
   - CI needs a job with R available for bridge tests; other jobs skip them via a
  -  pytest marker `r`.
  +  pytest marker `r`. Since slice 3C that is the `r-bridge` job in
  +  `.github/workflows/test.yaml` (R 4.5.3, Bioconductor 3.22, the image's CRAN
  +  snapshot), which blocks merges.

   [^spec]: Python Microbiome Toolkit development report, sections Risks and Open questions
  ```
  Log line:
  ```text
  - **Update**: [r-bridge-before-ports](decisions/r-bridge-before-ports.md) consequence: the `r-bridge` CI job runs the `r` tests. [phase-3-stats](roadmap/phase-3-stats.md) task 3.11 done.
  ```
- [x] **Step 8: Gate and commit** (tick 3.11)
  ```bash
  git add .github/workflows/test.yaml tests/test_ci.py tests/da/test_consensus.py docs/contributing.md \
    .knowledge/decisions/r-bridge-before-ports.md .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
  uvx prek run --all-files
  git commit -m "ci: run the R bridge tests in an r-bridge job

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  The job runs for the first time on the pull request (Checkpoint C); read its
  log for the rpy2 build line (`CFFI_MODE.API`), the package install time and
  the `r` count (`19 passed`).

### Checkpoint C - review slice 3C
- [x] Review the whole slice (superpowers:requesting-code-review, opus: the
  bridges are glue whose failure modes are silent sign and unit errors) against
  every contract, pure-by-default, the Phase 3 and slice 3C review focus and the
  measured parity; then a fix pass, one commit per finding, each with a test (an
  `r` test where only R shows the bug, plus a no-R test through `fake_rpy2`
  where the mapping can show it). Record the counts (Critical / Important /
  Minor) and the fix range here. Done: 0 Critical, 4 Important, 5 Minor and 2
  nits; the fix range is the commits after `39383e7` (M3 and M4 belong to the
  knowledge step and the PR).
  Recorded: 0 Critical / 4 Important / 5 Minor; fix pass `39383e7..64fe39d`;
  re-review 9/9.
- [x] Run the exit-gate check for 3C: in the R environment `uv run --group test
  --extra r pytest -m r -q` (expected `19 passed`; `24 passed` after the
  Checkpoint C fix pass), and on the host `uv run
  --group test pytest -q -W error::UserWarning` (expected `1224 passed, 48 deselected`;
  the gate gave `1225` after fix round 1's timeout test and `1230 passed, 53
  deselected` after the Checkpoint C fix pass) and
  `uv run --group test pytest -m "golden or network" -q` (expected `36
  passed`); record the three counts. The Phase 3 exit gate's "`da.aldex2`,
  `da.maaslin3` in the `r-bridge` job" box is ticked when that job is green on
  the pull request, and that job's log must show `CFFI_MODE.API` and `24 passed`
  (the apt packages, the pak install and rpy2's build cannot be checked locally).
  Recorded at `64fe39d`: host `1230 passed, 53 deselected`; host golden/network
  `36 passed, 1247 deselected`; container `-m r` `24 passed, 1259 deselected`.
- [x] Knowledge (codebase-map templates; R12.2-R12.4):
  - **Update `.knowledge/modules/da.md`** for the bridges: Responsibility gains
    `aldex2` and `maaslin3` (R-only methods through rpy2, extra `r`) and drops
    "the R bridges are later tasks"; Entry points `_aldex2.py:aldex2`,
    `_conditions`, `_maaslin3.py:maaslin3`, `_string_levels`,
    `_r.py:r_seed`, `r_function`, `call_r`; Invariants: rpy2 only through
    `import_optional` inside `_r.py`, one `set.seed(r_seed(seed))` per call
    inside the R closure, local converter only, BH recomputed (not `we.eBH`, not
    `qval_individual`), ALDEx2 labels `"0"`/`"1"`, MaAsLin 3 factors with string
    levels and plain names `x0..`, fit errors untested; Gotchas: ALDEx2 breaks on
    a factor `conds` and sorts labels by locale, `aldex.effect` needs two samples
    per level, maaslin3 stops on `evaluate_only` with `warn_prevalence = TRUE`,
    its raw coefficient is not what its p-value tests (hence `subtract_median`),
    `qval_individual` pools covariates, rpy2 sends non-string categories as
    sorted strings with only a UserWarning, a namespace loaded by `::` during an
    rpy2 call prints "stack imbalance", rpy2's silent ABI fallback without the
    link headers, R `NA` strings return as `NA_character_`; Dependencies: rpy2
    (extra `r`), R packages ALDEx2 and maaslin3 (user-installed);
    Verification adds `uv run --group test --extra r pytest -m r tests/da -q`.
    Description copied into `modules/index.md`.
  - **Update `.knowledge/decisions/optional-heavy-dependencies.md`** (R9.2), a
    sentence in the core-dependencies history: "Phase 3 task 3.6 added the extra
    `r` (`rpy2>=3.6.8`, approved 2026-10-05): ALDEx2 and MaAsLin 3 exist only in
    R. rpy2 (GPL-2.0-or-later) brings rpy2-rinterface, rpy2-robjects, cffi,
    jinja2 and tzlocal; rpy2-rinterface is an sdist on Linux that links against
    the user's R and needs R's link headers (`libpcre2-dev`, `libdeflate-dev`,
    `libzstd-dev` on Ubuntu), and fails to install without R. biotapy imports it
    only through `import_optional` in `da/_r.py` and never bundles it." The
    extras table's `r` row gains "(rpy2 3.6.8)".
  - **Refresh frontmatter** (`generated.at`, `commit`) of the concepts slice 3C
    edited: r-golden-parity, r-bridge-before-ports (its `verified` stays the
    human's), da, optional-heavy-dependencies, phase-3-stats.
  - **Update `.knowledge/roadmap/phase-3-stats.md`**: tick Checkpoint C's
    boxes; task 3.14's line drops "for the R bridges" (done here).
  - **Verification bump** for the concepts `bash scripts/knowledge_stale.sh
    --touched` lists whose statements still hold, as Checkpoints A and B did.
  - Log lines for each.
- [x] Push `phase-3c`, open the PR and merge (merge commit) when CI is green,
  including the new `r-bridge` job, docs, the network job and the knowledge
  report (push, PR and merge on green approved for Phase 3 slice branches,
  2026-10-05). If `r-bridge` fails only on numbers (CI's R build differs from
  the image), stop and report the measured differences instead of loosening a
  tolerance (R11.3, R11.5). Done 2026-10-06: PR #23 merged as `29efd9c`, Test run
  `37479113322` 21/21 green (`r-bridge` `24 passed`, `network` `36 passed`).
- [x] Ask the user to review slice 3C before slice 3D is expanded (R1.2a). Approved 2026-10-07.

### Slice 3C decisions for the user
Recommended answer first. Items 1-5 change what the approved plan said; 6-14
resolve the outline or an [UNVERIFIED]; 15-16 need fresh approval.

1. **Bridge tests carry the marker `r` only**, never `golden` or `network`, even
   when they read a golden file or download GlobalPatterns. Then `pytest -m
   "golden or network"` (contributing.md, the network job, this plan's gate)
   never needs R, and the `r-bridge` job gets the pooch cache. The outline gave
   the consensus test `r` + `network`, which the network job would select and
   fail. Alternative: keep the extra markers and change the network job to `-m
   "(network or golden) and not r"` (and every documented command).
2. **Bridge goldens are seed-matched and compared elementwise at `rtol=1e-7`**,
   not by signs and ranks. The R script calls `set.seed(1165433077)`, the
   integer `seed=20260927` becomes, so R and biotapy draw the same numbers;
   measured differences are 1e-14 to 1e-12. The outline said sign agreement and
   Spearman "numbers fixed when measured". Cost: the R script hard-codes a
   NumPy-derived integer, so a change in NumPy's `Generator.integers` stream
   would break the two goldens loudly (not silently).
3. **ALDEx2's `qvalue` is BH of `we.ep`**, computed by `_schema.result` as for
   every method, not ALDEx2's `we.eBH` (design note 7 said `qvalue = we.eBH`).
   `we.eBH` averages the BH values of the Monte Carlo draws, a different
   quantity; on the benchmark it calls 19 genera where BH of `we.ep` calls 11.
   Alternative: carry `we.eBH` (needs a `qvalue` argument in `result()` and
   breaks "BH over finite p" for one method).
4. **MaAsLin 3's `qvalue` is BH over the group's p-values only**, not
   `qval_individual`, which pools the covariates' p-values (and, without
   `evaluate_only`, the prevalence model's). With one covariate it calls 45
   genera against MaAsLin 3's 20; without covariates the two agree (52).
   Design note 5 said `qval_individual` is BH. Alternative: carry
   `qval_individual` (a covariate then changes the group's q-values).
5. **MaAsLin 3's `effect` is the coefficient minus the median coefficient**
   (`subtract_median = TRUE`), the quantity its default p-value tests, not the
   raw `coef` (design note 7). Same p-values; the sign differs from the raw
   coefficient for 109 of 600 tested genera (none of the 52 called), and agrees
   with LinDA's on 484 of 600 instead of 437. Alternative: the raw coefficient,
   whose sign can contradict the test.
6. **MaAsLin 3 results come from the return value** (`fit_data_abundance$results`),
   not `all_results.tsv`; MaAsLin 3 still writes its folder, an R `tempfile()`
   deleted on exit; pbapply's progress bar is off for the call.
7. **MaAsLin 3 gets factors, not a `reference` string**, and plain names
   `x0, x1, ...` with `"~ x0 + x1"` (as `da.ancombc2` does for patsy): its
   `"var,level"` syntax breaks on commas in level names and its formula on
   column names such as `"body site"`. Categorical levels go as strings: bool
   categories otherwise reach R as plain strings that maaslin3 sorts, which
   flipped `reference="True"` silently in the prototype.
8. **ALDEx2 receives the labels `"0"` (reference) and `"1"`**: `aldex.clr`
   fails on a factor, and the order of character labels depends on R's locale.
9. **ALDEx2 refuses a numeric group, a level in fewer than two samples, and
   `mc_samples` that is not an int >= 1**, with messages naming the argument;
   ALDEx2 itself would stop on the second deep inside R.
10. **Untested rows:** a feature ALDEx2 drops (no read) and a MaAsLin 3 row with
    an `error` are NaN and left out of BH (MaAsLin 3 also leaves them out of its
    own q-values; 36 genera on the benchmark, the same as ANCOM-BC2's).
11. **No mypy override for rpy2.** rpy2 3.6.8 ships `py.typed` only in
    `rinterface`, `rinterface_lib`, `rlike` and `situation`, not `robjects`;
    biotapy reaches rpy2 only through `import_optional` (typed `Any`), and the
    lint environment cannot install rpy2 without R. The outline planned a
    `follow_untyped_imports` override.
12. **`da/_r.py` is three functions** (`r_seed`, `r_function(code, *, package,
    func)`, `call_r(function, *args)`), not the outline's `r_packages(*names)`
    and `to_r_frame`/`from_r_frame`: each bridge needs one package and one call,
    and the local converter converts both ways. Packages load with `importr`
    because a namespace loaded by `::` during an rpy2 call prints R "stack
    imbalance" warnings.
13. **No-R tests replace the rpy2 modules** (`fake_rpy2` in `sys.modules`)
    rather than monkeypatching a private R call (design note 8), so every line
    of biotapy's bridge code runs without R; the canned tables are R's real
    output on `toy()`, checked by `r` tests.
14. **Examples are `# doctest: +SKIP`**, as the network examples in
    `datasets` are: the doctest run has no R. The `r` tests run the same calls.
15. **CI job details (fresh approval):** `runs-on: ubuntu-24.04`; CRAN at the
    image's snapshot `https://p3m.dev/cran/__linux__/noble/2026-04-23` with
    `use-public-rspm: false`, so CI's CRAN dependencies equal the image's and the
    seed-matched goldens hold (alternative: latest P3M, faster updates, digits
    may drift); `setup-r-dependencies` with `packages: bioc::ALDEx2,
    bioc::maaslin3` and `dependencies: '"hard"'` (its default `deps::.` needs a
    DESCRIPTION file); apt `libpcre2-dev libdeflate-dev libzstd-dev liblzma-dev
    libbz2-dev libicu-dev` on the runner; `RPY2_CFFI_MODE=API` so rpy2's silent
    ABI fallback fails the build instead of the import. Both actions at
    `f9a764fea8d5c63df6ef9a5c7795bf7deb5d7e05` (v2.14.0).
16. **`VERSIONS.txt` also records multcomp** (1.4.30), whose `glht` gives
    MaAsLin 3's median-test p-values. A record, not a pin.

Edits outside this section that the plan commit should make: design note 5's
"ALDEx2's `we.eBH` and MaAsLin 3's `qval_individual` (abundance model only) are
BH" and design note 7's "`qvalue = we.eBH`", "`effect = coef`", "`qvalue =
qval_individual`" and "read back from `all_results.tsv`" follow decisions 3-6 if
approved; the global constraint on mypy and rpy2 follows decision 11; design
note 8's "monkeypatching the private R call" follows decision 13; the outline
"Slice 3C - R bridges (outline)" is replaced by this section.

### Slice 3C self-review
1. **Spec coverage.** The extra `r` (3.6 Step 8) and its R9.2 line
   (Checkpoint C); `da/_r.py` with `import_optional`, `importr`, the install
   message, `set.seed` from `as_generator`, a local converter (3.6); ALDEx2's
   arguments, column mapping, two-group rule and seeding (3.6, decisions 3, 8,
   9); MaAsLin 3 abundance-only with `min_prevalence = 0`, the reference, the
   output folder, the mapping (3.7, decisions 4-7); every [UNVERIFIED] of the
   outline resolved: `evaluate_only` (needs `warn_prevalence = FALSE`), the
   `reference` syntax (not used), `name`/`metadata` values (`metadata == "x0"`),
   log2 (`LOG` is `log2`), rpy2's failure without R (`ImportError`), rpy2's
   `py.typed` (partial), pak's Bioconductor for R 4.5.3 (3.22), the job's time
   (measured locally); tests without R (canned output, seed derivation, missing
   rpy2, missing R package) and with R (`r`, registered already); goldens on the
   benchmark, seed-matched; the CI job pinned by SHA, in `check.needs`, with a
   `tests/test_ci.py` assertion; the four-method consensus with its counts;
   Checkpoint C with knowledge and push/PR/merge.
2. **Placeholder scan.** Every step carries the full file or the exact diff,
   rendered from the scratch clone's commits that passed every gate; no "TBD".
   The values left to the executor are the log heading's date and, if the `v2`
   tag moves, the action SHA (Step 2 says how to look it up).
3. **Type consistency.** `r_seed(seed) -> int`, `r_function(code, *, package,
   func) -> Callable[..., Any]`, `call_r(function, /, *args) -> pd.DataFrame`;
   `aldex2(adata, group, *, mc_samples=128, reference=None, seed=None)` and
   `maaslin3(adata, group, *, covariates=(), reference=None, seed=None)` as
   roadmap 3.6 and 3.7 give them; both return `_schema.result`'s table; the
   fixture's `output`, `installed`, `calls` as the tests use them.
4. **Review focus.** Every test named in the slice 3C review focus exists in
   the rendered code (checked by searching this section for each name).
5. **Known residual risks and what was not verified.**
   - CI was not run: the `r-bridge` job's first run is on the pull request.
     Unverified there: `setup-r`'s `cran` input with `use-public-rspm: false`
     giving P3M binaries on the runner, the job time (about 5 minutes cold was
     measured in a container, 2 min 40 s of it for 98 R packages), and whether
     GitHub's R 4.5.3 build gives the same digits as the rocker image (the
     goldens' `rtol=1e-7` leaves room for BLAS-level differences; the Monte Carlo
     draws use R's own RNG, which is platform independent).
   - The seed-matched goldens hard-code NumPy's `default_rng(20260927)
     .integers(2**31 - 1)`; NumPy does not promise that stream across releases.
     A change fails both goldens loudly; the fix is to re-derive the integer and
     regenerate.
   - [Superseded by Checkpoint C (I2): the bridge wrapper restores R's
     `.Random.seed`, or its absence, on exit, so `set.seed` leaves the user's R
     session as it found it (R3.4); the docstrings no longer disclose a changed
     state.] ~~`set.seed` changes the R session's global random state; the
     docstrings say so. Restoring the previous state would need more R glue (R2.3).~~
   - [Superseded by fix round 2 (eb5040b): R warnings are caught in R and
     re-emitted as `UserWarning`; the messages themselves are untouched. Since
     Checkpoint C (M1) an R error carries the warnings raised before it in its
     message.] ~~R's console messages reach Python as log records of
     `rpy2.rinterface_lib.callbacks`, not as `warnings`; biotapy does not
     translate them.~~
   - The bridges are Linux/macOS-first (rpy2 builds from source); Windows and
     macOS were not run.
   - A missing R after rpy2 was built surfaces as the "pip install
     'biotapy[r]'" message with the real cause chained (`libR.so` not found):
     accurate about the extra, indirect about R.
   - MaAsLin 3's median test makes its p-values move with the seed by about
     0.001; its `effect` and `se` do not.
   - The ALDEx2 image layer was rebuilt from the committed Dockerfile in the
     prototype and the ALDEx2 and ANCOMBC layers recompiled (a pruned build cache); that image regenerated all 41 files at `b0519fe` byte for byte, so the committed Dockerfile, not a cached layer, reproduces the goldens. The maaslin3 image was confirmed only by a cache hit on the committed file (same image id).

## Slice 3D - Docs and release

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax for
> tracking. This section replaces the slice 3D outline (rules.md R1.2a).

**Goal:** a user can read, per differential abundance method, what it fits, what biotapy changes
from the R defaults and how closely it matches R; runs the exit-gate analysis (GlobalPatterns
genera, human-associated samples against the rest) in an executed tutorial; finds the Phase 3
functions on the Coming-from-R page; sees asv baselines for PhILR, LinDA and ANCOM-BC2; and gets
biotapy 0.3.0 on PyPI, with every outward step (tag, GitHub release, PyPI upload) waiting for the
user's explicit approval.

**Architecture:** docs, tests, benchmarks, knowledge and the release only. No `src/` file
changes, except the two approved `docs(da)` docstring commits (`3b6239f`, `b9cb006`) and the Checkpoint D fix pass's docstring lines. New pages: `docs/methods/` (an index and five pages) and
`docs/tutorials/differential_abundance.md`; the DA guide is restructured (the per-method depth
moves to the method pages). The tutorial executes ANCOM-BC2 and LinDA; the two R bridges are shown
as a non-executed block with the numbers the `r-bridge` CI job measures, kept in step by a test
that shares constants with the `r` test. New asv module `benchmarks/benchmarks/da.py` and a
`Philr` class in `pp.py`.

**Tech stack:** sphinx 9.1.0 · myst-nb 1.4.0 (`nb_execution_mode = "cache"`,
`nb_execution_timeout = 300`, both already set) · MyST `dollarmath` (already enabled) with
sphinxcontrib-katex · asv 0.6.6 · hatchling/`uv build` · twine. No new dependency (R9.1).

**Spec:** this concept's design note 10 (exit-gate notebook), the slice 3D outline it replaces,
the data-model-slots contract's "DA results" (design notes 5, 7 and 8 are superseded for the
bridges), rules.md, [cut-a-release](/playbooks/cut-a-release.md) and Phase 2's slice 2D (how 0.2.0
was documented, benchmarked and released).

**How slice 3D was checked.** Every file below was written into a scratch clone of the
repository at `29efd9c` (master after PR #23, slice 3C) and replayed as one commit per step that
says "commit": six on `phase-3d` (3.10, 3.10b, 3.12, 3.13 suite, 3.13 page, 3.14) and two on
`release-0.3.0` (the release and its knowledge commit). Each committed state was gated with a
clean tree (`git status --short` empty), `BIOTAPY_DATA_DIR` set to an absolute scratch pooch
path, `docs/_build` and `docs/generated` deleted before every docs build:

| Commit (scratch hash) | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` (tutorial execution) | `asv check --python=same` |
|---|---|---|---|---|---|
| `29efd9c` (baseline) | passed | 1230 passed, 53 deselected | 36 passed, 1247 deselected | build succeeded, 52.9 s | No problems found |
| 3.10 `b0a7975` | passed | 1230 passed, 53 deselected | 36 passed, 1247 deselected | build succeeded, 39.8 s | No problems found |
| 3.10b `508de9b` | passed | 1231 passed, 53 deselected | 36 passed, 1248 deselected | build succeeded, 50.3 s (8.9 s) | No problems found |
| 3.12 `3588378` | passed | 1238 passed, 53 deselected | 36 passed, 1255 deselected | build succeeded, 50.1 s (8.2 s) | No problems found |
| 3.13 suite `7b15403` | passed | 1238 passed, 53 deselected | 36 passed, 1255 deselected | build succeeded, 49.8 s | No problems found |
| 3.13 page `d6d82ef` | passed | 1238 passed, 53 deselected | 36 passed, 1255 deselected | build succeeded, 60.2 s | No problems found |
| 3.14 `048fb2a` | passed | 1238 passed, 53 deselected | 36 passed, 1255 deselected | build succeeded, 47.2 s | No problems found |
| 3.15 release `93f47c3` | passed | 1238 passed, 53 deselected | 36 passed, 1255 deselected | build succeeded, 49.5 s | No problems found |
| 3.15 knowledge `cc371a6` | passed | 1238 passed, 53 deselected | 36 passed, 1255 deselected | build succeeded, 68.3 s | No problems found |

(The hashes are the scratch clone's. The tutorial executed in 8.0-9.1 s in every build; the
whole build stayed under 70 s with the pooch cache warm.)

- prek covers ruff 0.16.9 check and format, `mypy --strict`, import-linter, pyproject-fmt,
  zizmor and the file hooks (14 hooks, all `Passed`). `~/.cache/biotapy` was checked absent after
  every gate. A relative `BIOTAPY_DATA_DIR` resolves inside `docs/` when myst-nb runs a notebook
  (the prototype's first build wrote a `docs/pooch/` that way): always pass an absolute path.
- The RED results in each "expect failure" step were reproduced by running the new tests
  against the previous commit's code.
- **R side.** The `-m r` suite ran in the local `biotapy-r-py` image (R 4.5.3, ALDEx2 1.42.0,
  maaslin3 1.2.0, rpy2 3.6.8, Python 3.13; never pushed) on an export of `508de9b`, the commit that
  touches the `r` test: `24 passed, 1260 deselected`. The four-method consensus, run there on the
  exit-gate data with the tutorial's code, took 16.9 s and gave calls 208 (ANCOM-BC2), 118 (LinDA),
  13 (ALDEx2) and 52 (MaAsLin 3); `n_tested` 4 for 600 genera and 2 for 36; `n_significant`
  0/1/2/3/4 for 413/114/62/35/12 genera; consensus 223, 109, 47 and 12 genera at `min_methods`
  1, 2, 3 and 4; no conflict; of the 109, 106 are up in human-associated samples and 3 down; all 12
  four-way calls are up. These equal slice 3C's measurement.
- **CI evidence already on record (PR #23, merged as `29efd9c` on 2026-10-06).** Test workflow
  run `37479113322`, every one of its 21 checks `SUCCESS`: the `network` job (job
  `112322946830`) `36 passed, 1247 deselected`, which holds the `pp.clr`, `pp.philr`, `da.linda`
  and `da.ancombc2` goldens; the `r-bridge` job (job `112322946896`) with R 4.5.3, ALDEx2 1.42.0,
  maaslin3 1.2.0, rpy2 3.6.8 and `RPY2_CFFI_MODE: API`, `24 passed, 1259 deselected`, which holds
  the `da.aldex2` and `da.maaslin3` goldens; the `docs` job (job `112322946845`) executed the four
  existing notebooks. (The `r-bridge` log shows the env var, not a `CFFI_MODE.API` build line: uv
  builds rpy2 quietly. The import and the 24 tests are the proof the API build worked.)
- **The release state.** `uv build` and `uvx twine check --strict dist/*` PASSED for both
  files; the sdist holds no `benchmarks/` and does hold `docs/methods/` and the tutorial; the
  wheel's METADATA reads `Version: 0.3.0`, `Requires-Dist: mudata>=0.4`, `Provides-Extra: r` and
  `Requires-Dist: rpy2>=3.6.8; extra == 'r'`; `pytest` inside the extracted sdist: `1145 passed,
  53 deselected`; `import biotapy; biotapy.__version__` prints `0.3.0`.
- **Benchmarks.** asv 0.6.6 ran on `7b15403` (load average 2.21, 3.04, 2.75 before and 5.37,
  4.03, 3.12 after; 61 s): `time_philr` 3.04 s, `peakmem_philr` 2.14 GB, `time_linda` 469 ms,
  `peakmem_linda` 1.15 GB, `time_ancombc2` 3.92 s, `peakmem_ancombc2` 1.07 GB. asv's
  `branches: ["master"]` needs a local `master` branch: the scratch clone had none and asv stopped
  with `Unknown branch master in configuration` until `git branch master origin/master` (the real
  repository has one). Outside asv, the calls' own peaks above the loaded table: `pp.philr` 1.72 GB
  (4.3x its 400 MB dense table), `da.linda` 821 MB (5.1x its 160 MB table), `da.ancombc2` 723 MB
  (4.5x).
- **APIs and facts checked** (R2.2): myst-nb 1.4.0's `NbParserConfig` (printed by every build:
  `execution_mode='cache'`, `execution_timeout=300`, `execution_raise_on_error=True`); a fenced
  ` ```python ` block in a MyST notebook is not executed, a `{code-cell}` is;
  `tests/test_ci.py::test_no_page_overrides_the_notebook_execution_settings` forbids per-page
  `mystnb` metadata, so a skip-execution route is closed. Read the Docs runs `uvx hatch run
  docs:build` (`.readthedocs.yaml`), i.e. `sphinx-build -M html docs docs/_build -W`, executing
  every notebook itself with no pooch cache; the phyloseq vignette already downloads GlobalPatterns
  in the same build, so the tutorial adds no download. The R defaults in the method pages were read
  with `formals()` in the golden image: MicrobiomeStat 1.4 `linda`, ANCOMBC 2.12.0 `ancombc2`,
  ALDEx2 1.42.0 `aldex`, maaslin3 1.2.0 `maaslin3`. scikit-bio 0.7.4's ANCOM-BC2 steps (E-M
  mixture bias, sampling-fraction correction and refit, `(sqrt(v) + sqrt(v_delta))^2` plus the
  `var_quantile` offset, two-sided t with per-feature df) were read at source. Every citation was
  checked against Crossref (LinDA, ANCOM-BC2, ALDEx2, MaAsLin 3 at *Nature Methods* 23:554-564,
  Nearing, Pelto). The phyloseq analysis vignette's "Multiple Testing and Differential Abundance"
  section only points to DESeq2 (read at source on GitHub).

### Slice 3D design

Each settled question gives the answer and the reason. **(user)** marks the ones repeated under
"Slice 3D decisions for the user".

#### Where the work goes

| File | Holds |
|---|---|
| `docs/methods/{index,linda,ancombc2,aldex2,maaslin3,consensus}.md` (new, 3.10) | one page per method: model, units, the R defaults biotapy changes, agreement with R, reference |
| `docs/guide/differential_abundance.md` (3.10, 3.10b) | the workflow: one table, choosing methods (a table linking the pages), methods that run in R (label `da-methods-in-r`), consensus, plot |
| `docs/index.md` (3.10) | `methods/index.md` in the "User guide" toctree |
| `docs/tutorials/differential_abundance.md` (new, 3.10b), `docs/tutorials/index.md` | the exit-gate notebook |
| `tests/da/test_consensus.py` (3.10b) | the `FOUR_METHOD_*` constants shared by the `r` test and a new tutorial test |
| `tests/test_coming_from_r.py`, `docs/_data/r_idioms.toml`, `docs/coming_from_r.md`, `docs/tutorials/phyloseq_analysis.md` (3.12) | the Phase 3 rows pinned; "not in 0.3" |
| `benchmarks/benchmarks/{_data,pp,da}.py` (3.13), `docs/performance.md` | the suite and its baselines |
| `.knowledge/...` (3.14, 3.15) | `modules/da.md`, `playbooks/add-a-function.md`, `playbooks/cut-a-release.md`, verification bumps |
| `pyproject.toml`, `CHANGELOG.md`, `README.md` (3.15) | the release |

#### 1. Method pages and the guide (3.10)

- **Five pages and an index under `docs/methods/`, linked from the "User guide" toctree. (user)**
  Each page: a runnable call, `## Model` (MyST `dollarmath`, as the function guide's functional
  redundancy section), `## Units`, `## Compared with R` (a table of the R function's arguments,
  their R defaults as `formals()` gives them, and what biotapy does), `## Agreement with R` (the
  golden test's measured numbers), `## Choosing the reference` where the method is not
  antisymmetric (ANCOM-BC2, ALDEx2, MaAsLin 3), and `## Reference` with the DOI.
- **The guide keeps the workflow and drops the per-method depth.** Its "LinDA", "ANCOM-BC2",
  "ALDEx2" and "MaAsLin 3" sections move, text kept, into the pages; a new "Choosing methods"
  section holds a comparison table that links each page and the two cross-method warnings
  (untested features, `reference`). "Methods that run in R" keeps install, licence, the
  `ImportError`, densifying, and now also the generic bridge behaviour (seed and R's restored
  state, warnings, errors, repeated names) that used to sit in the ALDEx2 subsection. The guide
  stays the docstrings' `Guide:` target, so no `src/` file changes (bar the docstring-prose commits the ledger approved later). Duplicating the depth in both
  places would let them drift.
- **Citations are plain text with a DOI link**, as the function guide cites Tian et al.; no
  `references.bib` entries (the bibliography page would need `{cite}` roles across pages and a new
  convention). The consensus page cites Nearing et al. and Pelto et al.
- No new test: the docs build with `-W` (nitpicky, so every `{func}`, `{doc}` and `{ref}` must
  resolve) is the check, as for Phase 2's guide pages.

#### 2. The exit-gate tutorial (3.10b)

- **Data: the `benchmark` fixture's contrast exactly.** GlobalPatterns glommed to genus, one
  `pp.filter_features(min_prevalence=0.2)`, `host` = `"human"` for Feces, Skin and Tongue (9
  samples) and `"other"` for the 17 others, `reference="other"`. The 17 "others" are 14
  environmental samples and 3 mock communities; the tutorial's note says so, and that the contrast
  is a demonstration. Design note 10's "8 vs 18" was a miscount. Keeping the fixture's contrast
  makes the tutorial's numbers the goldens' and the `r` test's numbers. **(user)**
- **Live cells:** load, glom, filter, `host`; `ancombc2` and `linda` with their call counts; the
  36 genera ANCOM-BC2 cannot test; `da.consensus`; `n_significant` counts; consensus and conflict
  totals; a `pd.crosstab` of the two methods' calls; `pl.consensus` of the top 30 (genus labels),
  and `pl.consensus` of the genera only one method calls (hollow dots visible). No numeric claim
  in the prose that the outputs do not show, except data facts pinned by golden tests (636 genera,
  36 untested, 9 against 17).
- **The R bridges: a fenced, non-executed ` ```python ` block with the four-method code, then two
  static tables and one sentence of the numbers the `r-bridge` job measures. (user)** Read the
  Docs and the `docs` job have no R. The numbers cannot drift silently:
  `tests/da/test_consensus.py` gains `FOUR_METHOD_CALLS`, `FOUR_METHOD_N_SIGNIFICANT` and
  `FOUR_METHOD_CONSENSUS`; the existing `r` test `test_four_methods_on_the_exit_gate_data` asserts
  against them, and a new default-run test `test_the_tutorial_quotes_the_four_method_counts`
  asserts the tutorial's table rows and sentence render them. A change in any method's numbers
  fails the `r` test; fixing the constants then fails the tutorial test until the page matches.
  No silent fallback (R7.4): the page says the block is not run and why.
  Rejected: (a) "recorded outputs" through a `.ipynb` with stored outputs or a `skip-execution`
  tag: per-page execution metadata is forbidden by
  `test_no_page_overrides_the_notebook_execution_settings`, the format would differ from every
  other tutorial, and stored outputs are not checked; (b) a CSV written by the `r-bridge` job and
  committed for an executed cell to read: a new artifact path and a CI change for numbers that a
  constant pins as well; (c) dropping the mock samples so "other" is purely environmental: new R
  runs and constants, and the tutorial would no longer match the goldens.
- **No workflow change and no Read the Docs change. (user)** The `docs` job already caches pooch
  under the `_remote.py` key and executes every notebook; Read the Docs runs the same
  `hatch run docs:build` without a cache, as Phase 2 left it ("Read the Docs has no pooch cache":
  no workflow change). GlobalPatterns is already fetched in the same build by the phyloseq
  vignette, so the tutorial adds about 9 s and no download.
- `docs/methods/index.md` and the guide link the tutorial (added in 3.10b, because 3.10's build
  would fail on a `{doc}` link to a page that does not exist yet).

#### 3. Coming from R (3.12)

The seven Phase 3 rows already come from docstrings: `ANCOMBC::ancombc2`, `MicrobiomeStat::linda`,
`ALDEx2::aldex`, `maaslin3::maaslin3`, `philr::philr`, `vegan::decostand` and
`mia::transformAssay` (now `pp.clr`, `pp.philr`, `pp.relative`); the table has 58 rows. A
parametrized test pins them, as 2.13 pinned the mia importers. The five "not in 0.2" labels
become "not in 0.3" (cut-a-release step 2c), with the page sentence, the test and the phyloseq
vignette's three "0.2" statements. The vignette's bullet "multiple testing and differential
abundance (biotapy 0.3)" becomes true text: phyloseq's section points to DESeq2
(`phyloseq_to_deseq2`), which biotapy does not wrap; biotapy's own methods are in the new tutorial.
No `phyloseq_to_deseq2` row is added to `r_idioms.toml` (not asked for; R1.4).

#### 4. Benchmarks (3.13)

- **Synthetic data, measured only (R10.1). (user)** asv never downloads (Phases 1 and 2 used
  synthetic tables only; the lint job's `asv check` imports the suite on every PR). The outline's
  "pp.philr on GlobalPatterns" would need pooch inside asv; a synthetic table of chosen size says
  more about scaling anyway.
- `Philr` in `pp.py`: `synthetic(n_obs=1_000)`, the Phase 1 table's construction (50,000
  features, 2% non-zero, balanced tree whose one-child nodes PhILR skips) at 1,000 samples, so `X`
  densifies to 400 MB; the full 5,000 samples would densify 2 GB and peak near 10 GB.
- `DifferentialAbundance` in the new `da.py`: `synthetic_genera()`, 2,000 samples in two groups
  of 1,000 x 10,000 features, 30% non-zero counts 1-99 (a large genus-level study; zeros are what
  make ANCOM-BC2 treat values as missing), with `time_` and `peakmem_` for `linda` and
  `ancombc2`. ALDEx2 and MaAsLin 3 are not benchmarked: their time is R's, and asv has no R.
- Two commits, as 2.14: the suite (gated), then the page with the run's numbers, the commit hash
  and the load averages.
- **Fixed, approved 2026-10-07 (decision 8):** `da.ancombc2`'s docstring said peak memory is
  "about seven" dense arrays (measured on 400 x 1,000 and 200 x 2,000); on the asv table it is
  4.5x. The docstring now says 4.5x to 7.5x depending on the table's shape, in its own `docs(da)`
  commit; the page states both ends.

#### 5. Knowledge (3.14) and the release (3.15)

- **3.14 is mostly done:** `modules/da.md` was rewritten for the bridges at Checkpoint C (entry
  points, R bridge invariants, gotchas, the `-m r` verification). What remains: `da.md`
  Verification names the tutorial and method pages; its memory gotcha gains the asv number; one new
  gotcha (the tutorial's quoted counts and the shared constants); `add-a-function` step 8 (a `da`
  method also gets a `docs/methods/` page and a guide table row); verification bumps for the four
  concepts `knowledge_stale.sh --touched` lists (`phase-0-foundation`, `phase-1-core`,
  `engine-parity`, `r-golden-parity`) after reading each against the diff (nothing they state
  changed). `pl.md`, `pp.md`, `core.md` and `optional-heavy-dependencies.md` were re-read: no false
  statement (`pp.md`'s "philr peaks at about five arrays" agrees with the measured 4.3x plus the
  basis).
- **3.15** follows cut-a-release and Task 2.15. `## [Unreleased]` is empty again (no 3A-3D PR
  wrote an entry), so the 0.3.0 entries are written from `git log v0.2.0..master --no-merges`:
  `feat:` commits under Added; under Changed the `fix:` commits that change a 0.2.0 function
  (`dfa71fe`: readers infer counts only for non-negative whole numbers and count-needing functions
  refuse negatives; `c920c7b`: `pl.bar`/`pl.heatmap` refuse negative tables) and 3.12's labels.
  README: "0.3", the new Preprocessing, Differential abundance and Plots lines, the `r` extra with
  its GPL note under Installation, "Next, in 0.4" from Phase 4's description. **(user)**
- **Approvals.** Push, PR and merge-on-green of `phase-3d` are covered by the standing approval
  for Phase 3 slice branches (2026-10-05). `release-0.3.0`'s push/PR/merge, then the tag, the tag
  push, the GitHub release (the PyPI upload) and `close-phase-3` each need an explicit yes at the
  moment named (R13.3). **(user)**

#### 6. Which commit ticks which box

| Box | Ticked in | Evidence |
|---|---|---|
| Checkpoint C's push/PR/merge box | the plan commit | PR #23 merged as `29efd9c`, Test run `37479113322` 21/21 green |
| Checkpoint C's "ask the user to review 3C" box and the checklist's "Checkpoint C" | the plan commit | the user's approval of this slice 3D plan, which follows the 3C review |
| Exit gate 2, per-method agreement with R | the plan commit | PR #23: `network` job `36 passed` (CLR, PhILR, LinDA, ANCOM-BC2 goldens); `r-bridge` job `24 passed` (ALDEx2, MaAsLin 3 goldens) |
| 3.10, 3.10b, 3.12, 3.13, 3.14 checklist lines | each task's last commit | the gates above |
| Checkpoint D (checklist and its boxes) | the 3.15 knowledge commit (`release-0.3.0`) | the `phase-3d` PR merged green |
| Exit gate 1, executed notebook in the CI docs job | the 3.15 knowledge commit | that PR's `docs` job log: `tutorials/differential_abundance.md: Executed notebook` |
| Exit gate 3, Phase 1 and 2 gates green | the 3.15 knowledge commit | every job of that PR's Test run green |
| 3.15 Steps 1-8 | the 3.15 knowledge commit | |
| 3.15 Steps 9-13, `phase_state: done` | `close-phase-3`, after the upload | the release workflow run |

### Slice 3D global constraints (in addition to the Phase 3 list)

- No `src/` change (except the two approved `docs(da)` docstring commits (`3b6239f`, `b9cb006`) and the Checkpoint D fix pass's docstring lines) and no new dependency (runtime, extra or dev). New calls, each checked in the
  installed versions: none in Python beyond what the tutorial uses (`pd.crosstab`, already pandas
  3.0.6); docs syntax: MyST `{code-cell}`, `{doc}`, `{ref}`, `(label)=` targets and `$$` math.
- Every `pytest`, Python, asv and Sphinx run exports an **absolute**
  `BIOTAPY_DATA_DIR=<scratchpad>/pooch`; `~/.cache/biotapy` must not appear, nor `docs/pooch/`.
- Delete `docs/_build` and `docs/generated` before every docs build used as a gate.
- The tutorial executes on every docs build (CI and Read the Docs); it must never call an R
  bridge in a `{code-cell}`.
- Method pages state only what a golden test, a docstring or the R function's `formals()` backs;
  numbers are the golden tests' measured values.
- Commits stage explicit paths only (never `.claude/`, `.superpowers/`, `.worktrees/`,
  `notebooks/`, `build/`, `dist/`). Each task's last commit stages
  `.knowledge/roadmap/phase-3-stats.md` with its checklist line ticked and `.knowledge/log.md`
  with its line under `## <date of the commit> (Phase 3, slice 3D)`, newest line first.
- Diff blocks below are for reading, not `git apply`.
- Branch `phase-3d` from `master` (`29efd9c`). The release branch, the tag, the GitHub release and
  the phase-closing branch each need the user's explicit approval at that moment (R13.3).

### Slice 3D review focus

The five ways a user is most likely to get a wrong answer, or a broken release, from this slice
without an error, most likely first, each pinned by a test or a named check:

1. **The tutorial quotes R-bridge numbers that are no longer what the bridges give.** The docs
   cannot run R. Expected: the quoted counts and the `r` test read the same constants. Tests
   (3.10b): `tests/da/test_consensus.py::test_the_tutorial_quotes_the_four_method_counts` and the
   `r` test `test_four_methods_on_the_exit_gate_data`.
2. **A method page states an R default or setting biotapy does not use**, so a user reproduces it
   in R with different arguments. Expected: every row of each "Compared with R" table matches the
   R function's `formals()` (recorded above) and the call in the docstring and in
   `tests/r/export_golden.R`. Check: Checkpoint D's review reads each table against those three.
3. **A user misreads the contrast.** The tutorial's "other" includes three mock communities;
   `effect` is "human vs other". Expected: the note says what "other" holds, and the code passes
   `reference="other"` to every method (consensus refuses mixed contrasts:
   `tests/da/test_consensus.py::test_results_with_different_contrasts_raise`).
4. **The Coming-from-R table maps an R function to the wrong biotapy function.** Expected: the
   seven Phase 3 rows exactly. Test (3.12):
   `tests/test_coming_from_r.py::test_transforms_and_da_methods_map_to_their_r_functions`.
5. **The 0.3.0 wheel lacks the `r` extra or the README promises what 0.3 does not ship.**
   Expected: METADATA has `Provides-Extra: r` and `Requires-Dist: rpy2>=3.6.8; extra == 'r'`; the
   README lists only `__all__` members. Check (3.15 Step 6): the METADATA grep and, after the
   upload, Step 12's `bt.da.__all__` print.

Execution order: **plan commit -> 3.10 -> 3.10b -> 3.12 -> 3.13 -> Checkpoint D (review, fix pass,
re-review, 3.14 knowledge, gates, PR) -> 3.15.** The tutorial links the method pages and the
guide's `da-methods-in-r` label (3.10); 3.12 links the tutorial from the vignette; the benchmarks
share no file with the docs; the knowledge pass documents the final state; the release comes
last.

**Plan commit (controller, before Task 3.10).** With this section, the plan commit also updates,
outside it:
- the header note: "Slice 3A has full TDD steps; slices 3B-3D are outlines, expanded into TDD
  steps (rules.md R1.2a) when each is reached." becomes "Every slice has full TDD steps.";
- "# Tasks (checklist)": `- [ ] Checkpoint C (review, exit-gate check and knowledge done; push
  open)` becomes `- [x] Checkpoint C (PR #23 merged; the user approved slice 3C on 2026-10-07)`, and
  `- [ ] 3.14 Knowledge (updates ...)` becomes
  `- [ ] 3.14 Knowledge (the DA docs and benchmarks; da.md was updated for the bridges at Checkpoint C)`;
- Checkpoint C: tick the push box, adding "Done 2026-10-06: PR #23 merged as `29efd9c`, Test run
  `37479113322` 21/21 green (`r-bridge` `24 passed`, `network` `36 passed`)", and the review-ask
  box with the user's approval date;
- "# Exit gate": tick the second item, adding "Proven on PR #23 (run `37479113322`): `network`
  job `36 passed`, `r-bridge` job `24 passed`.";
- design note 10: "(Feces, Skin, Tongue: 8) vs environmental (18)" becomes "(Feces, Skin,
  Tongue: 9) vs the 17 others (14 environmental, 3 mock communities)", and "the ALDEx2 and MaAsLin
  3 calls are shown as a non-executed code block with their extra" gains "and the numbers the
  `r-bridge` job checks, pinned by a test";
- "# Decisions for the user", item 16 gains "(slice 3D decision 2)";
- `.knowledge/log.md`, under `## 2026-10-07 (Phase 3, slice 3D plan)`:
  `- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) expands slice 3D (3.10 method pages and guide, 3.10b tutorial, 3.12 Coming-from-R, 3.13 benchmarks, Checkpoint D with 3.14 knowledge, 3.15 release 0.3.0) into full TDD steps, prototyped and gated per commit; ticks Checkpoint C (PR #23 merged, slice 3C approved) and exit-gate item 2 (PR #23's network and r-bridge jobs).`

---
### Task 3.10: Method pages and the DA guide

**Files:**
- Create: `docs/methods/index.md`, `docs/methods/linda.md`, `docs/methods/ancombc2.md`,
  `docs/methods/aldex2.md`, `docs/methods/maaslin3.md`, `docs/methods/consensus.md`.
- Modify: `docs/guide/differential_abundance.md`, `docs/index.md`.

**Not touched:** any `src/` file (the docstrings' `Guide:` lines keep pointing at the guide; the
docstring-prose fixes the ledger approved are the exception),
`docs/api.md` (it already lists the five `da` functions and `pl.consensus`), `docs/guide/transforms.md`
(CLR and PhILR stay documented there; the outline's method pages are for DA), `docs/references.bib`,
`docs/conf.py`, the tests.

**Interfaces:**
- Consumes: the docstrings' Notes of `da.linda`, `da.ancombc2`, `da.aldex2`, `da.maaslin3`,
  `da.consensus`; the golden tests' measured tolerances (`tests/da/test_*_golden.py`,
  [r-golden-parity](/contracts/r-golden-parity.md)); the decision
  [da-consensus-agreement](/decisions/da-consensus-agreement.md); the R functions' `formals()`.
- Produces: the pages `methods/*`; the label `da-methods-in-r` on the guide's "Methods that run
  in R" section, which 3.10b's tutorial and the ALDEx2 and MaAsLin 3 pages reference.

- [x] **Step 1: The method pages.** Create `docs/methods/index.md`:

````markdown
# Methods

One page per differential abundance method: the model it fits, the units of its `effect`, how
biotapy's call differs from the R function's defaults, how closely the two agree on real data,
and the paper to cite. The {doc}`differential abundance guide </guide/differential_abundance>`
says how to run them.

```{toctree}
:maxdepth: 1

linda
ancombc2
aldex2
maaslin3
consensus
```
````

  `docs/methods/linda.md`:

````markdown
# LinDA

{func}`bt.da.linda <biotapy.da.linda>` fits one linear model per feature to the log2 centred
log-ratios of the counts and removes the bias that compositionality adds to every coefficient.
It runs in Python; no R is needed.

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.obs["age"] = [30, 41, 52, 38, 45, 60]
table = bt.da.linda(tdata, "group", covariates=["age"])
table[table["qvalue"] < 0.05]
```

## Model

For counts $y_{ij}$ of feature $j$ in sample $i$ ($n$ samples, $m$ features), 0.5 is added to
every count when the table holds a zero. Each sample's log2 counts are centred:

$$
w_{ij} = \log_2 y_{ij} - \frac{1}{m} \sum_{k=1}^{m} \log_2 y_{ik}
$$

Each feature gets one least-squares fit $w_{\cdot j} = X \beta_j + \varepsilon_j$ on the design
$X$ (an intercept, `group` and the `covariates`), with $p$ columns. The group coefficient
$\tilde\beta_j$ estimates the log2 fold change plus a bias $b$ that is the same for every feature,
because the centring divides by the geometric mean of a composition. LinDA estimates the bias as
the mode of the coefficients and subtracts it:

$$
\hat b = \frac{1}{\sqrt n}\, \operatorname{mode}\{\sqrt n\, \tilde\beta_1, \dots, \sqrt n\, \tilde\beta_m\},
\qquad \mathrm{effect}_j = \tilde\beta_j - \hat b
$$

The mode is found by a Gaussian mean shift started at the mean of the shortest half of the values
(`shorth`), with `bw.nrd0`'s bandwidth, as `modeest::mlv(method = "meanshift")` does. The
p-value is a two-sided t-test of $\mathrm{effect}_j / \mathrm{se}_j$ on $n - p$ degrees of
freedom, where $\mathrm{se}_j$ is the least-squares standard error of $\tilde\beta_j$, and
`qvalue` is the Benjamini-Hochberg correction over all features. Most features end up near zero;
the ones that changed stand out.

## Units

`effect` is a log2 fold change of the other level over `reference`. Numeric columns, a numeric
`group` included, are scaled to unit variance first, as LinDA does, so a numeric group's effect
is per standard deviation.

## Compared with R

biotapy equals `MicrobiomeStat::linda(t(counts), meta, "~group + covariates", feature.dat.type =
"count", is.winsor = FALSE)` (MicrobiomeStat 1.4) with fixed effects:

| `linda` argument | R default | biotapy |
|---|---|---|
| `formula` | required | `group` and `covariates`, fixed effects only: no `(1 \| subject)` |
| `feature.dat.type` | `"count"` | counts only: `X` must hold raw counts |
| `prev.filter`, `mean.abund.filter`, `max.abund.filter` | `0` | no filter: run `bt.pp.filter_features` once, before any method |
| `is.winsor`, `outlier.pct` | `TRUE`, `0.03` | no winsorisation |
| `adaptive`, `zero.handling`, `pseudo.cnt` | `TRUE`, `"pseudo-count"`, `0.5` | 0.5 added to every count when the table has a zero, which is what R runs (below) |
| `p.adj.method` | `"BH"` | Benjamini-Hochberg |
| `alpha` | `0.05` | not an argument: `bt.da.consensus` makes the calls |

MicrobiomeStat 1.4 announces an imputation of zeros when library size depends on the model, but
its switch compares `"Imputation"` with `"imputation"` and never fires, so it always adds the
pseudocount. biotapy computes what R returns.

## Agreement with R

The golden test compares `effect`, `se`, `pvalue` and `qvalue` with `MicrobiomeStat::linda` on
the 636 GlobalPatterns genera in at least 20% of samples, human hosts against the rest, with and
without a sequencing-depth covariate: they agree to a relative 1e-7.

## Reference

Zhou H, He K, Chen J, Zhang X (2022) LinDA: linear models for differential abundance analysis of
microbiome compositional data. *Genome Biology* 23:95.
[doi:10.1186/s13059-022-02655-5](https://doi.org/10.1186/s13059-022-02655-5)
````

  `docs/methods/ancombc2.md`:

````markdown
# ANCOM-BC2

{func}`bt.da.ancombc2 <biotapy.da.ancombc2>` runs scikit-bio's ANCOM-BC2, which corrects each
sample's log counts for the share of its ecosystem that was sequenced before it fits one linear
model per feature. It runs in Python; no R is needed.

```python
import biotapy as bt

table = bt.da.ancombc2(bt.datasets.toy(), "group")
```

## Model

The log of an observed count is the log of the feature's abundance in the ecosystem plus the log
of the sample's sampling fraction $\theta_i$, which differs between samples and biases every
comparison of raw counts:

$$
\log o_{ij} = \theta_i + \log a_{ij}, \qquad \log a_{ij} = x_i^\top \beta_j + \varepsilon_{ij}
$$

ANCOM-BC2 (scikit-bio's {func}`~skbio.stats.composition.ancombc2`):

1. treats zeros as missing (no pseudocount) and fits one least-squares model per feature to the
   log counts on the design $x_i$ (an intercept, `group` and the `covariates`);
2. estimates the bias the coefficients of each term share with an E-M algorithm that models them as
   a mixture of a null, a negative and a positive normal component, at most 100 iterations;
3. estimates each sample's $\theta_i$ from the bias-corrected fit, subtracts it from the sample's
   log counts and fits the models again;
4. widens every variance by the bias estimate's own variance and then by the 5% quantile of all
   variances, so that rare features with tiny standard errors do not dominate, and tests
   $W_j = \hat\beta_j / \mathrm{se}_j$ two-sided against a t distribution whose degrees of freedom
   are the feature's observed samples minus the model's rank.

`qvalue` is the Benjamini-Hochberg correction over the features ANCOM-BC2 could test. A feature
whose zeros leave one level of `group` without an observed value, or that has no more observed
samples than model terms, cannot be tested: its row is NaN.

## Units

ANCOM-BC2 reports natural logs; biotapy divides `effect` and `se` by ln 2, so they are log2 like
every other method's. A numeric `group`'s effect is per unit, where `da.linda`'s and
`da.maaslin3`'s are per standard deviation.

## Compared with R

scikit-bio's implementation matches `ANCOMBC::ancombc2` (ANCOMBC 2.12.0) with these settings:

| `ancombc2` argument | R default | biotapy |
|---|---|---|
| `fix_formula` | required | `group` and `covariates` |
| `rand_formula` | `NULL` | not offered: no random effects |
| `p_adj_method` | `"holm"` | Benjamini-Hochberg |
| `prv_cut`, `lib_cut` | `0.10`, `0` | `0`, `0`: run `bt.pp.filter_features` once, before any method |
| `pseudo` | `0` | `0`: zeros are missing |
| `pseudo_sens` | `TRUE` | not run, so the table has no `passed_ss` flag |
| `s0_perc` | `0.05` | `0.05` |
| `struc_zero`, `neg_lb` | `FALSE` | not run |
| `em_control` | `tol = 1e-5`, `max_iter = 100` | the same |
| `global`, `pairwise`, `dunnet`, `trend` | `FALSE` | not offered: `group` has two levels or is numeric |
| `alpha` | `0.05` | not an argument: `bt.da.consensus` makes the calls |

R reports a feature it cannot test with p = 1 and counts it in the correction; biotapy reports
NaN and leaves it out, so its q-values are BH over the features actually tested.

## Agreement with R

The golden test compares biotapy with `ANCOMBC::ancombc2` on the 636 GlobalPatterns genera in at
least 20% of samples, human hosts against the rest. With a sequencing-depth covariate, `effect`,
`se` and `pvalue` agree to 1e-6. Without it, the bias E-M stops at R's cap of 100 iterations
before it has converged, on a slightly different iterate in scikit-bio: effects agree within
0.015 log2, standard errors to a relative 2e-3 and p-values within 0.02, and the significant
genera are the same. With 1,000 iterations R moves every effect by about -0.28 log2 and calls 220
genera instead of 208; biotapy keeps R's default.

## Choosing the reference

Unlike LinDA, ANCOM-BC2 is not antisymmetric in `reference`: the bias-corrected fit is made
against the reference level, in R as here, so swapping it changes more than the sign of
`effect`. On the GlobalPatterns genera the calls at q < 0.05 go from 208 to 230, and an effect
plus its swap is about -0.38 log2, not 0; `da.consensus` with LinDA goes from 104 to 112 genera.
Choose `reference` on the biology (the control or baseline level), not to change the results.

## Reference

Lin H, Peddada SD (2024) Multigroup analysis of compositions of microbiomes with covariate
adjustments and repeated measures. *Nature Methods* 21:83-91.
[doi:10.1038/s41592-023-02092-7](https://doi.org/10.1038/s41592-023-02092-7)
````

  `docs/methods/aldex2.md`:

````markdown
# ALDEx2

{func}`bt.da.aldex2 <biotapy.da.aldex2>` runs `ALDEx2::aldex` in R through rpy2: it draws
plausible relative abundances for each sample, takes their log2 centred log-ratios and compares
the two groups with a Welch t-test on every draw. It needs R, the R package ALDEx2 and biotapy's
`r` extra ({ref}`Methods that run in R <da-methods-in-r>`).

```python
import biotapy as bt

table = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
```

## Model

For each sample $i$, ALDEx2 draws `mc_samples` proportion vectors from a Dirichlet distribution
whose parameters are the sample's counts plus 0.5, so a draw reflects how uncertain proportions
are at that sequencing depth. Each draw is turned into log2 centred log-ratios over all features
(`denom = "all"`):

$$
c^{(k)}_{ij} = \log_2 p^{(k)}_{ij} - \frac{1}{m} \sum_{l=1}^{m} \log_2 p^{(k)}_{il},
\qquad p^{(k)}_{i\cdot} \sim \operatorname{Dirichlet}(y_{i\cdot} + 0.5)
$$

For each draw $k$ and feature $j$, a Welch t-test compares the two groups' values $c^{(k)}_{ij}$;
`pvalue` is ALDEx2's `we.ep`: for each draw ALDEx2 doubles the one-sided Welch p-value in each
direction (capped at 1), averages each direction over the draws, and keeps the smaller average.
`effect` is ALDEx2's `diff.btw`, the median difference between values resampled at random from each
group's pooled draws, already log2. ALDEx2
reports no standard error, so `se` is NaN for every feature, and `qvalue` is the
Benjamini-Hochberg correction of `pvalue`, as for every method. A feature with no read in any
sample is dropped by ALDEx2 and is NaN in the table (not tested).

## Units

`effect` is the median log2 difference of the centred log-ratios between the other level and
`reference`. ALDEx2 compares two groups only: a numeric `group`, covariates and a level in fewer
than two samples are refused.

## Compared with R

biotapy calls `ALDEx2::aldex` (ALDEx2 1.42.0) with its defaults:

| `aldex` argument or output | R default | biotapy |
|---|---|---|
| `conditions` | required | the labels `"0"` (`reference`) and `"1"`: ALDEx2 fails on a factor and sorts character labels in the R session's locale |
| `mc.samples` | `128` | `mc_samples=128`; R warns below 128, and biotapy re-emits it as a `UserWarning` |
| `test`, `effect`, `denom` | `"t"`, `TRUE`, `"all"` | the same |
| `paired.test`, `iterate`, `gamma` | `FALSE`, `FALSE`, `NULL` | the same: no paired test, no scale model |
| `we.ep` | output | `pvalue` |
| `we.eBH` | output | not carried: it averages the draws' corrected values, a different quantity; `qvalue` is BH of `we.ep` |
| `diff.btw` | output | `effect` |
| `effect`, `overlap`, `wi.ep`, `wi.eBH`, `rab.*`, `diff.win` | outputs | not carried: ALDEx2's `effect` is a standardised size, not a fold change |
| random state | R's global seed | `seed` sets R's seed for the call, and the R session's random state is put back afterwards |

On the GlobalPatterns genera below, BH of `we.ep` calls 11 genera and ALDEx2's `we.eBH` 19.

## Agreement with R

ALDEx2 is random, but `seed` becomes one integer for R's `set.seed`, so biotapy and R draw the same
numbers. The golden test runs `set.seed(...); ALDEx2::aldex(...)` with that integer on the 636
GlobalPatterns genera in at least 20% of samples, human hosts against the rest: `effect` and
`diff.btw` agree to a relative 4e-15 and `pvalue` and `we.ep` to 8.4e-13 (checked at 1e-7). The
golden passes `reference="human"`, the level R sorts first, so both sides resample `diff.btw` in the same
label order. With other seeds, effects correlate
at 0.993 (Spearman) and 11 to 13 genera are called, depending on the seed (13 at `seed=0`, the
tutorial's).

## Choosing the reference

Swapping `reference` does more than flip the sign. ALDEx2 draws the same Monte Carlo instances, so
`pvalue` and the calls are unchanged, but `diff.btw` comes from a random resampling of each group's
values done in label order, so the same `seed` gives different effects, and the two runs are not
exact mirror images. They differ from exact antisymmetry by 0.1 to 0.3 log2 on `toy()`, depending
on the seed (its effects are about 4 log2 wide), and by up to 0.65 log2 on the GlobalPatterns
genera, 15 of 636 of which then do not change direction. The p-values and the calls are the same.

## Reference

Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB (2014) Unifying the
analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing
and selective growth experiments by compositional data analysis. *Microbiome* 2:15.
[doi:10.1186/2049-2618-2-15](https://doi.org/10.1186/2049-2618-2-15)
````

  `docs/methods/maaslin3.md`:

````markdown
# MaAsLin 3

{func}`bt.da.maaslin3 <biotapy.da.maaslin3>` runs the abundance model of `maaslin3::maaslin3` in R
through rpy2: one linear model of each feature's log2 relative abundance, on the samples where the
feature is present, tested against the median coefficient. It needs R, the R package maaslin3 and
biotapy's `r` extra ({ref}`Methods that run in R <da-methods-in-r>`).

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.obs["age"] = [30, 41, 52, 38, 45, 60]
table = bt.da.maaslin3(tdata, "group", covariates=["age"], seed=0)
```

## Model

Counts are divided by their sample's total (total sum scaling). For feature $j$, on the samples
where it is present ($y_{ij} > 0$):

$$
\log_2 \frac{y_{ij}}{\sum_k y_{ik}} = x_i^\top \beta_j + \varepsilon_{ij}
$$

with numeric columns of the design $x_i$ standardised. Compositionality shifts every coefficient
of a term by about the same amount, so MaAsLin 3 tests each coefficient against the median
coefficient $\tilde\beta$ of that term, taken over the features without a fit error whose own
p-value is below 0.95; the test draws 10,000 normal samples, which is why `seed` matters.
biotapy reports

$$
\mathrm{effect}_j = \hat\beta_j - \tilde\beta
$$

(`subtract_median = TRUE`), so the sign of `effect` is the side of the median the p-value is
about. `qvalue` is the Benjamini-Hochberg correction of the group's p-values. A feature whose fit
reports an error is NaN in the table (not tested), as MaAsLin 3 leaves it out of its own
correction.

## Units

`effect` is a log2 fold change of the relative abundance where the feature is present, minus the
median. Numeric columns are standardised, so a numeric `group`'s effect is per standard
deviation, as in `da.linda`.

## Compared with R

biotapy calls `maaslin3::maaslin3` (maaslin3 1.2.0) with:

| `maaslin3` argument | R default | biotapy |
|---|---|---|
| `formula` | `NULL` | `~ x0 + x1 + ...`, the `group` and `covariates` columns under plain names, so a column such as `"body site"` needs no quoting |
| `reference` | `NULL` | not passed: categorical columns arrive as factors with string levels, `reference` first |
| `random_effects`, `group_effects`, `ordered_effects`, `strata_effects` | `NULL` | not offered |
| `min_abundance`, `min_prevalence` | `0`, `0` | the same: no filter |
| `normalization`, `transform`, `standardize` | `"TSS"`, `"LOG"`, `TRUE` | the same |
| `median_comparison_abundance` | `TRUE` | the same |
| `subtract_median` | `FALSE` | `TRUE`: `effect` is what the p-value tests |
| `evaluate_only` | `NULL` (both models) | `"abundance"`: the prevalence model's log-odds cannot share a column with fold changes |
| `warn_prevalence` | `TRUE` | `FALSE`, which `evaluate_only` requires |
| `correction` | `"BH"` | BH over the group's p-values only; MaAsLin 3's `qval_individual` corrects them together with every covariate's |
| `plot_summary_plot`, `plot_associations` | `TRUE` | `FALSE`; the `output` folder is temporary and deleted |
| `cores`, `verbosity` | `1`, `"FINEST"` | `1`, `"ERROR"` |

With one numeric covariate on the GlobalPatterns genera below, `qval_individual` calls 20 genera
where biotapy's `qvalue` calls 45; without covariates the two agree (52).

## Agreement with R

MaAsLin 3's median test simulates, but `seed` becomes one integer for R's `set.seed`, so biotapy
and R draw the same numbers. The golden test runs `set.seed(...); maaslin3::maaslin3(...)` with
that integer on the 636 GlobalPatterns genera in at least 20% of samples, human hosts against the
rest, with and without a sequencing-depth covariate: the same 36 genera fail to fit, and `effect`,
`se` and `pvalue` agree with R's `coef`, `stderr` and `pval_individual` to a relative 7.2e-13 or
better (checked at 1e-7). With other seeds only the p-values move, by up to 0.0015.

## Choosing the reference

Swapping `reference` negates `effect` exactly, but the median test's simulation draws around the
coefficients rather than their negatives, so the same `seed` moves `pvalue` by up to 1e-3 and
`qvalue` by up to 2e-3 on the GlobalPatterns genera (the calls are the same).

## Reference

Nickols WA, Kuntz T, Shen J, Maharjan S, Mallick H, Franzosa EA, Thompson KN, Nearing JT,
Huttenhower C (2026) MaAsLin 3: refining and extending generalized multivariable linear models
for meta-omic association discovery. *Nature Methods* 23:554-564.
[doi:10.1038/s41592-025-02923-9](https://doi.org/10.1038/s41592-025-02923-9)
````

  `docs/methods/consensus.md`:

````markdown
# Consensus

{func}`bt.da.consensus <biotapy.da.consensus>` puts the result tables of several methods side by
side and counts, per feature, the methods that call it and whether they agree on its direction.
It runs no method itself: you choose and run the methods, and pass their tables.
{func}`bt.pl.consensus <biotapy.pl.consensus>` draws the table it returns.

```python
import biotapy as bt

tdata = bt.datasets.toy()
results = [bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")]
table = bt.da.consensus(results, alpha=0.05, min_methods=2)
```

## Rule

For feature $f$ and method $m$, with $q_{fm}$ the method's `qvalue` and $d_{fm}$ its `direction`:

| Column | Definition |
|---|---|
| `significant_<m>` | $q_{fm} < \alpha$, strictly: $q = \alpha$ is not a call |
| `n_tested` | methods whose `pvalue` for $f$ is finite; a method that could not test $f$, or whose table lacks it, does not count |
| `n_significant` | methods that call $f$ |
| `direction` | the sign $d_{fm}$ every calling method shares; 0 when none calls $f$ or their signs differ |
| `consensus` | `n_significant` $\ge$ `min_methods` and `direction` is not 0 |
| `conflict` | calling methods give $f$ opposite signs; never a consensus |

A call whose `effect` is exactly 0 has no direction: it counts in `n_significant`, but then
`direction` is 0, with no consensus and no conflict. The table also carries each method's
`effect_<m>` and `qvalue_<m>`, in the order of `results`.

## Why it is defined this way

- **One correction.** Every method's `qvalue` is Benjamini-Hochberg over the features it tested,
  so "significant" means the same false discovery rate in every column. R's defaults would mix
  Holm (ANCOM-BC2), BH over the group's p-values (LinDA), BH pooled with the covariates'
  (MaAsLin 3) and averaged per-draw BH (ALDEx2).
- **Strict calls, recomputed.** Calls come from `qvalue` with one rule, never from a method's own
  flag: LinDA's `reject` uses $q \le \alpha$, and ANCOM-BC2's `Signif` sits on Holm by default.
- **Untested is not "not significant".** ANCOM-BC2 and MaAsLin 3 cannot fit a feature that is
  absent from one group; counting that as a negative would let one method veto a feature it never
  tested. `min_methods=len(results)` still asks for every method.
- **Same contrast.** The tables must come from different methods and compare the same `contrast`;
  `bt.da.consensus` raises otherwise, so run every method with the same `group` and `reference`.
- **The user picks the methods.** A one-call consensus would hide which methods ran and invite
  trying methods until one agrees. Methods disagree a lot on real data: Nearing et al. (2022)
  recommend a consensus of several methods, Pelto et al. (2025) one simple method. Either way,
  choose the methods before you look at their results. Methods that share a model, such as
  ANCOM-BC and ANCOM-BC2, agree more often for that reason alone.

The rule is recorded, with the options rejected, in the knowledge bundle's decision
[da-consensus-agreement](https://github.com/pedrocr83/biotapy/blob/master/.knowledge/decisions/da-consensus-agreement.md).
There is no R equivalent.

## References

Nearing JT, Douglas GM, Hayes MG, MacDonald J, Desai DK, Allward N, Jones CMA, Wright RJ, Dhanani
AS, Comeau AM, Langille MGI (2022) Microbiome differential abundance methods produce different
results across 38 datasets. *Nature Communications* 13:342.
[doi:10.1038/s41467-022-28034-z](https://doi.org/10.1038/s41467-022-28034-z)

Pelto J, Auranen K, Kujala JV, Lahti L (2025) Elementary methods provide more replicable results
in microbial differential abundance analysis. *Briefings in Bioinformatics* 26:bbaf130.
[doi:10.1093/bib/bbaf130](https://doi.org/10.1093/bib/bbaf130)
````

  Where each statement comes from (re-check each against the code after any fix, R0.2):
  - the R defaults: `formals()` of `MicrobiomeStat::linda` (1.4), `ANCOMBC::ancombc2` (2.12.0),
    `ALDEx2::aldex` (1.42.0) and `maaslin3::maaslin3` (1.2.0) in the golden image
    (`docker run --rm --network none biotapy-golden Rscript -e 'formals(...)'`);
  - what biotapy passes: `_linda.py`, `_ancombc.py` (`p_adjust=None`, BH in `_schema.result`),
    `_aldex2.py:_ALDEX`, `_maaslin3.py:_MAASLIN`;
  - ANCOM-BC2's steps: scikit-bio 0.7.4 `skbio/stats/composition/_ancombc.py`
    (`_estimate_bias_em`, `_sample_fractions`, `_adjust_variance`, `_calc_pvalues`);
  - the agreement numbers: `tests/da/test_linda_golden.py` (`rtol=1e-7`),
    `test_ancombc_golden.py` (`MODELS`), `test_aldex2_golden.py`, `test_maaslin3_golden.py` and the
    slice 3B/3C measurements in this concept;
  - the reference-swap paragraphs: moved from the guide and the docstrings' Notes, unchanged.
- [x] **Step 2: The guide keeps the workflow.** `docs/guide/differential_abundance.md` (the four
  method sections move into the pages; a "Choosing methods" table links them; the generic bridge
  behaviour moves up into "Methods that run in R", which gets the label `da-methods-in-r`):

````diff
diff --git a/docs/guide/differential_abundance.md b/docs/guide/differential_abundance.md
index 9e8d2ad..7cb7117 100644
--- a/docs/guide/differential_abundance.md
+++ b/docs/guide/differential_abundance.md
@@ -42,51 +42,29 @@ filter features: filter once, with `bt.pp.filter_features`, before running any m
 method tests the same features. Methods need raw counts in `X` and a read in every sample
 (`bt.pp.filter_samples(tdata, min_depth=1)` drops empty ones).

-## LinDA
+## Choosing methods

-`bt.da.linda` fits one linear model per feature to the log2 centred log-ratios of the counts (with
-0.5 added to every count when the table has a zero). Compositionality biases every coefficient by
-the same amount; LinDA estimates that bias as the mode of all features' coefficients and subtracts
-it, so most features end up near zero and the ones that changed stand out:
+biotapy has four methods. They model compositional counts differently, so they call different
+features, and no one of them is right on every dataset:

-```python
-tdata.obs["age"] = [30, 41, 52, 38, 45, 60]
-table = bt.da.linda(tdata, "group", covariates=["age"])
-table[table["qvalue"] < 0.05]
-```
-
-It equals `MicrobiomeStat::linda(..., is.winsor = FALSE)` in R with fixed effects. Numeric columns
-are scaled to unit variance first, as LinDA does, so a numeric group's effect is per standard
-deviation. Not available: winsorisation (MicrobiomeStat's default), random effects such as
-`(1 | subject)`, and LinDA's own prevalence filters.
-
-## ANCOM-BC2
+| Method | Model | Runs in | Covariates | Numeric `group` | `seed` |
+|---|---|---|---|---|---|
+| {doc}`bt.da.linda </methods/linda>` | least squares of log2 centred log-ratios, minus the mode of all coefficients | Python | yes | per standard deviation | no |
+| {doc}`bt.da.ancombc2 </methods/ancombc2>` | least squares of log counts corrected for each sample's sampling fraction | Python (scikit-bio) | yes | per unit | no |
+| {doc}`bt.da.aldex2 </methods/aldex2>` | Welch t-tests on Monte Carlo draws of centred log-ratios | R | no | no | yes |
+| {doc}`bt.da.maaslin3 </methods/maaslin3>` | least squares of log2 relative abundance where present, tested against the median | R | yes | per standard deviation | yes |

-`bt.da.ancombc2` runs scikit-bio's ANCOM-BC2 (Lin and Peddada 2024): it estimates each sample's
-sampling fraction, corrects the log counts for it and the coefficients for their shared bias, and
-fits one linear model per feature. Zeros are treated as missing rather than given a pseudocount,
-so a feature with no read in one of the groups cannot be fitted: its row is NaN and it is left out
-of the Benjamini-Hochberg correction (R's `ANCOMBC::ancombc2` reports it with p = 1 and counts it).
-
-```python
-table = bt.da.ancombc2(tdata, "group")
-```
+Each method's page gives its model, the R defaults biotapy changes, and how closely biotapy
+matches the R package on real data. Two things to know before you choose:

-ANCOM-BC2 reports natural logs; biotapy divides `effect` and `se` by ln 2, so they are log2 like
-every other method's; a numeric `group`'s effect is per unit, where `da.linda`'s is per standard
-deviation. The settings are R's `ancombc2(..., p_adj_method = "BH", prv_cut = 0,
-pseudo_sens = FALSE)`: biotapy does not run R's pseudocount sensitivity analysis (its
-`passed_ss` flag) or its 10% prevalence filter. On the GlobalPatterns genera that the golden tests
-use, biotapy's effects are within 0.012 log2 of R's and the significant genera are the same; the
-small difference comes from the bias estimate, whose iterations stop at R's cap of 100 before they
-have converged on that data, in R as in scikit-bio.
-
-Unlike `da.linda`, ANCOM-BC2 is not antisymmetric in `reference`: swapping it changes more than the
-sign of `effect`, because the bias-corrected E-M is fitted against the reference level, in R as
-here. On the GlobalPatterns genera (`host`) the calls at q < 0.05 go from 208 to 230, and `effect`
-plus its swap is about -0.38 log2, not 0; `da.consensus` with LinDA goes from 104 to 112 genera.
-Choose `reference` on the biology (the control or baseline level), not to change the results.
+- ANCOM-BC2 and MaAsLin 3 leave zeros out, so a feature with no read in one group cannot be
+  fitted: its row is NaN. LinDA adds 0.5 to every count when the table has a zero, and ALDEx2
+  draws around every count, so both test it.
+- `reference` changes more than the sign in ANCOM-BC2 (its bias correction is fitted against the
+  reference level) and ALDEx2 (its effect's random resampling follows the label order), and
+  slightly in MaAsLin 3 (its p-values move by up to 1e-3 with the same `seed`, the calls do not).
+  Choose it on the biology, the control or baseline level, before you look at any result.

+(da-methods-in-r)=
 ## Methods that run in R

 ALDEx2 and MaAsLin 3 exist only in R, so biotapy calls them there through
@@ -102,57 +80,17 @@ rpy2 is GPL-2.0-or-later and the R packages have their own licences; biotapy its
 any of them. Without rpy2 or the R package, the call raises an `ImportError` that names what to
 install. Each call converts `X` to a dense table once, because rpy2 has no sparse converter.

-### ALDEx2
-
-`bt.da.aldex2` runs `ALDEx2::aldex` (Fernandes et al. 2014): Monte Carlo draws from each sample's
-Dirichlet posterior, their log2 centred log-ratios, and a Welch t-test per draw, whose p-values are
-averaged. `effect` is ALDEx2's `diff.btw`, the median log2 difference between the two groups:
+Both methods draw random numbers in R. `seed` seeds R for the call, and your R session's own
+random state is restored afterwards, so the same seed gives the same table. Each R warning raised
+during the call, such as ALDEx2's for fewer than 128 `mc_samples`, is re-emitted as a Python
+`UserWarning`; an R error is raised as a `RuntimeError`. Repeated `var_names` or `obs_names`
+raise: call `adata.var_names_make_unique()` first.

 ```python
 table = bt.da.aldex2(tdata, "group", seed=0)
+table = bt.da.maaslin3(tdata, "group", seed=0)
 ```

-ALDEx2 compares two groups without covariates, and each group needs two samples. It is random:
-`seed` seeds R for the call (your R session's own random state is restored afterwards), so the same seed gives the same table, and on the GlobalPatterns
-genera biotapy's numbers equal R's `set.seed(...); aldex(...)` to a relative 8.4e-13 (when `reference` is R's
-first sorted level, and R is seeded with the integer biotapy derives from `seed`). `qvalue` is the
-Benjamini-Hochberg correction of ALDEx2's expected p-value `we.ep`, as for every method; ALDEx2's own
-`we.eBH` averages the corrections of the draws instead and calls more features (19 against 11 on
-those genera).
-
-Swapping `reference` does more than flip the sign: ALDEx2 takes its Monte Carlo draws in label order,
-so the same `seed` gives different effects, and the two runs are not exact mirror images: the effects
-differ from exact antisymmetry by up to 0.25 log2 on `toy()` (where they are about 4 log2 wide) and up to
-0.65 log2 on the GlobalPatterns genera, 15 of 636 of which then do not change direction. The p-values and
-the calls are the same on both. Each R warning raised during the call, such as the one
-for fewer than 128 `mc_samples`, is re-emitted as a Python `UserWarning`; an R error is raised as a
-`RuntimeError`. Repeated `var_names` or `obs_names` raise: call `adata.var_names_make_unique()` first.
-
-### MaAsLin 3
-
-`bt.da.maaslin3` runs MaAsLin 3's abundance model (`maaslin3::maaslin3` with `evaluate_only =
-"abundance"`): counts become relative abundances, zeros are left out, and one linear model of the
-log2 abundance per feature is fitted on the samples where the feature is present. Each coefficient
-is tested against the median coefficient, MaAsLin 3's correction for
-compositionality (the median over the features without a fit error whose own p-value is below
-0.95, as MaAsLin 3 computes it); `effect` is the coefficient minus that median, so its sign says on which side
-of the median the feature moved (MaAsLin 3 reports the coefficient itself unless asked to subtract):
-
-```python
-table = bt.da.maaslin3(tdata, "group", covariates=["age"], seed=0)
-```
-
-Numeric columns are standardised, so a numeric group's effect is per standard deviation, as in
-`da.linda`. The test against the median simulates, so `seed` makes the p-values reproducible; on the
-GlobalPatterns genera biotapy's numbers equal R's to 1e-12. The prevalence model, whose effects
-are log-odds rather than fold changes, is not run. `qvalue` corrects the group's p-values only;
-MaAsLin 3's `qval_individual` corrects them together with every covariate's, which with one numeric
-covariate called 20 genera where biotapy calls 45. R warnings and errors surface as for ALDEx2.
-
-Swapping `reference` negates `effect` exactly, but the median test's simulation draws around the
-coefficients rather than their negatives, so the same `seed` moves `pvalue` by up to 1e-3 and
-`qvalue` by up to 2e-3 on those genera (the calls are the same); `seed` fixes the draws.
-
 ## Where methods agree

 `bt.da.consensus` puts the tables of several methods side by side and counts, for each feature,
@@ -169,7 +107,8 @@ A feature is a consensus hit when at least `min_methods` methods call it and all
 the same sign. Methods that call it in opposite directions mark a `conflict`, which is never a
 consensus. A method that could not test a feature does not count against it: `n_tested` says how
 many methods tested each feature. The tables must come from different methods and compare the same
-`contrast`, so run every method with the same `group` and `reference`.
+`contrast`, so run every method with the same `group` and `reference`. The
+{doc}`consensus page </methods/consensus>` defines every column.

 Methods disagree a lot on real data, and the literature is split on what to do about it: Nearing
 et al. (2022) recommend a consensus of several methods, Pelto et al. (2025) one simple method. Either
````

  and the toctree:

````diff
diff --git a/docs/index.md b/docs/index.md
index 34a172c..2dc0253 100644
--- a/docs/index.md
+++ b/docs/index.md
@@ -16,6 +16,7 @@ tutorials/index.md
 :caption: User guide

 guide/index.md
+methods/index.md
 coming_from_r.md
 ```

````

- [x] **Step 3: Build and read.**
  `rm -rf docs/_build docs/generated && BIOTAPY_DATA_DIR=<scratchpad>/pooch uv run --group doc sphinx-build -W -b html docs docs/_build/html`
  -> `build succeeded.` (39.8 s in the prototype). Open `docs/_build/html/methods/linda.html`:
  the math renders (KaTeX), the `linda` table's first row reads `(1 | subject)` (the escaped pipe);
  `methods/aldex2.html` links to `guide/differential_abundance.html#da-methods-in-r`; the guide's
  "Choosing methods" table links the four pages. Every `{func}` link resolves (the build is
  nitpicky).
- [x] **Step 4: Bookkeeping.**
  - `.knowledge/roadmap/phase-3-stats.md`: `- [ ] 3.10 Method pages in \`docs/methods/\` and the DA guide`
    becomes `- [x] ...`.
  - `.knowledge/log.md`, a new heading at the top:

```markdown
## <date> (Phase 3, slice 3D)
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.10 done: `docs/methods/` holds one page per `da` method and the consensus (model, units, the R defaults biotapy changes, agreement with R, reference); the DA guide keeps the workflow and links to them.
```

- [x] **Step 5: Gate and commit.**

```bash
git add docs/methods/index.md docs/methods/linda.md docs/methods/ancombc2.md docs/methods/aldex2.md \
  docs/methods/maaslin3.md docs/methods/consensus.md docs/guide/differential_abundance.md docs/index.md \
  .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
git commit -m "docs: add a page per differential abundance method and link them from the guide

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
git status --short                                   # empty
uvx prek run --all-files                             # all Passed
uv run --group test pytest -q -W error::UserWarning  # 1230 passed, 53 deselected
uv run --group test pytest -q -m "golden or network" # 36 passed, 1247 deselected
rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html  # build succeeded.
```

### Task 3.10b: The exit-gate tutorial

**Files:**
- Create: `docs/tutorials/differential_abundance.md`.
- Modify: `docs/tutorials/index.md`, `docs/methods/index.md`, `docs/guide/differential_abundance.md`,
  `tests/da/test_consensus.py`.

**Not touched:** `docs/conf.py` (`nb_execution_timeout = 300` and the cache mode already hold),
the workflows and `.readthedocs.yaml` (design 2: no change), `tests/da/conftest.py` (the
`benchmark` fixture is the tutorial's data, unchanged), any `src/` file.

**Interfaces:**
- Consumes: `bt.datasets.global_patterns()`, `bt.pp.tax_glom(adata, rank)`,
  `bt.pp.filter_features(adata, *, min_prevalence)`, `bt.da.ancombc2` / `bt.da.linda`
  `(adata, group, *, covariates=(), reference=None)`, `bt.da.consensus(results, *, alpha=0.05,
  min_methods=2)`, `bt.pl.consensus(table, *, top=30, ax=None)`; the guide label
  `da-methods-in-r` (3.10).
- Produces: the page `tutorials/differential_abundance` (exit gate item 1); the constants
  `FOUR_METHOD_CALLS`, `FOUR_METHOD_N_SIGNIFICANT`, `FOUR_METHOD_CONSENSUS` and `TUTORIAL` in
  `tests/da/test_consensus.py`.

- [x] **Step 1: Failing test.** In `tests/da/test_consensus.py`, the `r` test reads the shared
  constants and a default-run test checks the tutorial quotes them:

```diff
diff --git a/tests/da/test_consensus.py b/tests/da/test_consensus.py
index 6d53ae0..9dc61c5 100644
--- a/tests/da/test_consensus.py
+++ b/tests/da/test_consensus.py
@@ -1,3 +1,5 @@
+from pathlib import Path
+
 import numpy as np
 import pandas as pd
 import pytest
@@ -168,20 +170,34 @@ def test_numpy_scalars_are_valid_options():
     assert out["consensus"].tolist() == [True]


+# The Phase 3 exit-gate consensus: GlobalPatterns genera, human hosts against the rest, the two native methods and the
+# two R bridges, seeded. The r test measures these counts; the tutorial, whose docs build has no R, quotes them.
+FOUR_METHOD_CALLS = {"ancombc2": 208, "linda": 118, "aldex2": 13, "maaslin3": 52}
+FOUR_METHOD_N_SIGNIFICANT = {0: 413, 1: 114, 2: 62, 3: 35, 4: 12}
+FOUR_METHOD_CONSENSUS = 109
+FOUR_METHOD_CONFLICTS = 0
+TUTORIAL = Path(__file__).parents[2] / "docs" / "tutorials" / "differential_abundance.md"
+
+
 @pytest.mark.r
 def test_four_methods_on_the_exit_gate_data(benchmark):
-    # The Phase 3 exit-gate consensus: GlobalPatterns genera, human hosts against the rest, the two native methods and the
-    # two R bridges. Seeded, so the counts are fixed: calls 208 (ANCOM-BC2), 118 (LinDA), 13 (ALDEx2), 52 (MaAsLin 3).
     results = [
         bt.da.ancombc2(benchmark, "host", reference="other"),
         bt.da.linda(benchmark, "host", reference="other"),
         bt.da.aldex2(benchmark, "host", reference="other", seed=0),
         bt.da.maaslin3(benchmark, "host", reference="other", seed=0),
     ]
-    assert [int((table["qvalue"] < 0.05).sum()) for table in results] == [208, 118, 13, 52]
+    assert {table["method"].iloc[0]: int((table["qvalue"] < 0.05).sum()) for table in results} == FOUR_METHOD_CALLS
     table = bt.da.consensus(results)
     # ANCOM-BC2 and MaAsLin 3 cannot fit the 36 genera absent from one group; LinDA and ALDEx2 test all 636.
     assert table["n_tested"].value_counts().to_dict() == {4: 600, 2: 36}
-    assert table["n_significant"].value_counts().sort_index().to_dict() == {0: 413, 1: 114, 2: 62, 3: 35, 4: 12}
-    assert int(table["consensus"].sum()) == 109 and not table["conflict"].any()
+    assert table["n_significant"].value_counts().sort_index().to_dict() == FOUR_METHOD_N_SIGNIFICANT
+    assert int(table["consensus"].sum()) == FOUR_METHOD_CONSENSUS
+    assert int(table["conflict"].sum()) == FOUR_METHOD_CONFLICTS
     assert (table.loc[table["n_significant"] == 4, "direction"] == 1).all()
+
+
+def test_the_tutorial_quotes_the_four_method_counts():
+    page = TUTORIAL.read_text(encoding="utf-8")
+    for method, calls in FOUR_METHOD_CALLS.items():
+        assert f"| `{method}` | {calls} |" in page
+    assert "| Genera | " + " | ".join(map(str, FOUR_METHOD_N_SIGNIFICANT.values())) + " |" in page
+    sentence = (
+        f"With four methods, {FOUR_METHOD_CONSENSUS} genera are a consensus at `min_methods=2`, "
+        f"{_count_word(FOUR_METHOD_CONFLICTS)} is a conflict, and the {FOUR_METHOD_N_SIGNIFICANT[4]} genera all four "
+        "methods call are all more abundant in human-associated samples"
+    )
+    assert sentence in " ".join(page.split())
+
+
+def _count_word(count: int) -> str:
+    return "none" if count == 0 else str(count)
```

- [x] **Step 2: Run, expect failure.**
  `uv run --group test pytest tests/da/test_consensus.py -q` -> `1 failed, 25 passed, 1
  deselected`: `test_the_tutorial_quotes_the_four_method_counts` raises `FileNotFoundError` (the
  page does not exist). The `r` test is deselected here; Step 6 runs it with R.
- [x] **Step 3: The tutorial.** Create `docs/tutorials/differential_abundance.md`:

````markdown
---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Differential abundance in GlobalPatterns

This tutorial asks which genera differ between human-associated samples and the rest of
GlobalPatterns, runs two differential abundance methods on the same question, and reports where
they agree. It then shows the same report with the two methods that run in R.

:::{note}
The contrast is a demonstration of the workflow, not a biological finding. GlobalPatterns has 9
human-associated samples (feces, skin and tongue) and 17 others: 14 from soil, sediment,
freshwater and the ocean, and 3 mock communities. The data are downloaded from phyloseq's
repository on first use and cached; they are licensed to phyloseq's authors under AGPL-3.
:::

```{code-cell} ipython3
import numpy as np
import pandas as pd

import biotapy as bt
```

## The data

Counts are merged to genus, and genera present in fewer than 20% of the samples are dropped. This
is the one filter: the methods never filter on their own, so each of them is given the same 636
genera.

```{code-cell} ipython3
global_patterns = bt.datasets.global_patterns()
genus = bt.pp.tax_glom(global_patterns, "genus")
tdata = bt.pp.filter_features(genus, min_prevalence=0.2)
human = tdata.obs["SampleType"].isin(["Feces", "Skin", "Tongue"])
tdata.obs["host"] = pd.Categorical(np.where(human, "human", "other"))
tdata.obs["host"].value_counts()
```

## Two methods, chosen first

ANCOM-BC2 and LinDA both run in Python. They are chosen before looking at any result, and run with
the same `group` and `reference`, so their tables compare the same thing (`"human vs other"`): a
positive `effect` means more abundant in human-associated samples.

```{code-cell} ipython3
results = [
    bt.da.ancombc2(tdata, "host", reference="other"),
    bt.da.linda(tdata, "host", reference="other"),
]
{table["method"].iloc[0]: int((table["qvalue"] < 0.05).sum()) for table in results}
```

ANCOM-BC2 leaves zeros out, so it cannot fit the 36 genera with no read in one of the two groups:
their rows are NaN, and they count as not tested, not as "not significant".

```{code-cell} ipython3
int(results[0]["pvalue"].isna().sum())
```

## Where they agree

`bt.da.consensus` counts, per genus, the methods that call it at `qvalue < 0.05` and whether they
agree on its direction. A consensus genus is called by at least `min_methods` methods (2 by
default), all with the same sign.

```{code-cell} ipython3
table = bt.da.consensus(results)
table["n_significant"].value_counts().sort_index()
```

```{code-cell} ipython3
table[["consensus", "conflict"]].sum()
```

Which method calls what, genus by genus:

```{code-cell} ipython3
pd.crosstab(table["significant_ancombc2"], table["significant_linda"])
```

The dot plot shows the genera called by most methods first, then by the size of their effect.
Feature ids in GlobalPatterns are OTU numbers, so the rows are labelled with the genus name; a
name can appear more than once, because GlobalPatterns files some genera under more than one
family and `tax_glom` keeps the lineages apart.

```{code-cell} ipython3
genera = table.rename(index=tdata.var["genus"])
bt.pl.consensus(genera, top=30);
```

The genera only one method calls are where the two disagree. A hollow dot is a genus the other
method tested without calling it:

```{code-cell} ipython3
bt.pl.consensus(genera[genera["n_significant"] == 1], top=30);
```

## Adding the methods that run in R

ALDEx2 and MaAsLin 3 run in R ({ref}`Methods that run in R <da-methods-in-r>`), which this
documentation is built without, so the code below is shown but not run here. To run it, install
the extra with `pip install 'biotapy[r]'` and, in R, `BiocManager::install(c("ALDEx2", "maaslin3"))`:

```python
results = [
    bt.da.ancombc2(tdata, "host", reference="other"),
    bt.da.linda(tdata, "host", reference="other"),
    bt.da.aldex2(tdata, "host", reference="other", seed=0),
    bt.da.maaslin3(tdata, "host", reference="other", seed=0),
]
table = bt.da.consensus(results)
```

biotapy's continuous integration runs exactly this on every pull request, in R 4.5.3 with
Bioconductor 3.22's ALDEx2 and maaslin3, and checks that it gives these numbers:

| Method | Genera called |
|---|---|
| `ancombc2` | 208 |
| `linda` | 118 |
| `aldex2` | 13 |
| `maaslin3` | 52 |

| Methods calling a genus | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| Genera | 413 | 114 | 62 | 35 | 12 |

With four methods, 109 genera are a consensus at `min_methods=2`, none is a conflict, and the 12
genera all four methods call are all more abundant in human-associated samples. The methods differ
a lot in how many genera they call; the consensus table says which calls do not depend on the
choice of method.

## Reading the result

- A genus called by one method only is a result of that method, not of the data alone.
- `min_methods` is part of the analysis: decide it, like the methods, before looking.
- The {doc}`method pages </methods/index>` say what each `effect` measures and how each method's
  numbers compare with its R package.
````

  Link it from the tutorials toctree, the methods index and the guide:

````diff
diff --git a/docs/guide/differential_abundance.md b/docs/guide/differential_abundance.md
index 7cb7117..76fed0e 100644
--- a/docs/guide/differential_abundance.md
+++ b/docs/guide/differential_abundance.md
@@ -2,7 +2,8 @@

 Differential abundance (DA) asks which features are more abundant in one group of samples than in
 another. Every `bt.da` method takes the same arguments and returns the same table, so their answers
-can be put side by side.
+can be put side by side. The {doc}`tutorial </tutorials/differential_abundance>` runs them on
+GlobalPatterns.

 ## One question, one table

diff --git a/docs/methods/index.md b/docs/methods/index.md
index c16667e..93e9853 100644
--- a/docs/methods/index.md
+++ b/docs/methods/index.md
@@ -3,7 +3,8 @@
 One page per differential abundance method: the model it fits, the units of its `effect`, how
 biotapy's call differs from the R function's defaults, how closely the two agree on real data,
 and the paper to cite. The {doc}`differential abundance guide </guide/differential_abundance>`
-says how to run them.
+says how to run them, and the {doc}`tutorial </tutorials/differential_abundance>` runs them on
+GlobalPatterns.

 ```{toctree}
 :maxdepth: 1
diff --git a/docs/tutorials/index.md b/docs/tutorials/index.md
index 018f1e4..bcd81ef 100644
--- a/docs/tutorials/index.md
+++ b/docs/tutorials/index.md
@@ -7,4 +7,5 @@ getting_started
 quick_tour
 phyloseq_analysis
 function
+differential_abundance
 ```
````

  Run `uv run --group test pytest tests/da/test_consensus.py -q` -> `26 passed, 1 deselected`.
- [x] **Step 4: Build and read the page.**
  `rm -rf docs/_build docs/generated && BIOTAPY_DATA_DIR=<scratchpad>/pooch uv run --group doc sphinx-build -W -b html docs docs/_build/html`
  -> `build succeeded.` with `tutorials/differential_abundance.md: Executed notebook in <n>
  seconds` (8.0-9.1 s in the prototype, GlobalPatterns cached). Open
  `docs/_build/html/tutorials/differential_abundance.html` and check the outputs: `host` counts
  `other 17`, `human 9`; `{'ancombc2': 208, 'linda': 118}`; `36`; `n_significant` `0 414`,
  `1 118`, `2 104`; `consensus 104`, `conflict 0`; the crosstab (414 neither, 104 ANCOM-BC2 only,
  14 LinDA only, 104 both); two figures: 30 bold rows of red dots in both columns (Porphyromonas,
  Dialister, Finegoldia, ... at the top), then 30 rows with one filled dot and one hollow dot each.
  The R block is shown as code with no output.
- [x] **Step 5: Bookkeeping.**
  - The checklist line `- [ ] 3.10b \`docs/tutorials/differential_abundance.md\`, the exit-gate notebook`
    becomes `- [x] ...`.
  - Log line, first under the slice 3D heading:
    `- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.10b done: \`docs/tutorials/differential_abundance.md\` runs ANCOM-BC2, LinDA and their consensus on the GlobalPatterns genera on every docs build, and quotes the four-method counts the \`r-bridge\` job measures (a test keeps the two in step).`
- [x] **Step 6: Gate and commit.**

```bash
git add docs/tutorials/differential_abundance.md docs/tutorials/index.md docs/methods/index.md \
  docs/guide/differential_abundance.md tests/da/test_consensus.py .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
git commit -m "docs: add the differential abundance tutorial on GlobalPatterns

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
git status --short                                   # empty
uvx prek run --all-files                             # all Passed
uv run --group test pytest -q -W error::UserWarning  # 1231 passed, 53 deselected
uv run --group test pytest -q -m "golden or network" # 36 passed, 1248 deselected
rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html  # build succeeded.
```

  Where R and the packages are installed (the prototype used the local `biotapy-r-py` image on an
  export of the commit, with the pooch cache mounted):
  `BIOTAPY_DATA_DIR=<pooch> RPY2_CFFI_MODE=API uv run --group test --extra r pytest -m r -q` ->
  `24 passed, 1260 deselected`. Without a local R, the `r-bridge` job on the PR is this check.

### Task 3.12: Coming-from-R check

**Files:**
- Modify: `tests/test_coming_from_r.py`, `docs/_data/r_idioms.toml`, `docs/coming_from_r.md`,
  `docs/tutorials/phyloseq_analysis.md`.

**Not touched:** `docs/extensions/coming_from_r.py` (it already reads every `R equivalent:`
line), any docstring, the other idioms (no `phyloseq_to_deseq2` row: not asked for, R1.4).

**Interfaces:**
- Consumes: `docs/extensions/coming_from_r.py:rows()`, unchanged; the `R equivalent:` lines of
  `pp.relative`, `pp.clr`, `pp.philr` and the four `da` methods; the tutorial page (3.10b) the
  vignette now links.
- Produces: the label `"not in 0.3"`; a test pinning the seven Phase 3 rows.

- [x] **Step 1: Failing tests.** In `tests/test_coming_from_r.py`:

```diff
diff --git a/tests/test_coming_from_r.py b/tests/test_coming_from_r.py
index fc82dc0..af6009c 100644
--- a/tests/test_coming_from_r.py
+++ b/tests/test_coming_from_r.py
@@ -16,7 +16,7 @@ PHYLOSEQ_31 = """otu_table sample_data tax_table phy_tree refseq nsamples ntaxa
 sample_sums taxa_sums rank_names sample_variables get_taxa_unique prune_taxa prune_samples subset_taxa
 subset_samples filter_taxa transform_sample_counts rarefy_even_depth tax_glom estimate_richness distance
 UniFrac ordinate plot_bar plot_richness plot_ordination plot_heatmap import_biom""".split()
-NOT_IN_0_2 = ["tip_glom", "merge_samples", "psmelt", "plot_tree", "plot_net"]
+NOT_IN_0_3 = ["tip_glom", "merge_samples", "psmelt", "plot_tree", "plot_net"]


 def test_the_list_has_31_functions():
@@ -26,12 +26,12 @@ def test_the_list_has_31_functions():
 @pytest.mark.parametrize("name", PHYLOSEQ_31)
 def test_table_maps_each_of_the_31(name):
     cells = coming_from_r.rows()[f"phyloseq::{name}"]
-    assert cells and "not in 0.2" not in cells
+    assert cells and "not in 0.3" not in cells


-@pytest.mark.parametrize("name", NOT_IN_0_2)
+@pytest.mark.parametrize("name", NOT_IN_0_3)
 def test_uncovered_functions_are_marked(name):
-    assert coming_from_r.rows()[f"phyloseq::{name}"] == ["not in 0.2"]
+    assert coming_from_r.rows()[f"phyloseq::{name}"] == ["not in 0.3"]


 @pytest.mark.parametrize(
@@ -41,6 +41,22 @@ def test_mia_importers_map_to_the_function_readers(r_name, function):
     assert coming_from_r.rows()[f"mia::{r_name}"] == [f"{{func}}`bt.io.{function} <biotapy.io.{function}>`"]


+@pytest.mark.parametrize(
+    ("r_name", "functions"),
+    [
+        ("ANCOMBC::ancombc2", ["da.ancombc2"]),
+        ("MicrobiomeStat::linda", ["da.linda"]),
+        ("ALDEx2::aldex", ["da.aldex2"]),
+        ("maaslin3::maaslin3", ["da.maaslin3"]),
+        ("philr::philr", ["pp.philr"]),
+        ("vegan::decostand", ["pp.clr"]),
+        ("mia::transformAssay", ["pp.clr", "pp.philr", "pp.relative"]),
+    ],
+)
+def test_transforms_and_da_methods_map_to_their_r_functions(r_name, functions):
+    assert coming_from_r.rows()[r_name] == [f"{{func}}`bt.{name} <biotapy.{name}>`" for name in functions]
+
+
 def test_plot_functions_link_to_pl():
     assert coming_from_r.rows()["phyloseq::plot_bar"] == ["{func}`bt.pl.bar <biotapy.pl.bar>`"]

```

- [x] **Step 2: Run, expect failure.**
  `uv run --group test pytest tests/test_coming_from_r.py -q` -> `5 failed, 45 passed`: each
  `test_uncovered_functions_are_marked` case fails (`['not in 0.2'] == ['not in 0.3']`). The seven
  new cases pass at once: the rows already come from the 3A-3C docstrings; they are pins.
- [x] **Step 3: Implement.**

````diff
diff --git a/docs/_data/r_idioms.toml b/docs/_data/r_idioms.toml
index fe7a431..1bab3b9 100644
--- a/docs/_data/r_idioms.toml
+++ b/docs/_data/r_idioms.toml
@@ -23,8 +23,8 @@
 "phyloseq::prune_samples" = '`tdata[keep].copy()`'
 "phyloseq::subset_taxa" = '`tdata[:, tdata.var["phylum"] == "Chlamydiae"].copy()`'
 "phyloseq::subset_samples" = '`tdata[tdata.obs["SampleType"] == "Feces"].copy()`'
-"phyloseq::tip_glom" = "not in 0.2"
-"phyloseq::merge_samples" = "not in 0.2"
-"phyloseq::psmelt" = "not in 0.2"
-"phyloseq::plot_tree" = "not in 0.2"
-"phyloseq::plot_net" = "not in 0.2"
+"phyloseq::tip_glom" = "not in 0.3"
+"phyloseq::merge_samples" = "not in 0.3"
+"phyloseq::psmelt" = "not in 0.3"
+"phyloseq::plot_tree" = "not in 0.3"
+"phyloseq::plot_net" = "not in 0.3"
diff --git a/docs/coming_from_r.md b/docs/coming_from_r.md
index 6788ba4..1125af9 100644
--- a/docs/coming_from_r.md
+++ b/docs/coming_from_r.md
@@ -2,7 +2,7 @@

 Every public biotapy function names its R equivalent in its docstring. This table is generated
 from those lines each time the docs are built, plus a short list of phyloseq accessors that are
-plain AnnData/TreeData code (`docs/_data/r_idioms.toml`). Rows marked "not in 0.2" have no
+plain AnnData/TreeData code (`docs/_data/r_idioms.toml`). Rows marked "not in 0.3" have no
 biotapy equivalent yet.

 biotapy keeps samples as rows, so `tdata.X` is phyloseq's `otu_table` with
diff --git a/docs/tutorials/phyloseq_analysis.md b/docs/tutorials/phyloseq_analysis.md
index dba4dcf..6dcad6f 100644
--- a/docs/tutorials/phyloseq_analysis.md
+++ b/docs/tutorials/phyloseq_analysis.md
@@ -15,16 +15,18 @@ kernelspec:

 This notebook redoes the sections of phyloseq's
 [analysis vignette](https://github.com/joey711/phyloseq/blob/master/vignettes/phyloseq-analysis.Rmd)
-that biotapy 0.2 covers, on the same three datasets: GlobalPatterns, enterotype and esophagus.
+that biotapy 0.3 covers, on the same three datasets: GlobalPatterns, enterotype and esophagus.
 Each section names the R chunk it follows. Everything runs in biotapy; nothing is read from R.

-**Not in 0.2**, so left out:
+**Not in 0.3**, so left out:

 - `plot_tree` (exploratory tree plots) and `plot_net` (sample networks);
 - correspondence analysis (`ordinate(..., "CCA")`) and DPCoA, with their scree, species and biplot plots;
 - the `ACE` richness estimator, and `betadiver` distances such as `distance(esophagus, "g")`;
 - `hclust` dendrograms: SciPy's `scipy.cluster.hierarchy.linkage(..., method="average")` does it;
-- multiple testing and differential abundance (biotapy 0.3).
+- the "Multiple Testing and Differential Abundance" section, which points to DESeq2
+  (`phyloseq_to_deseq2`): biotapy does not wrap DESeq2, and runs its own methods in the
+  {doc}`differential abundance tutorial <differential_abundance>`.

 ```{code-cell} ipython3
 import matplotlib.pyplot as plt
@@ -168,7 +170,7 @@ bt.pl.ordination(global_patterns, basis="nmds", color="SampleType");

 `distance(esophagus, "bray")`, `"wunifrac"` and `"jaccard"`. phyloseq's `"jaccard"` is vegan's
 quantitative Jaccard; biotapy's is presence/absence, phyloseq's
-`distance(esophagus, "jaccard", binary = TRUE)`. The `betadiver` method `"g"` is not in 0.2.
+`distance(esophagus, "jaccard", binary = TRUE)`. The `betadiver` method `"g"` is not in 0.3.

 ```{code-cell} ipython3
 esophagus = bt.datasets.esophagus()
````

- [x] **Step 4: Run, expect pass.**
  `uv run --group test pytest tests/test_coming_from_r.py -q` -> `50 passed`.
- [x] **Step 5: Check the rendered page.** Build the docs, then read
  `docs/generated/coming_from_r_table.md`: 58 rows (`grep -c '^| \`' ...`); the seven Phase 3 rows
  link `bt.pp.*` and `bt.da.*`; five rows read "not in 0.3"; no "not in 0.2" anywhere
  (`grep -rn "not in 0.2" docs --include=*.md --include=*.toml` prints nothing outside
  `docs/_build`). The vignette's first bullet list links the tutorial.
- [x] **Step 6: Bookkeeping.**
  - The checklist line `- [ ] 3.12 Coming-from-R check` becomes `- [x] 3.12 Coming-from-R check`.
  - Log line: `- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.12 done: the Coming-from-R rows for \`pp.clr\`, \`pp.philr\` and the four \`da\` methods are pinned by a test; phyloseq calls without an equivalent read "not in 0.3", and the phyloseq vignette points its differential abundance section to the new tutorial.`
- [x] **Step 7: Gate and commit.**

```bash
git add docs/_data/r_idioms.toml docs/coming_from_r.md docs/tutorials/phyloseq_analysis.md tests/test_coming_from_r.py \
  .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
git commit -m "docs: pin the transform and DA rows of Coming from R, mark uncovered calls not in 0.3

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
git status --short                                   # empty
uvx prek run --all-files                             # all Passed
uv run --group test pytest -q -W error::UserWarning  # 1238 passed, 53 deselected
uv run --group test pytest -q -m "golden or network" # 36 passed, 1255 deselected
rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html  # build succeeded.
```

### Task 3.13: Benchmarks

**Files:**
- Create: `benchmarks/benchmarks/da.py`.
- Modify (the `docs(da)` commit): `src/biotapy/da/_ancombc.py` (docstring), `.knowledge/modules/da.md`.
- Modify: `benchmarks/benchmarks/_data.py`, `benchmarks/benchmarks/pp.py`, `docs/performance.md`.

**Not touched:** `Preprocessing`, `tl.py`, `fn.py`, `synthetic()`, `asv.conf.json`, the workflows
(the `lint` job's `asv check` already imports every module), any `src/` file but `da.ancombc2`'s
docstring. No optimisation (R10.1). The ANCOM-BC2 memory figure is reworded in its own
`docs(da)` commit (decision 8, approved 2026-10-07: a fix, not a report).

**Interfaces:**
- Consumes: `bt.pp.philr(tdata, *, pseudocount=0.5)`, `bt.da.linda`, `bt.da.ancombc2`;
  `_data.py:synthetic(n_obs, n_vars)`.
- Produces: `_data.py:N_PHILR_OBS`, `N_DA_OBS`, `N_DA_VARS`, `DA_DENSITY`,
  `synthetic_genera() -> AnnData`; the asv classes `pp.Philr` and `da.DifferentialAbundance`; the
  "Transforms and differential abundance" table in `docs/performance.md`.

- [x] **Step 1: The suite.**

```diff
diff --git a/benchmarks/benchmarks/_data.py b/benchmarks/benchmarks/_data.py
index f16e962..8bd5f49 100644
--- a/benchmarks/benchmarks/_data.py
+++ b/benchmarks/benchmarks/_data.py
@@ -1,4 +1,4 @@
-"""Synthetic benchmark data: a sparse count table with taxonomy and a tree, a function table, and taxon traits."""
+"""Synthetic benchmark data: a count table with taxonomy and a tree, a genus table, a function table, taxon traits."""

 import numpy as np
 import pandas as pd
@@ -16,6 +16,10 @@ _RANK_SIZES = {"phylum": 10_000, "class": 5_000, "order": 2_500, "family": 500,
 N_SAMPLES, N_FUNCTIONS, N_STRATA, N_GROUPS = 1_600, 500, 43, 50
 # Tian et al.'s functional redundancy at 2,000 taxa, as measured when fn.functional_redundancy was written.
 N_TAXA, N_GENES, N_ABUNDANCE_SAMPLES = 2_000, 2_500, 100
+# PhILR densifies X: 1,000 samples keep it at 400 MB on the 50,000-feature table.
+N_PHILR_OBS = 1_000
+# Differential abundance on a large genus-level study: 2,000 samples in two groups, 10,000 features, 30% non-zero.
+N_DA_OBS, N_DA_VARS, DA_DENSITY = 2_000, 10_000, 0.3


 def synthetic(n_obs: int = N_OBS, n_vars: int = N_VARS) -> TreeData:
@@ -72,3 +76,14 @@ def synthetic_traits() -> tuple[AnnData, pd.DataFrame]:
     )
     copies = rng.integers(1, 6, size=(N_TAXA, N_GENES)) * (rng.random((N_TAXA, N_GENES)) >= 0.7)
     return adata, pd.DataFrame(copies, index=taxa, columns=[f"g{k}" for k in range(N_GENES)])
+
+
+def synthetic_genera() -> AnnData:
+    """2,000 samples x 10,000 features, 30% non-zero counts from 1 to 99, ``obs["group"]`` a / b, seed 0."""
+    rng = np.random.default_rng(SEED)
+    X = sp.random(
+        N_DA_OBS, N_DA_VARS, density=DA_DENSITY, format="csr", rng=rng, data_rvs=lambda size: rng.integers(1, 100, size)
+    )
+    group = pd.Categorical(np.repeat(["a", "b"], N_DA_OBS // 2))
+    obs = pd.DataFrame({"group": group}, index=[f"s{i}" for i in range(N_DA_OBS)])
+    return AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"f{j}" for j in range(N_DA_VARS)]))
diff --git a/benchmarks/benchmarks/pp.py b/benchmarks/benchmarks/pp.py
index 4966725..a0aa8b6 100644
--- a/benchmarks/benchmarks/pp.py
+++ b/benchmarks/benchmarks/pp.py
@@ -1,9 +1,9 @@
-"""Preprocessing at 5,000 samples x 50,000 features."""
+"""Preprocessing at 5,000 samples x 50,000 features, and PhILR at 1,000 samples."""

 import biotapy as bt
 from biotapy._core import TreeData

-from ._data import synthetic
+from ._data import N_PHILR_OBS, synthetic


 class Preprocessing:
@@ -26,3 +26,21 @@ class Preprocessing:
     def time_rarefy(self, tdata: TreeData) -> None:
         """``pp.rarefy`` to the smallest sample depth."""
         bt.pp.rarefy(tdata, seed=0)
+
+
+class Philr:
+    """``pp.philr`` at 1,000 samples x 50,000 features: ``X`` densified once (400 MB), then 49,999 balances."""
+
+    number, repeat, rounds, timeout = 1, 3, 1, 600
+
+    def setup_cache(self) -> TreeData:
+        """Build the table once; asv pickles it for every benchmark."""
+        return synthetic(n_obs=N_PHILR_OBS)
+
+    def time_philr(self, tdata: TreeData) -> None:
+        """Centred log-ratios times the tree's sparse orthonormal basis."""
+        bt.pp.philr(tdata)
+
+    def peakmem_philr(self, tdata: TreeData) -> None:
+        """Peak memory: the dense table, its log-ratios and the balances."""
+        bt.pp.philr(tdata)
```

  Create `benchmarks/benchmarks/da.py`:

```python
"""Differential abundance at 2,000 samples x 10,000 features: the two methods that run without R."""

from anndata import AnnData

import biotapy as bt

from ._data import synthetic_genera


class DifferentialAbundance:
    """``da.linda`` and ``da.ancombc2``, two groups of 1,000 samples; each densifies ``X`` once (160 MB)."""

    number, repeat, rounds, timeout = 1, 3, 1, 600

    def setup_cache(self) -> AnnData:
        """Build the table once; asv pickles it for every benchmark."""
        return synthetic_genera()

    def time_linda(self, adata: AnnData) -> None:
        """One least-squares fit per feature, then the mode of the coefficients."""
        bt.da.linda(adata, "group")

    def peakmem_linda(self, adata: AnnData) -> None:
        """Peak memory; the docstring says about five dense copies of ``X``."""
        bt.da.linda(adata, "group")

    def time_ancombc2(self, adata: AnnData) -> None:
        """scikit-bio's ANCOM-BC2, its bias E-M capped at 100 iterations."""
        bt.da.ancombc2(adata, "group")

    def peakmem_ancombc2(self, adata: AnnData) -> None:
        """Peak memory; ``da.ancombc2``'s docstring gives 4.5x to 7.5x the dense ``X``, 4.5x on this table."""
        bt.da.ancombc2(adata, "group")
```

  Why the suite looks as it does: PhILR at 1,000 samples keeps the dense table at 400 MB (the
  full 5,000 would be 2 GB and peak near 10 GB); `synthetic(n_obs=...)` reuses Phase 1's
  construction, whose odd tree levels leave one-child nodes that `pp.philr` skips, as on a filtered
  real tree. The DA table is 30% non-zero because zeros drive ANCOM-BC2's missing-value path; 2,000
  x 10,000 is a large genus-level study. `AnnData` with counts inferred from the values (no `uns`),
  as `synthetic_traits` builds its table.
- [x] **Step 2: Check it imports.**
  `cd benchmarks && uv run --group dev asv check --python=same` -> `No problems found.`; then
  `uvx ruff check benchmarks` -> `All checks passed!` and `uvx ruff format --check benchmarks` ->
  `6 files already formatted`.
- [x] **Step 3: Gate and commit the suite.**

```bash
git add benchmarks/benchmarks/_data.py benchmarks/benchmarks/pp.py benchmarks/benchmarks/da.py
git commit -m "perf: add asv benchmarks for PhILR, LinDA and ANCOM-BC2

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
git status --short                                   # empty
uvx prek run --all-files                             # all Passed
uv run --group test pytest -q -W error::UserWarning  # 1238 passed, 53 deselected
```

- [x] **Step 3b: The docstring fix** (decision 8), its own commit after the baseline run, so the
  baseline is measured on the suite commit. In `src/biotapy/da/_ancombc.py`, `da.ancombc2`'s Notes:

```diff
-    samples x features. Peak memory is about seven such arrays (7.2x to 7.5x
-    measured on 400 x 1,000 and 200 x 2,000 tables).
+    samples x features. Peak memory is 4.5x to 7.5x such an array, depending on
+    the table's shape: 4.5x measured on a 2,000 x 10,000 table (30% non-zero),
+    7.2x to 7.5x on 400 x 1,000 and 200 x 2,000 tables.
```

  `.knowledge/modules/da.md` ("ANCOM-BC2 about 7x" becomes "4.5x to 7.5x, depending on the table's
  shape"; `generated` and `commit` bumped) and a `.knowledge/log.md` line go in the same commit:
  `docs(da): state ANCOM-BC2's measured peak-memory range`.
- [x] **Step 4: Run the baseline** on that commit, with `git status --short` empty and other heavy
  work closed (the numbers are the record). asv's `branches: ["master"]` needs a local `master`
  branch (`git branch --list master` must print it):

```bash
cat /proc/loadavg
cd benchmarks
uv run --group dev env HOME="$PWD/../.asv" asv machine --yes
uv run --group dev env HOME="$PWD/../.asv" asv run --python=same --bench "^(da\.|pp\.Philr)" --set-commit-hash "$(git rev-parse HEAD)" --show-stderr
uv run --group dev env HOME="$PWD/../.asv" asv show "$(git rev-parse HEAD)" --bench "^(da\.|pp\.Philr)"
cat /proc/loadavg
```

  Done 2026-10-07 on `51e0475` (load average 1.66, 1.30, 1.21 before, 4.64, 2.23, 1.54 after;
  55 s): `time_philr` 3.12 s, `peakmem_philr` 2.14 GB, `time_linda` 425 ms, `peakmem_linda` 1.15 GB,
  `time_ancombc2` 3.62 s, `peakmem_ancombc2` 1.07 GB. The calls' own peaks on this run:
  `pp.philr` 1.72 GB (4.31x its dense table), `da.linda` 822 MB (5.14x), `da.ancombc2` 726 MB
  (4.54x). None is more than twice the prototype's.

  Expect six results, each `ok`, in about 60 s. The prototype's (`7b15403`, the laptop of
  Phases 1 and 2: 11th Gen Intel Core i7-11800H, 16 threads, 62 GB, Linux, Python 3.13):
  `time_philr` 3.04 s, `peakmem_philr` 2.14 GB, `time_linda` 469 ms, `peakmem_linda` 1.15 GB,
  `time_ancombc2` 3.92 s, `peakmem_ancombc2` 1.07 GB. A result more than twice these, or a
  failure, is reported before it is recorded (R14.1).
  The calls' own peaks (Step 5's second bullet), measured outside asv with
  `/proc/self/statm` before the call and `resource.getrusage(...).ru_maxrss` after it, in a fresh
  process per function on the same tables: `pp.philr` 1.72 GB, `da.linda` 821 MB, `da.ancombc2`
  723 MB. Replace them with the run's values.
- [x] **Step 5: The performance page.** Add the section before "## Running the benchmarks",
  replacing the commit, date, load averages and every number with Step 4's run:

````diff
diff --git a/docs/performance.md b/docs/performance.md
index bd48d26..df4e2a0 100644
--- a/docs/performance.md
+++ b/docs/performance.md
@@ -65,6 +65,37 @@ The hierarchy puts every function in two of 50 groups. Functional redundancy run
 - `functional_redundancy` compares every pair of taxa, so its time grows with the square of the
   taxa times the genes, and its memory with the square of the taxa.

+## Transforms and differential abundance
+
+Measured on commit `51e0475`, 2026-10-07, on the same laptop and environment, with
+`asv run --python=same --bench "^(da\.|pp\.Philr)"` (the two classes took about 55 s). The load
+average was 1.66, 1.30, 1.21 before the run and 4.64, 2.23, 1.54 after it. Release 0.3 sets no
+speed target either.
+
+PhILR runs on the synthetic table above built with 1,000 samples (50,000 features, 2% non-zero,
+the same balanced tree). The two differential abundance methods that run without R run on 2,000
+samples in two groups of 1,000 and 10,000 features, 30% non-zero counts from 1 to 99, the size of
+a large genus-level study. All of it comes from seed 0 (`benchmarks/benchmarks/_data.py`). ALDEx2
+and MaAsLin 3 are not benchmarked: their time is R's.
+
+| Benchmark | Result |
+|---|---|
+| `pp.philr`, 1,000 x 50,000 | 3.12 s |
+| `pp.philr`, peak memory | 2.14 GB |
+| `da.linda`, 2,000 x 10,000 | 425 ms |
+| `da.linda`, peak memory | 1.15 GB |
+| `da.ancombc2`, 2,000 x 10,000 | 3.62 s |
+| `da.ancombc2`, peak memory | 1.07 GB |
+
+- Each of the three densifies `X` once (rules.md R6.2): 400 MB for PhILR, 160 MB for the
+  differential abundance table. The peaks also hold the Python process (418 MB with PhILR's table
+  loaded, 345 MB with the differential abundance table).
+- Measured outside asv on the same tables, the call itself added 1.72 GB for `pp.philr` (4.3 times
+  the dense table, plus scikit-bio's sparse basis, as its docstring says), 822 MB for `da.linda`
+  (5.1 times, as its docstring says) and 726 MB for `da.ancombc2` (4.5 times; its docstring gives
+  4.5 times here and 7.2 to 7.5 times on the 400 x 1,000 and 200 x 2,000 tables it was first
+  measured on, so the range depends on the table's shape). Nothing is being optimised
+  (rules.md R10.1).
+
 ## Running the benchmarks

 ```bash
@@ -77,6 +108,6 @@ uv run --group dev env HOME="$PWD/../.asv" asv show "$(git rev-parse HEAD)"
 ```

 `--python=same` runs in the current environment, and `--set-commit-hash` keeps the results, in
-`.asv/results`. Add `--bench "^fn\."` to `asv run` and `asv show` to run or show only the function
-benchmarks. The whole suite took about 13 minutes on this run and needs about 3 GB of free memory.
+`.asv/results`. Add `--bench "^fn\."` (or `"^(da\.|pp\.Philr)"`) to `asv run` and `asv show` to
+run or show only the function (or transform and differential abundance) benchmarks. The whole suite
+took about 13 minutes on this run and needs about 3 GB of free memory.
 `uv run --group dev asv check --python=same` imports the suite without running it; CI runs it.
````

- [x] **Step 6: Bookkeeping.**
  - The checklist line `- [ ] 3.13 asv benchmarks for \`pp.philr\`, \`da.linda\`, \`da.ancombc2\``
    becomes `- [x] ...`.
  - Log line: `- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.13 done: asv baselines for \`pp.philr\`, \`da.linda\` and \`da.ancombc2\` in docs/performance.md.`
- [x] **Step 7: Gate and commit the page.**

```bash
git add docs/performance.md .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
git commit -m "docs: record the asv baselines for PhILR, LinDA and ANCOM-BC2

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
git status --short                                   # empty
uvx prek run --all-files                             # all Passed
uv run --group test pytest -q -W error::UserWarning  # 1238 passed, 53 deselected
uv run --group test pytest -q -m "golden or network" # 36 passed, 1255 deselected
rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html  # build succeeded.
```

### Checkpoint D - review slice 3D

- [ ] **Review the whole slice** with superpowers:requesting-code-review (opus: the pages are
  claims about statistics a user will act on), against:
  - data-model-slots ("DA results"), function-shape, r-golden-parity, module-boundaries;
  - pure-by-default and da-consensus-agreement;
  - the slice 3D design, review focus and global constraints above.

  The reviewer reads each "Compared with R" table against the R `formals()` recorded in "How
  slice 3D was checked", the call in the docstring and `tests/r/export_golden.R` (review focus
  2), and each "Agreement with R" number against its golden test. Reviewers may download
  GlobalPatterns (always with an absolute `BIOTAPY_DATA_DIR`) and read the rendered pages. Then a
  fix pass, one commit per finding, each with a test where a test can show it (a docs-only finding
  is shown by the docs build or a `tests/test_ci.py`-style check); then a scoped re-review of the
  fix pass. Record: review counts (Critical / Important / Minor), fix pass range, re-review result.
- [ ] **Task 3.14** (below), after the fix pass.
- [ ] **Run the gates** on the committed tree (`git status --short` empty first), every command
  with an absolute `BIOTAPY_DATA_DIR`:
  - `uvx prek run --all-files`;
  - `uv run --group test pytest -q -W error::UserWarning` (1238 passed, 53 deselected before the
    fix pass);
  - `uv run --group test pytest -q -m "golden or network"` (36 passed);
  - `rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html`;
  - `cd benchmarks && uv run --group dev asv check --python=same` (`No problems found.`);
  - `bash scripts/knowledge_stale.sh --against HEAD` prints `25 current, 0 stale, 12 uncheckable`.

  Confirm `~/.cache/biotapy` and `docs/pooch/` do not exist.
- [ ] **Push** `phase-3d`, open the PR, and merge-commit on green (approved 2026-10-05 for Phase 3
  slice branches). CI must be green: every hatch-test job, `lint` (with `asv check`, which imports
  `da.py`), `import-without-extras`, `network`, `r-bridge` (its log shows `24 passed`: the `r`
  test now reads the constants the tutorial test checks), `docs` (its log shows five `Executed
  notebook` lines, one of them `tutorials/differential_abundance.md`). Record the PR number, the
  merge commit, the Test run id and the docs job id: they are exit-gate items 1 and 3's evidence,
  ticked in Task 3.15 Step 8. If `r-bridge` fails only on the four-method numbers, stop and report
  the measured differences; do not edit the constants to match without the user (R11.5).
- [ ] **After the merge**, once Read the Docs has built `master`:
  `curl -s -o /dev/null -w "%{http_code}\n" https://biotapy.readthedocs.io/en/latest/tutorials/differential_abundance.html`
  and the same for `methods/linda.html` print `200`, and the tutorial shows its two figures. Read
  the Docs keeps no pooch cache, so its build downloads GlobalPatterns (already the case for the
  phyloseq vignette). If a page is missing, read the build log on readthedocs.org and report it;
  do not work around it.
- [ ] **Ask the user, in one message:**
  1. to review slice 3D before the release;
  2. to approve pushing `release-0.3.0`, opening its PR and merge-committing it on green (Task
     3.15 Steps 3-9; R13.3: the standing approval names slice branches only);
  3. to confirm that on pypi.org project `biotapy` still lists the trusted publisher (owner
     `pedrocr83`, repository `biotapy`, workflow `release.yaml`, environment `pypi`) that published
     0.2.0. (Checked 2026-10-07 from GitHub's side: environment `pypi` exists with no protection
     rules.)

  Stop until the user answers.

### Task 3.14: Knowledge

Run as Checkpoint D's knowledge step, after the fix pass, so it documents the final code. Use the
codebase-map templates (R12.2-R12.4). Most of `modules/da.md` was written at Checkpoints B and C;
this task adds only what slice 3D changed.

**Files:**
- Modify: `.knowledge/modules/da.md`, `.knowledge/playbooks/add-a-function.md`,
  `.knowledge/roadmap/phase-3-stats.md`, `.knowledge/log.md`, and `commit`/`generated` only on
  the concepts Step 2 lists.

**Not touched:** `modules/pp.md`, `modules/pl.md`, `modules/core.md`,
`decisions/optional-heavy-dependencies.md` (re-read: nothing they state changed; `pp.md`'s "philr
peaks at about five arrays" agrees with the measured 4.3x plus the basis), any `src/` file, every
contract's body.

**Interfaces:**
- Consumes: the final slice 3D tree; `bash scripts/knowledge_stale.sh --touched --against master`.
- Produces: a bundle with 0 stale concepts against `HEAD`.

- [ ] **Step 1: Concepts the slice changes.** Re-check each sentence against the code after the
  fix pass (a fix may change a number). `modules/da.md` already has the ANCOM-BC2 memory gotcha
  (`b9cb006`: "4.5x to 7.5x, depending on the table's shape") and the MaAsLin 3 and ALDEx2
  bounds as measured; what is left is the frontmatter, the Verification paragraph and the
  tutorial gotcha:

````diff
diff --git a/.knowledge/modules/da.md b/.knowledge/modules/da.md
index 49b19ac..0fbbac9 100644
--- a/.knowledge/modules/da.md
+++ b/.knowledge/modules/da.md
@@ -5,8 +5,8 @@ description: Four differential abundance methods, native LinDA and ANCOM-BC2 and
 resource: /src/biotapy/da/
 paths: ["src/biotapy/da/**"]
 tags: [da, differential-abundance, linda, ancombc2, aldex2, maaslin3, rpy2]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-07T15:00:00Z }
-commit: df27399
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 status: stable
 ---

@@ -142,6 +142,9 @@ uv run --group test pytest -m "golden or network" tests/da -q
 uv run --group test --extra r pytest -m r tests/da -q
 ```

+The docs build also runs `da.linda`, `da.ancombc2` and `da.consensus` in
+`docs/tutorials/differential_abundance.md`; each method's model, R settings
+and measured agreement with R are written up in `docs/methods/`.
 The second needs the pooch cache of GlobalPatterns (`BIOTAPY_DATA_DIR`); at
 the 64fe39d check it gave `4 passed`. The third needs R with ALDEx2 and
 maaslin3 installed (the `r-bridge` CI job, or the golden image in
@@ -232,6 +235,13 @@ per-module command exists: unknown.
   table at peak; ANCOM-BC2 4.5x to 7.5x, depending on the table's shape). rules.md R6.2 allows it: a native
   method whose algorithm needs the full table may densify once (the user
   approved that wording after Checkpoint B, 2026-10-05).
+- **The tutorial quotes the four-method counts.** The docs build has no R, so
+  `docs/tutorials/differential_abundance.md` prints the counts the `r` test
+  `tests/da/test_consensus.py::test_four_methods_on_the_exit_gate_data`
+  measures; both read the constants `FOUR_METHOD_*` in that file, and
+  `test_the_tutorial_quotes_the_four_method_counts` fails until the tutorial
+  matches them. A change to any method's numbers on GlobalPatterns updates
+  the constants and the tutorial together.
 - Replicate rows within a group give `se = 0` in LinDA and p-values that are
   floating-point noise; R does the same.


diff --git a/.knowledge/playbooks/add-a-function.md b/.knowledge/playbooks/add-a-function.md
index 9425ed2..6350bad 100644
--- a/.knowledge/playbooks/add-a-function.md
+++ b/.knowledge/playbooks/add-a-function.md
@@ -5,8 +5,8 @@ description: The only sanctioned path from "we need X" to a merged public functi
 tags: [workflow, api, testing]
 status: stable
 paths: ["src/biotapy/**", "tests/**", "docs/**"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T14:12:26Z }
-commit: 64fe39d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 ---

 # When
@@ -39,7 +39,9 @@ Any new entry in a subpackage `__all__`.
    feeds the generated Coming-from-R page; `tests/test_docstrings.py` fails
    if it cannot be parsed.
 8. **Docs page**: add or extend the `docs/guide/<concept>.md` page the
-   docstring links to; add the function to `docs/api.md`.
+   docstring links to; add the function to `docs/api.md`. A `da` method also
+   gets a page in `docs/methods/` (model, units, the R defaults it changes,
+   agreement with R, reference) and a row in the guide's method table.
 9. **Knowledge**: update `.knowledge/` only if a contract, decision or module
    concept is now wrong or incomplete - see [maintain-knowledge](/playbooks/maintain-knowledge.md).
 10. **Tick** the task checkbox in the active phase concept.
````

- [ ] **Step 2: Concepts the diff only touches.** Run
  `bash scripts/knowledge_stale.sh --against HEAD` before Step 1's edits; on the fix pass's final
  commit it reports `16 current, 9 stale, 12 uncheckable`, the nine being `phase-0-foundation`,
  `phase-1-core`, `da`, `data-model-slots`, `engine-parity`, `function-shape`, `module-boundaries`,
  `r-golden-parity` and `add-a-function`. (`phase-3-stats` is current because the fix pass bumped
  it; Steps 3 and 4 edit it again, so bump it here too.) `da` and `add-a-function` are edited by
  Step 1; read each of the other seven against the slice diff, and nothing they state became false.
  Engine-parity's 5,000 x 50,000 table still gates compiled engines, and r-golden-parity's rows are
  unchanged (the `r` test only reads constants now). The three contracts that cover
  `src/biotapy/**` are stale because `3b6239f`, `b9cb006` and the fix pass edit docstring prose in
  `src/biotapy/da/` (`_aldex2.py`, `_ancombc.py`, `_maaslin3.py`) and nothing else: no signature,
  slot, result key, layer or import changed, so data-model-slots, function-shape and
  module-boundaries state the same things. Bump only `generated` and `commit` on
  each, and on `modules/da.md`, `add-a-function` and `phase-3-stats`:

```diff
diff --git a/.knowledge/contracts/engine-parity.md b/.knowledge/contracts/engine-parity.md
index a7d96ec..819d2d2 100644
--- a/.knowledge/contracts/engine-parity.md
+++ b/.knowledge/contracts/engine-parity.md
@@ -5,8 +5,8 @@ description: A compiled kernel is a drop-in behind an existing public function v
 tags: [performance, testing]
 status: stable
 paths: ["src/biotapy/**/*.py", "rust/**", "benchmarks/**"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T14:12:26Z }
-commit: 64fe39d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 sources:
   - id: spec
     resource: ../../plan.md
diff --git a/.knowledge/contracts/r-golden-parity.md b/.knowledge/contracts/r-golden-parity.md
index 17076d8..eb7afa6 100644
--- a/.knowledge/contracts/r-golden-parity.md
+++ b/.knowledge/contracts/r-golden-parity.md
@@ -5,8 +5,8 @@ description: Every computation with an R equivalent is tested against gzip CSV g
 tags: [testing, r, validation]
 status: stable
 paths: ["tests/r/**", "tests/golden/**", "tests/**/test_*.py"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T13:50:53Z }
-commit: 8a1f184
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 sources:
   - id: spec
     resource: ../../plan.md
diff --git a/.knowledge/contracts/data-model-slots.md b/.knowledge/contracts/data-model-slots.md
index 0c6fdb1..20817b8 100644
--- a/.knowledge/contracts/data-model-slots.md
+++ b/.knowledge/contracts/data-model-slots.md
@@ -5,8 +5,8 @@ description: Which AnnData/TreeData slot holds what, the exact result keys, the
 tags: [data-model, api]
 status: stable
 paths: ["src/biotapy/_core/**", "src/biotapy/io/**", "src/biotapy/pp/**", "src/biotapy/tl/**", "src/biotapy/fn/**", "src/biotapy/da/**"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T14:12:26Z }
-commit: 64fe39d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 sources:
   - id: spec
     resource: ../../plan.md
diff --git a/.knowledge/contracts/function-shape.md b/.knowledge/contracts/function-shape.md
index 23963bd..df5b987 100644
--- a/.knowledge/contracts/function-shape.md
+++ b/.knowledge/contracts/function-shape.md
@@ -5,8 +5,8 @@ description: One task = one public function `verb(data, required, *, options) ->
 tags: [api, conventions, docs]
 status: stable
 paths: ["src/biotapy/**/*.py"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T14:12:26Z }
-commit: 64fe39d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 sources:
   - id: spec
     resource: ../../plan.md
diff --git a/.knowledge/contracts/module-boundaries.md b/.knowledge/contracts/module-boundaries.md
index b1067bc..d9cdb45 100644
--- a/.knowledge/contracts/module-boundaries.md
+++ b/.knowledge/contracts/module-boundaries.md
@@ -5,8 +5,8 @@ description: Layered package (_core at the bottom, pl/ml/da at the top); public
 tags: [architecture, modularization]
 status: stable
 paths: ["src/biotapy/**", "pyproject.toml"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T14:12:26Z }
-commit: 64fe39d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 sources:
   - id: spec
     resource: ../../plan.md
diff --git a/.knowledge/roadmap/phase-0-foundation.md b/.knowledge/roadmap/phase-0-foundation.md
index de08f7c..1555cd0 100644
--- a/.knowledge/roadmap/phase-0-foundation.md
+++ b/.knowledge/roadmap/phase-0-foundation.md
@@ -9,8 +9,8 @@ phase_state: done
 effort: ~1 week part-time
 depends_on: []
 paths: ["pyproject.toml", ".pre-commit-config.yaml", ".github/**", "docs/**", "src/biotapy/__init__.py", "src/biotapy/_core/**", "tests/**", "scripts/**"]
-generated: { by: claude-code/claude-opus-5-5, at: 2026-10-06T14:27:16Z }
-commit: bade2ba
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 sources:
   - id: spec
     resource: ../../plan.md
diff --git a/.knowledge/roadmap/phase-1-core.md b/.knowledge/roadmap/phase-1-core.md
index fe70fd5..0f6a419 100644
--- a/.knowledge/roadmap/phase-1-core.md
+++ b/.knowledge/roadmap/phase-1-core.md
@@ -9,8 +9,8 @@ phase_state: done
 effort: 6-8 weeks part-time (spec); slices 1A-1D with checkpoints
 depends_on: [/roadmap/phase-0-foundation.md]
 paths: ["src/biotapy/**", "tests/**", "docs/**", "benchmarks/**"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T14:12:26Z }
-commit: 64fe39d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
+commit: d6d82ef
 sources:
   - id: spec
     resource: ../../plan.md
```

  Re-run `bash scripts/knowledge_stale.sh --against HEAD` after the commit: `25 current, 0 stale,
  12 uncheckable`.
- [ ] **Step 3: Roadmap.** The checklist line `- [ ] 3.14 Knowledge (...)` becomes `- [x] ...`;
  record Checkpoint D's review counts, fix range and re-review in its first box.
- [ ] **Step 4: Log**, first under the slice 3D heading:

```markdown
- **Update**: [da](modules/da.md): Verification names the tutorial and the method pages; new gotcha: the tutorial quotes the `r` test's four-method counts through shared constants and a test.
- **Update**: [add-a-function](playbooks/add-a-function.md) step 8: a `da` method also gets a page in `docs/methods/` and a row in the guide's method table.
- **Verification**: re-checked against `29efd9c..<HEAD>` and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [phase-1-core](roadmap/phase-1-core.md), [engine-parity](contracts/engine-parity.md), [r-golden-parity](contracts/r-golden-parity.md), [data-model-slots](contracts/data-model-slots.md), [function-shape](contracts/function-shape.md), [module-boundaries](contracts/module-boundaries.md).
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) task 3.14 done.
```

- [ ] **Step 5: Gate and commit.**

```bash
git add .knowledge/modules/da.md .knowledge/playbooks/add-a-function.md .knowledge/roadmap/phase-0-foundation.md \
  .knowledge/roadmap/phase-1-core.md .knowledge/contracts/engine-parity.md .knowledge/contracts/r-golden-parity.md \
  .knowledge/contracts/data-model-slots.md .knowledge/contracts/function-shape.md .knowledge/contracts/module-boundaries.md \
  .knowledge/roadmap/phase-3-stats.md .knowledge/log.md
git commit -m "docs(knowledge): document the DA docs and benchmarks and refresh concepts for slice 3D

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
git status --short                                   # empty
uvx prek run --all-files                             # all Passed
uv run --group test pytest -q -W error::UserWarning  # 1238 passed, 53 deselected (plus the fix pass's tests)
bash scripts/knowledge_stale.sh --against HEAD       # 25 current, 0 stale, 12 uncheckable
```

### Task 3.15: Release 0.3.0

Follows [cut-a-release](/playbooks/cut-a-release.md) and Task 2.15. Publishing to PyPI is
irreversible: a version can never be uploaded twice. Every outward step waits for the user's
explicit approval of that step (R13.3).

**Files:**
- Modify (release commit): `pyproject.toml` (version), `CHANGELOG.md`, `README.md`.
- Modify (knowledge commit): `.knowledge/playbooks/cut-a-release.md`, the concepts the version
  bump makes stale, `.knowledge/roadmap/phase-3-stats.md`, `.knowledge/log.md`.
- Modify (after the upload, branch `close-phase-3`): `.knowledge/roadmap/phase-3-stats.md`,
  `.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/roadmap/index.md`, `.knowledge/log.md`.

**Interfaces:**
- Consumes: `master` after Checkpoint D's merge; `.github/workflows/release.yaml` (on
  `release: published`: `uv build`, then `pypa/gh-action-pypi-publish` v1.14.2 in environment
  `pypi`).
- Produces: `biotapy 0.3.0` on PyPI with the extra `r`; Phase 3 `phase_state: done`; Phase 4
  active.

**Not touched:** any `src/`, `tests/` or `docs/` file (the version reaches the docs through
package metadata; `docs/changelog.md` includes `CHANGELOG.md`), the workflows. `uv.lock` is
git-ignored, so the bump leaves no lock diff.

- [ ] **Step 1: Approval in hand.** Checkpoint D's last box asked for the slice review, the
  release branch push/PR/merge and the PyPI publisher check. Do not start until the user has
  answered all three.
- [ ] **Step 2: Check the publish action.**
  `gh api repos/pypa/gh-action-pypi-publish/releases/latest --jq .tag_name` prints `v1.14.2`
  (checked 2026-10-07; `release.yaml` pins `v1.14.2`). If it prints anything newer, stop and
  report it: a stale pin failed 0.0.1 (the playbook's Common mistakes).
- [ ] **Step 3: Branch.**
  `git switch master && git pull --ff-only && git switch -c release-0.3.0`.
- [ ] **Step 4: Version and changelog.** `## [Unreleased]` is empty (no slice 3A-3D PR wrote an
  entry). Write the 0.3.0 entries from `git log v0.2.0..master --no-merges --oneline` and the diff
  below, with the date from `date -u +%F`; leave a new empty `## [Unreleased]` above:

```diff
diff --git a/CHANGELOG.md b/CHANGELOG.md
index 4f09d94..63c6bc2 100644
--- a/CHANGELOG.md
+++ b/CHANGELOG.md
@@ -10,6 +10,44 @@ and this project adheres to [Semantic Versioning][].

 ## [Unreleased]

+## [0.3.0] - 2026-10-07
+
+### Added
+
+- `bt.pp.clr`: centred log-ratio transform into `layers["clr"]`, equal to
+  `vegan::decostand(x, "clr", pseudocount = 0.5)`.
+- `bt.pp.philr`: PhILR balances along the phylogeny into `obsm["X_philr"]`,
+  equal to `philr::philr` with its default weights.
+- `bt.da.linda`, `bt.da.ancombc2`: differential abundance by LinDA (a port
+  of `MicrobiomeStat::linda`) and ANCOM-BC2 (through scikit-bio), both in
+  Python and compared with R in golden tests.
+- `bt.da.aldex2`, `bt.da.maaslin3`: ALDEx2 and MaAsLin 3's abundance model,
+  run in R through rpy2. They need R, the R package, and the new extra
+  `biotapy[r]` (`rpy2>=3.6.8`, GPL-2.0-or-later, built against your R).
+- One result table for every method: log2 `effect`, `se`, `pvalue`,
+  Benjamini-Hochberg `qvalue` over the features the method tested,
+  `direction`, `method` and `contrast`. Methods take `group`, `covariates`
+  and `reference` instead of a formula, and never filter features.
+- `bt.da.consensus` and `bt.pl.consensus`: count, per feature, the methods
+  that call it and whether they agree on its direction, as a table and as a
+  dot matrix.
+- A differential abundance guide, a page per method, and a tutorial on the
+  GlobalPatterns genera.
+- asv benchmarks for `pp.philr`, `da.linda` and `da.ancombc2`, with their
+  baselines in the docs.
+
+### Changed
+
+- The readers infer `x_kind` `"counts"` only for non-negative whole numbers,
+  and the functions that need counts (`pp.rarefy`, `tl.alpha`'s
+  `observed_features` and `chao1`, weighted `tl.unifrac`) refuse a table with
+  a negative value.
+- `bt.pl.bar` and `bt.pl.heatmap` raise `ValueError` for a table with a
+  negative value, such as `layers["clr"]`, instead of drawing it as
+  abundances.
+- The Coming-from-R page marks the phyloseq calls biotapy does not cover yet
+  "not in 0.3".
+
 ## [0.2.0] - 2026-10-05

 ### Added
diff --git a/pyproject.toml b/pyproject.toml
index 195dbd5..794071b 100644
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -4,7 +4,7 @@ requires = [ "hatchling>=1.27" ]

 [project]
 name = "biotapy"
-version = "0.2.0"
+version = "0.3.0"
 description = "mia-style microbiome toolkit for Python on AnnData/TreeData"
 readme = "README.md"
 license = "BSD-3-Clause"
```

  Each "Changed" line is a fix to a function 0.2.0 shipped: `dfa71fe` (`infer_x_kind` and
  `require_counts`: readers, `pp.rarefy`, `tl.alpha`, weighted `tl.unifrac`), `c920c7b`
  (`pl.bar`, `pl.heatmap`), and Task 3.12. If Checkpoint D's fix pass changed another 0.2.0
  function, add its line.
- [ ] **Step 5: README** (PyPI's project page). Change only the Status section and Installation:

````diff
diff --git a/README.md b/README.md
index cd201c8..f1ab12c 100644
--- a/README.md
+++ b/README.md
@@ -14,9 +14,9 @@ against R on real data.

 ## Status

-**biotapy 0.2 is an early release.** The API can still change between minor versions.
+**biotapy 0.3 is an early release.** The API can still change between minor versions.

-What 0.2 does (full signatures in the [API reference][api]):
+What 0.3 does (full signatures in the [API reference][api]):

 - **Readers**: `bt.io.read_biom` (BIOM 1.0/2.1), `bt.io.read_qiime2`
   (`.qza` artifacts, no QIIME 2 install needed), `bt.io.read_dada2`
@@ -34,17 +34,24 @@ What 0.2 does (full signatures in the [API reference][api]):
   and `bt.datasets.enzyme` (the ENZYME EC hierarchy), downloaded and cached
   on first use.
 - **Preprocessing**: `bt.pp.relative`, `bt.pp.tax_glom`,
-  `bt.pp.filter_features`, `bt.pp.filter_samples` and `bt.pp.rarefy`.
+  `bt.pp.filter_features`, `bt.pp.filter_samples` and `bt.pp.rarefy`, and
+  the compositional transforms `bt.pp.clr` and `bt.pp.philr` (checked against
+  vegan and philr).
 - **Function**: `bt.fn.load_hierarchy`, `bt.fn.func_glom` and `bt.fn.renorm`
   (checked against HUMAnN's own output), `bt.fn.contributions` and
   `bt.fn.functional_redundancy`.
 - **Tools**: `bt.tl.alpha`, `bt.tl.beta`, `bt.tl.unifrac`, `bt.tl.pcoa`,
   `bt.tl.nmds` and `bt.tl.permanova`, each checked against R on real data.
+- **Differential abundance**: `bt.da.linda` and `bt.da.ancombc2` in Python,
+  `bt.da.aldex2` and `bt.da.maaslin3` through R (the `r` extra below), all
+  returning one result table and checked against their R packages, and
+  `bt.da.consensus`, which reports where the methods agree.
 - **Plots**: `bt.pl.bar`, `bt.pl.richness`, `bt.pl.ordination`,
-  `bt.pl.scree`, `bt.pl.heatmap` and `bt.pl.contributions`.
+  `bt.pl.scree`, `bt.pl.heatmap`, `bt.pl.contributions` and
+  `bt.pl.consensus`.

-Next, in 0.3: compositional transforms (CLR, PhILR) and differential abundance
-(ANCOM-BC, LinDA, ALDEx2, MaAsLin 3) behind one result format. See the
+Next, in 0.4: multi-omics conventions on MuData, leak-free scikit-learn
+transformers, a PyTorch loader and an interface for embedding models. See the
 [roadmap][roadmap]; no dates are promised.

 ## Installation
@@ -71,6 +78,17 @@ pixi add --pypi biotapy
 The development version installs straight from GitHub:
 `pip install git+https://github.com/pedrocr83/biotapy.git`.

+ALDEx2 and MaAsLin 3 run in R. They need R with the R packages
+(`BiocManager::install(c("ALDEx2", "maaslin3"))`) and the `r` extra, which
+builds rpy2 against that R (Linux and macOS):
+
+```bash
+pip install 'biotapy[r]'
+```
+
+rpy2 is GPL-2.0-or-later and the R packages carry their own licences; biotapy
+does not ship any of them.
+
 On Python 3.14, the `biom-format` dependency has no wheels yet, so it is built
 from source and needs a C compiler until biom-format publishes 3.14 wheels
 ([biocore/biom-format#1004][biom-format-1004]).
````

- [ ] **Step 6: Build check.**
  - `rm -rf dist && uv build && uvx twine check --strict dist/*`: both files print `PASSED`.
  - `tar tzf dist/biotapy-0.3.0.tar.gz | grep -c benchmarks` prints `0`;
    `tar tzf dist/biotapy-0.3.0.tar.gz | grep -c "docs/methods\|tutorials/differential"` prints `7`.
  - `unzip -p dist/biotapy-0.3.0-py3-none-any.whl 'biotapy-0.3.0.dist-info/METADATA' | grep -E "^(Version|Requires-Dist: (rpy2|mudata)|Provides-Extra)"`
    prints `Version: 0.3.0`, `Requires-Dist: mudata>=0.4`, `Provides-Extra: r` and
    `Requires-Dist: rpy2>=3.6.8; extra == 'r'`.
  - The sdist's own tests (playbook step 3): extract it under `<scratchpad>/sdist`, then inside it
    run `BIOTAPY_DATA_DIR=<scratchpad>/pooch uv run --group test pytest -q -p no:cacheprovider`
    -> `1145 passed, 53 deselected` (the sdist leaves out `tests/test_ci.py`, the knowledge tests,
    `tests/humann`, `tests/r`).
  - `rm -rf dist <scratchpad>/sdist`.
- [ ] **Step 7: Gate and commit.**

```bash
git add pyproject.toml CHANGELOG.md README.md
git commit -m "chore: release 0.3.0

Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
git status --short                                   # empty
uvx prek run --all-files                             # all Passed
uv run --group test pytest -q -W error::UserWarning  # 1238 passed, 53 deselected (plus the fix pass's tests)
uv run python -c "import biotapy; print(biotapy.__version__)"  # 0.3.0
```

- [ ] **Step 8: Knowledge, second commit.**
  - `bash scripts/knowledge_stale.sh --against HEAD`. The bump re-stales the concepts whose
    `paths` hold `pyproject.toml` or `CHANGELOG.md`; in the prototype: `phase-0-foundation`,
    `module-boundaries`, `tree-access`, `cut-a-release`. Check each against the bump, which changes
    nothing they state, and bump `commit` and `generated`.
  - `playbooks/cut-a-release.md` step 3, the check this release added, and the roadmap's exit
    gate (items 1 and 3 ticked with Checkpoint D's PR number, Test run id and docs job id; item 2
    was ticked by the plan commit, and item 1 already carries the wording of slice 3D decision 2)
    and Checkpoint D boxes:

````diff
diff --git a/.knowledge/playbooks/cut-a-release.md b/.knowledge/playbooks/cut-a-release.md
index c8681a6..81a47eb 100644
--- a/.knowledge/playbooks/cut-a-release.md
+++ b/.knowledge/playbooks/cut-a-release.md
@@ -5,8 +5,8 @@ description: Bump the version, move the changelog entry, tag, publish a GitHub r
 tags: [release, workflow]
 status: stable
 paths: ["pyproject.toml", "CHANGELOG.md", ".github/workflows/release.yaml"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-06T14:12:26Z }
-commit: 64fe39d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:50:00Z }
+commit: 93f47c3
 sources:
   - id: trusted-publishing
     resource: https://docs.pypi.org/trusted-publishers/
@@ -37,7 +37,9 @@ A phase exit gate asks for a release (0.0.1 in Phase 0, 0.1 in Phase 1, ...).
 3. Commit `chore: release X.Y.Z` and merge it to `master` through a PR
    (merge commit, not squash). Before the PR, build and run `pytest` from the
    extracted sdist (`uv build --sdist`, `tar xzf`, then `pytest` inside it); it
-   must pass with no errors.
+   must pass with no errors. Check the wheel's `METADATA` too: `Version`, every
+   runtime `Requires-Dist`, and a `Provides-Extra` line for each extra (0.3.0
+   added `r`).
 4. With explicit user approval for each (rules.md R13.3), tag the merged commit
    and push only the tag:
    ```bash
diff --git a/.knowledge/roadmap/phase-3-stats.md b/.knowledge/roadmap/phase-3-stats.md
index e8c8124..58cb3ad 100644
--- a/.knowledge/roadmap/phase-3-stats.md
+++ b/.knowledge/roadmap/phase-3-stats.md
@@ -9,8 +9,8 @@ phase_state: in-progress
 effort: ~4 weeks part-time
 depends_on: [/roadmap/phase-2-function.md]
 paths: ["src/biotapy/da/**", "src/biotapy/pp/**", "src/biotapy/pl/**"]
-generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:40:00Z }
-commit: d6d82ef
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-07T08:50:00Z }
+commit: 93f47c3
 sources:
   - id: spec
     resource: ../../plan.md
@@ -610,14 +610,16 @@ columns onto it.
 - [ ] 3.15 Release 0.3.0

 # Exit gate
-- [ ] Consensus report on one benchmark dataset (GlobalPatterns genus, human
-  vs environmental), executed notebook in docs (CI docs job).
+- [x] Consensus report on one benchmark dataset (GlobalPatterns genus,
+  human-associated vs the other samples: environmental and mock communities),
+  executed notebook in docs (CI docs job). Proven on PR #<Checkpoint D PR>:
+  docs job run <id> executed `tutorials/differential_abundance.md`.
 - [x] Per-method agreement with the R reference, per
   [r-golden-parity](/contracts/r-golden-parity.md): `pp.clr`, `pp.philr`,
   `da.linda`, `da.ancombc2` in the network job; `da.aldex2`, `da.maaslin3` in
   the `r-bridge` job. Proven on PR #23 (run `37479113322`): `network`
   job `36 passed`, `r-bridge` job `24 passed`.
-- [ ] All Phase 1 and 2 gates still green.
+- [x] All Phase 1 and 2 gates still green. Proven on PR #<Checkpoint D PR>: every CI job green (run <id>).

 # Risks
 - **rpy2 on Linux builds from source and lists Python <= 3.13** -> the CI job
````

  - Also tick Checkpoint D's boxes and its checklist line, with the PR, merge commit and Read the
    Docs check recorded; tick Task 3.15 Steps 1-8.
  - Log, under `## <date> (release 0.3.0)` at the top:

```markdown
## <date> (release 0.3.0)
- **Update**: [cut-a-release](playbooks/cut-a-release.md) step 3: check the wheel's `METADATA` for the version, the runtime requirements and a `Provides-Extra` line per extra.
- **Verification**: re-checked against the 0.3.0 version bump and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [module-boundaries](contracts/module-boundaries.md), [tree-access](contracts/tree-access.md).
- **Update**: [phase-3-stats](roadmap/phase-3-stats.md) ticks the exit gate's notebook and "all gates green" items (docs job <id> and Test run <id> on PR #<n>), Checkpoint D and Task 3.15 Steps 1-8.
```

  - `bash scripts/knowledge_stale.sh --against HEAD` prints `25 current, 0 stale, 12 uncheckable`.
  - Commit `docs(knowledge): refresh concepts for the 0.3.0 release` (with the attribution line),
    staging `.knowledge/playbooks/cut-a-release.md`, `.knowledge/roadmap/phase-0-foundation.md`,
    `.knowledge/contracts/module-boundaries.md`, `.knowledge/contracts/tree-access.md`,
    `.knowledge/roadmap/phase-3-stats.md` and `.knowledge/log.md`; gate as Step 7.
- [ ] **Step 9: PR and merge** (approved at Checkpoint D): push `release-0.3.0`, open the PR,
  wait for green CI (every job, `r-bridge` and `docs` included), merge with a merge commit (not a
  squash).
- [ ] **Step 10: STOP AND ASK.** Tagging and releasing publishes to PyPI, and that cannot be
  undone. Ask the user, in one message, for explicit approval to:
  1. tag the merged `master` commit `v0.3.0`;
  2. push only that tag;
  3. create the GitHub release `0.3.0`, which triggers the upload;
  4. afterwards, push `close-phase-3` (Step 13), open its PR and merge-commit it on green.

  Wait for a yes that names this release. An earlier approval does not count.
- [ ] **Step 11: Tag and release**, only after Step 10's approval:

```bash
git switch master && git pull --ff-only
git log -1 --format=%s   # Merge pull request #<n> from pedrocr83/release-0.3.0
git tag v0.3.0
git push origin v0.3.0
awk '/^## \[0.3.0\]/{on=1;next} /^## \[/{on=0} on' CHANGELOG.md > <scratchpad>/notes-0.3.0.md
gh release create v0.3.0 --title "0.3.0" --notes-file <scratchpad>/notes-0.3.0.md
gh run list --workflow release.yaml --limit 3 --json databaseId,headBranch,status
```

  - Pick the run whose `headBranch` is `v0.3.0`, then `gh run watch <its databaseId>
    --exit-status`. If the `pypi` environment asks for a reviewer, the user approves it in the run.
  - If the run fails before the upload: fix `master`, then delete the release and tag (`gh release
    delete v0.3.0 --cleanup-tag`) and re-tag, each only with the user's approval (playbook, Common
    mistakes).
  - If it fails after the upload, report it: the version is spent, and the fix is 0.3.1.
- [ ] **Step 12: Verify** (playbook Verification). Both print `0.3.0`:

```bash
curl -s https://pypi.org/pypi/biotapy/json | python3 -c "import json,sys; print(json.load(sys.stdin)['info']['version'])"
uv run --no-project --with biotapy==0.3.0 python -c "import biotapy as bt; print(bt.__version__, bt.da.__all__, 'philr' in bt.pp.__all__)"
```

  The second also prints `['aldex2', 'ancombc2', 'consensus', 'linda', 'maaslin3'] True`, and
  `curl -s https://pypi.org/pypi/biotapy/0.3.0/json | python3 -c "import json,sys; print(json.load(sys.stdin)['info']['provides_extra'])"`
  prints `['r']`.
- [ ] **Step 13: Close Phase 3** (approved in Step 10), as `close-phase-2` did:
  - `git switch master && git pull --ff-only && git switch -c close-phase-3`;
  - `roadmap/phase-3-stats.md`: `phase_state: in-progress` becomes `phase_state: done`; tick Task
    3.15 Steps 9-13 and the checklist's `3.15 Release 0.3.0`; bump `generated` and `commit`;
  - `roadmap/phase-4-ml-multiomics.md`: `phase_state: not-started` becomes `phase_state:
    in-progress` (its entry criteria: Phase 3's exit gate and 0.3 released); bump `generated` and
    `commit`;
  - `roadmap/index.md`: Phase 4 moves to "# Active phase" with `**phase_state: in-progress**`;
    Phase 3 moves to "# Phases" after Phase 2 as `**phase_state: done** (0.3.0 on PyPI, <date>)`;
  - log, under `## <date> (release 0.3.0)`:
    `- **Update**: Phase 3 closed after biotapy 0.3.0 reached PyPI (tag v0.3.0, release workflow run <id>). [phase-3-stats](roadmap/phase-3-stats.md) is \`phase_state: done\` with every Task 3.15 step and exit-gate item ticked; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) is \`phase_state: in-progress\`; the [roadmap index](roadmap/index.md) lists Phase 4 as active.`
  - gates (`uvx prek run --all-files`; `uv run --group test pytest -q -W error::UserWarning`, whose
    `tests/test_knowledge_bundle.py` checks every concept's frontmatter and index entry); commit
    `docs(knowledge): close Phase 3 after the 0.3.0 release` (with the attribution line); push,
    PR, merge-commit on green.
  - Phase 4 then starts at its own expansion (rules.md R1.2a): ask the user before expanding it.

### Slice 3D decisions for the user

Recommended answer first. Items 1-3 change what the outline or design note 10 said; 4-9 settle
open questions; 10-13 are outward-facing or need fresh approval.

1. **Method pages in `docs/methods/`** (an index and pages for LinDA, ANCOM-BC2, ALDEx2,
   MaAsLin 3 and the consensus), in the "User guide" toctree, and **the DA guide restructured**:
   its four method sections move, text kept, into the pages; it gains a "Choosing methods" table;
   the generic bridge behaviour moves into "Methods that run in R". No `src/` change in the plan
   (two approved docstring-prose commits followed): docstrings keep linking the guide. Alternative: Phase 2's convention, method text only in the guide (no
   `docs/methods/`), which the outline departed from.
2. **The tutorial uses the `benchmark` fixture's contrast exactly**: 9 human-associated samples
   against the 17 others, which are 14 environmental samples and 3 mock communities (design note
   10's "8 vs 18" was a miscount), `reference="other"`; the note says what "other" holds and that
   the contrast is a demonstration. Its numbers then equal the goldens' and the `r` test's.
   Alternative: drop the mock samples so "other" is purely environmental (new R runs, new
   constants, and the tutorial no longer matches the goldens).
3. **The R bridges in the tutorial: a fenced, non-executed code block plus static tables of the
   counts the `r-bridge` job checks, kept in step by `FOUR_METHOD_*` constants shared by the `r`
   test and a new default-run test** (`test_the_tutorial_quotes_the_four_method_counts`). The page
   says the block is not run and why. Alternatives: recorded outputs (`.ipynb` with stored outputs
   or a `skip-execution` tag; per-page execution metadata is forbidden by an existing test, and
   stored outputs are unchecked); a CSV written by the `r-bridge` job and read by an executed cell
   (a new artifact and CI change).
4. **No workflow and no Read the Docs change**: Read the Docs keeps executing every notebook
   itself with no pooch cache, as Phase 2 decided; the tutorial adds about 9 s and no download
   (the phyloseq vignette already fetches GlobalPatterns in the same build).
5. **Citations as plain text with a DOI link**, as the function guide does; no
   `references.bib` entries or `{cite}` roles. Alternative: sphinxcontrib-bibtex citations
   collected on the References page (a new convention across pages).
6. **Coming from R:** the seven Phase 3 rows pinned by a test; "not in 0.2" becomes "not in 0.3"
   with the page sentence, the test and the phyloseq vignette, whose "multiple testing and
   differential abundance (biotapy 0.3)" bullet now says phyloseq points to DESeq2, which biotapy
   does not wrap, and links the new tutorial. No `phyloseq_to_deseq2` row. Alternative: a
   version-free "not yet in biotapy" label that never needs moving.
7. **Benchmarks on synthetic data, not GlobalPatterns**: `pp.philr` at 1,000 x 50,000 (Phase 1's
   table construction) and `da.linda`/`da.ancombc2` at 2,000 x 10,000, 30% non-zero, each timed
   and peak-memory measured; ALDEx2 and MaAsLin 3 not benchmarked (R's time; asv has no R). The
   outline's "pp.philr on GlobalPatterns" would put pooch downloads inside asv.
8. **Approved 2026-10-07 as a fix, not a report: Task 3.13 rewords `da.ancombc2`'s docstring to "4.5x to 7.5x depending on table shape" in its own `docs(da)` commit. Original finding:** `da.ancombc2`'s docstring says peak memory is about seven dense
   arrays; on the asv table the call adds 4.5x.** The performance page states both. Suggested:
   leave the docstring (its figure is measured on other shapes) or reword it to "five to seven" in
   a later `docs(da)` change.
9. **Task 3.14 is small**: `modules/da.md` was rewritten for the bridges at Checkpoint C; 3.14
   adds the tutorial and method pages to its Verification, the asv memory number, one gotcha (the
   tutorial's quoted counts), `add-a-function` step 8 (a `da` method gets a method page) and four
   verification bumps.
10. **Release 0.3.0 text:** the CHANGELOG entries in Task 3.15 Step 4, written at release because
    `## [Unreleased]` is empty again, with "Changed" lines for the 0.2.0 behaviours that changed
    (counts must be non-negative; `pl.bar`/`pl.heatmap` refuse negative tables; "not in 0.3");
    README: "biotapy 0.3", `pp.clr`/`pp.philr`, a Differential abundance line, `pl.consensus`, the
    `r` extra with its install line and GPL note under Installation, and "Next, in 0.4" from Phase
    4's description.
11. **Approvals, as for 0.2.0:** the `release-0.3.0` push/PR/merge with the slice 3D review
    (Checkpoint D); the tag, the tag push, the GitHub release (the PyPI upload) and the
    `close-phase-3` push/PR/merge in one message after the release PR merges (Task 3.15 Step 10).
    On close, Phase 3 `done` and Phase 4 `in-progress`.
12. **Box ticking:** the plan commit ticks Checkpoint C (PR #23 merged, run `37479113322` green)
    and exit-gate item 2 (PR #23's `network` `36 passed` and `r-bridge` `24 passed`); the 3.15
    knowledge commit ticks items 1 and 3 and Checkpoint D with the `phase-3d` PR's run ids.
13. **Checkpoint C's last box ("ask the user to review slice 3C") is ticked by approving this
    plan**, recorded with the approval date. Say so if the 3C review is still open.

### Slice 3D self-review

1. **Spec coverage.**

   | Brief or outline item | Where it lands |
   |---|---|
   | 3.10 method pages (model, units, R defaults and changes, references) and the DA guide (filter once, choose methods first, consensus table, plot, `r` extra and GPL note) | Task 3.10; design 1; decisions 1, 5 |
   | 3.10b exit-gate notebook, executed by myst-nb in the docs job with the pooch cache | Task 3.10b; design 2; decisions 2-4 |
   | The bridges in a docs build without R, no silent fallback | design 2; decision 3; review focus 1 |
   | Read the Docs executes or reuses the docs job's output | design 2 (it executes itself, no cache, unchanged); decision 4 |
   | 3.12 new rows from docstrings, pinned; "not in 0.3" (cut-a-release 2c) | Task 3.12; decision 6 |
   | 3.13 asv for `pp.philr`, `da.linda`, `da.ancombc2`, baselines with the machine, no optimisation | Task 3.13; design 4; decisions 7-8 |
   | Checkpoint D with 3.14 knowledge (mostly done at C) | Checkpoint D; Task 3.14; decision 9 |
   | 3.15 version, CHANGELOG, README, release PR; tag, push, GitHub release, PyPI as ask-the-user steps | Task 3.15 Steps 4-5, 9, 10-13; decisions 10-11 |
   | Exit-gate boxes, which commit ticks which | design 6; decision 12 |

2. **Placeholder scan.** Every page, test, benchmark and diff block is rendered by `git show`
   from the scratch commits that passed the gates, not retyped. Values left for run time, each
   named where it is used: the implementer's model id and UTC time in `generated`, `HEAD`'s short
   hash in `commit`, the dates in log headings, the CHANGELOG and the Checkpoint C record, Task
   3.13's benchmark numbers, commit and load averages (Step 4's run replaces the prototype's), PR
   numbers, the Test, docs-job and release workflow run ids, and `<scratchpad>`. The `<...>`
   strings in the 3.15 knowledge diff are those run-time values.
3. **Type consistency.** The tutorial calls each function with the signature its Interfaces line
   gives; the constants `FOUR_METHOD_CALLS` (method name -> calls), `FOUR_METHOD_N_SIGNIFICANT`
   (count -> genera) and `FOUR_METHOD_CONSENSUS` match the tutorial's table rows and sentence the
   test builds; `synthetic_genera`, `N_PHILR_OBS` match `da.py`'s and `pp.py`'s imports; the label
   `da-methods-in-r` is defined in the guide (3.10) before the tutorial uses it (3.10b).
4. **Review focus.** Items 1, 3 and 4 name tests that exist in the blocks above or in the
   repository (`test_the_tutorial_quotes_the_four_method_counts`,
   `test_four_methods_on_the_exit_gate_data`, `test_results_with_different_contrasts_raise`,
   `test_transforms_and_da_methods_map_to_their_r_functions`; checked by searching); items 2 and 5
   name the review check and the build check that cover them.
5. **Known residual risks and what was not verified.**
   - CI was not run on the slice: the `docs` job's execution of the tutorial and the `r-bridge`
     job's run of the edited `r` test happen first on the PR. Locally both passed (docs build;
     `-m r` in the R image, `24 passed`).
   - Read the Docs was not built; Checkpoint D checks the two pages after the merge.
   - The method pages' R-default tables are read from `formals()` of the image's versions
     (Bioconductor 3.22); a later Bioconductor release can change a default without any test
     noticing. The pages name the versions.
   - The benchmark numbers come from one run on one laptop under a load average of 2-5; they are a
     record, not a gate.
   - The tutorial's prose names data facts (636 genera, 36 untested, 9 against 17) that golden
     tests pin; the live outputs carry every other number.
   - Coverage was not re-measured: the slice changes no `src/` code, only docstring prose (the last measurement,
     slice 3C's, stands).

# Decisions for the user
Each changes a contract, rule, dependency, CI or a roadmap signature, or is
a judgement call. Recommended answer first.

1. **Slices and order:** 3A transforms -> 3B native DA, consensus and plot
   -> 3C R bridges and CI -> 3D docs and release (design note 1).
2. **PhILR output:** `obsm["X_philr"]` (a samples x nodes DataFrame) on a
   copy of the TreeData, returned as `TreeData`; not a new AnnData of
   balances. Contract rows in data-model-slots (design note 3).
3. **PhILR tree rules:** a node with more than two children raises, the root
   included (no invented root balance); one-child nodes are skipped in
   order; `_core.get_skbio_tree` gains `split_root=False`. Alternative:
   resolve a wide root as UniFrac does (first child vs the rest), which lets
   `pp.philr(bt.datasets.toy())` run but adds an arbitrary balance.
4. **PhILR weights:** uniform only in 0.3; balance signs follow philr.
5. **`pp.clr` signature:** no `layer` argument; the layer is dense; a
   pseudocount larger than the smallest non-zero value warns.
6. **Golden image additions, one commit each:** philr + `libuv1` (3A),
   ANCOMBC and MicrobiomeStat (3B), ALDEx2 and maaslin3 (3C); GPL R
   packages live only in the image.
7. **No formulas in 0.3:** every `da` method takes `group`, `covariates`,
   `reference` instead of `formula`; two-level or numeric `group`; missing
   values raise; no random effects or interactions (design note 4).
8. **The result schema:** columns `effect, se, pvalue, qvalue, direction,
   method, contrast` (the roadmap's six plus `contrast`); log2; BH
   everywhere; `q < alpha` strict; methods never filter features (design
   note 5).
9. **Consensus is two-step:** `da.consensus(results, *, alpha=0.05,
   min_methods=2)` combines tables the user computed, replacing
   `da.consensus(adata, formula, *, methods, alpha, min_methods, n_jobs)`;
   `n_jobs` and joblib dropped; the agreement rule in a Decision concept
   (design note 6). function-shape (and rules.md R3.2's "annotate `data`")
   gain one sentence: `da.consensus` takes result tables and `pl.consensus`
   a consensus table, not an AnnData.
10. **`da.ancombc` -> `da.ancombc2`**, wrapping scikit-bio's `ancombc2`;
    `alpha` dropped from every method.
11. **LinDA:** native, fixed effects only, no winsorisation (golden vs
    `MicrobiomeStat::linda(is.winsor = FALSE)`), no statsmodels.
12. **ALDEx2 through the rpy2 bridge**, not scikit-bio's `dirmult_ttest`;
    two groups, no covariates; `effect = diff.btw`, `qvalue = we.eBH`
    [superseded by slice 3C decision 3: `qvalue` is BH of `we.ep`].
13. **MaAsLin 3 through the bridge, abundance model only.**
14. **Extra `r = ["rpy2>=3.6.8"]`** and its licensing stance: biotapy (BSD-3)
    imports GPL rpy2 only when installed by the user, never bundles it;
    docs say the extra pulls GPL software.
15. **CI job `r-bridge`** (Linux, Python 3.13, R 4.5.3 / Bioconductor 3.22,
    `-m r`), added to `check.needs` with a `tests/test_ci.py` test.
    Alternative: no R in CI, bridge tests local only.
16. **Exit-gate notebook (slice 3D decision 2):** GlobalPatterns genus, human-associated vs the other samples (environmental and mock communities);
    ANCOM-BC2 and LinDA live, the bridges as a non-executed block (no R in
    the docs builds); the four-method consensus runs in the `r-bridge` job.
17. **Task changes:** new 3.0, 3.10b, 3.12, 3.13, 3.14, 3.15; 3.3 delivered
    inside 3.5; execution order 3.5 -> 3.4 -> 3.8 -> 3.9 in 3B.
18. **mypy config:** `skbio.stats.composition._base` joins
    `untyped_calls_exclude` (config only, no dependency).
19. **Frontmatter:** the new `description`, `paths` gains `src/biotapy/pl/**`,
    and the sources gain philr, LinDA, Nearing and Pelto.
20. **Branch pushes:** approve push, PR and merge on green for the Phase 3
    slice branches (R13.3); the standing approval covered Phase 2 only.

# Self-review
Run against the brief, the roadmap concept and the writing-plans checklist.

1. **Spec coverage.**

   | Roadmap item | Where it lands |
   |---|---|
   | Design note "agree" | Design note 6; Decision concept in 3.8 |
   | Design note "one result schema" | Design note 5; 3.5 (with 3.3) |
   | 3.1 `pp.clr` | Task 3.1 (golden 3.0) |
   | 3.2 `pp.philr` | Task 3.2 (golden 3.0) |
   | 3.3 schema | inside 3.5 (3B) |
   | 3.4 ANCOM-BC | 3.4 `da.ancombc2` (3B) |
   | 3.5 LinDA | 3.5 (3B) |
   | 3.6 ALDEx2, 3.7 MaAsLin 3 | 3C |
   | 3.8 consensus, 3.9 plot | 3B |
   | 3.10 docs, `da` Module concept | 3D (3.10, 3.14) |
   | 3.11 CI | 3C |
   | Exit gate: notebook | 3.10b; design note 10 |
   | Exit gate: R agreement | goldens in 3.0, 3.4, 3.5, 3.6, 3.7 |
   | Risk: rpy2 on Windows/macOS CI | Risks; CI job Linux-only |

   Brief questions 1-11 are answered in design notes 1-11. Gaps found and
   fixed while writing: the roadmap had no golden task for CLR/PhILR (3.0);
   the formula signatures could not say which coefficient a one-row table
   reports (design note 4); the one-call consensus broke R5's six-argument
   limit (design note 6); the docs builds have no R (design note 10).
2. **Placeholder scan.** Every slice 3A step carries the full file or the
   exact lines, rendered from the scratch clone's commits that passed every
   gate. Outline slices give a proposed answer to each open question and
   mark what is [UNVERIFIED].
3. **Type consistency.** Used throughout: `pseudocounted(adata, pseudocount,
   *, func)`; `get_skbio_tree(adata, *, split_root=True)`;
   `pp.philr -> TreeData` with `obsm["X_philr"]`; every `da` method
   `(adata, group, *, covariates=(), reference=None[, seed])`; the schema
   columns `effect, se, pvalue, qvalue, direction, method, contrast`;
   consensus columns `effect_<method>`, `qvalue_<method>`, `n_tested`,
   `n_significant`, `direction`, `consensus`, `conflict`.
4. **Review focus.** Items 1-3 name tests that exist in the rendered 3A
   code (checked by searching this file for each name); items 4-5 name tests
   that 3B's expansion must write under those names.
5. **Known residual risks.**
   - The 3.0 commit pair was gated only on its second state; the first
     changes only the Dockerfile and a contract sentence.
   - R parity of scikit-bio's `ancombc2`, the LinDA mode port and the
     bridges was not run (no R on the host; R ran only inside the image for
     3A's goldens).
   - [UNVERIFIED]: how rpy2 fails without R; rpy2 `py.typed`; maaslin3's
     `evaluate_only` and `reference` details; pak's Bioconductor choice for
     R 4.5.3; CI R job runtime; the GPL stance (a user ruling, not a fact).

[^spec]: Python Microbiome Toolkit development report, sections Positioning and Roadmap
[^maaslin3]: MaAsLin 3, Nature Methods
[^philr]: Silverman et al. 2017, A phylogenetic transform enhances analysis of compositional microbiota data, eLife
[^linda]: Zhou et al. 2022, LinDA, Genome Biology
[^nearing]: Nearing et al. 2022, Microbiome differential abundance methods produce different results across 38 datasets, Nature Communications
[^pelto]: Pelto et al. 2025, Elementary methods provide more replicable results in microbial differential abundance analysis, Briefings in Bioinformatics
