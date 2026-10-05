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
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T15:08:33Z }
commit: 2216894
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
> (recommended) or superpowers:executing-plans. Slice 3A has full TDD steps; slices
> 3B-3D are outlines, expanded into TDD steps (rules.md R1.2a) when each is reached.

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
     BH. Mixed Holm/BH q-values would make "significant" mean different
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
     we.ep`, `qvalue = we.eBH`. Seeding: `as_generator(seed)` draws one
     integer, passed to R's `set.seed` before the call (R's RNG is global to
     the embedded R; documented in `Notes`); `useMC = FALSE`.
   - **MaAsLin 3: `da.maaslin3(adata, group, *, covariates=(),
     reference=None, seed=None) -> pd.DataFrame`, rpy2 bridge.**[^maaslin3] Abundance
     model only (`evaluate_only = "abundance"`): the prevalence model reports
     log-odds, which cannot share an `effect` column with log2 fold changes.
     `min_prevalence = 0`, plots off, output to a temporary directory read
     back from `all_results.tsv`. `effect = coef`, `se = stderr`,
     `pvalue = pval_individual`, `qvalue = qval_individual`.
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
     coverage holds without R.
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
    human samples (Feces, Skin, Tongue: 8) vs environmental (18): research C
    ran three scikit-bio methods on exactly this in 1.3 s with 5-43 hits each
    and visible disagreement. It is already a pooch dataset (no new
    download, hash or licence) and holds counts, which every method needs.
    The notebook says the contrast is a demonstration, not a biological
    claim. It runs ANCOM-BC2 and LinDA live: the docs build (CI docs job,
    Read the Docs) has no R, so the ALDEx2 and MaAsLin 3 calls are shown as
    a non-executed code block with their extra, and the CI `r-bridge` job
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
  adds it.
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

Execution order inside 3B: **3.5 -> 3.4 -> 3.8 -> 3.9 -> Checkpoint B.**
LinDA first: it is native end to end, so the schema's every column (with a
real `se`) is fixed by code biotapy owns before a wrapper maps a library's
columns onto it.

# Tasks (checklist)
- [x] 3.0 CLR and PhILR golden files (philr in the R image)
- [x] 3.1 `pp.clr(adata, *, pseudocount=0.5) -> AnnData`
- [x] 3.2 `pp.philr(tdata, *, pseudocount=0.5) -> TreeData`
- [ ] Checkpoint A
- [ ] 3.3 Result schema `da/_schema.py` (delivered inside 3.5)
- [ ] 3.5 `da.linda(adata, group, *, covariates=(), reference=None) -> pd.DataFrame`
- [ ] 3.4 `da.ancombc2(adata, group, *, covariates=(), reference=None) -> pd.DataFrame`
- [ ] 3.8 `da.consensus(results, *, alpha=0.05, min_methods=2) -> pd.DataFrame` and the agreement decision
- [ ] 3.9 `pl.consensus(table, *, top=30, ax=None) -> Axes`
- [ ] Checkpoint B
- [ ] 3.6 `da.aldex2(adata, group, *, mc_samples=128, reference=None, seed=None) -> pd.DataFrame` and the extra `r`
- [ ] 3.7 `da.maaslin3(adata, group, *, covariates=(), reference=None, seed=None) -> pd.DataFrame`
- [ ] 3.11 CI job `r-bridge` for `-m r` tests
- [ ] Checkpoint C
- [ ] 3.10 Method pages in `docs/methods/` and the DA guide
- [ ] 3.10b `docs/tutorials/differential_abundance.md`, the exit-gate notebook
- [ ] 3.12 Coming-from-R check
- [ ] 3.13 asv benchmarks for `pp.philr`, `da.linda`, `da.ancombc2`
- [ ] Checkpoint D
- [ ] 3.14 Knowledge (Module concept `da`)
- [ ] 3.15 Release 0.3.0

# Exit gate
- [ ] Consensus report on one benchmark dataset (GlobalPatterns genus, human
  vs environmental), executed notebook in docs (CI docs job).
- [ ] Per-method agreement with the R reference, per
  [r-golden-parity](/contracts/r-golden-parity.md): `pp.clr`, `pp.philr`,
  `da.linda`, `da.ancombc2` in the network job; `da.aldex2`, `da.maaslin3` in
  the `r-bridge` job.
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
      8 bytes x samples x features. Peak memory is about three such arrays (the
      dense copy of ``X``, the transform's temporaries and its output; 3.4x measured
      on a 400 x 500 table), so budget for that on large tables.

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

  `layers["clr"]` is dense: CLR has no zeros, so it takes 8 bytes per sample and feature, and the call peaks at about three such arrays.
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
- [ ] Review the whole slice (superpowers:requesting-code-review) against
  every contract, pure-by-default, the Phase 3 and slice 3A review focus;
  then a fix pass, one commit per finding, each with a test.
- [ ] Run the exit-gate check for 3A: `uv run --group test pytest -m "golden or network" tests/pp -q`
  and the full `uv run --group test pytest`.
- [ ] Knowledge (codebase-map templates; R12.2-R12.4):
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
## Slice 3B - Native DA and consensus (outline)

**Goal:** without R, a user compares two groups with ANCOM-BC2 and LinDA,
gets two tables with the same columns and units, and sees in one table and
one plot where they agree.

### Slice 3B design (proposed; expanded into TDD steps when reached)
- **Where the code goes.** `da/__init__.py` (imports and `__all__` only),
  `da/_schema.py` (`COLUMNS`, `validate_result`, `make_result(features,
  effect, se, pvalue, *, method, contrast) -> pd.DataFrame`, which computes
  BH `qvalue` and `direction`), `da/_design.py` (`design(adata, group, *,
  covariates, reference) -> tuple[pd.DataFrame, str]`: validates the `obs`
  columns, returns the design columns and the `contrast` text; shared by the
  four methods, so it lives in `da`, R4.3), `da/_linda.py`, `da/_ancombc.py`,
  `da/_consensus.py`, `pl/_consensus.py`. Tests mirror them under
  `tests/da/` and `tests/pl/test_consensus.py`.
- **3.5 `da.linda` (with the schema, 3.3).** Design note 7's algorithm.
  - `require_counts` (LinDA's `count` type; proportions are out of 0.3).
  - The imputation branch: zeros become `N_i / N_max`-scaled values as in
    MicrobiomeStat's `linda` (lines read in research A section 5); the
    branch chosen is logged at DEBUG, not warned (R7.3).
  - The mean-shift mode: about 12 lines of numpy in `_linda.py`, a single-use
    helper allowed by R4.4 to keep `linda` under 30 statements.
  - Golden: the image gains MicrobiomeStat (own commit); `export_golden.R`
    writes `linda(t(counts), meta, "~host", feature.dat.type = "count",
    is.winsor = FALSE, prev.filter = 0)$output$hostHuman` for the benchmark
    data (GlobalPatterns genus, `min_prevalence = 0.2`; the R side filters
    with the same rule as `pp.filter_features`) and one numeric-covariate
    case. Compared elementwise (LinDA is deterministic), `rtol=1e-6` if the
    mode iteration's stopping rule needs it, with the reason in the test.
  - Tests named in the review focus: `test_reference_sets_the_sign`,
    `test_missing_group_value_raises`, `test_no_feature_is_dropped`; plus
    three-level group raises, numeric group is per SD, all-zero feature
    gives NaN or a finite row (decided by the R output), single-sample
    raises (no df), purity, Hypothesis property: `direction ==
    sign(effect)` and `qvalue >= pvalue`.
- **3.4 `da.ancombc2`.** `skbio.stats.composition.ancombc2(table,
  metadata, formula, p_adjust="bh")` with a DataFrame of dense `X`, the
  metadata from `design()` (group as a two-category `Categorical` with the
  reference first, so patsy's `[T.<level>]` row is the effect), formula
  built with `Q("...")` quoting; the covariate row selected from the
  `(FeatureID, Covariate)` MultiIndex; `Log(FC)` and `SE` divided by
  `ln 2`. Golden: ANCOMBC in the image (own commit), settings of design note
  7. Tests: `test_qvalue_is_benjamini_hochberg` (equals
  `false_discovery_control` on its `pvalue`), the reference/missing tests,
  zero handling on `toy()`, purity, the sparse-to-dense memory note.
- **3.8 `da.consensus` and `decisions/da-consensus-agreement.md`.** Design
  note 6. The decision concept (`type: Decision`) records the definition,
  the options weighed (k-of-n with direction, intersection, rank
  aggregation, "agree among all with an estimate"), why strict `<`, why one
  BH, the literature split, and that the result is a robustness report.
  Tests: `test_q_equal_to_alpha_is_not_called`,
  `test_results_with_different_contrasts_raise`, conflict flagged and never
  consensus, NaN is not tested, union of features, `min_methods` bounds,
  duplicate method names raise, a Hypothesis property (`consensus` implies
  `n_significant >= min_methods` and `not conflict`).
- **3.9 `pl.consensus`.** Design note 9. Tests in the `tests/pl` style
  (Agg backend): returns the given `ax`; marker counts per kind on a
  hand-built table; `top` respected; an empty table (nothing called)
  raises `ValueError` saying nothing to draw, rather than an empty figure;
  purity of the table.
- **Contracts at 3B:** data-model-slots gains a "DA results" section (the
  schema table) and says `da` writes no slot; r-golden-parity's DA row gets
  the per-method rules measured in 3.4/3.5; function-shape's `data` bullet
  allows a result sequence for `da.consensus` and a table for
  `pl.consensus`; module-boundaries is unchanged (`da` is already a layer).
- **Open for 3B expansion:** whether scikit-bio's `ancombc2` equals R's
  (measured first); LinDA's behaviour on an all-zero feature (follow R).

## Slice 3C - R bridges (outline)

**Goal:** with `pip install 'biotapy[r]'` and R, a user runs ALDEx2 and
MaAsLin 3 from Python and gets the same schema; CI runs them on Linux.

### Slice 3C design (proposed)
- **Where the code goes.** `da/_r.py` (`r_packages(*names) -> ...`: imports
  `rpy2.robjects` through `import_optional(..., extra="r")`, loads R
  packages with `importr`, turns a failure into `ImportError` naming the
  package and the `BiocManager::install` line; `r_seed(seed)`: derives an
  `int` from `as_generator(seed)` and calls `set.seed`; `to_r_frame` /
  `from_r_frame` inside a local converter context), `da/_aldex2.py`,
  `da/_maaslin3.py`. `pyproject.toml` gains `[project.optional-dependencies]
  r = ["rpy2>=3.6.8"]` and the mypy override for rpy2.
- **3.6 `da.aldex2`.** `aldex.clr(reads, conds, mc.samples, denom = "all",
  useMC = FALSE)`, `aldex.ttest`, `aldex.effect`; features as rows; a
  two-level `group` only (more levels or `covariates` raise). Golden: ALDEx2
  in the image (own commit), `set.seed(20260927)`, the benchmark data;
  compared by sign agreement and Spearman correlation of `effect` (contract
  DA rule; numbers fixed when measured). Tests without R: canned
  `aldex` output -> schema mapping (`diff.btw` -> `effect`, `we.ep`,
  `we.eBH`), seed derivation, missing-rpy2 `ImportError` naming the extra,
  missing-R-package `ImportError` naming `BiocManager::install("ALDEx2")`.
- **3.7 `da.maaslin3`.** `maaslin3(input_data, input_metadata, output =
  <tempdir>, formula = "~ <group> + <covariates>", reference = "<group>,<ref>",
  evaluate_only = "abundance", min_prevalence = 0, plot_summary_plot =
  FALSE, plot_associations = FALSE, cores = 1, verbosity = "ERROR")`, then
  `all_results.tsv` filtered to `metadata == group` and `model ==
  "abundance"`. Golden: maaslin3 in the image (own commit). Same no-R tests.
  [UNVERIFIED] until expanded: `evaluate_only`'s exact behaviour, the
  `reference` string syntax and `name` values for a two-level factor in
  1.2.0.
- **3.11 CI `r-bridge` job.** Design note 8; `tests/test_ci.py` gains
  `test_r_bridge_job_runs_the_r_marker` and `check.needs` includes it. The
  job also runs the four-method consensus test on the benchmark data
  (marker `r` + `network`, the pooch cache step copied from the network job).
- **Decision update:** `optional-heavy-dependencies.md` records rpy2 (R9.2);
  `r-bridge-before-ports.md`'s consequence line is met.

## Slice 3D - Docs and release (outline)

### Slice 3D design (proposed)
- **3.10 docs.** `docs/methods/` with one page per method (`ancombc2.md`,
  `linda.md`, `aldex2.md`, `maaslin3.md`, `consensus.md`: model, units,
  what the R defaults are and which biotapy changes, references) and
  `docs/guide/differential_abundance.md` (filter once, choose methods
  first, read the consensus table, the plot; the `r` extra and its GPL
  note). Both listed in the docs toctrees; `api.md` gains a "Differential
  abundance" section.
- **3.10b tutorial.** Design note 10, executed by myst-nb in the docs job
  with the pooch cache, as Phase 2's function tutorial.
- **3.12 Coming-from-R check.** New rows come from docstrings
  (`ANCOMBC::ancombc2`, `MicrobiomeStat::linda`, `ALDEx2::aldex`,
  `maaslin3::maaslin3`, `philr::philr`, `vegan::decostand`); a test pins
  them as `tests/test_coming_from_r.py` pins the mia importers; the
  "not in 0.2" labels become "not in 0.3" where still true (cut-a-release
  step 2c).
- **3.13 benchmarks.** asv entries for `pp.philr` (GlobalPatterns, 2,572
  taxa), `da.linda`, `da.ancombc2` (benchmark data); baselines in
  `docs/performance.md`. Measurements only; no optimisation (R10.1).
- **Checkpoint D and 3.14 knowledge.** Create `.knowledge/modules/da.md`
  (`type: Module`; responsibility, entry points, invariants: the schema,
  BH, no filtering, no formulas, R imports inside functions; verification
  `uv run --group test pytest tests/da -q`; gotchas: patsy's alphabetical
  reference, Holm defaults, rpy2 on Linux); update `pl.md`, `core.md`,
  `optional-heavy-dependencies.md`; index entries and log lines.
- **3.15 release 0.3.0** per the cut-a-release playbook, after the user
  approves the tag (R13.3).

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
    two groups, no covariates; `effect = diff.btw`, `qvalue = we.eBH`.
13. **MaAsLin 3 through the bridge, abundance model only.**
14. **Extra `r = ["rpy2>=3.6.8"]`** and its licensing stance: biotapy (BSD-3)
    imports GPL rpy2 only when installed by the user, never bundles it;
    docs say the extra pulls GPL software.
15. **CI job `r-bridge`** (Linux, Python 3.13, R 4.5.3 / Bioconductor 3.22,
    `-m r`), added to `check.needs` with a `tests/test_ci.py` test.
    Alternative: no R in CI, bridge tests local only.
16. **Exit-gate notebook:** GlobalPatterns genus, human vs environmental;
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
