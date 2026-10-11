---
type: Phase
title: Phase 4 - ML-ready and multi-omics (0.4)
description: io.to_mudata over shared samples, tl.mmvec through scikit-bio, leak-free scikit-learn transformers (PrevalenceFilter, CLR), a PyTorch Dataset behind the extra torch, an entry-point interface for embedding models with MGM as the reference plugin, a leak-free cross-validation notebook.
tags: [roadmap, ml, multiomics]
status: stable
release: "0.4"
phase_state: done
effort: ~4-6 weeks part-time
depends_on: [/roadmap/phase-3-stats.md]
paths: ["src/biotapy/ml/**", "src/biotapy/tl/**", "src/biotapy/io/**", "src/biotapy/_core/**", "src/biotapy/datasets/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-10T23:56:11Z }
commit: 3be7a23
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: mmvec
    resource: https://doi.org/10.1038/s41592-019-0616-3
    title: Morton et al. 2019, Learning representations of microbe-metabolite interactions, Nature Methods
  - id: mgm
    resource: https://pmc.ncbi.nlm.nih.gov/articles/PMC13116254/
    title: MGM foundation model (microformer-mgm), HUST-NingKang-Lab
  - id: biomegpt
    resource: https://www.biorxiv.org/content/10.64898/2026.01.05.697599v1.full
    title: BiomeGPT preprint
  - id: mgm-notebook
    resource: https://github.com/HUST-NingKang-Lab/MGM/blob/bee1469fe13116d1d09e74a6169b29f6cd5a3c26/MGM_Interpretability_Guideline.ipynb
    title: MGM's interpretability notebook (sample embeddings)
---

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Every slice has full
> TDD steps.

**Goal:** 0.4 makes biotapy ML-ready and multi-omics: several data types over
the same samples become one MuData, mmvec relates microbes to metabolites,
preprocessing that learns from samples runs inside scikit-learn pipelines so
cross-validation does not leak, a table becomes a PyTorch dataset, and
foundation-model embeddings land in `obsm` through plugins.[^spec]

**Architecture:**
- `io.to_mudata(modalities)` builds the MuData; modality names are `taxa`,
  `function`, `function_by_taxon`, `metabolites`, `host`. Functions that pair
  modalities (`tl.mmvec`) check the samples are aligned and name
  `io.to_mudata` as the fix.
- `tl.mmvec` wraps scikit-bio 0.7.4's own mmvec (NumPy/SciPy L-BFGS, no
  TensorFlow) and returns a microbes x metabolites table.
- `ml` is a new top-layer subpackage (beside `pl` and `da`). Its scikit-learn
  transformers are classes, the one place rules.md R3.6 allows them:
  `PrevalenceFilter` (learns which features to keep, so it must be fitted
  per fold) and `CLR` (stateless). `ml.to_torch` and `ml.embed` are
  functions; torch and transformers are imported inside them through
  `import_optional` (R4.6).
- Embedding models are plugins found through the entry-point group
  `biotapy.embeddings`; biotapy registers one itself, MGM, whose MIT weights
  are fetched by pooch from the published wheel and never bundled.

**Tech stack:** anndata 0.13.4 · treedata 0.3.1 · mudata 0.4.1 · scikit-bio
0.7.4 (`stats.ordination.mmvec`, `stats.composition.clr`) · scikit-learn 1.9.1
(`SelectorMixin`, `OneToOneFeatureMixin`, `TransformerMixin`,
`validate_data`, `parametrize_with_checks`) · numpy 2.5.3 · scipy 1.18.1 ·
pandas 3.0.6 · torch >= 2.9 (extra `torch`, slice 4B) · transformers >= 5
(5.19.0 measured; extra `mgm`, slice 4C) · pooch · pytest/hypothesis.

**Spec:** [plan.md](../../plan.md). Contracts that bind every task:
[function-shape](/contracts/function-shape.md),
[data-model-slots](/contracts/data-model-slots.md),
[module-boundaries](/contracts/module-boundaries.md),
[tree-access](/contracts/tree-access.md). Decisions:
[pure-by-default](/decisions/pure-by-default.md),
[optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md),
[function-tables-as-mudata](/decisions/function-tables-as-mudata.md),
[treedata-as-container](/decisions/treedata-as-container.md). Facts marked
[V] below were run or read at source on 2026-10-08; [UNVERIFIED] marks the
rest, which the slice outlines resolve before code.

**How slice 4A was checked.** Every file in slice 4A was written into a
scratch clone of the repository at `872c16d` (master, after 0.3.0 was
closed) and replayed as one commit per task. Each committed state was gated
in a separate worktree checked out at that commit:

| Task state | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` |
|---|---|---|---|---|
| base `872c16d` | passed | 1238 passed, 54 deselected | 36 passed | build succeeded |
| 4.1 `feat(io)` | passed | 1256 passed, 54 deselected | 36 passed | build succeeded |
| 4.2 `feat(tl)` | passed | 1272 passed, 54 deselected | 36 passed | build succeeded |
| 4.A0 `refactor(core)` | passed | 1282 passed, 54 deselected | 36 passed | build succeeded |
| 4.3 `feat(ml)` | passed | 1410 passed, 2 skipped, 54 deselected | 36 passed | build succeeded |

- prek covers ruff check and format, `mypy --strict`, import-linter,
  pyproject-fmt, biome and zizmor (14 hooks). Every run exported
  `BIOTAPY_DATA_DIR` to a scratch pooch cache; `~/.cache/biotapy` was checked
  absent after each.
- The 2 skips in 4.3 are scikit-learn's `check_array_api_input` for the two
  transformers, which raises `SkipTest` unless `SCIPY_ARRAY_API` is set.
- Coverage on the final state (`coverage run -m pytest`):
  `io/_mudata.py`, `tl/_mmvec.py`, `_core/_composition.py`,
  `ml/_transformers.py`, `pp/_transform.py` and `pp/_philr.py` are each 100%.
- The property tests also passed under Hypothesis seeds 1, 2 and 3.
- The RED result in each "expect failure" step was reproduced by running the
  task's tests against the previous commit.
- APIs checked in the installed versions: mudata 0.4.1 (`MuData(mapping)`
  keeps a TreeData modality's type in memory and through `.copy()`;
  `write_h5mu` + `read_h5mu` returns it as AnnData without `vart`; the global
  `obs` has no columns until `pull_obs()`; sample sets that only overlap give
  a union with `obsmap` gaps; `mod_names` is deprecated with a
  `FutureWarning`), scikit-bio 0.7.4 (`skbio.stats.ordination.mmvec(X, Y,
  dimensions=3, optimizer="lbfgs", max_iter=1000, ..., seed=)` returns
  `MMvecResult` with `ranks`, a row-centred microbes x metabolites
  DataFrame; refuses sparse and AnnData input with `TypeError`; raises
  `ValueError: X contains all-zero columns` on an all-zero feature or
  sample; deterministic for a fixed seed with L-BFGS), scikit-learn 1.9.1
  (`sklearn.utils.Tags`, `validate_data(estimator, X, accept_sparse=,
  reset=)`, `check_non_negative(X, whom)`, `check_is_fitted`,
  `parametrize_with_checks`; `check_estimator` emits `SkipTestWarning`, a
  `UserWarning` subclass, so it cannot run under `-W error::UserWarning`).

**Research behind slices 4B-4C** (2026-10-08, notes in the session
scratchpad `p4-notes.md`):
- **MGM** [V]: PyPI `microformer-mgm` 0.5.8 (2025-02-17, MIT) ships the
  pretrained general model inside its 33 MB wheel:
  `mgm/resources/general_model/pytorch_model.bin` (35.7 MB, a GPT-2: 8
  layers, 8 heads, 256 dimensions, 512 positions, vocabulary 9,669),
  `phylogeny.csv` (9,665 genera with a mean and standard deviation each) and
  a pickled tokenizer whose vocabulary is `<pad>, <mask>, <bos>, <eos>`
  followed by `phylogeny.csv`'s genera in order (read with `pickletools`,
  never unpickled). The package pins `numpy==1.24.3`, `pandas==2.0.3`,
  `torch==2.0.1`, `scikit-learn==1.3.1`, `transformers==4.33.3`, so it
  cannot be installed beside biotapy. The weights are 101 tensors; the 100
  under `transformer.` load into `transformers.GPT2Model` 5.19.0 with no
  missing or unexpected keys (`torch.load(..., weights_only=True)`,
  `lm_head.weight` left out), and embed 4 samples in 0.3 s on CPU.
- **MGM2** [V model card]: Hugging Face `LudensZhang/MGM2`, no licence on the
  card; needs the 650M-parameter NTv3 nucleotide model, sequences and the
  repository's scripts. Not a light reference.
- **BiomeGPT**: a bioRxiv preprint; no public code or weights found.
  **Waypoint** (`waypoint-bio`): gated Hugging Face checkpoints.
- **torch** [V PyPI]: 2.14.1 (2026-09-30) is newest; cp314 wheels exist from
  2.9.0 (2025-10-15). The PyPI Linux wheel is 555 MB and pulls CUDA 13
  packages; `https://download.pytorch.org/whl/cpu` has CPU-only wheels for
  CPython 3.12-3.14. A scratch venv with torch 2.13.0+cpu, transformers
  5.19.0, numpy and pandas took 946 MB.
- **PyTorch Dataset** [V torch 2.13.0+cpu]: rows densified one at a time in
  `__getitem__` batch with the default `collate_fn`; per-row sparse COO
  tensors make the default collate raise `RuntimeError: Batches of sparse
  tensors are not currently supported by the default collate_fn`, and
  `torch.sparse_coo_tensor` emits a `UserWarning` about invariant checks.
- **Leakage measured** [V]: HMP2 taxa (130 samples, IBD vs non-IBD),
  5-fold stratified CV, logistic regression on CLR, ROC AUC. The
  prevalence filter outside vs inside the pipeline: 0.475 vs 0.478 (shuffled
  labels 0.622 vs 0.624), no measurable leak. A supervised selection step
  (`SelectKBest(f_classif, k=20)`) outside vs inside: 0.722 vs 0.571
  (shuffled labels 0.754 vs 0.550). Design note 9 uses this; the notebook's
  5 x 5 repeated-fold numbers are in slice 4D decision 38.

# Goal
Differentiator 3: ML-ready by default, and multi-omics as MuData.[^spec]

# Entry criteria
- Phase 3 exit gate passed and 0.3 released. (Done: v0.3.0 on PyPI,
  2026-10-08; Phase 3 closed in `2332a73`.)

# Design notes (resolved)

Each note answers one design question from the planning brief, with the
reason. Notes marked **(user)** change a contract, a rule, a dependency, CI or
a roadmap signature and are repeated under "Decisions for the user".

1. **Slices and order: 4A no new dependency, 4B the extra `torch`, 4C
   embedding plugins and MGM, 4D docs and release. (user)**
   - 4A holds everything the installed packages already do: `io.to_mudata`
     (4.1), `tl.mmvec` (4.2), the pseudocount helper's move to `_core`
     (4.A0) and the transformers (4.3). The brief put mmvec in a later slice
     on the assumption that it needs TensorFlow; scikit-bio 0.7.4 ships it
     natively ([V], design note 3), so it costs no dependency and joins 4A.
   - 4B adds the first heavy extra, `torch`, with `ml.to_torch` (4.5) and the
     CI job that installs it. It comes before the plugins because the
     reference plugin needs torch too, and the extra, its CPU wheel source and
     its CI job are one approval.
   - 4C adds the plugin interface `ml.embed` (4.4) and MGM, the reference
     plugin, behind its own extra `mgm` (torch + transformers) and a pooch
     download of the weights.
   - 4D writes the ML and multi-omics docs, the leak-free cross-validation
     notebook (4.6, exit gate 1), the end-to-end embedding page (exit gate
     2), the Coming-from-R check, knowledge (4.7) and the 0.4.0 release.
     The notebook needs only 4A, but it is an executed docs page, and docs
     pages come last as in Phase 3.

2. **`io.to_mudata(modalities: Mapping[str, AnnData]) -> MuData`; samples
   are intersected. (user: signature)**
   - The roadmap's `to_mudata(**modalities)` is replaced by one mapping
     argument: the data comes first (R3.1), there is no `**kwargs` (R3.7),
     and a function table's two modalities go in as
     `{**table.mod, "taxa": tdata}`, as the
     [function-tables-as-mudata](/decisions/function-tables-as-mudata.md)
     forward note requires. A MuData value raises `TypeError` naming that fix
     (MuData's modalities must be AnnData; it cannot nest).
   - "Aligning samples" means the intersection, in the first modality's
     order, with one `UserWarning` naming how many samples each modality
     lost; no shared sample raises. MuData's own default is the union with
     `obsmap` gaps [V], which no paired method can use.
   - Each modality is a copy (purity, R6.5); a TreeData stays a TreeData in
     memory [V]. Modality names are documented, not enforced, so a
     `proteins` modality is not refused.
   - No provenance entry: `add_provenance` creates `uns["biotapy"]` with
     `x_kind="counts"` when it is missing, which would mislabel a
     metabolite table.
   - **Saving.** `write_h5mu` drops a TreeData modality's tree [V]. Of the
     roadmap's two options, biotapy documents it (docstring, guide) and pins
     it with a test (`test_h5mu_drops_a_treedata_modality_tree`), which fails
     the day mudata or treedata keeps the tree; a biotapy writer would be a
     second file format for one gotcha. Filing an upstream issue is a GitHub
     action outside this repository, so it is asked as decision 4, not done.
   - `mudata.to_mudata` exists with another meaning (split one AnnData by a
     column). biotapy's is always reached as `bt.io.to_mudata`; the roadmap
     name is kept (decision 3 offers `io.combine`).

3. **`tl.mmvec(mdata, *, microbes="taxa", metabolites="metabolites",
   seed=None) -> pd.DataFrame` wraps scikit-bio. (user: route)**
   - scikit-bio 0.7.4 has `skbio.stats.ordination.mmvec` (added in 0.7.3):
     the multinomial model of Morton et al.[^mmvec] fitted by SciPy's L-BFGS-B, with
     a seeded start [V]. biocore/mmvec (TensorFlow 1, unmaintained) and a
     native port are not needed (R2.1).
   - The wrapper returns `MMvecResult.ranks`: microbes x metabolites log
     conditional probabilities, each row centred to sum 0 (scikit-bio's
     documented output and the quantity the mmvec paper ranks). The
     embeddings, `probs`, `predict` and `score` stay one scikit-bio call
     away; `Notes` says so. `dimensions`, `max_iter` and the Adam options
     are not exposed (R2.3: no test or use case yet).
   - Checks before the fit, each naming the argument: the modality exists
     (`KeyError`), the two hold the same samples in the same order
     (`ValueError` naming `io.to_mudata`), values are finite and
     non-negative, and no feature or sample is all zero (scikit-bio raises
     `X contains all-zero columns`, which names neither the modality nor the
     feature).
   - scikit-bio needs dense input [V]: both modalities are densified once,
     with the comment and the memory note (R6.2).
   - No golden test: mmvec has no R implementation (the original is Python).
     Tests pin the row centring (property), seed determinism, and that a
     metabolite built from one microbe ranks highest for it.
   - `tl` purity: it returns a table that no slot of one AnnData holds, so
     it has no `inplace`, like `tl.permanova`
     ([pure-by-default](/decisions/pure-by-default.md) Consequences).

4. **The pseudocount step moves to `_core/_composition.py` (4.A0).**
   `pp/_transform.py:pseudocounted` (with its validation) is needed by
   `ml.CLR`, a second subpackage, so it moves to `_core` (R4.3,
   module-boundaries rule 2). It takes `X` instead of an AnnData;
   `check_pseudocount` is split out because `CLR.fit` validates the
   pseudocount without densifying. A `refactor(core)` commit of its own, so
   the `feat(ml)` diff stays one feature; behaviour and messages are
   unchanged (every Phase 3 `pp` test passes untouched).

5. **`ml.PrevalenceFilter` and `ml.CLR`; `ml.RelativeAbundance` is dropped.
   (user)**
   - `PrevalenceFilter(SelectorMixin, BaseEstimator)`: `fit` learns
     `prevalence_` (the fraction of training samples in which each feature
     is non-zero) and keeps `prevalence_ >= min_prevalence`, the rule and the
     division of `pp.filter_features`. `SelectorMixin` supplies `transform`,
     `get_support`, `get_feature_names_out` and keeps CSR input sparse
     (R2.1). Nothing kept raises in `fit`, as scikit-learn's
     `VarianceThreshold` does. `min_prevalence` defaults to 0.1 because
     scikit-learn needs every parameter to have a default; `min_total` is
     not offered (R2.3).
   - `CLR(OneToOneFeatureMixin, TransformerMixin, BaseEstimator)`: stateless;
     `fit` checks `X` (`check_non_negative`) and the pseudocount; `transform`
     returns `skbio.stats.composition.clr` of the `_core` pseudocount step,
     equal to `pp.clr` (tested to 1e-12) with `pp.clr`'s pseudocount warning,
     ending in `(ml.CLR)`.
   - `RelativeAbundance` is `sklearn.preprocessing.Normalizer(norm="l1")`
     for non-negative data (each row divided by its sum, sparse kept, an
     all-zero row left at zero, as `pp.relative`). A biotapy class would
     duplicate it (R2.1); the guide names it.
   - R3.6: the classes inherit only scikit-learn's own bases and mixins,
     one level, no biotapy base class.
   - Parameters are validated in `fit`, never in `__init__`: scikit-learn's
     convention, which `check_estimator` tests
     (`check_no_attributes_set_in_init`, `check_dont_overwrite_parameters`).
   - `parametrize_with_checks` runs scikit-learn's estimator checks as
     pytest cases (47 per class; `check_array_api_input` skips without
     `SCIPY_ARRAY_API`).
     `check_estimator` itself emits `SkipTestWarning`, a `UserWarning`
     subclass, which the gate's `-W error::UserWarning` turns into a
     failure [V]. The checks feed random floats below 0.5, so `CLR`'s
     pseudocount warning is filtered in that one test, with the reason in a
     comment; the warning has its own test.
   - Tags: `input_tags.sparse = True` on both; `CLR` adds
     `positive_only = True` (so the checks require `fit` to refuse negative
     input) and `requires_fit = False`.
   - Typing: subclassing scikit-learn's bases passes `mypy --strict`; the
     call `super().__sklearn_tags__()` is reported as an untyped call because
     `untyped_calls_exclude = ["sklearn", ...]` does not reach calls through
     `super()` [V], so each `__sklearn_tags__` carries one
     `# type: ignore[no-untyped-call]` with a comment.
   - Docs: the default autosummary class page lists inherited scikit-learn
     methods whose docstrings fail `nitpicky` (`array-like`, `n_samples`,
     `MetadataRequest`) [V]. A template
     `docs/_templates/autosummary/class.rst` documents only the class's own
     members, and `intersphinx_mapping` gains scikit-learn for
     `:class:~sklearn.pipeline.Pipeline`.
   - They take arrays, sparse matrices and DataFrames, not AnnData:
     scikit-learn's protocol. `tdata.X` is the input.

6. **`ml` purity. (user)** [pure-by-default](/decisions/pure-by-default.md)
   gains a row: `ml` estimators follow scikit-learn (`fit` stores what it
   learns on the estimator and returns it, `transform` returns a new array,
   the data is never changed). Slice 4C adds `ml.embed` to the `tl` row's
   convention (returns the embedding; `inplace=True` writes
   `obsm["X_<model>"]`), and slice 4B `ml.to_torch`, which returns a new
   object that references `X` without copying it.
   [function-shape](/contracts/function-shape.md) gains one sentence: `ml`'s
   transformers are classes whose options are constructor arguments that
   `fit` validates, and the class docstring carries the skeleton
   (`tests/test_docstrings.py` already reads class docstrings; it gains
   `ml` in its module set).

7. **`ml.to_torch(adata, *, label_key=None, layer=None) -> Dataset` (slice
   4B). (user: extra, CI)**
   - Map-style dataset whose `__getitem__(i)` densifies row `i` of the CSR
     table to a float32 tensor: the full table is never dense (R6.2) and the
     default `collate_fn` stacks the rows [V]. Per-row sparse tensors would
     need a custom `collate_fn` and trip `-W error::UserWarning` [V], so they
     are not used.
   - `label_key` names an `obs` column: categories become int64 codes in
     category order, numbers float32; a missing value raises. Without it,
     items are bare tensors.
   - The return annotation is the bare `Dataset` (slice 4B design).
   - torch must not be imported at module level, and the
     `import-without-extras` CI job imports every module, so the `Dataset`
     subclass is defined inside the function after
     `import_optional("torch", extra="torch")` [V: mypy does not follow torch
     (`follow_imports = "skip"`), so one `# type: ignore[misc]` holds with and
     without torch; slice 4B design].
   - The extra is `torch = ["torch>=2.9"]`: 2.9.0 is the first release with
     CPython 3.14 wheels [V]. CI installs CPU wheels through a uv index
     (`[[tool.uv.index]] pytorch-cpu`, `explicit = true`, with
     `[tool.uv.sources] torch`), so the job does not pull CUDA
     [V: `uv lock` resolves `2.14.1+cpu` for Linux and Windows and `2.14.1`
     for macOS arm64, py3.12-3.14; Intel macOS has no wheel; slice 4B design]. A new marker `torch` is excluded from the default
     run like `r`, and a CI job `ml-extras` (Linux, Python 3.13) runs
     `-m torch` and joins `check.needs`.

8. **Embedding plugins: `ml.embed(adata, model, *, inplace=False)`;
   MGM is the reference plugin (slice 4C).
   (user: interface, reference plugin, extra)**
   - Discovery: `importlib.metadata.entry_points(group="biotapy.embeddings")`
     at call time; an unknown `model` raises `KeyError` listing the installed
     names. A plugin entry point loads a callable
     `embed(adata: AnnData) -> np.ndarray`.
   - biotapy checks what the plugin returns: a 2-D finite float array with
     one row per sample, else `ValueError` naming the plugin. The result is
     returned, or with `inplace=True` written to `obsm[f"X_{model}"]`
     (data-model-slots already reserves `X_<plugin>`).
   - **Reference plugin: MGM**,[^mgm] the only microbiome foundation model found
     with public weights, a licence (MIT) and a size (36 MB) a test can
     download. It ships in biotapy as `ml/_mgm.py`, registered by biotapy's
     own `[project.entry-points."biotapy.embeddings"] mgm =
     "biotapy.ml._mgm:embed"`, so it exercises the same discovery path as a
     third-party plugin. The extra `mgm = ["torch>=2.9", "transformers>=5"]`;
     `microformer-mgm` itself is not a dependency (its pins conflict). The
     weights, `phylogeny.csv` and config are fetched by pooch from the
     0.5.8 wheel on files.pythonhosted.org, pinned by sha256
     (`210891685565022ea869e88a7452769dc6fbe3f699d2a133e15338d1e30eb92b`),
     with `pooch.Unzip(members=...)`.
   - MGM's input is genus abundances: the plugin reads `var["genus"]`
     (after `pp.tax_glom(tdata, "genus")`), maps each genus to the token
     `g__<genus>`, reads `var["genus"]` with MGM's regex, sums features of one
     genus, leaves out features without a known genus with one warning
     counting them, and repeats MGM's preprocessing: relative
     abundance, z-score per genus with `phylogeny.csv`'s mean and standard
     deviation, keep genera above the z-score of zero, sort descending, wrap
     in `<bos>` ... `<eos>`, pad or cut to 512 [V, read from
     `mgm/src/MicroCorpus.py`]. The sample embedding is the mean of the last
     hidden layer over the sample's tokens, as MGM's paper describes
     (Methods 4.5) and its notebook uses for sample embeddings (decision 28).
   - Parity: `microformer-mgm` 0.5.8 installs in an isolated Python 3.11
     environment; `tests/data/mgm/` holds its own embeddings, matched to 1.7e-6
     (decision 33).
   - Alternative (decision 9): a separate `biotapy-mgm` package with its own
     release cadence; it keeps transformers out of biotapy's extras.

9. **The leak-free notebook shows a leak that is real (4D, outline).
   (user: content)** The roadmap's `Pipeline(PrevalenceFilter, CLR,
   classifier)` beside the leaky version would show almost no difference:
   measured on HMP2, 0.475 vs 0.478 AUC [V], because a prevalence filter
   never sees the labels. The notebook (HMP2 taxa, IBD vs non-IBD, already
   downloaded by the docs job for the function tutorial) shows that result
   honestly, then adds a supervised selection step (`SelectKBest`) outside
   and inside the pipeline: 0.722 vs 0.571 on real labels and 0.754 vs
   0.550 on shuffled labels [V] (single 5-fold, one shuffle). The notebook uses 5 x 5 repeated folds and 10 shuffles:
   0.527 vs 0.526, 0.558 vs 0.705, 0.544 vs 0.767 (slice 4D decision 38), where the leak manufactures signal from
   noise. It is executed in the docs build (exit gate 1).

10. **Exit gate 2, "one foundation model plugged in end to end", is proven
    in CI, not in the docs build. (user)** The docs build has no torch; adding
    torch and transformers to the `doc` group would add about a gigabyte to
    every docs build [V: 946 MB venv]. A test marked `mgm` in the
    `ml-extras` job downloads the MGM weights, embeds GlobalPatterns' genera
    and checks shape, determinism and parity with MGM's own embeddings
    (`tests/data/mgm`); the
    docs page shows the same code as a non-executed block and quotes the
    test's numbers, the Phase 3 pattern for the R bridges.

11. **Python 3.14.** torch >= 2.9 has cp314 wheels for Linux, macOS and
    Windows [V]; transformers is pure Python; scikit-learn, mudata and
    scikit-bio already run there in the existing matrix. The `ml-extras`
    job pins 3.13, like `r-bridge` [V: tokenizers, safetensors and
    hf-xet ship abi3 wheels and regex cp314 wheels; MGM ran on CPython 3.14.8
    with torch 2.14.1+cpu and transformers 5.19.0].

# Global constraints
- Python >= 3.12. No new runtime dependency in Phase 4. Extras only after the
  user approves them at the start of the task named below.
- `X` stays CSR. Densify once, only where the delegated library needs dense
  input (scikit-bio `mmvec` and `clr`; one row at a time in `to_torch`),
  with a comment and the memory cost in `Notes` (R6.2).
- Results go only to `obsm["X_<model>"]` (with `inplace=True`); `tl.mmvec`
  and the transformers write no slot (R3.3, R6.3).
- torch, transformers and the MGM weights are reached only inside functions
  through `import_optional` and pooch (R4.6); `import biotapy` and the
  `import-without-extras` CI job stay green without them.
- Randomness: `tl.mmvec` and anything stochastic in `ml` take `seed` and
  convert it once with `as_generator` (R3.4); `to_torch` has no randomness
  (shuffling is `DataLoader`'s).
- Classes only in `ml`, only on scikit-learn's or torch's bases (R3.6).
- Commits stage explicit paths only. Never stage `.claude/`, `.superpowers/`,
  `.worktrees/`, `notebooks/` or `build/`. Every task's last commit also
  stages `.knowledge/roadmap/phase-4-ml-multiomics.md` with that task's box
  ticked and `.knowledge/log.md` with its dated line (R12.4).
- Every pytest, sphinx or Python run exports `BIOTAPY_DATA_DIR` to a
  scratch pooch cache; `~/.cache/biotapy` must not appear.
- prek's `--all-files` only sees tracked files: stage new files before
  running it.

# Dependencies to approve (ask at the start of the task named)
| Task | Group | Package | Reason |
|---|---|---|---|
| 4.3 | docs config only | intersphinx entry for scikit-learn; `docs/_templates/autosummary/class.rst` | link `Pipeline`; keep inherited scikit-learn docstrings out of `nitpicky` |
| 4.B0 | extra `torch` | `torch>=2.9` (BSD-3-Clause and bundled licences) | `ml.to_torch`; the roadmap's planned extra |
| 4.B0 | uv config | index `https://download.pytorch.org/whl/cpu` for torch (`explicit = true`) | CI and contributors get a 192 MB CPU wheel, not 555 MB plus CUDA 13 packages |
| 4.B1 | CI | job `ml-extras` (Linux, Python 3.13, `-m torch`), in `check.needs` | run what the default jobs skip |
| 4.5 | docs and type-check config only | intersphinx entry for PyTorch; mypy override `follow_imports = "skip"` for torch; root `conftest.py` doctest-marker hook | link `Dataset`; one mypy answer with or without torch; run the example only where torch is |
| 4.C0 | none | `_core.make_pooch` | the cache directory and variable defined once (R4.3) |
| 4.C1 | extra `mgm` | `torch>=2.9`, `transformers>=5` (Apache-2.0) | the MGM reference plugin |
| 4.4b | data, run time | MGM 0.5.8 weights from files.pythonhosted.org via pooch (MIT, 33 MB download, 67 MB in the cache with its extracted files) | never bundled (R6.6) |
| 4.6c | data, run time | mmvec's `examples/soils` BIOM files from GitHub via pooch (BSD-3-Clause repository, 135 KB, commit `88ca33b`, SHA-256 pinned) | the multi-omics tutorial (decision 37) |
| - | none | `microformer-mgm`, TensorFlow, biocore `mmvec`; `huggingface_hub` | not imported by biotapy (transformers imports `huggingface_hub`; nothing contacts the Hub) (design notes 3, 8) |

# Review focus
The five ways real users are most likely to get a wrong answer from Phase 4
without an error. Each line names the test that pins it.

1. **Preprocessing fitted on every sample before cross-validation.**
   Expected: the filter learns from the training rows only and is refitted
   per fold. Tests: 4.3 `test_prevalence_filter_learns_only_from_the_samples_it_is_fitted_on`,
   `test_a_pipeline_refits_the_filter_in_every_fold`; 4D's notebook shows
   the size of the leak a supervised step causes
   (`test_the_leak_free_cv_tutorial_quotes_its_numbers`, and the network test that runs
   the page).
2. **Modalities paired on different samples.** Expected: `io.to_mudata`
   keeps the shared samples in one order and warns; `tl.mmvec` refuses
   unaligned modalities. Tests: 4.1
   `test_keeps_the_samples_every_modality_has_in_the_first_ones_order`,
   `test_every_modality_holds_the_same_samples`; 4.2
   `test_unaligned_samples_raise`.
3. **The pseudocount on the wrong scale inside a pipeline.** Expected: the
   `pp.clr`'s warning, ending in `(ml.CLR)`. Test: 4.3
   `test_clr_pseudocount_above_the_smallest_value_warns`.
4. **A tree lost on save.** Expected: documented, and pinned so an upstream
   change is noticed. Test: 4.1 `test_h5mu_drops_a_treedata_modality_tree`.
5. **An embedding whose rows are not the samples.** Expected: a plugin's
   output with the wrong number of rows or non-finite values raises naming
   the plugin. Tests:
   4.4 `test_plugin_with_wrong_row_count_raises`, `test_plugin_with_nan_raises`,
   `test_a_refused_result_is_not_written`.

# Slices
| Slice | Delivers | Tasks | Ends with |
|---|---|---|---|
| **4A - No new dependency** | multi-omics MuData, mmvec, leak-free transformers | 4.1 `io.to_mudata` · 4.2 `tl.mmvec` · 4.A0 `refactor(core)` pseudocount step · 4.3 `ml.PrevalenceFilter`, `ml.CLR` | Checkpoint A |
| **4B - torch** | the extra `torch`, `ml.to_torch`, CPU wheels in CI | 4.B0 `build` extra `torch` and its CPU index · 4.5 `ml.to_torch` (+ marker, mypy override) · 4.B1 CI job `ml-extras` | Checkpoint B |
| **4C - Embeddings** | the plugin interface and MGM | 4.C0 `refactor(core)` `make_pooch` · 4.C1 `build` extra `mgm` · 4.4 `ml.embed` and the entry-point group (+ marker `mgm`) · 4.4b MGM plugin (+ reference embeddings, `ml-extras` step) | Checkpoint C |
| **4D - Docs and release** | the guide pages, the two exit-gate pages, 0.4 | 4.F4 · 4.F5 `fix(core)` sums · 4.6 leak-free CV notebook · 4.6b embedding page · 4.6c `datasets.biocrust` and multi-omics tutorial · 4.D1 Coming-from-R check · 4.7 knowledge · 4.D2 release 0.4.0 | exit gate |

Execution order inside 4A: **4.1 -> 4.2 -> 4.A0 -> 4.3 -> Checkpoint A.**
`tl.mmvec`'s tests build their MuData with `io.to_mudata`; the `refactor`
lands before the `feat(ml)` that needs it.

Execution order inside 4C: **4.C0 -> 4.C1 -> 4.4 -> 4.4b -> Checkpoint C.**
The cache refactor and the extra land before the features that need them.

Execution order inside 4D: **4.F4 -> 4.F5 -> 4.6 -> 4.6b -> 4.6c -> 4.D1 -> Checkpoint D
(with 4.7) -> 4.D2.**

# Tasks (checklist)
- [x] 4.1 `io.to_mudata(modalities) -> MuData` and the multi-omics decision
- [x] 4.2 `tl.mmvec(mdata, *, microbes="taxa", metabolites="metabolites", seed=None) -> pd.DataFrame`
- [x] 4.A0 `refactor(core)`: the pseudocount step moves to `_core/_composition.py`
- [x] 4.3 `ml.PrevalenceFilter(min_prevalence=0.1)`, `ml.CLR(pseudocount=0.5)`; scikit-learn's estimator checks pass
- [x] Checkpoint A
- [x] 4.F1 `da.ancombc2`: pin scikit-bio's bias E-M underflow; the schema property keeps to fittable designs
- [x] 4.F2 `pp.filter_features` names a wrongly typed threshold
- [x] 4.F3 `refactor(core)`: one finite, non-negative check
- [x] 4.B0 `build`: the extra `torch = ["torch>=2.9"]`, installed by uv from PyTorch's CPU index
- [x] 4.5 `ml.to_torch(adata, *, label_key=None, layer=None) -> torch.utils.data.Dataset` and the marker `torch`
- [x] 4.B1 CI job `ml-extras` for `-m torch` tests
- [x] Checkpoint B (PR #29 merged as `fd6a8be`)
- [x] 4.C0 `refactor(core)`: one download cache, `_core.make_pooch`
- [x] 4.C1 `build`: the extra `mgm = ["torch>=2.9", "transformers>=5"]`
- [x] 4.4 `ml.embed(adata, model, *, inplace=False)` and the entry-point group `biotapy.embeddings`
- [x] 4.4b MGM reference plugin, its reference embeddings and the marker `mgm`
- [x] Checkpoint C (PR #30 merged as `28de1ad`; the user approved slice 4C on 2026-10-09)
- [x] 4.F4 `_core.sum_by` sums narrow integer and bool tables in NumPy's sum dtype, so `pp.tax_glom` cannot wrap
- [x] 4.F5 `_core.sum_pairs` does the same, for `fn.func_glom` (decision 36)
- [x] 4.6 Leak-free cross-validation notebook (exit gate 1)
- [x] 4.6b End-to-end embedding page (exit gate 2's docs)
- [x] 4.6c `bt.datasets.biocrust()` and the multi-omics tutorial (decision 37)
- [x] 4.D1 Coming-from-R check
- [x] Checkpoint D
- [x] 4.7 Knowledge: `ml` and `tl` gotchas, `embedding-plugins` stable, verification bumps
- [x] 4.D2 Release 0.4.0

# Exit gate
- [x] Leak-free CV example executed in docs: `docs/tutorials/leak_free_cv.md`
  runs in the docs CI job (design note 9). Proven on PR #31: docs job 114304008718
  executed `tutorials/leak_free_cv.md`.
- [x] One foundation model plugged in end to end: `ml-extras`' `-m mgm` step runs
  `tests/ml/test_mgm.py::test_embeds_global_patterns_end_to_end`: GlobalPatterns'
  genera through `bt.ml.embed(..., "mgm", inplace=True)`, 26 x 256, equal
  across calls, every sample's nearest neighbour from its own environment
  (design note 10). Proven on PR #30 (Test run `37962120924`): `ml-extras`'
  `-m mgm` step `13 passed`.
- [x] All Phase 1-3 gates still green. Proven on PR #31: every job of Test run 38083140914 green.

# Risks
- **MGM's pooling is the paper's description, not code**: the paper names
  mean pooling for the pretrained model but no layer (decision 28); parity is
  pinned against MGM's own code.
- **oneMKL's VML CPU-type cache in torch's CPU wheels** races on a process's
  first parallel call (pytorch/pytorch#188792; 4 threads hit it too, rarely)
  -> biotapy fills the cache with one serial call and tests every call at
  1e-5 (decision 33).
- **The weights URL is a PyPI file** -> files.pythonhosted.org URLs are
  immutable and pinned by hash; if the release is yanked the download still
  works (yanking does not delete files), but a deletion would break the
  plugin, so the registry entry names the version.
- **torch install size and platform wheels** -> CPU index in CI; users install
  what suits their hardware; the docs say so. torch has no Intel macOS wheel,
  and 2.14's macOS wheel needs macOS 14 (`uv pip compile` for an older macOS
  target falls back to 2.11.0); every `uv run` CI job reads
  download.pytorch.org while locking.
- **mypy without torch** -> resolved in 4B: torch is not followed, so the
  result does not depend on the environment.
- **mmvec convergence is silent**: scikit-bio's L-BFGS prints its status
  only with `verbose=True` and `MMvecResult` has no flag [V] -> `Notes`
  points to `res.convergence` through scikit-bio; reconsider a warning if a
  user reports a non-converged fit.
- **scikit-learn's estimator checks change between releases** -> the check
  suite runs in every test job; a new check that fails is a real
  incompatibility, fixed in biotapy rather than skipped (R11.5).
- **Leakage misread as "the prevalence filter leaks a lot"** -> the
  notebook prints the measured near-zero difference first (design note 9), pinned by
  `test_the_leak_free_cv_tutorial_quotes_its_numbers`.
- **Tutorial numbers move with scikit-learn, torch or a dataset** -> the page-running `network`/`mgm`
  tests fail, and the constants and the page change together (decision 39).

---
## Slice 4A - No new dependency

**Goal:** with only what biotapy already installs, a user combines taxa,
a function table and metabolites into one MuData over the samples they
share, ranks metabolites by how strongly each microbe predicts them, and
runs a scikit-learn pipeline whose prevalence filter and CLR are fitted
inside each cross-validation fold.

### Slice 4A design
- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `io/_mudata.py` | `to_mudata` and `_check_modality` (4.1) |
  | `tl/_mmvec.py` | `mmvec` and `_table` (4.2) |
  | `_core/_composition.py` | `check_pseudocount`, `pseudocounted`, moved from `pp/_transform.py` (4.A0) |
  | `ml/__init__.py`, `ml/_transformers.py` | the new subpackage: `PrevalenceFilter`, `CLR` (4.3) |
  | `docs/guide/multiomics.md`, `docs/guide/machine_learning.md` | the two new guide pages |
  | `docs/_templates/autosummary/class.rst` | own-members-only class pages (4.3) |
  | `.knowledge/decisions/multiomics-as-mudata.md` | the 4.1 decision (`draft` until the user confirms it) |

  `ml` is created in 4.3, the task that fills it (R4.8). Its layer is
  already in `pyproject.toml`'s import-linter contract
  (`(pl) | (ml) | (da)`), so no config change is needed for it.
- **Facts the tasks rely on** (measured on the prototype; re-check each,
  R2.2): see "How slice 4A was checked" above, and:
  - `toy()`'s smallest stored count is 1, so `CLR()` on it does not warn;
    `[[0.2, 0.8], [0.5, 0.5]]` warns with "(0.2)".
  - In `toy()`, `f6` (Prevotella) is high in s4-s6 and `f3`
    (Faecalibacterium) in s1-s3; a metabolite equal to `f6 + 1` ranks
    first in `f6`'s mmvec row, and one equal to `f3 + 1` in `f3`'s, with
    `seed=0`.
  - `bt.pp.filter_features(toy(), min_prevalence=1.0)` keeps `f3, f4`;
    `min_prevalence=0.8` keeps the same features as `PrevalenceFilter(min_prevalence=0.8)`.
  - Building an AnnData with repeated `obs_names` warns "Observation names
    are not unique" (anndata 0.13.4); the 4.1 test expects it.
  - `toy_humann()`'s MuData has modalities `function`,
    `function_by_taxon`.

### Slice 4A global constraints (in addition to the Phase 4 list)
- Tests reach functions through `bt.io.to_mudata`, `bt.tl.mmvec`,
  `bt.ml.*`. `_core` unit tests (`tests/core/test_composition.py`) and test
  fixtures built with `_core` (`TreeData`, `get_tree`, `tree_tips`) are the
  exceptions, as in Phases 1-3.
- Run commands: `uv run --group test pytest <path> -q`. Gate before every
  commit: `uvx prek run --all-files` (new files staged first); with docs
  changes also `uv run --group doc sphinx-build -W -b html docs
  docs/_build/html` after deleting `docs/_build` and `docs/generated`.
- No `mod_names` (deprecated in mudata 0.4.1, a `FutureWarning`): use
  `list(mdata.mod)`.

### Slice 4A review focus
The Phase 4 review focus items 1-4 are slice 4A's. In addition:
- **Purity:** `to_mudata` copies every modality; `mmvec` leaves both
  modalities unchanged; `CLR` leaves its input unchanged. Tests: 4.1
  `test_returns_copies_and_keeps_the_input`, 4.2 `test_keeps_the_input`,
  4.3 `test_clr_keeps_its_input`.
- **Equality with the `pp` verbs:** `PrevalenceFilter` keeps what
  `pp.filter_features` keeps; `CLR` equals `pp.clr`. Tests: 4.3
  `test_prevalence_filter_keeps_what_pp_filter_features_keeps`,
  `test_clr_equals_pp_clr`.
- **The refactor changes nothing:** every `tests/pp` test passes unedited
  after 4.A0.

---

### Task 4.1: `io.to_mudata` and the multi-omics decision

**Files:** create `src/biotapy/io/_mudata.py`, `tests/io/test_mudata.py`,
`docs/guide/multiomics.md`, `.knowledge/decisions/multiomics-as-mudata.md`;
modify `src/biotapy/io/__init__.py`, `docs/guide/index.md`, `docs/api.md`,
`.knowledge/decisions/index.md`,
`.knowledge/decisions/function-tables-as-mudata.md`,
`.knowledge/contracts/data-model-slots.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `datasets/_hmp2.py` (it builds its MuData with `MuData(...,
obs=)` and `push_obs()`, already follows the names, and needs no change);
`_core/_function.py`; any reader.
**Interfaces:**
- Consumes: `_core.warn_user`; mudata's `MuData(mapping)`.
- Produces: `bt.io.to_mudata(modalities: Mapping[str, AnnData]) -> MuData`.

- [x] **Step 1: Failing tests.** Create `tests/io/test_mudata.py`:
  ```python
  from contextlib import nullcontext

  import mudata
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from anndata import AnnData
  from hypothesis import given
  from hypothesis import strategies as st
  from mudata import MuData

  import biotapy as bt


  def _metabolites(samples: list[str]) -> AnnData:
      values = np.arange(len(samples) * 3, dtype=np.float64).reshape(len(samples), 3)
      return AnnData(
          X=sp.csr_matrix(values),
          obs=pd.DataFrame({"batch": ["a"] * len(samples)}, index=samples),
          var=pd.DataFrame(index=["m1", "m2", "m3"]),
      )


  def test_aligned_modalities_become_one_mudata():
      tdata = bt.datasets.toy()
      mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata.obs_names.tolist())})
      assert isinstance(mdata, MuData)
      assert list(mdata.mod) == ["taxa", "metabolites"]
      assert mdata.obs_names.tolist() == tdata.obs_names.tolist()
      assert mdata["metabolites"].obs_names.equals(mdata["taxa"].obs_names)


  def test_keeps_the_samples_every_modality_has_in_the_first_ones_order():
      tdata = bt.datasets.toy()
      with pytest.warns(UserWarning, match=r"dropped: taxa 3 of 6, metabolites 1 of 4"):
          mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(["s4", "s3", "s2", "s9"])})
      assert mdata.obs_names.tolist() == ["s2", "s3", "s4"]
      for name in ("taxa", "metabolites"):
          assert mdata[name].obs_names.tolist() == ["s2", "s3", "s4"]
      np.testing.assert_array_equal(mdata["metabolites"].X.toarray()[:, 0], [6.0, 3.0, 0.0])


  def test_keeps_a_treedata_modality_and_its_tree():
      tdata = bt.datasets.toy()
      mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata.obs_names.tolist())})
      taxa = mdata["taxa"]
      assert type(taxa) is bt._core.TreeData
      assert set(bt._core.tree_tips(bt._core.get_tree(taxa))) == set(taxa.var_names)


  def test_function_table_modalities_go_side_by_side():
      table = bt.datasets.toy_humann()
      mdata = bt.io.to_mudata({**table.mod, "taxa": bt.datasets.toy()})
      assert list(mdata.mod) == ["function", "function_by_taxon", "taxa"]


  def test_returns_copies_and_keeps_the_input(assert_unchanged):
      tdata = bt.datasets.toy()
      metabolites = _metabolites(["s1", "s2"])
      before = tdata.copy()
      with pytest.warns(UserWarning, match=r"dropped: taxa 4 of 6$"):
          mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
      assert_unchanged(before, tdata)
      assert not mdata["taxa"].is_view and mdata["taxa"] is not tdata
      mdata["metabolites"].X[0, 1] = 99.0
      assert metabolites.X[0, 1] == 1.0


  def test_h5mu_drops_a_treedata_modality_tree(tmp_path):
      # The documented reason to save a tree-bearing modality with write_h5td too (contracts/tree-access).
      mdata = bt.io.to_mudata({"taxa": bt.datasets.toy()})
      mdata.write_h5mu(tmp_path / "study.h5mu")
      assert type(mudata.read_h5mu(tmp_path / "study.h5mu")["taxa"]) is AnnData


  def test_single_modality():
      mdata = bt.io.to_mudata({"taxa": bt.datasets.toy()})
      assert list(mdata.mod) == ["taxa"] and mdata.n_obs == 6


  def test_no_shared_sample_raises():
      with pytest.raises(ValueError, match="modalities share no sample"):
          bt.io.to_mudata({"taxa": bt.datasets.toy(), "metabolites": _metabolites(["x1", "x2"])})


  def test_empty_mapping_raises():
      with pytest.raises(ValueError, match="modalities is empty"):
          bt.io.to_mudata({})


  def test_a_mudata_entry_raises_naming_the_fix():
      with pytest.raises(TypeError, match=r"modalities\['function'\] must be an AnnData, got MuData"):
          bt.io.to_mudata({"function": bt.datasets.toy_humann()})


  def test_repeated_sample_ids_raise():
      with pytest.warns(UserWarning, match="Observation names are not unique"):
          metabolites = _metabolites(["s1", "s1"])
      with pytest.raises(ValueError, match=r"modalities\['metabolites'\] repeats sample ids: \['s1'\]"):
          bt.io.to_mudata({"taxa": bt.datasets.toy(), "metabolites": metabolites})


  @pytest.mark.parametrize("name", ["", 3])
  def test_bad_modality_name_raises(name):
      with pytest.raises(TypeError, match="modality names must be non-empty strings"):
          bt.io.to_mudata({name: bt.datasets.toy()})


  @given(st.lists(st.sampled_from([f"s{i}" for i in range(1, 9)]), min_size=1, max_size=8, unique=True))
  def test_every_modality_holds_the_same_samples(samples):
      tdata = bt.datasets.toy()
      shared = [name for name in tdata.obs_names if name in samples]
      metabolites = _metabolites(samples)
      if not shared:
          with pytest.raises(ValueError, match="share no sample"):
              bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
          return
      with (
          pytest.warns(UserWarning, match="samples missing from another modality are dropped: ")
          if len(shared) < max(6, len(samples))
          else nullcontext()
      ):
          mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
      assert mdata.obs_names.tolist() == shared
      assert mdata["metabolites"].obs_names.tolist() == shared
  ```
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/io/test_mudata.py -q` -> `14 failed`
  (`AttributeError: module 'biotapy.io' has no attribute 'to_mudata'`).
- [x] **Step 3: Implement.** Create `src/biotapy/io/_mudata.py`:
  ```python
  """Several data types over the same samples as one MuData (decisions/multiomics-as-mudata)."""

  from collections.abc import Mapping

  import pandas as pd
  from anndata import AnnData
  from mudata import MuData

  from biotapy._core import warn_user


  def to_mudata(modalities: Mapping[str, AnnData]) -> MuData:
      """Combine data types measured on the same samples into one MuData.

      Parameters
      ----------
      modalities
          Modality name -> samples x features table. biotapy's names are
          ``"taxa"``, ``"function"``, ``"function_by_taxon"``, ``"metabolites"``
          and ``"host"``; a TreeData stays a TreeData.

      Returns
      -------
      MuData
          One modality per entry, in the mapping's order, each a copy holding only
          the samples every modality has, in the first modality's order. The
          global ``obs`` has no columns; ``mdata.pull_obs()`` gathers them.

      Raises
      ------
      TypeError
          A name is not a non-empty string, or a value is not an AnnData (such as
          a function table, which is already a MuData).
      ValueError
          ``modalities`` is empty, a modality repeats a sample id, or the
          modalities share no sample.

      Warns
      -----
      UserWarning
          Some samples are missing from at least one modality; the warning names
          how many each modality loses.

      Notes
      -----
      R equivalent: ``MultiAssayExperiment::MultiAssayExperiment``, ``MultiAssayExperiment::intersectColumns``
      Guide: :doc:`/guide/multiomics`

      A function table from ``bt.io.read_humann`` or ``bt.io.read_picrust2`` is
      a MuData of two modalities; pass them as two entries, ``{**table.mod,
      "taxa": tdata}``.

      ``MuData.write_h5mu`` does not keep a TreeData's tree: the modality reads
      back as an AnnData without ``vart``. Save a tree-bearing modality with
      ``TreeData.write_h5td`` as well.

      Examples
      --------
      >>> import biotapy as bt
      >>> table = bt.datasets.toy_humann()
      >>> mdata = bt.io.to_mudata({**table.mod, "taxa": bt.datasets.toy()})
      >>> list(mdata.mod)
      ['function', 'function_by_taxon', 'taxa']
      """
      if not modalities:
          msg = "modalities is empty; pass at least one table, e.g. {'taxa': tdata}"
          raise ValueError(msg)
      shared: pd.Index[str] | None = None
      for name, mod in modalities.items():
          _check_modality(name, mod)
          shared = mod.obs_names if shared is None else shared.intersection(mod.obs_names, sort=False)
      if shared is None or shared.empty:
          firsts = {name: mod.obs_names[:3].tolist() for name, mod in modalities.items()}
          msg = f"modalities share no sample; their first sample ids: {firsts}"
          raise ValueError(msg)
      lost = [
          f"{name} {mod.n_obs - len(shared)} of {mod.n_obs}"
          for name, mod in modalities.items()
          if mod.n_obs > len(shared)
      ]
      if lost:
          warn_user(f"samples missing from another modality are dropped: {', '.join(lost)}")
      return MuData({name: mod[shared].copy() for name, mod in modalities.items()})


  def _check_modality(name: object, mod: object) -> None:
      """Raise unless ``name`` is a non-empty string and ``mod`` an AnnData with unique sample ids."""
      if not isinstance(name, str) or not name:
          msg = f"modality names must be non-empty strings, got {name!r}"
          raise TypeError(msg)
      if not isinstance(mod, AnnData):
          msg = (
              f"modalities[{name!r}] must be an AnnData, got {type(mod).__name__}; "
              "pass a function table's modalities as two entries: {**table.mod, ...}"
          )
          raise TypeError(msg)
      repeated = mod.obs_names[mod.obs_names.duplicated()].unique().tolist()
      if repeated:
          msg = f"modalities[{name!r}] repeats sample ids: {repeated[:5]}"
          raise ValueError(msg)
  ```
  `src/biotapy/io/__init__.py`:
  ```python
  from ._biom import read_biom, write_biom
  from ._dada2 import read_dada2
  from ._humann import read_humann
  from ._metaphlan import read_metaphlan
  from ._mudata import to_mudata
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
      "to_mudata",
      "write_biom",
  ]
  ```
- [x] **Step 4: Run, expect pass** -
  `uv run --group test pytest tests/io/test_mudata.py src/biotapy/io/_mudata.py -q -W error::UserWarning`
  -> `15 passed`; the property test also under `--hypothesis-seed=1`, `2`,
  `3`.
- [x] **Step 5: Docs.** Create `docs/guide/multiomics.md`:
  ````markdown
  # Multi-omics

  When several data types are measured on the same samples - taxa, a function
  table, metabolites, host data - biotapy keeps them in one
  [MuData](https://mudata.scverse.org/): one AnnData per data type, called a
  modality, over shared samples.

  ## Modality names

  | Modality | Holds |
  |---|---|
  | `"taxa"` | taxa, usually a TreeData with its phylogeny |
  | `"function"` | a function table's community rows |
  | `"function_by_taxon"` | a function table's per-taxon rows |
  | `"metabolites"` | metabolite intensities or counts |
  | `"host"` | host measurements, such as transcripts or clinical values |

  biotapy functions that read a modality take its name as an argument and
  default to these names. Other names work too.

  ## Building one

  `bt.io.to_mudata` takes a mapping of name to AnnData:

  ```python
  import biotapy as bt

  table = bt.io.read_humann("pathabundance.tsv")  # already a MuData
  mdata = bt.io.to_mudata({"taxa": tdata, **table.mod, "metabolites": metabolites})
  ```

  A function table is already a MuData of two modalities, and a MuData cannot
  hold another MuData, so `**table.mod` adds `function` and `function_by_taxon`
  as two entries.

  Every modality keeps only the samples that all of them have, in the first
  modality's order, and biotapy warns with how many samples each modality lost.
  Each modality is a copy, so the tables you passed are unchanged.

  The global `mdata.obs` starts without columns; `mdata.pull_obs()` gathers
  the modalities' sample metadata into it.

  ## Saving

  `mdata.write_h5mu("study.h5mu")` saves every modality, but a TreeData
  modality reads back as a plain AnnData without its tree. Save that modality
  with its tree as well:

  ```python
  mdata.write_h5mu("study.h5mu")
  mdata["taxa"].write_h5td("taxa.h5td")
  ```
  ````
  In `docs/guide/index.md` add `multiomics` after `function` in the
  toctree. In `docs/api.md` add `io.to_mudata` after `io.read_qiime2`.
- [x] **Step 6: Knowledge.** Create
  `.knowledge/decisions/multiomics-as-mudata.md` (`generated.at` is the
  commit time):
  ````markdown
  ---
  type: Decision
  title: Multi-omics data is one MuData with fixed modality names
  description: Several data types over the same samples are one MuData whose modalities are named taxa, function, function_by_taxon, metabolites and host; io.to_mudata keeps only the samples every modality has; a TreeData modality loses its tree in h5mu.
  tags: [io, mudata, multiomics]
  status: draft
  paths: ["src/biotapy/io/_mudata.py"]
  generated: { by: claude-code/claude-opus-5-5, at: 2026-10-08T17:54:32Z }
  commit: 872c16d
  sources:
    - id: spec
      resource: ../../plan.md
      title: Python Microbiome Toolkit development report
      author: human:pedrocr83
  ---

  # Context
  The spec stores several data types as MuData modalities `taxa`, `function`,
  `metabolites` and `host`.[^spec] Phase 2 already made a function table a
  MuData of two modalities, `function` and `function_by_taxon`
  ([function-tables-as-mudata](/decisions/function-tables-as-mudata.md)), and
  MuData cannot nest one MuData in another as a modality. MuData also accepts
  modalities whose samples only partly overlap: its global `obs` is the union and
  `obsmap` marks the gaps, which no paired method (mmvec, a joint model) can use.

  # Decision
  - Modality names: `taxa` (a TreeData or AnnData of taxa), `function` and
    `function_by_taxon` side by side (a function table's two modalities,
    unchanged), `metabolites` and `host`. Functions that read a modality take
    its name as a keyword defaulting to these. Other names are allowed;
    nothing checks them.
  - `io.to_mudata(modalities)` takes a mapping of name -> AnnData, keeps the
    samples every modality has in the first modality's order, copies each
    modality and warns once naming how many samples each loses. A MuData value
    raises `TypeError`: a function table goes in as `{**table.mod, ...}`.
  - The global `obs` stays as MuData builds it (no columns); `mdata.pull_obs()`
    gathers them.
  - Saving: `write_h5mu` drops a TreeData modality's `vart` and reads it back
    as AnnData ([tree-access](/contracts/tree-access.md), Gotchas). biotapy
    adds no writer; the docstring and the guide say to save that modality with
    `write_h5td` too.

  # Rejected
  - **Nested MuData** (`function` holding the function table): MuData's
    modalities must be AnnData.
  - **Union of samples** (MuData's default): a sample missing from one modality
    becomes NaN rows that a paired method would have to drop anyway.
  - **A biotapy `write_h5mu`** that writes trees beside the file: a second file
    format for one gotcha; reconsider when treedata or mudata supports it.
  - **Rejecting unknown modality names**: blocks a `proteins` or `viruses`
    modality for no benefit.

  # Consequences
  - Every function that pairs modalities checks that they hold the same samples
    in the same order and names `io.to_mudata` as the fix.
  - `datasets.hmp2` already follows the names (`function`,
    `function_by_taxon`, `taxa`).

  [^spec]: Python Microbiome Toolkit development report, section Data model
  ````
  Apply to the other concepts:
  ```diff
  diff --git a/.knowledge/contracts/data-model-slots.md b/.knowledge/contracts/data-model-slots.md
  index 7535610..57d5ae7 100644
  --- a/.knowledge/contracts/data-model-slots.md
  +++ b/.knowledge/contracts/data-model-slots.md
  @@ -120,6 +120,15 @@ sample metadata in the global `obs`, pushed into every modality
   (`datasets/_hmp2.py:hmp2`). `fn` verbs take the modality they need;
   `fn.renorm` keeps any other modality (`fn/_renorm.py:renorm`).

  +## Multi-omics
  +Several data types over the same samples are one `MuData`
  +([multiomics-as-mudata](/decisions/multiomics-as-mudata.md)). Modality names:
  +`taxa`, `function` and `function_by_taxon` (a function table's two, side by
  +side), `metabolites`, `host`. `io.to_mudata` builds it from a mapping of
  +AnnData, keeping the samples every modality has, in the first modality's
  +order, each modality a copy (`io/_mudata.py:to_mudata`). `write_h5mu` drops a
  +TreeData modality's tree.
  +
   ## Taxonomic profiles (MetaPhlAn)
   `io.read_metaphlan` keeps one feature per leaf clade: a row that no other
   row descends from through any ancestor (MetaPhlAn can omit an intermediate
  diff --git a/.knowledge/decisions/function-tables-as-mudata.md b/.knowledge/decisions/function-tables-as-mudata.md
  index ad11ba5..bb93778 100644
  --- a/.knowledge/decisions/function-tables-as-mudata.md
  +++ b/.knowledge/decisions/function-tables-as-mudata.md
  @@ -66,8 +66,9 @@ was brought forward from Phase 4, where multi-omics needs it anyway.
   - **Forward note for Phase 4 task 4.1** ([phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md)):
     a function table is itself a two-modality MuData, so a multi-omics MuData
     needs `function` and `function_by_taxon` side by side with `taxa`,
  -  `metabolites` and `host`, not a nested `function` MuData. Task 4.1's modality
  -  names must not reuse `function` for anything else, and `io.to_mudata` must
  -  accept a function table's two modalities as two entries. `datasets.hmp2`
  +  `metabolites` and `host`, not a nested `function` MuData. Task 4.1 does
  +  this: `io.to_mudata` takes a function table's two modalities as two entries
  +  and refuses a MuData value
  +  ([multiomics-as-mudata](/decisions/multiomics-as-mudata.md)). `datasets.hmp2`
     already holds `function`, `function_by_taxon` and `taxa` side by side
     (`datasets/_hmp2.py:hmp2`); `fn.renorm` keeps the extra modality.
  diff --git a/.knowledge/decisions/index.md b/.knowledge/decisions/index.md
  index 20c89d3..a834781 100644
  --- a/.knowledge/decisions/index.md
  +++ b/.knowledge/decisions/index.md
  @@ -3,6 +3,7 @@
   * [TreeData and MuData are the only containers](treedata-as-container.md) - biotapy owns no data class; one data type is a TreeData, several are a MuData, every feature is a function over them.
   * [Samples are rows](samples-as-rows.md) - Every matrix is samples x features (scverse orientation); importers transpose once, no other function checks orientation.
   * [Function tables are a two-modality MuData](function-tables-as-mudata.md) - A HUMAnN-style function table is a MuData with a community modality and a stratified modality, adopted in Phase 2 instead of Phase 4; mudata is a runtime dependency.
  +* [Multi-omics data is one MuData with fixed modality names](multiomics-as-mudata.md) - Several data types over the same samples are one MuData whose modalities are named taxa, function, function_by_taxon, metabolites and host; io.to_mudata keeps only the samples every modality has; a TreeData modality loses its tree in h5mu.
   * [Pure by default, one inplace convention for tl](pure-by-default.md) - io/pp return new objects and never mutate input; tl returns results, and inplace=True writes them to the documented slot; pl returns Axes.
   * [Python first, compiled code last](python-first-compiled-last.md) - Pure Python/NumPy by default; delegate to compiled libraries; Numba then Rust only for a benchmarked hotspot; never new C/C++.
   * [Heavy dependencies are optional extras](optional-heavy-dependencies.md) - torch, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily.
  ```
  Tick 4.1 in this concept and add the dated section to `.knowledge/log.md`
  (below the title, above the newest section):
  ```markdown
  ## <date> (Phase 4, slice 4A)
  - **Creation**: [multiomics-as-mudata](decisions/multiomics-as-mudata.md) (`draft`): modality names, `io.to_mudata` keeps the shared samples, no h5mu tree writer.
  - **Update**: [data-model-slots](contracts/data-model-slots.md) gains a Multi-omics section; [function-tables-as-mudata](decisions/function-tables-as-mudata.md)'s forward note says task 4.1 delivered it; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.1.
  ```
- [x] **Step 7: Gate and commit**
  ```bash
  git add src/biotapy/io/_mudata.py src/biotapy/io/__init__.py tests/io/test_mudata.py \
    docs/guide/multiomics.md docs/guide/index.md docs/api.md \
    .knowledge/decisions/multiomics-as-mudata.md .knowledge/decisions/index.md \
    .knowledge/decisions/function-tables-as-mudata.md .knowledge/contracts/data-model-slots.md \
    .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html
  uv run --group test pytest -q -W error::UserWarning
  git commit -m "feat(io): combine data types over shared samples into one MuData

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `build succeeded`; `1256 passed, 54 deselected`.

### Task 4.2: `tl.mmvec`

**Files:** create `src/biotapy/tl/_mmvec.py`, `tests/tl/test_mmvec.py`;
modify `src/biotapy/tl/__init__.py`, `docs/guide/multiomics.md`,
`docs/api.md`, `.knowledge/decisions/pure-by-default.md`,
`.knowledge/decisions/multiomics-as-mudata.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `tl/_ordination.py` (mmvec is not an ordination of samples
and writes no `obsm`); scikit-bio's other `MMvecResult` fields.
**Interfaces:**
- Consumes: `bt.io.to_mudata` (4.1, in tests); `_core.as_csr`,
  `as_generator`; `skbio.stats.ordination.mmvec`.
- Produces: `bt.tl.mmvec(mdata: MuData, *, microbes: str = "taxa",
  metabolites: str = "metabolites", seed: int | np.random.Generator | None =
  None) -> pd.DataFrame` (microbes x metabolites, rows sum to 0).

- [x] **Step 1: Failing tests.** Create `tests/tl/test_mmvec.py`:
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  import skbio
  from anndata import AnnData
  from hypothesis import given, settings
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays
  from mudata import MuData

  import biotapy as bt


  def _metabolites(counts: np.ndarray, samples) -> AnnData:
      """Three metabolites: one follows f6 (Prevotella), one f3 (Faecalibacterium), one is flat."""
      values = np.column_stack([counts[:, 5] + 1, counts[:, 2] + 1, np.full(counts.shape[0], 10)])
      return AnnData(
          X=sp.csr_matrix(values.astype(np.float64)),
          obs=pd.DataFrame(index=list(samples)),
          var=pd.DataFrame(index=["m_prev", "m_faec", "m_flat"]),
      )


  def _toy_mudata():
      tdata = bt.datasets.toy()
      return bt.io.to_mudata({"taxa": tdata, "metabolites": _metabolites(tdata.X.toarray(), tdata.obs_names)})


  def test_ranks_are_microbes_by_metabolites():
      ranks = bt.tl.mmvec(_toy_mudata(), seed=0)
      assert isinstance(ranks, pd.DataFrame)
      assert ranks.index.tolist() == [f"f{i}" for i in range(1, 9)]
      assert ranks.columns.tolist() == ["m_prev", "m_faec", "m_flat"]
      np.testing.assert_allclose(ranks.sum(axis=1), 0.0, atol=1e-10)


  def test_a_metabolite_ranks_highest_for_the_microbe_it_follows():
      ranks = bt.tl.mmvec(_toy_mudata(), seed=0)
      assert ranks.loc["f6"].idxmax() == "m_prev"
      assert ranks.loc["f3"].idxmax() == "m_faec"


  def test_same_seed_same_ranks():
      first = bt.tl.mmvec(_toy_mudata(), seed=1)
      again = bt.tl.mmvec(_toy_mudata(), seed=np.random.default_rng(1))
      pd.testing.assert_frame_equal(first, again)


  def test_modality_names_are_arguments():
      mdata = _toy_mudata()
      renamed = bt.io.to_mudata({"microbes": mdata["taxa"], "compounds": mdata["metabolites"]})
      pd.testing.assert_frame_equal(
          bt.tl.mmvec(renamed, microbes="microbes", metabolites="compounds", seed=0), bt.tl.mmvec(mdata, seed=0)
      )


  def test_ranks_stay_a_labelled_dataframe_whatever_scikit_bio_outputs():
      previous = skbio.get_config("table_output")
      skbio.set_config("table_output", "numpy")
      try:
          ranks = bt.tl.mmvec(_toy_mudata(), seed=0)
      finally:
          skbio.set_config("table_output", previous)
      assert isinstance(ranks, pd.DataFrame)
      assert ranks.columns.tolist() == ["m_prev", "m_faec", "m_flat"]


  def test_keeps_the_input(assert_unchanged):
      mdata = _toy_mudata()
      before = {name: mod.copy() for name, mod in mdata.mod.items()}
      bt.tl.mmvec(mdata, seed=0)
      for name, mod in mdata.mod.items():
          assert_unchanged(before[name], mod)


  def test_missing_modality_raises():
      with pytest.raises(KeyError, match=r"metabolites='metabolome' is not a modality; found \['taxa', 'metabolites'\]"):
          bt.tl.mmvec(_toy_mudata(), metabolites="metabolome")


  def test_unaligned_samples_raise():
      tdata = bt.datasets.toy()
      metabolites = _metabolites(tdata.X.toarray(), tdata.obs_names)[::-1].copy()
      with pytest.raises(ValueError, match=r"hold different samples.*bt.io.to_mudata"):
          bt.tl.mmvec(MuData({"taxa": tdata, "metabolites": metabolites}))


  def test_all_zero_feature_raises():
      mdata = _toy_mudata()
      mdata["taxa"].X = sp.csr_matrix(np.column_stack([mdata["taxa"].X.toarray()[:, :7], np.zeros(6)]))
      with pytest.raises(ValueError, match=r"microbes='taxa' has all-zero features: \['f8'\]"):
          bt.tl.mmvec(mdata)


  def test_all_zero_sample_raises():
      mdata = _toy_mudata()
      values = mdata["metabolites"].X.toarray()
      values[2] = 0.0
      mdata["metabolites"].X = sp.csr_matrix(values)
      with pytest.raises(ValueError, match=r"metabolites='metabolites' has all-zero samples: \['s3'\]"):
          bt.tl.mmvec(mdata)


  def test_single_sample():
      tdata = bt.pp.filter_features(bt.datasets.toy()[:1].copy(), min_total=1)
      metabolites = AnnData(
          X=sp.csr_matrix([[3.0, 4.0, 5.0]]), obs=pd.DataFrame(index=["s1"]), var=pd.DataFrame(index=["m1", "m2", "m3"])
      )
      assert bt.tl.mmvec(bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites}), seed=0).shape == (6, 3)


  @pytest.mark.parametrize("bad", [-1.0, np.nan])
  def test_negative_or_missing_value_raises(bad):
      mdata = _toy_mudata()
      values = mdata["metabolites"].X.toarray()
      values[0, 0] = bad
      mdata["metabolites"].X = sp.csr_matrix(values)
      with pytest.raises(ValueError, match="tl.mmvec needs finite, non-negative values in metabolites='metabolites'"):
          bt.tl.mmvec(mdata)


  @settings(max_examples=25, deadline=None)
  @given(arrays(np.int64, st.tuples(st.integers(2, 5), st.integers(2, 4)), elements=st.integers(1, 50)))
  def test_every_microbe_row_is_centred(counts):
      samples = [f"s{i}" for i in range(counts.shape[0])]
      taxa = AnnData(
          X=sp.csr_matrix(counts),
          obs=pd.DataFrame(index=samples),
          var=pd.DataFrame(index=[f"f{i}" for i in range(counts.shape[1])]),
      )
      metabolites = AnnData(
          X=sp.csr_matrix(counts[:, ::-1] + 1),
          obs=pd.DataFrame(index=samples),
          var=pd.DataFrame(index=[f"m{i}" for i in range(counts.shape[1])]),
      )
      ranks = bt.tl.mmvec(bt.io.to_mudata({"taxa": taxa, "metabolites": metabolites}), seed=0)
      np.testing.assert_allclose(ranks.sum(axis=1), 0.0, atol=1e-9)
  ```
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/tl/test_mmvec.py -q` -> `14 failed`
  (`AttributeError: module 'biotapy.tl' has no attribute 'mmvec'`).
- [x] **Step 3: Implement.** Create `src/biotapy/tl/_mmvec.py`:
  ```python
  """mmvec: which metabolites co-occur with which microbes, over samples of both."""

  from typing import cast

  import numpy as np
  import pandas as pd
  from anndata import AnnData
  from mudata import MuData
  from skbio.stats.ordination import mmvec as skbio_mmvec

  from biotapy._core import as_csr, as_generator


  def mmvec(
      mdata: MuData,
      *,
      microbes: str = "taxa",
      metabolites: str = "metabolites",
      seed: int | np.random.Generator | None = None,
  ) -> pd.DataFrame:
      """Learn how likely each metabolite is given each microbe (Morton et al. 2019).

      Parameters
      ----------
      mdata
          Modalities over the same samples, in the same order, such as
          :func:`biotapy.io.to_mudata` returns; ``X`` holds counts or another
          non-negative abundance.
      microbes
          The modality whose features condition the model.
      metabolites
          The modality whose features are predicted.
      seed
          Seed or generator for the starting values of the fit.

      Returns
      -------
      pandas.DataFrame
          Microbes (rows, ``mdata[microbes].var_names``) x metabolites (columns):
          the log conditional probability of each metabolite given each microbe,
          centred so that every row sums to 0. Larger means more likely to
          co-occur.

      Raises
      ------
      KeyError
          ``microbes`` or ``metabolites`` is not a modality.
      ValueError
          The two modalities hold different samples or a different order; a value
          is negative or not finite; or a feature or sample is all zero.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/multiomics`

      Wraps :func:`skbio.stats.ordination.mmvec` with its defaults: three
      latent dimensions and up to 1,000 L-BFGS iterations. For the embeddings,
      the fitted probabilities or a prediction on new samples, call scikit-bio
      with the same tables.

      scikit-bio needs dense tables, so both ``X`` are densified once:
      8 bytes x samples x features of each modality. The fit and the result also
      hold dense microbes x metabolites arrays: 8 bytes x microbes x metabolites
      each.

      References
      ----------
      Morton JT et al. (2019) Learning representations of microbe-metabolite
      interactions. Nature Methods 16:1306-1314.

      Examples
      --------
      >>> import anndata as ad
      >>> import biotapy as bt
      >>> tdata = bt.datasets.toy()
      >>> metabolites = ad.AnnData(tdata.X[:, [5, 2]].toarray() + 1.0, obs=tdata.obs[[]])
      >>> mdata = bt.io.to_mudata({"taxa": tdata, "metabolites": metabolites})
      >>> ranks = bt.tl.mmvec(mdata, seed=0)
      >>> ranks.shape
      (8, 2)
      """
      rng = as_generator(seed)
      x_mod = _modality(mdata, microbes, argument="microbes")
      y_mod = _modality(mdata, metabolites, argument="metabolites")
      if not x_mod.obs_names.equals(y_mod.obs_names):
          msg = (
              f"microbes={microbes!r} and metabolites={metabolites!r} hold different samples or a different order; "
              "align them with bt.io.to_mudata"
          )
          raise ValueError(msg)
      x_table = _table(x_mod, microbes, argument="microbes")
      y_table = _table(y_mod, metabolites, argument="metabolites")
      return cast("pd.DataFrame", skbio_mmvec(x_table, y_table, seed=rng, output_format="pandas").ranks)


  def _modality(mdata: MuData, key: str, *, argument: str) -> AnnData:
      if key not in mdata.mod:
          msg = f"{argument}={key!r} is not a modality; found {list(mdata.mod)}"
          raise KeyError(msg)
      return cast("AnnData", mdata.mod[key])


  def _table(mod: AnnData, key: str, *, argument: str) -> pd.DataFrame:
      """Modality ``key`` as a dense samples x features DataFrame, checked as mmvec needs it."""
      X = as_csr(mod.X).astype(np.float64)
      if not np.all(np.isfinite(X.data)) or np.any(X.data < 0):
          msg = f"tl.mmvec needs finite, non-negative values in {argument}={key!r}"
          raise ValueError(msg)
      for axis, names, what in ((0, mod.var_names, "features"), (1, mod.obs_names, "samples")):
          empty = names[np.asarray(X.sum(axis=axis)).ravel() == 0]
          if len(empty):
              msg = f"{argument}={key!r} has all-zero {what}: {empty[:5].tolist()}"
              raise ValueError(msg)
      # scikit-bio's mmvec needs dense input (rules.md R6.2): one dense copy per modality.
      return pd.DataFrame(X.toarray(), index=mod.obs_names, columns=mod.var_names)
  ```
  `src/biotapy/tl/__init__.py`:
  ```python
  """Diversity and ordination on AnnData/TreeData (contracts/module-boundaries)."""

  from ._alpha import alpha
  from ._beta import beta, unifrac
  from ._mmvec import mmvec
  from ._ordination import nmds, pcoa
  from ._permanova import permanova

  __all__ = ["alpha", "beta", "mmvec", "nmds", "pcoa", "permanova", "unifrac"]
  ```
- [x] **Step 4: Run, expect pass** -
  `uv run --group test pytest tests/tl/test_mmvec.py src/biotapy/tl/_mmvec.py -q -W error::UserWarning`
  -> `15 passed`; the property test also under `--hypothesis-seed=1`, `2`,
  `3`.
- [x] **Step 5: Docs and knowledge.**
  ```diff
  diff --git a/.knowledge/decisions/multiomics-as-mudata.md b/.knowledge/decisions/multiomics-as-mudata.md
  index 23c6b1c..b88d084 100644
  --- a/.knowledge/decisions/multiomics-as-mudata.md
  +++ b/.knowledge/decisions/multiomics-as-mudata.md
  @@ -4,7 +4,7 @@ title: Multi-omics data is one MuData with fixed modality names
   description: Several data types over the same samples are one MuData whose modalities are named taxa, function, function_by_taxon, metabolites and host; io.to_mudata keeps only the samples every modality has; a TreeData modality loses its tree in h5mu.
   tags: [io, mudata, multiomics]
   status: draft
  -paths: ["src/biotapy/io/_mudata.py"]
  +paths: ["src/biotapy/io/_mudata.py", "src/biotapy/tl/_mmvec.py"]
   generated: { by: claude-code/claude-opus-5-5, at: 2026-10-08T17:54:32Z }
   commit: 872c16d
   sources:
  @@ -27,8 +27,8 @@ modalities whose samples only partly overlap: its global `obs` is the union and
   - Modality names: `taxa` (a TreeData or AnnData of taxa), `function` and
     `function_by_taxon` side by side (a function table's two modalities,
     unchanged), `metabolites` and `host`. Functions that read a modality take
  -  its name as a keyword defaulting to these. Other names are allowed;
  -  nothing checks them.
  +  its name as a keyword defaulting to these (`tl.mmvec(microbes="taxa",
  +  metabolites="metabolites")`). Other names are allowed; nothing checks them.
   - `io.to_mudata(modalities)` takes a mapping of name -> AnnData, keeps the
     samples every modality has in the first modality's order, copies each
     modality and warns once naming how many samples each loses. A MuData value
  @@ -52,7 +52,8 @@ modalities whose samples only partly overlap: its global `obs` is the union and

   # Consequences
   - Every function that pairs modalities checks that they hold the same samples
  -  in the same order and names `io.to_mudata` as the fix.
  +  in the same order and names `io.to_mudata` as the fix (`tl.mmvec`,
  +  `tl/_mmvec.py`).
   - `datasets.hmp2` already follows the names (`function`,
     `function_by_taxon`, `taxa`).

  diff --git a/.knowledge/decisions/pure-by-default.md b/.knowledge/decisions/pure-by-default.md
  index 533b301..6692999 100644
  --- a/.knowledge/decisions/pure-by-default.md
  +++ b/.knowledge/decisions/pure-by-default.md
  @@ -31,7 +31,7 @@ behaviour. Confirmed by the user on 2026-09-26.
   | `pl` | `matplotlib.axes.Axes` | never |

   Every `tl` function that returns per-sample or per-pair values supports both modes with
  -identical semantics; `tl.permanova` is the exception (see Consequences).
  +identical semantics; `tl.permanova` and `tl.mmvec` are the exceptions (see Consequences).

   # Rejected
   - **scanpy default (`copy=False`, mutate)**: contradicts the spec's purity rule.
  @@ -41,7 +41,9 @@ identical semantics; `tl.permanova` is the exception (see Consequences).
   # Consequences
   - Tests assert the input object is unchanged after every `pp`/`tl` call.
   - `tl.permanova` returns a test result, not per-sample or per-pair values, so
  -  it has no slot and no `inplace`. No `tl` function takes `key_added` until a
  +  it has no slot and no `inplace`. `tl.mmvec` returns a microbes x metabolites
  +  table across two modalities, which no slot of one AnnData holds, so it has no
  +  `inplace` either. No `tl` function takes `key_added` until a
     use case needs one (rules.md R2.3).

   [^spec]: Python Microbiome Toolkit development report, section Function-level implementation
  diff --git a/docs/api.md b/docs/api.md
  index 8c37070..503bd3d 100644
  --- a/docs/api.md
  +++ b/docs/api.md
  @@ -86,6 +86,7 @@ Public functions are listed here as they ship, from Phase 1 onward.

       tl.alpha
       tl.beta
  +    tl.mmvec
       tl.nmds
       tl.pcoa
       tl.permanova
  diff --git a/docs/guide/multiomics.md b/docs/guide/multiomics.md
  index b5e079e..8f178b2 100644
  --- a/docs/guide/multiomics.md
  +++ b/docs/guide/multiomics.md
  @@ -40,6 +40,29 @@ Each modality is a copy, so the tables you passed are unchanged.
   The global `mdata.obs` starts without columns; `mdata.pull_obs()` gathers
   the modalities' sample metadata into it.

  +## Microbes and metabolites: mmvec
  +
  +`bt.tl.mmvec` fits mmvec (Morton et al. 2019) through scikit-bio: it learns,
  +from samples with both tables, how likely each metabolite is given each
  +microbe.
  +
  +```python
  +ranks = bt.tl.mmvec(mdata, seed=0)  # microbes="taxa", metabolites="metabolites"
  +ranks.loc["f6"].sort_values(ascending=False).head()  # metabolites most tied to f6
  +```
  +
  +`ranks` has one row per microbe and one column per metabolite, holding the log
  +probability of the metabolite given the microbe, centred so each row sums to
  +0. Compare values within a row: a large value means the metabolite tends to be
  +abundant in samples where that microbe is. A low value means no association,
  +not necessarily a negative one.
  +
  +The two modalities must hold the same samples in the same order, which
  +`bt.io.to_mudata` guarantees. A feature or sample that is zero everywhere has
  +nothing to learn from, so biotapy refuses it and names it. Drop an all-zero
  +feature with `bt.pp.filter_features(..., min_total=1)`; drop an all-zero sample
  +with `bt.pp.filter_samples` on that modality, then rebuild the MuData with
  +`bt.io.to_mudata`. The fit starts from random
  +values, so pass `seed` for the same ranks every time.
  +
   ## Saving

   `mdata.write_h5mu("study.h5mu")` saves every modality, but a TreeData
  ```
  Tick 4.2 here and add to the slice 4A log section:
  ```markdown
  - **Update**: [pure-by-default](decisions/pure-by-default.md): `tl.mmvec` returns a table and has no `inplace`, like `tl.permanova`; [multiomics-as-mudata](decisions/multiomics-as-mudata.md) names `tl.mmvec` as the first function that pairs modalities; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.2.
  ```
- [x] **Step 6: Gate and commit**
  ```bash
  git add src/biotapy/tl/_mmvec.py src/biotapy/tl/__init__.py tests/tl/test_mmvec.py \
    docs/guide/multiomics.md docs/api.md .knowledge/decisions/pure-by-default.md \
    .knowledge/decisions/multiomics-as-mudata.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html
  uv run --group test pytest -q -W error::UserWarning
  git commit -m "feat(tl): add mmvec over two modalities of a MuData

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `build succeeded`; `1272 passed, 54 deselected`.

### Task 4.A0: the pseudocount step moves to `_core` (`refactor(core)`)

**Files:** create `src/biotapy/_core/_composition.py`,
`tests/core/test_composition.py`; modify `src/biotapy/_core/__init__.py`,
`src/biotapy/pp/_transform.py`, `src/biotapy/pp/_philr.py`,
`.knowledge/modules/core.md`, `.knowledge/modules/pp.md`,
`.knowledge/log.md`.
**Not touched:** any `tests/pp` file (they must pass unedited); `pp.clr`'s
and `pp.philr`'s signatures, docstrings and messages.
**Interfaces:**
- Produces: `_core.check_pseudocount(pseudocount: object) -> None`;
  `_core.pseudocounted(X: object, pseudocount: float, *, func: str,
  columns: npt.NDArray[np.intp] | None = None) -> npt.NDArray[np.float64]`
  (takes `X`, where Phase 3's took the AnnData).

- [x] **Step 1: Failing tests.** Create `tests/core/test_composition.py`:
  ```python
  import numpy as np
  import pytest
  import scipy.sparse as sp

  from biotapy._core import check_pseudocount, pseudocounted


  def test_dense_and_sparse_x_give_the_same_values():
      dense = np.array([[0.0, 2.0], [3.0, 1.0]])
      expected = dense + 0.5
      np.testing.assert_array_equal(pseudocounted(dense, 0.5, func="f"), expected)
      np.testing.assert_array_equal(pseudocounted(sp.csr_matrix(dense), 0.5, func="f"), expected)


  def test_columns_are_taken_in_order():
      dense = np.array([[1.0, 2.0, 3.0]])
      np.testing.assert_array_equal(pseudocounted(dense, 1, func="f", columns=np.array([2, 0])), [[4.0, 2.0]])


  def test_does_not_change_its_input():
      dense = np.array([[1.0, 2.0]])
      pseudocounted(dense, 0.5, func="f")
      np.testing.assert_array_equal(dense, [[1.0, 2.0]])


  @pytest.mark.parametrize(
      ("value", "error"), [(True, TypeError), ("1", TypeError), (-1, ValueError), (np.inf, ValueError)]
  )
  def test_check_pseudocount_refuses(value, error):
      with pytest.raises(error, match="pseudocount must be"):
          check_pseudocount(value)


  def test_messages_name_the_caller():
      with pytest.raises(ValueError, match="ml.CLR needs finite, non-negative values in X"):
          pseudocounted(np.array([[-1.0, 2.0]]), 0.5, func="ml.CLR")


  def test_pseudocount_above_the_smallest_value_warns_naming_it():
      values = np.array([[0.0, 0.2], [0.8, 0.5]])
      with pytest.warns(UserWarning, match=r"pseudocount=0\.5 is larger than the smallest non-zero value in X \(0\.2\)"):
          pseudocounted(values, 0.5, func="ml.CLR")


  def test_pseudocount_warning_ends_with_the_calling_step():
      with pytest.warns(UserWarning, match=r"on their scale \(ml\.CLR\)$"):
          pseudocounted(np.array([[0.0, 0.2], [0.8, 0.5]]), 0.5, func="ml.CLR")


  def test_zero_pseudocount_with_a_zero_in_x_raises():
      with pytest.raises(ValueError, match="X holds zeros, whose logarithm is undefined; pass pseudocount > 0 to ml.CLR"):
          pseudocounted(np.array([[0.0, 2.0]]), 0, func="ml.CLR")
  ```
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/core/test_composition.py -q` ->
  `1 error` (`ImportError: cannot import name 'check_pseudocount' from
  'biotapy._core'`).
- [x] **Step 3: Implement.** Create `src/biotapy/_core/_composition.py`:
  ```python
  """The pseudocount step shared by pp's log-ratio transforms and ml's CLR."""

  import numpy as np
  import numpy.typing as npt

  from ._matrix import as_csr
  from ._warnings import warn_user


  def check_pseudocount(pseudocount: object) -> None:
      """Raise unless ``pseudocount`` is a finite real number >= 0; a bool is refused."""
      if isinstance(pseudocount, bool) or not isinstance(pseudocount, int | float | np.integer | np.floating):
          msg = f"pseudocount must be a real number, got {pseudocount!r}"
          raise TypeError(msg)
      if not np.isfinite(pseudocount) or pseudocount < 0:
          msg = f"pseudocount must be a finite number >= 0, got {pseudocount!r}"
          raise ValueError(msg)


  def pseudocounted(
      X: object, pseudocount: float, *, func: str, columns: npt.NDArray[np.intp] | None = None
  ) -> npt.NDArray[np.float64]:
      """``X`` (its ``columns``, in that order, if given) as a dense float64 array plus ``pseudocount``, checked > 0.

      ``X`` is anything :func:`as_csr` takes. ``func`` names the caller in messages.
      """
      check_pseudocount(pseudocount)
      matrix = as_csr(X).astype(np.float64)
      if columns is not None:
          # Reordered while sparse, so the dense copy below is the only one.
          matrix = matrix[:, columns]
      if not np.all(np.isfinite(matrix.data)) or np.any(matrix.data < 0):
          msg = f"{func} needs finite, non-negative values in X"
          raise ValueError(msg)
      positive = matrix.data[matrix.data > 0]
      if positive.size and pseudocount > positive.min():
          warn_user(
              f"pseudocount={pseudocount} is larger than the smallest non-zero value in X ({positive.min():.3g}), "
              f"so it swamps the rarest features; for relative abundances pass a pseudocount on their scale ({func})"
          )
      # scikit-bio's log-ratio functions need dense input (rules.md R6.2): one dense copy of X.
      values = matrix.toarray()
      values += pseudocount
      if np.any(values <= 0):
          msg = f"X holds zeros, whose logarithm is undefined; pass pseudocount > 0 to {func}"
          raise ValueError(msg)
      return values
  ```
  and move the callers:
  ```diff
  diff --git a/src/biotapy/_core/__init__.py b/src/biotapy/_core/__init__.py
  index 879dbf9..cc52e45 100644
  --- a/src/biotapy/_core/__init__.py
  +++ b/src/biotapy/_core/__init__.py
  @@ -1,5 +1,6 @@
   """Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

  +from ._composition import check_pseudocount, pseudocounted
   from ._function import (
       BY_TAXON_KEY,
       FUNCTION_KEY,
  @@ -53,6 +54,7 @@ __all__ = [
       "argmax_by",
       "as_csr",
       "as_generator",
  +    "check_pseudocount",
       "divide_rows",
       "feature_subset",
       "function_var",
  @@ -63,6 +65,7 @@ __all__ = [
       "make_function_mudata",
       "make_treedata",
       "normalize_ranks",
  +    "pseudocounted",
       "relabel_tips",
       "replace_features",
       "require_categorical",
  diff --git a/src/biotapy/pp/_philr.py b/src/biotapy/pp/_philr.py
  index e4ce0c1..b6f99a8 100644
  --- a/src/biotapy/pp/_philr.py
  +++ b/src/biotapy/pp/_philr.py
  @@ -4,9 +4,7 @@ import pandas as pd
   from skbio import TreeNode
   from skbio.stats.composition import clr, tree_basis

  -from biotapy._core import TreeData, add_provenance, get_skbio_tree
  -
  -from ._transform import pseudocounted
  +from biotapy._core import TreeData, add_provenance, get_skbio_tree, pseudocounted


   def philr(tdata: TreeData, *, pseudocount: float = 0.5) -> TreeData:
  @@ -91,7 +89,7 @@ def philr(tdata: TreeData, *, pseudocount: float = 0.5) -> TreeData:
       if len(tips) != tdata.n_vars:
           msg = f"pp.philr needs every feature to be a tip of the tree; {tdata.n_vars - len(tips)} feature(s) are not"
           raise ValueError(msg)
  -    values = pseudocounted(tdata, pseudocount, func="pp.philr", columns=tdata.var_names.get_indexer(tips))
  +    values = pseudocounted(tdata.X, pseudocount, func="pp.philr", columns=tdata.var_names.get_indexer(tips))
       basis, nodes = tree_basis(tree)
       # tree_basis puts a node's first child in the denominator; philr::philr puts it in the numerator.
       balances = pd.DataFrame(-(clr(values) @ basis.T), index=tdata.obs_names, columns=nodes)
  diff --git a/src/biotapy/pp/_transform.py b/src/biotapy/pp/_transform.py
  index 7c1a2b1..d054661 100644
  --- a/src/biotapy/pp/_transform.py
  +++ b/src/biotapy/pp/_transform.py
  @@ -1,11 +1,10 @@
   """Per-sample transforms: add one layer, keep everything else."""

   import numpy as np
  -import numpy.typing as npt
   from anndata import AnnData
   from skbio.stats.composition import clr as skbio_clr

  -from biotapy._core import add_provenance, as_csr, divide_rows, warn_user
  +from biotapy._core import add_provenance, as_csr, divide_rows, pseudocounted


   def relative(adata: AnnData) -> AnnData:
  @@ -104,40 +103,8 @@ def clr(adata: AnnData, *, pseudocount: float = 0.5) -> AnnData:
       >>> round(float(out.layers["clr"][0, 0]), 3)
       1.048
       """
  -    values = pseudocounted(adata, pseudocount, func="pp.clr")
  +    values = pseudocounted(adata.X, pseudocount, func="pp.clr")
       out = adata.copy()
       out.layers["clr"] = skbio_clr(values)
       add_provenance(out, "pp.clr", pseudocount=pseudocount)
       return out
  -
  -
  -def pseudocounted(
  -    adata: AnnData, pseudocount: float, *, func: str, columns: npt.NDArray[np.intp] | None = None
  -) -> npt.NDArray[np.float64]:
  -    """``X`` (its ``columns``, in that order, if given) as a dense float64 array plus ``pseudocount``, checked > 0."""
  -    if isinstance(pseudocount, bool) or not isinstance(pseudocount, int | float | np.integer | np.floating):
  -        msg = f"pseudocount must be a real number, got {pseudocount!r}"
  -        raise TypeError(msg)
  -    if not np.isfinite(pseudocount) or pseudocount < 0:
  -        msg = f"pseudocount must be a finite number >= 0, got {pseudocount!r}"
  -        raise ValueError(msg)
  -    X = as_csr(adata.X).astype(np.float64)
  -    if columns is not None:
  -        # Reordered while sparse, so the dense copy below is the only one.
  -        X = X[:, columns]
  -    if not np.all(np.isfinite(X.data)) or np.any(X.data < 0):
  -        msg = f"{func} needs finite, non-negative values in X"
  -        raise ValueError(msg)
  -    positive = X.data[X.data > 0]
  -    if positive.size and pseudocount > positive.min():
  -        warn_user(
  -            f"pseudocount={pseudocount} is larger than the smallest non-zero value in X ({positive.min():.3g}), "
  -            "so it swamps the rarest features; for relative abundances pass a pseudocount on their scale"
  -        )
  -    # scikit-bio's log-ratio functions need dense input (rules.md R6.2): one dense copy of X.
  -    values = X.toarray()
  -    values += pseudocount
  -    if np.any(values <= 0):
  -        msg = f"X holds zeros, whose logarithm is undefined; pass pseudocount > 0 to {func}"
  -        raise ValueError(msg)
  -    return values
  ```
- [x] **Step 4: Run, expect pass** -
  `uv run --group test pytest tests/core/test_composition.py tests/pp -q -W error::UserWarning`
  -> `118 passed, 8 deselected`.
- [x] **Step 5: Knowledge.**
  ```diff
  diff --git a/.knowledge/modules/core.md b/.knowledge/modules/core.md
  index c71c00c..d5d22b8 100644
  --- a/.knowledge/modules/core.md
  +++ b/.knowledge/modules/core.md
  @@ -12,7 +12,8 @@ status: stable

   # Responsibility

  -Owns: sparse matrix kernels (`_matrix.py`), the canonical taxonomic rank order
  +Owns: sparse matrix kernels (`_matrix.py`), the pseudocount check shared by
  +log-ratio transforms (`_composition.py`), the canonical taxonomic rank order
   and lineage/rank-column parsing (`_taxonomy.py`), the `x_kind`/provenance/
   `feature_subset`/`replace_features` slot rules (`_slots.py`), function-table
   construction and the HUMAnN special-row constants (`_function.py`), TreeData construction and the only
  @@ -41,6 +42,12 @@ none of them back.
     full toward each (many-to-many); the pairs are a set, a repeated pair counts
     once; used by `fn.func_glom`. Where `sum_by` assigns each feature one group,
     this does not.
  +- `_composition.py:pseudocounted` - validates `pseudocount`
  +  (`_composition.py:check_pseudocount`), rejects negative or non-finite `X`,
  +  warns when the pseudocount exceeds the smallest non-zero value, and returns
  +  `X + pseudocount` as one dense float64 array; `columns=` reorders features
  +  while still sparse. Used by `pp.clr`, `pp.philr` and `ml.CLR` (which also
  +  calls `check_pseudocount` in `fit`); moved from `pp/_transform.py` in Phase 4.
   - `_taxonomy.py:split_ranks` - split `var`'s taxonomy columns at a target
     rank, raising `KeyError` naming the rank when it is absent.
   - `_taxonomy.py:normalize_ranks` - canonicalize rank column names (aliases,
  diff --git a/.knowledge/modules/pp.md b/.knowledge/modules/pp.md
  index 611f615..01cd295 100644
  --- a/.knowledge/modules/pp.md
  +++ b/.knowledge/modules/pp.md
  @@ -27,11 +27,9 @@ does NOT own any diversity/ordination computation ([tl](/modules/tl.md)).
     overflows for a subnormal total); all-zero samples stay all-zero, which differs from phyloseq's
     `transform_sample_counts` (returns `NaN` there).
   - `_transform.py:clr` - centred log-ratio into `layers["clr"]`, through
  -  scikit-bio's `clr` on `_transform.py:pseudocounted`.
  -- `_transform.py:pseudocounted` - the check-and-densify step `clr` and `philr`
  -  share: validates `pseudocount`, rejects negative or non-finite `X`, adds the
  -  pseudocount to every value and returns one dense float64 array; `columns=`
  -  reorders features while still sparse.
  +  scikit-bio's `clr` on `_core._composition.py:pseudocounted`, the
  +  check-and-densify step `clr`, `philr` and `ml.CLR` share
  +  ([core](/modules/core.md)).
   - `_philr.py:philr` - PhILR balances into `obsm["X_philr"]`, through
     scikit-bio's `tree_basis` on the tree from `_core.get_skbio_tree`.
   - `_philr.py:_binary_tree` - copy of the tree without one-child nodes, children
  @@ -84,7 +82,7 @@ does NOT own any diversity/ordination computation ([tl](/modules/tl.md)).
     `tests/pp/test_transform_golden.py` and `tests/pp/test_philr_golden.py`). A
     pseudocount of 0 raises when `X` holds a zero. A pseudocount above the
     smallest non-zero value of `X` gives one `UserWarning` through `warn_user`
  -  (the default 0.5 meets relative abundances). `_transform.py:pseudocounted`
  +  (the default 0.5 meets relative abundances). `_core._composition.py:pseudocounted`
   - `layers["clr"]` is dense float64, because CLR has no zeros; an all-zero
     sample gives an all-zero row. `_transform.py:clr`
   - `obsm["X_philr"]` is a samples x balances `DataFrame`, one column per
  @@ -145,10 +143,9 @@ does NOT own any diversity/ordination computation ([tl](/modules/tl.md)).
     sparse basis, which holds one value per tip under each node (113 MB on
     GlobalPatterns' 26 x 19,216). Both state it in their `Notes`.
     `_transform.py:clr`, `_philr.py:philr`.
  -- `pseudocounted` lives in `_transform.py` and `_philr.py` imports it from
  -  there, which module-boundaries allows inside one subpackage. If slice 3B's
  -  `da.linda` reuses it, it moves to `_core` (two or more subpackages,
  -  [module-boundaries](/contracts/module-boundaries.md) rule 2).
  +- `pseudocounted` lives in `_core/_composition.py` since Phase 4, when
  +  `ml.CLR` became its second subpackage
  +  ([module-boundaries](/contracts/module-boundaries.md) rule 2).
   - `filter_features` keeps a feature when `present / n_obs >= min_prevalence`;
     phyloseq's `sum(x > 0) >= p * length(x)` can drop it at an exact boundary
     through floating point (7 of 25 samples at `p = 0.28`).
  ```
  Log line:
  ```markdown
  - **Update**: [core](modules/core.md) owns `_composition.py:pseudocounted` and `check_pseudocount`, moved from `pp/_transform.py` for `ml.CLR`; [pp](modules/pp.md) points to it.
  ```
- [x] **Step 6: Gate and commit**
  ```bash
  git add src/biotapy/_core/_composition.py src/biotapy/_core/__init__.py src/biotapy/pp/_transform.py \
    src/biotapy/pp/_philr.py tests/core/test_composition.py .knowledge/modules/core.md \
    .knowledge/modules/pp.md .knowledge/log.md
  uvx prek run --all-files
  uv run --group test pytest -q -W error::UserWarning
  git commit -m "refactor(core): share the pseudocount step between pp and ml

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1282 passed, 54 deselected`.

### Task 4.3: `ml.PrevalenceFilter` and `ml.CLR`

**Files:** create `src/biotapy/ml/__init__.py`,
`src/biotapy/ml/_transformers.py`, `tests/ml/test_transformers.py`,
`docs/guide/machine_learning.md`, `docs/_templates/autosummary/class.rst`;
modify `src/biotapy/__init__.py`, `tests/test_docstrings.py`,
`docs/guide/index.md`, `docs/api.md`, `docs/conf.py`,
`.knowledge/decisions/pure-by-default.md`,
`.knowledge/contracts/function-shape.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `pp.filter_features` and `pp.clr` (the transformers reuse
their rule and helper, not their code paths); `pyproject.toml` (the `ml`
layer already exists in import-linter, and mypy needs no new override);
`docs/extensions/coming_from_r.py` (the classes say `R equivalent: none`,
so they add no row).
**Interfaces:**
- Consumes: `_core.as_csr`, `check_pseudocount`, `pseudocounted` (4.A0);
  scikit-learn's `BaseEstimator`, `SelectorMixin`, `OneToOneFeatureMixin`,
  `TransformerMixin`, `Tags`, `validate_data`, `check_is_fitted`,
  `check_non_negative`; `skbio.stats.composition.clr`.
- Produces: `bt.ml.PrevalenceFilter(*, min_prevalence: float = 0.1)` with
  `prevalence_`; `bt.ml.CLR(*, pseudocount: float = 0.5)`.

- [x] **Step 1: Failing tests.** Create `tests/ml/test_transformers.py`:
  ```python
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from hypothesis import given
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays
  from sklearn.exceptions import NotFittedError
  from sklearn.linear_model import LogisticRegression
  from sklearn.model_selection import StratifiedKFold, cross_validate
  from sklearn.pipeline import make_pipeline
  from sklearn.utils.estimator_checks import parametrize_with_checks

  import biotapy as bt


  # The checks feed random floats below CLR's default pseudocount of 0.5, which CLR warns about
  # (tested below); the warning is not what they check.
  @pytest.mark.filterwarnings("ignore:pseudocount=")
  @parametrize_with_checks([bt.ml.PrevalenceFilter(), bt.ml.CLR()])
  def test_scikit_learn_estimator_checks(estimator, check):
      check(estimator)


  def _toy_x() -> sp.csr_matrix:
      return bt.datasets.toy().X


  def test_prevalence_filter_keeps_what_pp_filter_features_keeps():
      tdata = bt.datasets.toy()
      kept = bt.pp.filter_features(tdata, min_prevalence=0.8).var_names
      selector = bt.ml.PrevalenceFilter(min_prevalence=0.8).fit(tdata.X)
      assert tdata.var_names[selector.get_support()].tolist() == kept.tolist()


  def test_prevalence_filter_learns_only_from_the_samples_it_is_fitted_on():
      X = np.array([[1, 0, 3], [2, 0, 1], [0, 5, 2], [1, 7, 0]])
      selector = bt.ml.PrevalenceFilter(min_prevalence=0.5).fit(X[:2])
      np.testing.assert_array_equal(selector.prevalence_, [1.0, 0.0, 1.0])
      np.testing.assert_array_equal(selector.transform(X[2:]), [[0, 2], [1, 0]])


  def test_prevalence_filter_boundary_is_inclusive_and_exact():
      X = np.zeros((25, 2))
      X[:7, 0] = 1.0
      X[:, 1] = 1.0
      # 7 / 25 >= 0.28 holds, while 7 >= 0.28 * 25 does not.
      assert bt.ml.PrevalenceFilter(min_prevalence=0.28).fit(X).get_support().tolist() == [True, True]


  def test_prevalence_filter_keeps_csr_sparse():
      out = bt.ml.PrevalenceFilter(min_prevalence=1.0).fit_transform(_toy_x())
      assert sp.issparse(out) and out.format == "csr" and out.shape == (6, 2)


  def test_prevalence_filter_keeps_feature_names():
      tdata = bt.datasets.toy()
      frame = pd.DataFrame(tdata.X.toarray(), index=tdata.obs_names, columns=tdata.var_names)
      selector = bt.ml.PrevalenceFilter(min_prevalence=1.0).set_output(transform="pandas").fit(frame)
      assert selector.transform(frame).columns.tolist() == ["f3", "f4"]


  def test_prevalence_filter_all_zero_feature_is_dropped_and_all_zero_sample_counts():
      X = np.array([[0, 0], [0, 1], [0, 2], [0, 0]])
      selector = bt.ml.PrevalenceFilter(min_prevalence=0.5).fit(X)
      np.testing.assert_array_equal(selector.prevalence_, [0.0, 0.5])
      assert selector.get_support().tolist() == [False, True]


  def test_prevalence_filter_single_sample():
      selector = bt.ml.PrevalenceFilter(min_prevalence=1.0).fit(np.array([[0, 3, 1]]))
      assert selector.get_support().tolist() == [False, True, True]


  @pytest.mark.parametrize("value", [-0.1, 1.5])
  def test_prevalence_filter_out_of_range_raises(value):
      with pytest.raises(ValueError, match=f"min_prevalence must be between 0 and 1, got {value}"):
          bt.ml.PrevalenceFilter(min_prevalence=value).fit(_toy_x())


  def test_prevalence_filter_with_nothing_kept_raises():
      with pytest.raises(ValueError, match=r"no feature is non-zero in at least min_prevalence=0\.6 of the 2 samples"):
          bt.ml.PrevalenceFilter(min_prevalence=0.6).fit(np.array([[1, 0], [0, 0]]))


  def test_prevalence_filter_stays_unfitted_after_a_failed_first_fit():
      selector = bt.ml.PrevalenceFilter(min_prevalence=2.0)
      with pytest.raises(ValueError, match="min_prevalence must be between 0 and 1"):
          selector.fit(np.array([[1, 0], [1, 1]]))
      with pytest.raises(NotFittedError):
          selector.transform(np.array([[3, 4]]))


  def test_prevalence_filter_stays_unfitted_after_a_failed_refit():
      selector = bt.ml.PrevalenceFilter(min_prevalence=0.5).fit(np.array([[1, 0], [1, 1]]))
      selector.set_params(min_prevalence=0.9)
      with pytest.raises(ValueError, match="no feature is non-zero"):
          selector.fit(np.array([[1, 0], [0, 1]]))
      with pytest.raises(NotFittedError):
          selector.transform(np.array([[3, 4]]))


  def test_prevalence_filter_with_nothing_kept_prints_the_threshold_unrounded():
      with pytest.raises(ValueError, match=r"min_prevalence=0\.999 "):
          bt.ml.PrevalenceFilter(min_prevalence=0.999).fit(np.array([[1, 0], [0, 1]]))


  @pytest.mark.parametrize("value", [True, "0.5", None])
  def test_prevalence_filter_wrong_type_raises(value):
      with pytest.raises(TypeError, match="min_prevalence must be a real number"):
          bt.ml.PrevalenceFilter(min_prevalence=value).fit(_toy_x())


  def test_prevalence_filter_keeps_its_input():
      X = sp.csr_matrix((np.array([0.0, 2.0, 3.0, 1.0]), np.array([0, 1, 1, 2]), np.array([0, 2, 4])), shape=(2, 3))
      before = (X.data.copy(), X.indices.copy(), X.indptr.copy())
      bt.ml.PrevalenceFilter(min_prevalence=0.5).fit_transform(X)
      for kept, original in zip((X.data, X.indices, X.indptr), before, strict=True):
          np.testing.assert_array_equal(kept, original)


  def test_clr_equals_pp_clr():
      tdata = bt.datasets.toy()
      out = bt.ml.CLR().fit_transform(tdata.X)
      np.testing.assert_allclose(out, bt.pp.clr(tdata).layers["clr"], rtol=1e-12)


  def test_clr_dense_and_sparse_input_agree():
      X = _toy_x()
      np.testing.assert_allclose(bt.ml.CLR().fit_transform(X.toarray()), bt.ml.CLR().fit_transform(X), rtol=1e-12)


  def test_clr_all_zero_sample_is_all_zero():
      out = bt.ml.CLR().fit_transform(np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]))
      np.testing.assert_array_equal(out[0], 0.0)


  def test_clr_keeps_its_input():
      X = _toy_x()
      before = X.copy()
      bt.ml.CLR(pseudocount=1).fit_transform(X)
      assert (X != before).nnz == 0


  def test_clr_negative_value_raises_in_fit():
      with pytest.raises(ValueError, match="Negative values in data passed to CLR"):
          bt.ml.CLR().fit(np.array([[-1.0, 2.0]]))


  @pytest.mark.parametrize(("pseudocount", "error"), [(-1, ValueError), (True, TypeError)])
  def test_clr_bad_pseudocount_raises_in_fit(pseudocount, error):
      with pytest.raises(error, match="pseudocount must be"):
          bt.ml.CLR(pseudocount=pseudocount).fit(_toy_x())


  def test_clr_pseudocount_above_the_smallest_value_warns():
      relative = np.array([[0.2, 0.8], [0.5, 0.5]])
      with pytest.warns(UserWarning, match=r"pseudocount=0.5 is larger than the smallest non-zero value in X \(0.2\)"):
          bt.ml.CLR().fit_transform(relative)


  def test_clr_keeps_feature_names():
      tdata = bt.datasets.toy()
      frame = pd.DataFrame(tdata.X.toarray(), index=tdata.obs_names, columns=tdata.var_names)
      out = bt.ml.CLR().set_output(transform="pandas").fit_transform(frame)
      assert out.columns.tolist() == tdata.var_names.tolist() and out.index.tolist() == tdata.obs_names.tolist()


  def test_a_pipeline_refits_the_filter_in_every_fold():
      tdata = bt.datasets.toy()
      pipeline = make_pipeline(bt.ml.PrevalenceFilter(min_prevalence=0.9), bt.ml.CLR(), LogisticRegression())
      result = cross_validate(
          pipeline, tdata.X, tdata.obs["group"], cv=StratifiedKFold(3), return_estimator=True, return_indices=True
      )
      for estimator, train in zip(result["estimator"], result["indices"]["train"], strict=True):
          own = bt.ml.PrevalenceFilter(min_prevalence=0.9).fit(tdata.X[train])
          np.testing.assert_array_equal(estimator[0].prevalence_, own.prevalence_)
      supports = {tuple(estimator[0].get_support()) for estimator in result["estimator"]}
      assert len(supports) > 1


  @given(
      arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 6)), elements=st.integers(0, 3)),
      st.floats(0, 1),
  )
  def test_prevalence_filter_keeps_exactly_the_features_at_or_above_the_threshold(X, threshold):
      present = (X != 0).sum(axis=0) / X.shape[0]
      if not (present >= threshold).any():
          with pytest.raises(ValueError, match="no feature is non-zero"):
              bt.ml.PrevalenceFilter(min_prevalence=threshold).fit(X)
          return
      support = bt.ml.PrevalenceFilter(min_prevalence=threshold).fit(X).get_support()
      np.testing.assert_array_equal(support, present >= threshold)


  @given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 1000)))
  def test_clr_rows_sum_to_zero(X):
      np.testing.assert_allclose(bt.ml.CLR().fit_transform(X).sum(axis=1), 0.0, atol=1e-9)


  @pytest.mark.parametrize(("estimator", "value"), [(bt.ml.PrevalenceFilter, 0.8), (bt.ml.CLR, 1)])
  def test_transformer_options_are_keyword_only(estimator, value):
      with pytest.raises(TypeError):
          estimator(value)
  ```
  In `tests/test_docstrings.py` the module set becomes
  `{"da", "datasets", "fn", "io", "ml", "pl", "pp", "tl"}`.
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/ml -q` -> `1 error` during collection
  (`AttributeError: module 'biotapy' has no attribute 'ml'`).
- [x] **Step 3: Implement.** Create `src/biotapy/ml/_transformers.py`:
  ```python
  """scikit-learn transformers, so preprocessing is fitted inside each cross-validation fold."""

  from typing import Self

  import numpy as np
  import numpy.typing as npt
  import scipy.sparse as sp
  from skbio.stats.composition import clr as skbio_clr
  from sklearn.base import BaseEstimator, OneToOneFeatureMixin, TransformerMixin
  from sklearn.feature_selection import SelectorMixin
  from sklearn.utils import Tags
  from sklearn.utils.validation import check_is_fitted, check_non_negative, validate_data

  from biotapy._core import as_csr, check_pseudocount, pseudocounted

  # What fit and transform take: a samples x features array, sparse matrix or DataFrame.
  Table = npt.ArrayLike | sp.spmatrix


  class PrevalenceFilter(SelectorMixin, BaseEstimator):
      """Keep the features that are non-zero in enough of the training samples.

      Parameters
      ----------
      min_prevalence
          Keep features that are non-zero in at least this fraction of the samples
          given to ``fit``, from 0 to 1.

      Attributes
      ----------
      prevalence_
          Each feature's fraction of training samples in which it is non-zero.
      n_features_in_
          Number of features seen in ``fit``.
      feature_names_in_
          Feature names seen in ``fit``, when it was given a DataFrame.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/machine_learning`

      The scikit-learn form of :func:`biotapy.pp.filter_features` with
      ``min_prevalence``, keeping the same features. Which features pass depends
      on the samples, so filtering before splitting lets the test samples choose
      the features: inside a :class:`~sklearn.pipeline.Pipeline` the filter is
      fitted on each training fold only. A sparse ``X`` stays sparse.
      ``transform``, ``fit_transform``, ``get_support`` and
      ``get_feature_names_out`` come from
      :class:`~sklearn.feature_selection.SelectorMixin`.

      Examples
      --------
      >>> import biotapy as bt
      >>> X = bt.datasets.toy().X
      >>> bt.ml.PrevalenceFilter(min_prevalence=1.0).fit_transform(X).shape
      (6, 2)
      """

      def __init__(self, *, min_prevalence: float = 0.1) -> None:
          self.min_prevalence = min_prevalence

      def fit(self, X: Table, y: object = None) -> Self:
          """Learn each feature's prevalence in ``X``, samples x features."""
          # validate_data sets n_features_in_ before the checks below can raise, so prevalence_ alone marks a fit.
          self.__dict__.pop("prevalence_", None)
          X = validate_data(self, X, accept_sparse="csr")
          if isinstance(self.min_prevalence, bool) or not isinstance(
              self.min_prevalence, int | float | np.integer | np.floating
          ):
              msg = f"min_prevalence must be a real number, got {self.min_prevalence!r}"
              raise TypeError(msg)
          if not 0 <= self.min_prevalence <= 1:
              msg = f"min_prevalence must be between 0 and 1, got {self.min_prevalence}"
              raise ValueError(msg)
          matrix = as_csr(X)
          present = np.bincount(matrix.indices[matrix.data != 0], minlength=matrix.shape[1])
          prevalence = present / matrix.shape[0]
          if not (prevalence >= self.min_prevalence).any():
              msg = (
                  f"no feature is non-zero in at least min_prevalence={self.min_prevalence} of the "
                  f"{matrix.shape[0]} samples; lower min_prevalence"
              )
              raise ValueError(msg)
          self.prevalence_ = prevalence
          return self

      def _get_support_mask(self) -> npt.NDArray[np.bool_]:
          check_is_fitted(self, "prevalence_")
          # Divide rather than multiply, as pp.filter_features does: 7 / 25 >= 0.28 holds, 7 >= 0.28 * 25 does not.
          return np.asarray(self.prevalence_ >= self.min_prevalence)

      def __sklearn_tags__(self) -> Tags:
          # scikit-learn's tag methods are unannotated, and mypy's untyped_calls_exclude does not reach super().
          tags: Tags = super().__sklearn_tags__()  # type: ignore[no-untyped-call]
          tags.input_tags.sparse = True
          return tags


  class CLR(OneToOneFeatureMixin, TransformerMixin, BaseEstimator):
      """Centred log-ratio transform of each sample, as a scikit-learn transformer.

      Parameters
      ----------
      pseudocount
          Added to every value before the logarithm, so zeros have one.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/machine_learning`

      Gives the values of :func:`biotapy.pp.clr`. Each sample is transformed on
      its own, so ``fit`` learns nothing (it only checks ``X``); the class exists
      so the transform can sit in a :class:`~sklearn.pipeline.Pipeline`. The
      output is a dense float64 array: 8 bytes x samples x features.
      ``fit_transform``, ``get_feature_names_out`` and ``set_output`` come from
      scikit-learn's :class:`~sklearn.base.TransformerMixin` and
      :class:`~sklearn.base.OneToOneFeatureMixin`.

      Examples
      --------
      >>> import biotapy as bt
      >>> out = bt.ml.CLR().fit_transform(bt.datasets.toy().X)
      >>> round(float(out[0].sum()), 6)
      0.0
      """

      def __init__(self, *, pseudocount: float = 0.5) -> None:
          self.pseudocount = pseudocount

      def fit(self, X: Table, y: object = None) -> Self:
          """Check ``X`` (samples x features, non-negative) and ``pseudocount``; nothing is learned."""
          X = validate_data(self, X, accept_sparse="csr")
          check_pseudocount(self.pseudocount)
          check_non_negative(X, "CLR")
          return self

      def transform(self, X: Table) -> npt.NDArray[np.float64]:
          """Return the CLR of each sample of ``X`` as a dense float64 array."""
          check_is_fitted(self)
          X = validate_data(self, X, accept_sparse="csr", reset=False)
          return np.asarray(skbio_clr(pseudocounted(X, self.pseudocount, func="ml.CLR")), dtype=np.float64)

      def __sklearn_tags__(self) -> Tags:
          # scikit-learn's tag methods are unannotated, and mypy's untyped_calls_exclude does not reach super().
          tags: Tags = super().__sklearn_tags__()  # type: ignore[no-untyped-call]
          tags.input_tags.sparse = True
          tags.input_tags.positive_only = True
          tags.requires_fit = False
          return tags
  ```
  `src/biotapy/ml/__init__.py`:
  ```python
  from ._transformers import CLR, PrevalenceFilter

  __all__ = ["CLR", "PrevalenceFilter"]
  ```
  `src/biotapy/__init__.py`:
  ```python
  from importlib.metadata import version

  from . import da, datasets, fn, io, ml, pl, pp, tl

  __all__ = ["__version__", "da", "datasets", "fn", "io", "ml", "pl", "pp", "tl"]

  __version__ = version("biotapy")
  ```
- [x] **Step 4: Run, expect pass** -
  `uv run --group test pytest tests/ml src/biotapy/ml tests/test_docstrings.py -q -W error::UserWarning`
  -> `181 passed, 2 skipped`; `tests/ml` also under
  `--hypothesis-seed=1`, `2`, `3` (`123 passed, 2 skipped` each);
  `uv run --group dev --group doc mypy` -> `Success`.
- [x] **Step 5: Docs.** Create `docs/guide/machine_learning.md`:
  ````markdown
  # Machine learning

  `bt.ml` holds scikit-learn transformers, so microbiome preprocessing can sit
  inside a [Pipeline](https://scikit-learn.org/stable/modules/compose.html) and
  be fitted on the training samples of each cross-validation fold only.

  ## Which steps leak

  A step leaks when what it does to one sample depends on other samples. Run
  before the data is split, it lets the test samples shape the training data,
  and cross-validated scores come out better than they will be on new samples.

  | Step | Learns from the samples? | In a pipeline |
  |---|---|---|
  | Prevalence filter (`bt.ml.PrevalenceFilter`) | yes: which features are common enough | must be inside |
  | Relative abundance (`sklearn.preprocessing.Normalizer(norm="l1")`) | no: each sample on its own | either way |
  | CLR (`bt.ml.CLR`) | no: each sample on its own | either way |
  | Scaling (`sklearn.preprocessing.StandardScaler`) | yes: each feature's mean and spread | must be inside |

  `bt.pp.filter_features` on the whole table before cross-validation is the
  leaky version of `bt.ml.PrevalenceFilter`: the same rule, fitted on every
  sample at once.

  ## A leak-free pipeline

  ```python
  import biotapy as bt
  from sklearn.linear_model import LogisticRegression
  from sklearn.model_selection import cross_validate
  from sklearn.pipeline import make_pipeline

  tdata = bt.datasets.toy()
  pipeline = make_pipeline(
      bt.ml.PrevalenceFilter(min_prevalence=0.1),
      bt.ml.CLR(pseudocount=0.5),
      LogisticRegression(),
  )
  scores = cross_validate(pipeline, tdata.X, tdata.obs["group"], cv=3)
  ```

  The transformers take what scikit-learn takes: a samples x features array, a
  sparse matrix (`tdata.X`) or a DataFrame, whose column names they keep with
  `set_output(transform="pandas")`. They do not take an AnnData.

  `PrevalenceFilter` keeps a sparse `X` sparse. `CLR` returns a dense array,
  because a log-ratio has no zeros, and gives the values of `bt.pp.clr`; like
  `bt.pp.clr` it warns when the pseudocount is larger than the smallest non-zero
  value, as when the default 0.5 meets relative abundances.

  ## Relative abundance

  scikit-learn already has it: `Normalizer(norm="l1")` divides each sample by
  its total, keeps a sparse matrix sparse and leaves an all-zero sample at zero,
  as `bt.pp.relative` does. Before `CLR` it changes nothing but the scale the
  pseudocount is on.
  ````
  Create `docs/_templates/autosummary/class.rst`:
  ```rst
  {{ fullname | escape | underline }}

  .. currentmodule:: {{ module }}

  .. autoclass:: {{ objname }}
     :members:
  ```
  and:
  ```diff
  diff --git a/docs/api.md b/docs/api.md
  index 503bd3d..54ea9ed 100644
  --- a/docs/api.md
  +++ b/docs/api.md
  @@ -109,6 +109,19 @@ Public functions are listed here as they ship, from Phase 1 onward.
       da.maaslin3
   ```

  +## Machine learning
  +
  +```{eval-rst}
  +.. module:: biotapy.ml
  +.. currentmodule:: biotapy
  +
  +.. autosummary::
  +    :toctree: generated
  +
  +    ml.CLR
  +    ml.PrevalenceFilter
  +```
  +
   ## Plots

   ```{eval-rst}
  diff --git a/docs/conf.py b/docs/conf.py
  index 606b1e4..6fbcf48 100644
  --- a/docs/conf.py
  +++ b/docs/conf.py
  @@ -111,6 +111,7 @@ intersphinx_mapping = {
       "networkx": ("https://networkx.org/documentation/stable/", None),
       "matplotlib": ("https://matplotlib.org/stable/", None),
       "mudata": ("https://mudata.scverse.org/stable/", None),
  +    "sklearn": ("https://scikit-learn.org/stable/", None),
   }

   # List of patterns, relative to source directory, that match files and
  diff --git a/docs/guide/index.md b/docs/guide/index.md
  index 391d35a..5132c6a 100644
  --- a/docs/guide/index.md
  +++ b/docs/guide/index.md
  @@ -16,5 +16,6 @@ filtering
   diversity
   ordination
   differential_abundance
  +machine_learning
   plotting
   ```
  ```
- [x] **Step 6: Knowledge.**
  ```diff
  diff --git a/.knowledge/contracts/function-shape.md b/.knowledge/contracts/function-shape.md
  index 302816e..162ef54 100644
  --- a/.knowledge/contracts/function-shape.md
  +++ b/.knowledge/contracts/function-shape.md
  @@ -24,6 +24,10 @@ Every public function in `io`, `datasets`, `pp`, `tl`, `fn`, `da`, `ml`, `pl`:
        `da.consensus` takes a sequence of `da` result tables and `pl.consensus` the
        table `da.consensus` returns instead: they combine and draw results, not data.
      - Everything after the required arguments is keyword-only (`*`).
  +   - `ml`'s scikit-learn transformers are classes (rules.md R3.6), not
  +     functions: their options are constructor keyword arguments that `fit`
  +     validates (scikit-learn's convention, which `check_estimator` tests),
  +     and the class docstring carries the skeleton below.
      - No `**kwargs` pass-through, except a documented `plot_kwargs` in `pl`.
   2. **Return and mutation**: per [pure-by-default](/decisions/pure-by-default.md).
   3. **Randomness**: any stochastic function takes
  diff --git a/.knowledge/decisions/pure-by-default.md b/.knowledge/decisions/pure-by-default.md
  index 6692999..bd10470 100644
  --- a/.knowledge/decisions/pure-by-default.md
  +++ b/.knowledge/decisions/pure-by-default.md
  @@ -29,6 +29,7 @@ behaviour. Confirmed by the user on 2026-09-26.
   | `tl` with `inplace=True` | `None` | writes to the slot named in [data-model-slots](/contracts/data-model-slots.md) |
   | `fn`, `da` | new object or result `pd.DataFrame` | never |
   | `pl` | `matplotlib.axes.Axes` | never |
  +| `ml` estimators (`PrevalenceFilter`, `CLR`) | scikit-learn's protocol: `fit` stores what it learns on the estimator and returns it; `transform` returns a new array | never the data |

   Every `tl` function that returns per-sample or per-pair values supports both modes with
   identical semantics; `tl.permanova` and `tl.mmvec` are the exceptions (see Consequences).
  ```
  Tick 4.3 here; log line:
  ```markdown
  - **Update**: [pure-by-default](decisions/pure-by-default.md) gains the `ml` estimators' row; [function-shape](contracts/function-shape.md) says `ml`'s transformers are classes whose options `fit` validates; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.3.
  ```
- [x] **Step 7: Gate and commit**
  ```bash
  git add src/biotapy/ml/__init__.py src/biotapy/ml/_transformers.py src/biotapy/__init__.py \
    tests/ml/test_transformers.py tests/test_docstrings.py docs/guide/machine_learning.md \
    docs/guide/index.md docs/api.md docs/conf.py docs/_templates/autosummary/class.rst \
    .knowledge/decisions/pure-by-default.md .knowledge/contracts/function-shape.md \
    .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html
  uv run --group test pytest -q -W error::UserWarning
  uv run --group test pytest -q -m "golden or network"
  git commit -m "feat(ml): add PrevalenceFilter and CLR scikit-learn transformers

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `build succeeded`;
  `1410 passed, 2 skipped, 54 deselected`; `36 passed`.

### Checkpoint A - review slice 4A
- [x] Review the whole slice (superpowers:requesting-code-review) against
  every contract, pure-by-default, the Phase 4 and slice 4A review focus;
  then a fix pass, one commit per finding, each with a test. Record the
  counts and the fix range here.
  - Review: 0 Critical / 0 Important / 5 Minor (+3 nits). Fix pass
    `6ea27aa..9883786` (6 commits). Scoped re-review: all 6 addressed, no new
    findings.
- [x] Run the slice's checks: `uv run --group test pytest -q -W error::UserWarning`
  and `uv run --group test pytest -q -m "golden or network"`; record both
  counts at the last commit; `coverage run -m pytest` then `coverage
  report` for the six slice files (100% each on the prototype).
  - At `9883786`: `pytest -q -W error::UserWarning` 1410 passed, 2 skipped, 54
    deselected; `-m "golden or network"` 36 passed.
  - `coverage run -m pytest tests/io/test_mudata.py tests/tl/test_mmvec.py
    tests/core/test_composition.py tests/ml` (162 passed, 2 skipped), then
    `coverage report` on `io/_mudata.py`, `tl/_mmvec.py`,
    `_core/_composition.py`, `ml/_transformers.py`, `ml/__init__.py` and
    `_core/__init__.py`: 100% each (161 statements, 0 missed).
- [x] Knowledge (codebase-map templates; R12.2-R12.4):
  - **Create `.knowledge/modules/ml.md`** (`type: Module`, `paths:
    ["src/biotapy/ml/**"]`): Responsibility (scikit-learn transformers now;
    `to_torch` and `embed` from 4B/4C); Entry points
    `_transformers.py:PrevalenceFilter`, `_transformers.py:CLR`; Invariants:
    validation in `fit`, `PrevalenceFilter`'s rule equals
    `pp.filter_features`', `CLR` equals `pp.clr`, CSR kept sparse, no
    AnnData input; Gotchas: `check_estimator` cannot run under `-W
    error::UserWarning` (use `parametrize_with_checks`), the checks need
    the pseudocount warning filtered, `super().__sklearn_tags__()` needs a
    targeted type ignore, inherited methods stay out of the API pages
    (class template); Verification: `uv run --group test pytest tests/ml`.
    Add it to `modules/index.md` and the `.knowledge/index.md` Modules line.
  - **Update `.knowledge/modules/io.md`**: `to_mudata` (entry point,
    intersection rule, copies, no provenance and why).
  - **Update `.knowledge/modules/tl.md`**: `mmvec` (entry point, the
    pre-checks, dense inputs, no `inplace`).
  - Log lines for each.
- [x] Push the branch and open the PR only after the user approves that push
  (R13.3). CI green, including docs and the network job.
  - PR #27, Test run 37872650983 green (21 jobs, including docs, network and
    r-bridge); merged as `cdc3b07`.
- [x] Ask the user to review slice 4A, and to confirm
  [multiomics-as-mudata](/decisions/multiomics-as-mudata.md) (it stays
  `draft` until then), before slice 4B is expanded.
  - Approved 2026-10-09; the decision is `stable`.

## Fixes found in slice 4A (4.F1-4.F3)

Found during slice 4A's reviews; approved by the user on 2026-10-09 as one
small PR (branch `fix-4f`) before slice 4B. Each is one commit with its test.

### Task 4.F1: `da.ancombc2` on a two-sample reference group (scikit-bio defect)

**Root cause** (traced 2026-10-09 in scikit-bio and in R): when the reference group
has two samples whose centred log abundances of one feature nearly agree, the
intercept's `var_hat` is tiny (3.3e-6 on the pinned table). Nelder-Mead pushes
`kappa1`/`kappa2` to their bound 0 and all three Gaussian densities underflow,
so scikit-bio 0.7.4's `_estimate_bias_em` computes `resp /= resp.sum(0)` as
0/0 = NaN; `pi` turns NaN and `np.quantile` raises "Quantiles must be in the
range [0, 1]". R's `.bias_em` sets those responsibilities to 0
(`r0i[is.na(r0i)] = 0`), a line scikit-bio never ported; R fits the table
(log2 B vs A: f0 2.002, f1 0.744, f2 0.088, f3 -1.315, f4 -0.218, f5 -0.325,
f6 -0.605), and scikit-bio with R's line matches R to 6 decimals on all 16
failing tables found by sweep (one to 0.011, both runs at the E-M cap).
Unfixed on scikit-bio `main`; reported as [scikit-bio#2631](https://github.com/scikit-bio/scikit-bio/issues/2631). Rate: 0.16% of 6 x 7 tables
with a 2-sample reference group, 0 with 3 or more; the property test fails in
about 1.2% of runs, then every run once Hypothesis has saved the example.

biotapy already reports it clearly (`ValueError: da.ancombc2: scikit-bio could
not fit the model: ...`, R7.4). A biotapy-side fix would mean monkeypatching
scikit-bio's private function or reimplementing the E-M, so the fix is: pin the
behaviour, restrict the property test to fittable designs for this one defect,
document it, and report it upstream.

**Files:** modify `tests/da/test_ancombc.py`, `tests/da/test_schema.py`,
`src/biotapy/da/_ancombc.py` (Notes only).

- [x] **Step 1: regression test (pins current behaviour).** Append to
  `tests/da/test_ancombc.py`:
  ```python
  # Two reference samples whose centred log abundances of f5 nearly agree (variance 3.3e-6): scikit-bio 0.7.4's bias
  # E-M divides 0 by 0 where R's .bias_em sets the responsibility to 0. R fits this table (B vs A, log2): f0 2.002,
  # f1 0.744, f2 0.088, f3 -1.315, f4 -0.218, f5 -0.325, f6 -0.605. When scikit-bio fixes it, this test fails:
  # turn it into a check against those values.
  EM_UNDERFLOW = [
      [11, 39, 115, 65, 122, 90, 12],
      [81, 41, 194, 22, 174, 192, 186],
      [86, 68, 125, 1, 103, 14, 42],
      [21, 24, 43, 20, 38, 80, 6],
      [134, 61, 155, 8, 84, 169, 31],
      [71, 17, 64, 28, 63, 0, 10],
  ]


  def test_two_reference_samples_that_underflow_scikit_bio_raise_naming_it():
      counts = np.array(EM_UNDERFLOW)
      adata = ad.AnnData(
          X=sp.csr_matrix(counts),
          obs=pd.DataFrame({"g": ["a"] * 2 + ["b"] * 4}, index=[f"s{i}" for i in range(6)]),
          var=pd.DataFrame(index=[f"f{i}" for i in range(7)]),
      )
      with pytest.raises(ValueError, match=r"da\.ancombc2: scikit-bio could not fit the model: Quantiles must be in"):
          bt.da.ancombc2(adata, "g")
  ```
  (Use the file's existing imports; add only the missing ones.) Run:
  `uv run --group test pytest -q tests/da/test_ancombc.py -k underflow -W error::UserWarning`.
  Expected: `1 passed` (it pins existing, documented behaviour; R11.1's RED is
  Step 2).
- [x] **Step 2: RED.** In `tests/da/test_schema.py`, import `example` and
  `reject` from hypothesis and add, under the property test's `@given`:
  ```python
  @example(counts=np.array(EM_UNDERFLOW), split=2)
  ```
  with `EM_UNDERFLOW` copied as a module constant (tests do not import each
  other). Run `uv run --group test pytest -q tests/da/test_schema.py -k direction -W error::UserWarning`.
  Expected: `1 failed, 2 passed` (ancombc2 fails on the example).
- [x] **Step 3: GREEN.** Replace `out = method(adata, "g")` in that test with:
  ```python
  try:
      out = method(adata, "g")
  except ValueError as err:
      # Only scikit-bio's bias E-M underflow, pinned in test_ancombc.py; any other error still fails.
      if method is bt.da.ancombc2 and "Quantiles must be in the range [0, 1]" in str(err):
          reject()
      raise
  ```
  Expected: `3 passed`. This weakens no assertion (R11.5): the rejection covers
  one method and one upstream message, the excluded input is asserted by the
  Step 1 test, and it extends the test's existing restriction to fittable
  designs (`unique=True` excludes exact zero variance).
- [x] **Step 4: Notes.** In `da.ancombc2`'s docstring Notes, after the E-M
  sentence, add: "With two samples in the reference level, scikit-bio 0.7.4 can
  fail to estimate the bias (about 1 in 600 random small tables) where R
  returns results; biotapy raises naming scikit-bio rather than guess." Docs
  build with `-W`.
- [x] **Step 5: commit** `test(da): pin scikit-bio's ancombc2 bias E-M underflow and keep the schema property on fittable designs`.
- [x] **Upstream (needs separate approval, R13.3):** an issue on
  scikit-bio/scikit-bio; its draft text is in the 4.F PR description.
  Filed as [scikit-bio#2631](https://github.com/scikit-bio/scikit-bio/issues/2631).

### Task 4.F2: `pp.filter_features` names a wrongly typed threshold

**Files:** `src/biotapy/pp/_filter.py`, `tests/pp/test_filter.py`.
Same rule as `ml.PrevalenceFilter` (R3.5): a bool or a non-real
`min_prevalence` or `min_total` raises `TypeError` naming the argument.

- [x] **Step 1: RED.**
  ```python
  @pytest.mark.parametrize("value", [True, "0.5", [0.5]])
  @pytest.mark.parametrize("argument", ["min_prevalence", "min_total"])
  def test_wrongly_typed_threshold_raises_naming_it(argument, value):
      with pytest.raises(TypeError, match=f"{argument} must be a real number"):
          bt.pp.filter_features(bt.datasets.toy(), **{argument: value})
  ```
  Expected: 6 failed (`True` passes silently; the others raise an unnamed
  `TypeError` or numpy error).
- [x] **Step 2: GREEN.** Before the range checks, for each non-None
  threshold: `if isinstance(value, bool) or not isinstance(value, int | float | np.integer | np.floating): raise TypeError(f"{name} must be a real number, got {value!r}")`,
  written once as a loop over the two names (R5 limits). Same wording as
  `ml.PrevalenceFilter`. Expected: tests/pp passes.
- [x] **Step 3: commit** `fix(pp): name min_prevalence or min_total when its type is wrong`.

### Task 4.F3: one finite, non-negative check in `_core`

The predicate `np.isfinite(X.data).all() and not (X.data < 0).any()` is
written in `_core/_composition.py`, `tl/_mmvec.py` and `fn/_redundancy.py`
(R4.3). Move it to `_core/_matrix.py` as a predicate, so each caller keeps its
own message (no user-visible change).

**Files:** `src/biotapy/_core/_matrix.py`, `src/biotapy/_core/__init__.py`,
the three callers, `tests/core/test_matrix.py` (or the existing `_core`
matrix test file).

- [x] **Step 1: RED.**
  ```python
  @pytest.mark.parametrize(
      ("data", "expected"),
      [([0.0, 1.5], True), ([], True), ([1.0, -0.1], False), ([1.0, np.nan], False), ([np.inf], False)],
  )
  def test_finite_non_negative(data, expected):
      X = sp.csr_matrix((np.array(data, dtype=float), (np.zeros(len(data), int), np.arange(len(data)))), shape=(1, max(len(data), 1)))
      assert finite_non_negative(X) is expected
  ```
  Expected: ImportError.
- [x] **Step 2: GREEN.**
  ```python
  def finite_non_negative(X: sp.csr_matrix) -> bool:
      """True when every stored value of ``X`` is finite and >= 0."""
      return bool(np.isfinite(X.data).all() and not (X.data < 0).any())
  ```
  Export it from `_core/__init__.py`; replace the three inline predicates.
  `tests/pp`, `tests/tl`, `tests/fn`, `tests/core`, `tests/ml` pass unedited.
- [x] **Step 3:** update `modules/core.md` (entry point) and log line.
  Commit `refactor(core): share the finite, non-negative check`.

### Gates and merge
Done: PR #28, Test run 37881255064 green (21 jobs); merged as `ef82843`.
Commit, `git status --short` empty, then: prek; `pytest -q -W error::UserWarning`
(1412 + 1 + 6 + 5 passed, 2 skipped, 54 deselected); `-m "golden or network"` 36;
docs `-W`; `knowledge_stale.sh` 0 stale. Push `fix-4f`, PR, merge commit on
green (decision 20 covers Phase 4 branches).

---

## Slice 4B - torch

**Goal:** with `pip install 'biotapy[torch]'`, a user turns a table into a
PyTorch dataset whose rows stay sparse until they are read
(`bt.ml.to_torch`), and CI runs the PyTorch tests on CPU wheels in a job of
their own that blocks merges.

### Slice 4B design
- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `pyproject.toml` | the extra `torch`, `[tool.uv]` index and source (4.B0); the mypy override, the marker `torch`, `addopts` (4.5) |
  | `ml/_torch.py` | `to_torch`, `_table`, `_labels` (4.5) |
  | `conftest.py` (root) | the hook that gives `biotapy.ml._torch`'s doctests the marker `torch` (4.5) |
  | `tests/ml/test_torch.py` | 22 tests marked `torch`, 6 that need no torch (4.5) |
  | `.github/workflows/test.yaml` | the job `ml-extras`, in `check.needs` (4.B1) |
  | `docs/guide/machine_learning.md`, `docs/contributing.md` | "PyTorch" guide section; "PyTorch tests" (4.5, 4.B1) |

  `uv.lock` is git-ignored (`/uv.lock`, "resolve fresh in CI"), so no task
  commits it.
- **How slice 4B was checked.** Every file below was written into a scratch
  clone at `cdc3b07` (master after slice 4A) on branch `phase-4b` and
  committed one task at a time, after a stand-in `docs(roadmap)` commit that
  only adds the 4.B0 checklist line; the branch was then rebased onto
  `ef82843` (master after the 4.F fixes, which touch none of these files)
  and every commit gated again. The counts below are the rebased ones
  (`79c6b44` stand-in, `fada00b`, `875110f`, `5350e0e`). The gates ran on each committed tree
  (`git status --short` empty) with `.venv` holding no torch, as in CI's
  default jobs; the `-m torch` run used a second environment with the extra
  (`torch 2.14.1+cpu`, Python 3.13.2). Every run exported `BIOTAPY_DATA_DIR`
  to a scratch pooch cache; `~/.cache/biotapy` was absent after each.

  | Task state | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` | `--extra torch pytest -q -m torch -W error::UserWarning` |
  |---|---|---|---|---|---|
  | base `ef82843` | passed (14 hooks) | 1424 passed, 2 skipped, 54 deselected | 36 passed | build succeeded | - |
  | 4.B0 `build` (`fada00b`) | passed | 1425 passed, 2 skipped, 54 deselected | 36 passed | build succeeded | 0 selected, 1481 deselected (no marker yet) |
  | 4.5 `feat(ml)` (`875110f`) | passed | 1432 passed, 2 skipped, 77 deselected | 36 passed | build succeeded | 23 passed, 1488 deselected |
  | 4.B1 `ci` (`5350e0e`) | passed | 1434 passed, 2 skipped, 77 deselected | 36 passed | build succeeded | 23 passed, 1490 deselected |

  - The default run on `ef82843` and after also ends in "2 warnings" or "3
    warnings": scikit-bio's `RuntimeWarning: invalid value encountered in
    divide` from 4.F1's `ancombc2` tests (the count follows Hypothesis's
    draws), not from 4B. On `cdc3b07` the same commits gave 1413, 1420 and
    1422 passed with the same skips and deselections.
  - The 23 deselected by default in 4.5 are the 22 `torch` tests and
    `to_torch`'s doctest; the 7 new passes are the 6 tests that need no
    torch and `test_docstring_has_the_contract_sections[bt.ml.to_torch]`.
  - At the 4.5 commit (on `cdc3b07`), with torch installed in the
    environment, the default run gives the same `1420 passed, 2 skipped, 77
    deselected` as without it; `mypy --strict`
    passes with and without torch (`Success: no issues found in 65 source
    files`). Without the extra, `import biotapy` never imports torch: the
    `import-without-extras` command imports every module with `'torch' in
    sys.modules` False. With torch installed, scikit-bio 0.7.4 imports it
    (`skbio.util._testing` at module level), so the check means something only
    in a torch-free environment.
  - CI's `test` job, replayed through hatch at the 4.B1 commit (on `cdc3b07`)
    (`hatch run hatch-test.py3.13-stable:run-cov -n auto`): `1422 passed, 2
    skipped`, coverage total 99%; hatch installed no torch.
  - Coverage of `ml/_torch.py`: 100% (54 statements) with the extra
    (`coverage run -m pytest -m "torch or not torch" tests/ml/test_torch.py
    src/biotapy/ml/_torch.py`, 32 passed); 71% without it (the class and
    the import are not reached).
  - After the Checkpoint B fix pass (on `d33d661`, three tests added to the
    default run and two to `-m torch`): `1437 passed, 2 skipped, 79
    deselected`; `-m "golden or network"` `36 passed, 1482 deselected`;
    `-m torch` `25 passed, 1493 deselected` (24 tests and the doctest).
  - Hypothesis seeds 1, 2 and 3: `22 passed, 6 deselected` each.
  - Mutation check: densifying the whole table in `to_torch` fails
    `test_a_large_sparse_table_is_never_dense` (`assert peak < 50_000_000`)
    and `test_references_x_so_a_later_change_shows`.
  - The `ml-extras` job replayed locally with a cold uv cache and a new
    environment (Python 3.13): the sync step 19.4 s, the test step 14.8 s
    (pytest 6.2 s), environment 1.6 GB (torch 727 MB).
- **Resolved: mypy and the in-function `Dataset` subclass** (design note 7's
  [UNVERIFIED]). prek's mypy hook runs `uv run --group dev --group doc mypy`,
  an inexact sync of `.venv`: CI's `lint` job has no torch, while a
  contributor who ran `--extra torch` keeps it. Measured with mypy 2.4.0:

  | Code | torch absent | torch 2.14.1+cpu present |
  |---|---|---|
  | `torch = import_optional(...)`; `class _Rows(torch.utils.data.Dataset)` | `Name "torch.utils.data.Dataset" is not defined [name-defined]` | same |
  | no override; `TYPE_CHECKING` import of `Dataset` | `Cannot find implementation or library stub for module named "torch" [import-not-found]` | `Missing type arguments for generic type "Dataset" [type-arg]` |
  | `ignore_missing_imports` only; typed `Dataset["Tensor"]` | passes only with `# type: ignore[misc, unused-ignore]` and `[no-any-return, unused-ignore]`: the two environments need different ignores | same ignores; mypy analyses torch (20.4 s against 9.6 s on the same cache state) |
  | **chosen:** `torch: Any = import_optional("torch", extra="torch")`, `class AnnDataDataset(torch.utils.data.Dataset):  # type: ignore[misc]`, override `{ module = "torch" / "torch.*", ignore_missing_imports = true, follow_imports = "skip" }` | `Success` (cold 29 s) | `Success` (cold 25 s: torch is not analysed) |

  `follow_imports = "skip"` makes every torch name `Any` in both
  environments, so one `# type: ignore[misc]` (subclassing `Any`) is used in
  both and `warn_unused_ignores` never fires. R4.6 holds: torch is reached
  only through `import_optional` inside `to_torch`. R3.6 holds: one level of
  inheritance from `torch.utils.data.Dataset`, checked by
  `test_is_a_torch_dataset_with_one_item_per_sample`.
- **Resolved: the return annotation.** sphinx-autodoc-typehints 3.13.9,
  without torch in the docs environment, renders `-> "Dataset[Tensor]"` (and
  a union of them, and `"torch.utils.data.Dataset[torch.Tensor]"`) as the
  broken reference `torch.utils.data.Dataset.torch.Tensor`, which
  `nitpicky` fails. The bare `-> "Dataset"` with `if TYPE_CHECKING: from
  torch.utils.data import Dataset` renders as `Dataset` linked to
  `https://docs.pytorch.org/docs/stable/data.html#torch.utils.data.Dataset`
  through a new intersphinx entry (PyTorch's inventory has
  `torch.utils.data.Dataset`, `DataLoader` and `torch.Tensor`). With
  `follow_imports = "skip"` the bare name raises no `type-arg`.
- **Resolved: `uv lock` across platforms** (design note 7's [UNVERIFIED]).
  With `[tool.uv] sources.torch = { index = "pytorch-cpu" }` and the index
  `explicit = true`, `uv lock` (uv 0.6.13, and 0.12.24 gives the same
  versions) resolves 195 packages instead of 189: it adds fsspec 2026.9.0,
  mpmath 1.3.0, setuptools 84.0.0, sympy 1.14.0 and two torch entries, and
  moves no other version (lock diff +139/-9 lines). torch forks by platform:
  `2.14.1+cpu` (wheels for Linux x86_64, aarch64, s390x and Windows amd64,
  arm64; cp312, cp313, cp314, cp314t) and `2.14.1` for macOS (only
  `macosx_14_0_arm64`). `uv pip compile --extra torch --python-platform`:
  Windows py3.12 and py3.14 -> `torch==2.14.1+cpu`; Linux aarch64 ->
  `2.14.1+cpu`; macOS arm64 with `MACOSX_DEPLOYMENT_TARGET=14.0` ->
  `2.14.1` (an older macOS target resolves `2.11.0`); Intel macOS -> no
  solution (torch publishes no Intel macOS wheel). `uv sync --dry-run
  --extra torch` on Linux py3.12 and py3.14 -> `+ torch==2.14.1+cpu`.
  - **What a git-ignored `uv.lock` means for CI.** Every job that calls
    `uv run` (`lint`'s mypy, import-linter and asv steps,
    `import-without-extras`, `network`, `r-bridge`, `ml-extras`) locks
    afresh, so it now reads torch's versions from download.pytorch.org even
    when it installs no torch: `uv lock` takes 0.37 s instead of 0.10 s with
    a warm cache, 1.24 s instead of 1.20 s cold. CI therefore always gets
    the newest torch >= 2.9 (2.14.1 today), and an outage of
    download.pytorch.org fails those jobs at lock time. The hatch jobs
    (`test` matrix, `docs`, Read the Docs) do not read `[tool.uv]` for
    torch: `hatch env create hatch-test.py3.13-stable` took 5.7 s and
    installed no torch. The published wheel is unchanged apart from the
    extra: `uv build` writes `Provides-Extra: torch` and `Requires-Dist:
    torch>=2.9; extra == 'torch'`, and no index URL.
- **Resolved: how the default jobs skip the torch tests.** `addopts` becomes
  `-m "not network and not r and not torch"`, so every job but `ml-extras`
  deselects them; `ml-extras` passes `-m torch`, which overrides it, as the
  `network` and `r-bridge` jobs already do. No `pytest.importorskip`: it
  would show as skips in the default counts and turn a broken torch install
  in `ml-extras` into a green job. `tests/ml/test_torch.py` imports torch
  only inside a fixture, so collecting it works without torch and its 6
  tests that need no torch run in every job. `to_torch`'s doctest gets the
  marker from a root `conftest.py` hook (pytest 9.1.1 exports
  `pytest.DoctestItem`; the item is named `biotapy.ml._torch.to_torch`), so
  the example runs in `ml-extras` and is deselected elsewhere.
  `-W error::UserWarning`: `python -W error -c "import torch; import
  torch.utils.data"` prints nothing, and the 23 torch items pass under it.
- **Resolved: `to_torch`'s signature and items.** `to_torch(adata:
  AnnData, *, label_key: str | None = None, layer: str | None = None) ->
  "Dataset"` (R3.1; AnnData is the widest type, R3.2). Item `i` is a 1-D
  float32 tensor of sample `i`'s features, or `(features, label)` with
  `label_key`: a category, string or bool column -> int64 codes in
  `pd.Categorical(...).categories` order; a numeric column -> float32. The
  default `collate_fn` stacks them into `(batch, features)` float32 and
  `(batch,)` labels. A CSR table is referenced (`as_csr` does not copy a
  `csr_matrix`); another sparse format is converted to CSR once; a dense
  array (such as `layers["clr"]`) is referenced. `__getitem__` densifies one
  row with `table[i].toarray()[0]` and copies it with `np.array(...,
  dtype=float32)`, so an item never shares memory with the AnnData.
  Arguments are validated before torch is imported, so those errors are
  tested without the extra.
- **Facts the tasks rely on** (measured on the prototype; re-check each,
  R2.2):
  - `bt.datasets.toy()` is a 6 x 8 TreeData, `X` a CSR int64 matrix, `obs`
    one categorical column `group` (`A, A, A, B, B, B`); `bt.pp.clr` writes a
    dense `layers["clr"]`, `bt.pp.relative` a CSR `layers["relative"]`.
  - anndata 0.13.4 lists `X` as `layers[None]`, so `list(adata.layers)` is
    `[None]` on `toy()`; the missing-layer message therefore lists nothing.
  - `AnnData(obs=pd.DataFrame(index=["s1", "s2"]))` has `X is None`.
    Assigning a DataFrame to `layers` stores an ndarray, so no test can put a
    DataFrame there.
  - A view (`tdata[[0, 2, 4]]`) has a `SparseCSRMatrixView` `X`, a
    `csr_matrix` subclass, and works unchanged.
  - The first batch of 64 over a 2,000 x 50,000 CSR at 0.1% density peaks
    at 13.2 MB under `tracemalloc`; the dense table would be 800 MB.
  - `monkeypatch.setitem(sys.modules, "torch", None)` makes `import torch`
    raise `ImportError` whether or not torch is installed.

### Slice 4B global constraints (in addition to the Phase 4 list)
- The user approved the extra, the index and the CI job on 2026-10-09
  (decisions 10-12); 4.B0 does not ask again.
- Two environments, one `.venv`. `uv sync --all-groups` is an exact sync:
  it removes torch (0.6 s), giving CI's default-job environment. `uv run
  --group test --extra torch ...` adds it back from the uv cache (2.7 s
  warm). Run the default gates after `uv sync --all-groups`, the torch
  commands last.
- Gate before every commit, in this order (new files staged first):
  ```bash
  export BIOTAPY_DATA_DIR=<scratch pooch cache>
  uv sync --all-groups
  uvx prek run --all-files
  uv run --group test pytest -q -W error::UserWarning
  uv run --group test pytest -q -m "golden or network"
  rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html
  uv run --group test --extra torch pytest -q -m torch -W error::UserWarning
  ```
- Tests reach the function through `bt.ml.to_torch`; no test imports
  `biotapy.ml._torch`.
- No test file imports torch at module level (it would break collection in
  every job without the extra).

### Slice 4B review focus
1. **The whole table densified by accident** (a `.toarray()` on `X`).
   Expected: one row at a time. Test: 4.5
   `test_a_large_sparse_table_is_never_dense` (mutation-checked).
2. **A missing label turned into a code.** `pd.Categorical` codes a missing
   value as -1, which a loss function would take as a class index or reject
   far from the cause. Expected: `ValueError` naming the column. Test: 4.5
   `test_missing_label_raises`.
3. **Labels coded in an order the user did not expect.** Expected: the
   categorical's own order, sorted values for a string or bool column.
   Tests: 4.5 `test_labels_follow_the_category_order`,
   `test_string_and_bool_labels_become_sorted_codes`.
4. **A batch that writes back into the AnnData.** Expected: every item is a
   new tensor. Tests: 4.5 `test_keeps_the_input`,
   `test_an_item_does_not_share_memory_with_a_dense_float32_layer`.
5. **An integer class column read as a regression target.** Expected, and
   documented: a numeric column gives float32; class ids go in as a
   category. Test: 4.5 `test_numeric_labels_become_float32`. And the torch
   tests not running at all: `ml-extras` blocks merges (4.B1
   `test_ml_extras_job_blocks_merges`) and nothing uses `importorskip`.

---

### Task 4.B0: the extra `torch` and its CPU index (`build`)

**Files:** modify `pyproject.toml`, `tests/test_ci.py`,
`.knowledge/decisions/optional-heavy-dependencies.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `uv.lock` (git-ignored); the hatch environments (they
install no extra); `README.md` and `CHANGELOG.md` (4.D2 writes the release
notes and the install lines); the `import-without-extras` job (it must keep
installing no extra).
**Interfaces:**
- Consumes: nothing new.
- Produces: `pip install 'biotapy[torch]'`; `uv run --extra torch` installs
  torch from `https://download.pytorch.org/whl/cpu`.

- [x] **Step 1: Failing test.** In `tests/test_ci.py`, after
  `test_r_bridge_job_blocks_merges`:
  ```python
  def test_the_torch_extra_comes_from_the_cpu_index_and_nothing_else_does():
      pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
      assert pyproject["project"]["optional-dependencies"]["torch"] == ["torch>=2.9"]
      uv = pyproject["tool"]["uv"]
      assert uv["sources"] == {"torch": {"index": "pytorch-cpu"}}
      assert uv["index"] == [{"name": "pytorch-cpu", "url": "https://download.pytorch.org/whl/cpu", "explicit": True}]
  ```
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/test_ci.py -q` -> `1 failed, 16 passed`
  (`KeyError: 'torch'`).
- [x] **Step 3: Implement.**
  ```diff
  diff --git a/pyproject.toml b/pyproject.toml
  --- a/pyproject.toml
  +++ b/pyproject.toml
  @@ -46,6 +46,8 @@ dependencies = [
   ]
   # The ALDEx2 and MaAsLin 3 bridges (decisions/optional-heavy-dependencies); R and the R packages are the user's install.
   optional-dependencies.r = [ "rpy2>=3.6.8" ]
  +# ml.to_torch (decisions/optional-heavy-dependencies); 2.9.0 is torch's first release with CPython 3.14 wheels.
  +optional-dependencies.torch = [ "torch>=2.9" ]
   # https://docs.pypi.org/project_metadata/#project-urls
   urls.Documentation = "https://biotapy.readthedocs.io/"
   urls.Homepage = "https://github.com/pedrocr83/biotapy"
  @@ -124,6 +126,12 @@ envs.hatch-test.overrides.matrix.deps.env-vars = [
   ]
   envs.hatch-test.dependency-groups = [ "dev", "test" ]

  +[tool.uv]
  +# CI and contributors get torch's CPU wheels (196 MB on Linux), not PyPI's Linux wheel (555 MB plus CUDA packages).
  +# explicit = true: only torch, named in sources, comes from this index; the published metadata is unchanged.
  +sources.torch = { index = "pytorch-cpu" }
  +index = [ { name = "pytorch-cpu", url = "https://download.pytorch.org/whl/cpu", explicit = true } ]
  +
   [tool.ruff]
   line-length = 120
   src = [ "src" ]
  ```
- [x] **Step 4: Run, expect pass, and record the lock** -
  `uv run --group test pytest tests/test_ci.py -q` -> `17 passed`.
  `uv lock` -> `Resolved 195 packages` (189 at `cdc3b07`; the new entries
  are fsspec, mpmath, setuptools, sympy and torch twice, and no other
  version moves). `grep -A1 '^name = "torch"' uv.lock` -> `version =
  "2.14.1"` (macOS) and `version = "2.14.1+cpu"` (Linux, Windows); the torch
  version is the newest on the day, record it. `uv build --wheel`, then
  `unzip -p dist/*.whl '*/METADATA' | grep torch` -> `Provides-Extra:
  torch`, `Requires-Dist: torch>=2.9; extra == 'torch'`; delete `dist/`.
- [x] **Step 5: Knowledge.** In
  `.knowledge/decisions/optional-heavy-dependencies.md` (`generated` is
  `claude-code/<model>` at the commit time; `commit:` the parent's short
  sha):
  ```diff
  @@ -58,14 +58,26 @@ torch or an R installation into every install is unacceptable.[^spec]
     Ubuntu; without them rpy2 silently falls back to an ABI mode that cannot
     load R, so the `r-bridge` CI job sets `RPY2_CFFI_MODE=API` to fail the build
     instead, `.github/workflows/test.yaml`). biotapy imports it only through
  -  `import_optional` in `da/_r.py` and never bundles it.
  +  `import_optional` in `da/_r.py` and never bundles it. Phase 4 task 4.B0
  +  added the extra `torch` (`torch>=2.9`, approved 2026-10-09) for
  +  `ml.to_torch`: 2.9.0 is torch's first release with CPython 3.14 wheels.
  +  torch (BSD-3-Clause) brings filelock, fsspec, jinja2, networkx, setuptools,
  +  sympy (with mpmath) and typing-extensions; of these, fsspec, mpmath,
  +  setuptools and sympy are new to `uv.lock`, and no other version moves. uv
  +  installs it from PyTorch's CPU index (`[tool.uv]` in `pyproject.toml`:
  +  index `pytorch-cpu`, `https://download.pytorch.org/whl/cpu`, `explicit =
  +  true`, so no other package resolves there): `2.14.1+cpu` on Linux and
  +  Windows, `2.14.1` on macOS arm64 (no wheel exists for Intel macOS); the
  +  Linux wheel is 196 MB, against PyPI's 555 MB plus CUDA packages. The index
  +  is uv configuration only: the wheel's metadata says `torch>=2.9; extra ==
  +  'torch'`, and pip users get PyPI's torch.
   - Extras (names fixed now so docs never change), each added in the phase that first uses it:

     | Extra | Pulls | First used |
     |---|---|---|
     | `numba` | numba (>=0.67, supports up to Python 3.14) | perf track, only on benchmark evidence; also unlocks scikit-bio's `engine="numba"` |
     | `r` | rpy2 (3.6.8) | Phase 3 `da` bridges |
  -  | `torch` | torch | Phase 4 `ml` loaders |
  +  | `torch` | torch (>=2.9; under uv, CPU wheels from PyTorch's index) | Phase 4 `ml.to_torch` |
     | `plotnine` | plotnine | only if a `pl` function needs it |

  @@ -89,5 +101,8 @@ torch or an R installation into every install is unacceptable.[^spec]
   - CI imports every module with no extras installed (`import-without-extras`),
     so a lazy import leaking to module level fails fast. A job with all extras
     comes with the first extra.
  +- `uv.lock` is not committed, so every CI job that runs `uv run` resolves
  +  afresh, and reads torch's versions from download.pytorch.org even when it
  +  installs no torch (measured: `uv lock` 0.10 s -> 0.37 s with a warm cache).
  ```
  Use the versions and sizes Step 4 recorded if they differ. Tick 4.B0 in
  this concept. In `.knowledge/log.md`, below the title and above the
  newest section:
  ```markdown
  ## <date> (Phase 4, slice 4B)
  - **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md): the extra `torch` (`torch>=2.9`), what it brings, the CPU index uv installs it from and why, the extras table row, and that every CI `uv run` now reads torch's index; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.B0.
  ```
- [x] **Step 6: Gate and commit**
  ```bash
  git add pyproject.toml tests/test_ci.py .knowledge/decisions/optional-heavy-dependencies.md \
    .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  # the slice gate (Slice 4B global constraints) without its last command: the marker torch does not exist yet,
  # so `-m torch` selects nothing and pytest exits 5
  git commit -m "build: add the torch extra, installed by uv from PyTorch's CPU index

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1425 passed, 2 skipped, 54 deselected`; `36
  passed`; `build succeeded`.

### Task 4.5: `ml.to_torch`

**Files:** create `src/biotapy/ml/_torch.py`, `tests/ml/test_torch.py`;
modify `src/biotapy/ml/__init__.py`, `conftest.py`, `pyproject.toml`,
`docs/guide/machine_learning.md`, `docs/api.md`, `docs/conf.py`,
`docs/contributing.md`, `.knowledge/modules/ml.md`,
`.knowledge/modules/index.md`, `.knowledge/decisions/pure-by-default.md`,
`.knowledge/decisions/optional-heavy-dependencies.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `_core/_optional.py` (`import_optional` already names the
extra; `tests/core/test_optional.py` already uses `extra="torch"`);
`ml/_transformers.py`; `tests/conftest.py` (the doctest hook must also see
`src/biotapy`, which only the root `conftest.py` does);
`docs/extensions/coming_from_r.py` (`R equivalent: none` adds no row);
`.github/workflows/test.yaml` (4.B1).
**Interfaces:**
- Consumes: `_core.as_csr`, `_core.import_optional`; the extra `torch`
  (4.B0); fixtures `assert_unchanged`, `make_adata` (`tests/conftest.py`).
- Produces: `bt.ml.to_torch(adata: AnnData, *, label_key: str | None =
  None, layer: str | None = None) -> torch.utils.data.Dataset`; the pytest
  marker `torch`, deselected by default; `-m torch` runs the 22 tests and the
  doctest.

- [x] **Step 1: Failing tests.** Register the marker and deselect it by
  default:
  ```diff
  diff --git a/pyproject.toml b/pyproject.toml
  @@ -228,12 +228,13 @@
   [tool.pytest]
   addopts = [
  -  "--import-mode=importlib", # allow using test files with same name
  +  "--import-mode=importlib",             # allow using test files with same name
     "--doctest-modules",
     "-m",
  -  "not network and not r",
  +  "not network and not r and not torch",
   ]
   markers = [
     "golden: compares against R or HUMAnN golden files (contracts/r-golden-parity)",
     "network: downloads data; runs only in the dedicated CI job",
     "r: needs R and rpy2 (extra `r`)",
  +  "torch: needs PyTorch (extra `torch`); runs only in the ml-extras CI job",
   ]
  ```
  (The `--import-mode` line is pyproject-fmt's realignment, which prek
  applies when `addopts` changes.) Create `tests/ml/test_torch.py`:
  ```python
  import pickle
  import sys
  import tracemalloc

  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from anndata import AnnData
  from hypothesis import given
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays

  import biotapy as bt


  @pytest.fixture
  def torch():
      # Imported here, not at module level: without the extra, collecting this file must still work so the
      # tests below that need no torch run in every job.
      import torch

      return torch


  def _stacked(dataset, torch):
      return torch.stack([dataset[i] for i in range(len(dataset))]).numpy()


  @pytest.mark.torch
  def test_batches_hold_float32_rows_and_int64_labels(torch):
      loader = torch.utils.data.DataLoader(bt.ml.to_torch(bt.datasets.toy(), label_key="group"), batch_size=4)
      features, labels = next(iter(loader))
      assert features.shape == (4, 8) and features.dtype == torch.float32
      assert labels.shape == (4,) and labels.dtype == torch.int64
      assert [len(batch[0]) for batch in loader] == [4, 2]


  @pytest.mark.torch
  def test_is_a_torch_dataset_with_one_item_per_sample(torch):
      dataset = bt.ml.to_torch(bt.datasets.toy())
      assert isinstance(dataset, torch.utils.data.Dataset) and len(dataset) == 6


  @pytest.mark.torch
  def test_items_are_the_rows_of_x(torch):
      tdata = bt.datasets.toy()
      np.testing.assert_array_equal(_stacked(bt.ml.to_torch(tdata), torch), tdata.X.toarray().astype(np.float32))


  @pytest.mark.torch
  def test_labels_follow_the_category_order(torch):
      tdata = bt.datasets.toy()
      tdata.obs["group"] = tdata.obs["group"].cat.reorder_categories(["B", "A"])
      dataset = bt.ml.to_torch(tdata, label_key="group")
      assert [int(dataset[i][1]) for i in range(6)] == [1, 1, 1, 0, 0, 0]


  @pytest.mark.torch
  def test_a_subset_of_one_dataset_keeps_the_codes_that_per_split_datasets_shift(torch):
      tdata = bt.datasets.toy()
      test = [3, 5]  # both group B: the split lacks class A
      whole = bt.ml.to_torch(tdata, label_key="group")
      subset = torch.utils.data.Subset(whole, test)
      assert [int(subset[i][1]) for i in range(2)] == [1, 1]
      # anndata drops the unused category on subsetting, so a dataset built per split recodes B as 0.
      per_split = bt.ml.to_torch(tdata[test], label_key="group")
      assert [int(per_split[i][1]) for i in range(2)] == [0, 0]


  @pytest.mark.torch
  @pytest.mark.parametrize(
      ("values", "codes"), [(["b", "a", "c", "a", "b", "c"], [1, 0, 2, 0, 1, 2]), ([True, False] * 3, [1, 0] * 3)]
  )
  def test_string_and_bool_labels_become_sorted_codes(torch, values, codes):
      tdata = bt.datasets.toy()
      tdata.obs["label"] = values
      dataset = bt.ml.to_torch(tdata, label_key="label")
      assert [int(dataset[i][1]) for i in range(6)] == codes and dataset[0][1].dtype == torch.int64


  @pytest.mark.torch
  def test_numeric_labels_become_float32(torch):
      tdata = bt.datasets.toy()
      tdata.obs["age"] = [30, 41, 25, 60, 52, 47]
      dataset = bt.ml.to_torch(tdata, label_key="age")
      assert dataset[1][1].dtype == torch.float32 and float(dataset[1][1]) == 41.0


  @pytest.mark.torch
  @pytest.mark.parametrize(("make", "layer"), [(bt.pp.clr, "clr"), (bt.pp.relative, "relative")])
  def test_layer_is_read_instead_of_x(torch, make, layer):
      adata = make(bt.datasets.toy())
      expected = sp.csr_matrix(adata.layers[layer]).toarray().astype(np.float32)
      np.testing.assert_array_equal(_stacked(bt.ml.to_torch(adata, layer=layer), torch), expected)


  @pytest.mark.torch
  def test_a_view_of_some_samples_gives_those_samples(torch):
      tdata = bt.datasets.toy()
      dataset = bt.ml.to_torch(tdata[[0, 2, 4]], label_key="group")
      rows = torch.stack([dataset[i][0] for i in range(3)]).numpy()
      np.testing.assert_array_equal(rows, tdata.X[[0, 2, 4]].toarray())
      assert [int(dataset[i][1]) for i in range(3)] == [0, 0, 1]


  @pytest.mark.torch
  def test_keeps_the_input(torch, assert_unchanged):
      tdata = bt.datasets.toy()
      before = tdata.copy()
      dataset = bt.ml.to_torch(tdata, label_key="group")
      for features, _ in torch.utils.data.DataLoader(dataset, batch_size=4):
          features += 1
      dataset[0][0][:] = 99
      assert_unchanged(before, tdata)


  @pytest.mark.torch
  def test_an_item_does_not_share_memory_with_a_dense_float32_layer(torch):
      adata = bt.pp.clr(bt.datasets.toy())
      adata.layers["clr"] = adata.layers["clr"].astype(np.float32)
      before = adata.layers["clr"].copy()
      bt.ml.to_torch(adata, layer="clr")[0][:] = 99
      np.testing.assert_array_equal(adata.layers["clr"], before)


  @pytest.mark.torch
  def test_the_dataset_survives_pickling_and_spawned_workers(torch):
      # spawn and forkserver (macOS, Windows, Linux on Python 3.14) pickle the dataset to send it to each worker.
      dataset = bt.ml.to_torch(bt.datasets.toy(), label_key="group")
      clone = pickle.loads(pickle.dumps(dataset))
      assert len(clone) == 6 and torch.equal(clone[4][0], dataset[4][0]) and int(clone[4][1]) == int(dataset[4][1])
      workers = torch.utils.data.DataLoader(dataset, batch_size=2, num_workers=2, multiprocessing_context="spawn")
      alone = torch.utils.data.DataLoader(dataset, batch_size=2)
      for (rows, labels), (expected_rows, expected_labels) in zip(workers, alone, strict=True):
          assert torch.equal(rows, expected_rows) and torch.equal(labels, expected_labels)


  @pytest.mark.torch
  def test_editing_a_returned_label_does_not_change_the_next_fetch(torch):
      dataset = bt.ml.to_torch(bt.datasets.toy(), label_key="group")
      dataset[0][1].add_(3)
      assert int(dataset[0][1]) == 0


  @pytest.mark.torch
  @pytest.mark.parametrize("sparse", [True, False])
  def test_a_slice_index_raises_and_integer_indices_work(torch, sparse):
      tdata = bt.datasets.toy()
      adata = AnnData(X=tdata.X if sparse else tdata.X.toarray())
      dataset = bt.ml.to_torch(adata)
      with pytest.raises(TypeError, match="slice"):
          dataset[1:3]
      expected = tdata.X.toarray().astype(np.float32)
      np.testing.assert_array_equal(dataset[np.int64(2)].numpy(), expected[2])
      np.testing.assert_array_equal(dataset[-1].numpy(), expected[-1])
      with pytest.raises(IndexError):
          dataset[6]


  @pytest.mark.torch
  def test_a_bool_index_raises(torch):
      dataset = bt.ml.to_torch(bt.datasets.toy())
      with pytest.raises(TypeError, match="bool"):
          dataset[True]


  @pytest.mark.torch
  def test_references_x_so_a_later_change_shows(torch):
      tdata = bt.datasets.toy()
      dataset = bt.ml.to_torch(tdata)
      tdata.X.data[:] = 0
      assert float(dataset[0].sum()) == 0.0


  @pytest.mark.torch
  def test_a_large_sparse_table_is_never_dense(torch):
      # 2,000 x 50,000 at 0.1% density: 100,000 stored values; dense it would be 800 MB as float64.
      adata = AnnData(X=sp.random(2_000, 50_000, density=0.001, format="csr", random_state=0))
      tracemalloc.start()
      try:
          dataset = bt.ml.to_torch(adata)
          batch = next(iter(torch.utils.data.DataLoader(dataset, batch_size=64)))
          peak = tracemalloc.get_traced_memory()[1]
      finally:
          tracemalloc.stop()
      assert batch.shape == (64, 50_000)
      assert peak < 50_000_000


  @pytest.mark.torch
  def test_all_zero_sample_and_feature_are_zeros(torch, make_adata):
      adata = make_adata(np.array([[0, 0, 0], [3, 0, 1]]))
      rows = _stacked(bt.ml.to_torch(adata), torch)
      np.testing.assert_array_equal(rows, [[0, 0, 0], [3, 0, 1]])


  @pytest.mark.torch
  def test_single_sample(torch, make_adata):
      dataset = bt.ml.to_torch(make_adata(np.array([[2, 0, 5]])))
      batch = next(iter(torch.utils.data.DataLoader(dataset, batch_size=4)))
      assert len(dataset) == 1 and batch.tolist() == [[2.0, 0.0, 5.0]]


  @pytest.mark.torch
  def test_other_sparse_formats_are_read_as_csr(torch, make_adata):
      adata = make_adata(np.array([[1, 0], [0, 4]]))
      adata.X = sp.csc_matrix(adata.X)
      np.testing.assert_array_equal(_stacked(bt.ml.to_torch(adata), torch), [[1, 0], [0, 4]])


  @pytest.mark.torch
  @given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 1000)))
  def test_items_stack_back_to_the_table(dense):
      # Hypothesis rejects function-scoped fixtures, so this test imports torch itself.
      import torch

      rows = torch.stack(list(bt.ml.to_torch(AnnData(X=sp.csr_matrix(dense))))).numpy()
      np.testing.assert_array_equal(rows, dense.astype(np.float32))


  def test_missing_label_column_raises():
      with pytest.raises(KeyError, match="label_key='diet' is not a column of obs"):
          bt.ml.to_torch(bt.datasets.toy(), label_key="diet")


  def test_missing_label_raises():
      tdata = bt.datasets.toy()
      tdata.obs.loc["s2", "group"] = np.nan
      with pytest.raises(ValueError, match=r"label_key='group' has 1 missing value\(s\)"):
          bt.ml.to_torch(tdata, label_key="group")


  def test_missing_numeric_label_raises():
      tdata = bt.datasets.toy()
      tdata.obs["age"] = [30.0, np.nan, 25.0, 60.0, 52.0, 47.0]
      with pytest.raises(
          ValueError, match=r"label_key='age' has 1 missing value\(s\); drop those samples or fill them first"
      ):
          bt.ml.to_torch(tdata, label_key="age")


  def test_missing_layer_raises():
      with pytest.raises(KeyError, match=r"layer='clr' is not in adata.layers"):
          bt.ml.to_torch(bt.datasets.toy(), layer="clr")


  def test_a_table_that_is_not_an_array_raises():
      adata = AnnData(obs=pd.DataFrame(index=["s1", "s2"]))
      with pytest.raises(TypeError, match="adata.X is a NoneType; to_torch reads a NumPy array or a SciPy sparse matrix"):
          bt.ml.to_torch(adata)


  def test_options_are_keyword_only():
      with pytest.raises(TypeError):
          bt.ml.to_torch(bt.datasets.toy(), "group")


  def test_without_torch_names_the_extra(monkeypatch):
      # None in sys.modules makes `import torch` fail whether or not the extra is installed.
      monkeypatch.setitem(sys.modules, "torch", None)
      with pytest.raises(ImportError, match=r"pip install 'biotapy\[torch\]'"):
          bt.ml.to_torch(bt.datasets.toy())
  ```
  R11.2 coverage: happy path (batches, rows, labels, layers); edge cases
  (all-zero sample and feature, single sample; a NaN taxonomy rank does not
  apply, `to_torch` reads no taxonomy); purity (`test_keeps_the_input`, the
  shared-memory test); a Hypothesis property (items stack back to the
  table); no golden test (no R equivalent). The missing-torch test puts
  `None` in `sys.modules` instead of monkeypatching `import_optional`, so
  it reaches only the public API (R4.9) and holds with torch installed.
- [x] **Step 2: Run, expect failure** - `uv sync --all-groups && uv run
  --group test pytest tests/ml/test_torch.py -q` -> `6 failed, 22
  deselected`; `uv run --group test --extra torch pytest
  tests/ml/test_torch.py -q -m torch` -> `22 failed, 6 deselected`; every
  failure is `AttributeError: module 'biotapy.ml' has no attribute
  'to_torch'`.
- [x] **Step 3: Implement.** Create `src/biotapy/ml/_torch.py`:
  ```python
  """A samples x features table as a PyTorch dataset (extra ``torch``, decisions/optional-heavy-dependencies)."""

  import operator
  from typing import TYPE_CHECKING, Any, cast

  import numpy as np
  import numpy.typing as npt
  import pandas as pd
  import scipy.sparse as sp
  from anndata import AnnData

  from biotapy._core import as_csr, import_optional

  if TYPE_CHECKING:
      from torch import Tensor
      from torch.utils.data import Dataset

  # What a dataset item is read from: CSR referenced as is, or a dense array.
  Table = sp.csr_matrix | npt.NDArray[Any]


  def to_torch(adata: AnnData, *, label_key: str | None = None, layer: str | None = None) -> "Dataset":
      """A PyTorch dataset over the samples, one row of features per item.

      Parameters
      ----------
      adata
          Samples x features.
      label_key
          An ``obs`` column to pair with each row. A category, string or bool
          column becomes int64 codes in category order (sorted values for a
          string column), for a classifier; a numeric column becomes float32,
          for a regression. By default an item is the row alone.
      layer
          Read ``layers[layer]``, such as ``"clr"`` from :func:`biotapy.pp.clr`,
          instead of ``X``.

      Returns
      -------
      Dataset
          A map-style :class:`torch.utils.data.Dataset` of ``adata.n_obs`` items
          in ``obs`` order. Item ``i`` is sample ``i``'s features as a 1-D
          float32 tensor, or with ``label_key`` the pair ``(features, label)``.

      Raises
      ------
      ImportError
          torch is not installed: ``pip install 'biotapy[torch]'``.
      KeyError
          ``label_key`` is not a column of ``obs``, or ``layer`` is not a layer.
      ValueError
          The ``label_key`` column has a missing value.
      TypeError
          The table is neither a NumPy array nor a SciPy sparse matrix.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/machine_learning`

      A row is densified when its item is read, so a sparse table is never dense
      in full: an item costs 4 bytes x features, and a
      :class:`torch.utils.data.DataLoader` stacks items into batches, shuffling
      them if asked. The dataset holds the table it was given instead of
      copying it (a sparse table in another format than CSR is converted to CSR
      once), and reads the labels once, at construction. Do not modify ``adata``
      while using the dataset. Each item is a new tensor, so editing it leaves
      ``adata`` unchanged. Label codes follow
      ``pd.Categorical(adata.obs[label_key]).categories``. Label codes are per
      dataset: build one dataset and split it with
      :class:`torch.utils.data.Subset`, or splits that lack a class get
      different codes.

      Examples
      --------
      >>> import biotapy as bt
      >>> from torch.utils.data import DataLoader
      >>> dataset = bt.ml.to_torch(bt.datasets.toy(), label_key="group")
      >>> features, labels = next(iter(DataLoader(dataset, batch_size=4)))
      >>> features.shape, features.dtype, labels.tolist()
      (torch.Size([4, 8]), torch.float32, [0, 0, 0, 1])
      """
      table = _table(adata, layer)
      labels = None if label_key is None else _labels(adata, label_key)
      return _dataset(table, labels)


  def _dataset(table: Table, labels: npt.NDArray[np.int64] | npt.NDArray[np.float32] | None) -> "Dataset":
      """The dataset over a table and its labels; module-level so the dataset pickles for DataLoader workers."""
      # torch is the extra `torch`, so the Dataset subclass is defined only once it imports (rules.md R4.6, R3.6);
      # mypy treats torch as Any (pyproject.toml), so the subclassing needs the ignore with or without torch installed.
      torch: Any = import_optional("torch", extra="torch")
      targets = None if labels is None else torch.from_numpy(labels)

      class AnnDataDataset(torch.utils.data.Dataset):  # type: ignore[misc]
          def __len__(self) -> int:
              return int(table.shape[0])

          def __getitem__(self, index: int) -> "Tensor | tuple[Tensor, Tensor]":
              # One row at a time, so the full table is never dense (rules.md R6.2).
              if isinstance(index, bool):
                  msg = "a bool is not an index; use an integer"
                  raise TypeError(msg)
              position = operator.index(index)  # a slice would silently give one CSR row but several dense rows
              row = table[position].toarray()[0] if isinstance(table, sp.csr_matrix) else table[position]
              features = torch.from_numpy(np.array(row, dtype=np.float32))
              return features if targets is None else (features, targets[position].clone())

          def __reduce__(self) -> tuple[Any, tuple[Table, Any]]:
              # A class defined in a function cannot be pickled, which spawn and forkserver workers need.
              return (_dataset, (table, labels))

      return AnnDataDataset()


  def _table(adata: AnnData, layer: str | None) -> Table:
      """``X`` or ``layers[layer]``, referenced: CSR or dense as is, any other sparse format as CSR."""
      if layer is not None and layer not in adata.layers:
          msg = f"layer={layer!r} is not in adata.layers"
          raise KeyError(msg)
      values = adata.X if layer is None else adata.layers[layer]
      if isinstance(values, np.ndarray):
          return values
      if sp.issparse(values):
          return as_csr(values)
      name = "adata.X" if layer is None else f"layers[{layer!r}]"
      msg = f"{name} is a {type(values).__name__}; to_torch reads a NumPy array or a SciPy sparse matrix"
      raise TypeError(msg)


  def _labels(adata: AnnData, label_key: str) -> npt.NDArray[np.int64] | npt.NDArray[np.float32]:
      """``obs[label_key]`` as int64 codes in category order, or float32 for a numeric column."""
      if label_key not in adata.obs.columns:
          msg = f"label_key={label_key!r} is not a column of obs"
          raise KeyError(msg)
      # anndata types obs columns as Series | DataArray (its lazy variant); the data model guarantees a Series.
      values = cast("pd.Series", adata.obs[label_key])
      missing = int(values.isna().sum())
      if missing:
          msg = f"label_key={label_key!r} has {missing} missing value(s); drop those samples or fill them first"
          raise ValueError(msg)
      if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
          return values.to_numpy(dtype=np.float32, copy=True)
      return pd.Categorical(values).codes.astype(np.int64)
  ```
  `_table` and `_labels` are single-use helpers (R4.4): `to_torch` with
  their branches inline exceeds R5's 8 branches; no library call reads an
  AnnData slot or codes a label column this way (R2.1: `pd.Categorical` is
  the call `_labels` wraps). `src/biotapy/ml/__init__.py`:
  ```python
  from ._torch import to_torch
  from ._transformers import CLR, PrevalenceFilter

  __all__ = ["CLR", "PrevalenceFilter", "to_torch"]
  ```
  Append to the root `conftest.py`:
  ```python


  @pytest.hookimpl(tryfirst=True)
  def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
      """Give the doctests of a module that needs the extra `torch` its marker, before `-m` deselects."""
      for item in items:
          if isinstance(item, pytest.DoctestItem) and item.name.startswith("biotapy.ml._torch."):
              item.add_marker(pytest.mark.torch)
  ```
  and the mypy override:
  ```diff
  @@ -221,6 +221,10 @@ overrides = [
     { module = "sklearn.*", follow_untyped_imports = true, implicit_reexport = true },
     # threadpoolctl 3.7.0 is a single module, which cannot carry a py.typed marker; same treatment.
     { module = "threadpoolctl", follow_untyped_imports = true },
  +  # torch is the extra `torch`: absent where prek runs mypy, present where a contributor synced the extra. Not
  +  # following it makes its names Any in both, so mypy gives one answer, and skips analysing 700 MB of torch.
  +  { module = "torch", ignore_missing_imports = true, follow_imports = "skip" },
  +  { module = "torch.*", ignore_missing_imports = true, follow_imports = "skip" },
   ]
  ```
- [x] **Step 4: Run, expect pass** - `uv sync --all-groups && uv run --group
  test pytest tests/ml/test_torch.py src/biotapy/ml -q -W
  error::UserWarning` -> `8 passed, 23 deselected` (the 6 tests that need
  no torch and the two transformer doctests); `uv run --group dev --group
  doc mypy` -> `Success: no issues found in 65 source files`; `uv run
  --no-dev python -c "import importlib, pkgutil, sys, biotapy;
  [importlib.import_module(m.name) for m in
  pkgutil.walk_packages(biotapy.__path__, 'biotapy.')]; print('torch' in
  sys.modules)"` -> `False`. Then `uv run --group test --extra torch pytest
  tests/ml/test_torch.py src/biotapy/ml -q -W error::UserWarning -m torch`
  -> `23 passed, 8 deselected`; `tests/ml/test_torch.py -m torch` also
  under `--hypothesis-seed=1`, `2`, `3` (`22 passed, 6 deselected` each);
  `uv run --group dev --group doc --extra torch mypy` -> `Success: no
  issues found in 65 source files`.
- [x] **Step 5: Docs.**
  ````diff
  diff --git a/docs/guide/machine_learning.md b/docs/guide/machine_learning.md
  @@ -2,7 +2,8 @@

   `bt.ml` holds scikit-learn transformers, so microbiome preprocessing can sit
   inside a [Pipeline](https://scikit-learn.org/stable/modules/compose.html) and
  -be fitted on the training samples of each cross-validation fold only.
  +be fitted on the training samples of each cross-validation fold only, and
  +turns a table into a PyTorch dataset.

   ## Which steps leak

  @@ -53,3 +54,51 @@ scikit-learn already has it: `Normalizer(norm="l1")` divides each sample by
   its total, keeps a sparse matrix sparse and leaves an all-zero sample at zero,
   as `bt.pp.relative` does. Before `CLR` it changes nothing but the scale the
   pseudocount is on.
  +
  +## PyTorch
  +
  +`bt.ml.to_torch` turns a table into a PyTorch
  +[dataset](https://docs.pytorch.org/docs/stable/data.html) with one item per
  +sample, for a `DataLoader` to batch and shuffle. It needs the extra `torch`:
  +
  +```bash
  +pip install 'biotapy[torch]'
  +```
  +
  +```python
  +import biotapy as bt
  +from torch.utils.data import DataLoader
  +
  +tdata = bt.pp.clr(bt.datasets.toy())
  +dataset = bt.ml.to_torch(tdata, label_key="group", layer="clr")
  +for features, labels in DataLoader(dataset, batch_size=32, shuffle=True):
  +    ...  # features: float32, batch x features; labels: int64 codes
  +```
  +
  +An item is a sample's features as a float32 tensor, or with `label_key` the
  +pair `(features, label)`. A category, string or bool column gives int64 codes
  +in category order, for a classifier; a numeric column gives float32, for a
  +regression. A missing label raises: drop those samples first.
  +
  +Rows are densified one at a time, as they are read, so a sparse table is never
  +dense in full. The dataset holds the table without copying it and reads the
  +labels once, so do not modify the AnnData while you use the dataset.
  +
  +`to_torch` neither splits nor fits anything. Build one dataset over the whole
  +table and split it with `torch.utils.data.Subset`: label codes are per
  +dataset, and anndata drops a category a subset lacks, so a dataset built per
  +split recodes the classes when a split misses one. A step that learns from the
  +samples, like the prevalence filter, is fitted on the training samples only and
  +applied to both splits:
  +
  +```python
  +from torch.utils.data import Subset
  +
  +train, test = [0, 1, 3, 4], [2, 5]
  +keep = bt.ml.PrevalenceFilter(min_prevalence=0.1).fit(tdata[train].X).get_support()
  +dataset = bt.ml.to_torch(tdata[:, keep], label_key="group")
  +train_set, test_set = Subset(dataset, train), Subset(dataset, test)
  +```
  +
  +torch publishes no wheel for Intel macOS, so the extra does not install there.
  +
  +On Linux, pip installs PyPI's torch, which brings CUDA libraries. For a
  +CPU-only torch, install it from PyTorch's CPU index first:
  +
  +```bash
  +pip install torch --index-url https://download.pytorch.org/whl/cpu
  +pip install 'biotapy[torch]'
  +```
  diff --git a/docs/api.md b/docs/api.md
  @@ -120,6 +120,7 @@ Public functions are listed here as they ship, from Phase 1 onward.

       ml.CLR
       ml.PrevalenceFilter
  +    ml.to_torch
   ```

   ## Plots
  diff --git a/docs/conf.py b/docs/conf.py
  @@ -112,6 +112,7 @@ intersphinx_mapping = {
       "matplotlib": ("https://matplotlib.org/stable/", None),
       "mudata": ("https://mudata.scverse.org/stable/", None),
       "sklearn": ("https://scikit-learn.org/stable/", None),
  +    "torch": ("https://docs.pytorch.org/docs/stable/", None),
   }

   # List of patterns, relative to source directory, that match files and
  diff --git a/docs/contributing.md b/docs/contributing.md
  @@ -49,7 +49,7 @@ uv run --group test pytest
   ```

   Network and golden tests are excluded by default (`[tool.pytest]` in
  -`pyproject.toml` sets `-m "not network and not r"`). Run them explicitly:
  +`pyproject.toml` sets `-m "not network and not r and not torch"`). Run them explicitly:

   ```bash
   uv run --group test pytest -m "network or golden"
  @@ -80,6 +80,16 @@ BIOTAPY_DATA_DIR=.pooch uv run --group test --extra r pytest -m r
   CI runs them in the `r-bridge` job, with R 4.5.3 and the Bioconductor 3.22 packages the golden image
   pins.

  +### PyTorch tests
  +
  +Tests that need PyTorch (`bt.ml.to_torch`, and its docstring example) carry the marker `torch` and are
  +excluded from the runs above. The `torch` extra installs torch's CPU wheel from PyTorch's index
  +(`[tool.uv]` in `pyproject.toml`):
  +
  +```bash
  +uv run --group test --extra torch pytest -m torch
  +```
  +
   ### Regenerating the R golden files
  ````
  The two guide snippets were run with the extra (the split one with
  `train, test = [0, 1, 3, 4], [2, 5]` and `min_prevalence=0.9`); the API
  page shows `Return type: Dataset` linked to PyTorch's `Dataset`.
- [x] **Step 6: Knowledge.** (`generated` and `commit:` as in 4.B0.)
  ```diff
  diff --git a/.knowledge/modules/ml.md b/.knowledge/modules/ml.md
  @@ -1,21 +1,21 @@
   ---
   type: Module
   title: ml
  -description: scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold; they take arrays, sparse matrices and DataFrames, never AnnData.
  +description: scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold, taking arrays, sparse matrices and DataFrames; and to_torch, a PyTorch dataset over an AnnData's rows behind the extra torch.
   resource: /src/biotapy/ml/
   paths: ["src/biotapy/ml/**"]
  -tags: [ml, scikit-learn]
  +tags: [ml, scikit-learn, torch]
   status: stable
  -generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T01:55:36Z }
  -commit: 9883786
  +generated: { by: claude-code/<model>, at: <UTC commit time> }
  +commit: <parent short sha>
   ---

   # Responsibility

   Owns `bt.ml.*`, the top layer's machine-learning entry points. Today that is
  -two scikit-learn transformers (`_transformers.py`); `ml.to_torch` (slice 4B)
  -and `ml.embed` (slice 4C) are planned in
  -[phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md) and do not exist
  +two scikit-learn transformers (`_transformers.py`) and `to_torch`
  +(`_torch.py`, extra `torch`); `ml.embed` (slice 4C) is planned in
  +[phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md) and does not exist
   yet. Owns no reader and no table-level transform: `pp.filter_features` and
   `pp.clr` stay `pp`'s, the transformers are their fold-safe forms.
  @@ -26,6 +26,8 @@
   - `_transformers.py:CLR` - the centred log-ratio of each sample; stateless,
     so `fit` only validates.
  +- `_torch.py:to_torch` - a map-style `torch.utils.data.Dataset` over `X` or a
  +  layer, one float32 row per item, paired with an `obs` label when asked.
  @@ -33,6 +33,7 @@

  -- They take arrays, sparse matrices and DataFrames, not AnnData: inside a
  -  `Pipeline` the splitter hands over `X`, not an AnnData. This is why `ml`
  -  is the one place classes are allowed
  +- The transformers take arrays, sparse matrices and DataFrames, not AnnData:
  +  inside a `Pipeline` the splitter hands over `X`, not an AnnData (`to_torch`
  +  is the exception: it takes an AnnData, for its `obs` labels and layers).
  +  This is why `ml` is the one place classes are allowed
     ([function-shape](/contracts/function-shape.md), rules.md R3.6).
  @@ -52,6 +54,13 @@
   - CSR stays sparse: `PrevalenceFilter` reads `indices` and `data` of the CSR
     and its output stays sparse; `CLR` accepts sparse input and densifies once
     inside `pseudocounted` (rules.md R6.2).
  +- `to_torch` densifies one row when its item is read, never the table, and
  +  holds `X` (or the layer) instead of copying it, so the AnnData must not be
  +  modified while the dataset is used; every item is a new tensor (label
  +  included), so editing it leaves the AnnData unchanged. A non-integer index
  +  raises `TypeError`. Labels are converted once, at construction: category, string or bool -> int64 codes in category order, numeric
  +  -> float32; a missing label raises. Codes are per dataset (anndata drops a category a subset lacks), so build one dataset and split it with `torch.utils.data.Subset`; per-split datasets recode a class a split lacks. `_torch.py:to_torch`, `_torch.py:_labels`.
  +- `to_torch` validates its arguments before it imports torch, so its error
  +  tests run in every CI job, not only in `ml-extras`. `_torch.py:to_torch`.
  @@ -62,6 +71,8 @@
   - scikit-learn: `BaseEstimator`, `SelectorMixin`, `TransformerMixin`,
     `OneToOneFeatureMixin`, `validate_data`. scikit-bio: `clr`.
  +- torch (extra `torch`), only through `import_optional` inside `to_torch`
  +  ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).
  @@ -69,7 +80,9 @@
   `uv run --group test pytest tests/ml` (the scikit-learn estimator checks, the
  -`pp` parity tests, a `Pipeline` cross-validation test). The pseudocount
  +`pp` parity tests, a `Pipeline` cross-validation test, `to_torch`'s argument
  +errors). `uv run --group test --extra torch pytest -m torch` runs the
  +`to_torch` tests and its docstring example. The pseudocount
   warning's text is unit-tested in `tests/core/test_composition.py`.
  @@ -91,3 +104,17 @@
     class exists. Both transformers have no R equivalent (`R equivalent: none`)
     and so no golden test; parity is against `pp`.
  +- `to_torch`'s `Dataset` subclass is defined inside the module-level
  +  `_dataset(table, labels)`, after `import_optional`: a module-level class
  +  would import torch with biotapy. A class defined in a function cannot be
  +  pickled, so it defines `__reduce__` returning `(_dataset, (table, labels))`;
  +  that is what lets `DataLoader(num_workers>0)` work under spawn and
  +  forkserver (macOS, Windows, Linux on Python 3.14).
  +  mypy does not follow torch (`follow_imports = "skip"` in `pyproject.toml`),
  +  so it is `Any` with or without the extra installed and the class line
  +  carries `# type: ignore[misc]` (subclassing `Any`) in both. The return
  +  annotation is the bare `"Dataset"`: sphinx-autodoc-typehints renders a
  +  subscripted `Dataset[Tensor]` from a `TYPE_CHECKING` import as a broken
  +  cross-reference, which `nitpicky` fails. `_torch.py:_dataset`.
  +- The docstring example needs torch: the root `conftest.py` gives the
  +  doctests of `biotapy.ml._torch` the marker `torch`, so the default run
  +  deselects them and `-m torch` runs them. `conftest.py:pytest_collection_modifyitems`.
  +- Without the extra, `import biotapy` never imports torch. With torch
  +  installed, scikit-bio 0.7.4 imports it at module level
  +  (`skbio.util._testing`, reached through `_core._tree`), so `'torch' in
  +  sys.modules` after `import biotapy` tells nothing there; the
  +  `import-without-extras` check is meaningful only on a torch-free
  +  environment, as in CI.
  +- anndata 0.13 lists `X` as `layers[None]`, so `list(adata.layers)` holds
  +  `None` even when no layer was added; `to_torch`'s missing-layer error does
  +  not list the layers. `_torch.py:_table`.
  diff --git a/.knowledge/decisions/pure-by-default.md b/.knowledge/decisions/pure-by-default.md
  @@ -30,6 +30,7 @@
   | `ml` estimators (`PrevalenceFilter`, `CLR`) | scikit-learn's protocol: `fit` stores what it learns on the estimator and returns it; `transform` returns a new array | never the data |
  +| `ml.to_torch` | a new `torch.utils.data.Dataset` that holds `X` (or the layer) without copying it; each item a new tensor | never |
  diff --git a/.knowledge/decisions/optional-heavy-dependencies.md b/.knowledge/decisions/optional-heavy-dependencies.md
  @@ -89,6 +89,12 @@
   - Optional modules are imported inside the function through
     `biotapy._core.import_optional(name, extra)`, which raises `ImportError`
     naming the extra to install.
  +- A class that must inherit from an extra's base (rules.md R3.6), such as
  +  `ml.to_torch`'s torch `Dataset`, is defined inside the function after
  +  `import_optional` (`ml/_torch.py:_dataset`, built by `to_torch`). mypy does not follow the
  +  extra (`follow_imports = "skip"`, `ignore_missing_imports` in
  +  `pyproject.toml`), so the type check gives the same answer whether or not
  +  the extra is installed.
  ```
  Copy the new `ml` description into `.knowledge/modules/index.md`'s `ml`
  line (lowercase `scikit-learn`, as now). Tick 4.5 here; log line, first
  in the slice 4B section:
  ```markdown
  - **Update**: [ml](modules/ml.md) gains `to_torch` (entry point, per-row densify, references `X`, labels, validation before the torch import, the class-inside-the-function and mypy gotcha, the doctest marker, anndata's `layers[None]`), with the description copied into the [modules index](modules/index.md); [pure-by-default](decisions/pure-by-default.md) gains the `ml.to_torch` row; [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) says how a class inherits from an extra's base; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.5.
  ```
- [x] **Step 7: Gate and commit**
  ```bash
  git add src/biotapy/ml/_torch.py src/biotapy/ml/__init__.py tests/ml/test_torch.py conftest.py pyproject.toml \
    docs/guide/machine_learning.md docs/api.md docs/conf.py docs/contributing.md \
    .knowledge/modules/ml.md .knowledge/modules/index.md .knowledge/decisions/pure-by-default.md \
    .knowledge/decisions/optional-heavy-dependencies.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  # the slice gate (Slice 4B global constraints)
  git commit -m "feat(ml): add to_torch, a PyTorch dataset that densifies one row at a time

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1432 passed, 2 skipped, 77 deselected`; `36
  passed`; `build succeeded`; `23 passed, 1488 deselected`.

### Task 4.B1: CI job `ml-extras`

**Files:** modify `.github/workflows/test.yaml`, `tests/test_ci.py`,
`docs/contributing.md`, `.knowledge/decisions/optional-heavy-dependencies.md`,
`.knowledge/modules/ml.md`, `.knowledge/roadmap/phase-4-ml-multiomics.md`,
`.knowledge/log.md`.
**Not touched:** the `import-without-extras` job (it proves biotapy imports
without torch, so it installs no extra); the hatch `test` matrix and the
`docs` job (no torch, decision 17); a pooch cache in `ml-extras` (no torch
test downloads anything; 4.4b adds it with the MGM `network` test).
**Interfaces:**
- Consumes: the marker `torch` and `-m torch` (4.5); the extra (4.B0).
- Produces: the job `ml-extras`, required through `check.needs`.

- [x] **Step 1: Failing tests.** In `tests/test_ci.py`, before
  `test_the_torch_extra_comes_from_the_cpu_index_and_nothing_else_does`:
  ```python
  def test_ml_extras_job_runs_the_torch_marker_on_python_3_13():
      job = WORKFLOW["jobs"]["ml-extras"]
      setup = next(step for step in job["steps"] if step.get("uses", "").startswith("astral-sh/setup-uv@"))
      assert job["runs-on"] == "ubuntu-latest" and setup["with"]["python-version"] == "3.13"
      runs = [step.get("run", "").strip() for step in job["steps"]]
      assert "uv run --group test --extra torch pytest -m torch" in runs


  def test_the_sdist_ships_the_root_conftest_that_marks_the_torch_doctests():
      sdist = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["hatch"]["build"]["targets"][
          "sdist"
      ]
      assert "/conftest.py" in sdist["include"]


  def test_ml_extras_job_has_a_timeout():
      assert WORKFLOW["jobs"]["ml-extras"]["timeout-minutes"] == 30


  def test_ml_extras_job_blocks_merges():
      assert "ml-extras" in WORKFLOW["jobs"]["check"]["needs"]
  ```
  The Checkpoint B fix pass added the two tests before it: `/conftest.py`
  joins `build.targets.sdist.include` in `pyproject.toml` (its marker hook
  gives the `to_torch` doctest the marker `torch`, so an sdist test run needs
  it), and the job gets `timeout-minutes: 30`, as `r-bridge` has.
- [x] **Step 2: Run, expect failure** -
  `uv run --group test pytest tests/test_ci.py -q` -> `2 failed, 17 passed`
  (`KeyError: 'ml-extras'`; `AssertionError` on `check.needs`).
- [x] **Step 3: Implement.**
  ````diff
  diff --git a/.github/workflows/test.yaml b/.github/workflows/test.yaml
  @@ -212,6 +212,26 @@ jobs:
             RPY2_CFFI_MODE: API
           run: uv run --group test --extra r pytest -m r

  +  # Runs the tests that need PyTorch (marker torch, extra torch), which every other job deselects; uv installs torch's
  +  # CPU wheel from PyTorch's index ([tool.uv] in pyproject.toml, decisions/optional-heavy-dependencies).
  +  ml-extras:
  +    runs-on: ubuntu-latest
  +    timeout-minutes: 30
  +    steps:
  +      - uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1 # v7.0.1
  +        with:
  +          filter: blob:none
  +          fetch-depth: 0
  +          persist-credentials: false
  +      - name: Install uv
  +        uses: astral-sh/setup-uv@c18668ad3cf93ea998bef934396af7bb5c839dc7 # v10.2.0
  +        with:
  +          python-version: "3.13"
  +      # Syncs the environment, so the install shows apart from the tests in the job log.
  +      - name: Log the torch version
  +        run: uv run --group test --extra torch python -c "import torch; print('torch', torch.__version__)"
  +      - name: Run the PyTorch tests
  +        run: uv run --group test --extra torch pytest -m torch
  +
     # Builds the docs as Read the Docs does, executing every notebook (Phase 1 exit gate): the
     # phyloseq vignette downloads GlobalPatterns, enterotype and esophagus through the network
     # job's pooch cache.
  @@ -269,6 +289,7 @@ jobs:
         - import-without-extras
         - network
         - r-bridge
  +      - ml-extras
         - docs
       runs-on: ubuntu-latest
       steps:
  diff --git a/docs/contributing.md b/docs/contributing.md
  @@ -90,6 +90,8 @@ excluded from the runs above. The `torch` extra installs torch's CPU wheel from
   uv run --group test --extra torch pytest -m torch
   ```

  +CI runs them in the `ml-extras` job, on Linux with Python 3.13.
  +
   ### Regenerating the R golden files
  ````
- [x] **Step 4: Run, expect pass** -
  `uv run --group test pytest tests/test_ci.py -q` -> `19 passed`. Replay the
  job: `uv sync --all-groups` (no torch, as on a new runner), then its two
  commands -> `torch 2.14.1+cpu` (or the newest version on the day) and
  `23 passed, 1490 deselected`.
- [x] **Step 5: Knowledge.** (`generated` and `commit:` as in 4.B0.)
  ```diff
  diff --git a/.knowledge/decisions/optional-heavy-dependencies.md b/.knowledge/decisions/optional-heavy-dependencies.md
  @@ -101,8 +101,8 @@
   - CI imports every module with no extras installed (`import-without-extras`),
  -  so a lazy import leaking to module level fails fast. A job with all extras
  -  comes with the first extra.
  +  so a lazy import leaking to module level fails fast. Each extra has a job
  +  that installs it and runs its marker, deselected everywhere else: `r-bridge`
  +  (`r`) and `ml-extras` (`torch`, Linux, Python 3.13), both in `check.needs`.
  diff --git a/.knowledge/modules/ml.md b/.knowledge/modules/ml.md
  @@ -80,9 +80,9 @@
  -`to_torch` tests and its docstring example. The pseudocount
  +`to_torch` tests and its docstring example, as CI's `ml-extras` job does. The pseudocount
   warning's text is unit-tested in `tests/core/test_composition.py`.
  ```
  Tick 4.B1 here; log line, first in the slice 4B section:
  ```markdown
  - **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md): each extra has its CI job (`r-bridge`, `ml-extras`); [ml](modules/ml.md)'s Verification names the `ml-extras` job; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.B1.
  ```
- [x] **Step 6: Gate and commit**
  ```bash
  git add .github/workflows/test.yaml tests/test_ci.py docs/contributing.md \
    .knowledge/decisions/optional-heavy-dependencies.md .knowledge/modules/ml.md \
    .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  # the slice gate (Slice 4B global constraints); zizmor checks the new job
  git commit -m "ci: run the PyTorch tests in an ml-extras job

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1434 passed, 2 skipped, 77 deselected`; `36
  passed`; `build succeeded`; `23 passed, 1490 deselected`.

### Checkpoint B - review slice 4B
- [x] Review the whole slice (superpowers:requesting-code-review) against
  every contract, pure-by-default, optional-heavy-dependencies, the Phase 4
  and slice 4B review focus; then a fix pass, one commit per finding, each
  with a test. Record the counts and the fix range here.
  Review: 0 Critical, 1 Important, 5 Minor (+5 counted). Fix pass
  `d33d661..7a9c07a`, 7 commits. Scoped re-review: all 6 findings addressed,
  no new findings.
- [x] Run the slice's checks at the last commit and record them: the slice
  gate's five counts; `uv run --group test --extra torch coverage run -m
  pytest -m "torch or not torch" tests/ml/test_torch.py
  src/biotapy/ml/_torch.py` then `coverage report --include
  "src/biotapy/ml/_torch.py"` (exits 0; `ml/*` would add `_transformers.py`,
  which this run does not exercise, and `fail_under = 90` would fail the
  total). `_torch.py` was 100% on the prototype and is 100% (54 statements)
  after the fix pass.
  At `7a9c07a`: `pytest -q -W error::UserWarning` 1437 passed, 2 skipped, 79
  deselected; `-m "golden or network"` 36 passed; `--extra torch -m torch` 25
  passed, 1493 deselected; coverage `_torch.py` 54 statements, 100%.
- [x] Knowledge: [ml](/modules/ml.md) and
  [optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)
  already carry 4B (4.5, 4.B1); re-check them against the fix pass and
  bump only what changed, with log lines.
  Done at `7a9c07a`: [ml](/modules/ml.md) gains the bool index and the sdist
  `conftest.py`; [optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)
  the `ml-extras` timeout; the nine stale concepts re-checked (see the log).
- [x] Push the branch and open the PR only after the user approves that push
  (R13.3). The PR body carries the R9.2 reason: "Extra `torch` (`torch>=2.9`,
  BSD-3-Clause, approved 2026-10-09): `ml.to_torch` subclasses
  `torch.utils.data.Dataset`; an extra, never core (R9.3). uv installs the CPU
  wheel from `https://download.pytorch.org/whl/cpu` (`explicit = true`, torch
  only); the published metadata is unchanged." CI green, including
  `ml-extras`; record its runtime and the torch version it installed here
  (the prototype's local replay: sync 19.4 s cold, tests 14.8 s, 1.6 GB).
  - PR #29, Test run 37919311278 green (22 jobs); `ml-extras` ran 46 s with torch 2.14.1+cpu, 25 passed; merged as `fd6a8be`.
- [x] Ask the user to review slice 4B before slice 4C is expanded.
  - Approved 2026-10-09.

---

## Slice 4C - Embeddings

**Goal:** with `pip install 'biotapy[mgm]'`, `bt.ml.embed(tdata, "mgm")`
returns MGM's embedding of every sample (or writes `obsm["X_mgm"]`), equal to
MGM's own code to 2e-6, and any package adds a model by declaring an entry
point in the group `biotapy.embeddings`; CI's `ml-extras` job proves it end to
end on GlobalPatterns (exit gate 2).

### Slice 4C design
- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `_core/_download.py` | `make_pooch`, the one download cache, moved out of `datasets/_remote.py` (4.C0) |
  | `pyproject.toml` | the extra `mgm` (4.C1); the marker `mgm` and `addopts` (4.4); the entry point `mgm`, `pooch.processors.Unzip` in `untyped_calls_exclude`, `/tests/mgm` out of the sdist (4.4b) |
  | `ml/_embed.py` | `embed` and the group name `GROUP` (4.4) |
  | `conftest.py` (root) | `_EXTRA_DOCTESTS`: `biotapy.ml._embed`'s doctests get the marker `mgm` (4.4) |
  | `ml/_mgm.py` | the MGM plugin: `embed` (checks, imports, download, model, forward), `_sentences` (MGM's tokens), `_listed` (4.4b) |
  | `tests/ml/test_embed.py` | 19 tests with fake plugins, no extra needed (4.4) |
  | `tests/ml/test_mgm.py` | 4 tests that need no extra, 8 marked `mgm` (4.4b) |
  | `tests/data/mgm/` | `counts.csv` (synthetic genus table), `embeddings.csv` (MGM's own embeddings of it), `NOTICE.txt` (4.4b) |
  | `tests/mgm/export_reference.py` | writes `tests/data/mgm/` with MGM 0.5.8's own code, in a Python 3.11 environment of its own, never in CI (4.4b) |
  | `.github/workflows/test.yaml` | `ml-extras` gains a pooch cache and the `-m mgm` step (4.4b) |
  | `docs/guide/machine_learning.md` | "Embeddings" (4.4), "MGM" and "Filtering and leakage" (4.4b) |
  | `.knowledge/decisions/embedding-plugins.md` | the plugin decision (`draft` until the user confirms it) (4.4) |

  `uv.lock` stays git-ignored. No mypy override for transformers is needed:
  nothing imports it statically (`_mgm.py` reaches torch and transformers
  through `import_optional`, typed `Any`).
- **How slice 4C was checked.** Every file below was written into a scratch
  clone of the repository at `fd6a8be` (master after slice 4B) on branch
  `phase-4c` and committed one task at a time, after a stand-in
  `docs(roadmap)` commit that only adds the 4.C0 and 4.C1 checklist lines and
  the new 4.4 and 4.4b wording (`c55853d`). The commits are `82c6e6d` (4.C0),
  `5b73a1b` (4.C1), `0785751` (4.4), `a702d59` (4.4b). Each was gated on its
  committed tree (`git status --short` empty) with `.venv` synced without
  extras (`uv sync --all-groups`), as CI's default jobs are; `-m torch` ran
  with `--extra torch`, `-m mgm` with `--extra mgm` (torch 2.14.1+cpu,
  transformers 5.19.0, Python 3.13.2). Every run exported `BIOTAPY_DATA_DIR`
  to a scratch pooch cache, and `UV_CACHE_DIR`, `XDG_CACHE_HOME`, `HF_HOME`
  and `MPLCONFIGDIR` to scratch directories; `~/.cache/biotapy` never
  appeared, and no gate wrote under `~/.cache`.

  | Task state | `uvx prek run --all-files` | `uv run --group test pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` | `--extra torch pytest -q -m torch -W error::UserWarning` | `--extra mgm pytest -q -m mgm -W error::UserWarning` |
  |---|---|---|---|---|---|---|
  | base `fd6a8be` | passed (14 hooks) | 1437 passed, 2 skipped, 79 deselected | 36 passed, 1482 deselected | build succeeded | 25 passed, 1493 deselected | - |
  | stand-in `c55853d` | passed | 1437 passed, 2 skipped, 79 deselected | 36 passed, 1482 deselected | build succeeded | 25 passed, 1493 deselected | - |
  | 4.C0 `refactor(core)` (`82c6e6d`) | passed | 1439 passed, 2 skipped, 79 deselected | 36 passed, 1484 deselected | build succeeded | 25 passed, 1495 deselected | - (no marker yet) |
  | 4.C1 `build` (`5b73a1b`) | passed | 1440 passed, 2 skipped, 79 deselected | 36 passed, 1485 deselected | build succeeded | 25 passed, 1496 deselected | - (no marker yet) |
  | 4.4 `feat(ml)` (`0785751`) | passed | 1462 passed, 2 skipped, 80 deselected | 36 passed, 1508 deselected | build succeeded | 25 passed, 1519 deselected | 1 failed, 1543 deselected (see below) |
  | 4.4b `feat(ml)` (`a702d59`) | passed | 1470 passed, 2 skipped, 88 deselected | 36 passed, 1524 deselected | build succeeded | 25 passed, 1535 deselected | 9 passed, 1551 deselected |

  - Every default run also ends in "2 warnings": scikit-bio's
    `RuntimeWarning: invalid value encountered in divide` from 4.F1's
    `ancombc2` tests, as on master.
  - At 4.4, `-m mgm` selects only `ml.embed`'s docstring example, which runs
    MGM; MGM is registered in 4.4b, so it fails with `KeyError: "model='mgm'
    is not an installed embedding plugin; installed: []"`. No CI job runs
    `-m mgm` before 4.4b adds the step (decision 34).
  - The new default-run items: 4.C0 +2 (`tests/core/test_download.py`); 4.C1
    +1 (`tests/test_ci.py`); 4.4 +22 (19 in `tests/ml/test_embed.py`,
    `test_docstring_has_the_contract_sections[bt.ml.embed]`, and the
    knowledge-bundle checks of the new concept, 2); 4.4b +8 (4 in
    `tests/ml/test_mgm.py`, 1 in `tests/test_ci.py`, 3 size checks of
    `tests/data/mgm/*`). The 8 new deselections at 4.4b are the `mgm` tests.
  - Coverage with the extra (`uv run --group test --extra mgm coverage run -m
    pytest -m "mgm or not mgm" tests/ml/test_embed.py tests/ml/test_mgm.py
    tests/core/test_download.py src/biotapy/ml/_embed.py`, 34 passed):
    `_core/_download.py` 4 statements, `ml/_embed.py` 32, `ml/_mgm.py` 65,
    each 100%.
  - Hypothesis seeds 1, 2 and 3 on `tests/ml/test_embed.py`: `19 passed`
    each.
  - `uv run --group dev --group doc --extra mgm mypy` -> `Success: no issues
    found in 68 source files`, as without the extra. Without extras, the
    `import-without-extras` command leaves `'torch' in sys.modules` and
    `'transformers' in sys.modules` both `False`.
  - The `ml-extras` job replayed (Python 3.13, warm uv cache, cold pooch
    cache): the torch steps 11.8 s, then the transformers log and `-m mgm`
    21.2 s including the downloads (`9 passed`); the pooch cache holds 67 MB
    (the 33 MB wheel, its 36 MB of extracted files, GlobalPatterns).
  - Parity mutations (each applied to `_mgm.py` alone, then reverted): the
    relative-abundance denominator over all features, an ascending sort,
    last-token pooling, dropout left on, keeping `z > 0`, vocabulary order
    instead of the sort. Each fails `test_matches_mgm_s_own_embedding`; all
    but the first also fail `test_embeds_global_patterns_end_to_end`.
- **Resolved: MGM's sample embedding** (design note 8's [UNVERIFIED]).
  - MGM's paper, Methods 4.5 ("Sample Representation"): "we opted
    element-wise mean pooling to get the sample representation from our
    pretrained model"; for the fine-tuned model "the last token ('eos')".
    It names no layer and does not say whether `<bos>`/`<eos>` are in the
    mean.[^mgm]
  - MGM's only code that embeds samples is its notebook
    `MGM_Interpretability_Guideline.ipynb` (GitHub HUST-NingKang-Lab/MGM,
    commit `bee1469`, step 2 "Extract Sample Embeddings"), on a fine-tuned
    `GPT2ForSequenceClassification`: `model(**input,
    output_hidden_states=True).hidden_states[-1]`, then the row of the last
    token whose attention mask is 1, one sample at a time. The CLI's `predict`
    pools the same way through transformers' sequence classifier.
  - biotapy embeds with the pretrained general model, so it follows the
    paper: the mean of the last hidden layer (`last_hidden_state`, which is
    `hidden_states[-1]`, after GPT-2's final layer norm) over the sample's
    tokens, `<bos>`, genera and `<eos>`. `<bos> <eos>` keeps the mean defined
    for a sample with no known genus (decision 28).
- **Resolved: parity with MGM's own forward pass** (design note 8's
  [UNVERIFIED] Python 3.11 install). `microformer-mgm==0.5.8` with its pins
  installs on CPython 3.11.16 (numpy 1.24.3, pandas 2.0.3, torch 2.0.1+cpu
  from PyTorch's CPU index, transformers 4.33.3, scikit-learn 1.3.1,
  pytorch-lightning 2.0.6) through `uv run --no-project --python 3.11
  --with microformer-mgm==0.5.8 --with torch==2.0.1+cpu --extra-index-url
  https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match`,
  in 8 s with a warm cache. `tests/mgm/export_reference.py` builds a 7-sample
  genus table, runs MGM's `MicroCorpus` and `GPT2LMHeadModel` on it and writes
  `tests/data/mgm/embeddings.csv` (23,682 bytes) beside `counts.csv` (20,463
  bytes); two runs gave byte-identical files. biotapy's embedding of the same
  table differs from it by at most 1.67e-6 (torch 2.14.1+cpu, transformers
  5.19.0, 8 threads; 9.8e-7 with 4; 1.43e-6 on torch 2.9.0+cpu with
  transformers 5.0.0, the floors; 1.67e-6 on CPython 3.14.8).
- **Resolved: a torch first-call race.** On the prototype's CPU (16 logical
  cores, torch's default 8 threads), the first embedding in a fresh process
  differed from MGM's by up to 1.6e-4 in 5 of 8 processes (6 of 8 with eager
  attention), and later calls in the same process never (1.67e-6). Forward
  hooks put the first difference in layer 0's GELU (`h.0.mlp.act`,
  1.7e-4), and `torch.tanh` alone on a seeded 7 x 512 x 1024 tensor
  reproduces it: in 14 of 37 fresh processes the first call is off by up to
  9.08e-5 near `x = -5` (as if it saturated to -1; `tanh(-5)` is
  -0.9999092), and agrees with NumPy to 6e-8 from the second call on. With
  `OMP_NUM_THREADS` 2 or 4 it never happened (24 processes), with 8 in 3 of
  12; torch 2.13.0+cpu has it too (1 of 12), torch 2.0.1+cpu did not show it
  (0 of 8). Calling `torch.tanh(torch.zeros(1))` first removed it (16 of
  16). CI's runners have 4 vCPUs. (Amended 2026-10-10: the race is oneMKL's
  VML CPU-type cache, pytorch/pytorch#188792, and it also hits 4 threads,
  rarely; the claim above that 2 or 4 threads never showed it was a sample
  too small. biotapy now makes that serial call itself and tests at 1e-5;
  decision 33.)
- **Resolved: how MGM's vocabulary meets a biotapy table.**
  - MGM reads a genus from a column name with
    `str.extract(r'(g__[A-Za-z0-9_]+)')`, sums columns that give the same
    token, reindexes to `phylogeny.csv`'s 9,665 genera (dropping the rest),
    drops all-zero samples, divides by each sample's total over its known
    genera, standardises with the per-genus `mean` and `std`, keeps genera
    above the standardised zero, sorts them with pandas' `sort_values
    (ascending=False)`, and builds `<bos>` genera `<eos>` cut to 512
    (`mgm/src/MicroCorpus.py`). The vocabulary is `<pad>, <mask>, <bos>,
    <eos>` and then `phylogeny.csv`'s order (ids 0-3, genus `i` is `i + 4`);
    every `std` is positive (smallest 1.95e-7).
  - The plugin applies the same regex to `"g__" + var["genus"]`, so
    `Escherichia-Shigella` is `Escherichia` as in MGM. 3 of the 9,665 tokens
    can never be produced by MGM's regex (`g__Pseudo-nitzschia`,
    `g__Acrocirridae_gen._2_KJO-2009`, `g__Oxystominidae_gen._'Cricohalalaimus'`).
  - GlobalPatterns after `pp.tax_glom(..., "genus")`: 996 features (983
    genus names; `tax_glom` groups by lineage, so a name can repeat), 96 of
    them with a genus MGM does not know (`4-29`, `4041AA30`, `A17`,
    `Aquamonas`, `Arctic95A-2`, ...); 202 to 490 tokens per sample, none cut.
    enterotype: 551 genus names, 501 known. Neither dataset has a genus name
    with a space, so replacing spaces would change nothing (decision 30).
  - Unknown genera: one `UserWarning` per call, counting the features without
    a genus and those whose genus MGM does not know, naming up to five; a
    sample left with no known genus is embedded from `<bos> <eos>` (MGM would
    drop it) with a second warning naming up to five such samples.
- **Resolved: transformers' wheels and the lock** (design note 11's
  [UNVERIFIED]). transformers 5.19.0 (2026-10-06, Apache-2.0, `py3-none-any`,
  `Requires-Python >=3.10`) needs huggingface-hub (>=1.31,<3), numpy,
  packaging, pyyaml, regex, tokenizers (>=0.23.1,<0.24), typer,
  safetensors (>=0.8), tqdm. The binary ones ship `cp310-abi3` (tokenizers
  0.23.3, safetensors 0.8.0), `cp38-abi3` (hf-xet 1.7.0) or per-version
  wheels up to cp314 (regex 2026.9.29), plus `cp314t` for tokenizers, hf-xet
  and regex. `uv pip compile pyproject.toml --extra mgm --python-platform`
  for `x86_64-manylinux_2_28` and `aarch64-manylinux_2_28` on Python 3.12,
  3.13 and 3.14, `x86_64-pc-windows-msvc` and `aarch64-apple-darwin`
  (`MACOSX_DEPLOYMENT_TARGET=14.0`) on 3.12 and 3.14: `transformers==5.19.0`,
  `tokenizers==0.23.3`, torch `2.14.1+cpu` (Linux, Windows) or `2.14.1`
  (macOS); `x86_64-apple-darwin` has no solution (no torch wheel), as for the
  extra `torch`. The plain `x86_64-unknown-linux-gnu` target means
  manylinux_2_17, for which torch publishes nothing; that is the `torch`
  extra's behaviour already. `uv lock`: 211 packages against 195, the 16 new
  ones annotated-doc, anyio, h11, hf-xet, httpcore2, httpx2, httpx2-jsfetch,
  huggingface-hub 2.2.0, regex, safetensors, shellingham, tokenizers, tqdm,
  transformers, truststore, typer; no version moves. Their licences:
  Apache-2.0 (transformers, huggingface-hub, tokenizers, safetensors,
  hf-xet), Apache-2.0 AND CNRI-Python (regex), MPL-2.0 AND MIT (tqdm), MIT
  (typer, anyio, h11, truststore, annotated-doc), BSD-3-Clause (httpx2,
  httpcore2, httpx2-jsfetch), ISC (shellingham). Installed: about 97 MB
  beside torch's 742 MB (transformers 61 MB). transformers 5.0.0 was
  released 2026-01-26. On CPython 3.14.8, `pip install -e ".[mgm]"` (uv
  0.12.24, scratch environment) imports torch 2.14.1+cpu, transformers
  5.19.0 and tokenizers 0.23.3 and passes the parity check.
  `python -W error -c "import torch, transformers; from transformers import
  GPT2Config, GPT2Model"` prints nothing; building the model wrote nothing
  under a fresh `HOME` (transformers never contacts the Hugging Face Hub
  here: config and weights are local files).
- **Resolved: `batch_size`** (the outline's open question). MGM on the CPU,
  512 synthetic samples of 512 tokens, 8 threads: batch 1 took 11.5 s with
  no memory beyond the model's; 4, 8, 16, 32, 64 and 256 took 16.6-22.3 s
  and +89 MB, +214 MB, +428 MB, +778 MB, +1.2 GB and +4.8 GB of peak RSS.
  520 GlobalPatterns profiles (26 x 20): batch 1 8.8 s (8 threads) and 11.9 s
  (4 threads), batch 8 17.0 s / 25.1 s and +196 MB / +172 MB, batch 64
  18.5 s / 26.6 s and +1.6 GB. Dropping the attention mask changed nothing.
  The plugin runs one unpadded sample at a time, and `ml.embed` has no
  `batch_size` (decision 27).
- **Resolved: the plugin interface** (decision 15, as amended by decision
  27).
  - Group `biotapy.embeddings`; the entry point's name is the model's name;
    its value loads a callable `embed(adata: AnnData) -> numpy.ndarray`.
  - `ml.embed(adata: AnnData, model: str, *, inplace: bool = False) ->
    npt.NDArray[np.floating] | None` reads
    `importlib.metadata.entry_points(group="biotapy.embeddings")` on every
    call. Unknown name -> `KeyError: "model='nope' is not an installed
    embedding plugin; installed: [...]"`; a name two distributions register
    -> `ValueError` naming both (an `EntryPoints` lookup by name would pick
    one silently).
  - The result must be a `numpy.ndarray` (else `TypeError`), 2-D with
    `adata.n_obs` rows and at least one column, a float dtype, all finite
    (else `ValueError`); every message starts `plugin '<model>' returned`.
    It is returned unchanged, or with `inplace=True` written to
    `obsm[f"X_{model}"]` (data-model-slots; R3.3's `tl` convention) and
    `None` returned; a refused result is never written.
  - Tests install fake plugins as a real distribution: `<name>-1.0.dist-info`
    with `METADATA` and `entry_points.txt` under `tmp_path`,
    `monkeypatch.syspath_prepend`, and the callables on a module object put
    in `sys.modules`. `importlib.metadata` finds them through the same path
    as an installed package, so no test imports `ml/_embed.py` (R4.9).
- **Resolved: the weights through pooch.** URL
  `https://files.pythonhosted.org/packages/4b/9c/829a1e59d5e618756ce8e57553e15d13bb16c58d37896f03233d8697515a/microformer_mgm-0.5.8-py3-none-any.whl`,
  33,453,947 bytes, uploaded 2025-02-17, SHA-256
  `210891685565022ea869e88a7452769dc6fbe3f699d2a133e15338d1e30eb92b`
  (PyPI's published digest, and the download's). Extracted members (pooch
  1.9.0 `Unzip(members=...)`): `config.json` (760 bytes),
  `pytorch_model.bin` (35,731,997 bytes, SHA-256 `7a83b442...511b8`),
  `phylogeny.csv` (410,742 bytes, SHA-256 `016b8264...78db`) into
  `<cache>/microformer_mgm-0.5.8-py3-none-any.whl.unzip/mgm/resources/...`.
  `Unzip` returns paths in `os.walk` order, so `embed` keys them by file
  name. `torch.load(..., map_location="cpu", weights_only=True)` reads 101
  tensors; the 100 under `transformer.` load into `GPT2Model` with
  `strict=True` (`lm_head.weight` is left out). The cache is
  `_core.make_pooch` (4.C0): `BIOTAPY_DATA_DIR` or
  `pooch.os_cache("biotapy")`, the datasets' cache; R4.3 moves it to `_core`
  now that a second subpackage needs it (decision 32).
- **Resolved: markers and CI** (decision 31). The MGM tests that need the
  extra or the download carry the marker `mgm`, deselected by default
  (`addopts`: `-m "not network and not r and not torch and not mgm"`), and
  run in `ml-extras` after the torch tests with `BIOTAPY_DATA_DIR` on a pooch
  cache keyed on `src/biotapy/datasets/_remote.py` and `src/biotapy/ml/_mgm.py`
  (the end-to-end test also reads GlobalPatterns). They are not marked
  `network`: the `network` job runs `-m "network or golden"` without torch.
  The 4 MGM tests that need neither (registration, the two input errors, the
  missing extra) run in every job.
- **Facts the tasks rely on** (measured on the prototype; re-check each,
  R2.2):
  - `bt.datasets.toy()`: `f4` and `f5` are both `Bacteroides`, `f8` has no
    genus; `bt.pp.tax_glom(toy(), "genus")` (`dropna=True` by default) keeps
    6 features, `Blautia, Roseburia, Faecalibacterium, Bacteroides,
    Prevotella, Escherichia`, all in MGM's vocabulary, and embeds without a
    warning.
  - `pd.read_csv("tests/data/mgm/counts.csv", index_col=0)["genus"]` is
    pandas 3's `str` dtype with one NaN; `.astype("string")` turns it into
    `<NA>`, which `str.extract` keeps and `Index.get_indexer` maps to -1.
  - `_core.sum_by(X, codes, n)` drops negative codes and sums the rest;
    `_core.divide_rows` leaves a zero-total row at zero.
  - `importlib.metadata.EntryPoint.dist.name` is the `Name:` of the
    distribution's `METADATA`.
  - pytest's `--doctest-modules` imports every `.py` under `tests/`, so
    `tests/mgm/export_reference.py` imports MGM, torch and transformers
    inside its functions only.
  - `transformers.GPT2Config.from_json_file` on MGM's `config.json` emits no
    warning under `-W error`; `GPT2Model` defaults to the `sdpa` attention.
  - pyproject-fmt realigns the `--import-mode` comment in `addopts` when the
    marker expression grows, and moves `entry-points` after `urls`.

### Slice 4C global constraints (in addition to the Phase 4 list)
- The user approved the extra `mgm`, the in-tree plugin, the weights
  download and the interface on 2026-10-09 (decisions 9, 13-15); decisions
  27-35 below must be approved before 4.C0 starts. 4.C1 does not ask for the
  extra again.
- Three environments, one `.venv`: `uv sync --all-groups` (exact) gives CI's
  default environment; `uv run --group test --extra torch ...` adds torch;
  `uv run --group test --extra mgm ...` adds torch and transformers. Run the
  default gates after `uv sync --all-groups`, the extras last.
- Gate before every commit, in this order (new files staged first):
  ```bash
  export BIOTAPY_DATA_DIR=<scratch pooch cache>
  uv sync --all-groups
  uvx prek run --all-files
  uv run --group test pytest -q -W error::UserWarning
  uv run --group test pytest -q -m "golden or network"
  rm -rf docs/_build docs/generated && uv run --group doc sphinx-build -W -b html docs docs/_build/html
  uv run --group test --extra torch pytest -q -m torch -W error::UserWarning
  uv run --group test --extra mgm pytest -q -m mgm -W error::UserWarning   # from 4.4b on
  uv sync --all-groups
  ```
- Tests reach the code through `bt.ml.embed` and `importlib.metadata`; no
  test imports `biotapy.ml._embed` or `biotapy.ml._mgm`, and no test file
  imports torch or transformers at module level.
- Every `mgm` test may download (MGM's 33 MB wheel, GlobalPatterns); no other
  test does.
- `microformer-mgm` is never installed in biotapy's environment;
  `tests/mgm/export_reference.py` runs only through the `uv run --no-project
  --python 3.11 ...` command in its docstring.

### Slice 4C review focus
1. **An embedding whose rows are not the samples** (Phase 4 review focus 5).
   Expected: a wrong row count, NaN or infinity raises naming the plugin,
   and nothing is written. Tests: 4.4 `test_plugin_with_wrong_row_count_raises`,
   `test_plugin_with_nan_raises`, `test_a_refused_result_is_not_written`, and
   the property `test_a_finite_result_comes_back_unchanged_and_any_other_raises`.
2. **Two packages registering one model name.** Expected: an error naming
   both, never an arbitrary pick. Test: 4.4
   `test_two_packages_registering_one_name_raise`.
3. **MGM reading little of the table without saying so.** Expected: one
   warning counting the features left out (no genus, or a genus MGM does not
   know) and one naming the samples embedded from `<bos> <eos>` alone. Tests:
   4.4b `test_matches_mgm_s_own_embedding` (both messages, exactly),
   `test_embeds_global_patterns_end_to_end` (96 of 996).
4. **MGM's tokens built differently from MGM** (denominator, sort, kept
   genera, pooling, dropout). Expected: MGM's own embeddings within 1e-3.
   Test: 4.4b `test_matches_mgm_s_own_embedding` (mutation-checked, six
   mutations).
5. **A table that is not at genus level.** Expected: features of one genus
   summed, so an OTU table and its `tax_glom` give the same embedding; the
   order of features and the scale of a sample change nothing. Tests: 4.4b
   `test_features_of_one_genus_are_summed_so_tax_glom_changes_nothing`,
   `test_feature_order_does_not_matter`,
   `test_relative_abundances_embed_as_their_counts`.

---

### Task 4.C0: one download cache, `_core.make_pooch` (`refactor(core)`)

**Files:** create `src/biotapy/_core/_download.py`, `tests/core/test_download.py`;
modify `src/biotapy/_core/__init__.py`, `src/biotapy/datasets/_remote.py`,
`.knowledge/modules/core.md`, `.knowledge/roadmap/phase-4-ml-multiomics.md`,
`.knowledge/log.md`.
**Not touched:** every `tests/datasets` test (they pass unedited: the pooch
they get is the same); `datasets/_enzyme.py`'s docstring (it names
`pooch.os_cache("biotapy")` and `BIOTAPY_DATA_DIR`, still true);
`.knowledge/modules/datasets.md` (`_remote.py:_pooch` is still the `@cache`d
builder it describes); the CI cache keys (they hash `_remote.py`, so the
network job's cache is rebuilt once, which is harmless).
**Interfaces:**
- Consumes: pooch 1.9.0 (`pooch.create(path, base_url, registry=, urls=,
  env=)`, `pooch.os_cache`).
- Produces: `biotapy._core.make_pooch(base_url: str, registry: dict[str, str
  | None], *, urls: dict[str, str] | None = None) -> pooch.Pooch`, whose
  cache is `BIOTAPY_DATA_DIR` if set, else `pooch.os_cache("biotapy")`.
  4.4b's `ml/_mgm.py` calls it.

- [x] **Step 1: Failing test.** Create `tests/core/test_download.py` (a
  `_core` unit test, R4.9's exception):
  ```python
  from pathlib import Path

  import pooch

  from biotapy._core import make_pooch


  def test_the_cache_is_biotapy_data_dir_when_it_is_set(monkeypatch, tmp_path):
      monkeypatch.setenv("BIOTAPY_DATA_DIR", str(tmp_path))
      cache = make_pooch("https://example.org/data/", {"a.txt": None}, urls={"a.txt": "https://example.org/a.txt"})
      assert Path(cache.abspath) == tmp_path
      assert cache.registry == {"a.txt": None} and cache.get_url("a.txt") == "https://example.org/a.txt"


  def test_without_biotapy_data_dir_the_cache_is_pooch_s_per_user_directory(monkeypatch):
      monkeypatch.delenv("BIOTAPY_DATA_DIR", raising=False)
      cache = make_pooch("https://example.org/data/", {})
      assert Path(cache.abspath) == Path(pooch.os_cache("biotapy"))
  ```
  The second test creates no directory: `pooch.create` makes the cache
  directory only when `fetch` runs (`pooch/utils.py:make_local_storage`).
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest
  tests/core/test_download.py -q` -> `1 error` during collection: `E
  ImportError: cannot import name 'make_pooch' from 'biotapy._core'`.
- [x] **Step 3: Implement.** Create `src/biotapy/_core/_download.py`:
  ```python
  """The one download cache: every file biotapy fetches at run time goes through a pooch made here."""

  from typing import cast

  import pooch


  def make_pooch(base_url: str, registry: dict[str, str | None], *, urls: dict[str, str] | None = None) -> pooch.Pooch:
      """A pooch over biotapy's cache: ``BIOTAPY_DATA_DIR`` if set, else pooch's per-user cache directory."""
      # pooch ships no py.typed marker, so mypy --strict infers Any for the untyped `create`; cast it back to Pooch.
      return cast(
          pooch.Pooch,
          pooch.create(
              path=pooch.os_cache("biotapy"), base_url=base_url, registry=registry, urls=urls, env="BIOTAPY_DATA_DIR"
          ),
      )
  ```
  ```diff
  diff --git a/src/biotapy/_core/__init__.py b/src/biotapy/_core/__init__.py
  @@ -1,6 +1,7 @@
   """Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

   from ._composition import check_pseudocount, pseudocounted
  +from ._download import make_pooch
   from ._function import (
       BY_TAXON_KEY,
       FUNCTION_KEY,
  @@ -64,6 +65,7 @@ __all__ = [
       "import_optional",
       "infer_x_kind",
       "make_function_mudata",
  +    "make_pooch",
       "make_treedata",
       "normalize_ranks",
       "pseudocounted",
  diff --git a/src/biotapy/datasets/_remote.py b/src/biotapy/datasets/_remote.py
  @@ -1,11 +1,10 @@
   """Datasets downloaded once and cached with pooch: phyloseq's examples, the ENZYME files and the HMP2 tables."""

   from functools import cache
  -from typing import cast

   import pooch

  -from biotapy._core import TreeData
  +from biotapy._core import TreeData, make_pooch
   from biotapy.io import read_phyloseq

   # Pinned to one phyloseq commit so the SHA-256 hashes stay valid.
  @@ -35,14 +34,7 @@ _URLS = {

   @cache
   def _pooch() -> pooch.Pooch:
  -    # BIOTAPY_DATA_DIR overrides the per-user cache directory. pooch ships no py.typed
  -    # marker, so mypy --strict infers Any for the untyped `create`; cast it back to Pooch.
  -    return cast(
  -        pooch.Pooch,
  -        pooch.create(
  -            path=pooch.os_cache("biotapy"), base_url=_BASE_URL, registry=_REGISTRY, urls=_URLS, env="BIOTAPY_DATA_DIR"
  -        ),
  -    )
  +    return make_pooch(_BASE_URL, _REGISTRY, urls=_URLS)
  ```
  R2.1: pooch is the library; `make_pooch` only fixes biotapy's directory
  and variable in one place, which R4.3 requires once `ml` needs the cache
  too (4.4b).
- [x] **Step 4: Run, expect pass** - `uv run --group test pytest
  tests/core/test_download.py tests/datasets -q` -> `32 passed, 6
  deselected` (the 6 are the `network` tests, which the gate's `-m "golden
  or network"` run covers: every dataset is read through `make_pooch`
  there, `36 passed`).
- [x] **Step 5: Knowledge.** (`generated` is `claude-code/<model>` at the
  commit time; `commit:` the parent's short sha.)
  ```diff
  diff --git a/.knowledge/modules/core.md b/.knowledge/modules/core.md
  @@ -18,8 +18,9 @@ and lineage/rank-column parsing (`_taxonomy.py`), the `x_kind`/provenance/
   construction and the HUMAnN special-row constants (`_function.py`), TreeData construction and the only
   import of `treedata`/`networkx` in the package (`_tree.py`), the single
   user-facing warning entry point (`_warnings.py`), lazy optional-dependency
  -import (`_optional.py`), and the single RNG entry point (`_rng.py`).
  +import (`_optional.py`), the single RNG entry point (`_rng.py`), and the one
  +download cache (`_download.py`).

  @@ -135,6 +136,12 @@ none of them back.
   - `_rng.py:as_generator` - the single entry point that turns a seed into a
     `np.random.Generator` without touching global RNG state.
  +- `_download.py:make_pooch` - the `pooch.Pooch` every run-time download goes
  +  through: `BIOTAPY_DATA_DIR` if set, else `pooch.os_cache("biotapy")`.
  +  Used by `datasets/_remote.py` (phyloseq, ENZYME and HMP2 files); every
  +  later download goes through it too, so all of them land in one cache that
  +  CI and tests point elsewhere with one variable. Building the pooch
  +  downloads and creates nothing; `fetch` does.
  ```
  Tick 4.C0 in this concept. In `.knowledge/log.md`, below the title and
  above the newest section:
  ```markdown
  ## <date> (Phase 4, slice 4C)
  - **Update**: [core](modules/core.md) owns `_download.py:make_pooch`, the one download cache (`BIOTAPY_DATA_DIR` or pooch's per-user cache), moved out of `datasets/_remote.py` so `ml`'s MGM weights (4.4b) use the same cache; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.C0.
  ```
- [x] **Step 6: Gate and commit**
  ```bash
  git add src/biotapy/_core/_download.py src/biotapy/_core/__init__.py src/biotapy/datasets/_remote.py \
    tests/core/test_download.py .knowledge/modules/core.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  # the slice gate without its -m mgm line (no marker yet)
  git commit -m "refactor(core): one download cache, make_pooch, for datasets and ml

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1439 passed, 2 skipped, 79 deselected`; `36
  passed, 1484 deselected`; `build succeeded`; `25 passed, 1495
  deselected`.

### Task 4.C1: the extra `mgm` (`build`)

**Files:** modify `pyproject.toml`, `tests/test_ci.py`,
`.knowledge/decisions/optional-heavy-dependencies.md`,
`.knowledge/decisions/index.md`, `.knowledge/roadmap/phase-4-ml-multiomics.md`,
`.knowledge/log.md`.
**Not touched:** `[tool.uv]` (its `sources.torch` already sends torch to the
CPU index for any extra that needs it); `uv.lock` (git-ignored); the entry
point (4.4b, with the module it names); `README.md` and `CHANGELOG.md`
(4.D2); the `import-without-extras` job.
**Interfaces:**
- Consumes: the `pytorch-cpu` index (4.B0).
- Produces: `pip install 'biotapy[mgm]'` (torch >= 2.9, transformers >= 5);
  `uv run --extra mgm` installs torch from the CPU index and transformers
  from PyPI.

- [x] **Step 1: Failing test.** In `tests/test_ci.py`, before
  `test_coverage_below_90_percent_fails_the_test_job`:
  ```python
  def test_the_mgm_extra_is_torch_and_transformers_with_torch_from_the_cpu_index():
      pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
      assert pyproject["project"]["optional-dependencies"]["mgm"] == ["torch>=2.9", "transformers>=5"]
      assert pyproject["tool"]["uv"]["sources"] == {"torch": {"index": "pytorch-cpu"}}
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest
  tests/test_ci.py -q` -> `1 failed, 21 passed` (`KeyError: 'mgm'`).
- [x] **Step 3: Implement.**
  ```diff
  diff --git a/pyproject.toml b/pyproject.toml
  @@ -44,6 +44,8 @@ dependencies = [
     "treedata>=0.3.1,<0.4",
     "xarray",
   ]
  +# MGM, ml.embed's reference model (decisions/optional-heavy-dependencies): a GPT-2 that transformers runs on torch.
  +optional-dependencies.mgm = [ "torch>=2.9", "transformers>=5" ]
   # The ALDEx2 and MaAsLin 3 bridges (decisions/optional-heavy-dependencies); R and the R packages are the user's install.
   optional-dependencies.r = [ "rpy2>=3.6.8" ]
  ```
  (pyproject-fmt keeps the extras sorted, so `mgm` goes before `r`.)
- [x] **Step 4: Run, expect pass, and record the lock** - `uv run --group
  test pytest tests/test_ci.py -q` -> `22 passed`. `uv lock` -> 211 packages
  (195 before; `grep -c '^name = ' uv.lock`); the 16 new ones are
  annotated-doc, anyio, h11, hf-xet, httpcore2, httpx2, httpx2-jsfetch,
  huggingface-hub, regex, safetensors, shellingham, tokenizers, tqdm,
  transformers, truststore, typer, and no other version moves; `grep -A1
  '^name = "transformers"' uv.lock` -> the newest on the day (5.19.0 on the
  prototype), record it. `uv build --wheel`, then `unzip -p dist/*.whl
  '*/METADATA' | grep -E "Provides-Extra|Requires-Dist: (torch|transformers)"`
  -> `Provides-Extra: mgm`, `Requires-Dist: torch>=2.9; extra == 'mgm'`,
  `Requires-Dist: transformers>=5; extra == 'mgm'`, `Provides-Extra: r`,
  `Provides-Extra: torch`, `Requires-Dist: torch>=2.9; extra == 'torch'`;
  delete `dist/`.
- [x] **Step 5: Knowledge.** (`generated` and `commit:` as in 4.C0.)
  ```diff
  diff --git a/.knowledge/decisions/optional-heavy-dependencies.md b/.knowledge/decisions/optional-heavy-dependencies.md
  @@ -1,7 +1,7 @@
   ---
   type: Decision
   title: Heavy dependencies are optional extras
  -description: torch, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily; `pip install biotapy` stays light.
  +description: torch, transformers, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily; `pip install biotapy` stays light.
  @@ -74,13 +74,27 @@ torch or an R installation into every install is unacceptable.[^spec]
     Linux wheel is 196 MB, against PyPI's 555 MB plus CUDA packages. The index
     is uv configuration only: the wheel's metadata says `torch>=2.9; extra ==
  -  'torch'`, and pip users get PyPI's torch.
  +  'torch'`, and pip users get PyPI's torch. Phase 4 task 4.C1 added the extra
  +  `mgm` (`torch>=2.9`, `transformers>=5`, approved 2026-10-09) for MGM, the
  +  reference model of `ml.embed`: a GPT-2 that transformers (Apache-2.0) runs
  +  on torch, whose weights are downloaded at run time and never bundled.
  +  transformers brings huggingface-hub, tokenizers, safetensors, regex, tqdm
  +  and typer with their own dependencies; 16 packages are new to `uv.lock`
  +  (211 against 195) and no version moves. They take about 97 MB installed
  +  beside torch's 742 MB (Python 3.14, Linux). Every binary among them ships
  +  abi3 or per-version wheels for CPython 3.12-3.14 on Linux, macOS and
  +  Windows; torch still comes from the CPU index under uv. 5.0.0 is the first
  +  transformers 5 release (2026-01-26); MGM's plugin was run on 5.0.0 with
  +  torch 2.9.0 and on 5.19.0 with torch 2.14.1. `microformer-mgm` itself is
  +  not a dependency: it pins numpy 1.24, pandas 2.0, torch 2.0 and
  +  transformers 4.33.
   - Extras (names fixed now so docs never change), each added in the phase that first uses it:

     | Extra | Pulls | First used |
     |---|---|---|
     | `numba` | numba (>=0.67, supports up to Python 3.14) | perf track, only on benchmark evidence; also unlocks scikit-bio's `engine="numba"` |
  +  | `mgm` | torch (>=2.9), transformers (>=5) | Phase 4 `ml.embed`'s MGM plugin |
     | `r` | rpy2 (3.6.8) | Phase 3 `da` bridges |
  ```
  Copy the new description into
  `.knowledge/decisions/index.md`'s line (it keeps its shorter ending, "are
  imported lazily."). Tick 4.C1 here; log line, first in the slice 4C
  section:
  ```markdown
  - **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md): the extra `mgm` (`torch>=2.9`, `transformers>=5`), what transformers brings, its size and wheels, the versions MGM was run on, why `microformer-mgm` is not a dependency, the extras table row, and transformers in the description, copied into the [decisions index](decisions/index.md); [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.C1.
  ```
- [x] **Step 6: Gate and commit**
  ```bash
  git add pyproject.toml tests/test_ci.py .knowledge/decisions/optional-heavy-dependencies.md \
    .knowledge/decisions/index.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  # the slice gate without its -m mgm line (no marker yet)
  git commit -m "build: add the mgm extra, torch and transformers, for MGM embeddings

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1440 passed, 2 skipped, 79 deselected`; `36
  passed, 1485 deselected`; `build succeeded`; `25 passed, 1496
  deselected`.

### Task 4.4: `ml.embed` and the entry-point group `biotapy.embeddings`

**Files:** create `src/biotapy/ml/_embed.py`, `tests/ml/test_embed.py`,
`.knowledge/decisions/embedding-plugins.md`; modify
`src/biotapy/ml/__init__.py`, `conftest.py`, `pyproject.toml`,
`docs/guide/machine_learning.md`, `docs/api.md`, `docs/contributing.md`,
`.knowledge/decisions/index.md`, `.knowledge/contracts/data-model-slots.md`,
`.knowledge/decisions/pure-by-default.md`, `.knowledge/modules/ml.md`,
`.knowledge/modules/index.md`, `.knowledge/roadmap/phase-4-ml-multiomics.md`,
`.knowledge/log.md`.
**Not touched:** `ml/_torch.py`, `ml/_transformers.py`; `_core` (nothing here
is shared with another subpackage, R4.3); `tests/conftest.py`
(the doctest hook must see `src/biotapy`, which only the root `conftest.py`
does); `docs/extensions/coming_from_r.py` (`R equivalent: none` adds no row);
`.github/workflows/test.yaml` (4.4b adds the `-m mgm` step, once MGM exists);
`function-shape` (an `ml` function with the `tl` signature needs no new
sentence).
**Interfaces:**
- Consumes: `importlib.metadata.entry_points` (standard library); fixtures
  `assert_unchanged`, `make_adata` (`tests/conftest.py`); the extra `mgm`
  (4.C1), which the new marker names.
- Produces: `bt.ml.embed(adata: AnnData, model: str, *, inplace: bool =
  False) -> npt.NDArray[np.floating] | None`; the plugin protocol: an entry
  point in the group `biotapy.embeddings` whose value loads `embed(adata:
  AnnData) -> numpy.ndarray`; the pytest marker `mgm`, deselected by default,
  which the root `conftest.py` gives to `biotapy.ml._embed`'s doctests.
  4.4b registers `mgm`.

- [x] **Step 1: Failing tests.** Register the marker and deselect it by
  default, and give `ml.embed`'s doctest the marker (its example runs MGM,
  which 4.4b registers; decision 34):
  ```diff
  diff --git a/pyproject.toml b/pyproject.toml
  @@ -232,16 +232,17 @@ overrides = [

   [tool.pytest]
   addopts = [
  -  "--import-mode=importlib",             # allow using test files with same name
  +  "--import-mode=importlib",                         # allow using test files with same name
     "--doctest-modules",
     "-m",
  -  "not network and not r and not torch",
  +  "not network and not r and not torch and not mgm",
   ]
   markers = [
     "golden: compares against R or HUMAnN golden files (contracts/r-golden-parity)",
     "network: downloads data; runs only in the dedicated CI job",
     "r: needs R and rpy2 (extra `r`)",
     "torch: needs PyTorch (extra `torch`); runs only in the ml-extras CI job",
  +  "mgm: needs the extra `mgm` and downloads MGM's weights; runs only in the ml-extras CI job",
   ]
   strict = true
   testpaths = [ "tests", "src/biotapy" ]
  diff --git a/conftest.py b/conftest.py
  @@ -26,9 +26,15 @@ def _close_figures() -> Iterator[None]:
       plt.close("all")


  +# Doctests that need an extra, by module: their examples run only where the marker's CI job installs it.
  +_EXTRA_DOCTESTS = {"biotapy.ml._torch.": pytest.mark.torch, "biotapy.ml._embed.": pytest.mark.mgm}
  +
  +
   @pytest.hookimpl(tryfirst=True)
   def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
  -    """Give the doctests of a module that needs the extra `torch` its marker, before `-m` deselects."""
  +    """Give the doctests of a module that needs an extra its marker, before `-m` deselects."""
       for item in items:
  -        if isinstance(item, pytest.DoctestItem) and item.name.startswith("biotapy.ml._torch."):
  -            item.add_marker(pytest.mark.torch)
  +        if isinstance(item, pytest.DoctestItem):
  +            for prefix, marker in _EXTRA_DOCTESTS.items():
  +                if item.name.startswith(prefix):
  +                    item.add_marker(marker)
  ```
  (The `--import-mode` line is pyproject-fmt's realignment again.) Create
  `tests/ml/test_embed.py`:
  ```python
  import sys
  import types

  import numpy as np
  import pytest
  from hypothesis import given
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays

  import biotapy as bt

  GROUP = "biotapy.embeddings"


  def _installer(patch, root):
      """Install embedding plugins as a package would: a distribution's entry points under ``root`` on sys.path."""
      module = types.ModuleType(f"biotapy_test_plugins_{root.name}")
      patch.setitem(sys.modules, module.__name__, module)
      patch.syspath_prepend(str(root))
      registered = {}

      def install(name, plugin, *, distribution="biotapy-test-plugins"):
          attribute = f"plugin_{len(vars(module))}"
          setattr(module, attribute, plugin)
          registered.setdefault(distribution, []).append(f"{name} = {module.__name__}:{attribute}")
          info = root / f"{distribution.replace('-', '_')}-1.0.dist-info"
          info.mkdir(exist_ok=True)
          (info / "METADATA").write_text(f"Metadata-Version: 2.1\nName: {distribution}\nVersion: 1.0\n", encoding="utf-8")
          lines = "\n".join(registered[distribution])
          (info / "entry_points.txt").write_text(f"[{GROUP}]\n{lines}\n", encoding="utf-8")

      return install


  @pytest.fixture
  def install(tmp_path, monkeypatch):
      return _installer(monkeypatch, tmp_path)


  @pytest.fixture(scope="module")
  def echo(tmp_path_factory):
      """A plugin that returns ``uns["echo"]``, installed once per module: Hypothesis rejects function-scoped fixtures."""
      with pytest.MonkeyPatch.context() as patch:
          _installer(patch, tmp_path_factory.mktemp("echo"))("echo", lambda adata: adata.uns["echo"])
          yield


  def _ones(adata):
      return np.ones((adata.n_obs, 3))


  def test_returns_the_plugin_s_embedding(install):
      install("fake", lambda adata: np.arange(adata.n_obs * 2, dtype=float).reshape(-1, 2))
      tdata = bt.datasets.toy()
      result = bt.ml.embed(tdata, "fake")
      np.testing.assert_array_equal(result, np.arange(12, dtype=float).reshape(6, 2))
      assert "X_fake" not in tdata.obsm


  def test_inplace_writes_obsm_x_model_and_returns_none(install):
      install("fake", _ones)
      tdata = bt.datasets.toy()
      assert bt.ml.embed(tdata, "fake", inplace=True) is None
      np.testing.assert_array_equal(tdata.obsm["X_fake"], np.ones((6, 3)))


  def test_the_plugin_gets_the_data_itself(install):
      calls = []
      install("fake", lambda adata: calls.append(adata) or np.ones((adata.n_obs, 1)))
      tdata = bt.datasets.toy()
      bt.ml.embed(tdata, "fake")
      assert len(calls) == 1 and calls[0] is tdata


  def test_unknown_model_lists_the_installed_ones(install):
      install("fake", _ones)
      install("other", _ones)
      with pytest.raises(
          KeyError, match=r"model='nope' is not an installed embedding plugin; installed: \[.*'fake', .*'other'"
      ):
          bt.ml.embed(bt.datasets.toy(), "nope")


  @pytest.mark.parametrize("model", ["", "a/b", "with space", "x=y", "../up"])
  def test_a_model_name_that_is_not_an_obsm_key_raises(install, model):
      install("fake", _ones)
      with pytest.raises(ValueError, match=r"model=.* must be letters, digits, '_', '-' or '\.'"):
          bt.ml.embed(bt.datasets.toy(), model)


  @pytest.mark.parametrize("model", ["fake", "My.Model_2", "mgm-v2"])
  def test_a_model_name_of_letters_digits_and_separators_is_accepted(install, model):
      install(model, _ones)
      assert bt.ml.embed(bt.datasets.toy(), model).shape == (6, 3)


  def test_a_package_registering_one_name_twice_is_named_once(install):
      install("fake", _ones, distribution="only-package")
      install("fake", _ones, distribution="only-package")
      with pytest.raises(ValueError, match=r"several packages \(only-package\); uninstall all but one"):
          bt.ml.embed(bt.datasets.toy(), "fake")


  def test_two_packages_registering_one_name_raise(install):
      install("fake", _ones, distribution="first-package")
      install("fake", _ones, distribution="second-package")
      with pytest.raises(
          ValueError, match=r"model='fake' is registered by several packages \(first-package, second-package\)"
      ):
          bt.ml.embed(bt.datasets.toy(), "fake")


  def test_plugin_with_wrong_row_count_raises(install):
      install("fake", lambda adata: np.ones((adata.n_obs - 1, 3)))
      with pytest.raises(ValueError, match=r"plugin 'fake' returned shape \(5, 3\); expected \(6, dimensions\)"):
          bt.ml.embed(bt.datasets.toy(), "fake")


  @pytest.mark.parametrize("bad", [np.nan, np.inf])
  def test_plugin_with_nan_raises(install, bad):
      def plugin(adata):
          out = np.ones((adata.n_obs, 3))
          out[2, 1] = bad
          return out

      install("fake", plugin)
      with pytest.raises(ValueError, match="plugin 'fake' returned NaN or infinite values"):
          bt.ml.embed(bt.datasets.toy(), "fake")


  @pytest.mark.parametrize("shape", [(6,), (6, 0), (6, 2, 2)])
  def test_a_result_that_is_not_samples_by_dimensions_raises(install, shape):
      install("fake", lambda adata: np.ones(shape))
      with pytest.raises(ValueError, match=r"plugin 'fake' returned shape"):
          bt.ml.embed(bt.datasets.toy(), "fake")


  def test_an_integer_result_raises(install):
      install("fake", lambda adata: np.ones((adata.n_obs, 3), dtype=np.int64))
      with pytest.raises(ValueError, match="plugin 'fake' returned int64 values; expected floats"):
          bt.ml.embed(bt.datasets.toy(), "fake")


  def test_a_result_that_is_not_an_array_raises(install):
      install("fake", lambda adata: [[1.0, 2.0]] * adata.n_obs)
      with pytest.raises(TypeError, match="plugin 'fake' returned a list, not a NumPy array"):
          bt.ml.embed(bt.datasets.toy(), "fake")


  @pytest.mark.parametrize(
      ("make", "kind"),
      [
          (lambda n: np.ma.masked_invalid(np.full((n, 3), np.nan)), "MaskedArray"),
          (lambda n: np.ones((n, 3)).view(np.matrix), "matrix"),
      ],
  )
  def test_an_array_subclass_raises(install, make, kind):
      install("fake", lambda adata: make(adata.n_obs))
      with pytest.raises(TypeError, match=f"plugin 'fake' returned a {kind}, not a NumPy array"):
          bt.ml.embed(bt.datasets.toy(), "fake")


  @pytest.mark.parametrize("where", ["X", "layer", "obsm", "uns"])
  def test_a_result_that_shares_memory_with_the_input_is_copied(install, where):
      tdata = bt.datasets.toy()
      tdata.X = np.ones((tdata.n_obs, tdata.n_vars))
      tdata.layers["own"] = np.ones((tdata.n_obs, tdata.n_vars))
      tdata.obsm["own"] = np.ones((tdata.n_obs, 3))
      tdata.uns["own"] = np.ones((tdata.n_obs, 3))
      held = {
          "X": lambda a: a.X,
          "layer": lambda a: a.layers["own"],
          "obsm": lambda a: a.obsm["own"],
          "uns": lambda a: a.uns["own"],
      }
      install("fake", held[where])
      result = bt.ml.embed(tdata, "fake")
      assert not np.shares_memory(result, held[where](tdata))
      np.testing.assert_array_equal(result, held[where](tdata))


  def test_a_plugin_that_raises_keeps_its_error_and_is_named(install):
      def plugin(adata):
          msg = "model exploded"
          raise RuntimeError(msg)

      install("fake", plugin)
      with pytest.raises(RuntimeError, match="model exploded") as info:
          bt.ml.embed(bt.datasets.toy(), "fake")
      assert any("embedding plugin 'fake'" in note for note in info.value.__notes__)


  def test_a_plugin_that_is_not_callable_is_named(install):
      install("fake", 3)
      with pytest.raises(TypeError, match="not callable") as info:
          bt.ml.embed(bt.datasets.toy(), "fake")
      assert any("embedding plugin 'fake'" in note for note in info.value.__notes__)


  def test_a_plugin_that_cannot_be_loaded_is_named(install, tmp_path):
      install("fake", _ones)
      module = sys.modules[f"biotapy_test_plugins_{tmp_path.name}"]
      for attribute in [name for name in vars(module) if name.startswith("plugin_")]:
          delattr(module, attribute)
      with pytest.raises(AttributeError) as info:
          bt.ml.embed(bt.datasets.toy(), "fake")
      assert any("embedding plugin 'fake'" in note for note in info.value.__notes__)


  def test_a_refused_result_is_not_written(install):
      install("fake", lambda adata: np.full((adata.n_obs, 3), np.nan))
      tdata = bt.datasets.toy()
      with pytest.raises(ValueError, match="NaN"):
          bt.ml.embed(tdata, "fake", inplace=True)
      assert "X_fake" not in tdata.obsm


  def test_keeps_the_input(install, assert_unchanged):
      install("fake", _ones)
      tdata = bt.datasets.toy()
      before = tdata.copy()
      bt.ml.embed(tdata, "fake")
      assert_unchanged(before, tdata)


  def test_a_feature_change_drops_the_embedding(install):
      install("fake", _ones)
      tdata = bt.datasets.toy()
      bt.ml.embed(tdata, "fake", inplace=True)
      assert "X_fake" not in bt.pp.filter_features(tdata, min_prevalence=0.5).obsm


  def test_single_sample_and_all_zero_sample(install, make_adata):
      install("fake", lambda adata: np.asarray(adata.X.sum(axis=1), dtype=float))
      np.testing.assert_array_equal(bt.ml.embed(make_adata(np.array([[0, 0, 0]])), "fake"), [[0.0]])
      np.testing.assert_array_equal(bt.ml.embed(make_adata(np.array([[0, 0], [2, 3]])), "fake"), [[0.0], [5.0]])


  def test_options_are_keyword_only(install):
      install("fake", _ones)
      with pytest.raises(TypeError):
          bt.ml.embed(bt.datasets.toy(), "fake", True)


  @given(
      arrays(
          np.float64,
          st.tuples(st.just(6), st.integers(1, 5)),
          elements=st.floats(allow_nan=True, allow_infinity=True) | st.floats(-1e6, 1e6),
      )
  )
  def test_a_finite_result_comes_back_unchanged_and_any_other_raises(echo, result):
      tdata = bt.datasets.toy()
      tdata.uns["echo"] = result
      if np.isfinite(result).all():
          np.testing.assert_array_equal(bt.ml.embed(tdata, "echo"), result)
      else:
          with pytest.raises(ValueError, match="plugin 'echo' returned NaN or infinite values"):
              bt.ml.embed(tdata, "echo")
  ```
  R11.2 coverage: happy path (returned, `inplace=True`, the plugin receives
  the AnnData itself); edge cases (all-zero sample, single sample; NaN
  taxonomy does not apply to `embed`, which reads none, and is MGM's, 4.4b);
  purity (`test_keeps_the_input`); a Hypothesis property (a finite result
  comes back unchanged, any other raises); no golden test (no R
  equivalent). Plugins are installed as a real distribution on `sys.path`
  (`_installer`), so discovery runs the code path an installed package
  takes and no test imports `ml/_embed.py` (R4.9); the module-scoped `echo`
  fixture exists because Hypothesis rejects function-scoped fixtures.
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest
  tests/ml/test_embed.py -q` -> `19 failed`; every failure is
  `AttributeError: module 'biotapy.ml' has no attribute 'embed'`.
- [x] **Step 3: Implement.** Create `src/biotapy/ml/_embed.py`:
  ```python
  """Sample embeddings from models that plugins register in the entry-point group ``biotapy.embeddings``."""

  import re
  from collections.abc import Iterator
  from importlib.metadata import entry_points
  from typing import cast

  import numpy as np
  import numpy.typing as npt
  import scipy.sparse as sp
  from anndata import AnnData

  # The entry-point group a package registers an embedding model in (decisions/embedding-plugins).
  GROUP = "biotapy.embeddings"


  def embed(adata: AnnData, model: str, *, inplace: bool = False) -> npt.NDArray[np.floating] | None:
      """One embedding per sample from a model that a plugin provides.

      Parameters
      ----------
      adata
          Samples x features, in the form the model reads (``"mgm"``: counts or
          relative abundances with ``var["genus"]``).
      model
          The name a plugin registers in the entry-point group
          ``biotapy.embeddings``: letters, digits, ``_``, ``-`` or ``.``, so that
          ``obsm["X_<model>"]`` is a plain key. biotapy registers ``"mgm"``.
      inplace
          Write the embedding to ``obsm[f"X_{model}"]`` and return ``None``.

      Returns
      -------
      numpy.ndarray or None
          Samples x dimensions, a float array in ``obs`` order.

      Raises
      ------
      KeyError
          No installed plugin registers ``model``; the message lists those that do.
      TypeError
          The plugin returns something other than a plain NumPy array.
      ValueError
          ``model`` has other characters, two installed packages register
          ``model``, or the plugin's array is not 2-D, float and finite with one
          row per sample. Errors about the result name the plugin.
      Exception
          What the plugin itself raises passes through with its type and a note
          that names the plugin. ``"mgm"`` raises ``ImportError`` without the
          extra, ``KeyError`` without ``var["genus"]`` and ``ValueError`` for
          negative or non-finite ``X``.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/machine_learning`

      A plugin is a callable ``embed(adata)`` that returns a samples x
      dimensions NumPy array and leaves ``adata`` unchanged. A package registers
      it under a name in its ``pyproject.toml``::

          [project.entry-points."biotapy.embeddings"]
          mymodel = "mypackage.module:embed"

      biotapy finds it when ``embed`` is called, so installing the package is
      enough. It checks the array before returning or storing it.

      Examples
      --------
      >>> import biotapy as bt
      >>> genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      >>> bt.ml.embed(genera, "mgm").shape
      (6, 256)
      """
      if not re.fullmatch(r"[\w.-]+", model):
          msg = f"model={model!r} must be letters, digits, '_', '-' or '.', so that obsm['X_<model>'] is a plain key"
          raise ValueError(msg)
      found = [point for point in entry_points(group=GROUP) if point.name == model]
      if not found:
          installed = sorted({point.name for point in entry_points(group=GROUP)})
          msg = f"model={model!r} is not an installed embedding plugin; installed: {installed}"
          raise KeyError(msg)
      if len(found) > 1:
          packages = sorted({point.dist.name if point.dist else point.value for point in found})
          msg = f"model={model!r} is registered by several packages ({', '.join(packages)}); uninstall all but one"
          raise ValueError(msg)
      # A plugin is not trusted: a wrong row count would pair embeddings with the wrong samples (decisions/embedding-plugins).
      try:
          result = found[0].load()(adata)
      except Exception as err:
          err.add_note(f"while using the embedding plugin {model!r} ({found[0].value})")
          raise
      result = _checked(model, adata, result)
      if not inplace:
          return result
      adata.obsm[f"X_{model}"] = result
      return None


  def _checked(model: str, adata: AnnData, result: object) -> npt.NDArray[np.floating]:
      """``result`` if it is a float, finite, samples x dimensions array that shares no memory with ``adata``'s X, layers, obsm, varm, obsp, varp or top-level uns, else an error."""
      if type(result) is not np.ndarray:  # a masked array would hide NaN from the check below, a matrix is always 2-D
          msg = f"plugin {model!r} returned a {type(result).__name__}, not a NumPy array"
          raise TypeError(msg)
      if result.ndim != 2 or result.shape[0] != adata.n_obs or result.shape[1] == 0:
          msg = (
              f"plugin {model!r} returned shape {result.shape}; expected ({adata.n_obs}, dimensions), one row per sample"
          )
          raise ValueError(msg)
      if not np.issubdtype(result.dtype, np.floating):
          msg = f"plugin {model!r} returned {result.dtype} values; expected floats"
          raise ValueError(msg)
      if not np.isfinite(result).all():
          msg = f"plugin {model!r} returned NaN or infinite values"
          raise ValueError(msg)
      if any(np.shares_memory(result, held) for held in _arrays(adata)):
          return result.copy()  # editing the embedding must not edit the caller's table
      return result


  def _arrays(adata: AnnData) -> Iterator[npt.NDArray[np.generic]]:
      """The arrays adata holds: X, layers, obsm, varm, obsp, varp and the arrays in uns (a sparse matrix by its data)."""
      for held in (
          adata.X,
          *adata.layers.values(),
          *adata.obsm.values(),
          *adata.varm.values(),
          *adata.obsp.values(),
          *adata.varp.values(),
          *adata.uns.values(),
      ):
          if sp.issparse(held):
              yield cast("sp.csr_matrix", held).data
          elif isinstance(held, np.ndarray):
              yield held
  ```
  No private helper: with the lookup and the checks inline, `embed` has 7
  branches and 2 returns, within R5 (ruff passes), so R4.4 allows none. R2.1:
  `importlib.metadata.entry_points` is the discovery and NumPy the checks;
  nothing checks an embedding for biotapy. `src/biotapy/ml/__init__.py`:
  ```python
  from ._embed import embed
  from ._torch import to_torch
  from ._transformers import CLR, PrevalenceFilter

  __all__ = ["CLR", "PrevalenceFilter", "embed", "to_torch"]
  ```
- [x] **Step 4: Run, expect pass** - `uv run --group test pytest
  tests/ml/test_embed.py src/biotapy/ml -q -W error::UserWarning` -> `39
  passed, 2 deselected` (the 37 tests and the two transformer doctests; the
  deselected are `to_torch`'s and `embed`'s doctests); `--hypothesis-seed=1`,
  `2`, `3` -> `37 passed` each; `uv run --group dev --group doc mypy` ->
  `Success: no issues found in 67 source files`.
- [x] **Step 5: Docs.**
  ````diff
  diff --git a/docs/api.md b/docs/api.md
  @@ -120,6 +120,7 @@ Public functions are listed here as they ship, from Phase 1 onward.

       ml.CLR
       ml.PrevalenceFilter
  +    ml.embed
       ml.to_torch
   ```

  diff --git a/docs/contributing.md b/docs/contributing.md
  @@ -49,7 +49,7 @@ uv run --group test pytest
   ```

   Network and golden tests are excluded by default (`[tool.pytest]` in
  -`pyproject.toml` sets `-m "not network and not r and not torch"`). Run them explicitly:
  +`pyproject.toml` sets `-m "not network and not r and not torch and not mgm"`). Run them explicitly:

   ```bash
   uv run --group test pytest -m "network or golden"
  diff --git a/docs/guide/machine_learning.md b/docs/guide/machine_learning.md
  @@ -2,8 +2,8 @@

   `bt.ml` holds scikit-learn transformers, so microbiome preprocessing can sit
   inside a [Pipeline](https://scikit-learn.org/stable/modules/compose.html) and
  -be fitted on the training samples of each cross-validation fold only, and
  -turns a table into a PyTorch dataset.
  +be fitted on the training samples of each cross-validation fold only, turns
  +a table into a PyTorch dataset, and embeds samples with pretrained models.

   ## Which steps leak

  @@ -109,3 +109,42 @@ CPU-only torch, install it from PyTorch's CPU index first:
   pip install torch --index-url https://download.pytorch.org/whl/cpu
   pip install 'biotapy[torch]'
   ```
  +
  +## Embeddings
  +
  +`bt.ml.embed` runs a pretrained model over every sample and returns one
  +vector per sample, in `obs` order, or with `inplace=True` stores it in
  +`obsm["X_<model>"]`:
  +
  +```python
  +import biotapy as bt
  +
  +genera = bt.pp.tax_glom(bt.datasets.global_patterns(), "genus")
  +bt.ml.embed(genera, "mgm", inplace=True)
  +genera.obsm["X_mgm"].shape  # (26, 256)
  +```
  +
  +A model reads what it was trained on, so filter first: a later feature change
  +(`bt.pp.filter_features`, `bt.pp.tax_glom`) drops `obsm`, and the embedding
  +with it. MGM embeds each sample from that sample alone, so computing its
  +embedding before a cross-validation split leaks nothing. A plugin that
  +normalises across samples would leak; check how yours works.
  +
  +### Adding a model
  +
  +Models are plugins. A package provides a function that takes the AnnData and
  +returns a samples x dimensions NumPy array, leaving the AnnData unchanged, and
  +registers it under a name in the entry-point group `biotapy.embeddings` of its
  +`pyproject.toml`:
  +
  +```toml
  +[project.entry-points."biotapy.embeddings"]
  +mymodel = "mypackage.embedding:embed"
  +```
  +
  +Once the package is installed, `bt.ml.embed(adata, "mymodel")` finds it;
  +nothing is registered by hand. biotapy checks what the plugin returns - a
  +finite float array with one row per sample - and raises an error naming the
  +plugin otherwise. An unknown name raises `KeyError` listing the installed
  +models, and a name that two installed packages register raises instead of
  +picking one.
  ````
  The guide's code blocks are plain Markdown, not executed: the `mgm` model
  arrives in 4.4b, whose end-to-end test runs this very call.
- [x] **Step 6: Knowledge.** Create
  `.knowledge/decisions/embedding-plugins.md`:
  ```markdown
  ---
  type: Decision
  title: Embedding models are plugins found through entry points
  description: ml.embed(adata, model) loads the callable a package registers under model in the entry-point group biotapy.embeddings, checks that it returned a finite 2-D float array with one row per sample, and returns it or writes obsm["X_<model>"]; biotapy's own MGM is registered the same way, its weights downloaded, never bundled.
  tags: [ml, plugins, api]
  status: draft
  paths: ["src/biotapy/ml/_embed.py", "pyproject.toml"]
  generated: { by: claude-code/<model>, at: <UTC commit time> }
  commit: <parent short sha>
  sources:
    - id: spec
      resource: ../../plan.md
      title: Python Microbiome Toolkit development report
      author: human:pedrocr83
    - id: entry-points
      resource: https://packaging.python.org/en/latest/specifications/entry-points/
      title: Entry points specification (PyPA)
  ---

  # Context
  The spec plans foundation-model embeddings in `obsm` through plugins.[^spec]
  Microbiome foundation models come and go (MGM, MGM2, BiomeGPT), each with its
  own framework, weights and licence; biotapy cannot depend on all of them, and
  a model's authors may want to ship it themselves. Python's packaging already
  has a registry for this: entry points, which a package declares in its
  `pyproject.toml` and `importlib.metadata` lists at run time.[^entry-points]

  # Decision
  - **Group** `biotapy.embeddings`. An entry point's name is the model's name;
    its value loads a callable `embed(adata) -> numpy.ndarray` that returns
    samples x dimensions and leaves `adata` unchanged.
  - **Discovery at call time**: `ml.embed(adata, model)` reads
    `importlib.metadata.entry_points(group="biotapy.embeddings")` on every call,
    so installing a package is enough; nothing is imported until a model is
    asked for. An unknown name raises `KeyError` listing the installed names;
    a name two packages register raises `ValueError` naming both, never picks
    one (`ml/_embed.py:embed`).
  - **biotapy checks the result**, not the plugin: a plain `numpy.ndarray`, not
    a masked array or a matrix (else `TypeError`), 2-D with one row per sample and at least one column, float,
    finite (else `ValueError`), each message naming the plugin
    (`ml/_embed.py:_checked`). A refused result is never stored, and one that
    shares memory with `X`, a layer, `obsm`, `varm`, `obsp`, `varp` or a
    top-level `uns` array (nothing nested deeper) is copied. A plugin's own exception
    propagates with its type and a note naming the plugin. The model name is
    letters, digits, `_`, `-` or `.`, so `obsm["X_<model>"]` is a plain key.
  - **The plugin gets the caller's AnnData, not a copy**: a copy would double
    `X` in memory for every call. The protocol therefore forbids a plugin to
    modify it; biotapy does not enforce that (MGM's own test does:
    `tests/ml/test_mgm.py:test_keeps_the_input`).
  - **Where it goes**: returned, or with `inplace=True` written to
    `obsm[f"X_{model}"]` with `None` returned, the `tl` convention
    ([pure-by-default](/decisions/pure-by-default.md),
    [data-model-slots](/contracts/data-model-slots.md)).
  - **No options pass through**: the callable takes the AnnData alone (rules.md
    R3.7). The reference model, MGM, ran fastest one sample at a time on the CPU
    (phase-4 plan, task 4.4b), so no `batch_size` exists until a model needs it.
  - **Weights are never bundled** (rules.md R6.6): a plugin downloads them at run
    time; biotapy's own go through `_core.make_pooch`, pinned by SHA-256.
  - **biotapy's reference model registers itself the same way**:
    `[project.entry-points."biotapy.embeddings"] mgm = "biotapy.ml._mgm:embed"`
    in biotapy's `pyproject.toml`, behind the extra `mgm`
    ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).

  # Rejected
  - **A registry function** (`bt.ml.register(name, fn)`): needs an import with
    a side effect before every use, and two packages would race for a name.
  - **A base class plugins inherit from**: rules.md R3.6 allows classes only
    where an external protocol requires one; a callable is enough.
  - **Passing `**kwargs` to the plugin**: rules.md R3.7.
  - **Trusting the plugin's output**: a wrong row count would silently pair
    embeddings with the wrong samples.
  - **A separate `biotapy-mgm` package**: no release cadence of its own yet
    (Phase 4 decision 9); the entry point keeps that move cheap.

  # Consequences
  - Third-party models need no change to biotapy; the guide shows the
    `pyproject.toml` lines (`docs/guide/machine_learning.md`, Embeddings).
  - Tests install fake plugins as a distribution under `tmp_path` on
    `sys.path`, so they go through the same discovery
    (`tests/ml/test_embed.py:_installer`).
  - A feature-changing step drops `obsm`, and with it `X_<model>`
    ([data-model-slots](/contracts/data-model-slots.md), Propagation).

  [^spec]: Python Microbiome Toolkit development report, section Roadmap
  [^entry-points]: Entry points specification (PyPA)
  ```
  List it in `.knowledge/decisions/index.md`, above "Pure by default", with
  the description as its text. Then (`generated` and `commit:` as in 4.C0;
  the frontmatter hunks are not shown; data-model-slots also gains
  `"src/biotapy/ml/**"` in `paths`):
  ```diff
  diff --git a/.knowledge/contracts/data-model-slots.md b/.knowledge/contracts/data-model-slots.md
  @@ -32,7 +32,7 @@ Extends the spec's data-model table with exact keys.[^spec]
   | `obs` | sample metadata; `tl` per-sample results with `inplace=True` | `alpha_<metric>` (e.g. `alpha_shannon`) |
   | `var` | taxonomy, one lowercase column per rank; sequences; QIIME 2 assignment confidence | ranks from `kingdom, phylum, class, order, family, genus, species`; `sequence`; `confidence` (float, from a QIIME 2 `FeatureData[Taxonomy]` artifact's `Confidence` column); function tables: see Function tables |
   | `vart` | phylogeny as `networkx.DiGraph`, leaves = `var_names`, edge attribute `length` | `phylo` only |
  -| `obsm` | ordinations and embeddings | `X_pcoa`, `X_nmds`; `X_philr` from `pp.philr` (a samples x balances `DataFrame`, one column per internal tree node with two children, named after the node in `vart["phylo"]`, in preorder; R's names differ, so match balances across tools by their taxa); `X_<plugin>` |
  +| `obsm` | ordinations and embeddings | `X_pcoa`, `X_nmds`; `X_philr` from `pp.philr` (a samples x balances `DataFrame`, one column per internal tree node with two children, named after the node in `vart["phylo"]`, in preorder; R's names differ, so match balances across tools by their taxa); `X_<model>` from `ml.embed(adata, model, inplace=True)`, `model` being the plugin's entry-point name (`X_mgm`; [embedding-plugins](/decisions/embedding-plugins.md)) |
   | `obsp` | sample-sample distance matrices | metric name: `braycurtis`, `jaccard`, `unweighted_unifrac`, `weighted_unifrac` |
   | `uns["biotapy"]` | biotapy metadata, nothing else | `x_kind`, `provenance`, `pcoa` (`eigenvalues`, `proportion_explained`), `nmds` (`stress`) |

  @@ -167,7 +167,7 @@ raise, a categorical `group` has exactly two levels.
   | Feature-changing (`pp.filter_features`, `pp.tax_glom`, `fn.func_glom`, `fn.renorm`, `pp.rarefy`) | `obs`, `var` rows kept, `vart` (pruned by TreeData), `uns["biotapy"]["x_kind"]` and `["provenance"]` | all `layers`, `obsm`, `obsp`, `varm`, `varp`, `uns["biotapy"]["pcoa"]`, `["nmds"]`, other `uns` keys |
   | Sample-only (`pp.filter_samples`) | everything, subset by AnnData indexing; a kept `obsm` ordination and its `pcoa`/`nmds` summary still reflect the dropped samples, so recompute them | nothing |
   | Layer-adding (`pp.relative`, `pp.clr`) | everything | nothing; adds one layer |
  -| Embedding-adding (`pp.philr`) | everything | nothing; adds `obsm["X_philr"]`, which a later feature change drops |
  +| Embedding-adding (`pp.philr`; `ml.embed` with `inplace=True`) | everything | nothing; adds `obsm["X_philr"]` or `obsm["X_<model>"]`, which a later feature change drops |

   Feature-changing operations go through `_core.feature_subset`, or `_core.replace_features` when the new
   features are groups rather than a subset; both keep only `KEPT_META`, the
  diff --git a/.knowledge/decisions/pure-by-default.md b/.knowledge/decisions/pure-by-default.md
  @@ -31,6 +31,8 @@ behaviour. Confirmed by the user on 2026-09-26.
   | `pl` | `matplotlib.axes.Axes` | never |
   | `ml` estimators (`PrevalenceFilter`, `CLR`) | scikit-learn's protocol: `fit` stores what it learns on the estimator and returns it; `transform` returns a new array | never the data |
   | `ml.to_torch` | a new `torch.utils.data.Dataset` that holds `X` (or the layer) without copying it; each item a new tensor | never |
  +| `ml.embed` (default `inplace=False`) | the embedding (`np.ndarray`), as `tl` does | never |
  +| `ml.embed` with `inplace=True` | `None` | writes `obsm["X_<model>"]` |

   Every `tl` function that returns per-sample or per-pair values supports both modes with
   identical semantics; `tl.permanova` and `tl.mmvec` are the exceptions (see Consequences).
  diff --git a/.knowledge/modules/ml.md b/.knowledge/modules/ml.md
  @@ -1,22 +1,21 @@
   ---
   type: Module
   title: ml
  -description: scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold, taking arrays, sparse matrices and DataFrames; and to_torch, a PyTorch dataset over an AnnData's rows behind the extra torch.
  +description: scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold, taking arrays, sparse matrices and DataFrames; to_torch, a PyTorch dataset over an AnnData's rows behind the extra torch; and embed, one embedding per sample from a model a plugin registers.
   resource: /src/biotapy/ml/
   paths: ["src/biotapy/ml/**"]
  -tags: [ml, scikit-learn, torch]
  +tags: [ml, scikit-learn, torch, plugins]
   status: stable
  -generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T12:00:00Z }
  -commit: 7a9c07a
  +generated: { by: claude-code/<model>, at: <UTC commit time> }
  +commit: <parent short sha>
   ---

   # Responsibility

   Owns `bt.ml.*`, the top layer's machine-learning entry points. Today that is
  -two scikit-learn transformers (`_transformers.py`) and `to_torch`
  -(`_torch.py`, extra `torch`); `ml.embed` (slice 4C) is planned in
  -[phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md) and does not exist
  -yet. Owns no reader and no table-level transform: `pp.filter_features` and
  +two scikit-learn transformers (`_transformers.py`), `to_torch`
  +(`_torch.py`, extra `torch`) and `embed` (`_embed.py`), which runs a model
  +that a plugin registers ([embedding-plugins](/decisions/embedding-plugins.md)). Owns no reader and no table-level transform: `pp.filter_features` and
   `pp.clr` stay `pp`'s, the transformers are their fold-safe forms.

   # Entry points
  @@ -28,6 +27,10 @@ yet. Owns no reader and no table-level transform: `pp.filter_features` and
     so `fit` only validates.
   - `_torch.py:to_torch` - a map-style `torch.utils.data.Dataset` over `X` or a
     layer, one float32 row per item, paired with an `obs` label when asked.
  +- `_embed.py:embed` - loads the callable registered under `model` in the
  +  entry-point group `biotapy.embeddings` (`_embed.py:GROUP`), calls it with
  +  the AnnData and returns its samples x dimensions array, or writes
  +  `obsm["X_<model>"]` with `inplace=True`.

   # Invariants

  @@ -63,6 +66,18 @@ yet. Owns no reader and no table-level transform: `pp.filter_features` and
     -> float32; a missing label raises. Codes are per dataset (anndata drops a category a subset lacks), so build one dataset and split it with `torch.utils.data.Subset`; per-split datasets recode a class a split lacks. `_torch.py:to_torch`, `_torch.py:_labels`.
   - `to_torch` validates its arguments before it imports torch, so its error
     tests run in every CI job, not only in `ml-extras`. `_torch.py:to_torch`.
  +- `embed` reads the entry points on every call and loads only the one asked
  +  for; a name no package registers raises `KeyError` listing those installed,
  +  a name two packages register raises `ValueError` naming both.
  +  `_embed.py:embed`.
  +- `embed` trusts no plugin: the result must be a plain `numpy.ndarray` (not
  +  a masked array or a matrix), 2-D with one row per sample and at least one
  +  column, float and finite, or it raises naming the plugin, before anything
  +  is written; a result that shares memory with `X`, a layer, an `obsm`, `varm`,
  +  `obsp` or `varp` entry or a top-level `uns` array (not one nested deeper) is
  +  copied. A plugin's own exception keeps its type and gains a note naming the
  +  plugin. The model name must be letters, digits, `_`, `-` or `.`.
  +  `_embed.py:embed`, `_embed.py:_checked`.
   - Inherited scikit-learn methods (`transform`, `fit_transform`,
     `get_support`, `get_feature_names_out`, `set_output`) are named in each
     class's `Notes`, because the class template leaves inherited members off
  @@ -83,7 +98,7 @@ yet. Owns no reader and no table-level transform: `pp.filter_features` and

   `uv run --group test pytest tests/ml` (the scikit-learn estimator checks, the
   `pp` parity tests, a `Pipeline` cross-validation test, `to_torch`'s argument
  -errors). `uv run --group test --extra torch pytest -m torch` runs the
  +errors, `embed` with fake plugins). `uv run --group test --extra torch pytest -m torch` runs the
   `to_torch` tests and its docstring example, as CI's `ml-extras` job does (30-minute timeout,
   `.github/workflows/test.yaml`). The pseudocount
   warning's text is unit-tested in `tests/core/test_composition.py`.
  @@ -130,6 +145,15 @@ warning's text is unit-tested in `tests/core/test_composition.py`.
     sys.modules` after `import biotapy` tells nothing there; the
     `import-without-extras` check is meaningful only on a torch-free
     environment, as in CI.
  +- `embed`'s tests install fake plugins as a real distribution: a
  +  `<name>-1.0.dist-info` with `METADATA` and `entry_points.txt` under
  +  `tmp_path`, prepended to `sys.path`, and a module object in `sys.modules`
  +  holding the callables. Discovery is then the same code path as for an
  +  installed package, and no test imports `ml/_embed.py`.
  +  `tests/ml/test_embed.py:_installer`.
  +- `embed`'s docstring example runs MGM, so the root `conftest.py` gives
  +  `biotapy.ml._embed`'s doctests the marker `mgm`, deselected by default like
  +  `torch`. `conftest.py:_EXTRA_DOCTESTS`.
   - anndata 0.13 lists `X` as `layers[None]`, so `list(adata.layers)` holds
     `None` even when no layer was added; `to_torch`'s missing-layer error does
     not list the layers. `_torch.py:_table`.
  ```
  Copy the new `ml` description into `.knowledge/modules/index.md`'s `ml`
  line. Tick 4.4 here; log lines, first in the slice 4C section:
  ```markdown
  - **Creation**: [embedding-plugins](decisions/embedding-plugins.md) (`draft`): the entry-point group `biotapy.embeddings`, the callable `embed(adata)`, the output checks, `obsm["X_<model>"]`, weights never bundled; listed in the [decisions index](decisions/index.md).
  - **Update**: [ml](modules/ml.md) gains `embed` (entry point, discovery per call, the output checks, the fake-plugin distribution in tests, the doctest marker `mgm`), with the description copied into the [modules index](modules/index.md); [data-model-slots](contracts/data-model-slots.md): `X_<plugin>` becomes `X_<model>` and `ml.embed` joins the embedding-adding row; [pure-by-default](decisions/pure-by-default.md) gains the `ml.embed` rows; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.4.
  ```
- [x] **Step 7: Gate and commit**
  ```bash
  git add src/biotapy/ml/_embed.py src/biotapy/ml/__init__.py tests/ml/test_embed.py conftest.py pyproject.toml \
    docs/guide/machine_learning.md docs/api.md docs/contributing.md \
    .knowledge/decisions/embedding-plugins.md .knowledge/decisions/index.md .knowledge/contracts/data-model-slots.md \
    .knowledge/decisions/pure-by-default.md .knowledge/modules/ml.md .knowledge/modules/index.md \
    .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  # the slice gate; its -m mgm line selects only embed's doctest, which fails until 4.4b registers mgm (decision 34)
  git commit -m "feat(ml): add embed, sample embeddings from models that plugins register

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1480 passed, 2 skipped, 80 deselected`; `36
  passed, 1526 deselected`; `build succeeded`; `25 passed, 1537
  deselected`; `-m mgm`: `1 failed, 1561 deselected`
  (`src/biotapy/ml/_embed.py::biotapy.ml._embed.embed`, `KeyError:
  "model='mgm' is not an installed embedding plugin; installed: []"`).


### Task 4.4b: MGM, the reference plugin

**Files:** create `src/biotapy/ml/_mgm.py`, `tests/ml/test_mgm.py`,
`tests/mgm/export_reference.py`, `tests/data/mgm/counts.csv`,
`tests/data/mgm/embeddings.csv`, `tests/data/mgm/NOTICE.txt`; modify
`src/biotapy/ml/_embed.py` (docstring), `pyproject.toml`,
`.github/workflows/test.yaml`, `tests/test_ci.py`,
`docs/guide/machine_learning.md`, `docs/contributing.md`,
`.knowledge/modules/ml.md`, `.knowledge/contracts/r-golden-parity.md`,
`.knowledge/decisions/optional-heavy-dependencies.md`,
`.knowledge/modules/core.md`, `.knowledge/roadmap/phase-4-ml-multiomics.md`,
`.knowledge/log.md`.
**Not touched:** `ml/__init__.py` (the plugin is reached through its entry
point, never `bt.ml`); `docs/api.md` (no public function: MGM is documented
in `ml.embed`'s docstring and the guide); `conftest.py` (4.4's hook already
marks `embed`'s doctest); the `network` job (no `mgm` test is marked
`network`, decision 31); `.knowledge/decisions/embedding-plugins.md` (it
already names MGM's entry point); the end-to-end docs page (4.6b, slice 4D).
**Interfaces:**
- Consumes: `bt.ml.embed` and the group `biotapy.embeddings` (4.4);
  `_core.make_pooch` (4.C0), `_core.as_csr`, `divide_rows`,
  `finite_non_negative`, `import_optional`, `sum_by`, `warn_user`; the extra
  `mgm` (4.C1); the marker `mgm` (4.4); `bt.datasets.global_patterns`,
  `bt.pp.tax_glom`, `bt.pp.relative`.
- Produces: the entry point `mgm = "biotapy.ml._mgm:embed"`, a callable
  `embed(adata: AnnData) -> npt.NDArray[np.float32]` (samples x 256);
  `tests/data/mgm/` and the script that writes it; `ml-extras` runs
  `uv run --group test --extra mgm pytest -m mgm` with a pooch cache.

- [x] **Step 1: MGM's own embeddings.** Create
  `tests/mgm/export_reference.py`:
  ```python
  """Write tests/data/mgm/: a small genus table and MGM 0.5.8's own embedding of each of its samples.

  microformer-mgm pins numpy 1.24, pandas 2.0, torch 2.0 and transformers 4.33, so it cannot be installed
  beside biotapy. Run from the repository root, never in CI:

      HF_HOME=<scratch directory> \
      uv run --no-project --python 3.11 --with microformer-mgm==0.5.8 --with torch==2.0.1+cpu \
          --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match \
          python tests/mgm/export_reference.py

  The table is synthetic; its genus names are MGM's vocabulary, plus one name MGM lacks and one feature with no
  genus. MGM's own code builds the tokens (`MicroCorpus`) and runs the model (`GPT2LMHeadModel`, transformers 4.33),
  one sample at a time as MGM's notebook does. The embedding is the mean of the last hidden layer over the sample's
  tokens (<bos>, its genera, <eos>), the "element-wise mean pooling" MGM's paper uses for the pretrained model
  (Methods 4.5). MGM drops a sample with no count in its vocabulary; such a sample is embedded here from the tokens
  <bos> <eos>, as biotapy does. "Blautia" appears twice (f5 and f16): summing the two moves its rank in s2 from 5 to 2,
  so a reader that does not sum features of one genus fails the parity test.
  """

  import re
  from pathlib import Path

  import numpy as np
  import pandas as pd

  # MGM, torch and transformers are imported inside the functions: pytest imports this file for doctests in every run.
  OUT = Path("tests/data/mgm")
  GUT = ["Bacteroides", "Prevotella", "Faecalibacterium", "Bifidobacterium", "Akkermansia", "Blautia", "Roseburia",
         "Alistipes", "Streptococcus", "Lactobacillus", "Veillonella", "Ruminococcus"]  # fmt: skip


  def table() -> pd.DataFrame:
      """Features x (genus, s1..s7), counts."""
      from mgm.CLI.CLI_utils import find_pkg_resource

      vocabulary = pd.read_csv(find_pkg_resource("resources/phylogeny.csv"), index_col=0).index
      plain = [name[3:] for name in vocabulary if re.fullmatch(r"g__[A-Za-z0-9_]+", name) and name[3:] not in GUT]
      many = plain[::16][:600]
      genera = [*GUT, "Escherichia", "Escherichia-Shigella", "Notagenus", None, "Blautia", *many]
      counts = np.zeros((len(genera), 7), dtype=np.int64)
      gut = np.arange(1, len(GUT) + 1)
      counts[: len(GUT), 0] = gut * 10  # s1: twelve genera
      counts[: len(GUT), 1] = gut[::-1] * 7  # s2: the same, other counts, and the rest below
      # Escherichia twice (one token), a genus MGM lacks, no genus, Blautia again
      counts[12:17, 1] = [3, 5, 40, 25, 40]
      counts[2, 2] = 9  # s3: one genus
      counts[14:16, 3] = [8, 2]  # s4: only a genus MGM lacks and no genus; s5: all zero
      counts[17:, 5] = (np.arange(600) * 7919) % 600 + 1  # s6: 600 genera, more than the 510 tokens a sample can hold
      counts[[0, 3, 20, 400], 6] = [500, 1, 1, 30]  # s7: few genera, two of them once
      frame = pd.DataFrame(counts, columns=[f"s{i}" for i in range(1, 8)], index=[f"f{i}" for i in range(len(genera))])
      frame.insert(0, "genus", genera)
      frame.index.name = "feature"
      return frame


  def embed(input_ids, attention_mask, model) -> np.ndarray:
      """MGM's notebook `cal_embed` on CPU, with the mean over the sample's tokens instead of its last token."""
      model.eval()
      hidden = model(input_ids=input_ids, attention_mask=attention_mask, output_hidden_states=True).hidden_states[-1]
      return hidden.squeeze(0)[attention_mask.squeeze(0) == 1].mean(0).detach().numpy()


  def main() -> None:
      import torch
      from mgm.CLI.CLI_utils import find_pkg_resource
      from mgm.src.MicroCorpus import MicroCorpus
      from mgm.src.utils import CustomUnpickler
      from transformers import GPT2LMHeadModel

      frame = table()
      OUT.mkdir(parents=True, exist_ok=True)
      frame.to_csv(OUT / "counts.csv")
      columns = ["g__" + genus if isinstance(genus, str) else "unclassified" for genus in frame["genus"]]
      abundance = pd.DataFrame(frame.drop(columns="genus").T.to_numpy(), index=frame.columns[1:], columns=columns)
      with open(find_pkg_resource("resources/MicroTokenizer.pkl"), "rb") as file:
          tokenizer = CustomUnpickler(file).load()
      corpus = MicroCorpus(abu=abundance.astype(float), tokenizer=tokenizer, max_len=512, preprocess=True)
      model = GPT2LMHeadModel.from_pretrained(find_pkg_resource("resources/general_model"))
      kept = list(corpus.data.index)
      empty = torch.tensor([[2, 3] + [0] * 510]), torch.tensor([[1.0, 1.0] + [0.0] * 510])
      rows = {}
      with torch.no_grad():
          for sample in abundance.index:
              if sample in kept:
                  item = corpus[kept.index(sample)]
                  rows[sample] = embed(item["input_ids"][None], item["attention_mask"][None], model)
              else:
                  rows[sample] = embed(*empty, model)
      embeddings = pd.DataFrame(rows).T
      embeddings.columns = [f"e{i}" for i in range(embeddings.shape[1])]
      embeddings.index.name = "sample"
      embeddings.to_csv(OUT / "embeddings.csv", float_format="%.9g")


  if __name__ == "__main__":
      main()
  ```
  Run it from the repository root with the command in its docstring (uv
  downloads CPython 3.11 if none is installed; nothing touches `.venv`).
  Expected: MGM prints `2 samples are dropped for all zeroes` and `Total 5
  samples. Max length is 602. Average length is 128.0. Min length is 3.`;
  `md5sum tests/data/mgm/*.csv` on the prototype (x86-64) ->
  `bdacda5f68ba7184bbb16f99c22c6168  tests/data/mgm/counts.csv`,
  `4bdd9fe18a9b72396c5d4ddd3a70f477  tests/data/mgm/embeddings.csv`; a second
  run gives the same bytes. `counts.csv` is arithmetic and must match
  exactly; `embeddings.csv` may differ in its last digits on another CPU,
  which the test's tolerance absorbs (the prototype's files are also in the
  session scratchpad, `p4c/tests/data/mgm/`). Create
  `tests/data/mgm/NOTICE.txt`:
  ```text
  counts.csv is synthetic, written for biotapy (BSD-3-Clause): 7 samples x 617
  features. 614 of their genus names come from MGM's vocabulary
  (mgm/resources/phylogeny.csv); "Escherichia-Shigella" reads as MGM's
  Escherichia, "Notagenus" is not in it, one feature has no genus, and "Blautia"
  is the genus of two features.

  embeddings.csv holds MGM's own embedding of each sample of counts.csv (7 x 256
  float32 values), computed by tests/mgm/export_reference.py with
  microformer-mgm 0.5.8 (its wheel's SHA-256
  210891685565022ea869e88a7452769dc6fbe3f699d2a133e15338d1e30eb92b), torch
  2.0.1+cpu, transformers 4.33.3, numpy 1.24.3 and pandas 2.0.3 on CPython
  3.11.16: MGM's MicroCorpus builds the tokens and GPT2LMHeadModel runs its
  pretrained general model, one sample at a time; the embedding is the mean of
  the last hidden layer over the sample's tokens. Two runs gave byte-identical
  files. They are numbers derived from MGM's MIT-licensed model:

  MIT License

  Copyright (c) 2024 NingLab

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
  ```
- [x] **Step 2: Failing tests.** Create `tests/ml/test_mgm.py`:
  ```python
  import sys
  from importlib.metadata import entry_points
  from pathlib import Path

  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp
  from anndata import AnnData
  from sklearn.neighbors import NearestNeighbors

  import biotapy as bt

  DATA = Path(__file__).parents[1] / "data" / "mgm"
  # torch 2.13-2.14's first tanh in a process can saturate on CPUs running more than 4 threads (measured: the
  # embedding's first call off by up to 1.6e-4, later calls by 1.7e-6), so equality is checked to 1e-3.
  ATOL = 1e-3
  # Only a process's first call is affected by that race, so a later call matches MGM's own output to 1e-5 (measured 1.7e-6);
  # a changed last token (an <eos> kept at the cut) moves a row by 7.7e-4, which 1e-3 hides.
  WARM_ATOL = 1e-5


  def _reference_table():
      """tests/data/mgm/counts.csv as an AnnData: 7 samples x 617 features, var["genus"] (one missing)."""
      counts = pd.read_csv(DATA / "counts.csv", index_col=0)
      return AnnData(
          X=sp.csr_matrix(counts.drop(columns="genus").T.to_numpy()),
          obs=pd.DataFrame(index=counts.columns[1:]),
          var=counts[["genus"]],
      )


  def test_mgm_is_registered_as_a_plugin():
      points = [point for point in entry_points(group="biotapy.embeddings") if point.name == "mgm"]
      assert [point.value for point in points] == ["biotapy.ml._mgm:embed"]


  def test_without_a_genus_column_raises(make_adata):
      with pytest.raises(KeyError, match=r"mgm reads var\['genus'\], which this table does not have"):
          bt.ml.embed(make_adata(np.array([[1, 2]])), "mgm")


  def test_negative_values_raise(make_adata):
      adata = make_adata(np.array([[1.0, -2.0]]))
      adata.var["genus"] = ["Bacteroides", "Prevotella"]
      with pytest.raises(ValueError, match="mgm reads counts or relative abundances; adata.X holds negative"):
          bt.ml.embed(adata, "mgm")


  def test_without_the_extra_names_it(monkeypatch):
      # None in sys.modules makes the import fail whether or not the extra is installed.
      monkeypatch.setitem(sys.modules, "torch", None)
      monkeypatch.setitem(sys.modules, "transformers", None)
      with pytest.raises(ImportError, match=r"pip install 'biotapy\[mgm\]'"):
          bt.ml.embed(bt.pp.tax_glom(bt.datasets.toy(), "genus"), "mgm")


  @pytest.mark.mgm
  def test_matches_mgm_s_own_embedding():
      expected = pd.read_csv(DATA / "embeddings.csv", index_col=0)
      with pytest.warns(UserWarning) as record:
          result = bt.ml.embed(_reference_table(), "mgm")
      assert [str(warning.message) for warning in record] == [
          "mgm leaves out 2 of 617 features: 1 without a genus and 1 whose genus is not one of MGM's (Notagenus)",
          "mgm embeds 2 sample(s) from <bos> <eos> alone, none of their genera being MGM's: s4, s5",
      ]
      assert result.shape == (7, 256) and result.dtype == np.float32
      np.testing.assert_allclose(result, expected.to_numpy(), rtol=0, atol=ATOL)
      with pytest.warns(UserWarning):
          warm = bt.ml.embed(_reference_table(), "mgm")
      np.testing.assert_allclose(warm, expected.to_numpy(), rtol=0, atol=WARM_ATOL)


  @pytest.mark.mgm
  def test_features_of_one_genus_are_summed_so_tax_glom_changes_nothing():
      tdata = bt.datasets.toy()  # f4 and f5 are both Bacteroides; f8 has no genus
      with pytest.warns(UserWarning, match="1 without a genus"):
          raw = bt.ml.embed(tdata, "mgm")
      np.testing.assert_allclose(raw, bt.ml.embed(bt.pp.tax_glom(tdata, "genus"), "mgm"), rtol=0, atol=ATOL)


  @pytest.mark.mgm
  def test_the_leave_out_warning_lists_unknown_genera_only_when_there_are_some():
      with pytest.warns(UserWarning) as record:
          bt.ml.embed(bt.datasets.toy(), "mgm")
      assert [str(warning.message) for warning in record] == [
          "mgm leaves out 1 of 8 features: 1 without a genus and 0 whose genus is not one of MGM's"
      ]


  @pytest.mark.mgm
  def test_relative_abundances_embed_as_their_counts():
      genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      relative = genera.copy()
      relative.X = sp.csr_matrix(bt.pp.relative(genera).layers["relative"])
      np.testing.assert_allclose(bt.ml.embed(relative, "mgm"), bt.ml.embed(genera, "mgm"), rtol=0, atol=ATOL)


  @pytest.mark.mgm
  def test_feature_order_does_not_matter():
      genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      shuffled = genera[:, ::-1].copy()
      np.testing.assert_allclose(bt.ml.embed(shuffled, "mgm"), bt.ml.embed(genera, "mgm"), rtol=0, atol=ATOL)


  @pytest.mark.mgm
  def test_an_all_zero_feature_changes_nothing():
      genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      padded = AnnData(
          X=sp.hstack([genera.X, sp.csr_matrix((genera.n_obs, 1), dtype=genera.X.dtype)], format="csr"),
          obs=genera.obs,
          var=pd.DataFrame({"genus": [*genera.var["genus"], "Akkermansia"]}, index=[*genera.var_names, "zero"]),
      )
      np.testing.assert_allclose(bt.ml.embed(padded, "mgm"), bt.ml.embed(genera, "mgm"), rtol=0, atol=ATOL)


  @pytest.mark.mgm
  def test_single_sample():
      genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      alone = bt.ml.embed(genera[[3]].copy(), "mgm")
      np.testing.assert_allclose(alone, bt.ml.embed(genera, "mgm")[[3]], rtol=0, atol=ATOL)


  @pytest.mark.mgm
  def test_keeps_the_input(assert_unchanged):
      genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      before = genera.copy()
      bt.ml.embed(genera, "mgm")
      assert_unchanged(before, genera)


  @pytest.mark.mgm
  def test_embeds_global_patterns_end_to_end():
      # Phase 4 exit gate 2: one foundation model plugged in end to end (decisions 10, 17).
      tdata = bt.pp.tax_glom(bt.datasets.global_patterns(), "genus")
      with pytest.warns(UserWarning, match="mgm leaves out 96 of 996 features: 0 without a genus and 96 whose") as record:
          bt.ml.embed(tdata, "mgm", inplace=True)
      assert not [warning for warning in record if "<bos> <eos>" in str(warning.message)]  # no sample is empty
      embedding = tdata.obsm["X_mgm"]
      assert embedding.shape == (26, 256) and embedding.dtype == np.float32 and np.isfinite(embedding).all()
      with pytest.warns(UserWarning, match="mgm leaves out 96"):
          np.testing.assert_allclose(bt.ml.embed(tdata, "mgm"), embedding, rtol=0, atol=ATOL)
      # Every sample's nearest neighbour in the embedding comes from the same environment.
      neighbour = NearestNeighbors(n_neighbors=2, metric="cosine").fit(embedding).kneighbors(embedding)[1][:, 1]
      types = tdata.obs["SampleType"].to_numpy()
      assert (types[neighbour] == types).all()


  def _cached_file(name):
      """The file ``name`` that MGM's wheel was extracted to, in biotapy's data cache."""
      from biotapy.ml import _mgm

      cache = _mgm.make_pooch(_mgm._BASE_URL, {_mgm._WHEEL: _mgm._SHA256})
      paths = cache.fetch(_mgm._WHEEL, processor=_mgm.pooch.Unzip(members=_mgm._MEMBERS))
      return next(Path(path) for path in paths if Path(path).name == name)


  @pytest.mark.mgm
  @pytest.mark.parametrize("name", ["phylogeny.csv", "config.json"])
  def test_a_damaged_extracted_file_is_extracted_again(name):
      genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      expected = bt.ml.embed(genera, "mgm")
      path = _cached_file(name)
      original = path.read_bytes()
      path.write_bytes(original[: len(original) // 2])
      result = bt.ml.embed(genera, "mgm")
      assert path.read_bytes() == original
      np.testing.assert_allclose(result, expected, rtol=0, atol=ATOL)


  @pytest.mark.mgm
  def test_a_file_that_stays_damaged_raises(monkeypatch):
      from biotapy.ml import _mgm

      monkeypatch.setitem(_mgm._SHA256_OF, "phylogeny.csv", "0" * 64)
      genera = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      with pytest.raises(
          ValueError, match=r"after extracting it again; delete .*microformer_mgm-0\.5\.8-py3-none-any\.whl\.unzip\n"
      ):
          bt.ml.embed(genera, "mgm")
  ```
  In `tests/test_ci.py`, before `test_ml_extras_job_has_a_timeout`:
  ```diff
  diff --git a/tests/test_ci.py b/tests/test_ci.py
  @@ -78,6 +78,14 @@ def test_the_sdist_ships_the_root_conftest_that_marks_the_torch_doctests():
       assert "/conftest.py" in sdist["include"]


  +def test_ml_extras_job_runs_the_mgm_marker_with_a_pooch_cache():
  +    steps = WORKFLOW["jobs"]["ml-extras"]["steps"]
  +    run = [step for step in steps if step.get("run", "").strip() == "uv run --group test --extra mgm pytest -m mgm"]
  +    assert run and run[0]["env"]["BIOTAPY_DATA_DIR"] == "${{ github.workspace }}/.pooch"
  +    cache = next(step for step in steps if step.get("uses", "").startswith("actions/cache@"))["with"]
  +    assert cache["path"] == "${{ github.workspace }}/.pooch" and "src/biotapy/ml/_mgm.py" in cache["key"]
  +
  +
   def test_ml_extras_job_has_a_timeout():
       assert WORKFLOW["jobs"]["ml-extras"]["timeout-minutes"] == 30

  ```
  R11.2 coverage: happy path (MGM's own embeddings, GlobalPatterns end to
  end); edge cases (an all-zero sample and one with only unknown genera,
  `s5` and `s4`; a NaN genus, `f15`; a sample past 510 genera, `s6`; a
  single sample; an all-zero feature; a truncated extracted file, extracted again; a warm call at 1e-5); purity (`test_keeps_the_input`); invariants as tests rather
  than Hypothesis properties, because each MGM call runs a transformer
  (feature order, scale, summing per genus); a parity test against MGM's own
  code instead of an R golden (no R equivalent). The four tests without the
  marker need neither the extra nor the download.
- [x] **Step 3: Run, expect failure** - `uv sync --all-groups && uv run
  --group test pytest tests/ml/test_mgm.py tests/test_ci.py -q` -> `5
  failed, 22 passed, 12 deselected` (no entry point yet: the input tests get
  `KeyError: "model='mgm' is not an installed embedding plugin; installed:
  []"`; `test_ci` has no `-m mgm` step). `uv run --group test --extra mgm
  pytest tests/ml/test_mgm.py -q -m mgm` -> `12 failed, 4 deselected` (the
  same `KeyError`, or `DID NOT WARN` around it).
- [x] **Step 4: Implement.** Create `src/biotapy/ml/_mgm.py`:
  ```python
  """MGM, the reference embedding plugin: entry point ``mgm`` in ``biotapy.embeddings``, extra ``mgm``."""

  import hashlib
  import logging
  from collections.abc import Iterable
  from pathlib import Path
  from typing import Any, cast

  import numpy as np
  import numpy.typing as npt
  import pandas as pd
  import pooch
  from anndata import AnnData

  from biotapy._core import as_csr, divide_rows, finite_non_negative, import_optional, make_pooch, sum_by, warn_user

  # microformer-mgm 0.5.8's wheel (MIT) carries the pretrained general model; files.pythonhosted.org never changes a file.
  _WHEEL = "microformer_mgm-0.5.8-py3-none-any.whl"
  _BASE_URL = (
      "https://files.pythonhosted.org/packages/4b/9c/829a1e59d5e618756ce8e57553e15d13bb16c58d37896f03233d8697515a/"
  )
  _SHA256 = "sha256:210891685565022ea869e88a7452769dc6fbe3f699d2a133e15338d1e30eb92b"
  _MEMBERS = [
      "mgm/resources/general_model/config.json",
      "mgm/resources/general_model/pytorch_model.bin",
      "mgm/resources/phylogeny.csv",
  ]
  # SHA-256 of each extracted file: Unzip re-extracts only a missing member, so a truncated one would stay.
  _SHA256_OF = {
      "config.json": "88e4940036404ca90d2b1a79f38948efd419454484eb394196edeb9957fec499",
      "pytorch_model.bin": "7a83b442dcb0cbe2a4248ea9922de58e15edfa3a1ad71d07519dc8594f4511b8",
      "phylogeny.csv": "016b8264eed155b134ab55e387e8f21c7e46241bf987a83df49d6604c6b378db",
  }
  _LOG = logging.getLogger(__name__)
  # How MGM reads a genus from a column name (mgm/src/MicroCorpus.py, MicroCorpus._preprocess).
  _TOKEN = r"(g__[A-Za-z0-9_]+)"
  # MGM's vocabulary: <pad>, <mask>, <bos>, <eos>, then phylogeny.csv's genera in file order (mgm/resources/MicroTokenizer.pkl).
  _BOS, _EOS, _FIRST_GENUS = 2, 3, 4
  _MAX_TOKENS = 512


  def embed(adata: AnnData) -> npt.NDArray[np.float32]:
      """MGM's embedding of each sample, the mean of its last hidden layer over the sample's tokens (see ``ml.embed``)."""
      # Checked before torch, transformers or the download, so these errors need no extra (and are tested in every job).
      if "genus" not in adata.var.columns:
          msg = "mgm reads var['genus'], which this table does not have"
          raise KeyError(msg)
      if not finite_non_negative(as_csr(adata.X)):
          msg = "mgm reads counts or relative abundances; adata.X holds negative or non-finite values"
          raise ValueError(msg)
      torch: Any = import_optional("torch", extra="mgm")
      transformers: Any = import_optional("transformers", extra="mgm")
      files = _extracted_files()
      # anndata types var columns as Series | DataArray (its lazy variant); the data model guarantees a Series.
      genus = cast("pd.Series[str]", adata.var["genus"])
      sentences = _sentences(adata, genus, pd.read_csv(files["phylogeny.csv"], index_col=0))
      model = transformers.GPT2Model(transformers.GPT2Config.from_json_file(files["config.json"]))
      weights = torch.load(files["pytorch_model.bin"], map_location="cpu", weights_only=True)
      # The file holds a GPT2LMHeadModel: the body under "transformer.", and the head, which an embedding does not use.
      model.load_state_dict(
          {key.removeprefix("transformer."): value for key, value in weights.items() if key != "lm_head.weight"}
      )
      model.eval()
      out = np.empty((len(sentences), model.config.n_embd), dtype=np.float32)
      # One sample at a time, unpadded: on the CPU, batches were slower and grew memory (phase-4 plan, task 4.4b).
      with torch.inference_mode():
          for row, sentence in enumerate(sentences):
              out[row] = model(input_ids=torch.from_numpy(sentence)[None]).last_hidden_state[0].mean(0).numpy()
      return out


  def _extracted_files() -> dict[str, Path]:
      """MGM's three files from the cached wheel, extracted again once if one does not match its SHA-256."""
      cache = make_pooch(_BASE_URL, {_WHEEL: _SHA256})

      def fetch() -> tuple[dict[str, Path], list[str]]:
          paths = cache.fetch(_WHEEL, processor=pooch.Unzip(members=_MEMBERS))
          files = {Path(path).name: Path(path) for path in paths}  # Unzip lists them in os.walk order
          digest = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in files.items()}
          return files, [name for name, value in digest.items() if value != _SHA256_OF[name]]

      files, damaged = fetch()
      if damaged:
          _LOG.warning("extracting %s again: it does not match its SHA-256", ", ".join(damaged))
          for name in damaged:
              files[name].unlink()
          files, damaged = fetch()
      if damaged:
          msg = f"{', '.join(damaged)} does not match its SHA-256 after extracting it again; delete {Path(cache.abspath) / (_WHEEL + '.unzip')}"
          raise ValueError(msg)
      return files


  def _sentences(adata: AnnData, genus: "pd.Series[str]", phylogeny: pd.DataFrame) -> list[npt.NDArray[np.int64]]:
      """Each sample's tokens as MGM builds them: <bos>, its genera by standardised abundance, <eos>, cut to 512."""
      tokens = ("g__" + genus.astype("string")).str.extract(_TOKEN, expand=False)
      codes = phylogeny.index.get_indexer(pd.Index(tokens))
      missing = genus.isna().to_numpy()
      unknown = (codes < 0) & ~missing
      if missing.any() or unknown.any():
          msg = (
              f"mgm leaves out {int(missing.sum() + unknown.sum())} of {codes.size} features: {int(missing.sum())} "
              f"without a genus and {int(unknown.sum())} whose genus is not one of MGM's"
          )
          if unknown.any():
              msg += f" ({_listed(genus[unknown])})"
          warn_user(msg)
      # Features of one genus are summed, and relative abundance is taken over MGM's genera only, as MGM does.
      grouped = sum_by(as_csr(adata.X), codes, len(phylogeny))
      grouped.sort_indices()
      relative = divide_rows(grouped, np.asarray(grouped.sum(axis=1), dtype=np.float64).ravel())
      mean, std = phylogeny["mean"].to_numpy(), phylogeny["std"].to_numpy()
      sentences = []
      for row in range(relative.shape[0]):
          span = slice(relative.indptr[row], relative.indptr[row + 1])
          genera = relative.indices[span]
          standardised = (relative.data[span] - mean[genera]) / std[genera]
          # MGM keeps a genus whose standardised value is above that of zero abundance, and sorts with pandas.
          kept = pd.Series(standardised, index=genera)[standardised > (0 - mean[genera]) / std[genera]]
          order = kept.sort_values(ascending=False).index.to_numpy()
          sentences.append(np.r_[_BOS, order + _FIRST_GENUS, _EOS][:_MAX_TOKENS].astype(np.int64))
      empty = [name for name, sentence in zip(adata.obs_names, sentences, strict=True) if sentence.size == 2]
      if empty:
          warn_user(
              f"mgm embeds {len(empty)} sample(s) from <bos> <eos> alone, none of their genera being MGM's: {_listed(empty)}"
          )
      return sentences


  def _listed(names: "Iterable[str]") -> str:
      """Up to five distinct names, sorted, for a warning."""
      distinct = sorted({str(name) for name in names})
      return ", ".join(distinct[:5]) + (", ..." if len(distinct) > 5 else "")
  ```
  `_sentences` is the one single-use helper (R4.4): `embed` with it inline
  would hold over 40 statements, past R5's 30; `_listed` serves both
  warnings. R2.1: transformers runs the model (`GPT2Model`), torch loads the
  weights, pooch downloads and unzips, `_core.sum_by` and `divide_rows`
  group and scale, pandas sorts exactly as MGM does; no library builds MGM's
  sentence. Then:
  ```diff
  diff --git a/.github/workflows/test.yaml b/.github/workflows/test.yaml
  @@ -212,8 +212,10 @@ jobs:
             RPY2_CFFI_MODE: API
           run: uv run --group test --extra r pytest -m r

  -  # Runs the tests that need PyTorch (marker torch, extra torch), which every other job deselects; uv installs torch's
  -  # CPU wheel from PyTorch's index ([tool.uv] in pyproject.toml, decisions/optional-heavy-dependencies).
  +  # Runs the tests that need PyTorch (marker torch, extra torch) and MGM (marker mgm, extra mgm), which every other job
  +  # deselects; uv installs torch's CPU wheel from PyTorch's index ([tool.uv] in pyproject.toml,
  +  # decisions/optional-heavy-dependencies). The MGM tests download its weights and GlobalPatterns into a pooch cache
  +  # of their own, keyed on both loaders.
     ml-extras:
       runs-on: ubuntu-latest
       timeout-minutes: 30
  @@ -232,6 +234,17 @@ jobs:
           run: uv run --group test --extra torch python -c "import torch; print('torch', torch.__version__)"
         - name: Run the PyTorch tests
           run: uv run --group test --extra torch pytest -m torch
  +      - name: Cache pooch downloads
  +        uses: actions/cache@55cc8345863c7cc4c66a329aec7e433d2d1c52a9 # v6.1.0
  +        with:
  +          path: ${{ github.workspace }}/.pooch
  +          key: pooch-mgm-${{ hashFiles('src/biotapy/datasets/_remote.py', 'src/biotapy/ml/_mgm.py') }}
  +      - name: Log the transformers version
  +        run: uv run --group test --extra mgm python -c "import transformers; print('transformers', transformers.__version__)"
  +      - name: Run the MGM tests
  +        env:
  +          BIOTAPY_DATA_DIR: ${{ github.workspace }}/.pooch
  +        run: uv run --group test --extra mgm pytest -m mgm

     # Builds the docs as Read the Docs does, executing every notebook (Phase 1 exit gate): the
     # phyloseq vignette downloads GlobalPatterns, enterotype and esophagus through the network
  diff --git a/pyproject.toml b/pyproject.toml
  @@ -54,6 +54,8 @@ optional-dependencies.torch = [ "torch>=2.9" ]
   urls.Documentation = "https://biotapy.readthedocs.io/"
   urls.Homepage = "https://github.com/pedrocr83/biotapy"
   urls.Source = "https://github.com/pedrocr83/biotapy"
  +# ml.embed's models (decisions/embedding-plugins); biotapy registers its own reference model as any plugin would.
  +entry-points."biotapy.embeddings".mgm = "biotapy.ml._mgm:embed"

   [dependency-groups]
   dev = [
  @@ -105,6 +107,7 @@ build.targets.sdist.exclude = [
     "/docs/generated",
     "/docs/jupyter_execute",
     "/tests/humann",
  +  "/tests/mgm",
     "/tests/r",
     "/tests/test_ci.py",
     "/tests/test_knowledge_bundle.py",
  @@ -189,11 +192,13 @@ lint.pylint.max-statements = 30
   files = [ "docs/extensions", "src/biotapy" ]
   python_version = "3.12"
   # biom-format and scikit-learn ship no type annotations at all, nor do scikit-bio's
  -# subsample_counts, tree_basis and ancombc2, threadpoolctl's threadpool_limits and mudata's MuData.update; exempt only
  -# calls into them (not our own code) from strict's disallow_untyped_calls (mypy/checkexpr.py matches by callee fullname).
  +# subsample_counts, tree_basis and ancombc2, threadpoolctl's threadpool_limits, mudata's MuData.update and pooch's
  +# Unzip; exempt only calls into them (not our own code) from strict's disallow_untyped_calls (mypy/checkexpr.py
  +# matches by callee fullname).
   untyped_calls_exclude = [
     "biom",
     "mudata",
  +  "pooch.processors.Unzip",
     "skbio.stats._subsample",
     "skbio.stats.composition._ancombc.ancombc2",
     "skbio.stats.composition._base.tree_basis",
  diff --git a/src/biotapy/ml/_embed.py b/src/biotapy/ml/_embed.py
  @@ -64,6 +64,36 @@ def embed(adata: AnnData, model: str, *, inplace: bool = False) -> npt.NDArray[n

       biotapy finds it when ``embed`` is called, so installing the package is
       enough. It checks the array before returning or storing it.
  +
  +    ``"mgm"`` is MGM, the Microbial General Model (MIT licence), a GPT-2
  +    pretrained on genus profiles from MGnify. It needs the extra ``mgm``
  +    (``pip install 'biotapy[mgm]'``) and, at the first call, downloads the
  +    pretrained model (33 MB, microformer-mgm 0.5.8's wheel from PyPI, checked
  +    against its SHA-256) into biotapy's data cache: ``BIOTAPY_DATA_DIR`` if
  +    set, else pooch's per-user cache. Each feature's ``var["genus"]`` is read
  +    as MGM reads ``g__<genus>``: the name up to its first character other than
  +    a letter, digit or underscore (``Escherichia-Shigella`` is
  +    ``Escherichia``). Features whose genus is missing or outside MGM's 9,665
  +    genera are left out, with one warning that counts them, and features of
  +    the same genus are summed, so ``pp.tax_glom(tdata, "genus")`` first
  +    changes nothing. Then, as MGM's own preprocessing does, each sample becomes
  +    relative abundances, its genera are sorted by abundance standardised with
  +    MGM's per-genus mean and standard deviation, and the sentence ``<bos>``,
  +    genera, ``<eos>`` is cut to 512 tokens. A sample without a single known
  +    genus is embedded from ``<bos> <eos>``, with a warning naming it. The
  +    embedding is the mean of the model's last hidden layer over the sample's
  +    tokens (256 float32 values), the mean pooling MGM's authors use for the
  +    pretrained model; it matches MGM 0.5.8's own forward pass to 2e-6 (up to about 2e-4 on the
  +    first call of a process, on a CPU with more than four threads). The
  +    model runs on the CPU, one sample at a time: batches were slower there and
  +    needed up to 1.5 GB more memory. Cite Zhang et al. (2026) when you publish
  +    results that use it.
  +
  +    References
  +    ----------
  +    Zhang H, Zhang Y, Kang Z, Xiong J, Yang R, Ning K (2026) MGM as a
  +    large-scale pretrained foundation model for microbiome analyses in diverse
  +    contexts. Adv Sci 13:e13333.

       Examples
       --------
  ```
  pyproject-fmt moves the `entry-points` line after `urls`, as shown. The
  `ml-extras` cache key hashes `_mgm.py` (its URL and hash) and `_remote.py`
  (GlobalPatterns' pin).
- [x] **Step 5: Run, expect pass** - `uv sync --all-groups && uv run --group
  test pytest tests/ml/test_mgm.py tests/test_ci.py -q` -> `27 passed, 12
  deselected`; with `BIOTAPY_DATA_DIR` set, `uv run --group test --extra mgm
  pytest tests/ml/test_mgm.py src/biotapy/ml -q -m "mgm or not mgm" -W
  error::UserWarning` -> `20 passed` (the 16 tests, `embed`'s, `to_torch`'s
  and the two transformer doctests; the first run downloads 33 MB); `uv run
  --group dev --group doc --extra mgm mypy` -> `Success: no issues found in
  68 source files`; `uv run --no-dev python -c "import importlib, pkgutil,
  sys, biotapy; [importlib.import_module(m.name) for m in
  pkgutil.walk_packages(biotapy.__path__, 'biotapy.')]; print('torch' in
  sys.modules, 'transformers' in sys.modules)"` -> `False False`. Mutation
  check (each on `_mgm.py` alone, then `git checkout` it): relative
  abundance over `as_csr(adata.X)`'s row sums, `ascending=True`,
  `.last_hidden_state[0][-1]`, no `model.eval()`, `standardised > 0`,
  `kept.sort_index()`: `-m mgm -k "matches or global"` fails for each (the
  first fails only the parity test). `uv build --sdist` ships
  `tests/data/mgm/` and not `tests/mgm/`.
- [x] **Step 6: Docs.**
  ````diff
  diff --git a/docs/contributing.md b/docs/contributing.md
  @@ -92,6 +92,22 @@ uv run --group test --extra torch pytest -m torch

   CI runs them in the `ml-extras` job, on Linux with Python 3.13.

  +### MGM tests
  +
  +Tests that run MGM (`bt.ml.embed(..., "mgm")`, and `ml.embed`'s docstring example) carry the marker `mgm`
  +and are excluded by default too. They need the `mgm` extra and download MGM's weights (33 MB) and
  +GlobalPatterns through the pooch cache:
  +
  +```bash
  +BIOTAPY_DATA_DIR=.pooch uv run --group test --extra mgm pytest -m mgm
  +```
  +
  +`tests/data/mgm/embeddings.csv` is MGM's own embedding of `tests/data/mgm/counts.csv`. MGM's package pins
  +numpy 1.24 and torch 2.0, so it cannot share biotapy's environment; to regenerate the file, run
  +`tests/mgm/export_reference.py` as its docstring says (Python 3.11, never in CI).
  +
  +CI runs them in the `ml-extras` job, after the PyTorch tests, with a pooch cache of its own.
  +
   ### Regenerating the R golden files

   The golden CSVs under `tests/golden/` and the R-written fixtures under
  diff --git a/docs/guide/machine_learning.md b/docs/guide/machine_learning.md
  @@ -124,6 +124,50 @@ bt.ml.embed(genera, "mgm", inplace=True)
   genera.obsm["X_mgm"].shape  # (26, 256)
   ```

  +### MGM
  +
  +`"mgm"` is MGM, the Microbial General Model of Zhang et al. (2026): a GPT-2
  +with 8 layers and 256 dimensions, pretrained on genus profiles from MGnify,
  +released under the MIT licence. biotapy ships it as a plugin behind the extra
  +`mgm`, which installs torch and transformers:
  +
  +```bash
  +pip install 'biotapy[mgm]'
  +```
  +
  +The first call downloads the pretrained model once: 33 MB, the
  +`microformer-mgm` 0.5.8 wheel from PyPI, checked against its SHA-256, into the
  +cache `bt.datasets` uses (`BIOTAPY_DATA_DIR` if set, otherwise pooch's
  +per-user cache directory). `microformer-mgm` itself is never installed.
  +
  +MGM reads genera. Each feature's `var["genus"]` is read as MGM reads a
  +`g__<genus>` column: up to the first character that is not a letter, digit
  +or underscore, so `Escherichia-Shigella` is MGM's `Escherichia`. Features of
  +the same genus are summed, so a table at any level works, and
  +`bt.pp.tax_glom(tdata, "genus")` first gives the same embedding. Features
  +without a genus, or whose genus MGM never saw, are left out with one warning
  +that counts them; on GlobalPatterns that is 96 of 996 genus-level features.
  +Each sample then goes through MGM's own preprocessing: relative abundances
  +over the genera MGM knows, standardised with MGM's per-genus mean and standard
  +deviation, genera sorted by that value, and the sentence `<bos>`, genera,
  +`<eos>` cut to 512 tokens. A sample with none of MGM's genera is embedded
  +from `<bos> <eos>` alone, with a warning naming it.
  +
  +The embedding is the mean of the last hidden layer over the sample's tokens,
  +the "element-wise mean pooling" MGM's paper uses for the pretrained model:
  +256 float32 values per sample, within 2e-6 of MGM 0.5.8's own code. The model
  +runs on the CPU, one sample at a time; on GlobalPatterns' genus profiles that
  +is about 60 samples a second on 8 threads, with no memory beyond the model's.
  +On a CPU running more than four threads, the first call in a session can
  +differ from later ones by up to about 2e-4: torch 2.13 and 2.14 sometimes
  +compute their first `tanh` less precisely.
  +
  +Cite MGM when you publish results that use it: Zhang H, Zhang Y, Kang Z,
  +Xiong J, Yang R, Ning K (2026) MGM as a large-scale pretrained foundation
  +model for microbiome analyses in diverse contexts. *Adv Sci* 13:e13333.
  +
  +### Filtering and leakage
  +
   A model reads what it was trained on, so filter first: a later feature change
   (`bt.pp.filter_features`, `bt.pp.tax_glom`) drops `obsm`, and the embedding
   with it. MGM embeds each sample from that sample alone, so computing its
  ````
  The numbers are the prototype's: 96 of 996 from
  `test_embeds_global_patterns_end_to_end`, 520 GlobalPatterns profiles in
  8.8 s on 8 threads (about 60 a second), the 1.6e-4 first-call difference
  (design, "a torch first-call race").
- [x] **Step 7: Knowledge.** (`generated` and `commit:` as in 4.C0; the
  frontmatter hunks are not shown.)
  ```diff
  diff --git a/.knowledge/contracts/r-golden-parity.md b/.knowledge/contracts/r-golden-parity.md
  @@ -32,6 +32,13 @@ sources:
      database, so there is no container. Output:
      `tests/golden/humann/<name>.csv.gz` (samples as rows, row ids without
      their `": name"`) and `tests/golden/humann/VERSIONS.txt`.
  +1c. MGM's reference embeddings (`ml.embed(..., "mgm")`) are produced by
  +   `tests/mgm/export_reference.py` with MGM's own package, run with
  +   `uv run --no-project --python 3.11 --with microformer-mgm==0.5.8 --with
  +   torch==2.0.1+cpu` (and PyTorch's CPU index), never in CI: MGM pins numpy
  +   1.24 and torch 2.0, so it cannot share biotapy's environment. Output:
  +   `tests/data/mgm/counts.csv` and `embeddings.csv`, compared at `atol=1e-3`
  +   (`tests/ml/test_mgm.py`, marker `mgm`, which gives the reason).
   2. Output: `tests/golden/<dataset>/<function>.csv.gz`, samples as rows
      (see [samples-as-rows](/decisions/samples-as-rows.md)), plus
      `tests/golden/VERSIONS.txt` listing R and package versions. Golden files
  @@ -65,7 +72,8 @@ sources:
      stay synthetic, except small files copied under a permissive licence
      with a `NOTICE.txt` beside them (`tests/data/humann`: HUMAnN's MIT test
      data; `tests/data/metaphlan`: a MetaPhlAn 4.0.6 profile from HUMAnN's MIT
  -   test data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0). PICRUSt2
  +   test data; `tests/data/enzyme`: an ENZYME excerpt, CC BY 4.0;
  +   `tests/data/mgm`: MGM's embeddings of a synthetic genus table, MIT). PICRUSt2
      (GPL-3) fixtures are always synthetic, written from its documented
      column headers.
   7. `pl` functions have an R equivalent but no golden test. They draw numbers
  diff --git a/.knowledge/decisions/optional-heavy-dependencies.md b/.knowledge/decisions/optional-heavy-dependencies.md
  @@ -121,7 +121,8 @@ torch or an R installation into every install is unacceptable.[^spec]
   - CI imports every module with no extras installed (`import-without-extras`),
     so a lazy import leaking to module level fails fast. Each extra has a job
     that installs it and runs its marker, deselected everywhere else: `r-bridge`
  -  (`r`) and `ml-extras` (`torch`, Linux, Python 3.13, 30-minute timeout), both in `check.needs`.
  +  (`r`) and `ml-extras` (`torch`, then `mgm` with a pooch cache of its own for MGM's weights and
  +  GlobalPatterns; Linux, Python 3.13, 30-minute timeout), both in `check.needs`.
   - `uv.lock` is not committed, so every CI job that runs `uv run` resolves
     afresh, and reads torch's versions from download.pytorch.org even when it
     installs no torch (measured: `uv lock` 0.10 s -> 0.37 s with a warm cache).
  diff --git a/.knowledge/modules/core.md b/.knowledge/modules/core.md
  @@ -138,8 +138,8 @@ none of them back.
     `np.random.Generator` without touching global RNG state.
   - `_download.py:make_pooch` - the `pooch.Pooch` every run-time download goes
     through: `BIOTAPY_DATA_DIR` if set, else `pooch.os_cache("biotapy")`.
  -  Used by `datasets/_remote.py` (phyloseq, ENZYME and HMP2 files); every
  -  later download goes through it too, so all of them land in one cache that
  +  Used by `datasets/_remote.py` (phyloseq, ENZYME and HMP2 files) and
  +  `ml/_mgm.py` (MGM's weights), so every download lands in one cache that
     CI and tests point elsewhere with one variable. Building the pooch
     downloads and creates nothing; `fetch` does.

  diff --git a/.knowledge/modules/ml.md b/.knowledge/modules/ml.md
  @@ -31,6 +31,9 @@ that a plugin registers ([embedding-plugins](/decisions/embedding-plugins.md)).
     entry-point group `biotapy.embeddings` (`_embed.py:GROUP`), calls it with
     the AnnData and returns its samples x dimensions array, or writes
     `obsm["X_<model>"]` with `inplace=True`.
  +- `_mgm.py:embed` - MGM, the reference plugin, registered by biotapy's own
  +  `pyproject.toml` as `mgm` (extra `mgm`): genus tokens from `var["genus"]`,
  +  MGM's pretrained GPT-2, the mean of the last hidden layer per sample.

   # Invariants

  @@ -73,11 +76,40 @@ that a plugin registers ([embedding-plugins](/decisions/embedding-plugins.md)).
   - `embed` trusts no plugin: the result must be a plain `numpy.ndarray` (not
     a masked array or a matrix), 2-D with one row per sample and at least one
     column, float and finite, or it raises naming the plugin, before anything
     is written; a result that shares memory with `X`, a layer, an `obsm`, `varm`,
     `obsp` or `varp` entry or a top-level `uns` array (not one nested deeper) is
     copied. A plugin's own exception keeps its type and gains a note naming the
     plugin. The model name must be letters, digits, `_`, `-` or `.`.
     `_embed.py:embed`, `_embed.py:_checked`.
  +- MGM's tokens are MGM's own (`mgm/src/MicroCorpus.py`, 0.5.8): a genus is
  +  read with MGM's regex `g__[A-Za-z0-9_]+` on `"g__" + genus`; features of
  +  one token are summed (`_core.sum_by`); relative abundance is taken over
  +  MGM's genera only, as MGM drops the others before dividing; a genus is kept
  +  when its standardised value `(rel - mean) / std` exceeds that of zero
  +  abundance, and the kept ones are sorted by it with pandas' own
  +  `sort_values(ascending=False)` on the vocabulary-ordered series (the order
  +  of tied genera follows numpy's quicksort, which can differ from MGM's numpy
  +  1.24 environment; GlobalPatterns and enterotype have no such tie); `<bos>` ... `<eos>` is cut to 512 tokens (the `<eos>` is
  +  lost past 510 genera). Vocabulary ids are the position in `phylogeny.csv`
  +  plus 4 (`<pad>`, `<mask>`, `<bos>`, `<eos>`); the pickled tokenizer is
  +  never loaded. `_mgm.py:_sentences`.
  +- Where MGM would drop a sample (no count in its vocabulary), the plugin
  +  embeds `<bos> <eos>` and warns naming it, because `ml.embed` needs one row
  +  per sample. Left-out features get one warning counting them.
  +  `_mgm.py:_sentences`.
  +- MGM checks `var["genus"]` and `X` (finite, non-negative) before it imports
  +  torch or transformers or downloads anything, so those errors are tested in
  +  every CI job. `_mgm.py:embed`.
  +- The weights come from microformer-mgm 0.5.8's wheel on files.pythonhosted.org
  +  (SHA-256 pinned), fetched through `_core.make_pooch` and unzipped by
  +  `pooch.Unzip` to `<cache>/microformer_mgm-0.5.8-py3-none-any.whl.unzip/`;
  +  `torch.load(..., weights_only=True)` reads the GPT2LMHeadModel checkpoint
  +  and its `transformer.` keys load into `transformers.GPT2Model` strictly
  +  (`lm_head.weight` dropped). `_mgm.py:embed`.
  +- One sample at a time, unpadded, under `torch.inference_mode()`, model in
  +  `eval()` (dropout off): measured on 520 GlobalPatterns-like profiles, batch
  +  1 took 8.8 s, batches of 8 and 64 took 17.0 s and 18.5 s and needed 196 MB
  +  and 1.6 GB more (8 threads). `_mgm.py:embed`.
   - Inherited scikit-learn methods (`transform`, `fit_transform`,
     `get_support`, `get_feature_names_out`, `set_output`) are named in each
     class's `Notes`, because the class template leaves inherited members off
  @@ -90,6 +122,8 @@ that a plugin registers ([embedding-plugins](/decisions/embedding-plugins.md)).
     `OneToOneFeatureMixin`, `validate_data`. scikit-bio: `clr`.
   - torch (extra `torch`), only through `import_optional` inside `to_torch`
     ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).
  +- torch and transformers (extra `mgm`), only through `import_optional` inside
  +  `_mgm.py:embed`; pooch (`Unzip`) and `_core.make_pooch` for its weights.
   - `ml` imports no sibling top-layer module (`pl`, `da`) and not `pp`; the
     link to `pp` is through `_core` only
     ([module-boundaries](/contracts/module-boundaries.md)).
  @@ -100,7 +134,9 @@ that a plugin registers ([embedding-plugins](/decisions/embedding-plugins.md)).
   `pp` parity tests, a `Pipeline` cross-validation test, `to_torch`'s argument
   errors, `embed` with fake plugins). `uv run --group test --extra torch pytest -m torch` runs the
   `to_torch` tests and its docstring example, as CI's `ml-extras` job does (30-minute timeout,
  -`.github/workflows/test.yaml`). The pseudocount
  +`.github/workflows/test.yaml`). `BIOTAPY_DATA_DIR=.pooch uv run --group test --extra mgm pytest -m mgm`
  +runs MGM against its own embeddings (`tests/data/mgm`), the end-to-end GlobalPatterns test (exit gate 2)
  +and `embed`'s docstring example, as the job's last step does. The pseudocount
   warning's text is unit-tested in `tests/core/test_composition.py`.

   # Gotchas
  @@ -154,6 +190,19 @@ warning's text is unit-tested in `tests/core/test_composition.py`.
   - `embed`'s docstring example runs MGM, so the root `conftest.py` gives
     `biotapy.ml._embed`'s doctests the marker `mgm`, deselected by default like
     `torch`. `conftest.py:_EXTRA_DOCTESTS`.
  +- MGM's tests compare at `atol=1e-3`: torch 2.13 and 2.14 sometimes compute
  +  their first `tanh` in a process less precisely on a CPU running more than
  +  four threads (measured: 9e-5 on `tanh(-5)`, never with 4 threads or fewer,
  +  never after a first call), which moves that call's embedding by up to
  +  1.6e-4; otherwise biotapy and MGM 0.5.8 agree to 1.7e-6. The mutations the
  +  parity test was checked against (wrong denominator, ascending sort,
  +  last-token pooling, dropout on, `z > 0`, vocabulary order) each fail it.
  +  `tests/ml/test_mgm.py:ATOL`.
  +- `microformer-mgm` cannot be installed beside biotapy (its pins), so the
  +  reference embeddings are written by `tests/mgm/export_reference.py` in a
  +  Python 3.11 environment of its own (`uv run --no-project --python 3.11 --with
  +  microformer-mgm==0.5.8 ...`); the script imports MGM inside its functions,
  +  because pytest's `--doctest-modules` imports every file under `tests/`.
   - anndata 0.13 lists `X` as `layers[None]`, so `list(adata.layers)` holds
     `None` even when no layer was added; `to_torch`'s missing-layer error does
     not list the layers. `_torch.py:_table`.
  ```
  Tick 4.4b here; log line, first in the slice 4C section:
  ```markdown
  - **Update**: [ml](modules/ml.md) gains MGM (`_mgm.py`: MGM's own tokens, the `<bos> <eos>` sample, validation before the imports, the weights from the wheel, one sample at a time and why, the parity tolerance and the torch `tanh` race, the reference script's environment); [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md): `ml-extras` runs `mgm` with its own pooch cache; [r-golden-parity](contracts/r-golden-parity.md) gains MGM's reference embeddings (item 1c, the `tests/data/mgm` exception); [core](modules/core.md): `make_pooch` serves `ml/_mgm.py` too; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.4b.
  ```
- [x] **Step 8: Gate and commit**
  ```bash
  git add src/biotapy/ml/_mgm.py src/biotapy/ml/_embed.py pyproject.toml tests/ml/test_mgm.py \
    tests/data/mgm/counts.csv tests/data/mgm/embeddings.csv tests/data/mgm/NOTICE.txt tests/mgm/export_reference.py \
    .github/workflows/test.yaml tests/test_ci.py docs/guide/machine_learning.md docs/contributing.md \
    .knowledge/modules/ml.md .knowledge/contracts/r-golden-parity.md .knowledge/decisions/optional-heavy-dependencies.md \
    .knowledge/modules/core.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  # the slice gate, with its -m mgm line; zizmor checks the new cache step
  git commit -m "feat(ml): add MGM, the reference embedding model, behind the extra mgm

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1488 passed, 2 skipped, 92 deselected`; `36
  passed, 1546 deselected`; `build succeeded`; `25 passed, 1557
  deselected`; `13 passed, 1569 deselected`.

### Checkpoint C - review slice 4C
- [x] Review the whole slice (superpowers:requesting-code-review) against
  every contract, pure-by-default, optional-heavy-dependencies, the new
  embedding-plugins decision, the Phase 4 review focus (item 5) and the
  slice 4C review focus; then a fix pass, one commit per finding, each with
  a test. Record the counts and the fix range here.
  - Result: 0 Critical, 0 Important, 5 Minor (+1 counted); the controller
    re-graded Minor 1 (a damaged extracted MGM file) to Important. Fix pass
    `f9b3f50..c87fc13` (11 commits), then a follow-up round `c87fc13..31aa11d`
    (4 commits). Scoped re-review: all 8 findings addressed; one new minor and
    three nits fixed in the follow-up.
- [x] Run the slice's checks at the last commit and record them: the slice
  gate's six counts; `uv run --group test --extra mgm coverage run -m pytest
  -m "mgm or not mgm" tests/ml/test_embed.py tests/ml/test_mgm.py
  tests/core/test_download.py src/biotapy/ml/_embed.py` then `coverage
  report --include
  "src/biotapy/ml/_embed.py,src/biotapy/ml/_mgm.py,src/biotapy/_core/_download.py"`
  (exits 0; 100% each: `_embed.py` 54, `_mgm.py` 87 and `_download.py` 4 statements);
  `tests/ml/test_embed.py` under `--hypothesis-seed=1`, `2`, `3`; the six
  parity mutations of 4.4b Step 5, each failing.
  - At `31aa11d`: `pytest -q -W error::UserWarning` 1488 passed, 2 skipped,
    92 deselected; `-m "golden or network"` 36 passed; `--extra torch -m torch`
    25 passed; `--extra mgm -m mgm` 13 passed; coverage 100% on `_download.py`
    (4 statements), `_embed.py` (54) and `_mgm.py` (87). `tests/ml/test_embed.py`
    with `--hypothesis-seed` 1, 2 and 3: 37 passed each, at `55c5f2e`.
    Mutations (Checkpoint C review and re-review): all 6 on `_embed.py`
    fail a test; on `_mgm.py`, last-token pooling, dropped `<bos>`/`<eos>`,
    five token-order changes, a wrong denominator, dropout left on,
    genus ids off by one, no per-genus summing, and keeping `<eos>` at the
    512 cut each fail parity (the last only at the warm 1e-5 check).
- [x] Knowledge: [ml](/modules/ml.md), [core](/modules/core.md),
  [embedding-plugins](/decisions/embedding-plugins.md),
  [optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md),
  [data-model-slots](/contracts/data-model-slots.md),
  [r-golden-parity](/contracts/r-golden-parity.md) and
  [pure-by-default](/decisions/pure-by-default.md) already carry 4C;
  re-check them against the fix pass and bump only what changed, with log
  lines; run `scripts/knowledge_stale.sh` and re-check what it flags.
- [x] Push the branch and open the PR only after the user approves that push
  (R13.3). The PR body carries the R9.2 reasons: "Extra `mgm` (`torch>=2.9`,
  `transformers>=5`; transformers Apache-2.0; approved 2026-10-09,
  decision 13): MGM, the reference model of `ml.embed`, is a GPT-2 that
  transformers runs on torch; an extra, never core (R9.3). It adds 16
  packages to the lock (huggingface-hub, tokenizers, safetensors, regex,
  tqdm, typer and their dependencies; Apache-2.0, MIT, BSD-3-Clause,
  MPL-2.0, ISC licences), about 97 MB installed beside torch; under uv torch
  still comes from the CPU index. Run-time data: MGM 0.5.8's pretrained
  weights (MIT, Copyright (c) 2024 NingLab), downloaded once by pooch from
  `microformer_mgm-0.5.8-py3-none-any.whl` on files.pythonhosted.org,
  pinned by SHA-256, never bundled; `microformer-mgm` itself is not a
  dependency (its pins conflict with numpy 2 and pandas 3)." CI green,
  including `ml-extras`; record its runtime and the torch and transformers
  versions it installed here (the prototype's local replay: torch steps
  11.8 s, MGM steps 21.2 s with a cold pooch cache).
  - Done 2026-10-09: PR #30 merged as `28de1ad`, Test run `37962120924` 22/22 jobs green; `ml-extras` (job `113927302730`, 86 s) installed torch 2.14.1+cpu (`25 passed` in 12.96 s) and transformers 5.19.0 (MGM `13 passed` in 28.92 s on a cold cache).
- [x] Ask the user to review slice 4C, and to confirm
  [embedding-plugins](/decisions/embedding-plugins.md) (`draft` until then),
  before slice 4D is expanded.
  - Approved 2026-10-09; the user confirmed embedding-plugins.

---

## Slice 4D - Docs and release

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax for
> tracking. This section replaces the slice 4D outline (rules.md R1.2a).

**Goal:** two small integer-overflow fixes land first; then a user can read three executed or
CI-checked tutorials (leak-free cross-validation on HMP2, which is exit gate 1; mmvec on a real
pair of microbe and metabolite tables; MGM embeddings of GlobalPatterns, the docs side of exit
gate 2), finds `bt.io.to_mudata` on the Coming-from-R page, reads current knowledge concepts, and
gets biotapy 0.4.0 on PyPI with the extras `torch` and `mgm`, every outward step waiting for the
user's explicit approval.

**Architecture:** two `fix(core)` commits (4.F4 `sum_by`, 4.F5 `sum_pairs`); three tutorial
pages under `docs/tutorials/`, each with tests that keep what its prose quotes equal to what its
code gives (a shared `run_page` fixture runs a page's own code in the `network` or `mgm` test, a
default-run test checks the quoted text); one new public loader, `bt.datasets.biocrust()` (4.6c,
decision 37), because the multi-omics tutorial needs real paired data; the Coming-from-R pins
and the "not in 0.4" labels; knowledge; the release.

**Tech stack:** sphinx 9.1 · myst-nb 1.4 (`nb_execution_mode = "cache"`, `nb_execution_timeout
= 300`, unchanged) · scikit-learn 1.9.1 (`RepeatedStratifiedKFold`, `cross_val_score`,
`SelectKBest(f_classif)`, `LogisticRegression`) · scikit-bio 0.7.4 (`mmvec`) · pooch 1.9.0 ·
torch 2.14.1+cpu and transformers 5.19.0 (extra `mgm`, CI only) · hatchling/`uv build` · twine.
No new dependency (R9.1).

**Spec:** this concept's design notes 9 (leak-free notebook) and 10 (exit gate 2), decisions 16
to 18 and 31, the slice 4D outline it replaces, rules.md,
[cut-a-release](/playbooks/cut-a-release.md), and Phase 3's slice 3D (how 0.3.0 was documented
and released).

**How slice 4D was checked.** Every file below was written into a scratch clone of the
repository at `28de1ad` (master after PR #30, slice 4C) on branch `phase-4d`, after a stand-in
`docs(roadmap)` commit holding only the checklist lines 4.F4 and 4.F5, the Checkpoint C tick and
4.F1's upstream link (`b5654af`); the release commits sit on `release-0.4.0` branched from it.
Each commit was gated on its committed tree (`git status --short` empty) with `.venv` synced
without extras (`uv sync --all-groups`), as CI's default jobs are; `-m torch` ran with `--extra
torch`, `-m mgm` with `--extra mgm` (torch 2.14.1+cpu, transformers 5.19.0, CPython 3.13.2).
Every run exported absolute `BIOTAPY_DATA_DIR`, `HF_HOME`, `UV_CACHE_DIR`, `XDG_CACHE_HOME` and
`MPLCONFIGDIR` under the session scratchpad, starting from an empty pooch cache;
`~/.cache/biotapy` never appeared and no gate wrote under `~/.cache`.

| Commit (scratch hash) | `uvx prek run --all-files` | `pytest -q -W error::UserWarning` | `pytest -q -m "golden or network"` | `sphinx-build -W` (wall; notebooks) | `--extra torch -m torch` | `--extra mgm -m mgm` | `knowledge_stale.sh --against HEAD` |
|---|---|---|---|---|---|---|---|
| base `28de1ad` | passed (14 hooks) | 1488 passed, 2 skipped, 92 deselected | 36 passed, 1546 deselected | build succeeded (59.9 s; 5) | 25 passed, 1557 deselected | 13 passed, 1569 deselected | 28 current, 0 stale |
| 4.F4 `d5ca5aa` | passed | 1496 passed, 2 skipped, 92 deselected | 36 passed, 1554 deselected | build succeeded (50.2 s; 5) | 25 passed, 1565 deselected | 13 passed, 1577 deselected | 17 current, 11 stale |
| 4.F5 `03d217c` | passed | 1504 passed, 2 skipped, 92 deselected | 36 passed, 1562 deselected | build succeeded (50.1 s; 5) | - | - | 17 current, 11 stale |
| 4.6 `e46f065` | passed | 1510 passed, 2 skipped, 93 deselected | 37 passed, 1568 deselected | build succeeded (60.6 s; 6, `leak_free_cv` 11.2 s) | - | - | 17 current, 11 stale |
| 4.6b `f5fb966` | passed | 1512 passed, 2 skipped, 95 deselected | 37 passed, 1572 deselected | build succeeded (67.8 s; 6) | 25 passed, 1584 deselected | 15 passed, 1594 deselected | 17 current, 11 stale |
| 4.6c `55392ef` | passed | 1521 passed, 2 skipped, 98 deselected | 40 passed, 1581 deselected | build succeeded (69.5 s; 7, `multiomics` 3.1 s) | - | - | 16 current, 12 stale |
| 4.D1 `b904fa6` | passed | 1524 passed, 2 skipped, 98 deselected | 40 passed, 1584 deselected | build succeeded (76.4 s; 7) | - | - | 16 current, 12 stale |
| 4.7 `a9937e9` | passed | 1524 passed, 2 skipped, 98 deselected | 40 passed, 1584 deselected | build succeeded (77.0 s; 7) | - | - | 28 current, 0 stale |
| 4.D2 release `20151ac` | passed | 1524 passed, 2 skipped, 98 deselected | 40 passed, 1584 deselected | build succeeded (75.5 s; 7) | 25 passed, 1599 deselected | 15 passed, 1609 deselected | 23 current, 5 stale |
| 4.D2 knowledge `ceada96` | passed | 1524 passed, 2 skipped, 98 deselected | - | - | - | - | 28 current, 0 stale |

- Every default run also ends in "2 warnings": scikit-bio's `RuntimeWarning: invalid value
  encountered in divide` from 4.F1's `ancombc2` tests, as on master. The 2 skips are
  scikit-learn's `check_array_api_input`, as on master.
- "-" means the commit touches nothing that run covers (4.F5 is `fn` only; 4.6, 4.6c, 4.D1 and
  4.7 touch no `ml` code). The stand-in `b5654af` was checked with `tests/test_knowledge_bundle.py`
  only (`82 passed`): it changes three concepts and no code.
- The stale counts between 4.F4 and 4.D1 are the bundle's convention, not a defect: a task's
  concept edit writes `commit:` as its parent, so the code change in the same commit stales it;
  4.7 re-stamps every flagged concept after reading it against the slice (Phase 3's Task 3.14
  did the same). The release commit re-stales the five concepts whose `paths` hold
  `pyproject.toml` or `CHANGELOG.md`; its knowledge commit re-stamps them.
- New default-run items: 4.F4 +8 (7 `sum_by` dtype cases, 1 `tax_glom`); 4.F5 +8; 4.6 +6 (the
  prose test and 5 `run_page` tests; the network test is the +1 deselected); 4.6b +2 (two prose
  tests; the two `mgm` tests are the +2 deselected); 4.6c +9 (3 loader tests, 2 link tests,
  1 figures-prose test, 1 tutorial prose test, the docstring check of `bt.datasets.biocrust` and its doctest; the
  loader's, the figures' and the tutorial's network tests are the +3 deselected); 4.D1 +3
  (2 MultiAssayExperiment cases, 1 no-row check).
- The docs job's notebooks: on PR #30's run (`37962120924`) the five existing notebooks took
  2.0-11.7 s on GitHub's runner, close to the local 2.3-16.8 s; `leak_free_cv` (11.2-13.7 s
  locally, 85 s of CPU time across threads in its network test) and `multiomics` (3.0-3.6 s)
  add roughly 15-20 s to the build. Read the
  Docs runs the same `hatch run docs:build` with the `doc` group and no extras and no pooch
  cache: `leak_free_cv` reuses the HMP2 download the function tutorial already makes (23 MB), and
  `multiomics` adds 135 KB from GitHub.
- Coverage (`coverage run -m pytest tests/datasets/test_biocrust.py tests/core/test_matrix.py`):
  `datasets/_biocrust.py` 9 statements, `_core/_matrix.py` 34, each 100%. Hypothesis seeds 1, 2
  and 3 on `tests/core/test_matrix.py tests/pp/test_glom.py tests/fn/test_glom.py`: `85 passed`
  each.
- The release state (`release-0.4.0`, `20151ac`): `uv build` and `uvx twine check --strict
  dist/*` PASSED for both files; the wheel's METADATA reads `Version: 0.4.0`, `Provides-Extra:
  mgm`, `Requires-Dist: torch>=2.9; extra == 'mgm'`, `Requires-Dist: transformers>=5; extra ==
  'mgm'`, `Provides-Extra: r`, `Requires-Dist: rpy2>=3.6.8; extra == 'r'`, `Provides-Extra:
  torch`, `Requires-Dist: torch>=2.9; extra == 'torch'`, and its `entry_points.txt` is
  `[biotapy.embeddings]` / `mgm = biotapy.ml._mgm:embed`; the sdist holds the root
  `conftest.py`, `tests/conftest.py`, the three new tutorials and `tests/datasets/test_biocrust.py`,
  and no `tests/mgm` or `benchmarks`; `pytest` inside the extracted sdist: `1405 passed, 2 skipped,
  96 deselected` (with a relative `XDG_CACHE_HOME` the sdist run fails
  `test_without_biotapy_data_dir_the_cache_is_pooch_s_per_user_directory`: `pooch.os_cache`
  keeps the `..` unnormalised; with absolute paths it passes, and the playbook now says so);
  `biotapy.__version__` prints `0.4.0`.
- **CI evidence already on record (PR #30, merged as `28de1ad` on 2026-10-09).** Test run
  `37962120924`, all 21 jobs `SUCCESS` (and the Build workflow's `package`): `ml-extras` (job
  `113927302730`, 86 s) installed torch 2.14.1+cpu (`25 passed, 1557 deselected in 12.96s`) and
  transformers 5.19.0 (`13 passed, 1569 deselected in 28.92s`), which holds
  `tests/ml/test_mgm.py::test_embeds_global_patterns_end_to_end`; the `docs` job
  (`113927302392`) executed five notebooks.
- **APIs and facts checked** (R2.2): scikit-learn 1.9.1 `RepeatedStratifiedKFold(n_splits=,
  n_repeats=, random_state=)`, `cross_val_score(..., scoring="roc_auc")`, `SelectKBest(f_classif,
  k=)`; no `ConvergenceWarning` with `max_iter=5000` (the page's code ran under `python -W
  error`); numpy 2's `np.zeros(0, dtype=d).sum().dtype` is `int64` for `int8`, `int16`, `int32`
  and `bool`, `uint64` for `uint8`, and the dtype itself for floats (NumPy's documented `sum`
  accumulator); `scipy.sparse` multiplies in the wider of two dtypes (`int8 @ int64` gives
  `int64`); a MyST page without a `jupytext`/`kernelspec` header is not executed by myst-nb (the
  embedding page builds with no "Executed notebook" line); biocore/mmvec's `examples/soils`
  (`microbes.biom` 81,409 bytes, `metabolites.biom` 53,681 bytes, unchanged since `7457c87`
  of 2019-10-17; repository BSD-3-Clause; master `88ca33b`) and its `check_soils.ipynb`, which
  asserts that 13 listed metabolites rank above zero for `rplo 1 (Cyanobacteria)`; mmvec's README
  calls it the "desert biocrust experiment" and says studies of about 19 samples need careful
  tuning of the latent dimension and priors.

### Slice 4D design

Each settled question gives the answer and the reason. **(user)** marks the ones repeated under
"Decisions for the user (slice 4D)".

#### Where the work goes

| File | Holds |
|---|---|
| `src/biotapy/_core/_matrix.py` (4.F4, 4.F5) | `sum_by` and `sum_pairs` sum in `numpy.sum`'s dtype for `X` |
| `tests/core/test_matrix.py`, `tests/pp/test_glom.py`, `tests/fn/test_glom.py` (4.F4, 4.F5) | the dtype cases and one public-API test each |
| `tests/conftest.py` (4.6) | `run_page`: runs a docs page's `{code-cell}`s and `python` blocks, returns the names they define |
| `docs/tutorials/leak_free_cv.md` (new, 4.6), `tests/ml/test_transformers.py` | exit gate 1; `LEAK_FREE_CV_AUC` and its two tests |
| `docs/tutorials/embeddings.md` (new, 4.6b), `tests/ml/test_mgm.py` | exit gate 2's page (not executed); `GLOBAL_PATTERNS_LEFT_OUT`, `GLOBAL_PATTERNS_SHAPE` and their two tests |
| `src/biotapy/datasets/_biocrust.py` (new, 4.6c), `_remote.py`, `__init__.py`, `tests/datasets/test_biocrust.py` | `bt.datasets.biocrust()` |
| `docs/tutorials/multiomics.md` (new, 4.6c), `tests/tl/test_mmvec.py` | the multi-omics tutorial; `MICROCOLEUS_METABOLITES` and its two tests |
| `docs/tutorials/index.md`, `docs/guide/{machine_learning,multiomics,datasets}.md`, `docs/api.md` | links, the `biocrust` section, the API entry |
| `tests/test_coming_from_r.py`, `docs/_data/r_idioms.toml`, `docs/coming_from_r.md`, `docs/tutorials/phyloseq_analysis.md` (4.D1) | the MultiAssayExperiment pins; "not in 0.4" |
| `.knowledge/...` (each task, 4.7, 4.D2) | `core`, `datasets`, `ml`, `tl`, `embedding-plugins`, `cut-a-release`, verification bumps |
| `pyproject.toml`, `CHANGELOG.md`, `README.md` (4.D2) | the release |

#### 1. 4.F4 and 4.F5: sums that cannot wrap

- **The defect.** `_core.sum_by` builds its 0/1 indicator in `X.dtype`, and SciPy multiplies
  in the wider of the two dtypes, so an `int8` table sums in `int8`: three features of 100 into
  one genus give 44 (300 - 256); `uint8` 600 gives 88, `int16` 60,000 gives -5,536, and a bool
  table gives `True` (an OR, not a count). `_core.sum_pairs` has the same line. Measured on the
  prototype through `bt.pp.tax_glom` (44 for 300) and `bt.fn.func_glom` (-56 for 200).
- **Callers.** `sum_by`: `pp.tax_glom` (`pp/_glom.py:69`, keeps the dtype in `X`),
  `pl.bar`'s `fill` by a `var` column (`pl/_abundance.py:_segments`, converts the small result to
  float64 after summing) and MGM (`ml/_mgm.py:109`, divides rows to float64 after summing).
  `sum_pairs`: `fn.func_glom` (`fn/_glom.py:122`, keeps the dtype). Every caller wants exact sums;
  two of them return the integer dtype to the user; no contract fixes `X`'s dtype
  (data-model-slots says only which values are counts).
- **The fix.** The indicator takes the dtype `numpy.sum` gives `X`'s
  (`np.zeros(0, dtype=X.dtype).sum().dtype`): bool and signed integers narrower than 64 bits
  become `int64`, unsigned ones `uint64`, floats keep theirs, so a `float32` table stays
  `float32` and an `int64` table stays `int64` (`test_sum_by_keeps_integer_dtype` holds). This is
  the rule a user gets from `X.sum()`, with no new branch. One visible change: an `int32` table's
  sums become `int64` (CHANGELOG "Changed"). Rejected: always `float64` (turns `tax_glom`'s
  integer counts into floats, against `test_sum_by_keeps_integer_dtype`);
  always `int64` for integers and `X.dtype` for floats (sends `uint64` to float64 through SciPy's
  upcast).
- **4.F4 is the approved fix; 4.F5 is the same line in `sum_pairs` (decision 36).** They are two
  commits so 4.F5 can be dropped without touching 4.F4. Both run first, because the tutorials
  below use `tax_glom`.

#### 2. The leak-free notebook (4.6, exit gate 1)

- **Data and question:** `bt.datasets.hmp2()["taxa"]`, 130 participants' first stool samples,
  579 MetaPhlAn 3 species as relative abundances (each sample sums to 1; smallest non-zero value
  1.5e-6; 206 species absent from all samples), IBD (CD 65 + UC 38 = 103) against non-IBD (27).
- **Pipeline:** `PrevalenceFilter(min_prevalence=0.1)` (128 species kept on the whole table),
  `CLR(pseudocount=1e-6)` (below the smallest value, so no pseudocount warning),
  `LogisticRegression(max_iter=5000)` (no `ConvergenceWarning`), mean ROC AUC.
- **Cross-validation: `RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)`.
  (user, decision 38)** The outline's single 5-fold split (design note 9's numbers, reproduced
  exactly: honest 0.478, filter outside 0.475, `SelectKBest` inside 0.571, outside 0.722) puts the
  honest pipeline below 0.5, which a reader takes for a bug, and moves by about 0.05 with the
  split. Measured on the prototype (mean ROC AUC):

  | Setting | honest | filter outside | `SelectKBest` inside | `SelectKBest` outside | shuffled: inside | shuffled: outside |
  |---|---|---|---|---|---|---|
  | species, 5-fold (outline) | 0.478 | 0.475 | 0.571 | 0.722 | 0.531 | 0.764 |
  | species, 5 x 5 (chosen) | 0.527 | 0.526 | 0.558 | 0.705 | 0.544 | 0.767 |
  | genus, 5-fold | 0.535 | 0.517 | 0.577 | 0.736 | 0.530 | 0.699 |
  | genus, 5 x 5 | 0.586 | 0.579 | 0.629 | 0.769 | 0.550 | 0.709 |

  (Shuffled columns: the mean over 10 permutations, `np.random.default_rng(0)`; their standard
  deviation over permutations is 0.04-0.08.) Other contrasts give no stronger honest signal with
  logistic regression on 5 x 5 (CD vs non-IBD 0.580, UC vs non-IBD 0.497; a 300-tree random
  forest reaches 0.636 for IBD), so the page keeps the approved question and says plainly that
  the honest answer is near chance: that is the result the leak would have hidden.
- **Shuffled labels: 10 permutations, averaged** (decision 38), instead of design note 9's one
  shuffle, whose number (0.754 vs 0.550) depends on which permutation was drawn.
- **The six numbers the prose quotes**: honest 0.527; filter outside 0.526 ("no measurable
  leak"); `SelectKBest(f_classif, k=20)` inside 0.558, outside 0.705; on shuffled labels inside
  0.544, outside 0.767 ("higher than on the real labels"). The page's
  `mean_auc` returns the unrounded `float`, and the cells that show a number wrap it in
  `round(..., 3)`, so the outputs print `0.527`, not `np.float64(...)`.
- **How the numbers stay honest. (user, decision 39)** `tests/ml/test_transformers.py` holds
  `LEAK_FREE_CV_AUC`. A `network` test runs the page's own cells through the new `run_page`
  fixture and compares the six unrounded names it defines with the constants (`approx(abs=1e-3)`,
  because one ranking flip moves a 25-fold mean by about 3.6e-4 and the quoted figures have three
  decimals), pins the "coin toss" wording (|shuffled inside - 0.5| < 0.05, inside < 0.6) and the
  range of the ten shuffled scores (`LEAK_FREE_CV_SHUFFLED_RANGE`, with the honest scores inside
  it); a default-run test checks the
  four prose sentences quote the constants and that the words hold (|inside - outside| < 0.01,
  outside - inside > 0.1 for selection, shuffled outside > real outside). A scikit-learn change
  that moves a number fails the `network` job; the constants and the page then change together.
  Phase 3's `FOUR_METHOD_*` shared constants only between the `r` test and a prose test, and the
  `r` test repeated the tutorial's code; running the page itself removes that copy. Read the
  Docs and the `docs` job execute the page; neither compares numbers.
- Runtime: 11-14 s locally (load 5 s, real labels 1.5 s, shuffles 6 s); HMP2's first call
  downloads 23 MB when the cache is cold, which the function tutorial already does in the same
  build; `hmp2()` peaks at about 1 GB, as in the function tutorial.

#### 3. The embedding page (4.6b, exit gate 2's docs)

- **A tutorial page, `docs/tutorials/embeddings.md`, plain Markdown, not executed. (user,
  decision 40)** It has no `jupytext`/`kernelspec` header, so myst-nb does not run it; the docs
  build has no torch (decision 17). Its `python` blocks are what a reader runs; its `text`
  blocks quote the outputs of `test_embeds_global_patterns_end_to_end`: the leave-out warning
  (96 of 996, the five names), `((26, 256), dtype('float32'))`, and 26 samples whose nearest
  neighbour shares their `SampleType`. A note says the page is not run by the docs build and that
  CI's `ml-extras` job checks it. The guide keeps MGM's input rules, licence, citation and the
  plugin `pyproject.toml` lines; the page links there.
- **How it stays honest.** `tests/ml/test_mgm.py` gains `GLOBAL_PATTERNS_LEFT_OUT` (the whole
  warning) and `GLOBAL_PATTERNS_SHAPE`; the exit-gate test matches the whole warning through
  them (it matched a prefix before, so this tightens it); an `mgm` test runs the page's `python`
  blocks through `run_page` and checks the warning, the shape, the dtype and `same_type == 26`; a
  default-run test checks the three `text` blocks and the sentence.
- **Timings quoted:** on the prototype's laptop (8 threads) the first `embed` took 7 s (importing
  torch and loading the model) and a second 1 s for 26 samples; CI's `ml-extras` job on PR #30
  ran the 13 MGM tests in 28.9 s including the 33 MB download. The page quotes the laptop
  numbers as a measurement with its versions; no test pins times.
- **Found while writing it:** of the 96 genus-level features MGM leaves out of GlobalPatterns, 46
  are clone names (`4-29`, `BD2-13`, ...), 22 other genera, and 28 are *Candidatus* genera that
  Greengenes writes as one word (`CandidatusPelagibacter`); 23 of those 28 are in MGM's
  vocabulary as `Candidatus_<Name>`. MGM's own regex reads them as written, so biotapy matches
  MGM's code (decision 30). The page and `ml.md` say so; changing it is decision 43.

#### 4. The multi-omics tutorial and `bt.datasets.biocrust()` (4.6c)

- **Decision 18 said the mmvec tutorial waits for 0.5 unless 4D finds a small public paired
  dataset. It found one (user, decision 37):** the example in mmvec's own repository,
  `biocore/mmvec` `examples/soils`: a desert biological soil crust from four successional stages
  (Swenson et al. 2018) sampled at five times after wetting, 466 microbe counts in 20 samples and 85 metabolite intensities in 19 (sample
  `9hr_late` has no metabolite profile), 135 KB in two BIOM files, BSD-3-Clause repository.
  `bt.tl.mmvec(seed=0)` fits it in 0.84 s.
- **Validation it brings:** mmvec's `check_soils.ipynb` asserts that 13 metabolites it lists for
  `rplo 1 (Cyanobacteria)` (*Microcoleus* sp.) all rank above zero. biotapy's fit (scikit-bio's
  L-BFGS, 3 dimensions; the example used TensorFlow's Adam and 1 dimension) gives 13 of 13 for
  seeds 0-4, with the same top four (adenine, adenosine, N6-acetyl-lysine,
  4-guanidinobutanoate). `bt.tl.mmvec` had no check on real data (no R implementation); this is
  one against the original authors' own example.
- **A loader, not a download in the notebook.** A tutorial that calls `pooch.retrieve` itself
  either ignores `BIOTAPY_DATA_DIR` (writing under `~/.cache` in every gate and docs build) or
  spells the cache logic out in a tutorial. `bt.datasets.biocrust() -> MuData` goes through
  `_remote.py`'s pooch like `hmp2`: registry names `biocrust_microbes.biom` and
  `biocrust_metabolites.biom`, URLs pinned to commit `88ca33b`, SHA-256
  `7f634b8d...f4a7` and `cfea3e72...2c72`; it fetches both before reading either, reads each with
  `bt.io.read_biom` (microbes `x_kind="counts"`, metabolites `"abundance"`), intersects the
  samples itself (so `io.to_mudata` does not warn about `9hr_late`) and returns
  `to_mudata({"taxa": ..., "metabolites": ...})`. No metadata columns: the sample names hold time
  and successional stage, and nothing asks for parsed columns (R2.3).
- **Tests:** three offline tests with synthetic BIOM files written by `bt.io.write_biom` (no
  data copied; the fixture's values are invented), patching `_biocrust._fetch` as `test_hmp2.py`
  does; two default-run tests that the guide and the tutorial link mmvec's example at the pinned
  commit and one that the docstring, the guide and the tutorial quote 135 KB, 20 samples and the
  most abundant microbe; two `network` tests, one of the shapes and one of those three figures
  on the downloaded files (shared constants). The tutorial's `network` test runs the page through `run_page` and checks its
  `microcoleus` set equals `MICROCOLEUS_METABOLITES` and `above_zero == 13`; a default-run test
  checks the two sentences that quote 13.
- **Docs:** the guide's datasets page gains a `biocrust` section, its caching and licensing
  paragraphs a sentence each; `docs/api.md` an entry; the multi-omics guide a link. The tutorial
  executes in 3 s.

#### 5. Coming from R (4.D1)

`io.to_mudata`'s docstring names two R calls, `MultiAssayExperiment::MultiAssayExperiment` and
`MultiAssayExperiment::intersectColumns`; every other Phase 4 public name (`tl.mmvec`,
`ml.PrevalenceFilter`, `ml.CLR`, `ml.to_torch`, `ml.embed`, `datasets.biocrust`) says `R
equivalent: none`. The table grows from 58 to 60 rows. A parametrized test pins the two rows,
one test checks no `bt.ml.`, `bt.tl.mmvec` or `bt.datasets.biocrust` cell appears, and the five
"not in 0.3" labels become "not in 0.4" (cut-a-release step 2c) with the page sentence and the
phyloseq vignette's three "0.3" statements, as Task 3.12 did.

#### 6. Knowledge (4.7) and the release (4.D2)

- **4.7** after Checkpoint D's fix pass: `ml.md` (Verification names the two tutorials and their
  tests; three gotchas: tutorial numbers in test constants checked by running the page, the size
  of the leak on HMP2, the one-word *Candidatus* genera), `tl.md` (mmvec's real-data check),
  `embedding-plugins.md` `draft` -> `stable` (the user confirms it by approving this plan; Checkpoint
  C's last box), and `generated`/`commit` bumps for the 11 concepts `knowledge_stale.sh` flags
  after reading each against the slice. `verified` is the user's alone (R12.3; decision 44).
  `core.md` and `datasets.md` are edited by 4.F4/4.F5 and 4.6c themselves.
- **4.D2** follows cut-a-release and Task 3.15: `## [Unreleased]` is empty again, so the 0.4.0
  entries are written from `git log v0.3.0..master --no-merges` (decision 41): Added for the
  `feat:` commits; Changed for the 0.3.0 behaviours that moved (the `(pp.clr)`/`(pp.philr)`
  suffix of the pseudocount warning, `5978028`; `filter_features`'s `TypeError`, `8e1567e`;
  4.F4/4.F5's 64-bit sums; "not in 0.4"). README: "0.4", `to_mudata`, `biocrust`, `tl.mmvec`, a
  Machine learning line, the `torch` and `mgm` extras under Installation, and "Next" from Phase 5
  (an unplanned backlog). The wheel check adds `entry_points.txt`.
- **Approvals (R13.3; decision 20 gives no standing approval for Phase 4; decision 42):**
  `phase-4d`'s push, PR and merge at Checkpoint D with the slice review; `release-0.4.0`'s push,
  PR and merge in the same message; then the tag, the tag push, the GitHub release (the PyPI
  upload) and `close-phase-4`'s push, PR and merge in one message after the release PR merges.

#### 7. Which commit ticks which box

| Box | Ticked in | Evidence |
|---|---|---|
| Checkpoint C's push box and checklist line | the plan commit | PR #30 merged as `28de1ad`, Test run `37962120924` 21/21 jobs green; `ml-extras` 86 s, torch 2.14.1+cpu, transformers 5.19.0 |
| Checkpoint C's "ask the user to review 4C" box | the plan commit | the user's approval of this plan |
| 4.F1's "Upstream" box | the plan commit | [scikit-bio#2631](https://github.com/scikit-bio/scikit-bio/issues/2631) |
| Exit gate 2, one foundation model end to end | the plan commit | PR #30's `ml-extras` job ran `test_embeds_global_patterns_end_to_end` (`13 passed`) |
| 4.F4, 4.F5, 4.6, 4.6b, 4.6c, 4.D1, 4.7 checklist lines | each task's last commit | the gates above |
| Checkpoint D (checklist and boxes) | the 4.D2 knowledge commit | the `phase-4d` PR merged green |
| Exit gate 1, executed notebook in the CI docs job | the 4.D2 knowledge commit | that PR's `docs` job log: `tutorials/leak_free_cv.md: Executed notebook` |
| Exit gate 3, Phase 1-3 gates green | the 4.D2 knowledge commit | every job of that PR's Test run green |
| 4.D2 Steps 1-9 | the 4.D2 knowledge commit | |
| 4.D2 Steps 10-14, `phase_state: done` | `close-phase-4`, after the upload | the release workflow run |

### Slice 4D global constraints (in addition to the Phase 4 list)

- No new dependency (runtime, extra, dev or docs). `src/` changes only in 4.F4, 4.F5 and 4.6c
  (`datasets/_biocrust.py`, `_remote.py`, `__init__.py`); the release changes `pyproject.toml`'s
  version only.
- Every pytest, Python and Sphinx run exports **absolute** `BIOTAPY_DATA_DIR` and `HF_HOME` (and,
  for the sdist check, `XDG_CACHE_HOME`); `~/.cache/biotapy` must not appear, nor `docs/pooch/`.
- Delete `docs/_build`, `docs/generated` and `docs/jupyter_execute` before every docs build used
  as a gate.
- The executed tutorials must not need torch, transformers or R; `docs/tutorials/embeddings.md`
  stays without a `jupytext`/`kernelspec` header. No page sets `mystnb` metadata
  (`tests/test_ci.py::test_no_page_overrides_the_notebook_execution_settings`).
- A number in a tutorial's prose is either shown by an executed cell and pinned by a constant its
  tests share, or a data fact a `network` test pins (130 participants, 19 shared samples,
  96 of 996).
- Run `uvx prek run --all-files` after `git add` and before `git commit`; if ruff reformats a
  file, stage it and run prek again (the prototype hit this twice, in `tests/ml/test_mgm.py` and
  `tests/test_coming_from_r.py`; the blocks below are the formatted versions).
- Commits stage explicit paths only (never `.claude/`, `.superpowers/`, `.worktrees/`,
  `notebooks/`, `build/`, `dist/`). Each task's last commit stages
  `.knowledge/roadmap/phase-4-ml-multiomics.md` with its checklist line ticked and
  `.knowledge/log.md` with its line under `## <date> (Phase 4, slice 4D)`, newest first.
- Diff blocks are for reading, not `git apply`; `generated`/`commit` values in them are the
  prototype's (write the model, `date -u +%FT%TZ` and the parent's short hash).
- Branch `phase-4d` from `master` (`28de1ad`) after the plan commit; `release-0.4.0` from
  `master` after `phase-4d` merges; `close-phase-4` after the upload.

### Slice 4D review focus

The five ways a user is most likely to get a wrong answer, or a broken release, from this slice
without an error, most likely first:

1. **A tutorial quotes a number its code no longer gives** (scikit-learn, torch or a dataset
   moves). Expected: the network or `mgm` test that runs the page fails. Tests:
   `tests/ml/test_transformers.py::test_the_leak_free_cv_tutorial_gives_the_numbers_it_quotes`,
   `tests/ml/test_mgm.py::test_the_embedding_tutorial_gives_the_outputs_it_quotes`,
   `tests/tl/test_mmvec.py::test_the_multiomics_tutorial_ranks_mmvec_s_microcoleus_metabolites_above_zero`,
   with the three prose tests beside them.
2. **The leak-free page read as "the prevalence filter leaks a lot", or as a classifier that
   works.** Expected: the page prints the near-zero difference first and says the honest score
   is near chance. Test: `test_the_leak_free_cv_tutorial_quotes_its_numbers` (the sentences and
   the three inequalities behind their words).
3. **A small-integer table summed wrongly.** Expected: `tax_glom`, `func_glom`, `pl.bar` and MGM
   sum `int8`/`int16`/`uint8`/bool tables exactly. Tests: 4.F4
   `test_sum_by_sums_in_the_dtype_numpy_sums_in`,
   `tests/pp/test_glom.py::test_narrow_integer_counts_sum_without_wrapping`; 4.F5
   `test_sum_pairs_sums_in_the_dtype_numpy_sums_in`,
   `tests/fn/test_glom.py::test_narrow_integer_counts_sum_without_wrapping`.
4. **`biocrust()` pairing the wrong samples.** Expected: the 19 shared samples in the microbe
   table's order, no warning. Tests:
   `tests/datasets/test_biocrust.py::test_keeps_the_samples_both_tables_have_in_the_microbe_order_without_a_warning`,
   `test_biocrust_downloads_and_loads`.
5. **The 0.4.0 wheel without an extra or the entry point**, so `pip install 'biotapy[mgm]'` or
   `bt.ml.embed(..., "mgm")` fails after install. Expected: METADATA has `Provides-Extra` `mgm`,
   `r`, `torch` with their `Requires-Dist`, and `entry_points.txt` the `mgm` line. Check: 4.D2
   Step 6 and, after the upload, Step 12.

Execution order: **plan commit -> 4.F4 -> 4.F5 -> 4.6 -> 4.6b -> 4.6c -> 4.D1 -> Checkpoint D
(review, fix pass, re-review, 4.7 knowledge, gates, PR) -> 4.D2.** The fixes come first (the
tutorials use `tax_glom`); 4.6 adds `run_page`, which 4.6b and 4.6c use; 4.6c adds
`datasets.biocrust`, which 4.D1's no-row test names.

---

### Task 4.F4: `_core.sum_by` sums small integer and bool tables without wrapping (`fix(core)`)

**Files:** modify `src/biotapy/_core/_matrix.py` (`sum_by`), `tests/core/test_matrix.py`,
`tests/pp/test_glom.py`, `.knowledge/modules/core.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `sum_pairs` (the same line; Task 4.F5, decision 36); the callers
(`pp/_glom.py`, `pl/_abundance.py`, `ml/_mgm.py`), which need no change; `divide_rows`;
`.knowledge/contracts/data-model-slots.md` (it fixes no dtype for `X`).
**Interfaces:**
- Consumes: NumPy's `sum` accumulator rule, read as `np.zeros(0, dtype=X.dtype).sum().dtype` by
  one private helper, `_sum_dtype(X)`, which 4.F5's `sum_pairs` shares.
- Produces: `sum_by(X, codes, n_groups) -> sp.csr_matrix`, unchanged signature; its dtype is
  `int64` for bool and signed integers, `uint64` for unsigned ones, and `X.dtype` for floats.

- [x] **Step 1: Failing tests.** `_core` unit test (R4.9's exception) and one through the public
  API:
````diff
diff --git a/tests/core/test_matrix.py b/tests/core/test_matrix.py
--- a/tests/core/test_matrix.py
+++ b/tests/core/test_matrix.py
@@ -30,6 +30,25 @@ def test_sum_by_keeps_integer_dtype():
     assert sum_by(X, np.array([0, 0, 0]), 1).dtype == np.int64


+@pytest.mark.parametrize(
+    ("dtype", "value", "summed"),
+    [
+        (np.int8, 100, np.int64),
+        (np.uint8, 200, np.uint64),
+        (np.int16, 20_000, np.int64),
+        (np.int32, 1_000_000_000, np.int64),
+        (np.uint32, 2_000_000_000, np.uint64),
+        (np.bool_, True, np.int64),
+        (np.float32, 0.5, np.float32),
+    ],
+)
+def test_sum_by_sums_in_the_dtype_numpy_sums_in(dtype, value, summed):
+    # Three columns into one group: 3 * value overflows int8, uint8 and int16, and bool would saturate at True.
+    out = sum_by(sp.csr_matrix(np.full((2, 3), value, dtype=dtype)), np.array([0, 0, 0]), 1)
+    np.testing.assert_array_equal(out.toarray(), np.full((2, 1), 3 * value))
+    assert out.dtype == summed
+
+
 def test_argmax_by_one_index_per_group_in_code_order():
     np.testing.assert_array_equal(argmax_by(np.array([1.0, 9.0, 3.0, 7.0]), np.array([1, 0, 1, 0])), [1, 2])

diff --git a/tests/pp/test_glom.py b/tests/pp/test_glom.py
--- a/tests/pp/test_glom.py
+++ b/tests/pp/test_glom.py
@@ -66,6 +66,15 @@ def test_genus_sums_members_and_drops_unassigned():
     np.testing.assert_array_equal(out[:, "f4"].X.toarray().ravel(), [30, 25, 36, 10, 14, 9])


+def test_narrow_integer_counts_sum_without_wrapping():
+    # int8 holds up to 127: summing three features of 100 into one genus must give 300, not 44.
+    var = pd.DataFrame({"kingdom": ["K"] * 3, "genus": ["G"] * 3}, index=["a", "b", "c"])
+    adata = AnnData(X=sp.csr_matrix(np.full((2, 3), 100, dtype=np.int8)), var=var, obs=pd.DataFrame(index=["s1", "s2"]))
+    out = bt.pp.tax_glom(adata, "genus")
+    np.testing.assert_array_equal(out.X.toarray(), [[300], [300]])
+    assert out.X.dtype == np.int64
+
+
 def test_dropna_false_keeps_unassigned_as_own_group():
     assert bt.pp.tax_glom(bt.datasets.toy(), "genus", dropna=False).n_vars == 7

````
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_matrix.py
  tests/pp/test_glom.py -q` -> `7 failed, 42 passed`: the `int8`, `uint8`, `int16`, `int32`, `uint32`
  and `bool` cases give `44 (ACTUAL), 300 (DESIRED)`, `88 (ACTUAL), 600 (DESIRED)`, `-5536 (ACTUAL), 60000
  (DESIRED)`, `-1294967296 (ACTUAL), 3000000000 (DESIRED)`, `1705032704 (ACTUAL), 6000000000
  (DESIRED)` and `True (ACTUAL), 3 (DESIRED)`, and `test_narrow_integer_counts_sum_without_wrapping`
  `44 (ACTUAL), 300 (DESIRED)` with `dtype=int8`. The `float32` case passes already.
- [x] **Step 3: Implement.**
````diff
diff --git a/src/biotapy/_core/_matrix.py b/src/biotapy/_core/_matrix.py
--- a/src/biotapy/_core/_matrix.py
+++ b/src/biotapy/_core/_matrix.py
@@ -34,11 +34,22 @@ def divide_rows(X: sp.csr_matrix, totals: npt.NDArray[np.float64]) -> sp.csr_mat
     return out


+def _sum_dtype(X: sp.csr_matrix) -> np.dtype[np.generic]:
+    """The dtype ``numpy.sum`` gives ``X``'s: bool and narrower integers widen to 64 bits, floats keep theirs."""
+    return cast("np.dtype[np.generic]", np.zeros(0, dtype=X.dtype).sum().dtype)
+
+
 def sum_by(X: sp.csr_matrix, codes: npt.NDArray[np.intp], n_groups: int) -> sp.csr_matrix:
-    """Sum the columns of ``X`` that share a group code; negative codes are dropped."""
+    """Sum the columns of ``X`` that share a group code; negative codes are dropped.
+
+    The sums have the dtype ``numpy.sum`` gives ``X``'s: bool and integers
+    narrower than 64 bits widen to 64 bits, floats keep theirs.
+    """
     rows = np.flatnonzero(codes >= 0)
+    # SciPy multiplies in the wider of the two dtypes; an int8 or bool indicator would wrap or saturate the sums.
+    dtype = _sum_dtype(X)
     indicator = sp.csr_matrix(
-        (np.ones(rows.size, dtype=X.dtype), (rows, codes[rows])),
+        (np.ones(rows.size, dtype=dtype), (rows, codes[rows])),
         shape=(codes.size, n_groups),
     )
     return sp.csr_matrix(X @ indicator)
````
  R2.1: NumPy's own accumulator rule, no branch of biotapy's; SciPy then multiplies in the wider
  dtype. The helper holds the expression because 4.F5 uses it a second time (R4.3).
- [x] **Step 4: Run, expect pass** - the same command -> `49 passed`.
- [x] **Step 5: Knowledge.** `core.md` said "today's only caller is `pp.tax_glom`", false since
  slice 4A (`pl.bar`) and 4C (MGM); fixed here (R0.2). Tick 4.F4; the log section opens:
````diff
diff --git a/.knowledge/log.md b/.knowledge/log.md
--- a/.knowledge/log.md
+++ b/.knowledge/log.md
@@ -1,5 +1,8 @@
 # Knowledge bundle log

+## 2026-10-10 (Phase 4, slice 4D)
+- **Update**: [core](modules/core.md): `sum_by` sums in `numpy.sum`'s dtype for `X` (an int8 or bool table no longer wraps or saturates) and names its three callers (`pp.tax_glom`, `pl.bar`'s `fill`, MGM), not "today's only caller"; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.F4.
+
 ## 2026-10-10 (Phase 4, slice 4D plan)
 - **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) expands slice 4D (4.F4 and 4.F5 sums, 4.6 leak-free tutorial, 4.6b embedding page, 4.6c `datasets.biocrust` and the multi-omics tutorial, 4.D1 Coming-from-R, Checkpoint D with 4.7 knowledge, 4.D2 release 0.4.0) into full TDD steps, prototyped and gated per commit, and records decisions 36-45; ticks Checkpoint C (PR #30 merged, slice 4C approved), exit-gate item 2 (PR #30's `ml-extras`) and 4.F1's upstream box (scikit-bio#2631); [da](modules/da.md)'s ANCOM-BC2 gotcha links the issue.

diff --git a/.knowledge/modules/core.md b/.knowledge/modules/core.md
--- a/.knowledge/modules/core.md
+++ b/.knowledge/modules/core.md
@@ -5,8 +5,8 @@ description: Private kernel package - sparse group math, taxonomic rank order, f
 resource: /src/biotapy/_core/
 paths: ["src/biotapy/_core/**"]
 tags: [core, kernel]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T16:44:46Z }
-commit: 31aa11d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T16:40:30Z }
+commit: b5654af
 status: stable
 ---

@@ -39,8 +39,13 @@ none of them back.
   `pp.relative`, `fn.renorm` (`_renorm.py:_rescaled`) and
   `fn.functional_redundancy`.
 - `_matrix.py:sum_by` / `_matrix.py:argmax_by` - grouped column sum and
-  grouped argmax by integer group codes, negative codes dropped; today's only
-  caller is `pp.tax_glom`.
+  grouped argmax by integer group codes, negative codes dropped. `sum_by`
+  sums in the dtype `numpy.sum` gives `X`'s (bool and every integer narrower
+  than 64 bits sum to int64, or uint64 if unsigned; floats keep theirs), so an
+  int8 or int32 table cannot wrap. The callers of `sum_by`
+  are `pp.tax_glom`, `pl.bar`'s `fill` by a `var` column
+  (`pl/_abundance.py:_segments`) and MGM (`ml/_mgm.py`); `argmax_by`'s is
+  `pp.tax_glom`.
 - `_matrix.py:sum_pairs` - column sums into groups from `(feature, group)`
   membership pairs, where a feature may sit in several groups and counts in
   full toward each (many-to-many); the pairs are a set, a repeated pair counts
````
- [x] **Step 6: Gate and commit.**
  ```bash
  git add src/biotapy/_core/_matrix.py tests/core/test_matrix.py tests/pp/test_glom.py \
    .knowledge/modules/core.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  git commit -m "fix(core): sum int8, int16 and bool tables in sum_by without wrapping

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Then the full gate, with both extras (MGM calls `sum_by`). Expected: prek passed; `1496
  passed, 2 skipped, 92 deselected`; `36 passed, 1554 deselected`; `build succeeded`; `25 passed,
  1565 deselected`; `13 passed, 1577 deselected`.

### Task 4.F5: `_core.sum_pairs`, the same fix (`fix(core)`, decision 36)

Drop this task if the user declines decision 36; nothing later depends on it.

**Files:** modify `src/biotapy/_core/_matrix.py` (`sum_pairs`), `tests/core/test_matrix.py`,
`tests/fn/test_glom.py`, `.knowledge/modules/core.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `fn/_glom.py` (its `agg="mean"` path multiplies by a float64 diagonal after the
sum, so it was already float); `sum_by`'s body (4.F4 already holds the helper `_sum_dtype`).
**Interfaces:**
- Consumes: 4.F4's rule and its helper `_sum_dtype(X)`.
- Produces: `sum_pairs(X, features, groups, *, n_groups)`, unchanged signature, same dtype rule.

- [x] **Step 1: Failing tests.**
````diff
diff --git a/tests/core/test_matrix.py b/tests/core/test_matrix.py
--- a/tests/core/test_matrix.py
+++ b/tests/core/test_matrix.py
@@ -92,6 +92,25 @@ def test_sum_pairs_keeps_integer_dtype():
     assert sum_pairs(X, np.array([0]), np.array([0]), n_groups=1).dtype == np.int64


+@pytest.mark.parametrize(
+    ("dtype", "value", "summed"),
+    [
+        (np.int8, 100, np.int64),
+        (np.uint8, 200, np.uint64),
+        (np.int16, 20_000, np.int64),
+        (np.int32, 1_000_000_000, np.int64),
+        (np.uint32, 2_000_000_000, np.uint64),
+        (np.bool_, True, np.int64),
+        (np.float32, 0.5, np.float32),
+    ],
+)
+def test_sum_pairs_sums_in_the_dtype_numpy_sums_in(dtype, value, summed):
+    dense = np.full((2, 3), value, dtype=dtype)
+    out = sum_pairs(sp.csr_matrix(dense), np.array([0, 1, 2]), np.array([0, 0, 0]), n_groups=1)
+    np.testing.assert_array_equal(out.toarray(), np.full((2, 1), 3 * value))
+    assert out.dtype == summed
+
+
 @given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 50)), st.data())
 def test_sum_pairs_with_one_group_per_feature_equals_sum_by(dense, data):
     codes = np.array(data.draw(st.lists(st.integers(-1, 2), min_size=dense.shape[1], max_size=dense.shape[1])))
diff --git a/tests/fn/test_glom.py b/tests/fn/test_glom.py
--- a/tests/fn/test_glom.py
+++ b/tests/fn/test_glom.py
@@ -113,6 +113,14 @@ def test_counts_keep_their_label_when_each_feature_has_one_parent(make_adata):
     assert bt.fn.func_glom(_counts(make_adata), "p", hierarchy=hierarchy).uns["biotapy"]["x_kind"] == "counts"


+def test_narrow_integer_counts_sum_without_wrapping(make_adata):
+    # int8 holds up to 127: two features of 100 under one parent must give 200, not -56.
+    hierarchy = pd.DataFrame({"child": ["f0", "f1"], "parent": "P1", "level": "p"})
+    out = bt.fn.func_glom(make_adata(np.full((2, 2), 100, dtype=np.int8)), "p", hierarchy=hierarchy)
+    np.testing.assert_array_equal(out.X.toarray(), [[200], [200]])
+    assert out.X.dtype == np.int64
+
+
 def test_a_mean_is_labelled_abundance(make_adata):
     hierarchy = pd.DataFrame({"child": ["f0", "f1"], "parent": "P1", "level": "p"})
     out = bt.fn.func_glom(_counts(make_adata), "p", hierarchy=hierarchy, agg="mean")
````
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_matrix.py
  tests/fn/test_glom.py -q` -> `7 failed, 63 passed` (44, 88, -5536, -1294967296, 1705032704 and `True` as in 4.F4;
  `func_glom` gives `-56 (ACTUAL), 200 (DESIRED)`).
- [x] **Step 3: Implement.**
````diff
diff --git a/src/biotapy/_core/_matrix.py b/src/biotapy/_core/_matrix.py
--- a/src/biotapy/_core/_matrix.py
+++ b/src/biotapy/_core/_matrix.py
@@ -72,10 +72,13 @@ def sum_pairs(

     A feature listed with several groups counts in full toward each of them
     (many-to-many, as ``humann_regroup_table`` does). The pairs are a set: a
-    repeated pair counts once. A feature in no pair is left out.
+    repeated pair counts once. A feature in no pair is left out. The sums
+    have the dtype ``numpy.sum`` gives ``X``'s, as in ``sum_by``.
     """
+    # As in sum_by: an int8 or bool indicator would wrap or saturate the sums.
+    dtype = _sum_dtype(X)
     indicator = sp.csr_matrix(
-        (np.ones(features.size, dtype=X.dtype), (features, groups)),
+        (np.ones(features.size, dtype=dtype), (features, groups)),
         shape=(X.shape[1], n_groups),
     )
     # The constructor sums repeated (feature, group) entries; membership is a set, so reset them to 1.
````
- [x] **Step 4: Run, expect pass** - `70 passed`.
- [x] **Step 5: Knowledge.** Tick 4.F5; log line first under the slice heading:
````diff
diff --git a/.knowledge/log.md b/.knowledge/log.md
--- a/.knowledge/log.md
+++ b/.knowledge/log.md
@@ -1,6 +1,7 @@
 # Knowledge bundle log

 ## 2026-10-10 (Phase 4, slice 4D)
+- **Update**: [core](modules/core.md): `sum_pairs` sums in `numpy.sum`'s dtype for `X`, as `sum_by` does, so `fn.func_glom` cannot wrap an int8 table; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.F5.
 - **Update**: [core](modules/core.md): `sum_by` sums in `numpy.sum`'s dtype for `X` (an int8 or bool table no longer wraps or saturates) and names its three callers (`pp.tax_glom`, `pl.bar`'s `fill`, MGM), not "today's only caller"; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.F4.

 ## 2026-10-10 (Phase 4, slice 4D plan)
diff --git a/.knowledge/modules/core.md b/.knowledge/modules/core.md
--- a/.knowledge/modules/core.md
+++ b/.knowledge/modules/core.md
@@ -5,8 +5,8 @@ description: Private kernel package - sparse group math, taxonomic rank order, f
 resource: /src/biotapy/_core/
 paths: ["src/biotapy/_core/**"]
 tags: [core, kernel]
-generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T16:40:30Z }
-commit: b5654af
+generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-10T18:21:05Z }
+commit: 024d532
 status: stable
 ---

@@ -42,7 +42,8 @@ none of them back.
   grouped argmax by integer group codes, negative codes dropped. `sum_by`
   sums in the dtype `numpy.sum` gives `X`'s (bool and every integer narrower
   than 64 bits sum to int64, or uint64 if unsigned; floats keep theirs), so an
-  int8 or int32 table cannot wrap. The callers of `sum_by`
+  int8 or int32 table cannot wrap (`sum_by` and `sum_pairs` take the dtype
+  from one private helper, `_matrix.py:_sum_dtype`). The callers of `sum_by`
   are `pp.tax_glom`, `pl.bar`'s `fill` by a `var` column
   (`pl/_abundance.py:_segments`) and MGM (`ml/_mgm.py`); `argmax_by`'s is
   `pp.tax_glom`.
@@ -50,7 +51,7 @@ none of them back.
   membership pairs, where a feature may sit in several groups and counts in
   full toward each (many-to-many); the pairs are a set, a repeated pair counts
   once; used by `fn.func_glom`. Where `sum_by` assigns each feature one group,
-  this does not.
+  this does not. It sums in `numpy.sum`'s dtype for `X`, as `sum_by` does.
 - `_composition.py:pseudocounted` - validates `pseudocount`
   (`_composition.py:check_pseudocount`), rejects negative or non-finite `X`,
   warns when the pseudocount exceeds the smallest non-zero value (the message ends
````
- [x] **Step 6: Gate and commit.**
  ```bash
  git add src/biotapy/_core/_matrix.py tests/core/test_matrix.py tests/fn/test_glom.py \
    .knowledge/modules/core.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  git commit -m "fix(core): sum int8, int16 and bool tables in sum_pairs without wrapping

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1504 passed, 2 skipped, 92 deselected`; `36 passed, 1562
  deselected`; `build succeeded`.

### Task 4.6: The leak-free cross-validation tutorial (exit gate 1)

**Files:** create `docs/tutorials/leak_free_cv.md` and `tests/test_run_page.py`; modify `tests/conftest.py`,
`tests/ml/test_transformers.py`, `docs/tutorials/index.md`, `docs/guide/machine_learning.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `docs/conf.py` (the cache mode and the 300 s cell timeout hold), the workflows
and `.readthedocs.yaml` (the `docs` job's pooch cache already holds HMP2; Read the Docs
downloads it for the function tutorial anyway), `src/` (no API change), `.knowledge/modules/ml.md`
(4.7).
**Interfaces:**
- Consumes: `bt.datasets.hmp2()`, `bt.pp.filter_features(adata, *, min_prevalence)`,
  `bt.ml.PrevalenceFilter(*, min_prevalence)`, `bt.ml.CLR(*, pseudocount)`; scikit-learn's
  `make_pipeline`, `RepeatedStratifiedKFold`, `cross_val_score`, `SelectKBest`, `f_classif`,
  `LogisticRegression`.
- Produces: the page `tutorials/leak_free_cv` (exit gate 1); the fixture `run_page(path: Path)
  -> dict[str, object]` in `tests/conftest.py` (4.6b and 4.6c use it; it skips a cell's MyST
  options and raises a `ValueError` on IPython magics and shell escapes); `LEAK_FREE_CV`,
  `LEAK_FREE_CV_AUC` and `LEAK_FREE_CV_SHUFFLED_RANGE` in `tests/ml/test_transformers.py`; the
  page's names `inside`, `outside`, `select_inside`, `select_outside`, `shuffled_inside`,
  `shuffled_outside` (unrounded floats; the page rounds where it shows them) and `shuffled_scores`.

- [x] **Step 1: Failing tests.** The fixture (with its own tests in `tests/test_run_page.py`),
  the constants, a `network` test that runs the page and a default-run test of its prose:
````diff
diff --git a/tests/conftest.py b/tests/conftest.py
--- a/tests/conftest.py
+++ b/tests/conftest.py
@@ -1,4 +1,6 @@
+import re
 from collections.abc import Callable
+from pathlib import Path

 import matplotlib
 import numpy as np
@@ -51,3 +53,31 @@ def _make_adata(dense: np.ndarray) -> AnnData:
 def make_adata() -> Callable[[np.ndarray], AnnData]:
     """AnnData from a dense samples x features array, with samples ``s0..`` and features ``f0..``."""
     return _make_adata
+
+
+# A MyST notebook's code cells and a page's fenced Python blocks, in page order.
+_PAGE_CODE = re.compile(r"^```(?:\{code-cell\} ipython3|python)\n(.*?)^```$", re.MULTILINE | re.DOTALL)
+
+
+# A cell's MyST options come first: ":key: value" lines, or a YAML block between two "---" lines.
+_CELL_OPTIONS = re.compile(r"\A(?:---\n.*?^---\n|(?::[\w-]+:[^\n]*\n)+)", re.MULTILINE | re.DOTALL)
+_MAGIC = re.compile(r"^\s*[%!]", re.MULTILINE)
+
+
+def _cell_source(cell: str) -> str:
+    source = _CELL_OPTIONS.sub("", cell, count=1)
+    if _MAGIC.search(source):
+        raise ValueError(f"run_page cannot run an IPython magic or shell escape:\n{source}")
+    return source
+
+
+def _run_page(path: Path) -> dict[str, object]:
+    namespace: dict[str, object] = {}
+    exec("\n".join(_cell_source(cell) for cell in _PAGE_CODE.findall(path.read_text(encoding="utf-8"))), namespace)
+    return namespace
+
+
+@pytest.fixture(scope="session")
+def run_page() -> Callable[[Path], dict[str, object]]:
+    """Run a docs page's code, as a reader would, and return the names it defines (tests quoting a page's numbers)."""
+    return _run_page
diff --git a/tests/ml/test_transformers.py b/tests/ml/test_transformers.py
--- a/tests/ml/test_transformers.py
+++ b/tests/ml/test_transformers.py
@@ -1,3 +1,5 @@
+from pathlib import Path
+
 import numpy as np
 import pandas as pd
 import pytest
@@ -202,3 +204,63 @@ def test_clr_rows_sum_to_zero(X):
 def test_transformer_options_are_keyword_only(estimator, value):
     with pytest.raises(TypeError):
         estimator(value)
+
+
+LEAK_FREE_CV = Path(__file__).parents[2] / "docs" / "tutorials" / "leak_free_cv.md"
+# The mean ROC AUCs of the leak-free cross-validation tutorial on HMP2 (Phase 4 exit gate 1), which its prose quotes.
+# The network test runs the page's cells and checks they still give these; the next test checks the prose quotes them.
+# The page keeps the AUCs unrounded and shows three decimals.
+LEAK_FREE_CV_AUC = {
+    "inside": 0.527,
+    "outside": 0.526,
+    "select_inside": 0.558,
+    "select_outside": 0.705,
+    "shuffled_inside": 0.544,
+    "shuffled_outside": 0.767,
+}
+# The range of the ten shuffled-label scores of the pipeline, which the prose quotes.
+LEAK_FREE_CV_SHUFFLED_RANGE = (0.484, 0.669)
+
+
+@pytest.mark.network
+def test_the_leak_free_cv_tutorial_gives_the_numbers_it_quotes(run_page):
+    namespace = run_page(LEAK_FREE_CV)
+    for name, quoted in LEAK_FREE_CV_AUC.items():
+        # One ranking flip moves a mean over 25 folds by about 3.6e-4, so the 3-decimal figures can slip by a unit.
+        assert namespace[name] == pytest.approx(quoted, abs=1e-3), name
+    # "Near a coin toss" in the prose.
+    assert abs(namespace["shuffled_inside"] - 0.5) < 0.05
+    assert namespace["inside"] < 0.6
+    low, high = LEAK_FREE_CV_SHUFFLED_RANGE
+    assert min(namespace["shuffled_scores"]) == pytest.approx(low, abs=1e-3)
+    assert max(namespace["shuffled_scores"]) == pytest.approx(high, abs=1e-3)
+    # "The honest scores fall inside that range."
+    assert (
+        min(namespace["shuffled_scores"])
+        < namespace["inside"]
+        < namespace["select_inside"]
+        < max(namespace["shuffled_scores"])
+    )
+
+
+def test_the_leak_free_cv_tutorial_quotes_its_numbers():
+    page = " ".join(LEAK_FREE_CV.read_text(encoding="utf-8").split())
+    auc = LEAK_FREE_CV_AUC
+    for sentence in [
+        f"The prevalence filter fitted inside the pipeline scores {auc['inside']} and the one fitted on every sample "
+        f"{auc['outside']}: no measurable leak.",
+        "here, fitting it on the test samples moves the score by 0.001",
+        f"Selected inside the pipeline, the score is {auc['select_inside']}; selected on every sample, "
+        f"{auc['select_outside']}.",
+        f"On shuffled labels the pipeline scores {auc['shuffled_inside']} on average, near a coin toss",
+        f"The ten shuffled-label scores of the pipeline range from {LEAK_FREE_CV_SHUFFLED_RANGE[0]} to "
+        f"{LEAK_FREE_CV_SHUFFLED_RANGE[1]}, and the honest scores of {auc['inside']} and {auc['select_inside']} fall "
+        "inside that range",
+        f"The version that selects species on every sample scores {auc['shuffled_outside']}: higher than it scored on "
+        "the real labels.",
+    ]:
+        assert sentence in page
+    # What the prose says about the numbers: no measurable leak, a large one, and a larger one on noise.
+    assert abs(auc["inside"] - auc["outside"]) < 0.01
+    assert auc["select_outside"] - auc["select_inside"] > 0.1
+    assert auc["shuffled_outside"] > auc["select_outside"]
diff --git a/tests/test_run_page.py b/tests/test_run_page.py
new file mode 100644
--- /dev/null
+++ b/tests/test_run_page.py
@@ -0,0 +1,30 @@
+import pytest
+
+PAGE = """# A page
+
+```{{code-cell}} ipython3
+{options}x = 1
+```
+
+```python
+y = x + 1
+```
+"""
+
+
+@pytest.mark.parametrize(
+    "options",
+    ["", ":tags: [hide-input]\n:mystnb-key: value\n", "---\ntags: [hide-input]\nraises-exception: false\n---\n"],
+)
+def test_run_page_skips_myst_cell_options(run_page, tmp_path, options):
+    page = tmp_path / "page.md"
+    page.write_text(PAGE.format(options=options), encoding="utf-8")
+    assert {name: run_page(page)[name] for name in ("x", "y")} == {"x": 1, "y": 2}
+
+
+@pytest.mark.parametrize("magic", ["%matplotlib inline", "!pip list"])
+def test_run_page_fails_clearly_on_ipython_magics(run_page, tmp_path, magic):
+    page = tmp_path / "page.md"
+    page.write_text(PAGE.format(options=f"{magic}\n"), encoding="utf-8")
+    with pytest.raises(ValueError, match="IPython magic"):
+        run_page(page)
````
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/ml/test_transformers.py
  -q -m "network or not network" -k leak_free` -> `2 failed, 125 deselected`, both
  `FileNotFoundError: ... docs/tutorials/leak_free_cv.md`.
- [x] **Step 3: The tutorial.** Create `docs/tutorials/leak_free_cv.md`:
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

# Leak-free cross-validation

Cross-validation estimates how well a model does on samples it has not seen. A preprocessing
step that learns from the samples, such as choosing which features to keep, leaks when it is
fitted on every sample before the split: the test samples then shape what the model is trained
on, and the score comes out better than the model will do on new samples. This tutorial measures
that on a real cohort, with `bt.ml`'s transformers inside a scikit-learn pipeline, so that every
step is fitted on the training samples of each fold only.

:::{note}
The HMP2 tables are downloaded from the [IBDMDB](https://ibdmdb.org/) on first use (23 MB) and
cached. The IBDMDB states no licence for them; biotapy ships none of them. Cite the study when
you use them: Lloyd-Price J et al. (2019) Multi-omics of the gut microbial ecosystem in
inflammatory bowel diseases. *Nature* 569:655-662.
:::

```{code-cell} ipython3
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline

import biotapy as bt
```

## The question

Does the species profile of a participant's first stool sample tell inflammatory bowel disease
(Crohn's disease or ulcerative colitis) from non-IBD? `bt.datasets.hmp2()` holds the HMP2
cohort's MetaPhlAn 3 species for 130 participants.

```{code-cell} ipython3
taxa = bt.datasets.hmp2()["taxa"]
taxa.obs["diagnosis"].value_counts()
```

```{code-cell} ipython3
ibd = (taxa.obs["diagnosis"] != "nonIBD").to_numpy()
```

MetaPhlAn's values are relative abundances: each sample sums to 1. `bt.ml.CLR`'s default
pseudocount of 0.5 would swamp them, so it is set on their scale, below the smallest non-zero
value:

```{code-cell} ipython3
f"{taxa.X.data.min():.2g}"
```

## The pipeline

The prevalence filter keeps the species present in at least 10% of the samples, the centred
log-ratio puts them on a log scale, and a logistic regression predicts IBD. Five-fold
cross-validation, stratified so every fold keeps the share of non-IBD participants, is repeated
five times with different splits, and the score is the mean ROC AUC over the 25 test folds:
0.5 is a coin toss, 1 a perfect ranking.

```{code-cell} ipython3
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=0)


def mean_auc(model, X, y):
    return float(cross_val_score(model, X, y, cv=cv, scoring="roc_auc").mean())


pipeline = make_pipeline(
    bt.ml.PrevalenceFilter(min_prevalence=0.1),
    bt.ml.CLR(pseudocount=1e-6),
    LogisticRegression(max_iter=5000),
)
inside = mean_auc(pipeline, taxa.X, ibd)
round(inside, 3)
```

Close to a coin toss: with this model, a first stool sample's species hardly separate IBD from
non-IBD participants in this cohort. That is the honest answer, and the one the leaky versions
below would have hidden.

## The prevalence filter outside the pipeline

The common mistake is to filter the whole table first, with `bt.pp.filter_features`, and
cross-validate only the model:

```{code-cell} ipython3
filtered = bt.pp.filter_features(taxa, min_prevalence=0.1)
outside = mean_auc(make_pipeline(bt.ml.CLR(pseudocount=1e-6), LogisticRegression(max_iter=5000)), filtered.X, ibd)
round(outside, 3)
```

The prevalence filter fitted inside the pipeline scores 0.527 and the one fitted on every
sample 0.526: no measurable leak. The filter never looks at the labels, so here, fitting it on the
test samples moves the score by 0.001; it only changes which rare species the model gets. A step
fitted on all samples without the labels can still bias cross-validation in general, so keep it
inside the pipeline all the same: it costs nothing, and then no step is fitted on a test sample
without anyone having to argue that this one is harmless.

## A step that reads the labels

Choosing features by how well they separate the groups is different. `SelectKBest` keeps the 20
species whose CLR values differ most between IBD and non-IBD (an F-test). Inside the pipeline,
it chooses them from the training samples of each fold; outside, from all 130 samples, test
samples included:

```{code-cell} ipython3
selecting = make_pipeline(
    bt.ml.PrevalenceFilter(min_prevalence=0.1),
    bt.ml.CLR(pseudocount=1e-6),
    SelectKBest(f_classif, k=20),
    LogisticRegression(max_iter=5000),
)
select_inside = mean_auc(selecting, taxa.X, ibd)
clr = bt.ml.CLR(pseudocount=1e-6).fit_transform(filtered.X)
select_outside = mean_auc(LogisticRegression(max_iter=5000), SelectKBest(f_classif, k=20).fit_transform(clr, ibd), ibd)
round(select_inside, 3), round(select_outside, 3)
```

Selected inside the pipeline, the score is 0.558; selected on every sample, 0.705. The second
number describes a classifier that does not exist: its 20 species were chosen because they
separate the very samples it is then tested on.

## Labels that mean nothing

The plainest check is to shuffle the labels, so that no species can predict them, and run both
versions again. Ten shuffles, each cross-validated as above:

```{code-cell} ipython3
rng = np.random.default_rng(0)
shuffles = [rng.permutation(ibd) for _ in range(10)]
shuffled_scores = [mean_auc(selecting, taxa.X, labels) for labels in shuffles]
shuffled_inside = float(np.mean(shuffled_scores))
shuffled_outside = float(
    np.mean(
        [
            mean_auc(LogisticRegression(max_iter=5000), SelectKBest(f_classif, k=20).fit_transform(clr, labels), labels)
            for labels in shuffles
        ]
    )
)
round(shuffled_inside, 3), round(shuffled_outside, 3), round(min(shuffled_scores), 3), round(max(shuffled_scores), 3)
```

On shuffled labels the pipeline scores 0.544 on average, near a coin toss, as it should. The
version that selects species on every sample scores 0.767: higher than it scored on the real
labels. The leak manufactures a signal from noise.

The ten shuffled-label scores of the pipeline range from 0.484 to 0.669, and the honest scores of
0.527 and 0.558 fall inside that range: on this cohort they cannot be told apart from scores on
shuffled labels.

```{code-cell} ipython3
pd.DataFrame(
    {"inside the pipeline": [inside, select_inside, shuffled_inside], "on every sample": [outside, select_outside, shuffled_outside]},
    index=["prevalence filter", "SelectKBest, real labels", "SelectKBest, shuffled labels"],
).round(3)
```

## What goes inside the pipeline

Every step that learns from more than one sample: which features to keep, by prevalence or by
their relation to the labels, and how to scale them. Steps that transform each sample on its own,
such as relative abundance and CLR, give the same values either way; putting them inside keeps the
rule simple. The [machine learning guide](../guide/machine_learning.md#which-steps-leak) lists
which `bt.ml` and scikit-learn steps learn from the samples.
````
  Link it from the tutorials toctree and the guide's "Which steps leak":
````diff
diff --git a/docs/guide/machine_learning.md b/docs/guide/machine_learning.md
--- a/docs/guide/machine_learning.md
+++ b/docs/guide/machine_learning.md
@@ -20,7 +20,10 @@ and cross-validated scores come out better than they will be on new samples.

 `bt.pp.filter_features` on the whole table before cross-validation is the
 leaky version of `bt.ml.PrevalenceFilter`: the same rule, fitted on every
-sample at once.
+sample at once. The {doc}`leak-free cross-validation tutorial
+</tutorials/leak_free_cv>` measures both versions on the HMP2 cohort: there
+the prevalence filter leaks almost nothing, while a step that chooses features
+by the labels inflates the score even on shuffled labels.

 ## A leak-free pipeline

diff --git a/docs/tutorials/index.md b/docs/tutorials/index.md
--- a/docs/tutorials/index.md
+++ b/docs/tutorials/index.md
@@ -8,4 +8,5 @@ quick_tour
 phyloseq_analysis
 function
 differential_abundance
+leak_free_cv
 ```
````
- [x] **Step 4: Run, expect pass** - the Step 2 command with `-W error::UserWarning` -> `2
  passed, 125 deselected` (about 11 s; the first run downloads HMP2, 23 MB).
- [x] **Step 5: Build and read the page.** `rm -rf docs/_build docs/generated
  docs/jupyter_execute && uv run --group doc sphinx-build -W -b html docs docs/_build/html` ->
  `build succeeded.` with `tutorials/leak_free_cv.md: Executed notebook in <n> seconds`
  (11.2-13.7 s on the prototype). In `docs/_build/html/tutorials/leak_free_cv.html` the outputs
  read, in order: `diagnosis` counts `CD 65`, `UC 38`, `nonIBD 27`; `'1.5e-06'`; `0.527`;
  `0.526`; `(0.558, 0.705)`; `(0.544, 0.767, 0.484, 0.669)`; and the table (prevalence filter 0.527 / 0.526,
  real labels 0.558 / 0.705, shuffled 0.544 / 0.767). No stderr block.
- [x] **Step 6: Bookkeeping.** Tick 4.6 (the exit-gate box waits for the CI docs job, 4.D2);
  log line:
  `- **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.6 done: \`docs/tutorials/leak_free_cv.md\` cross-validates a prevalence filter, CLR and logistic regression on HMP2's species on every docs build, inside and outside the pipeline, with a label-reading \`SelectKBest\` on real and shuffled labels; a network test runs the page's cells against the numbers its prose quotes (\`LEAK_FREE_CV_AUC\`).`
- [x] **Step 7: Gate and commit.**
  ```bash
  git add docs/tutorials/leak_free_cv.md docs/tutorials/index.md docs/guide/machine_learning.md \
    tests/conftest.py tests/test_run_page.py tests/ml/test_transformers.py \
    .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  git commit -m "docs: add the leak-free cross-validation tutorial on HMP2

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1510 passed, 2 skipped, 93 deselected`; `37 passed, 1568
  deselected`; `build succeeded` with six executed notebooks.

### Task 4.6b: The MGM embedding tutorial (exit gate 2's page)

**Files:** create `docs/tutorials/embeddings.md`; modify `tests/ml/test_mgm.py`,
`docs/tutorials/index.md`, `docs/guide/machine_learning.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `src/biotapy/ml/` (no code change); the `ml-extras` job (its `-m mgm` step
already runs every `mgm` test); `docs/conf.py`; the guide's MGM, licence and plugin sections
(the page links them).
**Interfaces:**
- Consumes: `run_page` (4.6); `bt.pp.tax_glom`, `bt.datasets.global_patterns`, `bt.ml.embed(adata,
  model, *, inplace=False)`; scikit-learn's `NearestNeighbors`.
- Produces: the page `tutorials/embeddings`; `EMBEDDINGS`, `GLOBAL_PATTERNS_LEFT_OUT`,
  `GLOBAL_PATTERNS_SHAPE`, `MGM_GENERA` in `tests/ml/test_mgm.py`; the page's names `genera`, `embedding`,
  `same_type`.

- [x] **Step 1: Failing tests.** The exit-gate test matches the whole warning through the new
  constant; an `mgm` test runs the page; a default-run test checks what it quotes:
````diff
diff --git a/tests/ml/test_mgm.py b/tests/ml/test_mgm.py
--- a/tests/ml/test_mgm.py
+++ b/tests/ml/test_mgm.py
@@ -1,3 +1,4 @@
+import re
 import sys
 from importlib.metadata import entry_points
 from pathlib import Path
@@ -12,6 +13,16 @@ from sklearn.neighbors import NearestNeighbors
 import biotapy as bt

 DATA = Path(__file__).parents[1] / "data" / "mgm"
+EMBEDDINGS = Path(__file__).parents[2] / "docs" / "tutorials" / "embeddings.md"
+# Phase 4 exit gate 2 on GlobalPatterns' genera. docs/tutorials/embeddings.md quotes these outputs; its build has no torch,
+# so an mgm test runs the page's code and a default-run test checks the page quotes them.
+GLOBAL_PATTERNS_LEFT_OUT = (
+    "mgm leaves out 96 of 996 features: 0 without a genus and 96 whose genus is not one of MGM's "
+    "(4-29, 4041AA30, A17, Aquamonas, Arctic95A-2, ...)"
+)
+GLOBAL_PATTERNS_SHAPE = (26, 256)
+# The size of MGM's vocabulary (phylogeny.csv's genera), which the tutorial quotes.
+MGM_GENERA = 9665
 # torch 2.13-2.14's first tanh in a process can saturate on CPUs running more than 4 threads (measured: the
 # embedding's first call off by up to 1.6e-4, later calls by 1.7e-6), so equality is checked to 1e-3.
 ATOL = 1e-3
@@ -133,11 +144,11 @@ def test_keeps_the_input(assert_unchanged):
 def test_embeds_global_patterns_end_to_end():
     # Phase 4 exit gate 2: one foundation model plugged in end to end (decisions 10, 17).
     tdata = bt.pp.tax_glom(bt.datasets.global_patterns(), "genus")
-    with pytest.warns(UserWarning, match="mgm leaves out 96 of 996 features: 0 without a genus and 96 whose") as record:
+    with pytest.warns(UserWarning, match=re.escape(GLOBAL_PATTERNS_LEFT_OUT)) as record:
         bt.ml.embed(tdata, "mgm", inplace=True)
     assert not [warning for warning in record if "<bos> <eos>" in str(warning.message)]  # no sample is empty
     embedding = tdata.obsm["X_mgm"]
-    assert embedding.shape == (26, 256) and embedding.dtype == np.float32 and np.isfinite(embedding).all()
+    assert embedding.shape == GLOBAL_PATTERNS_SHAPE and embedding.dtype == np.float32 and np.isfinite(embedding).all()
     with pytest.warns(UserWarning, match="mgm leaves out 96"):
         np.testing.assert_allclose(bt.ml.embed(tdata, "mgm"), embedding, rtol=0, atol=ATOL)
     # Every sample's nearest neighbour in the embedding comes from the same environment.
@@ -146,6 +157,25 @@ def test_embeds_global_patterns_end_to_end():
     assert (types[neighbour] == types).all()


+@pytest.mark.mgm
+def test_the_embedding_tutorial_gives_the_outputs_it_quotes(run_page):
+    with pytest.warns(UserWarning, match=re.escape(GLOBAL_PATTERNS_LEFT_OUT)):
+        namespace = run_page(EMBEDDINGS)
+    assert (namespace["embedding"].shape, namespace["embedding"].dtype) == (GLOBAL_PATTERNS_SHAPE, np.float32)
+    assert namespace["same_type"] == GLOBAL_PATTERNS_SHAPE[0]
+
+
+def test_the_embedding_tutorial_quotes_the_end_to_end_outputs():
+    page = EMBEDDINGS.read_text(encoding="utf-8")
+    assert f"```text\nUserWarning: {GLOBAL_PATTERNS_LEFT_OUT}\n```" in page
+    assert f"```text\n({GLOBAL_PATTERNS_SHAPE}, dtype('float32'))\n```" in page
+    assert f"```text\n{GLOBAL_PATTERNS_SHAPE[0]}\n```" in page
+    assert (
+        f"For all {GLOBAL_PATTERNS_SHAPE[0]} samples, the nearest neighbour is a sample of the same type"
+        in " ".join(page.split())
+    )
+
+
 def _cached_file(name):
     """The file ``name`` that MGM's wheel was extracted to, in biotapy's data cache."""
     from biotapy.ml import _mgm
@@ -155,6 +185,28 @@ def _cached_file(name):
     return next(Path(path) for path in paths if Path(path).name == name)


+@pytest.mark.mgm
+def test_the_embedding_tutorial_s_vocabulary_and_example_genera_hold():
+    vocabulary = pd.read_csv(_cached_file("phylogeny.csv"), index_col=0).index
+    genera = set(bt.pp.tax_glom(bt.datasets.global_patterns(), "genus").var["genus"].dropna())
+    assert len(vocabulary) == MGM_GENERA
+    # The clone name and the one-word Candidatus genus the page names: biotapy reads "BD2-13" as MGM's "BD2".
+    assert {"BD2-13", "CandidatusPelagibacter"} <= genera
+    assert "g__BD2" not in vocabulary and "g__CandidatusPelagibacter" not in vocabulary
+    assert "g__Candidatus_Pelagibacter" in vocabulary
+
+
+def test_the_embedding_tutorial_quotes_its_vocabulary_and_labels_its_timings():
+    page = " ".join(EMBEDDINGS.read_text(encoding="utf-8").split())
+    samples, dimensions = GLOBAL_PATTERNS_SHAPE
+    assert f"MGM's vocabulary holds {MGM_GENERA:,} genera" in page
+    assert f"each sample comes back as a vector of {dimensions} numbers" in page
+    assert f"{samples} samples" in page
+    left_out, features = re.match(r"mgm leaves out (\d+) of (\d+) features", GLOBAL_PATTERNS_LEFT_OUT).groups()
+    assert f"{left_out} of GlobalPatterns' {features} genus-level features name a genus outside it" in page
+    assert "one unpinned laptop run, not checked by CI" in page
+
+
 @pytest.mark.mgm
 @pytest.mark.parametrize("name", ["phylogeny.csv", "config.json"])
 def test_a_damaged_extracted_file_is_extracted_again(name):
````
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/ml/test_mgm.py -q -k
  tutorial` -> `2 failed, 18 deselected` (`FileNotFoundError`); `uv run --group test --extra mgm
  pytest tests/ml/test_mgm.py -q -m mgm -k tutorial` -> `1 failed, 1 passed, 18 deselected` (the
  failure is `Failed: DID NOT WARN`, from the `FileNotFoundError` inside `pytest.warns`; the
  vocabulary test needs no page).
- [x] **Step 3: The page.** Create `docs/tutorials/embeddings.md` (no notebook header, so the
  docs build does not run it):
````markdown
# Embedding samples with MGM

This tutorial runs one pretrained microbiome model end to end: GlobalPatterns' genus profiles go
through MGM, each sample comes back as a vector of 256 numbers, and a nearest-neighbour check
shows that the vectors keep samples from the same environment together.

:::{note}
This page is not run when the documentation is built. MGM needs PyTorch and transformers (the
extra `mgm`), which the documentation build does not install. The code below is run on every
pull request instead, in the `ml-extras` job of biotapy's continuous integration, which checks
that it gives the outputs quoted here (`tests/ml/test_mgm.py`).
:::

```bash
pip install 'biotapy[mgm]'
```

## Embedding the samples

MGM reads genera, so the operational taxonomic units are merged to genus first. The first call
downloads MGM's pretrained model once (33 MB) into the cache `bt.datasets` uses.

```python
import biotapy as bt

genera = bt.pp.tax_glom(bt.datasets.global_patterns(), "genus")
bt.ml.embed(genera, "mgm", inplace=True)
embedding = genera.obsm["X_mgm"]
```

```text
UserWarning: mgm leaves out 96 of 996 features: 0 without a genus and 96 whose genus is not one of MGM's (4-29, 4041AA30, A17, Aquamonas, Arctic95A-2, ...)
```

MGM's vocabulary holds 9,665 genera; 96 of GlobalPatterns' 996 genus-level features name a genus
outside it and are left out. Most are environmental lineages known only by a clone name (`4-29`,
`BD2-13`) or *Candidatus* genera, which GlobalPatterns' Greengenes taxonomy writes as one word
(`CandidatusPelagibacter`) where MGM's vocabulary has `Candidatus_Pelagibacter`; biotapy reads
each name as MGM's own code does. Every sample keeps genera MGM knows, so none is embedded from
an empty profile.

```python
embedding.shape, embedding.dtype
```

```text
((26, 256), dtype('float32'))
```

## What the embedding keeps

MGM never saw GlobalPatterns' sample types. For each sample, the nearest other sample in the
embedding, by cosine distance, should still come from the same environment:

```python
from sklearn.neighbors import NearestNeighbors

nearest = NearestNeighbors(n_neighbors=2, metric="cosine").fit(embedding).kneighbors(embedding)[1][:, 1]
sample_type = genera.obs["SampleType"].to_numpy()
same_type = int((sample_type[nearest] == sample_type).sum())
same_type
```

```text
26
```

For all 26 samples, the nearest neighbour is a sample of the same type: feces next to feces, soil
next to soil, the three mock communities together. The embedding is a starting point for any
model that takes a fixed-length vector per sample; `genera.obsm["X_mgm"]` travels with the
table.

## Time and memory

In one unpinned laptop run, not checked by CI, on a CPU with 8 threads (torch 2.14.1+cpu,
transformers 5.19.0), the first call took 7 s, most of it importing torch and loading the model,
and a second call about 1 s for the 26 samples (the timing is rounded, so the rate is in the tens
of samples a second, as the guide says). MGM runs one sample at a time and needs no memory beyond the
model's.

## More

The {doc}`machine learning guide </guide/machine_learning>` says how MGM reads a table's genera,
which features it leaves out, how the embedding is pooled, how to cite MGM and its licence, and
how a package adds its own model as a plugin.
````
````diff
diff --git a/docs/guide/machine_learning.md b/docs/guide/machine_learning.md
--- a/docs/guide/machine_learning.md
+++ b/docs/guide/machine_learning.md
@@ -127,6 +127,9 @@ bt.ml.embed(genera, "mgm", inplace=True)
 genera.obsm["X_mgm"].shape  # (26, 256)
 ```

+The {doc}`embedding tutorial </tutorials/embeddings>` runs this on
+GlobalPatterns and checks what the embedding keeps.
+
 ### MGM

 `"mgm"` is MGM, the Microbial General Model of Zhang et al. (2026): a GPT-2
@@ -160,7 +163,8 @@ The embedding is the mean of the last hidden layer over the sample's tokens,
 the "element-wise mean pooling" MGM's paper uses for the pretrained model:
 256 float32 values per sample, within 2e-6 of MGM 0.5.8's own code. The model
 runs on the CPU, one sample at a time; on GlobalPatterns' genus profiles that
-is about 60 samples a second on 8 threads, with no memory beyond the model's.
+is a few tens of samples a second on 8 threads (one unpinned laptop run, not checked by
+CI; the [MGM tutorial](../tutorials/embeddings.md) gives its timings), with no memory beyond the model's.
 On a CPU running more than four threads, the first call in a session can
 differ from later ones by up to about 2e-4: torch 2.13 and 2.14 sometimes
 compute their first `tanh` less precisely.
diff --git a/docs/tutorials/index.md b/docs/tutorials/index.md
--- a/docs/tutorials/index.md
+++ b/docs/tutorials/index.md
@@ -9,4 +9,5 @@ phyloseq_analysis
 function
 differential_abundance
 leak_free_cv
+embeddings
 ```
````
- [x] **Step 4: Run, expect pass** - `uv run --group test --extra mgm pytest tests/ml/test_mgm.py
  -q -m "mgm or not mgm" -W error::UserWarning` -> `20 passed` (the end-to-end test included);
  then `uv sync --all-groups`.
- [x] **Step 5: Build.** The docs build succeeds; `tutorials/embeddings.html` exists and the log
  has no "Executed notebook" line for it.
- [x] **Step 6: Bookkeeping.** Tick 4.6b; log line:
  `- **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.6b done: \`docs/tutorials/embeddings.md\` shows MGM on GlobalPatterns' genera end to end, not run by the docs build (no torch); an \`mgm\` test runs its code and a default-run test checks it quotes the outputs of exit gate 2's test (\`GLOBAL_PATTERNS_LEFT_OUT\`, \`GLOBAL_PATTERNS_SHAPE\`).`
- [x] **Step 7: Gate and commit.**
  ```bash
  git add docs/tutorials/embeddings.md docs/tutorials/index.md docs/guide/machine_learning.md \
    tests/ml/test_mgm.py .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  git commit -m "docs: add the MGM embedding tutorial on GlobalPatterns

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1512 passed, 2 skipped, 95 deselected`; `37 passed, 1572
  deselected`; `build succeeded`; `25 passed, 1584 deselected`; `15 passed, 1594 deselected`.

### Task 4.6c: `bt.datasets.biocrust()` and the multi-omics tutorial (decision 37)

Drop this task if the user keeps the mmvec tutorial for 0.5; 4.D1's no-row test then leaves out
`bt.datasets.biocrust`, and the CHANGELOG and README lose their `biocrust` lines.

**Files:** create `src/biotapy/datasets/_biocrust.py`, `tests/datasets/test_biocrust.py`,
`docs/tutorials/multiomics.md`; modify `src/biotapy/datasets/__init__.py`,
`src/biotapy/datasets/_remote.py`, `tests/tl/test_mmvec.py`, `docs/api.md`,
`docs/guide/datasets.md`, `docs/guide/multiomics.md`, `docs/tutorials/index.md`,
`.knowledge/modules/datasets.md`, `.knowledge/modules/index.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `bt.tl.mmvec` (no option exposed: the tutorial's check holds at scikit-bio's
defaults, R2.3); `io.read_biom`, `io.to_mudata`; the CI cache keys (they hash `_remote.py`, so
the `network`, `r-bridge` and `docs` caches are rebuilt once, as in 4.C0).
**Interfaces:**
- Consumes: `_remote._fetch`, `bt.io.read_biom`, `bt.io.to_mudata`, `bt.tl.mmvec(mdata, *,
  microbes="taxa", metabolites="metabolites", seed=None)`, `run_page`.
- Produces: `bt.datasets.biocrust() -> MuData` (modalities `"taxa"`, `"metabolites"`, 19
  samples); `MULTIOMICS` and `MICROCOLEUS_METABOLITES` in `tests/tl/test_mmvec.py`; the page's
  names `mdata`, `ranks`, `cyanobacterium`, `microcoleus`, `above_zero`.

- [x] **Step 1: Failing loader tests.** Create `tests/datasets/test_biocrust.py`:
````python
import warnings
from pathlib import Path

import mudata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

import biotapy as bt
from biotapy.datasets import _biocrust, _remote

# Synthetic tables in the layout of mmvec's soil example: the metabolite table lacks one microbe sample and lists the
# others in another order. No data of the example is copied.
_MICROBES = AnnData(
    X=sp.csr_matrix(np.array([[10, 0, 3], [4, 6, 0], [0, 2, 8]])),
    obs=pd.DataFrame(index=["3min_early", "3min_late", "9hr_early"]),
    var=pd.DataFrame(index=["rplo 1 (Cyanobacteria)", "rplo 2 (Firmicutes)", "rplo 3 (Proteobacteria)"]),
)
_METABOLITES = AnnData(
    X=sp.csr_matrix(np.array([[1500.5, 2.0], [800.25, 31.75]])),
    obs=pd.DataFrame(index=["9hr_early", "3min_early"]),
    var=pd.DataFrame(index=["adenine", "uracil"]),
)

# The figures the docstring, the guide and the tutorial quote; the network test checks them on the downloaded files and
# the default-run test checks the prose quotes them.
DOWNLOAD_KB = 135
MICROBE_SAMPLES = 20
MOST_ABUNDANT = "rplo 1 (Cyanobacteria)"
ROOT = Path(__file__).parents[2]


@pytest.fixture
def fetched(tmp_path, monkeypatch):
    # As in test_hmp2.py: the only offline route to the loader is its private _fetch (R11.4).
    bt.io.write_biom(_MICROBES, tmp_path / "biocrust_microbes.biom")
    bt.io.write_biom(_METABOLITES, tmp_path / "biocrust_metabolites.biom")
    names = []
    monkeypatch.setattr(_biocrust, "_fetch", lambda name: names.append(name) or str(tmp_path / name))
    return names


def test_keeps_the_samples_both_tables_have_in_the_microbe_order_without_a_warning(fetched):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mdata = bt.datasets.biocrust()
    assert isinstance(mdata, mudata.MuData) and list(mdata.mod) == ["taxa", "metabolites"]
    for mod in mdata.mod.values():
        assert mod.obs_names.tolist() == ["3min_early", "9hr_early"]


def test_values_and_units_come_from_the_reader(fetched):
    mdata = bt.datasets.biocrust()
    assert mdata["taxa"].uns["biotapy"]["x_kind"] == "counts"
    assert mdata["metabolites"].uns["biotapy"]["x_kind"] == "abundance"
    np.testing.assert_array_equal(mdata["taxa"].X.toarray(), [[10, 0, 3], [0, 2, 8]])
    np.testing.assert_array_equal(mdata["metabolites"].X.toarray(), [[800.25, 31.75], [1500.5, 2.0]])


def test_fetches_both_pinned_files_before_reading(fetched):
    bt.datasets.biocrust()
    assert fetched == ["biocrust_microbes.biom", "biocrust_metabolites.biom"]


@pytest.mark.parametrize("page", ["guide/datasets.md", "tutorials/multiomics.md"])
def test_the_pages_link_mmvec_s_example_at_the_pinned_commit(page):
    commit = _remote._MMVEC.split("/mmvec/")[1].split("/")[0]
    text = (ROOT / "docs" / page).read_text(encoding="utf-8")
    assert f"https://github.com/biocore/mmvec/tree/{commit}/examples/soils" in text


def test_the_prose_quotes_the_figures_the_network_test_checks():
    def prose(path):
        return " ".join((ROOT / path).read_text(encoding="utf-8").split())

    docstring, guide, tutorial = (
        prose(path)
        for path in ["src/biotapy/datasets/_biocrust.py", "docs/guide/datasets.md", "docs/tutorials/multiomics.md"]
    )
    for text in (docstring, guide, tutorial):
        assert f"({DOWNLOAD_KB} KB)" in text
    assert f"counted in {MICROBE_SAMPLES} samples" in guide
    assert f"{MICROBE_SAMPLES} samples" in docstring
    assert f"The microbe table has {MICROBE_SAMPLES} samples" in tutorial
    assert f"The most abundant microbe is `{MOST_ABUNDANT}`" in tutorial


@pytest.mark.network
def test_biocrust_downloads_and_loads():
    mdata = bt.datasets.biocrust()
    assert {key: mod.shape for key, mod in mdata.mod.items()} == {"taxa": (19, 466), "metabolites": (19, 85)}
    assert "9hr_late" not in mdata.obs_names


@pytest.mark.network
def test_biocrust_prose_figures_hold():
    # The figures the guide, the tutorial and the docstring quote about the downloaded files.
    paths = [Path(_biocrust._fetch(name)) for name in _biocrust.FILES]
    assert round(sum(path.stat().st_size for path in paths) / 1000) == DOWNLOAD_KB
    microbes = bt.io.read_biom(paths[0])
    assert microbes.n_obs == MICROBE_SAMPLES
    assert microbes.var_names[np.asarray(microbes.X.sum(axis=0)).ravel().argmax()] == MOST_ABUNDANT
````
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/datasets/test_biocrust.py
  -q` -> `1 error` during collection: `ImportError: cannot import name '_biocrust' from
  'biotapy.datasets'`.
- [x] **Step 3: The loader.** Create `src/biotapy/datasets/_biocrust.py`:
````python
"""A desert biocrust wetting experiment: microbes and metabolites over the same samples, from mmvec's repository."""

from mudata import MuData

from biotapy.io import read_biom, to_mudata

from ._remote import _fetch

# Fetched together, before either is parsed, so a download error comes first.
FILES = ["biocrust_microbes.biom", "biocrust_metabolites.biom"]


def biocrust() -> MuData:
    """Microbes and metabolites of a desert biological soil crust after wetting, over the samples both tables have.

    Downloaded once (135 KB) from the example of mmvec's repository
    (``examples/soils``), pinned to one commit: 466 microbes counted in 20
    samples and 85 metabolites measured in 19 of them. The 20 samples are
    four biocrust successional stages (``early``, ``earlymid``, ``latemid``,
    ``late``) at five times after wetting.

    Returns
    -------
    MuData
        The 19 samples both tables have, in the microbe table's order, in
        two modalities: ``"taxa"``, the microbes' counts
        (``x_kind == "counts"``, no taxonomy columns and no tree; each
        feature name ends with its phylum), and ``"metabolites"``, the
        metabolites' intensities (``x_kind == "abundance"``), both as
        ``bt.io.read_biom`` reads them. The tables carry no sample metadata;
        a sample's name holds its time after wetting and its crust's
        successional stage.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/datasets`

    The microbe table's sample ``9hr_late`` has no metabolite profile and
    is left out without a warning; the files themselves are unchanged.

    The files are distributed in mmvec's repository under its BSD-3-Clause
    licence. biotapy ships none of them; cite the mmvec paper when you use
    them.

    References
    ----------
    Swenson TL et al. (2018) Linking soil biology and chemistry in biological soil crust
    using isolate exometabolomics. Nature Communications 9:19.

    Morton JT et al. (2019) Learning representations of microbe-metabolite interactions.
    Nature Methods 16:1306-1314.

    Examples
    --------
    >>> import biotapy as bt
    >>> mdata = bt.datasets.biocrust()  # doctest: +SKIP
    >>> mdata["taxa"].shape, mdata["metabolites"].shape  # doctest: +SKIP
    ((19, 466), (19, 85))
    """
    paths = [_fetch(name) for name in FILES]
    microbes, metabolites = (read_biom(path) for path in paths)
    shared = microbes.obs_names.intersection(metabolites.obs_names, sort=False)
    return to_mudata({"taxa": microbes[shared], "metabolites": metabolites[shared]})
````
````diff
diff --git a/docs/api.md b/docs/api.md
--- a/docs/api.md
+++ b/docs/api.md
@@ -32,6 +32,7 @@ Public functions are listed here as they ship, from Phase 1 onward.
 .. autosummary::
     :toctree: generated

+    datasets.biocrust
     datasets.enterotype
     datasets.enzyme
     datasets.esophagus
diff --git a/docs/guide/datasets.md b/docs/guide/datasets.md
--- a/docs/guide/datasets.md
+++ b/docs/guide/datasets.md
@@ -1,9 +1,10 @@
 # Example datasets

-`biotapy.datasets` ships seven datasets: six examples and the ENZYME hierarchy
+`biotapy.datasets` ships eight datasets: seven examples and the ENZYME hierarchy
 for `bt.fn.func_glom`. Four of the examples return the [data model](data_model.md)
 every biotapy function relies on; `toy_humann` returns a `MuData` function
-table, and `hmp2` a `MuData` cohort of function tables and taxa.
+table, `hmp2` a `MuData` cohort of function tables and taxa, and `biocrust` a
+`MuData` of microbes and metabolites.

 ## `toy`

@@ -61,7 +62,8 @@ Each of `global_patterns`, `enterotype` and `esophagus` is downloaded once from
 [pooch](https://www.fatiando.org/pooch/). A later call re-hashes the cached
 file and, as long as the hash still matches, reads it straight from disk
 with no network access at all. `hmp2`'s three files are pinned the same way,
-to the SHA-256 hashes of the IBDMDB's dated releases.
+to the SHA-256 hashes of the IBDMDB's dated releases, and `biocrust`'s two to
+one commit of mmvec's repository.

 Set `BIOTAPY_DATA_DIR` to change the cache directory; the default is a
 per-user cache directory (`pooch.os_cache("biotapy")`):
@@ -111,7 +113,32 @@ mdata = bt.datasets.hmp2()
 mdata.obs["diagnosis"].value_counts()  # CD 65, UC 38, nonIBD 27
 ```

-The [function tutorial](../tutorials/function.md) analyses it.
+The [function](../tutorials/function.md) and [leak-free cross-validation](../tutorials/leak_free_cv.md)
+tutorials analyse it.
+
+## `biocrust`
+
+`bt.datasets.biocrust()` is the example mmvec's own repository uses: a desert
+biological soil crust from four successional stages (`early`, `earlymid`,
+`latemid`, `late`; Swenson et al. 2018), sampled at five times after wetting,
+with 466 microbes counted in 20 samples and 85 metabolites measured in 19 of
+them. It downloads
+the two BIOM files once (135 KB) and returns a `MuData` over the 19 samples
+both have:
+
+| Modality | Holds | `x_kind` |
+|---|---|---|
+| `"taxa"` | microbe counts, no taxonomy columns and no tree; each name ends with its phylum | `counts` |
+| `"metabolites"` | metabolite intensities | `abundance` |
+
+```python
+import biotapy as bt
+
+mdata = bt.datasets.biocrust()
+ranks = bt.tl.mmvec(mdata, seed=0)  # 466 microbes x 85 metabolites
+```
+
+The [multi-omics tutorial](../tutorials/multiomics.md) analyses it.

 ## Licensing

@@ -124,7 +151,13 @@ you publish results that use it. `hmp2` downloads the HMP2 tables from the
 [IBDMDB](https://ibdmdb.org/); the IBDMDB states no licence for them, so biotapy
 ships none of them, and you should cite the study when you use them: Lloyd-Price J
 et al. (2019) Multi-omics of the gut microbial ecosystem in inflammatory bowel
-diseases. *Nature* 569:655-662. biotapy itself is
+diseases. *Nature* 569:655-662. `biocrust` downloads two files from
+[mmvec's repository](https://github.com/biocore/mmvec/tree/88ca33b408a85b6bf90fae06982936247b860272/examples/soils),
+distributed under its BSD-3-Clause licence; cite Morton JT et al. (2019) Learning
+representations of microbe-metabolite interactions. *Nature Methods*
+16:1306-1314. The data come from Swenson TL et al. (2018) Linking soil biology and
+chemistry in biological soil crust using isolate exometabolomics. *Nature
+Communications* 9:19. biotapy itself is
 [BSD-3-Clause](https://github.com/pedrocr83/biotapy/blob/master/LICENSE).

 [phyloseq-data]: https://github.com/joey711/phyloseq/tree/master/data
diff --git a/src/biotapy/datasets/__init__.py b/src/biotapy/datasets/__init__.py
--- a/src/biotapy/datasets/__init__.py
+++ b/src/biotapy/datasets/__init__.py
@@ -1,6 +1,7 @@
+from ._biocrust import biocrust
 from ._enzyme import enzyme
 from ._hmp2 import hmp2
 from ._remote import enterotype, esophagus, global_patterns
 from ._toy import toy, toy_humann

-__all__ = ["enterotype", "enzyme", "esophagus", "global_patterns", "hmp2", "toy", "toy_humann"]
+__all__ = ["biocrust", "enterotype", "enzyme", "esophagus", "global_patterns", "hmp2", "toy", "toy_humann"]
diff --git a/src/biotapy/datasets/_remote.py b/src/biotapy/datasets/_remote.py
--- a/src/biotapy/datasets/_remote.py
+++ b/src/biotapy/datasets/_remote.py
@@ -1,4 +1,4 @@
-"""Datasets downloaded once and cached with pooch: phyloseq's examples, the ENZYME files and the HMP2 tables."""
+"""Datasets downloaded once and cached with pooch: phyloseq's examples, the ENZYME files, the HMP2 and biocrust tables."""

 from functools import cache

@@ -21,14 +21,20 @@ _REGISTRY = {
     "pathabundances_3.tsv.gz": "sha256:dd983871b0e155255844b91ec10d50fb09230d2f4e915464ab680fa3a9c9ddb3",
     "taxonomic_profiles_3.tsv.gz": "sha256:d790ff15e46d61ca0cadc55d9f918de4e3415d7f97c992ac37610aaee02117ed",
     "hmp2_metadata_2018-08-20.csv": "sha256:656b7bd97660ddb875548805e30bede31f2d1208293f7170d2d5755e33862ec9",
+    # mmvec's soil biocrust example, pinned to one commit of biocore/mmvec.
+    "biocrust_microbes.biom": "sha256:7f634b8dcdbe3ad0c25a97223e2e3b130cfc7ae69100c21415bffbc3d86ef4a7",
+    "biocrust_metabolites.biom": "sha256:cfea3e7233102680f628e71c75229a2a55640fa6c704583ce89b0b9b77f32c72",
 }
 _IBDMDB = "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/"
+_MMVEC = "https://raw.githubusercontent.com/biocore/mmvec/88ca33b408a85b6bf90fae06982936247b860272/examples/soils/"
 _URLS = {
     "enzyme.dat": "https://ftp.expasy.org/databases/enzyme/enzyme.dat",
     "enzclass.txt": "https://ftp.expasy.org/databases/enzyme/enzclass.txt",
     "pathabundances_3.tsv.gz": f"{_IBDMDB}products/HMP2/MGX/2018-05-04/pathabundances_3.tsv.gz",
     "taxonomic_profiles_3.tsv.gz": f"{_IBDMDB}products/HMP2/MGX/2018-05-04/taxonomic_profiles_3.tsv.gz",
     "hmp2_metadata_2018-08-20.csv": f"{_IBDMDB}metadata/hmp2_metadata_2018-08-20.csv",
+    "biocrust_microbes.biom": f"{_MMVEC}microbes.biom",
+    "biocrust_metabolites.biom": f"{_MMVEC}metabolites.biom",
 }


````
  The hashes are the SHA-256 of the two files at commit `88ca33b` (checked by downloading
  both, `sha256sum`); R6.6: nothing is committed but synthetic fixtures.
- [x] **Step 4: Run, expect pass** - `uv run --group test pytest tests/datasets/test_biocrust.py
  -q -m "network or not network" -W error::UserWarning` -> `8 passed` (the network tests
  download 135 KB).
- [x] **Step 5: Failing tutorial tests.** In `tests/tl/test_mmvec.py`:
````diff
diff --git a/tests/tl/test_mmvec.py b/tests/tl/test_mmvec.py
--- a/tests/tl/test_mmvec.py
+++ b/tests/tl/test_mmvec.py
@@ -1,3 +1,5 @@
+from pathlib import Path
+
 import numpy as np
 import pandas as pd
 import pytest
@@ -136,3 +138,38 @@ def test_every_microbe_row_is_centred(counts):
     )
     ranks = bt.tl.mmvec(bt.io.to_mudata({"taxa": taxa, "metabolites": metabolites}), seed=0)
     np.testing.assert_allclose(ranks.sum(axis=1), 0.0, atol=1e-9)
+
+
+MULTIOMICS = Path(__file__).parents[2] / "docs" / "tutorials" / "multiomics.md"
+# The 13 metabolites mmvec's own soil example checks for rplo 1, Microcoleus vaginatus (biocore/mmvec
+# examples/soils/check_soils.ipynb at 88ca33b): each ranks above zero in that microbe's row. The multi-omics tutorial
+# runs the check on bt.datasets.biocrust() and quotes the count.
+MICROCOLEUS_METABOLITES = {
+    "(3-methyladenine)",
+    "7-methyladenine",
+    "4-guanidinobutanoate",
+    "uracil",
+    "xanthine",
+    "hypoxanthine",
+    "(N6-acetyl-lysine)",
+    "cytosine",
+    "N-acetylornithine",
+    "succinate",
+    "adenosine",
+    "guanine",
+    "adenine",
+}
+
+
+@pytest.mark.network
+def test_the_multiomics_tutorial_ranks_mmvec_s_microcoleus_metabolites_above_zero(run_page):
+    namespace = run_page(MULTIOMICS)
+    assert namespace["microcoleus"] == MICROCOLEUS_METABOLITES
+    assert namespace["above_zero"] == len(MICROCOLEUS_METABOLITES)
+
+
+def test_the_multiomics_tutorial_quotes_its_check():
+    page = " ".join(MULTIOMICS.read_text(encoding="utf-8").split())
+    count = len(MICROCOLEUS_METABOLITES)
+    assert f"checks its fit of these data against {count} metabolites it lists for this microbe" in page
+    assert f"All {count} of the metabolites mmvec's example lists rank above zero for the cyanobacterium" in page
````
  Run `uv run --group test pytest tests/tl/test_mmvec.py -q -m "network or not network" -k
  tutorial` -> `2 failed, 14 deselected` (`FileNotFoundError`).
- [x] **Step 6: The tutorial.** Create `docs/tutorials/multiomics.md`:
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

# Microbes and metabolites in a desert biocrust

This tutorial puts two data types measured on the same samples into one `MuData`, then asks
which metabolites go with which microbe, with mmvec (Morton et al. 2019). The data are the
example from mmvec's own repository: a desert biological soil crust from four successional stages
(`early`, `earlymid`, `latemid`, `late`), sampled at five times after wetting, with its microbes
counted and its metabolites measured (Swenson et al. 2018).

:::{note}
The two tables are downloaded from [mmvec's repository](https://github.com/biocore/mmvec/tree/88ca33b408a85b6bf90fae06982936247b860272/examples/soils)
on first use (135 KB) and cached; they are distributed there under its BSD-3-Clause licence.
Cite the mmvec paper when you use them: Morton JT et al. (2019) Learning representations of
microbe-metabolite interactions. *Nature Methods* 16:1306-1314. The data come from Swenson TL et
al. (2018) Linking soil biology and chemistry in biological soil crust using isolate
exometabolomics. *Nature Communications* 9:19.
:::

```{code-cell} ipython3
import biotapy as bt
```

## Two tables, one MuData

The microbe table has 20 samples and the metabolite table 19; `bt.datasets.biocrust()` keeps the
19 both have, as `bt.io.to_mudata` does for your own tables. Each data type is a modality, named
as biotapy's {doc}`multi-omics guide </guide/multiomics>` names them.

```{code-cell} ipython3
mdata = bt.datasets.biocrust()
mdata
```

## Which metabolites go with which microbe

`bt.tl.mmvec` learns, from the samples, how likely each metabolite is given each microbe. Its
result has one row per microbe and one column per metabolite, holding log probabilities centred
so each row sums to 0: within a row, a large value is a metabolite that tends to be abundant
where that microbe is. The fit starts from random values, so `seed` fixes it.

```{code-cell} ipython3
ranks = bt.tl.mmvec(mdata, seed=0)
ranks.shape
```

The most abundant microbe is `rplo 1 (Cyanobacteria)`, which mmvec's example treats as the
cyanobacterium *Microcoleus* sp. (closest isolate *M. vaginatus* PCC 9802). Its ten highest-ranked metabolites:

```{code-cell} ipython3
cyanobacterium = ranks.loc["rplo 1 (Cyanobacteria)"].sort_values(ascending=False)
cyanobacterium.head(10)
```

## A check against mmvec's own example

mmvec's repository checks its fit of these data against 13 metabolites it lists for this
microbe: every one should rank above zero in its row. biotapy fits mmvec through scikit-bio
rather than mmvec's own TensorFlow code, so the same check says whether the answer holds:

```{code-cell} ipython3
microcoleus = {
    "(3-methyladenine)", "7-methyladenine", "4-guanidinobutanoate", "uracil", "xanthine", "hypoxanthine",
    "(N6-acetyl-lysine)", "cytosine", "N-acetylornithine", "succinate", "adenosine", "guanine", "adenine",
}
above_zero = int((cyanobacterium[list(microcoleus)] > 0).sum())
above_zero
```

All 13 of the metabolites mmvec's example lists rank above zero for the cyanobacterium, as in
mmvec's own fit, although biotapy's fit differs in its optimiser and in its number of
dimensions (scikit-bio's default of 3, where the example used 1).

## Reading the result

- A rank compares metabolites within one microbe's row. It is not a correlation, and a low rank
  means no association rather than a negative one.
- The page runs scikit-bio's default of 3 latent dimensions, while mmvec's README advises at least
  10 samples per dimension, which allows 1 dimension for 19 samples. Read the ranks with that in mind.
- 19 samples is few: mmvec's README says studies this small need careful tuning of the number
  of dimensions and the priors, which scikit-bio's `mmvec` exposes and `bt.tl.mmvec` leaves at
  their defaults.
- The {doc}`multi-omics guide </guide/multiomics>` says how to build the `MuData` from your own
  tables and what `bt.tl.mmvec` refuses.
````
````diff
diff --git a/docs/guide/multiomics.md b/docs/guide/multiomics.md
--- a/docs/guide/multiomics.md
+++ b/docs/guide/multiomics.md
@@ -3,7 +3,8 @@
 When several data types are measured on the same samples - taxa, a function
 table, metabolites, host data - biotapy keeps them in one
 [MuData](https://mudata.scverse.org/): one AnnData per data type, called a
-modality, over shared samples.
+modality, over shared samples. The {doc}`multi-omics tutorial
+</tutorials/multiomics>` runs mmvec on a real pair of tables.

 ## Modality names

diff --git a/docs/tutorials/index.md b/docs/tutorials/index.md
--- a/docs/tutorials/index.md
+++ b/docs/tutorials/index.md
@@ -8,6 +8,7 @@ quick_tour
 phyloseq_analysis
 function
 differential_abundance
+multiomics
 leak_free_cv
 embeddings
 ```
````
  Run `uv run --group test pytest tests/tl/test_mmvec.py tests/datasets/test_biocrust.py -q -m
  "network or not network" -W error::UserWarning` -> `24 passed`. Build the docs:
  `tutorials/multiomics.md: Executed notebook` (3.0-3.6 s); its outputs read `MuData object with
  n_obs × n_vars = 19 × 551` with `taxa: 19 × 466` and `metabolites: 19 × 85`; `(466, 85)`;
  adenine 5.05, adenosine 5.02, (N6-acetyl-lysine) 4.48, 4-guanidinobutanoate 3.98, isoleucine
  3.78 at the top of the ten; `13`.
- [x] **Step 7: Knowledge.** Tick the new checklist line `4.6c`; the module concept and its
  index description:
````diff
diff --git a/.knowledge/log.md b/.knowledge/log.md
--- a/.knowledge/log.md
+++ b/.knowledge/log.md
@@ -1,6 +1,7 @@
 # Knowledge bundle log

 ## 2026-10-10 (Phase 4, slice 4D)
+- **Update**: [datasets](modules/datasets.md) gains `biocrust()`, mmvec's soil example of microbes and metabolites as a MuData over the 19 shared samples, pinned to one mmvec commit and two SHA-256 hashes (description copied into the [modules index](modules/index.md)); [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.6c done: `docs/tutorials/multiomics.md` runs `bt.tl.mmvec` on it on every docs build, and a network test checks the page still ranks mmvec's 13 Microcoleus metabolites above zero.
 - **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.6b done: `docs/tutorials/embeddings.md` shows MGM on GlobalPatterns' genera end to end, not run by the docs build (no torch); an `mgm` test runs its code and a default-run test checks it quotes the outputs of exit gate 2's test (`GLOBAL_PATTERNS_LEFT_OUT`, `GLOBAL_PATTERNS_SHAPE`).
 - **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.6 done: `docs/tutorials/leak_free_cv.md` cross-validates a prevalence filter, CLR and logistic regression on HMP2's species on every docs build, inside and outside the pipeline, with a label-reading `SelectKBest` on real and shuffled labels; a network test runs the page's cells against the numbers its prose quotes (`LEAK_FREE_CV_AUC`).
 - **Update**: [core](modules/core.md): `sum_pairs` sums in `numpy.sum`'s dtype for `X`, as `sum_by` does, so `fn.func_glom` cannot wrap an int8 table; [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks 4.F5.
diff --git a/.knowledge/modules/datasets.md b/.knowledge/modules/datasets.md
--- a/.knowledge/modules/datasets.md
+++ b/.knowledge/modules/datasets.md
@@ -1,12 +1,12 @@
 ---
 type: Module
 title: datasets
-description: In-memory and pooch-cached example data for docs, doctests and tests - TreeData objects, a HUMAnN-style function MuData, the HMP2 cohort as a three-modality MuData, and the ENZYME hierarchy as an edge table.
+description: In-memory and pooch-cached example data for docs, doctests and tests - TreeData objects, a HUMAnN-style function MuData, the HMP2 cohort as a three-modality MuData, mmvec's soil biocrust example as a microbes-and-metabolites MuData, and the ENZYME hierarchy as an edge table.
 resource: /src/biotapy/datasets/
 paths: ["src/biotapy/datasets/**"]
 tags: [datasets]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T16:44:46Z }
-commit: 31aa11d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T17:01:49Z }
+commit: f5fb966
 status: stable
 ---

@@ -17,8 +17,9 @@ entirely in memory; `global_patterns()`/`enterotype()`/`esophagus()`,
 downloaded once from phyloseq's repository and cached with pooch (Task 1.11,
 esophagus Task 1.15b); `enzyme()`, the ENZYME EC hierarchy downloaded from
 ExPASy the same way; and `hmp2()`, the HMP2 (IBDMDB) cohort's pathway and
-taxon tables (Task 2.10). The return type varies: a TreeData, a MuData
-(`toy_humann`, `hmp2`) or a `pandas.DataFrame` (`enzyme`).
+taxon tables (Task 2.10); and `biocrust()`, mmvec's soil example of microbes
+and metabolites (Phase 4 Task 4.6c). The return type varies: a TreeData, a MuData
+(`toy_humann`, `hmp2`, `biocrust`) or a `pandas.DataFrame` (`enzyme`).

 # Entry points

@@ -44,6 +45,12 @@ taxon tables (Task 2.10). The return type varies: a TreeData, a MuData
   `visit_num` only orders and is not in `obs`, and a missing one sorts last.
   Seven metadata columns (`_hmp2.py:COLUMNS`) sit in the global `obs` and are
   pushed into every modality; `diagnosis` is categorical `nonIBD`, `UC`, `CD`.
+- `_biocrust.py:biocrust` - mmvec's soil example (`biocore/mmvec`,
+  `examples/soils`): 466 microbe counts in 20 samples and 85 metabolite
+  intensities in 19, read by `io.read_biom` and combined by `io.to_mudata` as
+  `taxa` and `metabolites` over the 19 shared samples, in the microbe table's
+  order. The loader intersects the samples itself, so `to_mudata` does not
+  warn about the one microbe sample (`9hr_late`) without metabolites.
 - `_remote.py:global_patterns` - GlobalPatterns: 26 samples x 19,216 OTUs,
   with taxonomy and a tree, read through `bt.io.read_phyloseq`.
 - `_remote.py:enterotype` - enterotype: 280 samples x 553 genera, as relative
@@ -78,6 +85,11 @@ taxon tables (Task 2.10). The return type varies: a TreeData, a MuData
   sample metadata; only its 1,638 stool metagenomes are used, and they share
   their ids with both tables (checked 2026-10-05), so the selected samples
   index every table. `_hmp2.py:hmp2`
+- `biocrust()`'s two files are pinned to one commit of `biocore/mmvec`
+  (`88ca33b`, `_remote.py:_MMVEC`) and to SHA-256 hashes; both are fetched
+  before either is parsed (`_biocrust.py:FILES`). Their registry names carry a
+  `biocrust_` prefix because the repository calls them `microbes.biom` and
+  `metabolites.biom`.
 - `enzyme()` is the one download that carries no pinned hash: ENZYME keeps only
   its current release online, so `_remote.py:_REGISTRY` lists `enzyme.dat` and
   `enzclass.txt` with `None`. The first download is cached for good, and
@@ -90,13 +102,14 @@ taxon tables (Task 2.10). The return type varies: a TreeData, a MuData
 - [core](/modules/core.md): `make_treedata`, `make_function_mudata`,
   `tree_from_edges`, `TreeData`.
 - [io](/modules/io.md): `read_phyloseq`, used by `_remote.py`'s three loaders;
-  `read_humann` and `read_metaphlan`, used by `_hmp2.py:hmp2`.
+  `read_humann` and `read_metaphlan`, used by `_hmp2.py:hmp2`; `read_biom` and
+  `to_mudata`, used by `_biocrust.py:biocrust`.
 - `mudata`: `hmp2` builds its MuData and pushes the global `obs` into the
   modalities (`MuData.push_obs`).
 - `pooch` (runtime, Task 1.11): fetches and caches `GlobalPatterns.RData`/
   `enterotype.RData`/`esophagus.RData` and, from the ExPASy FTP site,
   `enzyme.dat`/`enzclass.txt`, and from the IBDMDB's Globus endpoint the three
-  HMP2 files (`_remote.py:_URLS`); `BIOTAPY_DATA_DIR` overrides its
+  HMP2 files, and from GitHub mmvec's two BIOM files (`_remote.py:_URLS`); `BIOTAPY_DATA_DIR` overrides its
   cache directory (the pooch is built by `_core._download.py:make_pooch`,
   [core](/modules/core.md), which `ml/_mgm.py` shares).

@@ -106,8 +119,9 @@ taxon tables (Task 2.10). The return type varies: a TreeData, a MuData
 for `toy()` and `_remote.py`'s offline test; add
 `BIOTAPY_DATA_DIR=<dir> uv run --group test pytest -m network tests/datasets -q`
 to actually exercise the downloads (CI's dedicated `network` job runs
-`-m "network or golden"`). The docs job runs `docs/tutorials/function.md`,
-which calls `hmp2()`.
+`-m "network or golden"`). The docs job runs `docs/tutorials/function.md`
+and `docs/tutorials/leak_free_cv.md`, which call `hmp2()`, and
+`docs/tutorials/multiomics.md`, which calls `biocrust()`.

 # Gotchas

diff --git a/.knowledge/modules/index.md b/.knowledge/modules/index.md
--- a/.knowledge/modules/index.md
+++ b/.knowledge/modules/index.md
@@ -6,6 +6,6 @@
 * [tl](tl.md) - Diversity, ordination, PERMANOVA and mmvec over AnnData/TreeData/MuData - alpha, beta, UniFrac, PCoA, NMDS and PERMANOVA through scikit-bio and scikit-learn, and microbe-metabolite co-occurrence between two modalities, returning results or writing the data-model-slots keys.
 * [pl](pl.md) - Plots of what tl, pp, fn and da give - stacked bars, heatmap, a function's contributions per taxon, richness, ordination, scree and the da consensus dots - drawn with matplotlib on the given or a new Axes, computing nothing.
 * [fn](fn.md) - Function hierarchies, aggregation along them, HUMAnN-style renormalisation, per-taxon contributions and functional redundancy (Tian 2020) over function tables; owns no reader and no download.
-* [datasets](datasets.md) - In-memory and pooch-cached example data for docs, doctests and tests - TreeData objects, a HUMAnN-style function MuData, the HMP2 cohort as a three-modality MuData, and the ENZYME hierarchy as an edge table.
+* [datasets](datasets.md) - In-memory and pooch-cached example data for docs, doctests and tests - TreeData objects, a HUMAnN-style function MuData, the HMP2 cohort as a three-modality MuData, mmvec's soil biocrust example as a microbes-and-metabolites MuData, and the ENZYME hierarchy as an edge table.
 * [ml](ml.md) - scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold, taking arrays, sparse matrices and DataFrames; to_torch, a PyTorch dataset over an AnnData's rows behind the extra torch; and embed, one embedding per sample from a model a plugin registers.
 * [da](da.md) - Four differential abundance methods, native LinDA and ANCOM-BC2 and the R bridges ALDEx2 and MaAsLin 3 (rpy2, extra `r`), that return one result table schema, and a consensus table counting where the methods agree; da writes no slot and filters nothing.
````
- [x] **Step 8: Gate and commit.**
  ```bash
  git add src/biotapy/datasets/_biocrust.py src/biotapy/datasets/__init__.py src/biotapy/datasets/_remote.py \
    tests/datasets/test_biocrust.py tests/tl/test_mmvec.py docs/tutorials/multiomics.md docs/tutorials/index.md \
    docs/guide/multiomics.md docs/guide/datasets.md docs/api.md .knowledge/modules/datasets.md \
    .knowledge/modules/index.md .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  git commit -m "feat(datasets): add biocrust, mmvec's soil example, and the multi-omics tutorial

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1521 passed, 2 skipped, 98 deselected`; `40 passed, 1581
  deselected`; `build succeeded` with seven executed notebooks.

### Task 4.D1: Coming-from-R check

**Files:** modify `tests/test_coming_from_r.py`, `docs/_data/r_idioms.toml`,
`docs/coming_from_r.md`, `docs/tutorials/phyloseq_analysis.md`,
`.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
**Not touched:** `docs/extensions/coming_from_r.py` (it reads every `R equivalent:` line,
classes included); any docstring; no new idiom rows (R1.4).
**Interfaces:**
- Consumes: `coming_from_r.rows()`; `io.to_mudata`'s `R equivalent:` line.
- Produces: the label `"not in 0.4"`; pins for the two MultiAssayExperiment rows.

- [x] **Step 1: Failing tests.**
````diff
diff --git a/tests/test_coming_from_r.py b/tests/test_coming_from_r.py
--- a/tests/test_coming_from_r.py
+++ b/tests/test_coming_from_r.py
@@ -16,7 +16,7 @@ PHYLOSEQ_31 = """otu_table sample_data tax_table phy_tree refseq nsamples ntaxa
 sample_sums taxa_sums rank_names sample_variables get_taxa_unique prune_taxa prune_samples subset_taxa
 subset_samples filter_taxa transform_sample_counts rarefy_even_depth tax_glom estimate_richness distance
 UniFrac ordinate plot_bar plot_richness plot_ordination plot_heatmap import_biom""".split()
-NOT_IN_0_3 = ["tip_glom", "merge_samples", "psmelt", "plot_tree", "plot_net"]
+NOT_IN_0_4 = ["tip_glom", "merge_samples", "psmelt", "plot_tree", "plot_net"]


 def test_the_list_has_31_functions():
@@ -26,12 +26,12 @@ def test_the_list_has_31_functions():
 @pytest.mark.parametrize("name", PHYLOSEQ_31)
 def test_table_maps_each_of_the_31(name):
     cells = coming_from_r.rows()[f"phyloseq::{name}"]
-    assert cells and "not in 0.3" not in cells
+    assert cells and "not in 0.4" not in cells


-@pytest.mark.parametrize("name", NOT_IN_0_3)
+@pytest.mark.parametrize("name", NOT_IN_0_4)
 def test_uncovered_functions_are_marked(name):
-    assert coming_from_r.rows()[f"phyloseq::{name}"] == ["not in 0.3"]
+    assert coming_from_r.rows()[f"phyloseq::{name}"] == ["not in 0.4"]


 @pytest.mark.parametrize(
@@ -57,6 +57,22 @@ def test_transforms_and_da_methods_map_to_their_r_functions(r_name, functions):
     assert coming_from_r.rows()[r_name] == [f"{{func}}`bt.{name} <biotapy.{name}>`" for name in functions]


+@pytest.mark.parametrize(
+    "r_name", ["MultiAssayExperiment::MultiAssayExperiment", "MultiAssayExperiment::intersectColumns"]
+)
+def test_multi_assay_experiment_maps_to_to_mudata(r_name):
+    assert coming_from_r.rows()[r_name] == ["{func}`bt.io.to_mudata <biotapy.io.to_mudata>`"]
+
+
+def test_ml_mmvec_and_biocrust_add_no_row():
+    # Phase 4: the ml transformers, to_torch, embed, tl.mmvec and datasets.biocrust have no R equivalent; to_mudata is
+    # the only new row.
+    cells = [cell for cells in coming_from_r.rows().values() for cell in cells]
+    assert [
+        cell for cell in cells if any(name in cell for name in ("bt.ml.", "bt.tl.mmvec", "bt.datasets.biocrust"))
+    ] == []
+
+
 def test_plot_functions_link_to_pl():
     assert coming_from_r.rows()["phyloseq::plot_bar"] == ["{func}`bt.pl.bar <biotapy.pl.bar>`"]

````
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/test_coming_from_r.py
  -q` -> `5 failed, 48 passed`: each `test_uncovered_functions_are_marked` case (`['not in 0.3']
  == ['not in 0.4']`). The three new tests pass at once: they pin rows the 4A-4C docstrings
  already give.
- [x] **Step 3: Implement.**
````diff
diff --git a/docs/_data/r_idioms.toml b/docs/_data/r_idioms.toml
--- a/docs/_data/r_idioms.toml
+++ b/docs/_data/r_idioms.toml
@@ -23,8 +23,8 @@
 "phyloseq::prune_samples" = '`tdata[keep].copy()`'
 "phyloseq::subset_taxa" = '`tdata[:, tdata.var["phylum"] == "Chlamydiae"].copy()`'
 "phyloseq::subset_samples" = '`tdata[tdata.obs["SampleType"] == "Feces"].copy()`'
-"phyloseq::tip_glom" = "not in 0.3"
-"phyloseq::merge_samples" = "not in 0.3"
-"phyloseq::psmelt" = "not in 0.3"
-"phyloseq::plot_tree" = "not in 0.3"
-"phyloseq::plot_net" = "not in 0.3"
+"phyloseq::tip_glom" = "not in 0.4"
+"phyloseq::merge_samples" = "not in 0.4"
+"phyloseq::psmelt" = "not in 0.4"
+"phyloseq::plot_tree" = "not in 0.4"
+"phyloseq::plot_net" = "not in 0.4"
diff --git a/docs/coming_from_r.md b/docs/coming_from_r.md
--- a/docs/coming_from_r.md
+++ b/docs/coming_from_r.md
@@ -2,7 +2,7 @@

 Every public biotapy function names its R equivalent in its docstring. This table is generated
 from those lines each time the docs are built, plus a short list of phyloseq accessors that are
-plain AnnData/TreeData code (`docs/_data/r_idioms.toml`). Rows marked "not in 0.3" have no
+plain AnnData/TreeData code (`docs/_data/r_idioms.toml`). Rows marked "not in 0.4" have no
 biotapy equivalent yet.

 biotapy keeps samples as rows, so `tdata.X` is phyloseq's `otu_table` with
diff --git a/docs/tutorials/phyloseq_analysis.md b/docs/tutorials/phyloseq_analysis.md
--- a/docs/tutorials/phyloseq_analysis.md
+++ b/docs/tutorials/phyloseq_analysis.md
@@ -15,10 +15,10 @@ kernelspec:

 This notebook redoes the sections of phyloseq's
 [analysis vignette](https://github.com/joey711/phyloseq/blob/master/vignettes/phyloseq-analysis.Rmd)
-that biotapy 0.3 covers, on the same three datasets: GlobalPatterns, enterotype and esophagus.
+that biotapy 0.4 covers, on the same three datasets: GlobalPatterns, enterotype and esophagus.
 Each section names the R chunk it follows. Everything runs in biotapy; nothing is read from R.

-**Not in 0.3**, so left out:
+**Not in 0.4**, so left out:

 - `plot_tree` (exploratory tree plots) and `plot_net` (sample networks);
 - correspondence analysis (`ordinate(..., "CCA")`) and DPCoA, with their scree, species and biplot plots;
@@ -170,7 +170,7 @@ bt.pl.ordination(global_patterns, basis="nmds", color="SampleType");

 `distance(esophagus, "bray")`, `"wunifrac"` and `"jaccard"`. phyloseq's `"jaccard"` is vegan's
 quantitative Jaccard; biotapy's is presence/absence, phyloseq's
-`distance(esophagus, "jaccard", binary = TRUE)`. The `betadiver` method `"g"` is not in 0.3.
+`distance(esophagus, "jaccard", binary = TRUE)`. The `betadiver` method `"g"` is not in 0.4.

 ```{code-cell} ipython3
 esophagus = bt.datasets.esophagus()
````
- [x] **Step 4: Run, expect pass** - `53 passed`. `coming_from_r.render()` has 60 rows (58 in
  0.3.0); `grep -rn "not in 0\.3\|biotapy 0\.3\|Not in 0\.3" docs --include=*.md
  --include=*.toml` prints nothing outside `docs/_build` and `docs/generated`.
- [x] **Step 5: Bookkeeping.** Tick 4.D1; log line:
  `- **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.D1 done: the Coming-from-R rows \`MultiAssayExperiment::MultiAssayExperiment\` and \`::intersectColumns\` (to \`bt.io.to_mudata\`) are pinned, and a test checks \`ml\`, \`tl.mmvec\` and \`datasets.biocrust\` add no row; phyloseq calls without an equivalent read "not in 0.4", as do the page sentence and the phyloseq vignette.`
- [x] **Step 6: Gate and commit.**
  ```bash
  git add docs/_data/r_idioms.toml docs/coming_from_r.md docs/tutorials/phyloseq_analysis.md \
    tests/test_coming_from_r.py .knowledge/roadmap/phase-4-ml-multiomics.md .knowledge/log.md
  uvx prek run --all-files
  git commit -m "docs: pin the MultiAssayExperiment rows of Coming from R, mark uncovered calls not in 0.4

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  ```
  Expected: prek passed; `1524 passed, 2 skipped, 98 deselected`; `40 passed, 1584
  deselected`; `build succeeded`.

### Checkpoint D - review slice 4D

- [x] **Review the whole slice** with superpowers:requesting-code-review (opus: the pages make
  claims a user will act on), against data-model-slots, function-shape, module-boundaries,
  pure-by-default, optional-heavy-dependencies, the Phase 4 review focus, and the slice 4D
  design, review focus and global constraints. The reviewer reads each tutorial's prose against
  its rendered outputs and against the source it cites (mmvec's `check_soils.ipynb` and README
  at `88ca33b`; the HMP2 counts; the MGM warning), runs the three page-running tests
  (`-m network` with an absolute `BIOTAPY_DATA_DIR`; `--extra mgm -m mgm`), and checks 4.F4/4.F5
  against every caller. Then a fix pass, one commit per finding, each with a test where a test
  can show it; then a scoped re-review. Record: counts (Critical / Important / Minor), fix range,
  re-review result.
  Recorded: review 0 Critical / 0 Important / 4 Minor new (the per-task reviews' 4.6c Important
  I1, biocrust's successional stages, and their minors carried in); fix pass `4976b14..1842824`
  (11 commits) plus follow-up `1842824..38f9379` (3 commits); scoped re-review: 10 of 11
  addressed, item 2 partial and two stale design paragraphs, all fixed in the follow-up.
- [x] **Task 4.7** (below), after the fix pass.
- [x] **Run the gates** on the committed tree (`git status --short` empty), every command with
  absolute `BIOTAPY_DATA_DIR` and `HF_HOME`: the slice gate's seven counts (default, `golden or
  network`, docs, `-m torch`, `-m mgm`, prek, `knowledge_stale.sh --against HEAD` at `28
  current, 0 stale`); `coverage report` on `datasets/_biocrust.py` and `_core/_matrix.py`
  (100%); Hypothesis seeds 1-3 on the 4.F4/4.F5 test files. Confirm `~/.cache/biotapy` and
  `docs/pooch/` do not exist.
  - At `bc65712` (clean tree): prek passed; `pytest -q -W error::UserWarning` 1524 passed, 2
    skipped, 98 deselected; `-m "golden or network"` 40 passed; sphinx `-W` succeeded;
    `-m torch` 25 passed; `-m mgm` 15 passed; asv check no problems; coverage of
    `_core/_matrix.py` and `datasets/_biocrust.py` 100%; Hypothesis seeds 1-3 on
    `tests/core/test_matrix.py` 41 passed each; `knowledge_stale.sh` 28 current, 0 stale;
    `~/.cache/biotapy` and `docs/pooch/` absent. CI on PR #31 failed the three
    `hatch-test.py3.14-pre` jobs: pre-release pandas' Pandas4Warning on `future.infer_string`
    (read by anndata) tripped a blanket `simplefilter("error")` in a biocrust test, narrowed to
    `UserWarning` in `8353fe3`; then 22/22 green (run 38083140914), merged as `b0066c3`.
- [x] **Push** `phase-4d`, open the PR and merge-commit it on green, only after the user approves
  that push (decision 42; decision 20 gives no standing approval for Phase 4). CI must be green:
  every hatch-test job, `lint`, `import-without-extras`, `network` (its log shows `40 passed`:
  the leak-free and multi-omics pages and `biocrust`), `r-bridge`, `ml-extras` (its `-m mgm`
  step `15 passed`, the embedding page among them), `docs` (seven `Executed notebook` lines,
  `tutorials/leak_free_cv.md` and `tutorials/multiomics.md` among them). Record the PR number,
  merge commit, Test run id and the docs job id: exit-gate items 1 and 3's evidence, ticked in
  4.D2 Step 8. If a page-running test fails only on a number, stop and report the measured
  values; do not edit the constants without the user (R11.5).
- [x] **After the merge**, once Read the Docs has built `master`: `curl -s -o /dev/null -w
  "%{http_code}\n" https://biotapy.readthedocs.io/en/latest/tutorials/leak_free_cv.html` and the
  same for `tutorials/multiomics.html` and `tutorials/embeddings.html` print `200`, and the
  leak-free page shows its table. If a page is missing, read the build log on readthedocs.org
  and report it.
- [x] **Ask the user, in one message:** to review slice 4D; to approve pushing
  `release-0.4.0`, opening its PR and merge-committing it on green (4.D2 Steps 3-9); and to
  confirm that pypi.org project `biotapy` still lists the trusted publisher (owner `pedrocr83`,
  repository `biotapy`, workflow `release.yaml`, environment `pypi`) that published 0.3.0. Stop
  until the user answers.

### Task 4.7: Knowledge

Run as Checkpoint D's knowledge step, after the fix pass, so it documents the final code. Most
concept text was written by the tasks (`core.md` in 4.F4/4.F5, `datasets.md` in 4.6c); this task
adds what spans the slice and re-stamps what the slice touched.

**Files:** modify `.knowledge/modules/ml.md`, `.knowledge/modules/tl.md`,
`.knowledge/decisions/embedding-plugins.md`, `.knowledge/roadmap/phase-4-ml-multiomics.md`,
`.knowledge/log.md`, and `generated`/`commit` only on the concepts Step 2 lists.
**Not touched:** every contract's body (nothing they state changed: no slot, key, signature
shape, layer or import rule moved; `datasets` importing `io` is the existing layer order);
`optional-heavy-dependencies.md` (no extra changed); `multiomics-as-mudata.md` (the names and
the intersection rule hold, and `biocrust` follows them); `add-a-function.md`'s steps; `verified`
anywhere (R12.3).
**Interfaces:**
- Consumes: the final slice tree; `bash scripts/knowledge_stale.sh --against HEAD`.
- Produces: a bundle with `28 current, 0 stale, 12 uncheckable`.

- [x] **Step 1: Concepts the slice changes.** Re-check each sentence against the code after the
  fix pass (a fix may move a number). `embedding-plugins` becomes `stable` only because the user
  confirmed it with slice 4C (Checkpoint C's last box, ticked by approving this plan); if they
  have not, leave it `draft` and say so:
````diff
diff --git a/.knowledge/decisions/embedding-plugins.md b/.knowledge/decisions/embedding-plugins.md
--- a/.knowledge/decisions/embedding-plugins.md
+++ b/.knowledge/decisions/embedding-plugins.md
@@ -3,10 +3,10 @@ type: Decision
 title: Embedding models are plugins found through entry points
 description: ml.embed(adata, model) loads the callable a package registers under model in the entry-point group biotapy.embeddings, checks that it returned a finite 2-D float array with one row per sample, and returns it or writes obsm["X_<model>"]; biotapy's own MGM is registered the same way, its weights downloaded, never bundled.
 tags: [ml, plugins, api]
-status: draft
+status: stable
 paths: ["src/biotapy/ml/_embed.py", "pyproject.toml"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T16:44:46Z }
-commit: 31aa11d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T17:11:09Z }
+commit: b904fa6
 sources:
   - id: spec
     resource: ../../plan.md
diff --git a/.knowledge/modules/ml.md b/.knowledge/modules/ml.md
--- a/.knowledge/modules/ml.md
+++ b/.knowledge/modules/ml.md
@@ -6,8 +6,8 @@ resource: /src/biotapy/ml/
 paths: ["src/biotapy/ml/**"]
 tags: [ml, scikit-learn, torch, plugins]
 status: stable
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T16:44:46Z }
-commit: 31aa11d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T17:11:09Z }
+commit: b904fa6
 ---

 # Responsibility
@@ -143,6 +143,13 @@ runs MGM against its own embeddings (`tests/data/mgm`), the end-to-end GlobalPat
 and `embed`'s docstring example, as the job's last step does. The pseudocount
 warning's text is unit-tested in `tests/core/test_composition.py`.

+The docs build runs `docs/tutorials/leak_free_cv.md` on HMP2;
+`-m network` runs its cells again and compares them with `LEAK_FREE_CV_AUC`
+(`tests/ml/test_transformers.py`). `docs/tutorials/embeddings.md` is not run by
+the docs build (no torch): `-m mgm` runs its code, and a default-run test checks
+that it quotes `GLOBAL_PATTERNS_LEFT_OUT` and `GLOBAL_PATTERNS_SHAPE`
+(`tests/ml/test_mgm.py`).
+
 # Gotchas

 - `check_estimator` cannot run under `-W error::UserWarning`, the project's
@@ -216,6 +223,25 @@ warning's text is unit-tested in `tests/core/test_composition.py`.
   features on purpose: summing them moves its rank in one sample, so a reader
   that does not sum features of one genus fails the parity test.
   `tests/mgm/export_reference.py:table`.
+- **The two tutorials quote numbers that tests pin.** The prose of
+  `docs/tutorials/leak_free_cv.md` and the quoted outputs of
+  `docs/tutorials/embeddings.md` are constants in `tests/ml/test_transformers.py`
+  and `tests/ml/test_mgm.py`; the root `tests/conftest.py:run_page` fixture runs
+  a page's code cells and Python blocks, so the network and `mgm` tests check
+  the page's own code. A scikit-learn, torch or data change that moves a number
+  fails one of them: update the constant and the page together.
+- **The leak on HMP2 is small for the prevalence filter and large for a step
+  that reads the labels**: the filter outside the pipeline moves the mean AUC
+  by 0.001 (0.527 against 0.526), `SelectKBest(k=20)` by 0.147 on the real
+  labels and by 0.223 on shuffled ones (0.544 against 0.767). The honest
+  pipeline is near chance: a first stool sample's species barely separate IBD
+  from non-IBD there. The tutorial is about leakage, not a classifier.
+- **GlobalPatterns writes *Candidatus* genera as one word** (Greengenes:
+  `CandidatusPelagibacter`), MGM's vocabulary as `Candidatus_Pelagibacter`. Of
+  the 96 GlobalPatterns genus-level features MGM leaves out, 28 are such names,
+  23 of which MGM's vocabulary holds with the underscore. MGM's own regex reads
+  the name as written, so biotapy leaves them out as MGM would (Phase 4
+  decision 30); `test_embeds_global_patterns_end_to_end` pins the 96.
 - anndata 0.13 lists `X` as `layers[None]`, so `list(adata.layers)` holds
   `None` even when no layer was added; `to_torch`'s missing-layer error does
   not list the layers. `_torch.py:_table`.
diff --git a/.knowledge/modules/tl.md b/.knowledge/modules/tl.md
--- a/.knowledge/modules/tl.md
+++ b/.knowledge/modules/tl.md
@@ -5,8 +5,8 @@ description: Diversity, ordination, PERMANOVA and mmvec over AnnData/TreeData/Mu
 resource: /src/biotapy/tl/
 paths: ["src/biotapy/tl/**"]
 tags: [tl, diversity, ordination]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T03:45:07Z }
-commit: 1c5d1ae
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T17:11:09Z }
+commit: b904fa6
 status: stable
 ---

@@ -116,7 +116,11 @@ cache: `BIOTAPY_DATA_DIR=<cache> uv run --group test pytest -m golden tests/tl -
 - `mmvec` has no R golden (`R equivalent: none`); its tests check shape, row
   centring, that a metabolite ranks highest for the microbe it follows, and
   the checks above, not scikit-bio's learned values, which depend on the seed
-  (`tests/tl/test_mmvec.py`). The fit holds dense microbes x metabolites
+  (`tests/tl/test_mmvec.py`). On real data, a network test runs
+  `docs/tutorials/multiomics.md` on `bt.datasets.biocrust()` and repeats the
+  check mmvec's own soil example makes: all 13 metabolites it lists for
+  `rplo 1 (Cyanobacteria)` rank above zero (`MICROCOLEUS_METABOLITES`), with
+  scikit-bio's L-BFGS fit and 3 dimensions where the example used 1. The fit holds dense microbes x metabolites
   arrays, 8 bytes each, on top of the dense inputs. An all-zero sample or
   feature in either modality raises rather than being dropped.
 - All-zero samples give scikit-bio's values, with no custom mapping: two of
````
- [x] **Step 2: Concepts the diff only touches.** `bash scripts/knowledge_stale.sh --against
  HEAD` on the fix pass's last commit lists 12 stale on the prototype: `phase-0-foundation`,
  `phase-1-core`, `phase-2-function`, `phase-4-ml-multiomics`, `core`, `datasets`,
  `data-model-slots`, `engine-parity`, `function-shape`, `module-boundaries`, `r-golden-parity`,
  `add-a-function`. Read each against `git diff 28de1ad..HEAD -- <its paths>`: the `_core`
  change is a dtype rule, `datasets` gains a loader that follows the existing ones, `tests/`
  gains tests and a fixture, `docs/` gains pages; nothing any of them states became false
  (r-golden-parity's rows are unchanged: the mmvec check is against mmvec's example, not an R
  golden, and `tl.md` records it). Bump only `generated` and `commit` on each (and on `ml`, `tl`,
  `embedding-plugins` with Step 1): the model that runs the task, `date -u +%FT%TZ`, and the
  parent's short hash.
- [x] **Step 3: Roadmap.** Tick `4.7 Knowledge` in the checklist; record Checkpoint D's review
  counts, fix range and re-review in its first box.
- [x] **Step 4: Log**, first under the slice heading:
  ```markdown
  - **Update**: [ml](modules/ml.md): Verification names the two tutorials and the tests that pin what they quote; new gotchas: the tutorials' numbers live in test constants checked by running the pages' code (`tests/conftest.py:run_page`), the size of the leak on HMP2, and GlobalPatterns' one-word *Candidatus* genera that MGM leaves out.
  - **Update**: [embedding-plugins](decisions/embedding-plugins.md) is `stable`: the user confirmed it with slice 4C.
  - **Update**: [tl](modules/tl.md): `mmvec`'s gotcha names the network test that runs the multi-omics tutorial against mmvec's own soil-example check.
  - **Verification**: re-checked against the slice 4D diff and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [phase-1-core](roadmap/phase-1-core.md), [phase-2-function](roadmap/phase-2-function.md), [core](modules/core.md), [datasets](modules/datasets.md), [data-model-slots](contracts/data-model-slots.md), [engine-parity](contracts/engine-parity.md), [function-shape](contracts/function-shape.md), [module-boundaries](contracts/module-boundaries.md), [r-golden-parity](contracts/r-golden-parity.md), [add-a-function](playbooks/add-a-function.md).
  - **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) task 4.7 done.
  ```
- [x] **Step 5: Gate and commit.**
  ```bash
  git add .knowledge
  uvx prek run --all-files
  git commit -m "docs(knowledge): document the slice 4D tutorials and refresh concepts

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  bash scripts/knowledge_stale.sh --against HEAD       # 28 current, 0 stale, 12 uncheckable
  ```
  Expected: prek passed; `1524 passed, 2 skipped, 98 deselected`;
  `40 passed, 1584 deselected`; `build succeeded`.

### Task 4.D2: Release 0.4.0

Follows [cut-a-release](/playbooks/cut-a-release.md) and Task 3.15. Publishing to PyPI is
irreversible. Every outward step waits for the user's explicit approval of that step (R13.3).

**Files:**
- Modify (release commit): `pyproject.toml` (version), `CHANGELOG.md`, `README.md`.
- Modify (knowledge commit): `.knowledge/playbooks/cut-a-release.md`, the concepts the version
  bump makes stale, `.knowledge/roadmap/phase-4-ml-multiomics.md`, `.knowledge/log.md`.
- Modify (after the upload, branch `close-phase-4`): `.knowledge/roadmap/phase-4-ml-multiomics.md`,
  `.knowledge/roadmap/phase-5-beyond.md`, `.knowledge/roadmap/index.md`, `.knowledge/log.md`.
**Not touched:** any `src/`, `tests/` or `docs/` file (the version reaches the docs through
package metadata; `docs/changelog.md` includes `CHANGELOG.md`); the workflows; `uv.lock`
(git-ignored).
**Interfaces:**
- Consumes: `master` after Checkpoint D's merge; `.github/workflows/release.yaml` (on `release:
  published`: `uv build`, then `pypa/gh-action-pypi-publish` v1.14.2 in environment `pypi`).
- Produces: biotapy 0.4.0 on PyPI with the extras `mgm`, `r`, `torch` and the entry point
  `biotapy.embeddings:mgm`; Phase 4 `phase_state: done`.

- [x] **Step 1: Approval in hand.** Checkpoint D's last box asked for the slice review, the
  release branch push/PR/merge and the PyPI publisher check. Do not start until all three are
  answered.
- [x] **Step 2: Check the publish action.** `gh api
  repos/pypa/gh-action-pypi-publish/releases/latest --jq .tag_name` prints `v1.14.2` (checked
  2026-10-10; `release.yaml` pins `v1.14.2`). If newer, stop and report it.
- [x] **Step 3: Branch.** `git switch master && git pull --ff-only && git switch -c
  release-0.4.0`.
- [x] **Step 4: Version and changelog.** `## [Unreleased]` is empty (no 4A-4D PR wrote an
  entry). Write the 0.4.0 entries from `git log v0.3.0..master --no-merges --oneline` and the
  diff below, with `date -u +%F`; leave a new empty `## [Unreleased]` above. Each "Changed" line
  is a change to a function 0.3.0 shipped: `5978028` (the pseudocount warning's suffix), `8e1567e`
  (`filter_features`' `TypeError`), 4.F4/4.F5, 4.D1. If the fix pass changed another 0.3.0
  function, add its line; if decision 36 or 37 was declined, drop `bt.fn.func_glom` or the
  `biocrust` lines.
````diff
diff --git a/CHANGELOG.md b/CHANGELOG.md
--- a/CHANGELOG.md
+++ b/CHANGELOG.md
@@ -10,6 +10,53 @@ and this project adheres to [Semantic Versioning][].

 ## [Unreleased]

+## [0.4.0] - 2026-10-10
+
+### Added
+
+- `bt.io.to_mudata`: combine data types measured on the same samples into one
+  `MuData`, keeping the samples every modality has (in the first one's order,
+  with a warning naming how many each lost). Modality names: `taxa`,
+  `function`, `function_by_taxon`, `metabolites`, `host`.
+- `bt.tl.mmvec`: which metabolites go with which microbe (mmvec, Morton et
+  al. 2019), fitted by scikit-bio, as a microbes x metabolites table of
+  row-centred log probabilities.
+- `bt.ml.PrevalenceFilter` and `bt.ml.CLR`: scikit-learn transformers, so a
+  prevalence filter and CLR are fitted inside each cross-validation fold.
+  They pass scikit-learn's estimator checks and give `bt.pp.filter_features`'s
+  and `bt.pp.clr`'s results. Relative abundance is scikit-learn's own
+  `Normalizer(norm="l1")`.
+- `bt.ml.to_torch`: a PyTorch dataset over a table's samples, densifying one
+  row at a time, with integer or float labels from an `obs` column. It needs
+  the new extra `biotapy[torch]` (`torch>=2.9`).
+- `bt.ml.embed`: one embedding per sample from a pretrained model, returned or
+  stored in `obsm["X_<model>"]`. Models are plugins: a package registers one
+  in the entry-point group `biotapy.embeddings`, and biotapy checks what it
+  returns.
+- MGM, the Microbial General Model (Zhang et al. 2026), as the first such
+  model, behind the new extra `biotapy[mgm]` (`torch>=2.9`,
+  `transformers>=5`). Its pretrained weights (MIT) are downloaded once from
+  the `microformer-mgm` 0.5.8 wheel on PyPI and checked against its SHA-256;
+  biotapy's embeddings are within 2e-6 of MGM's own code.
+- `bt.datasets.biocrust`: mmvec's soil example, microbes and metabolites of
+  a desert biocrust after wetting, downloaded and cached on first use.
+- Guide pages for multi-omics and machine learning, and three tutorials:
+  multi-omics with mmvec, leak-free cross-validation on the HMP2 cohort, and
+  MGM embeddings of GlobalPatterns.
+
+### Changed
+
+- The pseudocount warning of `bt.pp.clr` and `bt.pp.philr` ends with the
+  step that gave it, `(pp.clr)` or `(pp.philr)`.
+- `bt.pp.filter_features` raises `TypeError` naming `min_prevalence` or
+  `min_total` when the threshold is a bool or not a number.
+- `bt.pp.tax_glom`, `bt.fn.func_glom` and `bt.pl.bar`'s `fill` sum a table
+  of small integers or bools in 64 bits, as NumPy's `sum` does: an `int8`
+  table no longer wraps around, and a bool table is counted instead of
+  saturating at `True`. An `int32` table now gives `int64` sums.
+- The Coming-from-R page marks the phyloseq calls biotapy does not cover yet
+  "not in 0.4".
+
 ## [0.3.0] - 2026-10-07

 ### Added
diff --git a/pyproject.toml b/pyproject.toml
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -4,7 +4,7 @@ requires = [ "hatchling>=1.27" ]

 [project]
 name = "biotapy"
-version = "0.3.0"
+version = "0.4.0"
 description = "mia-style microbiome toolkit for Python on AnnData/TreeData"
 readme = "README.md"
 license = "BSD-3-Clause"
````
- [x] **Step 5: README** (PyPI's project page): Status, Next and Installation only.
````diff
diff --git a/README.md b/README.md
--- a/README.md
+++ b/README.md
@@ -14,9 +14,9 @@ against R on real data.

 ## Status

-**biotapy 0.3 is an early release.** The API can still change between minor versions.
+**biotapy 0.4 is an early release.** The API can still change between minor versions.

-What 0.3 does (full signatures in the [API reference][api]):
+What 0.4 does (full signatures in the [API reference][api]):

 - **Readers**: `bt.io.read_biom` (BIOM 1.0/2.1), `bt.io.read_qiime2`
   (`.qza` artifacts, no QIIME 2 install needed), `bt.io.read_dada2`
@@ -25,14 +25,16 @@ What 0.3 does (full signatures in the [API reference][api]):
   and 4 profiles) all read into one `TreeData`. `bt.io.read_humann` and
   `bt.io.read_picrust2` read function tables into a `MuData` with a community
   and a per-taxon modality; `bt.io.read_picrust2_traits` reads PICRUSt2's
-  per-ASV gene copy numbers.
+  per-ASV gene copy numbers. `bt.io.to_mudata` combines data types measured
+  on the same samples into one `MuData`.
 - **Writer**: `bt.io.write_biom` writes a BIOM 2.1 or 1.0 table back out.
 - **Datasets**: `bt.datasets.toy` and `bt.datasets.toy_humann` (in-memory,
   for examples and tests), `bt.datasets.global_patterns`,
   `bt.datasets.enterotype` and `bt.datasets.esophagus` (phyloseq's example
-  datasets), `bt.datasets.hmp2` (the HMP2 inflammatory bowel disease cohort)
-  and `bt.datasets.enzyme` (the ENZYME EC hierarchy), downloaded and cached
-  on first use.
+  datasets), `bt.datasets.hmp2` (the HMP2 inflammatory bowel disease cohort),
+  `bt.datasets.biocrust` (mmvec's microbes and metabolites example) and
+  `bt.datasets.enzyme` (the ENZYME EC hierarchy), downloaded and cached on
+  first use.
 - **Preprocessing**: `bt.pp.relative`, `bt.pp.tax_glom`,
   `bt.pp.filter_features`, `bt.pp.filter_samples` and `bt.pp.rarefy`, and
   the compositional transforms `bt.pp.clr` and `bt.pp.philr` (checked against
@@ -41,18 +43,23 @@ What 0.3 does (full signatures in the [API reference][api]):
   (checked against HUMAnN's own output), `bt.fn.contributions` and
   `bt.fn.functional_redundancy`.
 - **Tools**: `bt.tl.alpha`, `bt.tl.beta`, `bt.tl.unifrac`, `bt.tl.pcoa`,
-  `bt.tl.nmds` and `bt.tl.permanova`, each checked against R on real data.
+  `bt.tl.nmds` and `bt.tl.permanova`, each checked against R on real data,
+  and `bt.tl.mmvec`, which relates microbes to metabolites.
 - **Differential abundance**: `bt.da.linda` and `bt.da.ancombc2` in Python,
   `bt.da.aldex2` and `bt.da.maaslin3` through R (the `r` extra below), all
   returning one result table and checked against their R packages, and
   `bt.da.consensus`, which reports where the methods agree.
+- **Machine learning**: `bt.ml.PrevalenceFilter` and `bt.ml.CLR`,
+  scikit-learn transformers fitted inside each cross-validation fold;
+  `bt.ml.to_torch`, a PyTorch dataset (the `torch` extra below); and
+  `bt.ml.embed`, sample embeddings from pretrained models that plugins
+  register, with MGM included (the `mgm` extra below).
 - **Plots**: `bt.pl.bar`, `bt.pl.richness`, `bt.pl.ordination`,
   `bt.pl.scree`, `bt.pl.heatmap`, `bt.pl.contributions` and
   `bt.pl.consensus`.

-Next, in 0.4: multi-omics conventions on MuData, leak-free scikit-learn
-transformers, a PyTorch loader and an interface for embedding models. See the
-[roadmap][roadmap]; no dates are promised.
+Next: what users and issues ask for. See the [roadmap][roadmap]; no dates are
+promised.

 ## Installation

@@ -89,6 +96,19 @@ pip install 'biotapy[r]'
 rpy2 is GPL-2.0-or-later and the R packages carry their own licences; biotapy
 does not ship any of them.

+`bt.ml.to_torch` needs PyTorch, and MGM embeddings need PyTorch and
+transformers, through two extras:
+
+```bash
+pip install 'biotapy[torch]'   # bt.ml.to_torch
+pip install 'biotapy[mgm]'     # bt.ml.embed(..., "mgm")
+```
+
+On Linux, pip installs PyPI's torch, which brings CUDA libraries; for a
+CPU-only torch, run `pip install torch --index-url
+https://download.pytorch.org/whl/cpu` first. torch publishes no wheel for
+Intel macOS.
+
 On Python 3.14, the `biom-format` dependency has no wheels yet, so it is built
 from source and needs a C compiler until biom-format publishes 3.14 wheels
 ([biocore/biom-format#1004][biom-format-1004]).
````
- [x] **Step 6: Build check.**
  - `rm -rf dist && uv build && uvx twine check --strict dist/*`: both files `PASSED`.
  - `unzip -p dist/biotapy-0.4.0-py3-none-any.whl 'biotapy-0.4.0.dist-info/METADATA' | grep -E
    "^(Version|Provides-Extra|Requires-Dist: (torch|transformers|rpy2))"` prints `Version:
    0.4.0`, `Provides-Extra: mgm`, `Requires-Dist: torch>=2.9; extra == 'mgm'`, `Requires-Dist:
    transformers>=5; extra == 'mgm'`, `Provides-Extra: r`, `Requires-Dist: rpy2>=3.6.8; extra ==
    'r'`, `Provides-Extra: torch`, `Requires-Dist: torch>=2.9; extra == 'torch'`.
  - `unzip -p dist/biotapy-0.4.0-py3-none-any.whl 'biotapy-0.4.0.dist-info/entry_points.txt'`
    prints `[biotapy.embeddings]` and `mgm = biotapy.ml._mgm:embed`.
  - `tar tzf dist/biotapy-0.4.0.tar.gz | grep -E
    "biotapy-0.4.0/conftest.py|tests/conftest.py|docs/tutorials/(leak_free_cv|embeddings|multiomics).md|tests/datasets/test_biocrust.py"`
    prints those six paths; `tar tzf dist/biotapy-0.4.0.tar.gz | grep -c "tests/mgm\|benchmarks"`
    prints `0`.
  - The sdist's own tests: extract it under `<scratchpad>/sdist`, then inside it, with absolute
    `BIOTAPY_DATA_DIR`, `XDG_CACHE_HOME`, `HF_HOME` and `MPLCONFIGDIR`, `uv run --group test
    pytest -q -p no:cacheprovider` -> `1405 passed, 2 skipped, 96 deselected` (plus the fix
    pass's tests; the sdist leaves out `tests/test_ci.py`, the knowledge tests, `tests/humann`,
    `tests/mgm`, `tests/r`). A relative `XDG_CACHE_HOME` fails
    `test_without_biotapy_data_dir_the_cache_is_pooch_s_per_user_directory` (an unnormalised
    `..` in `pooch.os_cache`); that is the command, not the package.
  - `rm -rf dist <scratchpad>/sdist`.
- [x] **Step 7: Gate and commit.**
  ```bash
  git add pyproject.toml CHANGELOG.md README.md
  uvx prek run --all-files
  git commit -m "chore: release 0.4.0

  Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>"
  uv run python -c "import biotapy; print(biotapy.__version__)"   # 0.4.0
  ```
  Expected: prek passed; `1524 passed, 2 skipped, 98 deselected`; `40 passed, 1584
  deselected`; `build succeeded`; `25 passed, 1599 deselected`; `15 passed, 1609 deselected`;
  `knowledge_stale.sh` `23 current, 5 stale` (Step 8 re-stamps them).
- [x] **Step 8: Knowledge, second commit.**
  - `bash scripts/knowledge_stale.sh --against HEAD` lists `phase-0-foundation`,
    `embedding-plugins`, `module-boundaries`, `tree-access`, `cut-a-release` (their `paths` hold
    `pyproject.toml` or `CHANGELOG.md`). Check each against the bump, which changes nothing they
    state, and bump `generated` and `commit`.
  - `cut-a-release` step 3 gains the check this release added:
````diff
diff --git a/.knowledge/playbooks/cut-a-release.md b/.knowledge/playbooks/cut-a-release.md
--- a/.knowledge/playbooks/cut-a-release.md
+++ b/.knowledge/playbooks/cut-a-release.md
@@ -5,8 +5,8 @@ description: Bump the version, move the changelog entry, tag, publish a GitHub r
 tags: [release, workflow]
 status: stable
 paths: ["pyproject.toml", "CHANGELOG.md", ".github/workflows/release.yaml"]
-generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T16:44:46Z }
-commit: 31aa11d
+generated: { by: claude-code/claude-opus-5-5, at: 2026-10-10T17:24:15Z }
+commit: 20151ac
 sources:
   - id: trusted-publishing
     resource: https://docs.pypi.org/trusted-publishers/
@@ -39,8 +39,11 @@ A phase exit gate asks for a release (0.0.1 in Phase 0, 0.1 in Phase 1, ...).
    extracted sdist (`uv build --sdist`, `tar xzf`, then `pytest` inside it); it
    must pass with no errors. Check the wheel's `METADATA` too: `Version`, every
    runtime `Requires-Dist`, and a `Provides-Extra` line for each extra (0.3.0
-   added `r`, slice 4B `torch`, slice 4C `mgm`). The sdist must ship the root `conftest.py`
-   (`pyproject.toml` `build.targets.sdist.include`): its tests need the marker hook.
+   added `r`, slice 4B `torch`, slice 4C `mgm`), and the wheel's `entry_points.txt`
+   (0.4.0 added `[biotapy.embeddings] mgm = biotapy.ml._mgm:embed`). The sdist must ship the root `conftest.py`
+   (`pyproject.toml` `build.targets.sdist.include`): its tests need the marker hook. Run the sdist's tests with
+   absolute cache paths (`BIOTAPY_DATA_DIR`, `XDG_CACHE_HOME`): `tests/core/test_download.py` compares
+   `pooch.os_cache` paths, which a relative `XDG_CACHE_HOME` leaves unnormalised.
 4. With explicit user approval for each (rules.md R13.3), tag the merged commit
    and push only the tag:
    ```bash
````
  - The roadmap: tick the exit gate's items 1 and 3 with Checkpoint D's PR, Test run and docs job
    (item 2 was ticked by the plan commit), Checkpoint D's boxes and checklist line, and Steps 1-8
    here; `<...>` are the run's values:
    ```diff
    -- [ ] Leak-free CV example executed in docs: `docs/tutorials/leak_free_cv.md`
    -  runs in the docs CI job (design note 9).
    +- [x] Leak-free CV example executed in docs: `docs/tutorials/leak_free_cv.md`
    +  runs in the docs CI job (design note 9). Proven on PR #<n>: docs job <id>
    +  executed `tutorials/leak_free_cv.md`.
    @@
    -- [ ] All Phase 1-3 gates still green.
    +- [x] All Phase 1-3 gates still green. Proven on PR #<n>: every job of Test run <id> green.
    ```
  - Log, under `## <date> (release 0.4.0)` at the top:
    ```markdown
    ## <date> (release 0.4.0)
    - **Update**: [cut-a-release](playbooks/cut-a-release.md) step 3: check the wheel's `entry_points.txt` too, and run the sdist's tests with absolute cache paths.
    - **Verification**: re-checked against the 0.4.0 version bump and bumped only: [phase-0-foundation](roadmap/phase-0-foundation.md), [embedding-plugins](decisions/embedding-plugins.md), [module-boundaries](contracts/module-boundaries.md), [tree-access](contracts/tree-access.md).
    - **Update**: [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md) ticks the exit gate's notebook and "all gates green" items (docs job <id> and Test run <id> on PR #<n>), Checkpoint D and Task 4.D2 Steps 1-8.
    ```
  - Commit `docs(knowledge): refresh concepts for the 0.4.0 release` (with the attribution
    line), staging the five concepts, `phase-4-ml-multiomics.md` and `log.md`; prek; the default
    run (`1524 passed, 2 skipped, 98 deselected`); `knowledge_stale.sh --against HEAD` -> `28
    current, 0 stale, 12 uncheckable`.
- [x] **Step 9: PR and merge** (approved at Checkpoint D): push `release-0.4.0`, open the PR,
  wait for green CI (every job, `ml-extras`, `r-bridge` and `docs` included), merge with a merge
  commit.
- [x] **Step 10: STOP AND ASK.** Tagging and releasing publishes to PyPI, which cannot be undone.
  Ask the user, in one message, for explicit approval to: tag the merged `master` commit
  `v0.4.0`; push only that tag; create the GitHub release `0.4.0`, which triggers the upload; and
  afterwards push `close-phase-4` (Step 13), open its PR and merge-commit it on green. Wait for a
  yes that names this release.
- [x] **Step 11: Tag and release**, only after Step 10's approval:
  ```bash
  git switch master && git pull --ff-only
  git log -1 --format=%s   # Merge pull request #<n> from pedrocr83/release-0.4.0
  git tag v0.4.0
  git push origin v0.4.0
  awk '/^## \[0.4.0\]/{on=1;next} /^## \[/{on=0} on' CHANGELOG.md > <scratchpad>/notes-0.4.0.md
  gh release create v0.4.0 --title "0.4.0" --notes-file <scratchpad>/notes-0.4.0.md
  gh run list --workflow release.yaml --limit 3 --json databaseId,headBranch,status
  ```
  Pick the run whose `headBranch` is `v0.4.0`, then `gh run watch <id> --exit-status`. If the
  `pypi` environment asks for a reviewer, the user approves it in the run. If it fails before
  the upload: fix `master`, then delete the release and tag (`gh release delete v0.4.0
  --cleanup-tag`) and re-tag, each only with the user's approval. If it fails after the upload,
  report it: the version is spent, and the fix is 0.4.1.
- [x] **Step 12: Verify** (playbook Verification):
  ```bash
  curl -s https://pypi.org/pypi/biotapy/json | python3 -c "import json,sys; print(json.load(sys.stdin)['info']['version'])"
  curl -s https://pypi.org/pypi/biotapy/0.4.0/json | python3 -c "import json,sys; print(json.load(sys.stdin)['info']['provides_extra'])"
  uv run --no-project --with biotapy==0.4.0 python -c "import biotapy as bt; from importlib.metadata import entry_points; print(bt.__version__, bt.ml.__all__, [e.value for e in entry_points(group='biotapy.embeddings')], 'biocrust' in bt.datasets.__all__)"
  ```
  Expected: `0.4.0`; `['mgm', 'r', 'torch']`; `0.4.0 ['CLR', 'PrevalenceFilter', 'embed',
  'to_torch'] ['biotapy.ml._mgm:embed'] True` (run with `BIOTAPY_DATA_DIR` set; nothing is
  downloaded).
  - Release 0.4.0: PR #32 (22/22 green, run 38084175167) merged as `5cdd6ee`; tag `v0.4.0`;
    GitHub release https://github.com/pedrocr83/biotapy/releases/tag/v0.4.0; release workflow
    run 38084530873 uploaded to PyPI. PyPI shows `0.4.0`, `provides_extra` `['mgm', 'r', 'torch']`;
    a fresh `uv run --with biotapy==0.4.0` prints `0.4.0 ['CLR', 'PrevalenceFilter', 'embed',
    'to_torch'] ['biotapy.ml._mgm:embed'] True`.
- [x] **Step 13: Close Phase 4** (approved in Step 10), as `close-phase-3` did:
  - `git switch master && git pull --ff-only && git switch -c close-phase-4`;
  - `roadmap/phase-4-ml-multiomics.md`: `phase_state: in-progress` becomes `done`; tick Steps
    9-13 and the checklist's `4.D2 Release 0.4.0`; bump `generated` and `commit`;
  - `roadmap/phase-5-beyond.md`: `phase_state: not-started` becomes `in-progress` (decision 45:
    its own rule still gives no item a task until an issue or the user asks); bump `generated`
    and `commit`;
  - `roadmap/index.md`: Phase 5 moves to "# Active phase" with `**phase_state: in-progress**`;
    Phase 4 moves to "# Phases" after Phase 3 as `**phase_state: done** (0.4.0 on PyPI, <date>)`;
  - log, under `## <date> (release 0.4.0)`: `- **Update**: Phase 4 closed after biotapy 0.4.0
    reached PyPI (tag v0.4.0, release workflow run <id>). [phase-4-ml-multiomics](roadmap/phase-4-ml-multiomics.md)
    is \`phase_state: done\` with every Task 4.D2 step and exit-gate item ticked;
    [phase-5-beyond](roadmap/phase-5-beyond.md) is \`phase_state: in-progress\`, its items still
    waiting for an issue or a user request; the [roadmap index](roadmap/index.md) lists it as active.`
  - gates (`uvx prek run --all-files`; `uv run --group test pytest -q -W error::UserWarning`);
    commit `docs(knowledge): close Phase 4 after the 0.4.0 release` (with the attribution line);
    push, PR, merge-commit on green.
- [ ] **Step 14: After 0.4.** Two upstream reports wait for the user's approval of their text
  (GitHub actions on other projects): the h5mu tree loss (decision 4) and the torch first-`tanh`
  race (decision 35, superseded: the issue exists upstream as pytorch/pytorch#188792). Ask once,
  after the release; do nothing without a yes.

# Decisions for the user
Each changes a contract, rule, dependency, CI or a roadmap signature, or is
a judgement call. Recommended answer first.

1. **Slices and order:** 4A no new dependency (4.1, 4.2, 4.A0, 4.3) -> 4B
   extra `torch` and `ml.to_torch` -> 4C plugin interface and MGM -> 4D
   docs and release. mmvec moves into 4A because scikit-bio 0.7.4 ships it
   (design note 1).
2. **`io.to_mudata(modalities: Mapping[str, AnnData]) -> MuData`**, not
   `to_mudata(**modalities)`; keeps the samples every modality has (first
   modality's order, one warning), copies each modality, writes no
   provenance, documents but does not enforce the names `taxa`, `function`,
   `function_by_taxon`, `metabolites`, `host`. Confirm the new decision
   concept [multiomics-as-mudata](/decisions/multiomics-as-mudata.md)
   (`draft` until you do) (design note 2).
3. **Name:** keep `bt.io.to_mudata` although `mudata.to_mudata` means
   something else (splitting one AnnData). Alternative: `bt.io.combine`.
4. **Trees in h5mu:** document the loss and pin it with a test; no biotapy
   writer (reversed 2026-10-11 at the user's request: `bt.io.write_h5mu`,
   task 5.1 in [phase-5-beyond](phase-5-beyond.md)). Filing an issue with mudata or treedata is a GitHub action on
   another project: recommended after 0.4, with your approval of the text.
5. **mmvec route:** wrap `skbio.stats.ordination.mmvec` (scikit-bio 0.7.4,
   NumPy/SciPy, no TensorFlow); return the row-centred ranks table only; no
   `dimensions`/`max_iter` arguments; no golden test (no R implementation).
   Rejected: biocore/mmvec (TensorFlow 1, unmaintained), a native port,
   dropping the task (design note 3).
6. **Transformers:** `ml.PrevalenceFilter(min_prevalence=0.1)` (raises when
   nothing passes; no `min_total`) and `ml.CLR(pseudocount=0.5)`;
   **`ml.RelativeAbundance` dropped** in favour of scikit-learn's
   `Normalizer(norm="l1")`, which the guide names (design note 5).
7. **Contracts for `ml`:** [pure-by-default](/decisions/pure-by-default.md)
   gains an `ml` estimators row (later `embed` on the `tl` convention,
   `to_torch` returning a new object over `X`);
   [function-shape](/contracts/function-shape.md) says the transformers are
   classes whose options `fit` validates and whose class docstring carries
   the skeleton (design note 6).
8. **Tooling for 4.3 (config only, no dependency):** an intersphinx entry for
   scikit-learn; `docs/_templates/autosummary/class.rst` listing only a
   class's own members; one `# type: ignore[no-untyped-call]` per
   `__sklearn_tags__` (mypy's `untyped_calls_exclude` does not reach
   `super()` calls).
9. **Reference plugin: MGM, shipped in biotapy** (`ml/_mgm.py`, registered by
   biotapy's own entry point). The only microbiome foundation model found
   with public weights, a licence (MIT) and a testable size (36 MB);
   BiomeGPT[^biomegpt] has no public weights, MGM2 needs a 650M-parameter nucleotide
   model and states no licence. Alternative: a separate `biotapy-mgm`
   package (design note 8).
10. **Extra `torch = ["torch>=2.9"]`** (2.9.0 is the first release with
    CPython 3.14 wheels).
11. **uv index for torch:** `https://download.pytorch.org/whl/cpu`,
    `explicit = true`, used only for torch, so CI and contributors install a
    192 MB CPU wheel (722 MB installed) instead of 555 MB plus CUDA 13
    packages. It changes `uv.lock`, not the published metadata.
    Alternative: a CI-only `uv pip install` step.
12. **CI:** a pytest marker `torch`, excluded by default like `r`, and a job
    `ml-extras` (Linux, Python 3.13, `-m torch`) in `check.needs`, with
    `tests/test_ci.py` tests. Alternative: torch tests local only.
13. **Extra `mgm = ["torch>=2.9", "transformers>=5"]`** and the entry point
    `[project.entry-points."biotapy.embeddings"] mgm =
    "biotapy.ml._mgm:embed"` in `pyproject.toml`. `microformer-mgm` is not a
    dependency (its pins conflict with numpy 2 and pandas 3).
14. **MGM weights at run time:** pooch downloads the 33 MB
    `microformer_mgm-0.5.8` wheel from files.pythonhosted.org, pinned by
    sha256, and extracts three files; nothing is bundled; the docs credit MGM
    and its MIT licence.
15. **Plugin interface:** entry-point group `biotapy.embeddings`; a plugin is
    a callable `embed(adata: AnnData, *, batch_size: int) -> np.ndarray`;
    biotapy refuses output that is not 2-D, finite, one row per sample;
    `ml.embed(adata, model, *, batch_size=64, inplace=False)` returns the
    array or writes `obsm["X_<model>"]` (design note 8) (amended by decision 27: no `batch_size`). A decision concept
    records it in 4C.
16. **Notebook content:** the leak-free CV page shows the prevalence filter's
    near-zero leak as measured and adds a supervised selection step to show
    a real one, on HMP2 (design note 9). The roadmap's literal version
    would show almost no difference (numbers: slice 4D decision 38).
17. **Exit gate 2 in CI, not in the docs build:** the `ml-extras` job's
    end-to-end MGM test is the proof; the docs page is a non-executed block
    quoting it (design note 10). Alternative: torch and transformers in the
    `doc` group (about 1 GB per docs build).
18. **Task changes:** new 4.A0 (`refactor(core)`), 4.B1 (CI), 4.4b (MGM),
    4.6b (embedding page), 4.D1 (Coming-from-R), 4.D2 (release); 4.2 moves
    to slice 4A; the mmvec tutorial on real paired metabolomics is deferred
    to 0.5 unless 4D finds a small public dataset (dataset found:
    slice 4D decision 37).
19. **Frontmatter:** the new `description`; `paths` gains
    `src/biotapy/_core/**`; sources gain the mmvec paper, MGM and the
    BiomeGPT preprint.
20. **Branch pushes:** approve push, PR and merge on green for the Phase 4
    slice branches (R13.3); no standing approval covers Phase 4.

**Slice 4B (approved by the user on 2026-10-09):**


21. **mypy does not follow torch** (`follow_imports = "skip"` and
   `ignore_missing_imports` for `torch`, `torch.*`), and the in-function
   `Dataset` subclass carries one `# type: ignore[misc]`. torch's names are
   `Any` to mypy with or without the extra, so the gate gives one answer and
   cold mypy skips 700 MB of torch. Alternative: follow torch when installed,
   which needs `unused-ignore` codes on lines whose error depends on the
   environment.
22. **Return annotation: the bare `"Dataset"`** (no `[Tensor]` parameter), plus
   an intersphinx entry for PyTorch's docs: the subscripted form renders as a
   broken reference that `nitpicky` fails in the torch-free docs build.
23. **The docstring example runs only in `ml-extras`**: a hook in the root
   `conftest.py` marks `biotapy.ml._torch`'s doctests `torch`. Alternative:
   `# doctest: +SKIP`, as the R bridges do, so the example never runs.
24. **A separate `build` commit (task 4.B0)** for the extra and the index,
   before `feat(ml)`, so the dependency change can be reviewed and reverted on
   its own. Alternative: fold it into 4.5, as the outline had it.
25. **No pooch cache in `ml-extras` until 4.4b**: no 4B test downloads
   anything (R2.3) (the cache is `ml-extras`' own, for the `mgm` step, decision 31). Alternative: add it now, as the outline said.
26. **Accept that every `uv run` CI job contacts download.pytorch.org while
   locking** (`uv.lock` is not committed; +0.27 s warm), and that CI always
   takes the newest torch >= 2.9. Alternative: keep the index out of
   `pyproject.toml` and install torch in `ml-extras` alone with `uv pip install
   --index-url` (decision 11's alternative), or commit `uv.lock` (reverses the
   "resolve fresh" choice in `.gitignore`).

**Slice 4C (approved by the user on 2026-10-09):**


27. **No `batch_size`: `ml.embed(adata, model, *, inplace=False)`, and a
    plugin is `embed(adata) -> np.ndarray`** (amends decision 15 and design
    note 8). MGM on the CPU is fastest one unpadded sample at a time and
    needs no memory beyond the model's: 520 GlobalPatterns profiles took
    8.8 s at batch 1, 17.0 s at 8 (+196 MB), 18.5 s at 64 (+1.6 GB), on 8
    threads (11.9 s, 25.1 s, 26.6 s on 4). With no plugin that benefits,
    the parameter would be speculative (R2.3); adding it later is a
    keyword-only addition. Alternative: keep `batch_size=64` (MGM pads
    batches: twice as slow and +1.6 GB), or keep it with default 1.
28. **MGM's sample embedding is the mean of the last hidden layer over the
    sample's tokens, `<bos>` and `<eos>` included**: MGM's paper uses
    "element-wise mean pooling" for the pretrained model (Methods 4.5), and
    the mean stays defined for a sample with no known genus. Alternative:
    the `<eos>` token's state, as MGM's notebook does for a fine-tuned
    model (one line in `_mgm.py` and in `tests/mgm/export_reference.py`,
    then regenerate `embeddings.csv`), or the mean over genus tokens only.
29. **Features of one MGM genus are summed, as MGM does; no "one feature per
    genus" error** (the outline proposed one, naming `pp.tax_glom`).
    GlobalPatterns after `tax_glom(..., "genus")` has 996 features for 983
    genus names (`tax_glom` groups by lineage), so that error would refuse
    exit gate 2's own input; summing makes an OTU table and its `tax_glom`
    give the same embedding (tested). Alternative: raise when a genus
    repeats.
30. **MGM's genus rule as is**: `g__[A-Za-z0-9_]+` on `"g__" + genus`
    (`Escherichia-Shigella` reads as `Escherichia`, `[Ruminococcus]` as
    nothing), with no rewriting of spaces or brackets: no genus name in
    GlobalPatterns or enterotype has a space, so a rewrite would change
    nothing measured (R2.3). Alternative: spaces to underscores before the
    regex (MGM's vocabulary spells `Candidatus_Accumulibacter`).
31. **A marker `mgm` for the tests that need the extra and the download,
    run by `ml-extras` with a pooch cache of its own** (amends the outline's
    `torch` + `network`, decision 25's cache and exit gate 2's wording). The
    `network` job runs `-m "network or golden"` without torch, so a
    `network` mark would fail there; the R bridges set the precedent (marker
    `r` only, the network job's cache). Exit gate 2 becomes
    `tests/ml/test_mgm.py::test_embeds_global_patterns_end_to_end` in
    `ml-extras`. Alternative: `torch` + `network` marks and the network job
    changed to `-m "(network or golden) and not torch"`.
32. **A new task 4.C0, `refactor(core)`: the pooch cache moves to
    `_core.make_pooch`** (R4.3: `ml` is the second subpackage that downloads;
    the cache directory and `BIOTAPY_DATA_DIR` stay defined once).
    Alternative: a second `pooch.create(path=pooch.os_cache("biotapy"),
    ..., env="BIOTAPY_DATA_DIR")` in `ml/_mgm.py`.
33. **Parity fixture and tolerance**: `tests/data/mgm/` (a synthetic table
    and MGM 0.5.8's own embeddings of it, MIT, `NOTICE.txt`), written by
    `tests/mgm/export_reference.py` in an isolated Python 3.11 environment,
    compared at `atol=1e-3`. biotapy and MGM agree to 1.7e-6, but torch 2.13
    and 2.14 sometimes compute their first `tanh` in a process less
    precisely on a CPU running more than four threads, which moves that
    call's embedding by up to 1.6e-4; biotapy does not work around torch.
    Alternative: call `torch.tanh(torch.zeros(1))` once before the forward
    (removed it in 16 of 16 runs on the prototype) and compare at 1e-5.
    (Amended 2026-10-10 at the user's request: the race is oneMKL's VML
    CPU-type cache, pytorch/pytorch#188792; biotapy fills the cache with one
    serial call and tests every call at 1e-5.)
34. **`ml.embed`'s docstring example runs MGM** (`bt.ml.embed(bt.pp.tax_glom(
    bt.datasets.toy(), "genus"), "mgm").shape` -> `(6, 256)`), marked `mgm`
    from 4.4 on; between 4.4 and 4.4b, `-m mgm` fails on it, and no CI job
    runs `-m mgm` until 4.4b adds the step. Alternative: `# doctest: +SKIP`
    in 4.4, removed in 4.4b; or one commit for 4.4 and 4.4b.
35. **Report the torch race upstream after 0.4**, with your approval of the
    text (a GitHub action outside this repository): a ten-line reproducer
    (`torch.tanh` on a seeded 7 x 512 x 1024 tensor, first call against the
    second, 8+ threads) shows it on torch 2.13.0+cpu and 2.14.1+cpu.
    Alternative: do not report. (Superseded: the issue exists upstream as
    pytorch/pytorch#188792.)

**Slice 4D (approved by the user on 2026-10-10):**

36. **Task 4.F5: fix `_core.sum_pairs` too**, the same one-line defect as `sum_by` (an `int8`
    table's `func_glom` gives -56 for 200), in its own `fix(core)` commit after 4.F4, so it can be
    dropped alone. Alternative: report it and fix it later (R1.4: the approval named `sum_by`).
37. **Task 4.6c: add `bt.datasets.biocrust()` and an executed multi-omics tutorial now**, as
    decision 18 allowed once a small public paired dataset was found: mmvec's own example
    (`biocore/mmvec` `examples/soils` at `88ca33b`, 135 KB, BSD-3-Clause repository; 19 shared
    samples, 466 microbes, 85 metabolites), downloaded through `_remote.py`'s pooch with pinned
    URLs and SHA-256, never committed. `bt.tl.mmvec` reproduces the check mmvec's example makes (13
    of 13 listed metabolites rank above zero for `rplo 1`, seeds 0-4), its first test on real
    data. It adds one public function (`datasets.biocrust() -> MuData`, `R equivalent: none`) and a
    run-time download. Alternative: defer the tutorial and the loader to 0.5 with the dataset now
    named, and keep the guide's toy example for 0.4.
38. **The leak-free notebook's protocol**: `RepeatedStratifiedKFold(n_splits=5, n_repeats=5,
    random_state=0)`, `CLR(pseudocount=1e-6)`, `LogisticRegression(max_iter=5000)`,
    `SelectKBest(f_classif, k=20)`, and the shuffled labels averaged over 10 permutations of
    `np.random.default_rng(0)`; quoted numbers 0.527 / 0.526 (filter inside / outside), 0.558 /
    0.705 (selection), 0.544 / 0.767 (selection on shuffled labels). Alternative: design note 9's
    single 5-fold split and one shuffle (0.478 / 0.475, 0.571 / 0.722), whose honest score sits
    below 0.5 and moves by about 0.05 with the split.
39. **Tutorial numbers are kept honest by running the page**: a `run_page` fixture in
    `tests/conftest.py` executes a page's code cells and Python blocks; a `network` test (4.6,
    4.6c) or an `mgm` test (4.6b) compares what the page computes with constants, and a
    default-run test checks the prose quotes the same constants. Alternative: Phase 3's pattern
    only (constants shared by a test that repeats the page's code and a prose test).
40. **Exit gate 2's page is a separate tutorial, `docs/tutorials/embeddings.md`, in plain
    Markdown that the docs build does not run**, quoting the end-to-end test's outputs and one
    laptop timing (not pinned). Alternative: a subsection of the ML guide's "Embeddings".
41. **Release text** (4.D2 Steps 4-5): Added for the seven new public names, MGM, the two extras
    and the docs; Changed for the `(pp.clr)`/`(pp.philr)` warning suffix, `filter_features`'
    `TypeError`, the 64-bit sums (an `int32` table now sums to `int64`) and "not in 0.4"; README
    "0.4", a Machine learning line, the `torch` and `mgm` install lines with the CPU-index tip,
    and "Next: what users and issues ask for" (Phase 5 is an unplanned backlog). Alternative:
    list Phase 5's candidates as "Next".
42. **Approvals:** at Checkpoint D, one message for the slice review, `phase-4d`'s and
    `release-0.4.0`'s push/PR/merge, and the PyPI publisher check; after the release PR merges,
    one message for the tag, the tag push, the GitHub release (the PyPI upload) and
    `close-phase-4`'s push/PR/merge (4.D2 Step 10). Nothing is pushed before.
43. **MGM keeps GlobalPatterns' one-word *Candidatus* genera out, as MGM's own code does**
    (decision 30's rule): 28 of the 96 features it leaves out are Greengenes names such as
    `CandidatusPelagibacter`, and MGM's vocabulary holds 23 of them as `Candidatus_<Name>`. The
    page and `ml.md` say so. Alternative: rewrite `Candidatus<Name>` to `Candidatus_<Name>` before
    MGM's regex, which would leave out 73 instead of 96 (counted from the vocabulary, not run), move
    exit gate 2's numbers, and differ from MGM's code on the same input; a 0.5 candidate.
44. **`verified` on [embedding-plugins](/decisions/embedding-plugins.md) is yours to add**
    (R12.3); 4.7 only changes it from `draft` to `stable`, on the strength of your approval of
    slice 4C. Say if you have not confirmed it, and it stays `draft`.
45. **After the release, Phase 5 becomes `in-progress`** so that rules.md R1.2's "exactly one
    phase in progress" holds; its own rule still gives no item a task until an issue or you ask.
    Alternative: no phase in progress, which needs R1.2 reworded ("at most one"); rules.md is
    yours to change.

# Self-review
Run against the brief, the roadmap outline and the writing-plans checklist.

1. **Spec coverage.**

   | Roadmap item | Where it lands |
   |---|---|
   | Design note "classes only for protocols" | Design notes 5, 7; Global constraints |
   | Design note "leakage from stateful steps" | Design note 5 (measured), 9 (notebook); 4.3 tests |
   | Design note "embedding plugins" | Design note 8; 4C |
   | Design note "saving loses trees" | Design note 2; 4.1 test and docs; decision 4 |
   | 4.1 MuData conventions | Task 4.1 (signature changed, decision 2) |
   | 4.2 `tl.mmvec` | Task 4.2 (scikit-bio route, decision 5) |
   | 4.3 transformers, `check_estimator` | Tasks 4.A0, 4.3 (`parametrize_with_checks`; `RelativeAbundance` dropped, decision 6) |
   | 4.4 `ml.embed`, reference plugin | 4C: 4.4, 4.4b |
   | 4.5 `ml.to_torch`, extra `torch` | 4B: 4.5, 4.B1 |
   | 4.6 leak-free CV notebook | 4D: 4.6 |
   | 4.7 knowledge | Checkpoint A (`ml` Module concept), 4C (plugin decision), 4D (4.7) |
   | Exit gate 1 | 4.6, executed in the docs job |
   | Exit gate 2 | 4.4b's end-to-end test in `ml-extras`; 4.6b page |
   | Spec: multi-omics tutorial | 4.6c (decision 37) |

   Brief research questions: mmvec (design note 3, decision 5); scikit-learn
   protocol and leakage (design note 5); plugins and weights (design note 8,
   decisions 9, 13-15); torch Dataset and collate (design note 7); Python
   3.14 (design note 11). Gaps found and fixed while writing: the roadmap's
   `**modalities` broke R3.1/R3.7 (decision 2); `check_estimator` cannot run
   under the gate's `-W error::UserWarning` (design note 5); the roadmap's
   notebook would show no leak (design note 9); `tl.mmvec`'s purity needed
   a pure-by-default line (design note 3).
2. **Placeholder scan.** Every slice 4A step carries the full file or the
   exact diff, rendered from the scratch clone's commits that passed every
   gate (`1ef640e`, `0c1c5f0`, `64fb759`, `c9536aa` on branch `phase-4a`).
   Slice 4B carries full steps, rendered from the scratch clone's commits
   `fada00b`, `875110f`, `5350e0e` on branch `phase-4b` (rebased onto
   `ef82843`).
   Slice 4C carries full steps, rendered from the scratch clone's commits
   `82c6e6d`, `5b73a1b`, `0785751`, `a702d59` on branch `phase-4c`.
   Slice 4D carries full steps, rendered from the scratch clone's commits `b5654af`..`ceada96`
   on branches `phase-4d` and `release-0.4.0`.
   The only open values are the log section's `<date>` and each concept's
   `generated.at`, filled at commit time.
3. **Type consistency.** Used throughout: `to_mudata(modalities:
   Mapping[str, AnnData]) -> MuData`; `mmvec(mdata, *, microbes="taxa",
   metabolites="metabolites", seed=None) -> pd.DataFrame`;
   `check_pseudocount(pseudocount: object) -> None`; `pseudocounted(X,
   pseudocount, *, func, columns=None)`; `PrevalenceFilter(*,
   min_prevalence: float = 0.1)` with `prevalence_`; `CLR(*, pseudocount: float = 0.5)`;
   `to_torch(adata, *, label_key=None, layer=None)`; `embed(adata, model, *,
   inplace=False)`; plugin `embed(adata)`; `make_pooch(base_url, registry, *,
   urls=None)`; `obsm["X_<model>"]`.
4. **Review focus.** Items 1-4 name tests that exist in the rendered 4A
   code (checked by searching this file for each name); item 5 names tests
   4C must write under those names.
5. **Known residual risks.**
   - The gate counts were measured on Python 3.13 only; the CI matrix adds
     3.12, 3.14 and pre-release dependencies. scikit-learn's estimator
     checks are the likeliest to differ under pre-releases.
   - That the pooling follows the paper's prose (decision 28). (The torch
     first-`tanh` race, decision 35, is fixed by a serial call; see decision 33.)
   - Upstream: the weights live in a PyPI wheel the MGM authors control;
     MGM2 may supersede MGM before 0.4 ships.

**Slice 4D:**

1. **Spec coverage.**

   | Outline or brief item | Where it lands |
   |---|---|
   | 4.F4 `sum_by` overflow (user-approved fix) | Task 4.F4; design 1; 4.F5 for the same line in `sum_pairs` (decision 36) |
   | 4.6 leak-free CV notebook, executed in docs (exit gate 1) | Task 4.6; design 2; decisions 38-39 |
   | Near-zero prevalence leak and a large supervised leak (decision 16) | Task 4.6's sections and prose test |
   | A test pins the quoted numbers (Phase 3 precedent) | `LEAK_FREE_CV_AUC`, `GLOBAL_PATTERNS_*`, `MICROCOLEUS_METABOLITES` with `run_page` (decision 39) |
   | 4.6b end-to-end embedding page, static outputs, kept honest | Task 4.6b; design 3; decision 40 |
   | Multi-omics tutorial, deferred unless a dataset is found (decision 18) | Task 4.6c; design 4; decision 37 |
   | 4.D1 Coming-from-R (`to_mudata` rows, a check) | Task 4.D1; design 5 |
   | 4.7 knowledge, `verified` only by the user | Checkpoint D; Task 4.7; decision 44 |
   | 4.D2 version, CHANGELOG (warning suffix, extras, `ml`), playbook, METADATA and entry point, approvals, close Phase 4 | Task 4.D2 Steps 1-14; decisions 41, 42, 45 |
   | 4.F1 upstream issue scikit-bio#2631 | "Other places" (4.F1 box, `da.md`) |
   | Benchmarks | none (outline: no measurement asks for one, R10.1) |

2. **Placeholder scan.** Every test, page, loader and diff block is rendered by `git show`/`git
   diff` from the scratch commits that passed the gates, not retyped. Values left for run time,
   each named where it is used: the implementer's model and UTC time in `generated`, the parent's
   short hash in `commit`, the dates in log headings and the CHANGELOG, PR numbers, Test run,
   docs job and release workflow ids, and `<scratchpad>`.
3. **Type consistency.** `run_page(path: Path) -> dict[str, object]` is defined in 4.6 and used
   by 4.6b and 4.6c; the names each test reads (`inside` ... `shuffled_outside`; `embedding`,
   `same_type`; `microcoleus`, `above_zero`) are defined by the pages shown; `biocrust() ->
   MuData` with `"taxa"` and `"metabolites"` matches `tl.mmvec`'s defaults; `sum_by`/`sum_pairs`
   keep their signatures.
4. **Review focus.** Each item names tests that exist in the blocks above (checked by searching
   the rendered draft) or the 4.D2 build check.
5. **Known residual risks and what was not verified.**
   - CI did not run the slice: the `network` job running the leak-free page (85 s of CPU on the
     prototype's 16 threads; GitHub's runner has 4 vCPUs), the `docs` job, and Read the Docs
     executing two more notebooks happen first on the PR.
   - The gates ran on CPython 3.13 on Linux only; the tutorials execute only in Linux jobs, and
     their AUCs were not measured on another BLAS.
   - mmvec's example data: the repository's licence and README were read; where its authors took
     the tables from was not traced, so the pages cite only the mmvec paper and the repository.
   - Decision 43's "73 instead of 96" is counted from MGM's vocabulary, not run through `embed`.
   - 4.D2 Steps 9-14 (push, CI, tag, upload, PyPI checks, close) cannot run before approval;
     Step 12's expected output was checked against the clone's installed 0.4.0, not PyPI.

[^spec]: Python Microbiome Toolkit development report, sections Positioning and Roadmap
[^mmvec]: Morton et al. 2019, Learning representations of microbe-metabolite interactions, Nature Methods
[^mgm]: MGM foundation model (microformer-mgm), HUST-NingKang-Lab
[^biomegpt]: BiomeGPT preprint
