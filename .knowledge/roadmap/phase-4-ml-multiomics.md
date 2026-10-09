---
type: Phase
title: Phase 4 - ML-ready and multi-omics (0.4)
description: io.to_mudata over shared samples, tl.mmvec through scikit-bio, leak-free scikit-learn transformers (PrevalenceFilter, CLR), a PyTorch Dataset behind the extra torch, an entry-point interface for embedding models with MGM as the reference plugin, a leak-free cross-validation notebook.
tags: [roadmap, ml, multiomics]
status: stable
release: "0.4"
phase_state: in-progress
effort: ~4-6 weeks part-time
depends_on: [/roadmap/phase-3-stats.md]
paths: ["src/biotapy/ml/**", "src/biotapy/tl/**", "src/biotapy/io/**", "src/biotapy/_core/**"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-09T03:37:48Z }
commit: cdc3b07
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
---

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Slice 4A has full TDD steps;
> slices 4B-4D are outlines, expanded (superpowers:writing-plans) and approved
> when reached (rules.md R1.2a).

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
(extra `mgm`, slice 4C) · pooch · pytest/hypothesis.

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
  cannot be installed beside biotapy. The weights load into
  `transformers.GPT2Model` 5.19.0 with no missing or unexpected keys
  (`torch.load(..., weights_only=True)`, keys stripped of `transformer.`),
  and embed 4 samples in 0.3 s on CPU.
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
  (shuffled labels 0.754 vs 0.550). Design note 9 uses this.

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
   4B, outline). (user: extra, CI)**
   - Map-style dataset whose `__getitem__(i)` densifies row `i` of the CSR
     table to a float32 tensor: the full table is never dense (R6.2) and the
     default `collate_fn` stacks the rows [V]. Per-row sparse tensors would
     need a custom `collate_fn` and trip `-W error::UserWarning` [V], so they
     are not used.
   - `label_key` names an `obs` column: categories become int64 codes in
     category order, numbers float32; a missing value raises. Without it,
     items are bare tensors.
   - torch must not be imported at module level, and the
     `import-without-extras` CI job imports every module, so the `Dataset`
     subclass is defined inside the function after
     `import_optional("torch", extra="torch")` [UNVERIFIED: mypy treats the
     base as `Any` when torch is absent from the type-check environment and
     reports `disallow_subclassing_any`; 4B measures this before choosing a
     targeted ignore or a plain map-style class, which `DataLoader` also
     accepts].
   - The extra is `torch = ["torch>=2.9"]`: 2.9.0 is the first release with
     CPython 3.14 wheels [V]. CI installs CPU wheels through a uv index
     (`[[tool.uv.index]] pytorch-cpu`, `explicit = true`, with
     `[tool.uv.sources] torch`), so the job does not pull CUDA
     [UNVERIFIED: uv's lock across platforms with that source; macOS wheels
     on the CPU index]. A new marker `torch` is excluded from the default
     run like `r`, and a CI job `ml-extras` (Linux, Python 3.13) runs
     `-m torch` and joins `check.needs`.

8. **Embedding plugins: `ml.embed(adata, model, *, batch_size=64,
   inplace=False)`; MGM is the reference plugin (slice 4C, outline).
   (user: interface, reference plugin, extra)**
   - Discovery: `importlib.metadata.entry_points(group="biotapy.embeddings")`
     at call time; an unknown `model` raises `KeyError` listing the installed
     names. A plugin entry point loads a callable
     `embed(adata: AnnData, *, batch_size: int) -> np.ndarray`.
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
     `g__<genus>`, drops genera outside the 9,665-token vocabulary with one
     warning naming how many, and repeats MGM's preprocessing: relative
     abundance, z-score per genus with `phylogeny.csv`'s mean and standard
     deviation, keep genera above the z-score of zero, sort descending, wrap
     in `<bos>` ... `<eos>`, pad or cut to 512 [V, read from
     `mgm/src/MicroCorpus.py`]. The sample embedding is the mean of the last
     hidden layer over non-padding tokens [UNVERIFIED: MGM's example notebook
     may use another pooling; 4C reads it first].
   - Parity: 4C tries to generate a golden fixture by running
     `microformer-mgm` 0.5.8 in a Python 3.11 environment with its own pins
     on a small genus table [UNVERIFIED that it installs].
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
   0.550 on shuffled labels [V], where the leak manufactures signal from
   noise. It is executed in the docs build (exit gate 1).

10. **Exit gate 2, "one foundation model plugged in end to end", is proven
    in CI, not in the docs build. (user)** The docs build has no torch; adding
    torch and transformers to the `doc` group would add about a gigabyte to
    every docs build [V: 946 MB venv]. A `torch`+`network` test in the
    `ml-extras` job downloads the MGM weights, embeds GlobalPatterns' genera
    and checks shape, determinism and (if 4C's fixture exists) parity; the
    docs page shows the same code as a non-executed block and quotes the
    test's numbers, the Phase 3 pattern for the R bridges.

11. **Python 3.14.** torch >= 2.9 has cp314 wheels for Linux, macOS and
    Windows [V]; transformers is pure Python; scikit-learn, mudata and
    scikit-bio already run there in the existing matrix. The `ml-extras`
    job pins 3.13, like `r-bridge` [UNVERIFIED: transformers' tokenizers
    wheel on 3.14, which the MGM plugin does not call].

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
| 4.5 | extra `torch` | `torch>=2.9` (BSD-3-Clause and bundled licences) | `ml.to_torch`; the roadmap's planned extra |
| 4.5 | uv config | index `https://download.pytorch.org/whl/cpu` for torch (`explicit = true`) | CI and contributors get a 192 MB CPU wheel, not 555 MB plus CUDA 13 packages |
| 4.5 | CI | job `ml-extras` (Linux, Python 3.13, `-m torch`), in `check.needs` | run what the default jobs skip |
| 4.4 | extra `mgm` | `torch>=2.9`, `transformers>=5` (Apache-2.0) | the MGM reference plugin |
| 4.4 | data, run time | MGM 0.5.8 weights from files.pythonhosted.org via pooch (MIT, 33 MB download) | never bundled (R6.6) |
| - | none | `microformer-mgm`, `huggingface_hub`, TensorFlow, biocore `mmvec` | not imported (design notes 3, 8) |

# Review focus
The five ways real users are most likely to get a wrong answer from Phase 4
without an error. Each line names the test that pins it.

1. **Preprocessing fitted on every sample before cross-validation.**
   Expected: the filter learns from the training rows only and is refitted
   per fold. Tests: 4.3 `test_prevalence_filter_learns_only_from_the_samples_it_is_fitted_on`,
   `test_a_pipeline_refits_the_filter_in_every_fold`; 4D's notebook shows
   the size of the leak a supervised step causes.
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
   the plugin. Tests (4C, names to be written there):
   `test_plugin_with_wrong_row_count_raises`,
   `test_plugin_with_nan_raises`.

# Slices
| Slice | Delivers | Tasks | Ends with |
|---|---|---|---|
| **4A - No new dependency** | multi-omics MuData, mmvec, leak-free transformers | 4.1 `io.to_mudata` · 4.2 `tl.mmvec` · 4.A0 `refactor(core)` pseudocount step · 4.3 `ml.PrevalenceFilter`, `ml.CLR` | Checkpoint A |
| **4B - torch** | the extra `torch`, `ml.to_torch`, CPU wheels in CI | 4.5 `ml.to_torch` (+ extra, uv index, marker) · 4.B1 CI job `ml-extras` | Checkpoint B |
| **4C - Embeddings** | the plugin interface and MGM | 4.4 `ml.embed` and the entry-point group · 4.4b MGM plugin (+ extra `mgm`, weights via pooch, parity fixture) | Checkpoint C |
| **4D - Docs and release** | the guide pages, the two exit-gate pages, 0.4 | 4.6 leak-free CV notebook · 4.6b embedding page · 4.D1 Coming-from-R check · 4.7 knowledge · 4.D2 release 0.4.0 | exit gate |

Execution order inside 4A: **4.1 -> 4.2 -> 4.A0 -> 4.3 -> Checkpoint A.**
`tl.mmvec`'s tests build their MuData with `io.to_mudata`; the `refactor`
lands before the `feat(ml)` that needs it.

# Tasks (checklist)
- [x] 4.1 `io.to_mudata(modalities) -> MuData` and the multi-omics decision
- [x] 4.2 `tl.mmvec(mdata, *, microbes="taxa", metabolites="metabolites", seed=None) -> pd.DataFrame`
- [x] 4.A0 `refactor(core)`: the pseudocount step moves to `_core/_composition.py`
- [x] 4.3 `ml.PrevalenceFilter(min_prevalence=0.1)`, `ml.CLR(pseudocount=0.5)`; scikit-learn's estimator checks pass
- [x] Checkpoint A
- [ ] 4.F1 `da.ancombc2`: pin scikit-bio's bias E-M underflow; the schema property keeps to fittable designs
- [ ] 4.F2 `pp.filter_features` names a wrongly typed threshold
- [ ] 4.F3 `refactor(core)`: one finite, non-negative check
- [ ] 4.5 `ml.to_torch(adata, *, label_key=None, layer=None) -> torch.utils.data.Dataset` and the extra `torch`
- [ ] 4.B1 CI job `ml-extras` for `-m torch` tests
- [ ] Checkpoint B
- [ ] 4.4 `ml.embed(adata, model, *, batch_size=64, inplace=False)` and the entry-point group `biotapy.embeddings`
- [ ] 4.4b MGM reference plugin and the extra `mgm`
- [ ] Checkpoint C
- [ ] 4.6 Leak-free cross-validation notebook (exit gate 1)
- [ ] 4.6b End-to-end embedding page (exit gate 2's docs)
- [ ] 4.D1 Coming-from-R check
- [ ] Checkpoint D
- [ ] 4.7 Knowledge: `ml` Module concept, the plugin-interface decision
- [ ] 4.D2 Release 0.4.0

# Exit gate
- [ ] Leak-free CV example executed in docs: `docs/tutorials/leak_free_cv.md`
  runs in the docs CI job (design note 9).
- [ ] One foundation model plugged in end to end: the `ml-extras` job's
  `torch`+`network` test embeds GlobalPatterns' genera with MGM through
  `bt.ml.embed(..., "mgm")` (design note 10).
- [ ] All Phase 1-3 gates still green.

# Risks
- **MGM's pooling and parity are unverified** -> 4C reads MGM's example
  notebook before writing `_mgm.py`, and tries a Python 3.11 golden fixture;
  without one, the docstring says the embedding is biotapy's reading of
  MGM's preprocessing, checked only for shape and determinism.
- **The weights URL is a PyPI file** -> files.pythonhosted.org URLs are
  immutable and pinned by hash; if the release is yanked the download still
  works (yanking does not delete files), but a deletion would break the
  plugin, so the registry entry names the version.
- **torch install size and platform wheels** -> CPU index in CI; users install
  what suits their hardware; the docs say so.
- **mypy without torch** -> design note 7's unverified point; measured
  first in 4B.
- **mmvec convergence is silent**: scikit-bio's L-BFGS prints its status
  only with `verbose=True` and `MMvecResult` has no flag [V] -> `Notes`
  points to `res.convergence` through scikit-bio; reconsider a warning if a
  user reports a non-converged fit.
- **scikit-learn's estimator checks change between releases** -> the check
  suite runs in every test job; a new check that fails is a real
  incompatibility, fixed in biotapy rather than skipped (R11.5).
- **Leakage misread as "the prevalence filter leaks a lot"** -> the
  notebook prints the measured near-zero difference first (design note 9).

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
Unfixed on scikit-bio `main`; no upstream issue. Rate: 0.16% of 6 x 7 tables
with a 2-sample reference group, 0 with 3 or more; the property test fails in
about 1.2% of runs, then every run once Hypothesis has saved the example.

biotapy already reports it clearly (`ValueError: da.ancombc2: scikit-bio could
not fit the model: ...`, R7.4). A biotapy-side fix would mean monkeypatching
scikit-bio's private function or reimplementing the E-M, so the fix is: pin the
behaviour, restrict the property test to fittable designs for this one defect,
document it, and report it upstream.

**Files:** modify `tests/da/test_ancombc.py`, `tests/da/test_schema.py`,
`src/biotapy/da/_ancombc.py` (Notes only).

- [ ] **Step 1: regression test (pins current behaviour).** Append to
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
- [ ] **Step 2: RED.** In `tests/da/test_schema.py`, import `example` and
  `reject` from hypothesis and add, under the property test's `@given`:
  ```python
  @example(counts=np.array(EM_UNDERFLOW), split=2)
  ```
  with `EM_UNDERFLOW` copied as a module constant (tests do not import each
  other). Run `uv run --group test pytest -q tests/da/test_schema.py -k direction -W error::UserWarning`.
  Expected: `1 failed, 2 passed` (ancombc2 fails on the example).
- [ ] **Step 3: GREEN.** Replace `out = method(adata, "g")` in that test with:
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
- [ ] **Step 4: Notes.** In `da.ancombc2`'s docstring Notes, after the E-M
  sentence, add: "With two samples in the reference level, scikit-bio 0.7.4 can
  fail to estimate the bias (about 1 in 600 random small tables) where R
  returns results; biotapy raises naming scikit-bio rather than guess." Docs
  build with `-W`.
- [ ] **Step 5: commit** `test(da): pin scikit-bio's ancombc2 bias E-M underflow and keep the schema property on fittable designs`.
- [ ] **Upstream (needs separate approval, R13.3):** an issue on
  scikit-bio/scikit-bio; its draft text is in the 4.F PR description.

### Task 4.F2: `pp.filter_features` names a wrongly typed threshold

**Files:** `src/biotapy/pp/_filter.py`, `tests/pp/test_filter.py`.
Same rule as `ml.PrevalenceFilter` (R3.5): a bool or a non-real
`min_prevalence` or `min_total` raises `TypeError` naming the argument.

- [ ] **Step 1: RED.**
  ```python
  @pytest.mark.parametrize("value", [True, "0.5", [0.5]])
  @pytest.mark.parametrize("argument", ["min_prevalence", "min_total"])
  def test_wrongly_typed_threshold_raises_naming_it(argument, value):
      with pytest.raises(TypeError, match=f"{argument} must be a real number"):
          bt.pp.filter_features(bt.datasets.toy(), **{argument: value})
  ```
  Expected: 6 failed (`True` passes silently; the others raise an unnamed
  `TypeError` or numpy error).
- [ ] **Step 2: GREEN.** Before the range checks, for each non-None
  threshold: `if isinstance(value, bool) or not isinstance(value, int | float | np.integer | np.floating): raise TypeError(f"{name} must be a real number, got {value!r}")`,
  written once as a loop over the two names (R5 limits). Same wording as
  `ml.PrevalenceFilter`. Expected: tests/pp passes.
- [ ] **Step 3: commit** `fix(pp): name min_prevalence or min_total when its type is wrong`.

### Task 4.F3: one finite, non-negative check in `_core`

The predicate `np.isfinite(X.data).all() and not (X.data < 0).any()` is
written in `_core/_composition.py`, `tl/_mmvec.py` and `fn/_redundancy.py`
(R4.3). Move it to `_core/_matrix.py` as a predicate, so each caller keeps its
own message (no user-visible change).

**Files:** `src/biotapy/_core/_matrix.py`, `src/biotapy/_core/__init__.py`,
the three callers, `tests/core/test_matrix.py` (or the existing `_core`
matrix test file).

- [ ] **Step 1: RED.**
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
- [ ] **Step 2: GREEN.**
  ```python
  def finite_non_negative(X: sp.csr_matrix) -> bool:
      """True when every stored value of ``X`` is finite and >= 0."""
      return bool(np.isfinite(X.data).all() and not (X.data < 0).any())
  ```
  Export it from `_core/__init__.py`; replace the three inline predicates.
  `tests/pp`, `tests/tl`, `tests/fn`, `tests/core`, `tests/ml` pass unedited.
- [ ] **Step 3:** update `modules/core.md` (entry point) and log line.
  Commit `refactor(core): share the finite, non-negative check`.

### Gates and merge
Commit, `git status --short` empty, then: prek; `pytest -q -W error::UserWarning`
(1412 + 1 + 6 + 5 passed, 2 skipped, 54 deselected); `-m "golden or network"` 36;
docs `-W`; `knowledge_stale.sh` 0 stale. Push `fix-4f`, PR, merge commit on
green (decision 20 covers Phase 4 branches).

---
## Slice 4B - torch (outline; expand with superpowers:writing-plans when reached)

**Goal:** `pip install 'biotapy[torch]'` gives `bt.ml.to_torch`, a PyTorch
dataset over a table whose rows stay sparse until they are fetched, and CI
runs the torch tests on CPU wheels.

**Tasks.**
- **4.5 `ml.to_torch(adata, *, label_key=None, layer=None) -> Dataset`**
  (`ml/_torch.py`, `tests/ml/test_torch.py` marked `torch`).
  - Start of task: ask for the extra `torch = ["torch>=2.9"]`, the uv index
    and the marker (decisions 10-12). Then `uv lock` and record the lock diff
    size.
  - Design as in design note 7: rows densified one at a time to float32;
    `label_key` -> int64 category codes (categorical or string column) or
    float32 (numeric); missing labels raise naming the column; `layer` reads
    `adata.layers[layer]` (CSR or dense; `layers["clr"]` is dense and is
    taken as is); `X`/the layer is referenced, not copied, and `Notes` says
    that editing `adata` later changes what the dataset returns.
  - First measurement, before the tests are written: run `mypy --strict`
    with the `Dataset` subclass defined inside the function and torch absent
    from the type-check environment. If it reports
    `disallow_subclassing_any`, choose between one
    `# type: ignore[misc]` with a comment and a plain map-style class (no
    torch base), and record the choice under "Slice 4B decisions".
  - Tests: batch shapes and dtypes through `DataLoader(batch_size=4)`; rows
    equal `X.toarray()` rows; labels in category order; NaN label raises;
    `layer` read; purity (`adata` unchanged; the full table never densified,
    checked by a CSR `X` of 2,000 x 50,000 at 0.1% density whose dense form
    would be 800 MB); without torch, `ImportError` naming
    `biotapy[torch]` (monkeypatched `import_optional`).
  - Docs: the ML guide gains "PyTorch"; API entry.
- **4.B1 CI job `ml-extras`**: Linux, Python 3.13, `uv run --group test
  --extra torch pytest -m torch`, pooch cache like `network`, added to
  `check.needs`; `tests/test_ci.py` gains the job's tests (runs `-m
  torch`, blocks merges, installs from the CPU index). `addopts` excludes
  `torch` like `r`. The job's runtime and the installed size go into the
  slice record.

**Open questions, each with a proposed answer.**
- uv index scope: per-package source for torch only, `explicit = true`, so
  no other package resolves from the PyTorch index. [UNVERIFIED: `uv lock`
  across Linux, macOS and Windows markers with that source.]
- Does `torch` belong in the `import-without-extras` job's install? No: that
  job exists to prove biotapy imports without it.

---
## Slice 4C - Embeddings (outline)

**Goal:** `bt.ml.embed(tdata, "mgm")` returns one MGM embedding per sample
(or writes `obsm["X_mgm"]`), and any package can add a model by declaring an
entry point in `biotapy.embeddings`.

**Tasks.**
- **4.4 `ml.embed(adata, model, *, batch_size=64, inplace=False) ->
  np.ndarray | None`** (`ml/_embed.py`, `tests/ml/test_embed.py`; no torch
  needed for these tests).
  - Discovery through `importlib.metadata.entry_points(group=
    "biotapy.embeddings")`; tests register a fake plugin by monkeypatching
    the `entry_points` that `_embed.py` imports (no installed test
    distribution).
  - Tests: an installed plugin's array is returned and, with
    `inplace=True`, written to `obsm["X_<model>"]` with `None` returned;
    unknown model -> `KeyError` listing installed names;
    `test_plugin_with_wrong_row_count_raises`,
    `test_plugin_with_nan_raises`, a 1-D result raises; `batch_size < 1`
    raises; purity with `inplace=False`; a feature change drops
    `obsm["X_<model>"]` (data-model-slots Propagation, already true).
  - Knowledge: the decision concept "embedding plugins" (group name, the
    callable's signature, the output checks, weights never bundled);
    data-model-slots' `X_<plugin>` key becomes `X_<model>` with the
    plugin rule; pure-by-default's `ml` row gains `embed`.
- **4.4b MGM plugin** (`ml/_mgm.py`, extra `mgm`, `tests/ml/test_mgm.py`
  marked `torch` and, for the download, `network`).
  - Start of task: ask for the extra `mgm = ["torch>=2.9",
    "transformers>=5"]` and the entry point in `pyproject.toml`
    (decisions 13-14).
  - Read MGM's example notebook for its sample-embedding pooling; use that
    pooling, or mean pooling over non-padding tokens if it names none, and
    say which in `Notes`.
  - pooch registry entry: the 0.5.8 wheel URL and sha256 (design note 8),
    `pooch.Unzip(members=["mgm/resources/general_model/config.json",
    "mgm/resources/general_model/pytorch_model.bin",
    "mgm/resources/phylogeny.csv"])`; `torch.load(..., weights_only=True)`;
    the tokenizer pickle is never loaded (the vocabulary is
    `phylogeny.csv`'s order after four special tokens, [V]).
  - Input rules: needs `var["genus"]`; one feature per genus (raise
    otherwise, naming `pp.tax_glom(tdata, "genus")`); unknown genera dropped
    with one warning naming the count; an all-zero sample, or one whose
    genera are all unknown, embeds `<bos><eos>` and warns.
  - Golden: try `microformer-mgm==0.5.8` in a Python 3.11 venv on a small
    genus table, store its embeddings under `tests/data/mgm/` (< 1 MB) with
    a NOTICE (MIT); compare to 1e-5. If it does not install, record that and
    test shape, determinism and the token sequence only.
  - The end-to-end test for exit gate 2 (design note 10): GlobalPatterns at
    genus level -> `bt.ml.embed(..., "mgm", inplace=True)` ->
    `obsm["X_mgm"]` of shape (26, 256), finite, equal across two calls.

**Open questions.**
- In-tree plugin vs a `biotapy-mgm` package (decision 9). Proposed:
  in-tree for 0.4; split out if MGM needs its own release cadence.
- `batch_size` meaning for MGM: samples per forward pass; default 64 keeps
  peak memory at 64 x 512 x 256 float32 activations per layer (about 34 MB)
  [UNVERIFIED for the whole forward].

---
## Slice 4D - Docs and release (outline)

**Goal:** the two exit-gate pages exist and run, the guide covers ML and
multi-omics, knowledge is current, and 0.4.0 is on PyPI.

**Tasks.**
- **4.6 `docs/tutorials/leak_free_cv.md`** (executed; exit gate 1): HMP2
  taxa (`bt.datasets.hmp2()["taxa"]`, relative abundances, so `CLR`'s
  pseudocount is set on their scale), IBD vs non-IBD, 5-fold stratified CV
  with a fixed `random_state`, logistic regression, ROC AUC. Sections: the
  pipeline; the prevalence filter outside vs inside (the measured near-zero
  difference, explained: it never sees labels); a supervised selection step
  outside vs inside, on real and shuffled labels (design note 9); what to
  put inside the pipeline (the guide's table). A test pins the quoted
  numbers through shared constants, as Phase 3 did for its tutorial.
- **4.6b embedding page** in the ML guide or `docs/tutorials/embeddings.md`:
  the MGM call as a non-executed block, the numbers from 4.4b's end-to-end
  test, how to write a plugin (the entry-point snippet a third-party
  `pyproject.toml` needs).
- **Multi-omics tutorial** (spec's Tutorials row): HMP2's `taxa` and
  `function` modalities through `io.to_mudata`; mmvec needs metabolites,
  and HMP2's metabolomics are not in `datasets.hmp2`. Proposed: defer the
  mmvec tutorial to 0.5 unless a small public paired dataset is found in
  4D [UNVERIFIED]; the guide's toy example stands for 0.4.
- **4.D1 Coming-from-R check**: `io.to_mudata` maps to
  `MultiAssayExperiment::MultiAssayExperiment`; nothing else in Phase 4 has
  an R equivalent; a test pins the row.
- **4.7 Knowledge**: `ml` Module concept refreshed for `to_torch`, `embed`
  and MGM; the plugin decision concept marked for the user's verification;
  roadmap index; log.
- **4.D2 Release 0.4.0** per [cut-a-release](/playbooks/cut-a-release.md):
  version bump, CHANGELOG, the wheel's `Provides-Extra` lines for `r`,
  `torch` and `mgm`; push, tag and publish only with explicit approval
  (R13.3).
- Benchmarks: none planned. `tl.mmvec` and the transformers delegate to
  scikit-bio and scikit-learn (R10.1: no measurement asks for one);
  `to_torch`'s per-row densify gets an asv benchmark only if 4B's large
  test shows a problem.

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
   writer. Filing an issue with mudata or treedata is a GitHub action on
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
    array or writes `obsm["X_<model>"]` (design note 8). A decision concept
    records it in 4C.
16. **Notebook content:** the leak-free CV page shows the prevalence filter's
    near-zero leak as measured and adds a supervised selection step to show
    a real one, on HMP2 (design note 9). The roadmap's literal version
    would show almost no difference.
17. **Exit gate 2 in CI, not in the docs build:** the `ml-extras` job's
    end-to-end MGM test is the proof; the docs page is a non-executed block
    quoting it (design note 10). Alternative: torch and transformers in the
    `doc` group (about 1 GB per docs build).
18. **Task changes:** new 4.A0 (`refactor(core)`), 4.B1 (CI), 4.4b (MGM),
    4.6b (embedding page), 4.D1 (Coming-from-R), 4.D2 (release); 4.2 moves
    to slice 4A; the mmvec tutorial on real paired metabolomics is deferred
    to 0.5 unless 4D finds a small public dataset.
19. **Frontmatter:** the new `description`; `paths` gains
    `src/biotapy/_core/**`; sources gain the mmvec paper, MGM and the
    BiomeGPT preprint.
20. **Branch pushes:** approve push, PR and merge on green for the Phase 4
    slice branches (R13.3); no standing approval covers Phase 4.

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
   | Spec: multi-omics tutorial | 4D, deferred in part (decision 18) |

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
   The only open values are the log section's `<date>` and each concept's
   `generated.at`, filled at commit time. Outline slices give a proposed
   answer to each open question and mark what is [UNVERIFIED].
3. **Type consistency.** Used throughout: `to_mudata(modalities:
   Mapping[str, AnnData]) -> MuData`; `mmvec(mdata, *, microbes="taxa",
   metabolites="metabolites", seed=None) -> pd.DataFrame`;
   `check_pseudocount(pseudocount: object) -> None`; `pseudocounted(X,
   pseudocount, *, func, columns=None)`; `PrevalenceFilter(*,
   min_prevalence: float = 0.1)` with `prevalence_`; `CLR(*, pseudocount: float = 0.5)`;
   `to_torch(adata, *, label_key=None, layer=None)`; `embed(adata, model, *,
   batch_size=64, inplace=False)`; plugin `embed(adata, *, batch_size)`;
   `obsm["X_<model>"]`.
4. **Review focus.** Items 1-4 name tests that exist in the rendered 4A
   code (checked by searching this file for each name); item 5 names tests
   4C must write under those names.
5. **Known residual risks.**
   - The gate counts were measured on Python 3.13 only; the CI matrix adds
     3.12, 3.14 and pre-release dependencies. scikit-learn's estimator
     checks are the likeliest to differ under pre-releases.
   - [UNVERIFIED]: mypy on a torch `Dataset` subclass without torch
     installed; uv's lock with a per-package torch index on every platform;
     MGM's sample pooling and a 3.11 parity fixture; transformers' wheels on
     3.14; `batch_size` memory for MGM.
   - Upstream: the weights live in a PyPI wheel the MGM authors control;
     MGM2 may supersede MGM before 0.4 ships.

[^spec]: Python Microbiome Toolkit development report, sections Positioning and Roadmap
[^mmvec]: Morton et al. 2019, Learning representations of microbe-metabolite interactions, Nature Methods
[^mgm]: MGM foundation model (microformer-mgm), HUST-NingKang-Lab
[^biomegpt]: BiomeGPT preprint
