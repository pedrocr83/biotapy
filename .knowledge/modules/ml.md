---
type: Module
title: ml
description: scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold, taking arrays, sparse matrices and DataFrames; to_torch, a PyTorch dataset over an AnnData's rows behind the extra torch; and embed, one embedding per sample from a model a plugin registers.
resource: /src/biotapy/ml/
paths: ["src/biotapy/ml/**"]
tags: [ml, scikit-learn, torch, plugins]
status: stable
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-10T21:29:09Z }
commit: 5cc503f
---

# Responsibility

Owns `bt.ml.*`, the top layer's machine-learning entry points. Today that is
two scikit-learn transformers (`_transformers.py`), `to_torch`
(`_torch.py`, extra `torch`) and `embed` (`_embed.py`), which runs a model
that a plugin registers ([embedding-plugins](/decisions/embedding-plugins.md)). Owns no reader and no table-level transform: `pp.filter_features` and
`pp.clr` stay `pp`'s, the transformers are their fold-safe forms.

# Entry points

- `_transformers.py:PrevalenceFilter` - keeps the features non-zero in at
  least `min_prevalence` of the training samples; `SelectorMixin` supplies
  `transform`, `get_support` and `get_feature_names_out`.
- `_transformers.py:CLR` - the centred log-ratio of each sample; stateless,
  so `fit` only validates.
- `_torch.py:to_torch` - a map-style `torch.utils.data.Dataset` over `X` or a
  layer, one float32 row per item, paired with an `obs` label when asked.
- `_embed.py:embed` - loads the callable registered under `model` in the
  entry-point group `biotapy.embeddings` (`_embed.py:GROUP`), calls it with
  the AnnData and returns its samples x dimensions array, or writes
  `obsm["X_<model>"]` with `inplace=True`.
- `_mgm.py:embed` - MGM, the reference plugin, registered by biotapy's own
  `pyproject.toml` as `mgm` (extra `mgm`): genus tokens from `var["genus"]`,
  MGM's pretrained GPT-2, the mean of the last hidden layer per sample.

# Invariants

- The transformers take arrays, sparse matrices and DataFrames, not AnnData:
  inside a `Pipeline` the splitter hands over `X`, not an AnnData (`to_torch`
  is the exception: it takes an AnnData, for its `obs` labels and layers).
  This is why `ml` is the one place classes are allowed
  ([function-shape](/contracts/function-shape.md), rules.md R3.6).
- Constructor options are keyword-only and stored unchanged; `fit` validates
  them (scikit-learn's convention), not `__init__`.
  `_transformers.py:PrevalenceFilter`, `_transformers.py:CLR`.
- `PrevalenceFilter`'s fitted state is `prevalence_` alone. A `fit` or refit
  that raises (bad option, no feature passing) leaves it unfitted, so
  `transform` raises `NotFittedError`, never a stale mask. `n_features_in_`
  is set by `validate_data` before the checks and is not the marker.
  `_transformers.py:PrevalenceFilter.fit`, `_transformers.py:PrevalenceFilter._get_support_mask`.
- `PrevalenceFilter` keeps the features `pp.filter_features(min_prevalence=)`
  keeps: it divides the count by the sample number (`prevalence >=
  min_prevalence`) rather than multiplying, because `7 / 25 >= 0.28` holds
  and `7 >= 0.28 * 25` does not. `_transformers.py:PrevalenceFilter._get_support_mask`.
- `CLR` gives `pp.clr`'s values (pinned by `tests/ml/test_transformers.py:test_clr_equals_pp_clr`);
  both build on `_core._composition.py:pseudocounted`
  ([core](/modules/core.md)). Output is dense float64, 8 bytes x samples x
  features. `_transformers.py:CLR.transform`.
- CSR stays sparse: `PrevalenceFilter` reads `indices` and `data` of the CSR
  and its output stays sparse; `CLR` accepts sparse input and densifies once
  inside `pseudocounted` (rules.md R6.2).
- `to_torch` densifies one row when its item is read, never the table, and
  holds `X` (or the layer) instead of copying it, so the AnnData must not be
  modified while the dataset is used; every item is a new tensor (label
  included), so editing it leaves the AnnData unchanged. A non-integer or bool index
  raises `TypeError`. Labels are converted once, at construction: category, string or bool -> int64 codes in category order, numeric
  -> float32; a missing label raises. Codes are per dataset (anndata drops a category a subset lacks), so build one dataset and split it with `torch.utils.data.Subset`; per-split datasets recode a class a split lacks. `_torch.py:to_torch`, `_torch.py:_labels`.
- `to_torch` validates its arguments before it imports torch, so its error
  tests run in every CI job, not only in `ml-extras`. `_torch.py:to_torch`.
- `embed` reads the entry points on every call and loads only the one asked
  for; a name no package registers raises `KeyError` listing those installed,
  a name two packages register raises `ValueError` naming both.
  `_embed.py:embed`.
- `embed` trusts no plugin: the result must be a plain `numpy.ndarray` (not
  a masked array or a matrix), 2-D with one row per sample and at least one
  column, float and finite, or it raises naming the plugin, before anything
  is written; a result that shares memory with `X`, a layer, an `obsm`, `varm`,
  `obsp` or `varp` entry or a top-level `uns` array (not one nested deeper) is
  copied. A plugin's own exception keeps its type and gains a note naming the
  plugin. The model name must match `[\w.-]+`, so that `obsm["X_<model>"]` is a plain key.
  `_embed.py:embed`, `_embed.py:_checked`.
- MGM's tokens are MGM's own (`mgm/src/MicroCorpus.py`, 0.5.8): a genus is
  read with MGM's regex `g__[A-Za-z0-9_]+` on `"g__" + genus`; features of
  one token are summed (`_core.sum_by`); relative abundance is taken over
  MGM's genera only, as MGM drops the others before dividing; a genus is kept
  when its standardised value `(rel - mean) / std` exceeds that of zero
  abundance, and the kept ones are sorted by it with pandas' own
  `sort_values(ascending=False)` on the vocabulary-ordered series (the order
  of tied genera follows numpy's quicksort, which can differ from MGM's numpy
  1.24 environment; GlobalPatterns and enterotype have no such tie); `<bos>` ... `<eos>` is cut to 512 tokens (the `<eos>` is
  lost past 510 genera). Vocabulary ids are the position in `phylogeny.csv`
  plus 4 (`<pad>`, `<mask>`, `<bos>`, `<eos>`); the pickled tokenizer is
  never loaded. `_mgm.py:_sentences`.
- Where MGM would drop a sample (no count in its vocabulary), the plugin
  embeds `<bos> <eos>` and warns naming it, because `ml.embed` needs one row
  per sample. Left-out features get one warning counting them.
  `_mgm.py:_sentences`.
- MGM checks `var["genus"]` and `X` (finite, non-negative) before it imports
  torch or transformers or downloads anything, so those errors are tested in
  every CI job. `_mgm.py:embed`.
- The weights come from microformer-mgm 0.5.8's wheel on files.pythonhosted.org
  (SHA-256 pinned), fetched through `_core.make_pooch` and unzipped by
  `pooch.Unzip` to `<cache>/microformer_mgm-0.5.8-py3-none-any.whl.unzip/`.
  `Unzip` re-extracts only a missing member, so each extracted file is also
  checked against its own SHA-256 (`_mgm.py:_SHA256_OF`); one that differs is
  deleted and extracted again once (logged), and one that still differs raises
  `ValueError` naming the `.unzip` folder to delete. `_mgm.py:_extracted_files`.
  `torch.load(..., weights_only=True)` reads the GPT2LMHeadModel checkpoint
  and its `transformer.` keys load into `transformers.GPT2Model` strictly
  (`lm_head.weight` dropped). `_mgm.py:embed`.
- One sample at a time, unpadded, under `torch.inference_mode()`, model in
  `eval()` (dropout off): measured on 520 GlobalPatterns-like profiles, batch
  1 took 8.8 s, batches of 8 and 64 took 17.0 s and 18.5 s and needed 196 MB
  and 1.6 GB more (8 threads). `_mgm.py:embed`.
- Inherited scikit-learn methods (`transform`, `fit_transform`,
  `get_support`, `get_feature_names_out`, `set_output`) are named in each
  class's `Notes`, because the class template leaves inherited members off
  the API page (rules.md R8.2).

# Dependencies

- [core](/modules/core.md): `as_csr`, `check_pseudocount`, `pseudocounted`.
- scikit-learn: `BaseEstimator`, `SelectorMixin`, `TransformerMixin`,
  `OneToOneFeatureMixin`, `validate_data`. scikit-bio: `clr`.
- torch (extra `torch`), only through `import_optional` inside `to_torch`
  ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).
- torch and transformers (extra `mgm`), only through `import_optional` inside
  `_mgm.py:embed`; pooch (`Unzip`) and `_core.make_pooch` for its weights.
- `ml` imports no sibling top-layer module (`pl`, `da`) and not `pp`; the
  link to `pp` is through `_core` only
  ([module-boundaries](/contracts/module-boundaries.md)).

# Verification

`uv run --group test pytest tests/ml` (the scikit-learn estimator checks, the
`pp` parity tests, a `Pipeline` cross-validation test, `to_torch`'s argument
errors, `embed` with fake plugins). `uv run --group test --extra torch pytest -m torch` runs the
`to_torch` tests and its docstring example, as CI's `ml-extras` job does (30-minute timeout,
`.github/workflows/test.yaml`). `BIOTAPY_DATA_DIR=.pooch uv run --group test --extra mgm pytest -m mgm` (`.pooch/` is git-ignored)
runs MGM against its own embeddings (`tests/data/mgm`), the end-to-end GlobalPatterns test (exit gate 2)
and `embed`'s docstring example, as the job's last step does. The pseudocount
warning's text is unit-tested in `tests/core/test_composition.py`.

The docs build runs `docs/tutorials/leak_free_cv.md` on HMP2;
`-m network` runs its cells again and compares them with `LEAK_FREE_CV_AUC`
(`tests/ml/test_transformers.py`). `docs/tutorials/embeddings.md` is not run by
the docs build (no torch): `-m mgm` runs its code, and a default-run test checks
that it quotes `GLOBAL_PATTERNS_LEFT_OUT` and `GLOBAL_PATTERNS_SHAPE`
(`tests/ml/test_mgm.py`).

# Gotchas

- `check_estimator` cannot run under `-W error::UserWarning`, the project's
  test mode; the tests use `parametrize_with_checks` with a
  `filterwarnings("ignore:pseudocount=")` mark, because the checks feed `CLR`
  data the pseudocount swamps. `tests/ml/test_transformers.py`.
- `CLR`'s pseudocount warning ends with `(ml.CLR)` (the `func=` name), but its
  reported location is scikit-learn's `sklearn/utils/_set_output.py`, not the
  caller: `transform` is wrapped by `_set_output`, and `warn_user` attributes
  a warning to the first frame outside biotapy. Match the message text, not
  the file. `_core/_warnings.py:warn_user`, `_transformers.py:CLR.transform`.
- `super().__sklearn_tags__()` returns an unannotated call; mypy's
  `untyped_calls_exclude` does not reach `super()`, so each class carries a
  targeted `# type: ignore[no-untyped-call]` rather than a global relaxation.
  `_transformers.py:PrevalenceFilter.__sklearn_tags__`.
- Prevalence is learned from the training fold only. Filtering before
  splitting lets the test samples pick the features; that is the reason the
  class exists. Both transformers have no R equivalent (`R equivalent: none`)
  and so no golden test; parity is against `pp`.
- `to_torch`'s `Dataset` subclass is defined inside the module-level
  `_dataset(table, labels)`, after `import_optional`: a module-level class
  would import torch with biotapy. A class defined in a function cannot be
  pickled, so it defines `__reduce__` returning `(_dataset, (table, labels))`;
  that is what lets `DataLoader(num_workers>0)` work under spawn and
  forkserver (macOS, Windows, Linux on Python 3.14).
  mypy does not follow torch (`follow_imports = "skip"` in `pyproject.toml`),
  so it is `Any` with or without the extra installed and the class line
  carries `# type: ignore[misc]` (subclassing `Any`) in both. The return
  annotation is the bare `"Dataset"`: sphinx-autodoc-typehints renders a
  subscripted `Dataset[Tensor]` from a `TYPE_CHECKING` import as a broken
  cross-reference, which `nitpicky` fails. `_torch.py:_dataset`.
- The docstring example needs torch: the root `conftest.py` gives the
  doctests of `biotapy.ml._torch` the marker `torch`, so the default run
  deselects them and `-m torch` runs them. The sdist ships the root `conftest.py` (`pyproject.toml`
  `build.targets.sdist.include`), because the sdist's tests need that hook; the
  wheel does not. `conftest.py:pytest_collection_modifyitems`.
- Without the extra, `import biotapy` never imports torch. With torch
  installed, scikit-bio 0.7.4 imports it at module level
  (`skbio.util._testing`, reached through `_core._tree`), so `'torch' in
  sys.modules` after `import biotapy` tells nothing there; the
  `import-without-extras` check is meaningful only on a torch-free
  environment, as in CI.
- `embed`'s tests install fake plugins as a real distribution: a
  `<name>-1.0.dist-info` with `METADATA` and `entry_points.txt` under
  `tmp_path`, prepended to `sys.path`, and a module object in `sys.modules`
  holding the callables. Discovery is then the same code path as for an
  installed package, and no test imports `ml/_embed.py`.
  `tests/ml/test_embed.py:_installer`.
- `embed`'s docstring example runs MGM, so the root `conftest.py` gives
  `biotapy.ml._embed`'s doctests the marker `mgm`, deselected by default like
  `torch`. `conftest.py:_EXTRA_DOCTESTS`.
- MGM's tests compare every call, the first of a process included, to MGM's own
  embeddings at `atol=1e-5` (biotapy and MGM 0.5.8 agree to 1.7e-6; 1e-3 would
  hide a changed last token, 7.7e-4). torch's CPU wheels bundle oneMKL 2024.2,
  whose VML caches the CPU type on its first call without a lock and publishes
  an unmapped value on the way (pytorch/pytorch#188792); a thread of a parallel
  first call that reads it runs AVX2's low-accuracy kernel on its chunk, moving
  an embedding by up to 2e-4 (about 15% of fresh processes at 32 threads, rarely
  at 4). `_mgm.py` fills the cache with one serial `torch.tanh(torch.zeros(1))`
  before the forward pass, and
  `tests/ml/test_mgm.py:test_a_fresh_process_s_first_embedding_is_its_second`
  runs 24 fresh 32-thread processes (Phase 4 decision 33).
  The mutations the
  parity test was checked against (wrong denominator, ascending sort,
  last-token pooling, dropout on, `z > 0`, vocabulary order) each fail it.
  `tests/ml/test_mgm.py:ATOL`.
- `microformer-mgm` cannot be installed beside biotapy (its pins), so the
  reference embeddings are written by `tests/mgm/export_reference.py` in a
  Python 3.11 environment of its own (`uv run --no-project --python 3.11 --with
  microformer-mgm==0.5.8 ...`); the script imports MGM inside its functions,
  because pytest's `--doctest-modules` imports every file under `tests/`.
  Set `HF_HOME` to a scratch directory in the command, else the old
  transformers writes under `~/.cache`. The fixture repeats "Blautia" in two
  features on purpose: summing them moves its rank in one sample, so a reader
  that does not sum features of one genus fails the parity test.
  `tests/mgm/export_reference.py:table`.
- **The two tutorials quote numbers that tests pin.** The prose of
  `docs/tutorials/leak_free_cv.md` and the quoted outputs of
  `docs/tutorials/embeddings.md` are constants in `tests/ml/test_transformers.py`
  and `tests/ml/test_mgm.py`; the root `tests/conftest.py:run_page` fixture runs
  a page's code cells and Python blocks (it strips MyST cell options and raises
  on an IPython magic), so the network and `mgm` tests check the page's own
  code. The leak-free page keeps its AUCs unrounded and shows three decimals, so
  its test compares them with `approx(abs=1e-3)`. A scikit-learn, torch or data
  change that moves a number fails one of them: update the constant and the page
  together.
- **The leak on HMP2 is small for the prevalence filter and large for a step
  that reads the labels**: the filter outside the pipeline moves the mean AUC
  by 0.001 (0.527 against 0.526), `SelectKBest(k=20)` by 0.147 on the real
  labels and by 0.223 on shuffled ones (0.544 against 0.767). The honest
  pipeline is near chance: a first stool sample's species barely separate IBD
  from non-IBD there. The tutorial is about leakage, not a classifier.
  `tests/ml/test_transformers.py:LEAK_FREE_CV_AUC`.
- **GlobalPatterns writes *Candidatus* genera as one word** (Greengenes:
  `CandidatusPelagibacter`), MGM's vocabulary as `Candidatus_Pelagibacter`. Of
  the 96 GlobalPatterns genus-level features MGM leaves out, 28 are such names,
  23 of which MGM's vocabulary holds with the underscore (counted against
  `phylogeny.csv` of mgm 0.5.8; no test pins the 28 and 23). MGM's own regex
  reads the name as written, so biotapy leaves them out as MGM would (Phase 4
  decision 30); `tests/ml/test_mgm.py:test_embeds_global_patterns_end_to_end`
  pins the 96.
- anndata 0.13 lists `X` as `layers[None]`, so `list(adata.layers)` holds
  `None` even when no layer was added; `to_torch`'s missing-layer error does
  not list the layers. `_torch.py:_table`.
