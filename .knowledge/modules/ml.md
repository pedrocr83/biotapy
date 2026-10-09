---
type: Module
title: ml
description: scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold, taking arrays, sparse matrices and DataFrames; and to_torch, a PyTorch dataset over an AnnData's rows behind the extra torch.
resource: /src/biotapy/ml/
paths: ["src/biotapy/ml/**"]
tags: [ml, scikit-learn, torch]
status: stable
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T12:00:00Z }
commit: 7a9c07a
---

# Responsibility

Owns `bt.ml.*`, the top layer's machine-learning entry points. Today that is
two scikit-learn transformers (`_transformers.py`) and `to_torch`
(`_torch.py`, extra `torch`); `ml.embed` (slice 4C) is planned in
[phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md) and does not exist
yet. Owns no reader and no table-level transform: `pp.filter_features` and
`pp.clr` stay `pp`'s, the transformers are their fold-safe forms.

# Entry points

- `_transformers.py:PrevalenceFilter` - keeps the features non-zero in at
  least `min_prevalence` of the training samples; `SelectorMixin` supplies
  `transform`, `get_support` and `get_feature_names_out`.
- `_transformers.py:CLR` - the centred log-ratio of each sample; stateless,
  so `fit` only validates.
- `_torch.py:to_torch` - a map-style `torch.utils.data.Dataset` over `X` or a
  layer, one float32 row per item, paired with an `obs` label when asked.

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
- `ml` imports no sibling top-layer module (`pl`, `da`) and not `pp`; the
  link to `pp` is through `_core` only
  ([module-boundaries](/contracts/module-boundaries.md)).

# Verification

`uv run --group test pytest tests/ml` (the scikit-learn estimator checks, the
`pp` parity tests, a `Pipeline` cross-validation test, `to_torch`'s argument
errors). `uv run --group test --extra torch pytest -m torch` runs the
`to_torch` tests and its docstring example, as CI's `ml-extras` job does (30-minute timeout,
`.github/workflows/test.yaml`). The pseudocount
warning's text is unit-tested in `tests/core/test_composition.py`.

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
- anndata 0.13 lists `X` as `layers[None]`, so `list(adata.layers)` holds
  `None` even when no layer was added; `to_torch`'s missing-layer error does
  not list the layers. `_torch.py:_table`.
