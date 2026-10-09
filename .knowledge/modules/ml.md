---
type: Module
title: ml
description: scikit-learn transformers over a samples x features table - PrevalenceFilter and CLR - so preprocessing is fitted inside each cross-validation fold; they take arrays, sparse matrices and DataFrames, never AnnData.
resource: /src/biotapy/ml/
paths: ["src/biotapy/ml/**"]
tags: [ml, scikit-learn]
status: stable
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T01:55:36Z }
commit: 9883786
---

# Responsibility

Owns `bt.ml.*`, the top layer's machine-learning entry points. Today that is
two scikit-learn transformers (`_transformers.py`); `ml.to_torch` (slice 4B)
and `ml.embed` (slice 4C) are planned in
[phase-4-ml-multiomics](/roadmap/phase-4-ml-multiomics.md) and do not exist
yet. Owns no reader and no table-level transform: `pp.filter_features` and
`pp.clr` stay `pp`'s, the transformers are their fold-safe forms.

# Entry points

- `_transformers.py:PrevalenceFilter` - keeps the features non-zero in at
  least `min_prevalence` of the training samples; `SelectorMixin` supplies
  `transform`, `get_support` and `get_feature_names_out`.
- `_transformers.py:CLR` - the centred log-ratio of each sample; stateless,
  so `fit` only validates.

# Invariants

- They take arrays, sparse matrices and DataFrames, not AnnData: inside a
  `Pipeline` the splitter hands over `X`, not an AnnData. This is why `ml`
  is the one place classes are allowed
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
- Inherited scikit-learn methods (`transform`, `fit_transform`,
  `get_support`, `get_feature_names_out`, `set_output`) are named in each
  class's `Notes`, because the class template leaves inherited members off
  the API page (rules.md R8.2).

# Dependencies

- [core](/modules/core.md): `as_csr`, `check_pseudocount`, `pseudocounted`.
- scikit-learn: `BaseEstimator`, `SelectorMixin`, `TransformerMixin`,
  `OneToOneFeatureMixin`, `validate_data`. scikit-bio: `clr`.
- `ml` imports no sibling top-layer module (`pl`, `da`) and not `pp`; the
  link to `pp` is through `_core` only
  ([module-boundaries](/contracts/module-boundaries.md)).

# Verification

`uv run --group test pytest tests/ml` (the scikit-learn estimator checks, the
`pp` parity tests, a `Pipeline` cross-validation test). The pseudocount
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
