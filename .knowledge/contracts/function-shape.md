---
type: Contract
title: Public function shape
description: One task = one public function `verb(data, required, *, options) -> result`, fully typed, keyword-only options, seeded randomness, NumPy docstring with a parseable R-equivalent line.
tags: [api, conventions, docs]
status: stable
paths: ["src/biotapy/**/*.py"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T16:44:46Z }
commit: 31aa11d
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Statement

Every public function in `io`, `datasets`, `pp`, `tl`, `fn`, `da`, `ml`, `pl`:

1. **Signature**: `verb(data, <required args>, *, <options>) -> <result>`.
   - `data` is annotated with the widest type that works: `AnnData` when no tree
     is needed, `TreeData` when `vart` is read, `MuData` for multi-modal.
     `da.consensus` takes a sequence of `da` result tables and `pl.consensus` the
     table `da.consensus` returns instead: they combine and draw results, not data.
   - Everything after the required arguments is keyword-only (`*`).
   - `ml`'s scikit-learn transformers are classes (rules.md R3.6), not
     functions: their options are keyword-only constructor arguments that `fit`
     validates (scikit-learn's convention, which `check_estimator` tests),
     and the class docstring carries the skeleton below. `ml.to_torch` is a
     function: its `Dataset` subclass is the one class, built inside
     `ml/_torch.py:_dataset` (rules.md R3.6).
   - No `**kwargs` pass-through, except a documented `plot_kwargs` in `pl`.
2. **Return and mutation**: per [pure-by-default](/decisions/pure-by-default.md).
3. **Randomness**: any stochastic function takes
   `seed: int | np.random.Generator | None = None`, converted once with
   `biotapy._core.as_generator(seed)`. No global RNG state is read or set.
4. **Types**: full hints on parameters and return; no `Any`, no untyped
   `dict`/`list` in public signatures; `Literal[...]` for string options.
5. **Validation**: inputs are checked at the public boundary with
   `biotapy._core` validators; errors are `ValueError`/`KeyError`/`TypeError`
   naming the offending argument. Private helpers do not re-validate.
6. **Docstring** (NumPy style, enforced by ruff `D` rules):

   ```text
   One-line summary ending with a period.

   Parameters / Returns sections.

   Notes
   -----
   R equivalent: ``phyloseq::tax_glom``, ``mia::agglomerateByRank``
   Guide: :doc:`/guide/aggregation`

   Examples
   --------
   >>> import biotapy as bt
   >>> tdata = bt.datasets.toy()
   >>> bt.pp.tax_glom(tdata, "phylum").n_vars
   3
   ```

   - Exactly one line starting `R equivalent:`; comma-separated
     ``pkg::fn`` items, or `none`. The "Coming from R" page is generated from it.
   - Examples use `bt.datasets.toy()`, or `bt.datasets.toy_humann()` for
     function tables (both built in memory, no download; a function that
     reads a file may write a small temp file in its example), and run under
     doctest in CI.

# Why
- Keyword-only options let parameters be added without breaking callers.
- A seeded generator is the only way stochastic results are reproducible and testable.
- A parseable `R equivalent:` line keeps the migration table generated, never hand-written.

# Enforced by
- ruff `D`, `PLR0913`, `PLR0917` and mypy strict (Phase 0, task 0.4).
- `tests/test_docstrings.py` checks every public function's docstring. It must have exactly one
  `R equivalent:` line, parsed by `docs/extensions/coming_from_r.py`, with `Guide:` on the next
  line, and an `Examples` section. The same parser writes the Coming-from-R table at each docs
  build, with `docs/_data/r_idioms.toml` for the R calls that are plain AnnData code (task 1.19).
- Doctests in CI (`--doctest-modules` in pytest config, Phase 0 task 0.4).
- Purity: each public `pp`, `tl` and `pl` function has a test asserting its input is unchanged
  (`assert_unchanged`, `tests/conftest.py`); `pl.scree` shares `pl.ordination`'s.

# Binds
- [module-boundaries](/contracts/module-boundaries.md)
- [data-model-slots](/contracts/data-model-slots.md)
