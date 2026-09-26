---
type: Contract
title: Module boundaries
description: Layered package (_core at the bottom, pl/ml/da at the top); public API only through subpackage __all__; private topic files; shared helpers only in _core.
tags: [architecture, modularization]
status: stable
paths: ["src/biotapy/**", "pyproject.toml"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
---

# Statement

## Layout
```text
src/biotapy/
  __init__.py          # imports subpackages, __version__; nothing else
  _core/               # private kernel: tree access, groupby, validation, rng, optional imports
  io/  datasets/  pp/  tl/  fn/  da/  ml/  pl/
    __init__.py        # public API: explicit imports + __all__, no logic
    _<topic>.py        # one topic per file, e.g. pp/_glom.py, tl/_alpha.py
```
A subpackage is created in the phase that first fills it, never empty in advance.

## Layers (higher may import lower, never the reverse)
```text
pl | ml | da        (top, mutually independent)
tl | fn             (independent)
pp
datasets
io
_core               (imports only third-party)
```

## Rules
1. Users and tests reach functions only through `bt.<module>.<fn>`; private
   `_topic.py` files are imported only by their own subpackage `__init__.py`
   and by `_core` unit tests.
2. A private helper lives in the topic file if one public function uses it;
   in `_core` if two or more subpackages use it. Never copy-paste between modules.
3. Topic file: <= 300 lines, <= 6 public functions. Over that, split by topic.
4. Function body: <= 30 statements (ruff `PLR0915`), complexity <= 8 (`C901`).
5. No module-level side effects beyond constant definitions; no module-level
   import of an optional dependency.
6. Tests mirror source: `src/biotapy/pp/_glom.py` -> `tests/pp/test_glom.py`.

# Why
Layering keeps `_core` small and stable and stops `pl` from computing things
(`pl` reads slots; `tl` computes). Explicit `__all__` makes the public surface
reviewable in one file per module.

# Enforced by
- import-linter `layers` contract in `pyproject.toml` (Phase 0, task 0.4),
  run in CI as `lint-imports`.
- ruff `C901`, `PLR0915`, `TID252` (no relative imports beyond siblings).
- Not automated: file-length limit, checked in review.

# Binds
- [function-shape](/contracts/function-shape.md)
- [tree-access](/contracts/tree-access.md)
