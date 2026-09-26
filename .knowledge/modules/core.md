---
type: Module
title: Core (`_core`)
description: Private kernel package - sparse group math, taxonomic rank order, x_kind/provenance/slot rules, and the sole gateway to TreeData and networkx.
resource: /src/biotapy/_core/
paths: ["src/biotapy/_core/**"]
tags: [core, kernel]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T14:08:28Z }
commit: 0fdbd4d
status: stable
---

# Responsibility

Owns: sparse matrix kernels (`_matrix.py`), the canonical taxonomic rank order
(`_taxonomy.py`), the `x_kind`/provenance/`feature_subset` slot rules
(`_slots.py`), TreeData construction and the only import of `treedata`/
`networkx` in the package (`_tree.py`), lazy optional-dependency import
(`_optional.py`), and the single RNG entry point (`_rng.py`).

Does NOT own any public verb (`bt.pp.*`, `bt.tl.*`, ...) - those live one
layer up, per [module-boundaries](/contracts/module-boundaries.md). `_core`
is the bottom layer: every other biotapy module may import it; it imports
none of them back.

# Entry points

- `_matrix.py:as_csr` - normalize any array-like or sparse input (including
  `AnnData.X`) to CSR, without copying one that already is.
- `_matrix.py:sum_by` / `_matrix.py:argmax_by` - grouped column sum and
  grouped argmax by integer group codes, negative codes dropped; today's only
  caller is `pp.tax_glom`.
- `_taxonomy.py:split_ranks` - split `var`'s taxonomy columns at a target
  rank, raising `KeyError` naming the rank when it is absent.
- `_slots.py:x_kind` / `require_counts` - read, or enforce, what `X` holds.
- `_slots.py:add_provenance` - append one provenance entry.
- `_slots.py:feature_subset` - the only place that implements the
  Propagation table in [data-model-slots](/contracts/data-model-slots.md).
- `_tree.py:make_treedata` / `get_tree` / `tree_from_edges` - build or read a
  TreeData's phylogeny; the only functions allowed to touch `treedata` or
  `networkx` ([tree-access](/contracts/tree-access.md)).
- `_optional.py:import_optional` - lazy import for a heavy extra, raising an
  `ImportError` that names the extra to install.
- `_rng.py:as_generator` - the single entry point that turns a seed into a
  `np.random.Generator` without touching global RNG state.

# Invariants

- `as_csr` may return an object sharing buffers with its input; a caller must
  never mutate the result in place (rules.md R6.2/R3.3 - matters for the
  future `pp.rarefy`). `_matrix.py:as_csr`
- `add_provenance` mutates its `adata` argument in place by design. Call it
  on a function's output copy, never on the caller's input.
  `_slots.py:add_provenance`
- `feature_subset` drops `layers`, `obsm`, `obsp`, `varm`, `varp` and every
  `uns` key except `biotapy`, but does not itself touch `vart`.
  `_slots.py:feature_subset`, `_slots.py:DERIVED_SLOTS`. TreeData's own
  subsetting prunes the tree to the kept leaves plus their ancestors, so
  unary internal nodes survive with their original edge lengths
  (`treedata._utils.subset_tree`, exercised by
  `tests/pp/test_glom.py:test_phylum_keeps_most_abundant_member_in_original_order`).
- Provenance entries are JSON strings, not dicts, because h5ad/h5td cannot
  store a list of dicts; after an `.h5td` round-trip the list comes back as
  an `ndarray` of `str`. `_slots.py:add_provenance`,
  `tests/conftest.py:_assert_unchanged`.
- Only `_tree.py` imports `treedata` or `networkx`; enforced by ruff
  `TID251` with a per-file exemption for that module
  ([tree-access](/contracts/tree-access.md)).

# Dependencies

None inside biotapy. Imports only third-party packages: `numpy`, `scipy`,
`pandas`, `anndata`, and, in `_tree.py` only, `treedata`/`networkx`.

# Verification

`uv run --group test pytest tests/core -q` plus `uvx prek run --all-files`.

# Gotchas

- `mypy --strict` type-checks `treedata` through `follow_untyped_imports`
  (`pyproject.toml` `[[tool.mypy.overrides]]`), not `ignore_missing_imports`:
  treedata ships no `py.typed`, and `ignore_missing_imports` would silently
  make the `TreeData` type `Any` everywhere it is used.
- `as_csr` takes `X: object`, not a sparse/array-like union, specifically so
  callers can pass `adata.X` directly without a cast - anndata types
  `AnnData.X` as a private union the old, narrower parameter type could not
  accept (commit fc26baa). Do not narrow the annotation back.
