---
type: Module
title: Core (`_core`)
description: Private kernel package - sparse group math, taxonomic rank order, x_kind/provenance/slot rules, and the sole gateway to TreeData and networkx.
resource: /src/biotapy/_core/
paths: ["src/biotapy/_core/**"]
tags: [core, kernel]
generated: { by: claude-code/claude-sonnet-5, at: 2026-09-27T10:23:10Z }
commit: 43d6efb
status: stable
---

# Responsibility

Owns: sparse matrix kernels (`_matrix.py`), the canonical taxonomic rank order
and lineage/rank-column parsing (`_taxonomy.py`), the `x_kind`/provenance/
`feature_subset` slot rules (`_slots.py`), TreeData construction and the only
import of `treedata`/`networkx` in the package (`_tree.py`), the single
user-facing warning entry point (`_warnings.py`), lazy optional-dependency
import (`_optional.py`), and the single RNG entry point (`_rng.py`).

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
- `_taxonomy.py:normalize_ranks` - canonicalize rank column names (aliases,
  lowercasing) and values (strip `k__`-style prefixes, map `""`/`"NA"`/a bare
  prefix to NaN) to the pandas `str` dtype; used directly by `io.read_dada2`
  and indirectly, via `split_lineage`, by `io.read_biom`/`read_qiime2`.
- `_taxonomy.py:split_lineage` - parse a `;`-separated lineage string
  (Greengenes `k__`, RESCRIPT/SILVA `d__`, SILVA `D_0__`, or unprefixed) into
  rank columns, then run them through `normalize_ranks`.
- `_slots.py:x_kind` / `require_counts` - read, or enforce, what `X` holds.
- `_slots.py:infer_x_kind` - classify a freshly read matrix as `"counts"`
  (every value a whole number), `"relative"` (every nonzero row sums to 1
  within `RELATIVE_TOLERANCE`), or `"abundance"`, for readers whose file
  format does not record `x_kind` itself.
- `_slots.py:add_provenance` - append one provenance entry.
- `_slots.py:feature_subset` - the only place that implements the
  Propagation table in [data-model-slots](/contracts/data-model-slots.md).
- `_tree.py:make_treedata` / `get_tree` / `tree_from_edges` - build or read a
  TreeData's phylogeny; the only functions allowed to touch `treedata` or
  `networkx` ([tree-access](/contracts/tree-access.md)). `make_treedata`
  also casts `obs`/`var` ids to unique, non-missing strings and, given a
  tree, aligns it with the table (see Invariants).
- `_tree.py:tree_from_newick` - parse one Newick string via scikit-bio into
  the same tip-named, uniquely-labelled-internal-node graph shape as
  `tree_from_edges`.
- `_tree.py:tree_from_phylo` - build a tree from an ape `phylo` edge matrix
  (R's `phy_tree` slot: 1-based, tips numbered `1..len(tips)`), naming
  internal nodes with the same collision-free `n<i>` scheme as
  `tree_from_newick`; used by `io.read_phyloseq`.
- `_tree.py:tree_tips` - the tree's leaf names (nodes with no children).
- `_tree.py:relabel_tips` - rename a subset of a tree's nodes (e.g. sequence
  -> ASV id), raising rather than silently merging nodes on a name collision.
- `_warnings.py:warn_user` - the single `UserWarning` entry point, attributed
  to the first stack frame outside biotapy; shared by `_tree.py` (tree/table
  mismatch) and `io/_join.py` (partial join).
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
- `make_treedata` requires every `obs`/`var` id to be present and, once cast
  to `str`, unique; a missing (NaN/None) id or a duplicate raises
  `ValueError` naming the axis and (for duplicates) up to 5 ids.
  `_tree.py:_with_str_ids`.
- Given a tree, `make_treedata` keeps only features that are tree tips and
  prunes tips outside the table (ancestors kept); a partial overlap gives one
  `UserWarning` naming both counts, and no overlap at all raises `ValueError`
  with example ids from each side. `_tree.py:_align_tree`.
- Rank columns (`_taxonomy.py:normalize_ranks`) are the pandas `str` dtype,
  never `object`: anndata's h5ad/h5td writer raises on an all-NaN `object`
  column, which a reader whose `species` rank is entirely missing would
  otherwise produce.

# Dependencies

None inside biotapy. Imports only third-party packages: `numpy`, `scipy`,
`pandas`, `anndata`, and, in `_tree.py` only, `treedata`/`networkx` and
scikit-bio (Newick parsing, `_tree.py:tree_from_newick`). scikit-bio is a
real cost at import time: measured at commit 43d6efb, `import biotapy` takes
~1.0-1.1s, of which roughly half (~0.5s) is scikit-bio, found by diffing
against importing biotapy's other runtime dependencies alone. A lazy
(function-local) import is deferred until Task 1.21's asv benchmark measures
it properly (rules.md R10.1: no optimization without a measurement).

# Verification

`uv run --group test pytest tests/core -q` plus `uvx prek run --all-files`.

# Gotchas

- `mypy --strict` type-checks `treedata` and `skbio` through
  `follow_untyped_imports` (`pyproject.toml` `[[tool.mypy.overrides]]`), not
  `ignore_missing_imports`: neither ships `py.typed`, and
  `ignore_missing_imports` would silently make their types `Any` everywhere
  they are used.
- `as_csr` takes `X: object`, not a sparse/array-like union, specifically so
  callers can pass `adata.X` directly without a cast - anndata types
  `AnnData.X` as a private union the old, narrower parameter type could not
  accept (commit fc26baa). Do not narrow the annotation back.
