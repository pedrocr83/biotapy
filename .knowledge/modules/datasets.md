---
type: Module
title: datasets
description: In-memory and (future) cached example TreeData objects for docs, doctests and tests.
resource: /src/biotapy/datasets/
paths: ["src/biotapy/datasets/**"]
tags: [datasets]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T14:08:28Z }
commit: 0fdbd4d
status: stable
---

# Responsibility

Owns `bt.datasets.*` example-data loaders. Today: `toy()`, built entirely in
memory. Downloaded reference datasets (`global_patterns`, `enterotype`, via
pooch) are Slice 1B (Task 1.11), not yet written.

# Entry points

- `_toy.py:toy` - six samples (`obs["group"]` A/B, three each) x eight
  features, with kingdom..genus taxonomy (`f8` unassigned at genus) and a
  phylogeny in `vart["phylo"]`.

# Invariants

- `toy()` never touches the network or disk, so docs and tests never
  download anything. `_toy.py:toy`
- Its counts are laid out as contiguous taxonomic blocks per phylum
  (`_toy.py:_TAXONOMY`, `_toy.py:_COUNTS`); it cannot exercise interleaved
  taxa, which is why `tests/pp/test_glom.py` builds its own fixture for that
  case rather than reusing `toy()`.
- Round-trips through `TreeData.write_h5td`/`treedata.read_h5td` unchanged in
  `X`, tree edge lengths and provenance
  (`tests/datasets/test_toy.py:test_toy_round_trips_through_h5td`); the
  provenance list comes back as an `ndarray` of `str`, not a `list`.

# Dependencies

- [core](/modules/core.md): `make_treedata`, `tree_from_edges`, `TreeData`.

# Verification

`uv run --group test pytest tests/datasets -q` plus `uvx prek run --all-files`.

# Gotchas

- `toy()` is deliberately small and hand-built so exact values (which
  feature is each phylum's archetype, which samples share a genus, ...) can
  be pinned in tests and docstring examples across `pp` and (later) `tl`/
  `pl`. Changing its numbers silently breaks those pinned expectations
  without touching `toy()`'s own tests - check `tests/pp/test_glom.py` and
  every doctest calling `bt.datasets.toy()` before changing it.
