---
type: Module
title: datasets
description: In-memory and pooch-cached example TreeData objects for docs, doctests and tests.
resource: /src/biotapy/datasets/
paths: ["src/biotapy/datasets/**"]
tags: [datasets]
generated: { by: claude-code/claude-sonnet-5, at: 2026-09-27T17:28:27Z }
commit: 54b8ef2
status: stable
---

# Responsibility

Owns `bt.datasets.*` example-data loaders: `toy()`, built entirely in memory,
and `global_patterns()`/`enterotype()`/`esophagus()`, downloaded once from
phyloseq's repository and cached with pooch (Task 1.11, esophagus Task
1.15b).

# Entry points

- `_toy.py:toy` - six samples (`obs["group"]` A/B, three each) x eight
  features, with kingdom..genus taxonomy (`f8` unassigned at genus) and a
  phylogeny in `vart["phylo"]`.
- `_remote.py:global_patterns` - GlobalPatterns: 26 samples x 19,216 OTUs,
  with taxonomy and a tree, read through `bt.io.read_phyloseq`.
- `_remote.py:enterotype` - enterotype: 280 samples x 553 genera, as relative
  abundances, no tree.
- `_remote.py:esophagus` - esophagus: 3 samples x 58 OTUs, with a tree, no
  taxonomy or sample data.

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
- `global_patterns()`/`enterotype()`/`esophagus()` never construct their
  `pooch.Pooch` (`_remote.py:_pooch`, `@cache`) until first called, so
  importing `biotapy.datasets` does no network or disk work (R4.7). No
  third-party data file is committed (stage2-constraints): all three are
  downloaded at run time, pinned to one phyloseq commit and a SHA-256 in
  `_remote.py:_REGISTRY`.

# Dependencies

- [core](/modules/core.md): `make_treedata`, `tree_from_edges`, `TreeData`.
- [io](/modules/io.md): `read_phyloseq`, used by `_remote.py`'s three loaders.
- `pooch` (runtime, Task 1.11): fetches and caches `GlobalPatterns.RData`/
  `enterotype.RData`/`esophagus.RData`; `BIOTAPY_DATA_DIR` overrides its
  cache directory (`pooch.create(..., env="BIOTAPY_DATA_DIR")`).

# Verification

`uv run --group test pytest tests/datasets -q` plus `uvx prek run --all-files`
for `toy()` and `_remote.py`'s offline test; add
`BIOTAPY_DATA_DIR=<dir> uv run --group test pytest -m network tests/datasets -q`
to actually exercise the three downloads (CI's dedicated `network` job runs
`-m "network or golden"`).

# Gotchas

- `toy()` is deliberately small and hand-built so exact values (which
  feature is each phylum's archetype, which samples share a genus, ...) can
  be pinned in tests and docstring examples across `pp` and (later) `tl`/
  `pl`. Changing its numbers silently breaks those pinned expectations
  without touching `toy()`'s own tests - check `tests/pp/test_glom.py` and
  every doctest calling `bt.datasets.toy()` before changing it.
- `global_patterns()`/`enterotype()`/`esophagus()`'s doctest examples are
  `# doctest: +SKIP`: running them for real would download phyloseq's (AGPL-3) data
  during `pytest --doctest-modules`, which R11.4 and stage2-constraints both
  forbid outside the dedicated network job. `tests/datasets/test_remote.py`
  covers the loaders offline instead, by monkeypatching the private
  `_remote._fetch` to return a local fixture - the one place R4.9's "tests
  call the public API" rule is deliberately broken, for network isolation.
