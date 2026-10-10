---
type: Module
title: datasets
description: In-memory and pooch-cached example data for docs, doctests and tests - TreeData objects, a HUMAnN-style function MuData, the HMP2 cohort as a three-modality MuData, mmvec's soil biocrust example as a microbes-and-metabolites MuData, and the ENZYME hierarchy as an edge table.
resource: /src/biotapy/datasets/
paths: ["src/biotapy/datasets/**"]
tags: [datasets]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-10T18:52:47Z }
commit: 38f9379
status: stable
---

# Responsibility

Owns `bt.datasets.*` example-data loaders: `toy()` and `toy_humann()`, built
entirely in memory; `global_patterns()`/`enterotype()`/`esophagus()`,
downloaded once from phyloseq's repository and cached with pooch (Task 1.11,
esophagus Task 1.15b); `enzyme()`, the ENZYME EC hierarchy downloaded from
ExPASy the same way; and `hmp2()`, the HMP2 (IBDMDB) cohort's pathway and
taxon tables (Task 2.10); and `biocrust()`, mmvec's soil example of microbes
and metabolites (Phase 4 Task 4.6c). The return type varies: a TreeData, a MuData
(`toy_humann`, `hmp2`, `biocrust`) or a `pandas.DataFrame` (`enzyme`).

# Entry points

- `_toy.py:toy` - six samples (`obs["group"]` A/B, three each) x eight
  features, with kingdom..genus taxonomy (`f8` unassigned at genus) and a
  phylogeny in `vart["phylo"]`.
- `_toy.py:toy_humann` - the toy samples' gene families regrouped to EC numbers
  as HUMAnN writes them: a MuData with `function` and `function_by_taxon`
  modalities, `x_kind` `rpk`, and `obs["group"]` as in `toy()`. Built through
  `_core.make_function_mudata`; the doctest and test fixture for `fn`.
- `_enzyme.py:enzyme` - the ENZYME hierarchy as an edge table for
  `fn.func_glom`: columns `child`, `parent`, `level` (`class`, `subclass`,
  `subsubclass`), `parent_name`, one row per EC number and ancestor, with
  `attrs["source"]` (the release read) and `attrs["license"]` (`"CC BY 4.0"`).
  Its text columns are the `str` dtype, the same schema `fn.load_hierarchy`
  gives (`_enzyme.py:_ancestors`).
- `_hmp2.py:hmp2` - the HMP2 inflammatory bowel disease cohort: each of the
  130 participants' first stool metagenome, as a MuData with `function` and
  `function_by_taxon` (HUMAnN 3 pathways in CPM, read by `io.read_humann`) and
  `taxa` (MetaPhlAn 3 species, read by `io.read_metaphlan`, no tree). A
  participant's metagenomes are ordered by `Participant ID`, `week_num`,
  `visit_num`, then `External ID`, and the first is kept (`_hmp2.py:ORDER`);
  `visit_num` only orders and is not in `obs`, and a missing one sorts last.
  Seven metadata columns (`_hmp2.py:COLUMNS`) sit in the global `obs` and are
  pushed into every modality; `diagnosis` is categorical `nonIBD`, `UC`, `CD`.
- `_biocrust.py:biocrust` - mmvec's soil example (`biocore/mmvec`,
  `examples/soils`): 466 microbe counts in 20 samples and 85 metabolite
  intensities in 19, read by `io.read_biom` and combined by `io.to_mudata` as
  `taxa` and `metabolites` over the 19 shared samples, in the microbe table's
  order. The loader intersects the samples itself, so `to_mudata` does not
  warn about the one microbe sample (`9hr_late`) without metabolites.
- `_remote.py:global_patterns` - GlobalPatterns: 26 samples x 19,216 OTUs,
  with taxonomy and a tree, read through `bt.io.read_phyloseq`.
- `_remote.py:enterotype` - enterotype: 280 samples x 553 genera, as relative
  abundances, no tree.
- `_remote.py:esophagus` - esophagus: 3 samples x 58 OTUs, with a tree, no
  taxonomy or sample data.

# Invariants

- `toy()` and `toy_humann()` never touch the network or disk, so docs and
  tests never download anything. `_toy.py:toy`, `_toy.py:toy_humann`
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

- `hmp2()`'s three files (`pathabundances_3.tsv.gz`,
  `taxonomic_profiles_3.tsv.gz`, `hmp2_metadata_2018-08-20.csv`) are pinned
  to SHA-256 hashes in `_remote.py:_REGISTRY`, with absolute IBDMDB URLs in
  `_URLS`. It fetches all three before parsing any, so a download error comes
  before a parse error (`_hmp2.py:FILES`). The metadata is the study's whole
  sample metadata; only its 1,638 stool metagenomes are used, and they share
  their ids with both tables (checked 2026-10-05), so the selected samples
  index every table. `_hmp2.py:hmp2`
- `biocrust()`'s two files are pinned to one commit of `biocore/mmvec`
  (`88ca33b`, `_remote.py:_MMVEC`) and to SHA-256 hashes; both are fetched
  before either is parsed (`_biocrust.py:FILES`). Their registry names carry a
  `biocrust_` prefix because the repository calls them `microbes.biom` and
  `metabolites.biom`.
- `enzyme()` is the one download that carries no pinned hash: ENZYME keeps only
  its current release online, so `_remote.py:_REGISTRY` lists `enzyme.dat` and
  `enzclass.txt` with `None`. The first download is cached for good, and
  `enzyme()` records the release in `attrs["source"]`. It raises `ValueError`
  (not a silent fallback) when either file lacks its `Release` line or the two
  releases differ. `_enzyme.py:enzyme`.

# Dependencies

- [core](/modules/core.md): `make_treedata`, `make_function_mudata`,
  `tree_from_edges`, `TreeData`.
- [io](/modules/io.md): `read_phyloseq`, used by `_remote.py`'s three loaders;
  `read_humann` and `read_metaphlan`, used by `_hmp2.py:hmp2`; `read_biom` and
  `to_mudata`, used by `_biocrust.py:biocrust`.
- `mudata`: `hmp2` builds its MuData and pushes the global `obs` into the
  modalities (`MuData.push_obs`).
- `pooch` (runtime, Task 1.11): fetches and caches `GlobalPatterns.RData`/
  `enterotype.RData`/`esophagus.RData` and, from the ExPASy FTP site,
  `enzyme.dat`/`enzclass.txt`, and from the IBDMDB's Globus endpoint the three
  HMP2 files, and from GitHub mmvec's two BIOM files (`_remote.py:_URLS`); `BIOTAPY_DATA_DIR` overrides its
  cache directory (the pooch is built by `_core._download.py:make_pooch`,
  [core](/modules/core.md), which `ml/_mgm.py` shares).

# Verification

`uv run --group test pytest tests/datasets -q` plus `uvx prek run --all-files`
for `toy()` and `_remote.py`'s offline test; add
`BIOTAPY_DATA_DIR=<dir> uv run --group test pytest -m network tests/datasets -q`
to actually exercise the downloads (CI's dedicated `network` job runs
`-m "network or golden"`). The docs job runs `docs/tutorials/function.md`
and `docs/tutorials/leak_free_cv.md`, which call `hmp2()`, and
`docs/tutorials/multiomics.md`, which calls `biocrust()`.

# Gotchas

- Delete both cached ENZYME files (not one) to take a newer release: with one
  stale and one fresh, `enzyme()` raises its two-releases `ValueError`.
  `_enzyme.py:enzyme`. A new ENZYME release can change the table silently
  because nothing pins it; the file format checks are the only guard.
- `enzyme()`'s doctest is `# doctest: +SKIP`, like the phyloseq loaders'.
  `tests/datasets/test_enzyme.py` tests the loader against local excerpts under
  `tests/data/enzyme/` through an `offline` fixture, and the real download is
  covered by the `network` marker.

- `hmp2()` reads the whole pathway table before selecting samples:
  `read_humann` builds a dense 22,113 x 1,638 `float64` array (about 290 MB),
  and the call peaks at about 1 GB resident (measured 1,029-1,033 MB), so it
  needs that much free memory. A cold cache downloads 23 MB, measured at
  10-25 s, which is why `docs/conf.py` sets `nb_execution_timeout = 300`.
- MuData converts the global `obs` with pandas' `convert_dtypes`, so
  `hmp2()`'s `week_num` is `int64`, `consent_age` the nullable `Int64` and
  text columns pandas' `string` dtype, not the `str` dtype the readers use
  (mudata 0.4.1). The h5mu round trip works
  (`tests/datasets/test_hmp2.py:test_round_trips_through_h5mu`).
- A MuData does not survive `pickle` (mudata 0.4.1 with anndata 0.13.4:
  `TypeError: cannot create weak reference to 'NoneType' object`); its
  AnnData modalities do. The asv function benchmarks cache a dict of
  modalities for that reason (`benchmarks/benchmarks/fn.py:FuncGlom`).
- `hmp2()`'s offline tests write synthetic files in HMP2's layout and patch
  `_hmp2._fetch`: the IBDMDB states no licence, so no HMP2 data is committed
  (`tests/datasets/test_hmp2.py`).
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
