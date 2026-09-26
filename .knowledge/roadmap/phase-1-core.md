---
type: Phase
title: Phase 1 - Core, a credible phyloseq replacement (0.1)
description: TreeData conventions in _core; io for phyloseq, BIOM, QIIME 2 and DADA2; pp filter, rarefy, relative, tax_glom; tl alpha, beta, UniFrac, PCoA, NMDS, PERMANOVA; pl bar, richness, ordination, heatmap; R golden tests; Coming-from-R table.
tags: [roadmap, core, io, pp, tl, pl]
status: stable
release: "0.1"
phase_state: in-progress
effort: 6-8 weeks part-time (spec); slices 1A-1D with checkpoints
depends_on: [/roadmap/phase-0-foundation.md]
paths: ["src/biotapy/**", "tests/**", "docs/**", "benchmarks/**"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 0fdbd4d
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: phyloseq-glom
    resource: https://github.com/joey711/phyloseq/blob/master/R/transform_filter-methods.R
    title: phyloseq tax_glom and merge_taxa source
  - id: phyloseq-vignette
    resource: https://github.com/joey711/phyloseq/blob/master/vignettes/phyloseq-analysis.Rmd
    title: phyloseq-analysis vignette
  - id: skbio
    resource: https://scikit.bio/docs/latest/
    title: scikit-bio 0.7.4 documentation
---

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Tasks 1.1-1.5 and slice 1B stage 1
> (1.6-1.9) have full TDD steps. Every later task lists files, interface, tests and done-when; expand
> it with superpowers:writing-plans and get approval before starting (rules.md R1.2a).

**Goal:** 0.1 is a credible phyloseq replacement on TreeData.[^spec]

**Architecture:** all shared mechanics (sparse group sums, slot rules, tree
construction) live in `_core`; every public function is a thin, pure verb over
TreeData that delegates math to SciPy/scikit-bio. Slices build bottom-up:
kernel and first verbs (1A), readers and golden infrastructure (1B), filtering
and diversity (1C), plots, docs and release (1D).

**Tech stack:** anndata · treedata 0.3.x · networkx · numpy · scipy · pandas ·
scikit-bio 0.7.4 · matplotlib · pooch · rdata · biom-format · scikit-learn (NMDS) · pytest/hypothesis · asv.

**Spec:** [plan.md](../../plan.md). Contracts that bind every task:
[function-shape](/contracts/function-shape.md),
[data-model-slots](/contracts/data-model-slots.md),
[module-boundaries](/contracts/module-boundaries.md),
[tree-access](/contracts/tree-access.md),
[r-golden-parity](/contracts/r-golden-parity.md).

# Global constraints
- Python >= 3.12; `treedata>=0.3.1,<0.4`; `scikit-bio>=0.7.4,<0.8`.
- `X` is CSR, samples x features; readers transpose once ([samples-as-rows](/decisions/samples-as-rows.md)).
- Function body <= 30 statements, complexity <= 8 (rules.md R5).
- scikit-bio's AnnData dispatch fails on sparse `X`: wrappers pass dense arrays plus ids (rules.md R6.2).
- UniFrac runs through scikit-bio; the `unifrac` package has no pip wheels ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).

# Dependencies to approve (ask at the start of the task named)
| Task | Group | Package | Reason |
|---|---|---|---|
| 1.1 | runtime | scipy, pandas (anndata already present) | sparse kernels, slots - approved 2026-09-26 |
| 1.1 | dev | pandas-stubs, scipy-stubs | `mypy --strict` cannot type untyped scipy/pandas - approved 2026-09-26 |
| 1.3 | runtime | treedata `>=0.3.1,<0.4`, networkx | container and tree |
| 1.3 | dev | types-networkx | networkx ships no type information - approved 2026-09-26 |
| 1.6 | spike only | rdata (+ xarray), throwaway `uv run --with` env | a runtime dependency only if the 1.6 decision picks the native route - spike env approved 2026-09-26 |
| 1.7c | runtime | biom-format `>=2.1.16` | BIOM 1.0 JSON and 2.1 HDF5; no CPython 3.14 wheels yet, builds from source - approved 2026-09-26 |
| 1.7a | runtime | scikit-bio `>=0.7.4,<0.8` | Newick parsing, diversity, ordination - approved 2026-09-26 |
| 1.11 | runtime | pooch | cached dataset downloads |
| 1.12 | test | pyarrow | read parquet golden files |
| 1.17 | runtime | scikit-learn | non-metric MDS (scikit-bio has none) |
| 1.18 | runtime | matplotlib | `pl` |
| 1.21 | dev | asv | benchmarks |

# Review focus
1. **Non-string or duplicated sample/feature ids** from readers (BIOM ids can be ints) -> readers cast to `str` and fail on duplicates naming them. Tests in 1.7a, 1.7c-1.9.
2. **Tree tips and table features disagree** -> readers keep the intersection and warn with both counts, never error deep inside TreeData. Tests in 1.7a, 1.7c, 1.9.
3. **Same genus name in different lineages** ("uncultured") -> `tax_glom` groups by lineage. Test in 1.5.
4. **All-zero samples** after filtering -> `relative` keeps zeros; `rarefy` drops them; `alpha` returns NaN, never raises. Tests in 1.4, 1.14, 1.15.
5. **Memory on realistic data** (5,000 x 50,000) -> `alpha` densifies in bounded row chunks; `beta` documents its one dense copy. Tests in 1.15-1.16 assert chunking; asv in 1.21 measures it.

---

## Slice 1A - Kernel and first verbs

### Task 1.1: `_core` sparse kernels

**Files:** create `src/biotapy/_core/_matrix.py`, `tests/core/test_matrix.py`;
modify `src/biotapy/_core/__init__.py`, `pyproject.toml` (runtime scipy, pandas; dev
pandas-stubs, scipy-stubs).
**Interfaces (produces):**
- `as_csr(X: object) -> sp.csr_matrix` (widened from `sp.spmatrix | sp.sparray | npt.ArrayLike`
  in commit fc26baa so callers pass `AnnData.X` directly without a cast)
- `sum_by(X: sp.csr_matrix, codes: npt.NDArray[np.intp], n_groups: int) -> sp.csr_matrix`
- `argmax_by(values: npt.NDArray[np.float64], codes: npt.NDArray[np.intp]) -> npt.NDArray[np.intp]`

- [x] **Step 1: Failing tests**
  ```python
  # tests/core/test_matrix.py
  import numpy as np
  import scipy.sparse as sp
  from hypothesis import given
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays

  from biotapy._core import argmax_by, as_csr, sum_by

  X = sp.csr_matrix(np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int64))


  def test_as_csr_converts_dense():
      assert isinstance(as_csr(np.eye(2)), sp.csr_matrix)


  def test_as_csr_does_not_copy_csr():
      assert as_csr(X) is X


  def test_sum_by_sums_columns_per_group():
      np.testing.assert_array_equal(sum_by(X, np.array([0, 1, 0]), 2).toarray(), [[4, 2], [10, 5]])


  def test_sum_by_drops_negative_codes():
      np.testing.assert_array_equal(sum_by(X, np.array([0, -1, 0]), 1).toarray(), [[4], [10]])


  def test_sum_by_keeps_integer_dtype():
      assert sum_by(X, np.array([0, 0, 0]), 1).dtype == np.int64


  def test_argmax_by_one_index_per_group_in_code_order():
      np.testing.assert_array_equal(argmax_by(np.array([1.0, 9.0, 3.0, 7.0]), np.array([1, 0, 1, 0])), [1, 2])


  def test_argmax_by_takes_first_on_ties():
      np.testing.assert_array_equal(argmax_by(np.array([5.0, 5.0]), np.array([0, 0])), [0])


  def test_argmax_by_skips_negative_codes():
      np.testing.assert_array_equal(argmax_by(np.array([9.0, 1.0, 2.0]), np.array([-1, 0, 0])), [2])


  def test_argmax_by_with_no_valid_codes_is_empty():
      assert argmax_by(np.array([1.0]), np.array([-1])).size == 0


  @given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 50)), st.data())
  def test_sum_by_preserves_sample_totals(dense, data):
      codes = np.array(data.draw(st.lists(st.integers(0, 2), min_size=dense.shape[1], max_size=dense.shape[1])))
      out = sum_by(sp.csr_matrix(dense), codes, 3)
      np.testing.assert_array_equal(np.asarray(out.sum(axis=1)).ravel(), dense.sum(axis=1))
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_matrix.py -q` -> `ImportError: cannot import name 'argmax_by'`.
- [x] **Step 3: Implement**
  ```python
  # src/biotapy/_core/_matrix.py
  """Sparse kernels shared by pp, fn and tl."""

  import numpy as np
  import numpy.typing as npt
  import scipy.sparse as sp


  def as_csr(X: sp.spmatrix | sp.sparray | npt.ArrayLike) -> sp.csr_matrix:
      """Return ``X`` as a CSR matrix, without copying one that already is."""
      if isinstance(X, sp.csr_matrix):
          return X
      return sp.csr_matrix(X)


  def sum_by(X: sp.csr_matrix, codes: npt.NDArray[np.intp], n_groups: int) -> sp.csr_matrix:
      """Sum the columns of ``X`` that share a group code; negative codes are dropped."""
      rows = np.flatnonzero(codes >= 0)
      indicator = sp.csr_matrix(
          (np.ones(rows.size, dtype=X.dtype), (rows, codes[rows])),
          shape=(codes.size, n_groups),
      )
      return sp.csr_matrix(X @ indicator)


  def argmax_by(values: npt.NDArray[np.float64], codes: npt.NDArray[np.intp]) -> npt.NDArray[np.intp]:
      """Index of the largest value per group, first on ties, ordered by group code."""
      valid = np.flatnonzero(codes >= 0)
      if valid.size == 0:
          return valid
      order = valid[np.lexsort((-values[valid], codes[valid]))]
      first = np.r_[True, np.diff(codes[order]) != 0]
      return order[first]
  ```
  Export all three from `src/biotapy/_core/__init__.py` (import + `__all__`).
- [x] **Step 4: Run, expect pass** - same command -> 10 passed.
- [x] **Step 5: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(core): add sparse group-sum and group-argmax kernels"`

### Task 1.2: `_core` taxonomy and slot rules

**Files:** create `src/biotapy/_core/_taxonomy.py`, `src/biotapy/_core/_slots.py`,
`tests/core/test_taxonomy.py`, `tests/core/test_slots.py`; modify `_core/__init__.py`.
**Interfaces (produces):**
- `RANKS: tuple[str, ...]` = kingdom .. species
- `split_ranks(adata: AnnData, rank: str) -> tuple[list[str], list[str]]` (up to and including `rank`, below it)
- `XKind = Literal["counts", "relative", "rpk", "cpm", "abundance"]`
- `x_kind(adata) -> XKind`; `require_counts(adata, *, func: str) -> None`
- `add_provenance(adata, step: str, **params: str | int | float | bool | None) -> None`
- `feature_subset(adata, index: npt.NDArray[np.intp]) -> AnnData`

- [x] **Step 1: Failing tests**
  ```python
  # tests/core/test_taxonomy.py
  import anndata as ad
  import numpy as np
  import pandas as pd
  import pytest

  from biotapy._core import split_ranks


  def _adata(columns: list[str]) -> ad.AnnData:
      var = pd.DataFrame({c: ["x"] for c in columns}, index=["f1"])
      return ad.AnnData(X=np.zeros((1, 1)), obs=pd.DataFrame(index=["s1"]), var=var)


  def test_split_ranks_uses_canonical_order():
      assert split_ranks(_adata(["genus", "kingdom", "phylum"]), "phylum") == (["kingdom", "phylum"], ["genus"])


  def test_split_ranks_ignores_non_rank_columns():
      assert split_ranks(_adata(["kingdom", "sequence"]), "kingdom") == (["kingdom"], [])


  def test_split_ranks_names_the_missing_rank():
      with pytest.raises(KeyError, match="genus"):
          split_ranks(_adata(["kingdom"]), "genus")
  ```
  ```python
  # tests/core/test_slots.py
  import json

  import anndata as ad
  import numpy as np
  import pandas as pd
  import pytest
  import scipy.sparse as sp

  from biotapy._core import add_provenance, feature_subset, require_counts, x_kind


  def _adata() -> ad.AnnData:
      adata = ad.AnnData(
          X=sp.csr_matrix(np.arange(6).reshape(2, 3)),
          obs=pd.DataFrame(index=["s1", "s2"]),
          var=pd.DataFrame(index=["f1", "f2", "f3"]),
      )
      adata.layers["relative"] = adata.X.copy()
      adata.obsm["X_pcoa"] = np.zeros((2, 2))
      adata.obsp["braycurtis"] = sp.csr_matrix((2, 2))
      adata.uns["other"] = 1
      return adata


  def test_x_kind_defaults_to_counts():
      assert x_kind(_adata()) == "counts"


  def test_require_counts_rejects_relative():
      adata = _adata()
      adata.uns["biotapy"] = {"x_kind": "relative"}
      with pytest.raises(ValueError, match="pp.rarefy needs raw counts"):
          require_counts(adata, func="pp.rarefy")


  def test_add_provenance_appends_json_entries():
      adata = _adata()
      add_provenance(adata, "pp.a", rank="genus")
      add_provenance(adata, "pp.b")
      entries = [json.loads(e) for e in adata.uns["biotapy"]["provenance"]]
      assert [e["step"] for e in entries] == ["pp.a", "pp.b"]
      assert entries[0]["params"] == {"rank": "genus"}


  def test_feature_subset_drops_derived_slots():
      out = feature_subset(_adata(), np.array([0, 2]))
      assert list(out.var_names) == ["f1", "f3"]
      assert not out.layers and not out.obsm and not out.obsp
      assert set(out.uns) == {"biotapy"}


  def test_feature_subset_leaves_input_alone():
      adata = _adata()
      feature_subset(adata, np.array([0]))
      assert adata.n_vars == 3 and "relative" in adata.layers and "other" in adata.uns
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core -q` -> ImportError for `split_ranks`.
- [x] **Step 3: Implement**
  ```python
  # src/biotapy/_core/_taxonomy.py
  """Canonical taxonomic ranks (contracts/data-model-slots)."""

  from anndata import AnnData

  RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")


  def split_ranks(adata: AnnData, rank: str) -> tuple[list[str], list[str]]:
      """Split the rank columns of ``var`` into up-to-and-including ``rank`` and below it."""
      present = [r for r in RANKS if r in adata.var.columns]
      if rank not in present:
          msg = f"rank={rank!r} is not a taxonomy column; available ranks: {present}"
          raise KeyError(msg)
      cut = present.index(rank) + 1
      return present[:cut], present[cut:]
  ```
  ```python
  # src/biotapy/_core/_slots.py
  """x_kind, provenance and feature-changing subsets (contracts/data-model-slots)."""

  import json
  from importlib.metadata import version
  from typing import Literal, cast

  import numpy as np
  import numpy.typing as npt
  from anndata import AnnData

  XKind = Literal["counts", "relative", "rpk", "cpm", "abundance"]
  ParamValue = str | int | float | bool | None
  DERIVED_SLOTS = ("layers", "obsm", "obsp", "varm", "varp")


  def x_kind(adata: AnnData) -> XKind:
      """What ``X`` holds; a missing key means raw counts."""
      return cast(XKind, adata.uns.get("biotapy", {}).get("x_kind", "counts"))


  def require_counts(adata: AnnData, *, func: str) -> None:
      """Raise unless ``X`` holds raw counts."""
      kind = x_kind(adata)
      if kind != "counts":
          msg = f"{func} needs raw counts in X, but uns['biotapy']['x_kind'] is {kind!r}"
          raise ValueError(msg)


  def add_provenance(adata: AnnData, step: str, **params: ParamValue) -> None:
      """Append one JSON entry to ``uns['biotapy']['provenance']`` (h5ad cannot store a list of dicts)."""
      meta = adata.uns.setdefault("biotapy", {"x_kind": "counts"})
      entry = json.dumps({"step": step, "version": version("biotapy"), "params": params})
      meta["provenance"] = [*meta.get("provenance", []), entry]


  def feature_subset(adata: AnnData, index: npt.NDArray[np.intp]) -> AnnData:
      """Subset features and drop every slot derived from the old feature set."""
      out = adata[:, index].copy()
      for slot in DERIVED_SLOTS:
          mapping = getattr(out, slot)
          for key in list(mapping.keys()):
              del mapping[key]
      out.uns = {"biotapy": out.uns.get("biotapy", {"x_kind": "counts"})}
      return out
  ```
  Export `RANKS`, `split_ranks`, `XKind`, `x_kind`, `require_counts`,
  `add_provenance`, `feature_subset` from `_core/__init__.py`.
- [x] **Step 4: Run, expect pass** -> 8 new tests pass.
- [x] **Step 5: Gate and commit** - `uvx prek run --all-files`; `git commit -am "feat(core): add rank splitting, x_kind and provenance slot rules"`

### Task 1.3: `_core` tree helpers and `datasets.toy()`

**Files:** create `src/biotapy/_core/_tree.py`, `src/biotapy/datasets/__init__.py`,
`src/biotapy/datasets/_toy.py`, `tests/core/test_tree.py`,
`tests/datasets/test_toy.py`, `docs/guide/index.md`, `docs/guide/data_model.md`;
modify `src/biotapy/__init__.py`, `_core/__init__.py`, `docs/index.md`, `docs/api.md`,
`pyproject.toml` (treedata, networkx; dev types-networkx).
**Interfaces (produces):**
- `PHYLO_KEY = "phylo"`; `TreeData` (re-exported type)
- `tree_from_edges(edges: Iterable[tuple[str, str, float]]) -> nx.DiGraph` (edge attribute `length`)
- `get_tree(tdata: TreeData) -> nx.DiGraph`
- `make_treedata(X, *, obs, var, tree: nx.DiGraph | None, x_kind: XKind, source: str) -> TreeData`
- `bt.datasets.toy() -> TreeData`: 6 samples x 8 features. `obs["group"]` A (s1-s3) / B (s4-s6);
  kingdom..genus with `f8` genus missing; phylum totals give archetypes f3, f6, f7.

- [x] **Step 1: Failing tests**
  ```python
  # tests/core/test_tree.py
  import numpy as np
  import pandas as pd
  import pytest

  from biotapy._core import get_tree, make_treedata, tree_from_edges


  def _make(tree):
      return make_treedata(
          np.ones((1, 2)), obs=pd.DataFrame(index=["s1"]), var=pd.DataFrame(index=["a", "b"]),
          tree=tree, x_kind="counts", source="test",
      )


  def test_tree_from_edges_stores_branch_length():
      assert tree_from_edges([("r", "a", 0.5)]).edges["r", "a"]["length"] == 0.5


  def test_make_treedata_adds_no_label_column():
      assert "tree" not in _make(tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)])).var.columns


  def test_get_tree_without_phylogeny_names_the_key():
      with pytest.raises(KeyError, match="phylo"):
          get_tree(_make(None))
  ```
  ```python
  # tests/datasets/test_toy.py
  import treedata as td

  import biotapy as bt
  from biotapy._core import get_tree


  def test_toy_shape_taxonomy_and_tree():
      tdata = bt.datasets.toy()
      assert tdata.shape == (6, 8)
      assert tdata.var["genus"].isna().tolist() == [False] * 7 + [True]
      assert set(get_tree(tdata).successors("n5")) == {"f4", "f5"}


  def test_toy_is_counts_with_one_provenance_entry():
      meta = bt.datasets.toy().uns["biotapy"]
      assert meta["x_kind"] == "counts" and len(meta["provenance"]) == 1


  def test_toy_round_trips_through_h5td(tmp_path):
      path = tmp_path / "toy.h5td"
      bt.datasets.toy().write_h5td(path)
      back = td.read_h5td(path)
      assert (back.X != bt.datasets.toy().X).nnz == 0
      assert back.vart["phylo"].edges["n4", "f1"]["length"] == 0.1
      assert list(back.uns["biotapy"]["provenance"]) == list(bt.datasets.toy().uns["biotapy"]["provenance"])
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_tree.py tests/datasets -q` -> ImportError.
- [x] **Step 3: Implement the tree helpers**
  ```python
  # src/biotapy/_core/_tree.py
  """The only module that imports treedata or networkx (contracts/tree-access)."""

  from collections.abc import Iterable

  import networkx as nx
  import numpy.typing as npt
  import pandas as pd
  import scipy.sparse as sp
  from treedata import TreeData

  from ._matrix import as_csr
  from ._slots import XKind, add_provenance

  PHYLO_KEY = "phylo"


  def tree_from_edges(edges: Iterable[tuple[str, str, float]]) -> nx.DiGraph:
      """Build a rooted tree from ``(parent, child, branch_length)`` triples."""
      tree = nx.DiGraph()
      tree.add_weighted_edges_from(edges, weight="length")
      return tree


  def get_tree(tdata: TreeData) -> nx.DiGraph:
      """Return the phylogeny in ``vart['phylo']``."""
      if PHYLO_KEY not in tdata.vart:
          msg = f"no phylogeny in vart[{PHYLO_KEY!r}]"
          raise KeyError(msg)
      return tdata.vart[PHYLO_KEY]


  def make_treedata(
      X: sp.spmatrix | npt.ArrayLike,
      *,
      obs: pd.DataFrame,
      var: pd.DataFrame,
      tree: nx.DiGraph | None,
      x_kind: XKind,
      source: str,
  ) -> TreeData:
      """Construct a TreeData that follows contracts/data-model-slots."""
      vart = None if tree is None else {PHYLO_KEY: tree}
      tdata = TreeData(X=as_csr(X), obs=obs, var=var, vart=vart, label=None)
      tdata.uns["biotapy"] = {"x_kind": x_kind}
      add_provenance(tdata, source)
      return tdata
  ```
  Export `PHYLO_KEY`, `TreeData`, `tree_from_edges`, `get_tree`, `make_treedata` from `_core/__init__.py`.
- [x] **Step 4: Implement the toy dataset**
  ```python
  # src/biotapy/datasets/_toy.py
  """Tiny in-memory dataset for docstring examples and tests."""

  import numpy as np
  import pandas as pd

  from biotapy._core import TreeData, make_treedata, tree_from_edges

  _COUNTS = np.array(
      [
          [10, 5, 20, 30, 0, 2, 1, 0],
          [8, 7, 25, 22, 3, 0, 0, 1],
          [12, 4, 18, 35, 1, 5, 2, 0],
          [2, 1, 5, 10, 0, 40, 15, 3],
          [0, 2, 3, 12, 2, 38, 20, 5],
          [1, 0, 4, 8, 1, 45, 12, 2],
      ],
      dtype=np.int64,
  )
  _FIRMICUTES = ("Bacteria", "Firmicutes", "Clostridia")
  _BACTEROIDOTA = ("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales")
  _PROTEOBACTERIA = ("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae")
  _TAXONOMY = [
      (*_FIRMICUTES, "Lachnospirales", "Lachnospiraceae", "Blautia"),
      (*_FIRMICUTES, "Lachnospirales", "Lachnospiraceae", "Roseburia"),
      (*_FIRMICUTES, "Oscillospirales", "Ruminococcaceae", "Faecalibacterium"),
      (*_BACTEROIDOTA, "Bacteroidaceae", "Bacteroides"),
      (*_BACTEROIDOTA, "Bacteroidaceae", "Bacteroides"),
      (*_BACTEROIDOTA, "Prevotellaceae", "Prevotella"),
      (*_PROTEOBACTERIA, "Escherichia"),
      (*_PROTEOBACTERIA, None),
  ]
  _EDGES = [
      ("root", "n1", 0.1), ("root", "n2", 0.1), ("root", "n3", 0.2),
      ("n1", "n4", 0.05), ("n4", "f1", 0.1), ("n4", "f2", 0.12), ("n1", "f3", 0.2),
      ("n2", "n5", 0.05), ("n5", "f4", 0.02), ("n5", "f5", 0.03), ("n2", "f6", 0.15),
      ("n3", "f7", 0.1), ("n3", "f8", 0.12),
  ]
  _RANK_COLUMNS = ["kingdom", "phylum", "class", "order", "family", "genus"]


  def toy() -> TreeData:
      """Six gut samples x eight features with taxonomy and a phylogeny.

      Built in memory, so examples and tests never download anything.

      Returns
      -------
      TreeData
          Counts in ``X``; ``obs['group']`` is ``A`` (s1-s3) or ``B`` (s4-s6);
          ``var`` holds kingdom to genus, with ``f8`` unassigned at genus;
          ``vart['phylo']`` carries branch lengths.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/data_model`

      Examples
      --------
      >>> import biotapy as bt
      >>> bt.datasets.toy().shape
      (6, 8)
      """
      obs = pd.DataFrame({"group": pd.Categorical(["A"] * 3 + ["B"] * 3)}, index=[f"s{i}" for i in range(1, 7)])
      var = pd.DataFrame(_TAXONOMY, index=[f"f{i}" for i in range(1, 9)], columns=_RANK_COLUMNS)
      return make_treedata(_COUNTS, obs=obs, var=var, tree=tree_from_edges(_EDGES), x_kind="counts", source="datasets.toy")
  ```
  `src/biotapy/datasets/__init__.py`: `from ._toy import toy` and `__all__ = ["toy"]`.
  `src/biotapy/__init__.py`: add `from . import datasets` and `"datasets"` to `__all__`.
- [x] **Step 5: Docs.** `docs/guide/index.md` (toctree of guide pages) linked
  from `docs/index.md`; `docs/guide/data_model.md`: the slot table and the
  samples-as-rows rule, written for users (from [data-model-slots](/contracts/data-model-slots.md)).
  Add `datasets.toy` to `docs/api.md`.
- [x] **Step 6: Run, expect pass** - `uv run --group test pytest -q` (includes the doctest) and
  `uv run --group doc sphinx-build -W -b html docs docs/_build/html`.
- [x] **Step 7: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(datasets): add in-memory toy TreeData and core tree helpers"`

### Task 1.4: `pp.relative`

**Files:** create `src/biotapy/pp/__init__.py`, `src/biotapy/pp/_transform.py`,
`tests/pp/test_transform.py`, `docs/guide/transforms.md`; modify
`tests/conftest.py`, `src/biotapy/__init__.py`, `docs/api.md`, `docs/guide/index.md`.
**Interfaces:** consumes `as_csr`, `add_provenance`; produces
`bt.pp.relative(adata: AnnData) -> AnnData` adding `layers["relative"]`; test fixture `assert_unchanged`.

- [x] **Step 1: Purity fixture** - append to `tests/conftest.py`:
  ```python
  from collections.abc import Callable

  import pandas as pd
  import pytest
  from anndata import AnnData

  from biotapy._core import as_csr


  def _assert_unchanged(before: AnnData, after: AnnData) -> None:
      assert (as_csr(before.X) != as_csr(after.X)).nnz == 0
      pd.testing.assert_frame_equal(before.obs, after.obs)
      pd.testing.assert_frame_equal(before.var, after.var)
      for slot in ("layers", "obsm", "obsp", "uns"):
          assert set(getattr(before, slot).keys()) == set(getattr(after, slot).keys()), slot


  @pytest.fixture
  def assert_unchanged() -> Callable[[AnnData, AnnData], None]:
      """Fail if a biotapy call mutated its input (rules.md R3.3)."""
      return _assert_unchanged
  ```
- [x] **Step 2: Failing tests**
  ```python
  # tests/pp/test_transform.py
  import json

  import anndata as ad
  import numpy as np
  import pandas as pd
  import scipy.sparse as sp
  from hypothesis import given
  from hypothesis import strategies as st
  from hypothesis.extra.numpy import arrays

  import biotapy as bt


  def _row_sums(M) -> np.ndarray:
      return np.asarray(M.sum(axis=1)).ravel()


  def test_relative_rows_sum_to_one():
      np.testing.assert_allclose(_row_sums(bt.pp.relative(bt.datasets.toy()).layers["relative"]), 1.0)


  def test_relative_keeps_x_and_input(assert_unchanged):
      tdata = bt.datasets.toy()
      before = tdata.copy()
      out = bt.pp.relative(tdata)
      assert_unchanged(before, tdata)
      assert (out.X != tdata.X).nnz == 0


  def test_relative_all_zero_sample_stays_zero():
      tdata = bt.datasets.toy()
      dense = tdata.X.toarray()
      dense[0] = 0
      tdata.X = sp.csr_matrix(dense)
      assert bt.pp.relative(tdata).layers["relative"][0].nnz == 0


  def test_relative_records_provenance():
      entries = bt.pp.relative(bt.datasets.toy()).uns["biotapy"]["provenance"]
      assert json.loads(entries[-1])["step"] == "pp.relative"


  @given(arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 8)), elements=st.integers(0, 1000)))
  def test_relative_rows_sum_to_one_or_zero(dense):
      adata = ad.AnnData(
          X=sp.csr_matrix(dense),
          obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
          var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
      )
      expected = np.where(dense.sum(axis=1) > 0, 1.0, 0.0)
      np.testing.assert_allclose(_row_sums(bt.pp.relative(adata).layers["relative"]), expected)
  ```
- [x] **Step 3: Run, expect failure** - `uv run --group test pytest tests/pp -q` -> `AttributeError: module 'biotapy' has no attribute 'pp'`.
- [x] **Step 4: Implement**
  ```python
  # src/biotapy/pp/_transform.py
  """Per-sample transforms: add one layer, keep everything else."""

  import numpy as np
  import scipy.sparse as sp
  from anndata import AnnData

  from biotapy._core import add_provenance, as_csr


  def relative(adata: AnnData) -> AnnData:
      """Add per-sample relative abundance as ``layers['relative']``.

      Parameters
      ----------
      adata
          Samples x features; ``X`` holds counts or another non-negative abundance.

      Returns
      -------
      AnnData
          A copy of ``adata`` (a TreeData stays a TreeData) with
          ``layers['relative']``; ``X`` is unchanged.

      Notes
      -----
      R equivalent: ``phyloseq::transform_sample_counts``, ``mia::transformAssay``
      Guide: :doc:`/guide/transforms`

      All-zero samples stay all-zero, where phyloseq returns ``NaN``.

      Examples
      --------
      >>> import biotapy as bt
      >>> out = bt.pp.relative(bt.datasets.toy())
      >>> round(float(out.layers["relative"][0].sum()), 6)
      1.0
      """
      X = as_csr(adata.X)
      sums = np.asarray(X.sum(axis=1), dtype=np.float64).ravel()
      scale = np.divide(1.0, sums, out=np.zeros_like(sums), where=sums > 0)
      out = adata.copy()
      out.layers["relative"] = sp.csr_matrix(sp.diags(scale) @ X)
      add_provenance(out, "pp.relative")
      return out
  ```
  `src/biotapy/pp/__init__.py`: `from ._transform import relative`, `__all__ = ["relative"]`.
  `src/biotapy/__init__.py`: add `pp`.
- [x] **Step 5: Docs** - `docs/guide/transforms.md` (what `relative` does, the
  zero-sample difference from phyloseq); add to guide toctree and `docs/api.md`.
- [x] **Step 6: Run, expect pass** - tests, doctest and `sphinx-build -W`.
- [x] **Step 7: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(pp): add relative abundance transform"`

### Task 1.5: `pp.tax_glom`

**Files:** create `src/biotapy/pp/_glom.py`, `tests/pp/test_glom.py`,
`docs/guide/aggregation.md`; modify `src/biotapy/pp/__init__.py`, `docs/api.md`, `docs/guide/index.md`.
**Interfaces:** consumes `split_ranks`, `as_csr`, `sum_by`, `argmax_by`,
`feature_subset`, `add_provenance`; produces
`bt.pp.tax_glom(adata: AnnData, rank: str, *, dropna: bool = True) -> AnnData`.
Semantics follow phyloseq exactly:[^phyloseq-glom] group by the lineage string
joined with `";_;"` (missing -> `"NA"`), archetype = most abundant member (first
on ties), ranks below `rank` set to `NaN`.

- [x] **Step 1: Failing tests**
  ```python
  # tests/pp/test_glom.py
  import numpy as np
  import pytest
  from hypothesis import given
  from hypothesis import strategies as st

  import biotapy as bt
  from biotapy._core import get_tree

  RANKS = ["kingdom", "phylum", "class", "order", "family", "genus"]


  def _leaves(tdata) -> set[str]:
      tree = get_tree(tdata)
      return {n for n in tree.nodes if tree.out_degree(n) == 0}


  def test_phylum_keeps_most_abundant_member_in_original_order():
      out = bt.pp.tax_glom(bt.datasets.toy(), "phylum")
      assert list(out.var_names) == ["f3", "f6", "f7"]
      assert out.var[["class", "order", "family", "genus"]].isna().all().all()
      assert _leaves(out) == {"f3", "f6", "f7"}


  def test_genus_sums_members_and_drops_unassigned():
      out = bt.pp.tax_glom(bt.datasets.toy(), "genus")
      assert list(out.var_names) == ["f1", "f2", "f3", "f4", "f6", "f7"]
      np.testing.assert_array_equal(out[:, "f4"].X.toarray().ravel(), [30, 25, 36, 10, 14, 9])


  def test_dropna_false_keeps_unassigned_as_own_group():
      assert bt.pp.tax_glom(bt.datasets.toy(), "genus", dropna=False).n_vars == 7


  def test_groups_by_lineage_not_label():
      tdata = bt.datasets.toy()
      tdata.var.loc["f7", "genus"] = "Blautia"
      assert bt.pp.tax_glom(tdata, "genus").n_vars == 6


  def test_drops_derived_slots():
      assert "relative" not in bt.pp.tax_glom(bt.pp.relative(bt.datasets.toy()), "phylum").layers


  def test_input_unchanged(assert_unchanged):
      tdata = bt.datasets.toy()
      before = tdata.copy()
      bt.pp.tax_glom(tdata, "genus")
      assert_unchanged(before, tdata)


  def test_unknown_rank_is_named():
      with pytest.raises(KeyError, match="species"):
          bt.pp.tax_glom(bt.datasets.toy(), "species")


  def test_rank_with_no_values_raises():
      tdata = bt.datasets.toy()
      tdata.var["genus"] = None
      with pytest.raises(ValueError, match="genus"):
          bt.pp.tax_glom(tdata, "genus")


  @given(st.sampled_from(RANKS))
  def test_sample_totals_preserved_without_dropna(rank):
      tdata = bt.datasets.toy()
      out = bt.pp.tax_glom(tdata, rank, dropna=False)
      np.testing.assert_array_equal(np.asarray(out.X.sum(axis=1)).ravel(), np.asarray(tdata.X.sum(axis=1)).ravel())
  ```
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/pp/test_glom.py -q` -> `AttributeError: ... 'tax_glom'`.
- [x] **Step 3: Implement**
  ```python
  # src/biotapy/pp/_glom.py
  """Aggregation along the taxonomy."""

  import numpy as np
  import pandas as pd
  from anndata import AnnData

  from biotapy._core import add_provenance, argmax_by, as_csr, feature_subset, split_ranks, sum_by


  def tax_glom(adata: AnnData, rank: str, *, dropna: bool = True) -> AnnData:
      """Aggregate features to a taxonomic rank.

      Features that share a lineage down to ``rank`` are summed into one
      feature, represented by the group's most abundant member.

      Parameters
      ----------
      adata
          Samples x features with taxonomy columns in ``var``.
      rank
          Canonical rank name, for example ``"genus"``.
      dropna
          Drop features with no value at ``rank`` before aggregating.

      Returns
      -------
      AnnData
          Same type as ``adata``; a TreeData keeps the representatives' subtree.
          Ranks below ``rank`` are ``NaN``. ``layers``, ``obsm`` and ``obsp`` are
          dropped because they described the old features.

      Raises
      ------
      KeyError
          ``rank`` is not a taxonomy column.
      ValueError
          No feature has a value at ``rank``.

      Notes
      -----
      R equivalent: ``phyloseq::tax_glom``, ``mia::agglomerateByRank``
      Guide: :doc:`/guide/aggregation`

      Examples
      --------
      >>> import biotapy as bt
      >>> bt.pp.tax_glom(bt.datasets.toy(), "phylum").n_vars
      3
      """
      upto, below = split_ranks(adata, rank)
      # phyloseq keys groups by paste(lineage, collapse=";_;"), NA included
      lineage = adata.var[upto].astype("string").fillna("NA").agg(";_;".join, axis=1)
      if dropna:
          lineage = lineage.where(adata.var[rank].notna())
      codes, groups = pd.factorize(lineage)
      if groups.size == 0:
          msg = f"no feature has a value at rank={rank!r}"
          raise ValueError(msg)
      X = as_csr(adata.X)
      keep = argmax_by(np.asarray(X.sum(axis=0), dtype=np.float64).ravel(), codes)
      order = np.argsort(keep)
      out = feature_subset(adata, keep[order])
      out.X = sum_by(X, codes, groups.size)[:, order]
      out.var[below] = np.nan
      add_provenance(out, "pp.tax_glom", rank=rank, dropna=dropna)
      return out
  ```
  Add `tax_glom` to `pp/__init__.py` imports and `__all__`.
- [x] **Step 4: Docs** - `docs/guide/aggregation.md`: lineage grouping, the
  archetype rule, what happens to the tree and to derived slots; add to toctree and `docs/api.md`.
- [x] **Step 5: Run, expect pass** - tests, doctests, `sphinx-build -W`.
- [x] **Step 6: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(pp): add tax_glom with phyloseq archetype semantics"`

### Checkpoint A
- [x] Review slice 1A against every contract (superpowers:requesting-code-review).
- [x] Write `Module` concepts `.knowledge/modules/core.md` and `.knowledge/modules/pp.md`
  (codebase-map templates), replace the "modules - not yet documented" line in
  `.knowledge/index.md` with `* [modules](modules/index.md) - ...`, create `modules/index.md`, log it.
- [x] Ask the user to review before slice 1B. (2026-09-26: the user told us to proceed to
  slice 1B; recorded as a go-ahead, not a line-by-line review, so no `verified` was added.)

---

## Slice 1B - Readers, datasets, golden infrastructure

Slice 1B is planned in two stages, because the 1.6 spike decides how phyloseq
objects are read:

- **Stage 1 (tasks 1.6-1.9)** has full TDD steps below. It covers the spike,
  the shared `_core` support for readers, and the BIOM, QIIME 2 and DADA2
  readers plus the BIOM writer.
- **Stage 2 (tasks 1.10-1.12)** is expanded with superpowers:writing-plans once
  the 1.6 decision is approved (rules.md R1.2a).

### Stage 1 design
- **One constructor.** `_core.make_treedata` is the only way readers build a
  TreeData. It casts ids to `str` and raises on duplicates, naming them. It
  also aligns the tree with the table: features that are not tips are dropped,
  tips that are not features are pruned, and one `UserWarning` gives both
  counts. So every reader and `datasets.toy()` get the same rules from one place.
- **One taxonomy normalizer.** `_core/_taxonomy.py` sits next to `RANKS`, so
  the rank vocabulary and its aliases have one home and unit tests.
  - `split_lineage` parses lineage strings: Greengenes `k__`, RESCRIPT/SILVA
    `d__`, SILVA `D_0__`, Greengenes2 truncated, or unprefixed.
  - `normalize_ranks` cleans rank columns (DADA2 `Kingdom`, `NA`).
- **Test fixtures are built at test time** in `tests/io/conftest.py`, with
  biom-format and `zipfile`. No binary fixture is committed, and fixtures never
  use biotapy's own writer, so tests never check biotapy against itself.
- **Not in stage 1:**
  - `.rds` input for `read_dada2`. `rdata` needs `xarray`, and it flattens
    character matrices, dropping `dim` and `dimnames`, so `.rds` support waits
    for the 1.6 decision on `rdata`.
  - Validating QIIME 2 `metadata.yaml` types. That would need a YAML runtime
    dependency; readers check the payload file name instead.

### Stage 1 global constraints (in addition to the Phase 1 list)
- `scikit-bio>=0.7.4,<0.8` (added in 1.7a) and `biom-format>=2.1.16` (added in
  1.7c; scikit-bio requires it anyway, and `io` imports it directly).
- biom-format 2.1.17 has no CPython 3.14 wheels (upstream fix merged, not yet
  released: biocore/biom-format#1004). On 3.14 it builds from the sdist
  (`setuptools`, `numpy`, `cython`). The user chose to keep 3.14 (2026-09-26):
  the stage-1 PR's CI must show the 3.14 jobs green on Linux, macOS and Windows
  before merge, and the docs say 3.14 users need a C compiler until wheels ship.
- Newick parsing always passes `convert_underscores=False`. skbio's default
  turns `ASV_1` into `ASV 1`.
- mypy: a library without `py.typed` gets `follow_untyped_imports = true`,
  never `ignore_missing_imports` (Checkpoint A finding F2).
- Readers transpose at most once, keep `X` CSR, set `x_kind="counts"`, and
  record provenance `io.<reader>`.
- Commits stage explicit paths only. Never commit `.claude/` or `.superpowers/`.

### Stage 1 review focus
1. **Integer or duplicated ids.** BIOM JSON can hold unquoted integer ids. Ids
   become `str`, and duplicates raise naming them: `make_treedata` does this, and
   for BIOM input biom-format itself raises `Duplicate sample IDs!` first. Tests
   in 1.7a and 1.7c.
2. **Tree and table disagree.** Keep the intersection and give one warning with
   both counts; if nothing is shared, raise. Tests in 1.7a, 1.7c and 1.9.
3. **Taxonomy dialects** (Greengenes, RESCRIPT, SILVA `D_n__`, Greengenes2,
   unprefixed, DADA2 capitalised columns with `NA`, QIIME `Unassigned`)
   become canonical columns with NaN for missing values. Tests in 1.7b, 1.8
   and 1.9.
4. **Newick quirks.** Underscores and support values such as `0.95` on
   internal nodes: tip names stay intact and internal nodes get unique names.
   Tests in 1.7a.
5. **Missing ranks surviving a round trip.** BIOM HDF5 drops empty taxonomy
   entries on read, so `write_biom` writes prefixed values (`g__`) and
   positions survive. Tests in 1.7d.

### Task 1.6: Spike - phyloseq objects through `rdata` (timebox: 1 day)

**Question:** can `rdata.read_rda` with a `constructor_dict` turn phyloseq's
`GlobalPatterns.RData` into counts, taxonomy, sample data and a tree?
GlobalPatterns is an S4 `phyloseq` object with slots `otu_table`, `tax_table`,
`sam_data`, `phy_tree` and `refseq`.

**Output:** the Decision concept `.knowledge/decisions/phyloseq-import-route.md`.
It chooses between the native `rdata` route and an R export script shipped in
the docs (write BIOM + Newick + TSV).

**Files:**
- Create `.knowledge/decisions/phyloseq-import-route.md`.
- Modify `.knowledge/decisions/index.md` and `.knowledge/log.md`.
- Spike code lives only in the scratchpad and is never committed.

**Known facts** (research 2026-09-26; re-check each in the spike):
- `rdata` 1.1.0 needs `numpy`, `xarray`, `pandas` and `typing_extensions`.
- `read_rda` returns a dict of the objects in the file. S4 objects become a
  `SimpleNamespace` with one attribute per slot plus `class`, and emit a
  "Missing constructor" `UserWarning` unless `constructor_dict` maps the class.
- Numeric matrices with dimnames become `xarray.DataArray`.
- **Character matrices come back flat** (1-D, column-major, `dim`/`dimnames`
  dropped). phyloseq's `taxonomyTable` is a character matrix, so a
  constructor must reshape it from `attrs["dim"]` and `attrs["dimnames"]`.
- ape `phylo` is an S3 list with `edge` (Nedge x 2, 1-based node ids;
  tips `1..Ntip`, internal nodes after), `edge.length`, `tip.label`, `Nnode`
  and an optional `node.label`.
- Files (raw GitHub, `joey711/phyloseq/master/data/`):
  - `GlobalPatterns.RData` (435,652 B)
  - `enterotype.RData` (195,260 B)
  - `esophagus.RData` (1,840 B)

- [x] **Step 1: Fetch.** Download the three files into the scratchpad (never into the repo).
- [x] **Step 2: Raw read.** Run
  `uv run --no-project --with rdata python spike.py` from the scratchpad. It
  uses a throwaway environment, not a project dependency; the user approved it
  with the stage-1 plan. `spike.py` calls `rdata.read_rda(path)` and prints, for
  each object and slot: its type, shape, the `class` attr and any warnings.
- [x] **Step 3: Constructors.** Add `constructor_dict` entries for:
  - `otu_table`: numeric matrix plus `taxa_are_rows`;
  - `taxonomyTable`: reshape the character matrix with `order="F"`;
  - `sample_data`;
  - `phylo`: build `(parent, child, length)` triples, naming tips from
    `tip.label` and internal nodes `n<i>`;
  - `phyloseq`.
  Record lines of code, run time and remaining warnings.
- [x] **Step 4: Check against known shapes.**
  - GlobalPatterns: 26 samples x 19,216 taxa, 7 rank columns, a tree with
    19,216 tips.
  - enterotype: holds relative abundances.
  - esophagus: 3 samples, with a tree.
  Note every mismatch.
- [x] **Step 5: Decide.** Write `decisions/phyloseq-import-route.md`
  (`type: Decision`, `status: draft`) with:
  - Context;
  - Options: native `rdata` route (new runtime deps `rdata` + `xarray`, the
    constructor code, fragility if phyloseq's S4 layout changes), or the R
    export script (no deps, needs an R install);
  - Evidence (numbers from Steps 2-4);
  - Recommendation;
  - Consequences for 1.9 `.rds`, 1.10, 1.11 and 1.12.
  Add its line to `decisions/index.md` and a `log.md` entry. Commit
  `docs(knowledge): record phyloseq import route spike`.
- [ ] **Step 6: Done when** the user approves the decision (status -> `stable`)
  and the scratchpad spike files are deleted.

### Task 1.7a: `_core` support for readers - Newick parsing and id/tree rules in `make_treedata`

**Files:**
- Modify:
  - `src/biotapy/_core/_tree.py`
  - `src/biotapy/_core/__init__.py`
  - `tests/core/test_tree.py`
  - `pyproject.toml`: runtime `scikit-bio>=0.7.4,<0.8`, plus a mypy override
    if needed
  - `.knowledge/contracts/tree-access.md`: Newick parsing now exists; gotchas
    `convert_underscores=False`, dropped internal labels, NaN lengths,
    alignment warning
  - `.knowledge/contracts/data-model-slots.md`: new convention 5, ids are
    unique `str`, enforced by `make_treedata`
  - `.knowledge/decisions/optional-heavy-dependencies.md`: one sentence
  - `.knowledge/log.md`
**Interfaces (produces):**
- `tree_from_newick(text: str) -> nx.DiGraph[str]`: tips keep their names. Internal nodes are named
  `n0, n1, ...` in preorder, skipping any name that is a tip. A missing
  branch length becomes `nan`. Unnamed or repeated tips raise `ValueError`,
  naming them.
- `make_treedata(X, *, obs, var, tree, x_kind, source) -> TreeData`: the
  signature is unchanged. New behaviour:
  - obs/var ids are cast to `str`; duplicates raise `ValueError("duplicate <axis> ids: [...]")`;
  - with a tree, only features that are tips are kept, and tips outside the
    table are pruned (ancestors kept, as TreeData subsetting does);
  - any mismatch gives one `UserWarning` naming both counts, attributed to the
    first frame outside biotapy (`skip_file_prefixes`, Python 3.12+);
  - no shared feature raises `ValueError`;
  - the input frames are never mutated.

- [ ] **Step 1: Dependency.** Add `"scikit-bio>=0.7.4,<0.8"` to `[project] dependencies`, then run
  `uv sync --group dev --group test --group doc`. On Python 3.14 this builds
  biom-format from source (see stage-1 constraints).
- [ ] **Step 2: Failing tests** - append to `tests/core/test_tree.py`:
  ```python
  import math

  from biotapy._core import tree_from_newick


  def _leaves(tree) -> set[str]:
      return {n for n in tree.nodes if tree.out_degree(n) == 0}


  def _frames(samples, features):
      return pd.DataFrame(index=samples), pd.DataFrame(index=features)


  def test_tree_from_newick_keeps_underscores_in_tip_names():
      assert _leaves(tree_from_newick("((ASV_1:0.1,ASV_2:0.2):0.05,ASV_3:0.3);")) == {"ASV_1", "ASV_2", "ASV_3"}


  def test_tree_from_newick_stores_branch_lengths():
      tree = tree_from_newick("((a:0.1,b:0.2):0.05,c:0.3);")
      parent = next(iter(tree.predecessors("a")))
      assert tree.edges[parent, "a"]["length"] == 0.1


  def test_tree_from_newick_gives_internal_nodes_unique_names():
      tree = tree_from_newick("((a:1,b:1)0.95:1,(c:1,d:1)0.95:1);")
      internal = set(tree.nodes) - {"a", "b", "c", "d"}
      assert len(internal) == 3 and "0.95" not in internal


  def test_tree_from_newick_internal_names_skip_tip_names():
      tree = tree_from_newick("((n0:1,n1:1):1,n2:1);")
      assert tree.number_of_nodes() == 5 and _leaves(tree) == {"n0", "n1", "n2"}


  def test_tree_from_newick_missing_length_is_nan():
      tree = tree_from_newick("(a,b:2);")
      assert math.isnan(tree.edges[next(iter(tree.predecessors("a"))), "a"]["length"])


  def test_tree_from_newick_repeated_tip_names_raise():
      with pytest.raises(ValueError, match="'a'"):
          tree_from_newick("(a:1,a:1);")


  def test_make_treedata_casts_ids_to_str():
      obs, var = _frames([1, 2], [10, 20])
      tdata = make_treedata(np.ones((2, 2)), obs=obs, var=var, tree=None, x_kind="counts", source="test")
      assert list(tdata.obs_names) == ["1", "2"] and list(tdata.var_names) == ["10", "20"]


  def test_make_treedata_names_duplicate_ids():
      obs, var = _frames(["s1", "s1"], ["a", "b"])
      with pytest.raises(ValueError, match="duplicate obs ids.*s1"):
          make_treedata(np.ones((2, 2)), obs=obs, var=var, tree=None, x_kind="counts", source="test")


  def test_make_treedata_keeps_shared_features_and_warns_once():
      obs, var = _frames(["s1"], ["a", "b", "c"])
      tree = tree_from_edges([("r", "a", 1.0), ("r", "b", 1.0), ("r", "d", 1.0)])
      with pytest.warns(UserWarning, match=r"1 feature\(s\) not in the tree and 1 tree tip\(s\)") as record:
          tdata = make_treedata(np.array([[1, 2, 3]]), obs=obs, var=var, tree=tree, x_kind="counts", source="test")
      assert len(record) == 1
      assert list(tdata.var_names) == ["a", "b"]
      assert tdata.X.toarray().tolist() == [[1, 2]]
      assert _leaves(get_tree(tdata)) == {"a", "b"}


  def test_make_treedata_without_shared_features_raises():
      obs, var = _frames(["s1"], ["a"])
      with pytest.raises(ValueError, match="no feature"):
          make_treedata(np.ones((1, 1)), obs=obs, var=var, tree=tree_from_edges([("r", "z", 1.0)]), x_kind="counts", source="test")


  def test_make_treedata_leaves_input_frames_alone():
      obs, var = _frames([1], ["a"])
      make_treedata(np.ones((1, 1)), obs=obs, var=var, tree=None, x_kind="counts", source="test")
      assert list(obs.index) == [1]
  ```
- [ ] **Step 3: Run, expect failure** - `uv run --group test pytest tests/core/test_tree.py -q` -> `ImportError: cannot import name 'tree_from_newick'`.
- [ ] **Step 4: Implement.** In `src/biotapy/_core/_tree.py`:
  - add imports `import itertools`, `import math`, `import warnings`,
    `from collections import Counter`, `from pathlib import Path`,
    `import numpy as np` and `from skbio import TreeNode`;
  - add the constant and the functions below;
  - replace `make_treedata`'s body.
  ```python
  # Warnings point at the first frame outside biotapy, however deep the call.
  _PACKAGE_DIR = str(Path(__file__).resolve().parents[1])


  def tree_from_newick(text: str) -> nx.DiGraph[str]:
      """Parse one Newick tree; tips keep their names, internal nodes get unique ones.

      Internal labels (often support values such as ``0.95``) repeat, so they
      cannot name graph nodes and are dropped. A missing branch length is NaN.
      """
      # skbio turns "_" into " " in unquoted names by default; ASV ids need them intact.
      root = TreeNode.read([text], convert_underscores=False)
      tips = [tip.name for tip in root.tips()]
      _require_unique_names(tips)
      names = _node_names(root, set(tips))
      return tree_from_edges(
          (names[id(node.parent)], names[id(node)], math.nan if node.length is None else float(node.length))
          for node in root.preorder(include_self=False)
      )


  def _require_unique_names(tips: list[str | None]) -> None:
      bad = [name for name, count in Counter(tips).items() if name is None or count > 1]
      if bad:
          msg = f"Newick tips need unique names; unnamed or repeated: {bad[:5]}"
          raise ValueError(msg)


  def _node_names(root: TreeNode, tips: set[str]) -> dict[int, str]:
      fresh = (name for name in (f"n{i}" for i in itertools.count()) if name not in tips)
      return {id(node): node.name if node.is_tip() else next(fresh) for node in root.preorder()}


  def make_treedata(
      X: object,
      *,
      obs: pd.DataFrame,
      var: pd.DataFrame,
      tree: nx.DiGraph[str] | None,
      x_kind: XKind,
      source: str,
  ) -> TreeData:
      """Construct a TreeData that follows contracts/data-model-slots.

      Ids become unique strings. With a tree, only features that are its tips
      are kept and tips outside the table are pruned; a mismatch warns once.
      """
      obs, var = _with_str_ids(obs, "obs"), _with_str_ids(var, "var")
      matrix = as_csr(X)
      if tree is not None:
          matrix, var, tree = _align_tree(matrix, var, tree)
      vart = None if tree is None else {PHYLO_KEY: tree}
      tdata = TreeData(X=matrix, obs=obs, var=var, vart=vart, label=None)
      tdata.uns["biotapy"] = {"x_kind": x_kind}
      add_provenance(tdata, source)
      return tdata


  def _with_str_ids(frame: pd.DataFrame, axis: str) -> pd.DataFrame:
      out = frame.copy()
      out.index = out.index.astype(str)
      duplicated = out.index[out.index.duplicated()].unique().tolist()
      if duplicated:
          msg = f"duplicate {axis} ids: {duplicated[:5]}"
          raise ValueError(msg)
      return out


  def _align_tree(
      X: sp.csr_matrix, var: pd.DataFrame, tree: nx.DiGraph[str]
  ) -> tuple[sp.csr_matrix, pd.DataFrame, nx.DiGraph[str]]:
      tips = {node for node in tree.nodes if tree.out_degree(node) == 0}
      shared = var.index.isin(tips)
      n_extra = len(tips - set(var.index))
      if shared.all() and n_extra == 0:
          return X, var, tree
      if not shared.any():
          msg = "no feature of the table is a tip of the tree"
          raise ValueError(msg)
      msg = (
          f"tree and table disagree: {int((~shared).sum())} feature(s) not in the tree and "
          f"{n_extra} tree tip(s) not in the table; keeping the {int(shared.sum())} shared features"
      )
      warnings.warn(msg, UserWarning, skip_file_prefixes=(_PACKAGE_DIR,))
      kept = var.index[shared]
      keep_nodes = set(kept).union(*(nx.ancestors(tree, tip) for tip in kept))
      return X[:, np.flatnonzero(shared)], var.loc[kept], tree.subgraph(keep_nodes).copy()
  ```
  Export `tree_from_newick` from `_core/__init__.py`. The `X` annotation widens to
  `object`, matching `as_csr` since fc26baa; `numpy.typing` may become unused, so
  drop it if ruff says so.
  - **mypy:** if it reports `import-untyped` for `skbio`, add
    `{ module = "skbio", follow_untyped_imports = true, implicit_reexport = true }` and the
    same for `"skbio.*"` next to the treedata overrides, with a one-line comment.
  - **Annotation-only fixes** follow ruling P1.
- [ ] **Step 5: Run, expect pass** - `uv run --group test pytest tests/core -q`. The existing `toy()`
  tests still pass with no warning, because toy's tree tips equal its features.
  Also record `uv run python -X importtime -c "import biotapy" 2>&1 | tail -1` in
  the report: rules.md R10.1 wants a measurement before anyone makes skbio
  import lazily.
- [ ] **Step 6: Knowledge and gate.**
  - Update the two contracts, the dependency decision and the log as listed
    under Files.
  - Run `uvx prek run --all-files` and `uv run --group test pytest`.
  - Commit `feat(core): parse Newick and align trees with tables in make_treedata`.

### Task 1.7b: `_core` taxonomy normalization

**Files:** modify `src/biotapy/_core/_taxonomy.py`, `src/biotapy/_core/__init__.py`,
`tests/core/test_taxonomy.py`.
**Interfaces (produces):**
- `normalize_ranks(frame: pd.DataFrame) -> pd.DataFrame`:
  - rank-like columns are renamed to canonical lowercase (case-insensitive;
    `domain` -> `kingdom`);
  - in rank columns, `k__`/`D_0__` prefixes are stripped, and `""`,
    whitespace, `"NA"` and bare prefixes become NaN;
  - other columns are untouched, and a new frame is returned.
- `split_lineage(lineage: pd.Series) -> pd.DataFrame`:
  - splits `;`-separated lineages;
  - a prefixed part (`k__`, `d__`, `p__` ... `s__`, `D_<n>__`) goes to its
    rank; an unprefixed part goes by position;
  - columns run from kingdom down to the deepest rank seen;
  - the index is kept; a NaN lineage gives an all-NaN row.

- [ ] **Step 1: Failing tests** - append to `tests/core/test_taxonomy.py`:
  ```python
  from biotapy._core import normalize_ranks, split_lineage

  LINEAGES = pd.Series(
      [
          "k__Bacteria; p__Firmicutes; c__; o__; f__; g__; s__",
          "d__Bacteria; p__Proteobacteria",
          "D_0__Archaea;D_1__Euryarchaeota",
          "Bacteria;Firmicutes",
          "p__Firmicutes_A;c__Clostridia_258483",
          np.nan,
          "Unassigned",
      ],
      index=[f"f{i}" for i in range(7)],
  )


  def test_split_lineage_reads_greengenes_prefixes():
      out = split_lineage(LINEAGES)
      assert list(out.columns) == ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
      assert out.loc["f0", ["kingdom", "phylum"]].tolist() == ["Bacteria", "Firmicutes"]
      assert out.loc["f0", "class":].isna().all()


  @pytest.mark.parametrize(
      ("feature", "expected"),
      [("f1", ["Bacteria", "Proteobacteria"]), ("f2", ["Archaea", "Euryarchaeota"]), ("f3", ["Bacteria", "Firmicutes"])],
  )
  def test_split_lineage_reads_silva_and_unprefixed(feature, expected):
      assert split_lineage(LINEAGES).loc[feature, ["kingdom", "phylum"]].tolist() == expected


  def test_split_lineage_places_truncated_lineage_by_prefix():
      row = split_lineage(LINEAGES).loc["f4"]
      assert pd.isna(row["kingdom"]) and row["phylum"] == "Firmicutes_A" and row["class"] == "Clostridia_258483"


  def test_split_lineage_missing_lineage_is_all_nan():
      assert split_lineage(LINEAGES).loc["f5"].isna().all()


  def test_split_lineage_keeps_unassigned_as_kingdom():
      assert split_lineage(LINEAGES).loc["f6", "kingdom"] == "Unassigned"


  def test_split_lineage_columns_stop_at_deepest_rank():
      assert list(split_lineage(pd.Series(["k__A; p__B"], index=["x"])).columns) == ["kingdom", "phylum"]


  def test_normalize_ranks_canonicalizes_names_and_missing_values():
      frame = pd.DataFrame(
          {"Domain": ["Bacteria", "NA", " "], "Genus": ["g__", None, "g__Blautia"], "sequence": ["AC", "GT", "TT"]},
          index=["a", "b", "c"],
      )
      out = normalize_ranks(frame)
      assert list(out.columns) == ["kingdom", "genus", "sequence"]
      assert out.loc["a", "kingdom"] == "Bacteria" and out["kingdom"].iloc[1:].isna().all()
      assert out["genus"].iloc[:2].isna().all() and out.loc["c", "genus"] == "Blautia"
      assert out["sequence"].tolist() == ["AC", "GT", "TT"]


  def test_normalize_ranks_leaves_input_alone():
      frame = pd.DataFrame({"Genus": ["g__Blautia"]})
      normalize_ranks(frame)
      assert list(frame.columns) == ["Genus"] and frame.iloc[0, 0] == "g__Blautia"
  ```
- [ ] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_taxonomy.py -q` -> ImportError.
- [ ] **Step 3: Implement** in `src/biotapy/_core/_taxonomy.py`. Add imports `re`,
  `numpy as np` and `pandas as pd`; change the module docstring to
  "Canonical taxonomic ranks and their normalization (contracts/data-model-slots).";
  then add:
  ```python
  RANK_ALIASES = {"domain": "kingdom"}
  # Greengenes/RESCRIPT "k__"-style or SILVA "D_0__"-style rank prefixes.
  _PREFIX = re.compile(r"^(?:(?P<letter>[kdpcofgs])__|D_(?P<level>\d+)__)")
  _PREFIX_RANK = {
      "k": "kingdom", "d": "kingdom", "p": "phylum", "c": "class",
      "o": "order", "f": "family", "g": "genus", "s": "species",
  }
  _MISSING = ("", "NA")


  def _canonical(column: object) -> object:
      name = RANK_ALIASES.get(str(column).strip().lower(), str(column).strip().lower())
      return name if name in RANKS else column


  def normalize_ranks(frame: pd.DataFrame) -> pd.DataFrame:
      """Canonical lowercase rank columns with missing values as NaN (contracts/data-model-slots)."""
      out = frame.rename(columns=_canonical)
      for rank in [column for column in out.columns if column in RANKS]:
          values = out[rank].astype("string").str.strip().str.replace(_PREFIX, "", regex=True)
          out[rank] = values.astype(object).mask(values.isna() | values.isin(_MISSING), np.nan)
      return out


  def _rank_and_value(part: str, position: int) -> tuple[str | None, str]:
      match = _PREFIX.match(part)
      if match is None:
          return (RANKS[position] if position < len(RANKS) else None), part
      if match["letter"] is not None:
          return _PREFIX_RANK[match["letter"]], part[match.end() :]
      level = int(match["level"])
      return (RANKS[level] if level < len(RANKS) else None), part[match.end() :]


  def _parse_lineage(text: str) -> dict[str, str]:
      ranks: dict[str, str] = {}
      for position, part in enumerate(text.split(";")):
          rank, value = _rank_and_value(part.strip(), position)
          if rank is not None:
              ranks[rank] = value
      return ranks


  def split_lineage(lineage: pd.Series) -> pd.DataFrame:
      """Split ``;``-separated lineages (``k__Bacteria; p__Firmicutes``) into rank columns."""
      frame = pd.DataFrame([_parse_lineage(t) if isinstance(t, str) else {} for t in lineage], index=lineage.index)
      depth = max((RANKS.index(column) + 1 for column in frame.columns), default=0)
      return normalize_ranks(frame.reindex(columns=list(RANKS[:depth])))
  ```
  Export `normalize_ranks` and `split_lineage` from `_core/__init__.py`. The
  controller prototyped this code on pandas 3.0.6 on 2026-09-26; mypy fixes
  follow ruling P1.
- [ ] **Step 4: Run, expect pass** -> 10 new tests pass.
- [ ] **Step 5: Gate and commit** - `uvx prek run --all-files`; commit `feat(core): normalize taxonomy ranks and split lineage strings`.

### Task 1.7c: `io.read_biom`

**Files:**
- Create:
  - `src/biotapy/io/__init__.py`
  - `src/biotapy/io/_biom.py`
  - `tests/io/conftest.py`
  - `tests/io/test_biom.py`
  - `docs/guide/reading_data.md`
- Modify:
  - `src/biotapy/__init__.py` (add `io`)
  - `pyproject.toml` (runtime `biom-format>=2.1.16`, plus a mypy override if needed)
  - `docs/api.md`, `docs/guide/index.md`
  - `.knowledge/decisions/optional-heavy-dependencies.md`, `.knowledge/log.md`
**Interfaces:**
- **Consumes:** `make_treedata`, `split_lineage`, `tree_from_newick`, `as_csr`.
- **Produces:**
  - `bt.io.read_biom(path: str | Path, *, tree: str | Path | None = None) -> TreeData`;
  - the private `_biom_parts(table: biom.Table) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]`
    (X as samples x features, obs, var), reused by 1.8.

- [ ] **Step 1: Dependency.** Add `"biom-format>=2.1.16"` to `[project] dependencies`, then run `uv sync`.
- [ ] **Step 2: Fixtures** - `tests/io/conftest.py`:
  ```python
  import json

  import biom
  import numpy as np
  import pytest
  import scipy.sparse as sp
  from biom.util import biom_open

  NEWICK = "((OTU_1:0.1,OTU_2:0.2):0.05,(OTU_3:0.3,OTU_4:0.4):0.1);"
  TAXONOMY = [
      ["k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__", "f__", "g__", "s__"],
      ["k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Lactobacillales", "f__", "g__", "s__"],
      ["k__Bacteria", "p__Bacteroidetes", "c__", "o__", "f__", "g__", "s__"],
      ["k__Archaea", "p__Euryarchaeota", "c__", "o__", "f__", "g__", "s__"],
  ]


  @pytest.fixture
  def biom_table() -> biom.Table:
      """4 features x 3 samples, laid out as BIOM stores them."""
      counts = sp.csr_matrix(np.array([[5, 0, 3], [1, 2, 0], [0, 4, 6], [7, 0, 0]]))
      return biom.Table(
          counts,
          ["OTU_1", "OTU_2", "OTU_3", "OTU_4"],
          ["S1", "S2", "S3"],
          observation_metadata=[{"taxonomy": t} for t in TAXONOMY],
          sample_metadata=[{"group": g} for g in ["A", "A", "B"]],
      )


  @pytest.fixture
  def biom_hdf5(tmp_path, biom_table):
      path = tmp_path / "table.biom"
      with biom_open(str(path), "w") as handle:
          biom_table.to_hdf5(handle, "biotapy tests")
      return path


  @pytest.fixture
  def biom_json(tmp_path, biom_table):
      path = tmp_path / "table.json.biom"
      path.write_text(biom_table.to_json("biotapy tests"))
      return path


  @pytest.fixture
  def newick(tmp_path):
      path = tmp_path / "tree.nwk"
      path.write_text(NEWICK)
      return path


  @pytest.fixture
  def biom_json_ids(tmp_path):
      """Write a minimal BIOM 1.0 JSON table with ids exactly as given (ints allowed)."""

      def write(sample_ids, observation_ids):
          path = tmp_path / "ids.biom"
          document = {
              "id": None, "format": "Biological Observation Matrix 1.0.0",
              "format_url": "http://biom-format.org", "type": "OTU table",
              "generated_by": "biotapy tests", "date": "2026-09-26T00:00:00",
              "rows": [{"id": i, "metadata": None} for i in observation_ids],
              "columns": [{"id": i, "metadata": None} for i in sample_ids],
              "matrix_type": "sparse", "matrix_element_type": "int",
              "shape": [len(observation_ids), len(sample_ids)], "data": [[0, 0, 5], [1, 1, 3]],
          }
          path.write_text(json.dumps(document))
          return path

      return write
  ```
- [ ] **Step 3: Failing tests** - `tests/io/test_biom.py`:
  ```python
  import biom
  import numpy as np
  import pandas as pd
  import pytest
  from biom.exception import TableException
  from biom.util import biom_open

  import biotapy as bt


  def test_read_biom_hdf5_puts_samples_in_rows(biom_hdf5, biom_table):
      tdata = bt.io.read_biom(biom_hdf5)
      assert tdata.shape == (3, 4)
      assert list(tdata.obs_names) == ["S1", "S2", "S3"]
      assert list(tdata.var_names) == ["OTU_1", "OTU_2", "OTU_3", "OTU_4"]
      np.testing.assert_array_equal(tdata.X.toarray(), biom_table.matrix_data.T.toarray())


  def test_read_biom_json_matches_hdf5(biom_json, biom_hdf5):
      from_json, from_hdf5 = bt.io.read_biom(biom_json), bt.io.read_biom(biom_hdf5)
      assert (from_json.X != from_hdf5.X).nnz == 0
      pd.testing.assert_frame_equal(from_json.var, from_hdf5.var)


  def test_read_biom_taxonomy_becomes_rank_columns(biom_hdf5):
      var = bt.io.read_biom(biom_hdf5).var
      assert list(var.columns) == ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
      assert var.loc["OTU_2", "order"] == "Lactobacillales"
      assert var.loc["OTU_1", ["order", "family", "genus", "species"]].isna().all()


  def test_read_biom_sample_metadata_goes_to_obs(biom_hdf5):
      assert bt.io.read_biom(biom_hdf5).obs["group"].tolist() == ["A", "A", "B"]


  def test_read_biom_integer_ids_become_strings(biom_json_ids):
      tdata = bt.io.read_biom(biom_json_ids([1, 2], [10, 20]))
      assert list(tdata.obs_names) == ["1", "2"] and list(tdata.var_names) == ["10", "20"]


  def test_read_biom_duplicate_ids_are_rejected_by_biom(biom_json_ids):
      # biom-format validates ids itself (biom/err.py SAMPDUP, default state "raise").
      with pytest.raises(TableException, match="Duplicate sample IDs"):
          bt.io.read_biom(biom_json_ids(["S1", "S1"], ["a", "b"]))


  def test_read_biom_attaches_tree(biom_hdf5, newick):
      tree = bt.io.read_biom(biom_hdf5, tree=newick).vart["phylo"]
      assert {n for n in tree.nodes if tree.out_degree(n) == 0} == {"OTU_1", "OTU_2", "OTU_3", "OTU_4"}


  def test_read_biom_tree_mismatch_keeps_shared_features(biom_hdf5, tmp_path):
      path = tmp_path / "partial.nwk"
      path.write_text("((OTU_1:1,OTU_2:1):1,OTU_9:1);")
      with pytest.warns(UserWarning, match="2 feature"):
          tdata = bt.io.read_biom(biom_hdf5, tree=path)
      assert list(tdata.var_names) == ["OTU_1", "OTU_2"]


  def test_read_biom_without_taxonomy_has_no_rank_columns(tmp_path, biom_table):
      path = tmp_path / "bare.biom"
      bare = biom.Table(biom_table.matrix_data, biom_table.ids("observation"), biom_table.ids())
      with biom_open(str(path), "w") as handle:
          bare.to_hdf5(handle, "biotapy tests")
      assert not set(bt.io.read_biom(path).var.columns) & {"kingdom", "phylum"}


  def test_read_biom_records_counts_and_provenance(biom_hdf5):
      meta = bt.io.read_biom(biom_hdf5).uns["biotapy"]
      assert meta["x_kind"] == "counts" and '"io.read_biom"' in meta["provenance"][-1]
  ```
- [ ] **Step 4: Run, expect failure** - `uv run --group test pytest tests/io -q` -> `AttributeError: module 'biotapy' has no attribute 'io'`.
- [ ] **Step 5: Implement** `src/biotapy/io/_biom.py`:
  ```python
  """BIOM tables (JSON 1.0 and HDF5 2.1) through biom-format."""

  from collections.abc import Iterable, Mapping, Sequence
  from pathlib import Path

  import biom
  import numpy as np
  import pandas as pd
  import scipy.sparse as sp

  from biotapy._core import TreeData, as_csr, make_treedata, split_lineage, tree_from_newick


  def read_biom(path: str | Path, *, tree: str | Path | None = None) -> TreeData:
      """Read a BIOM table (JSON 1.0 or HDF5 2.1) with samples as rows.

      Parameters
      ----------
      path
          BIOM file; JSON or HDF5 is detected from the content.
      tree
          Newick file whose tips are the table's observation ids.

      Returns
      -------
      TreeData
          Counts in ``X``; observation ``taxonomy`` metadata as rank columns in
          ``var``; sample metadata in ``obs``; the tree in ``vart['phylo']``.

      Raises
      ------
      biom.exception.TableException
          The file repeats a sample or observation id (raised by biom-format).

      Warns
      -----
      UserWarning
          Tree tips and table features differ; only shared features are kept.

      Notes
      -----
      R equivalent: ``phyloseq::import_biom``
      Guide: :doc:`/guide/reading_data`

      ``X`` is read as counts; BIOM does not record whether it holds counts.

      Examples
      --------
      >>> import tempfile
      >>> from pathlib import Path
      >>> import biom
      >>> from biom.util import biom_open
      >>> import biotapy as bt
      >>> toy = bt.datasets.toy()
      >>> path = Path(tempfile.mkdtemp()) / "toy.biom"
      >>> table = biom.Table(toy.X.T, list(toy.var_names), list(toy.obs_names))
      >>> with biom_open(str(path), "w") as handle:
      ...     table.to_hdf5(handle, "example")
      >>> bt.io.read_biom(path).shape
      (6, 8)
      """
      X, obs, var = _biom_parts(biom.load_table(str(path)))
      phylo = None if tree is None else tree_from_newick(Path(tree).read_text())
      return make_treedata(X, obs=obs, var=var, tree=phylo, x_kind="counts", source="io.read_biom")


  def _biom_parts(table: biom.Table) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]:
      """Counts as samples x features, plus obs and var frames, from a biom Table."""
      samples = [str(i) for i in table.ids(axis="sample")]
      features = [str(i) for i in table.ids(axis="observation")]
      # BIOM stores features x samples; transpose once so samples are rows (R6.1).
      X = as_csr(table.matrix_data.T)
      obs = _metadata_frame(table.metadata(axis="sample"), samples)
      var = _taxonomy_frame(table.metadata(axis="observation"), features)
      return X, obs, var


  def _metadata_frame(metadata: Sequence[Mapping[str, object]] | None, ids: list[str]) -> pd.DataFrame:
      if metadata is None:
          return pd.DataFrame(index=ids)
      return pd.DataFrame([dict(entry) for entry in metadata], index=ids)


  def _taxonomy_frame(metadata: Sequence[Mapping[str, object]] | None, ids: list[str]) -> pd.DataFrame:
      if metadata is None:
          return pd.DataFrame(index=ids)
      return split_lineage(pd.Series([_lineage(entry) for entry in metadata], index=ids, dtype=object))


  def _lineage(entry: Mapping[str, object]) -> str | float:
      # HDF5 gives a list of ranks; JSON gives whatever the producer wrote (list or string).
      value = entry.get("taxonomy") or entry.get("Taxonomy")
      if isinstance(value, str):
          return value
      if isinstance(value, Iterable):
          return "; ".join(str(part) for part in value)
      return np.nan
  ```
  - `src/biotapy/io/__init__.py`: `from ._biom import read_biom` and `__all__ = ["read_biom"]`.
  - `src/biotapy/__init__.py`: add `io`.
  - **mypy:** if `biom` is untyped, add `follow_untyped_imports` overrides for
    `biom` and `biom.*`, as in 1.7a.
- [ ] **Step 6: Docs.**
  - `docs/guide/reading_data.md`, titled "Reading and writing data", with a
    "BIOM" section covering: samples become rows (one transpose); taxonomy
    dialects become rank columns; ids are strings; how a tree/table mismatch is
    handled; duplicate ids are rejected by biom-format.
  - Add it to the guide toctree.
  - `docs/api.md` gets an "Input and output" autosummary block (`.. module:: biotapy.io`) listing `io.read_biom`.
- [ ] **Step 7: Run, expect pass** - the tests, the doctest and `sphinx-build -W`.
- [ ] **Step 8: Knowledge and gate.**
  - Add the dependency sentence to the dependency decision, and a log line.
  - Run `uvx prek run --all-files` and `uv run --group test pytest`.
  - Commit `feat(io): read BIOM tables`.

### Task 1.7d: `io.write_biom`
Added 2026-09-26 at the user's request. Renumbered from 1.7b when stage 1 was
expanded.

**Files:** modify `src/biotapy/io/_biom.py`, `src/biotapy/io/__init__.py`,
`tests/io/test_biom.py`, `docs/guide/reading_data.md`, `docs/api.md`.
**Interfaces:**
- **Consumes:** `read_biom` (for the round trip), `RANKS`, `as_csr`.
- **Produces:** `bt.io.write_biom(adata: AnnData, path: str | Path, *, fmt: Literal["hdf5", "json"] = "hdf5") -> None`.

- [ ] **Step 1: Failing tests** - append to `tests/io/test_biom.py`:
  ```python
  import anndata as ad
  import scipy.sparse as sp

  RANK_COLUMNS = ["kingdom", "phylum", "class", "order", "family", "genus"]


  @pytest.mark.parametrize("fmt", ["hdf5", "json"])
  def test_write_biom_round_trips_toy(tmp_path, fmt):
      toy = bt.datasets.toy()
      path = tmp_path / "toy.biom"
      bt.io.write_biom(toy, path, fmt=fmt)
      back = bt.io.read_biom(path)
      assert list(back.obs_names) == list(toy.obs_names) and list(back.var_names) == list(toy.var_names)
      np.testing.assert_array_equal(back.X.toarray(), toy.X.toarray())
      pd.testing.assert_frame_equal(
          back.var[RANK_COLUMNS].fillna("-"), toy.var[RANK_COLUMNS].fillna("-"), check_dtype=False
      )
      assert back.obs["group"].tolist() == toy.obs["group"].astype(str).tolist()


  def test_write_biom_stores_features_by_samples(tmp_path):
      path = tmp_path / "toy.biom"
      bt.io.write_biom(bt.datasets.toy(), path)
      table = biom.load_table(str(path))
      assert table.shape == (8, 6) and list(table.ids("observation"))[:2] == ["f1", "f2"]


  def test_write_biom_prefixes_ranks_in_canonical_order(tmp_path):
      path = tmp_path / "toy.biom"
      bt.io.write_biom(bt.datasets.toy(), path)
      table = biom.load_table(str(path))
      assert list(table.metadata("f1", "observation")["taxonomy"]) == [
          "k__Bacteria", "p__Firmicutes", "c__Clostridia", "o__Lachnospirales", "f__Lachnospiraceae", "g__Blautia",
      ]
      assert list(table.metadata("f8", "observation")["taxonomy"])[-1] == "g__"


  def test_write_biom_without_taxonomy_writes_no_metadata(tmp_path):
      adata = ad.AnnData(
          X=sp.csr_matrix(np.eye(2)), obs=pd.DataFrame(index=["s1", "s2"]), var=pd.DataFrame(index=["a", "b"])
      )
      path = tmp_path / "bare.biom"
      bt.io.write_biom(adata, path)
      assert biom.load_table(str(path)).metadata(axis="observation") is None


  def test_write_biom_leaves_input_alone(tmp_path, assert_unchanged):
      toy = bt.datasets.toy()
      before = toy.copy()
      bt.io.write_biom(toy, tmp_path / "toy.biom")
      assert_unchanged(before, toy)


  def test_write_biom_does_not_write_the_tree(tmp_path):
      path = tmp_path / "toy.biom"
      bt.io.write_biom(bt.datasets.toy(), path)
      assert "phylo" not in bt.io.read_biom(path).vart
  ```
- [ ] **Step 2: Run, expect failure** -> `AttributeError: ... 'write_biom'`.
- [ ] **Step 3: Implement.** Append to `src/biotapy/io/_biom.py`. Add imports
  `from importlib.metadata import version`, `from typing import Literal`,
  `from anndata import AnnData`, `from biom.util import biom_open`, and `RANKS`
  from `biotapy._core`.
  ```python
  def write_biom(adata: AnnData, path: str | Path, *, fmt: Literal["hdf5", "json"] = "hdf5") -> None:
      """Write ``X`` with taxonomy and sample metadata as a BIOM table.

      Parameters
      ----------
      adata
          Samples x features; rank columns in ``var`` become observation
          ``taxonomy`` metadata, ``obs`` columns become sample metadata.
      path
          Output file.
      fmt
          ``"hdf5"`` (BIOM 2.1) or ``"json"`` (BIOM 1.0).

      Notes
      -----
      R equivalent: ``biomformat::write_biom``
      Guide: :doc:`/guide/reading_data`

      BIOM has no slot for a tree, layers or embeddings: a TreeData's tree and
      everything outside ``X``, rank columns and ``obs`` are not written.
      Sample metadata is written as text. Missing ranks are written as bare
      prefixes (``g__``) so every rank keeps its place.

      Examples
      --------
      >>> import tempfile
      >>> from pathlib import Path
      >>> import biotapy as bt
      >>> path = Path(tempfile.mkdtemp()) / "toy.biom"
      >>> bt.io.write_biom(bt.datasets.toy(), path)
      >>> bt.io.read_biom(path).shape
      (6, 8)
      """
      table = biom.Table(
          # biotapy keeps samples as rows; BIOM stores features x samples.
          as_csr(adata.X).T,
          [str(i) for i in adata.var_names],
          [str(i) for i in adata.obs_names],
          observation_metadata=_taxonomy_metadata(adata.var),
          sample_metadata=_sample_metadata(adata.obs),
      )
      generated_by = f"biotapy {version('biotapy')}"
      if fmt == "json":
          Path(path).write_text(table.to_json(generated_by))
          return
      with biom_open(str(path), "w") as handle:
          table.to_hdf5(handle, generated_by)


  def _taxonomy_metadata(var: pd.DataFrame) -> list[dict[str, list[str]]] | None:
      ranks = [rank for rank in RANKS if rank in var.columns]
      if not ranks:
          return None
      # Prefixed values keep their rank: BIOM HDF5 drops empty list entries on read.
      values = var[ranks].astype("string").fillna("")
      return [
          {"taxonomy": [f"{rank[0]}__{value}" for rank, value in zip(ranks, row, strict=True)]}
          for row in values.itertuples(index=False)
      ]


  def _sample_metadata(obs: pd.DataFrame) -> list[dict[str, str]] | None:
      if obs.columns.empty:
          return None
      values = obs.astype("string").fillna("")
      return [dict(zip(map(str, values.columns), row, strict=True)) for row in values.itertuples(index=False)]
  ```
  - Export `write_biom`.
  - The first letters of `RANKS` (`k p c o f g s`) are unique; that is what makes
    `rank[0]` a safe prefix.
  - If biom-format's HDF5 writer rejects the text sample metadata, confirm the
    cause in `biom/table.py` `general_formatter` before changing anything (R2.2).
- [ ] **Step 4: Docs.**
  - Add a "Writing BIOM" subsection to `docs/guide/reading_data.md`: what is not
    written, and why ranks carry prefixes.
  - Add `io.write_biom` to `docs/api.md`.
- [ ] **Step 5: Run, expect pass**; gate; commit `feat(io): write BIOM tables`.

### Task 1.8: `io.read_qiime2`

**Files:**
- Create `src/biotapy/io/_qiime2.py` and `tests/io/test_qiime2.py`.
- Modify `src/biotapy/io/__init__.py`, `tests/io/conftest.py`,
  `docs/guide/reading_data.md` and `docs/api.md`.

**Interfaces:**
- **Consumes:** `_biom_parts` (from `io/_biom.py`), `make_treedata`,
  `split_lineage`, `tree_from_newick`.
- **Produces:** `bt.io.read_qiime2(table: str | Path, *, taxonomy: str | Path | None = None, tree: str | Path | None = None, metadata: str | Path | None = None) -> TreeData`.

**Formats** (q2-types and qiime2/rachis sources, 2026-09-26):
- **Artifact layout:** a `.qza` is a zip with one top-level `<uuid>/`
  directory. Payloads:
  - `<uuid>/data/feature-table.biom`: BIOM 2.1 HDF5, `FeatureTable[Frequency]`;
  - `<uuid>/data/taxonomy.tsv`: header starts `Feature ID\tTaxon`, optional
    extra columns such as `Confidence`;
  - `<uuid>/data/tree.nwk`: `Phylogeny[Rooted|Unrooted]`.
- **Metadata TSV:**
  - ID headers `id, sampleid, sample id, sample-id, featureid, feature id,
    feature-id` match case-insensitively; `#SampleID, #Sample ID, #OTUID,
    #OTU ID, sample_name` match exactly.
  - Leading `#` comments and blank rows are skipped.
  - An optional `#q2:types` row gives `categorical|numeric` (case-insensitive).
  - Encoding is `utf-8-sig`.

- [ ] **Step 1: Fixture** - append to `tests/io/conftest.py`:
  ```python
  import uuid
  import zipfile


  @pytest.fixture
  def make_qza(tmp_path):
      """Build a minimal .qza: <uuid>/metadata.yaml, <uuid>/VERSION, <uuid>/data/<payload>."""

      def make(name, payload, content, semantic_type):
          uid = str(uuid.uuid5(uuid.NAMESPACE_URL, name))  # deterministic (R11.4)
          path = tmp_path / f"{name}.qza"
          with zipfile.ZipFile(path, "w") as archive:
              archive.writestr(f"{uid}/metadata.yaml", f"uuid: {uid}\ntype: {semantic_type}\nformat: null\n")
              archive.writestr(f"{uid}/VERSION", "QIIME 2\narchive: 5\nframework: 2024.10\n")
              archive.writestr(f"{uid}/data/{payload}", content)
          return path

      return make
  ```
- [ ] **Step 2: Failing tests** - `tests/io/test_qiime2.py`:
  ```python
  import numpy as np
  import pytest

  import biotapy as bt

  NEWICK = "((OTU_1:0.1,OTU_2:0.2):0.05,(OTU_3:0.3,OTU_4:0.4):0.1);"
  TAXONOMY = (
      "Feature ID\tTaxon\tConfidence\n"
      "OTU_1\td__Bacteria; p__Firmicutes; c__Clostridia\t0.98\n"
      "OTU_2\td__Bacteria; p__Firmicutes\t0.7\n"
      "OTU_3\tUnassigned\t0.5\n"
  )
  METADATA = (
      "# written by hand\n"
      "sample-id\tdepth\tsite\n"
      "#q2:types\tnumeric\tcategorical\n"
      "S1\t1000\tgut\n"
      "S2\t\tskin\n"
      "\n"
      "S3\t500\tgut\n"
      "S9\t1\tgut\n"
  )


  @pytest.fixture
  def table_qza(make_qza, biom_hdf5):
      return make_qza("table", "feature-table.biom", biom_hdf5.read_bytes(), "FeatureTable[Frequency]")


  def test_read_qiime2_table(table_qza):
      tdata = bt.io.read_qiime2(table_qza)
      assert tdata.shape == (3, 4) and list(tdata.obs_names) == ["S1", "S2", "S3"]
      assert tdata.uns["biotapy"]["x_kind"] == "counts"


  def test_read_qiime2_taxonomy(table_qza, make_qza):
      taxonomy = make_qza("taxonomy", "taxonomy.tsv", TAXONOMY.encode(), "FeatureData[Taxonomy]")
      var = bt.io.read_qiime2(table_qza, taxonomy=taxonomy).var
      assert var.loc["OTU_1", ["kingdom", "phylum", "class"]].tolist() == ["Bacteria", "Firmicutes", "Clostridia"]
      assert var.loc["OTU_3", "kingdom"] == "Unassigned"
      assert var.loc["OTU_4"].isna().all()
      assert var["confidence"].tolist()[:3] == [0.98, 0.7, 0.5]


  def test_read_qiime2_taxonomy_without_confidence(table_qza, make_qza):
      taxonomy = make_qza("taxonomy", "taxonomy.tsv", b"Feature ID\tTaxon\nOTU_1\tk__Bacteria\n", "FeatureData[Taxonomy]")
      var = bt.io.read_qiime2(table_qza, taxonomy=taxonomy).var
      assert "confidence" not in var.columns and var.loc["OTU_1", "kingdom"] == "Bacteria"


  def test_read_qiime2_tree(table_qza, make_qza):
      tree = make_qza("tree", "tree.nwk", NEWICK.encode(), "Phylogeny[Rooted]")
      phylo = bt.io.read_qiime2(table_qza, tree=tree).vart["phylo"]
      assert {n for n in phylo.nodes if phylo.out_degree(n) == 0} == {"OTU_1", "OTU_2", "OTU_3", "OTU_4"}


  def test_read_qiime2_wrong_artifact_names_the_argument(make_qza):
      taxonomy = make_qza("taxonomy", "taxonomy.tsv", TAXONOMY.encode(), "FeatureData[Taxonomy]")
      with pytest.raises(ValueError, match=r"table=.*feature-table\.biom"):
          bt.io.read_qiime2(taxonomy)


  def test_read_qiime2_metadata(table_qza, tmp_path):
      path = tmp_path / "metadata.tsv"
      path.write_text(METADATA)
      obs = bt.io.read_qiime2(table_qza, metadata=path).obs
      assert list(obs.index) == ["S1", "S2", "S3"]
      assert obs["site"].tolist() == ["gut", "skin", "gut"]
      assert obs["depth"].dtype.kind == "f" and np.isnan(obs.loc["S2", "depth"])


  @pytest.mark.parametrize("header", ["id", "SampleID", "Sample-ID", "#SampleID", "sample_name"])
  def test_read_qiime2_metadata_id_headers(table_qza, tmp_path, header):
      path = tmp_path / "metadata.tsv"
      path.write_text(f"{header}\tdepth\nS1\t10\nS2\t20\nS3\t30\n")
      assert bt.io.read_qiime2(table_qza, metadata=path).obs["depth"].tolist() == [10, 20, 30]


  def test_read_qiime2_metadata_legacy_header_is_case_sensitive(table_qza, tmp_path):
      path = tmp_path / "metadata.tsv"
      path.write_text("#sampleid\tdepth\nS1\t10\n")
      with pytest.raises(ValueError, match="metadata="):
          bt.io.read_qiime2(table_qza, metadata=path)
  ```
- [ ] **Step 3: Run, expect failure** -> `AttributeError: ... 'read_qiime2'`.
- [ ] **Step 4: Implement** `src/biotapy/io/_qiime2.py`:
  ```python
  """QIIME 2 artifacts (.qza) and metadata files, read without a QIIME 2 install."""

  import csv
  import tempfile
  import zipfile
  from pathlib import Path

  import biom
  import numpy as np
  import pandas as pd

  from biotapy._core import TreeData, make_treedata, split_lineage, tree_from_newick

  from ._biom import _biom_parts

  _ID_HEADERS_ANY_CASE = frozenset({"id", "sampleid", "sample id", "sample-id", "featureid", "feature id", "feature-id"})
  _ID_HEADERS_EXACT = frozenset({"#SampleID", "#Sample ID", "#OTUID", "#OTU ID", "sample_name"})


  def read_qiime2(
      table: str | Path,
      *,
      taxonomy: str | Path | None = None,
      tree: str | Path | None = None,
      metadata: str | Path | None = None,
  ) -> TreeData:
      """Read QIIME 2 artifacts into one TreeData, without QIIME 2 installed.

      Parameters
      ----------
      table
          ``FeatureTable[Frequency]`` artifact (``.qza``).
      taxonomy
          ``FeatureData[Taxonomy]`` artifact; replaces any taxonomy in the table.
      tree
          ``Phylogeny[Rooted]`` or ``Phylogeny[Unrooted]`` artifact.
      metadata
          QIIME 2 sample metadata file (TSV); samples not in the table are ignored.

      Returns
      -------
      TreeData
          Counts in ``X``; rank columns and ``confidence`` in ``var``; metadata
          in ``obs``; the tree in ``vart['phylo']``.

      Raises
      ------
      ValueError
          An artifact lacks the expected payload, or the metadata has no ID header.

      Notes
      -----
      R equivalent: ``qiime2R::qza_to_phyloseq``
      Guide: :doc:`/guide/reading_data`

      Artifacts are recognized by their payload file, not by ``metadata.yaml``.

      Examples
      --------
      >>> import tempfile
      >>> import zipfile
      >>> from pathlib import Path
      >>> import biom
      >>> from biom.util import biom_open
      >>> import biotapy as bt
      >>> toy = bt.datasets.toy()
      >>> folder = Path(tempfile.mkdtemp())
      >>> with biom_open(str(folder / "t.biom"), "w") as handle:
      ...     biom.Table(toy.X.T, list(toy.var_names), list(toy.obs_names)).to_hdf5(handle, "example")
      >>> with zipfile.ZipFile(folder / "table.qza", "w") as archive:
      ...     archive.write(folder / "t.biom", "0000/data/feature-table.biom")
      >>> bt.io.read_qiime2(folder / "table.qza").shape
      (6, 8)
      """
      with tempfile.TemporaryDirectory() as directory:
          workdir = Path(directory)
          X, obs, var = _biom_parts(biom.load_table(str(_payload(table, "feature-table.biom", workdir, "table"))))
          if taxonomy is not None:
              var = _taxonomy(_payload(taxonomy, "taxonomy.tsv", workdir, "taxonomy")).reindex(var.index)
          phylo = None if tree is None else tree_from_newick(_payload(tree, "tree.nwk", workdir, "tree").read_text())
      if metadata is not None:
          obs = _metadata(Path(metadata)).reindex(obs.index)
      return make_treedata(X, obs=obs, var=var, tree=phylo, x_kind="counts", source="io.read_qiime2")


  def _payload(artifact: str | Path, filename: str, directory: Path, argument: str) -> Path:
      """Extract ``<uuid>/data/<filename>`` from a .qza into ``directory``."""
      with zipfile.ZipFile(artifact) as archive:
          members = [name for name in archive.namelist() if Path(name).parts[1:] == ("data", filename)]
          if len(members) != 1:
              msg = f"{argument}={str(artifact)!r} is not a QIIME 2 artifact holding data/{filename}"
              raise ValueError(msg)
          return Path(archive.extract(members[0], directory))


  def _taxonomy(path: Path) -> pd.DataFrame:
      frame = pd.read_csv(path, sep="\t", index_col=0, dtype=str)
      ranks = split_lineage(frame["Taxon"])
      if "Confidence" in frame.columns:
          ranks["confidence"] = pd.to_numeric(frame["Confidence"], errors="coerce")
      return ranks


  def _is_id_header(cell: str) -> bool:
      return cell in _ID_HEADERS_EXACT or cell.lower() in _ID_HEADERS_ANY_CASE


  def _metadata_rows(path: Path) -> tuple[list[str], list[list[str]], list[str] | None]:
      with path.open(newline="", encoding="utf-8-sig") as handle:
          rows = [[cell.strip() for cell in row] for row in csv.reader(handle, delimiter="\t")]
      rows = [row for row in rows if any(row)]
      while rows and rows[0][0].startswith("#") and not _is_id_header(rows[0][0]):
          rows = rows[1:]
      if not rows or not _is_id_header(rows[0][0]):
          msg = f"metadata={str(path)!r} has no QIIME 2 ID header (e.g. 'sample-id', 'id', '#SampleID')"
          raise ValueError(msg)
      header, body = rows[0], rows[1:]
      types = next((row for row in body if row[0] == "#q2:types"), None)
      data = [row + [""] * (len(header) - len(row)) for row in body if not row[0].startswith("#")]
      return header, data, types


  def _typed(values: pd.Series, declared: str) -> pd.Series:
      values = values.mask(values == "", np.nan)
      if declared == "categorical":
          return values
      numeric = pd.to_numeric(values, errors="coerce")
      # QIIME 2 infers numeric when every present value parses as a number.
      if declared == "numeric" or numeric.notna().sum() == values.notna().sum():
          return numeric
      return values


  def _metadata(path: Path) -> pd.DataFrame:
      header, data, types = _metadata_rows(path)
      frame = pd.DataFrame([row[1:] for row in data], index=[row[0] for row in data], columns=header[1:], dtype=object)
      declared = [t.lower() for t in types[1:]] if types else [""] * len(frame.columns)
      return pd.DataFrame({c: _typed(frame[c], d) for c, d in zip(frame.columns, declared, strict=False)}, index=frame.index)
  ```
  - The controller prototyped the metadata functions on 2026-09-26.
  - Export `read_qiime2`.
  - The docstring example writes a bare `0000/` artifact: the reader needs only
    the payload path.
- [ ] **Step 5: Docs** - add a "QIIME 2" section to `docs/guide/reading_data.md`
  (which artifacts are read, the metadata rules, no QIIME 2 install needed), and
  add `io.read_qiime2` to `docs/api.md`.
- [ ] **Step 6: Run, expect pass**; gate; commit `feat(io): read QIIME 2 artifacts and metadata`.

### Task 1.9: `io.read_dada2`

**Files:**
- Create `src/biotapy/io/_dada2.py` and `tests/io/test_dada2.py`.
- Modify:
  - `src/biotapy/_core/_tree.py` (add `relabel_tips`);
  - `src/biotapy/_core/__init__.py`;
  - `tests/core/test_tree.py`;
  - `src/biotapy/io/__init__.py`;
  - `docs/guide/reading_data.md`;
  - `docs/api.md`.

**Interfaces:**
- **Consumes:** `make_treedata`, `normalize_ranks`, `tree_from_newick`.
- **Produces:**
  - `relabel_tips(tree: nx.DiGraph[str], names: Mapping[str, str]) -> nx.DiGraph[str]` in `_core`;
  - `bt.io.read_dada2(seqtab: str | Path, taxa: str | Path | None = None, *, tree: str | Path | None = None) -> TreeData`.

**Formats** (DADA2 tutorial, R `write.table` docs):
- `seqtab` is samples x sequences, and its column names are the sequences.
- `assignTaxonomy` gives sequences x `Kingdom..Species`, with `NA` where a
  rank is unassigned.
- `write.csv` leaves the first header cell empty and writes row names.
- `.rds` input waits for the 1.6 decision (stage 2).

- [ ] **Step 1: Failing tests.** Add to `tests/core/test_tree.py`:
  ```python
  from biotapy._core import relabel_tips


  def test_relabel_tips_renames_listed_tips_only():
      tree = relabel_tips(tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)]), {"a": "x"})
      assert _leaves(tree) == {"x", "b"} and tree.edges["r", "x"]["length"] == 1.0
  ```
  Create `tests/io/test_dada2.py`:
  ```python
  import numpy as np
  import pytest

  import biotapy as bt

  SEQTAB = '"","ACGTACGT","TTGACCAA","GGGCCCAA","CCCCAAAA"\n"S1",10,0,5,0\n"S2",0,0,0,0\n"S3",3,7,1,0\n'
  TAXA = '"","Kingdom","Phylum","Genus"\n"ACGTACGT","Bacteria","Firmicutes","Blautia"\n"TTGACCAA","Bacteria","Bacteroidota",NA\n'


  @pytest.fixture
  def seqtab(tmp_path):
      path = tmp_path / "seqtab.csv"
      path.write_text(SEQTAB)
      return path


  @pytest.fixture
  def taxa(tmp_path):
      path = tmp_path / "taxa.csv"
      path.write_text(TAXA)
      return path


  def test_read_dada2_names_asvs_and_keeps_sequences(seqtab):
      tdata = bt.io.read_dada2(seqtab)
      assert tdata.shape == (3, 4) and list(tdata.var_names) == ["ASV1", "ASV2", "ASV3", "ASV4"]
      assert tdata.var["sequence"].tolist() == ["ACGTACGT", "TTGACCAA", "GGGCCCAA", "CCCCAAAA"]
      np.testing.assert_array_equal(tdata.X.toarray()[0], [10, 0, 5, 0])


  def test_read_dada2_keeps_all_zero_sample_and_feature(seqtab):
      tdata = bt.io.read_dada2(seqtab)
      assert tdata.X[1].nnz == 0 and tdata.X[:, 3].nnz == 0


  def test_read_dada2_normalizes_taxa(seqtab, taxa):
      var = bt.io.read_dada2(seqtab, taxa).var
      assert {"kingdom", "phylum", "genus", "sequence"} <= set(var.columns)
      assert var.loc["ASV1", "genus"] == "Blautia" and np.isnan(var.loc["ASV2", "genus"])
      assert var.loc["ASV3", ["kingdom", "phylum", "genus"]].isna().all()


  def test_read_dada2_reads_tsv(tmp_path):
      path = tmp_path / "seqtab.tsv"
      path.write_text(SEQTAB.replace(",", "\t"))
      assert bt.io.read_dada2(path).shape == (3, 4)


  def test_read_dada2_single_sample(tmp_path):
      path = tmp_path / "one.csv"
      path.write_text('"","ACGTACGT"\n"S1",4\n')
      assert bt.io.read_dada2(path).shape == (1, 1)


  def test_read_dada2_rejects_a_transposed_table(tmp_path):
      path = tmp_path / "asvs_by_samples.csv"
      path.write_text('"","S1","S2"\n"ACGTACGT",1,2\n')
      with pytest.raises(ValueError, match="samples x sequences"):
          bt.io.read_dada2(path)


  def test_read_dada2_tree_tips_named_by_sequence(seqtab, tmp_path):
      path = tmp_path / "tree.nwk"
      path.write_text("(((ACGTACGT:1,TTGACCAA:1):1,GGGCCCAA:1):1,CCCCAAAA:1);")
      phylo = bt.io.read_dada2(seqtab, tree=path).vart["phylo"]
      assert {n for n in phylo.nodes if phylo.out_degree(n) == 0} == {"ASV1", "ASV2", "ASV3", "ASV4"}
  ```
- [ ] **Step 2: Run, expect failure** -> ImportError for `relabel_tips`; `AttributeError` for `read_dada2`.
- [ ] **Step 3: Implement.** In `_core/_tree.py` (add `from collections.abc import Mapping`):
  ```python
  def relabel_tips(tree: nx.DiGraph[str], names: Mapping[str, str]) -> nx.DiGraph[str]:
      """Rename the nodes listed in ``names`` (e.g. sequence -> ASV id); others keep theirs."""
      return nx.relabel_nodes(tree, dict(names), copy=True)
  ```
  `src/biotapy/io/_dada2.py`:
  ```python
  """DADA2 sequence tables and taxonomy, as written by R's write.csv / write.table."""

  import re
  from pathlib import Path

  import pandas as pd

  from biotapy._core import TreeData, make_treedata, normalize_ranks, relabel_tips, tree_from_newick

  _SEQUENCE = re.compile(r"^[ACGTN]+$")


  def read_dada2(seqtab: str | Path, taxa: str | Path | None = None, *, tree: str | Path | None = None) -> TreeData:
      r"""Read a DADA2 sequence table, with optional taxonomy and tree.

      Parameters
      ----------
      seqtab
          CSV or TSV of ``seqtab``/``seqtab.nochim``: samples x sequences, row
          names in the first column (``write.csv(seqtab.nochim, ...)``).
      taxa
          CSV or TSV of ``assignTaxonomy``/``addSpecies`` output: sequences x ranks.
      tree
          Newick file whose tips are sequences or ASV ids.

      Returns
      -------
      TreeData
          Counts in ``X``; features named ``ASV1..n`` with the sequence in
          ``var['sequence']``; rank columns in ``var``; the tree in ``vart['phylo']``.

      Raises
      ------
      ValueError
          ``seqtab``'s column names are not DNA sequences (the table is transposed).

      Notes
      -----
      R equivalent: ``phyloseq::phyloseq``
      Guide: :doc:`/guide/reading_data`

      The text table is read densely once and stored sparse; DADA2 tables are
      small enough for this. ``.rds`` input is not supported yet.

      Examples
      --------
      >>> import tempfile
      >>> from pathlib import Path
      >>> import biotapy as bt
      >>> path = Path(tempfile.mkdtemp()) / "seqtab.csv"
      >>> _ = path.write_text('"","ACGT","TTGA"\n"S1",3,0\n"S2",1,4\n')
      >>> bt.io.read_dada2(path).var_names.tolist()
      ['ASV1', 'ASV2']
      """
      counts = _read_table(Path(seqtab))
      sequences = [str(column) for column in counts.columns]
      if not all(_SEQUENCE.match(sequence) for sequence in sequences):
          msg = f"seqtab={str(seqtab)!r} must be samples x sequences like DADA2's seqtab; its columns are not sequences"
          raise ValueError(msg)
      names = [f"ASV{i}" for i in range(1, len(sequences) + 1)]
      var = pd.DataFrame({"sequence": sequences}, index=names)
      if taxa is not None:
          ranks = normalize_ranks(_read_table(Path(taxa))).reindex(sequences)
          var = ranks.set_axis(names).join(var)
      phylo = None
      if tree is not None:
          phylo = relabel_tips(tree_from_newick(Path(tree).read_text()), dict(zip(sequences, names, strict=True)))
      obs = pd.DataFrame(index=counts.index)
      return make_treedata(counts.to_numpy(), obs=obs, var=var, tree=phylo, x_kind="counts", source="io.read_dada2")


  def _read_table(path: Path) -> pd.DataFrame:
      """A table whose first column holds row names, as R's write.csv writes it."""
      # write.csv leaves the corner cell empty; pandas would name the index "Unnamed: 0".
      return pd.read_csv(path, sep="," if path.suffix.lower() == ".csv" else "\t", index_col=0).rename_axis(index=None)
  ```
  - Export `relabel_tips` from `_core` and `read_dada2` from `io`.
  - The ASV naming and the sequence column follow
    [data-model-slots](/contracts/data-model-slots.md) (`var["sequence"]`).
- [ ] **Step 4: Docs** - add a "DADA2" section to `docs/guide/reading_data.md`
  (the expected orientation, ASV naming, tree tips by sequence), and add
  `io.read_dada2` to `docs/api.md`.
- [ ] **Step 5: Run, expect pass**; gate; commit `feat(io): read DADA2 sequence tables`.

### Checkpoint B1 - stage 1 review
- [ ] Whole-branch review of stage 1 against every contract and the stage-1
  review focus (superpowers:requesting-code-review); fix pass.
- [ ] Knowledge:
  - add a `Module` concept `.knowledge/modules/io.md`;
  - update `.knowledge/modules/core.md` (Newick, alignment, taxonomy normalizer);
  - update `modules/index.md` and the log.
- [ ] The PR's CI is green, including the Python 3.14 jobs that build
  biom-format from source.
- [ ] Ask the user to review stage 1 and approve the stage-2 expansion.

### Stage 2 - expanded after the 1.6 decision

### Task 1.10: `io.read_phyloseq`
- Per the 1.6 decision. **Interface:** `read_phyloseq(path: str | Path) -> TreeData`.
  Transposes when `taxa_are_rows`; `refseq` -> `var["sequence"]`.
- **Done when:** GlobalPatterns loads with 26 samples x 19,216 features and a tree.

### Task 1.11: `datasets.global_patterns()`, `datasets.enterotype()`
- pooch with SHA-256 hashes; source per 1.6 (phyloseq's `data/*.RData`, or
  converted files attached to a biotapy GitHub release).
- **Tests:** marker `network`; one dedicated CI job runs them with a cached pooch dir.
- **Done when:** both load; `x_kind` set (enterotype holds relative abundances).

### Task 1.12: R golden infrastructure
- **Files:** `tests/r/Dockerfile` (pinned `rocker/r-ver` tag + Bioconductor
  phyloseq, mia, vegan, arrow), `tests/r/export_golden.R`,
  `tests/golden/<dataset>/<function>.parquet`, `tests/golden/VERSIONS.txt`,
  `tests/data/esophagus/*` (tiny fixture exported by the same script),
  `.knowledge/playbooks/regenerate-golden-files.md`.
- **Tests:** golden tests for `relative` and `tax_glom` (phylum, genus) on
  GlobalPatterns, marker `golden`, per [r-golden-parity](/contracts/r-golden-parity.md).
- **Done when:** golden tests pass; the playbook regenerates files bit-identically.

### Checkpoint B - review slice 1B; `Module` concepts for `io` and `datasets`.

---

## Slice 1C - Filtering and diversity

### Task 1.13: `pp.filter_features`, `pp.filter_samples`
- **Interfaces:** `filter_features(adata, *, min_prevalence: float | None = None, min_total: float | None = None) -> AnnData`
  (prevalence = fraction of samples with a non-zero count; via `feature_subset`);
  `filter_samples(adata, *, min_depth: float | None = None) -> AnnData` (AnnData indexing, keeps all slots).
- phyloseq's `prune_*`/`subset_*` map to AnnData indexing in the Coming-from-R table; no wrappers (R2.1).
- **Tests:** thresholds inclusive; both `None` raises `ValueError`; golden vs `filter_taxa`.

### Task 1.14: `pp.rarefy`
- **Interface:** `rarefy(adata, depth: int | None = None, *, seed: int | np.random.Generator | None = None) -> AnnData`.
  `require_counts`; default depth = minimum sample depth (phyloseq default);
  samples below depth dropped; per-row `Generator.multivariate_hypergeometric`
  on the row's non-zeros (sampling without replacement; phyloseq defaults to
  `replace=TRUE`, stated in Notes); all-zero features removed via `feature_subset`.
- **Tests:** row sums == depth; no count exceeds the original; same seed ->
  same result; golden invariants only.

### Task 1.15: `tl.alpha`
- **Interface:** `alpha(adata, metrics: Sequence[str] = ("observed_features", "shannon", "simpson", "chao1"), *, inplace: bool = False) -> pd.DataFrame | None`;
  `faith_pd` requires a TreeData (tree passed as `skbio.TreeNode`, conversion in `_core/_tree.py`).
  Dense conversion in bounded row chunks (R6.2). `inplace=True` writes `obs["alpha_<metric>"]`.
- **Check first:** scikit-bio's Shannon log base and Simpson definition vs
  phyloseq `estimate_richness` (natural log, Gini-Simpson); pass parameters explicitly.
- **Tests:** all-zero sample -> NaN, no exception; chunked == unchunked; golden vs `estimate_richness`.

### Task 1.16: `tl.beta`, `tl.unifrac`
- **Interfaces:** `beta(adata, metric: Literal["braycurtis", "jaccard"] = "braycurtis", *, inplace: bool = False) -> pd.DataFrame | None`;
  `unifrac(tdata, *, weighted: bool = False, normalized: bool = True, inplace: bool = False) -> pd.DataFrame | None`
  via `skbio.diversity.beta_diversity` with `tree=` and `taxa=`. `inplace=True` writes `obsp[<metric>]`.
- **Gotchas:** vegan's `jaccard` is quantitative unless `binary=TRUE`, scikit-bio's is presence/absence:
  golden compares against `vegdist(binary=TRUE)` and the docstring says so.
  phyloseq's weighted UniFrac defaults to `normalized=TRUE`; scikit-bio's to `False`; pass it explicitly.
- **Tests:** symmetric, zero diagonal (Hypothesis); golden vs phyloseq `distance` on esophagus and GlobalPatterns.

### Task 1.17: `tl.pcoa`, `tl.nmds`, `tl.permanova`
- **Interfaces:** `pcoa(adata, *, distance: str = "braycurtis", n_components: int = 10, inplace: bool = False)`
  reads `obsp[distance]` (error message says which `tl.beta` call to run),
  writes `obsm["X_pcoa"]` and `uns["biotapy"]["pcoa"]`;
  `nmds(adata, *, distance: str = "braycurtis", n_components: int = 2, seed=None, inplace: bool = False)`
  via `sklearn.manifold.MDS` non-metric on the precomputed distances (check the current signature with Context7 first);
  `permanova(adata, grouping: str, *, distance: str = "braycurtis", permutations: int = 999, seed=None) -> pd.Series`
  via `skbio.stats.distance.permanova`.
- **Tests:** golden per [r-golden-parity](/contracts/r-golden-parity.md) (PCoA up to sign, NMDS stress and Procrustes, PERMANOVA statistic exact, p-value within 0.02 at 9,999 permutations).

### Checkpoint C - review slice 1C; `Module` concept for `tl`.

---

## Slice 1D - Plots, validation, docs, release

### Task 1.18: `pl.bar`, `pl.richness`, `pl.ordination`, `pl.heatmap`
- Each takes `ax: matplotlib.axes.Axes | None = None`, returns `Axes`, reads
  slots only (never computes a diversity or ordination; errors say which `tl` call to run).
  `pl.ordination` also plots the scree values stored by `tl.pcoa`.
- **Tests:** Agg backend; plotted data equals the source slot; input unchanged.

### Task 1.19: Docstring contract test and Coming-from-R table
- `tests/test_docstrings.py`: for every name in every subpackage `__all__`,
  the docstring has exactly one parseable `R equivalent:` line, a `Guide:` line
  and an `Examples` section.
- Generator (Sphinx build hook or pre-build script) writes the table from those
  lines plus `docs/_data/r_idioms.yaml` (idioms such as `nsamples` -> `tdata.n_obs`); output not committed.

### Task 1.20: Docs and the vignette notebook
- Getting started notebook (first ordination in ~10 lines); guide pages for
  diversity and ordination; `docs/tutorials/phyloseq_analysis.ipynb` reproducing
  the in-scope sections of phyloseq's `phyloseq-analysis.Rmd`;[^phyloseq-vignette]
  myst-nb executes notebooks in the docs build (`nb_execution_mode = "cache"`,
  `nb_execution_raise_on_error = True`).

### Task 1.21: asv benchmarks
- `benchmarks/` with `asv.conf.json` and a synthetic generator (5,000 x 50,000,
  density 0.02, seed 0); benchmarks for `relative`, `tax_glom`, `rarefy`, `beta`, `alpha`.
- Baseline recorded in the docs perf page; no speed gate in Phase 1.

### Task 1.22: Knowledge
- `Module` concepts for `pl`; refresh every Phase 1 concept's `commit`;
  `bash scripts/knowledge_stale.sh` reports 0 stale.

### Task 1.23: Release 0.1
- Follow `playbooks/cut-a-release.md`. After release, propose (do not submit
  without approval) a conda-forge recipe and the scverse ecosystem listing.

---

# Exit gate
- [ ] `docs/tutorials/phyloseq_analysis.ipynb` executes in CI and reproduces the
  vignette's in-scope sections on GlobalPatterns, enterotype and esophagus:
  bar plots, richness, UniFrac + PCoA with scree, NMDS, heatmap, `distance()`
  examples. Out-of-scope sections are listed in the notebook as "not in 0.1":
  `plot_tree`, `plot_net`, CCA, DPCoA.
- [ ] All `golden` tests pass per [r-golden-parity](/contracts/r-golden-parity.md).
- [ ] The generated Coming-from-R table maps all 31 functions below.
- [ ] asv baselines recorded.
- [ ] Coverage >= 90% on public functions (rules.md R11.6).
- [ ] `biotapy 0.1.0` on PyPI; Phase 1 `phase_state: done`; Phase 2 active.

## The 31 phyloseq functions the table must cover
| phyloseq | biotapy |
|---|---|
| `otu_table`, `sample_data`, `tax_table`, `phy_tree`, `refseq` | `tdata.X`, `tdata.obs`, `tdata.var[ranks]`, `tdata.vart["phylo"]`, `tdata.var["sequence"]` |
| `nsamples`, `ntaxa`, `sample_names`, `taxa_names` | `n_obs`, `n_vars`, `obs_names`, `var_names` |
| `sample_sums`, `taxa_sums`, `rank_names`, `sample_variables`, `get_taxa_unique` | `X.sum(axis=1)`, `X.sum(axis=0)`, `var.columns`, `obs.columns`, `var[rank].unique()` |
| `prune_taxa`, `prune_samples`, `subset_taxa`, `subset_samples` | AnnData indexing + `.copy()` |
| `filter_taxa` | `bt.pp.filter_features` |
| `transform_sample_counts` | `bt.pp.relative` |
| `rarefy_even_depth` | `bt.pp.rarefy` |
| `tax_glom` | `bt.pp.tax_glom` |
| `estimate_richness` | `bt.tl.alpha` |
| `distance`, `UniFrac` | `bt.tl.beta`, `bt.tl.unifrac` |
| `ordinate` | `bt.tl.pcoa`, `bt.tl.nmds` |
| `plot_bar`, `plot_richness`, `plot_ordination`, `plot_heatmap` | `bt.pl.bar`, `bt.pl.richness`, `bt.pl.ordination`, `bt.pl.heatmap` |
| `import_biom` | `bt.io.read_biom` |

Rows marked "not in 0.1" in the generated table: `tip_glom`, `merge_samples`, `psmelt`, `plot_tree`, `plot_net`.

[^spec]: Python Microbiome Toolkit development report, sections Roadmap and Testing
[^phyloseq-glom]: phyloseq tax_glom and merge_taxa source
[^phyloseq-vignette]: phyloseq-analysis vignette
