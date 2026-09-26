---
type: Phase
title: Phase 1 - Core, a credible phyloseq replacement (0.1)
description: TreeData conventions in _core; io for phyloseq, BIOM, QIIME 2 and DADA2; pp filter, rarefy, relative, tax_glom; tl alpha, beta, UniFrac, PCoA, NMDS, PERMANOVA; pl bar, richness, ordination, heatmap; R golden tests; Coming-from-R table.
tags: [roadmap, core, io, pp, tl, pl]
status: stable
release: "0.1"
phase_state: not-started
effort: 6-8 weeks part-time (spec); slices 1A-1D with checkpoints
depends_on: [/roadmap/phase-0-foundation.md]
paths: ["src/biotapy/**", "tests/**", "docs/**", "benchmarks/**"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: b77a226
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
> (recommended) or superpowers:executing-plans. Tasks 1.1-1.5 have full TDD
> steps. Every later task lists files, interface, tests and done-when; expand
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
| 1.1 | runtime | scipy, pandas, anndata | sparse kernels, slots |
| 1.3 | runtime | treedata `>=0.3.1,<0.4`, networkx | container and tree |
| 1.6 | runtime | rdata | read phyloseq `.rds`/`.RData` |
| 1.7 | runtime | biom-format | BIOM 1.0 JSON and 2.1 HDF5 |
| 1.7 | runtime | scikit-bio `>=0.7.4,<0.8` | Newick parsing, diversity, ordination |
| 1.11 | runtime | pooch | cached dataset downloads |
| 1.12 | test | pyarrow | read parquet golden files |
| 1.17 | runtime | scikit-learn | non-metric MDS (scikit-bio has none) |
| 1.18 | runtime | matplotlib | `pl` |
| 1.21 | dev | asv | benchmarks |

# Review focus
1. **Non-string or duplicated sample/feature ids** from readers (BIOM ids can be ints) -> readers cast to `str` and fail on duplicates naming them. Tests in 1.7-1.9.
2. **Tree tips and table features disagree** -> readers keep the intersection and warn with both counts, never error deep inside TreeData. Tests in 1.7, 1.8.
3. **Same genus name in different lineages** ("uncultured") -> `tax_glom` groups by lineage. Test in 1.5.
4. **All-zero samples** after filtering -> `relative` keeps zeros; `rarefy` drops them; `alpha` returns NaN, never raises. Tests in 1.4, 1.14, 1.15.
5. **Memory on realistic data** (5,000 x 50,000) -> `alpha` densifies in bounded row chunks; `beta` documents its one dense copy. Tests in 1.15-1.16 assert chunking; asv in 1.21 measures it.

---

## Slice 1A - Kernel and first verbs

### Task 1.1: `_core` sparse kernels

**Files:** create `src/biotapy/_core/_matrix.py`, `tests/core/test_matrix.py`;
modify `src/biotapy/_core/__init__.py`, `pyproject.toml` (scipy, pandas, anndata).
**Interfaces (produces):**
- `as_csr(X: sp.spmatrix | sp.sparray | npt.ArrayLike) -> sp.csr_matrix`
- `sum_by(X: sp.csr_matrix, codes: npt.NDArray[np.intp], n_groups: int) -> sp.csr_matrix`
- `argmax_by(values: npt.NDArray[np.float64], codes: npt.NDArray[np.intp]) -> npt.NDArray[np.intp]`

- [ ] **Step 1: Failing tests**
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
- [ ] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_matrix.py -q` -> `ImportError: cannot import name 'argmax_by'`.
- [ ] **Step 3: Implement**
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
- [ ] **Step 4: Run, expect pass** - same command -> 10 passed.
- [ ] **Step 5: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(core): add sparse group-sum and group-argmax kernels"`

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

- [ ] **Step 1: Failing tests**
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
- [ ] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core -q` -> ImportError for `split_ranks`.
- [ ] **Step 3: Implement**
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
- [ ] **Step 4: Run, expect pass** -> 8 new tests pass.
- [ ] **Step 5: Gate and commit** - `uvx prek run --all-files`; `git commit -am "feat(core): add rank splitting, x_kind and provenance slot rules"`

### Task 1.3: `_core` tree helpers and `datasets.toy()`

**Files:** create `src/biotapy/_core/_tree.py`, `src/biotapy/datasets/__init__.py`,
`src/biotapy/datasets/_toy.py`, `tests/core/test_tree.py`,
`tests/datasets/test_toy.py`, `docs/guide/index.md`, `docs/guide/data_model.md`;
modify `src/biotapy/__init__.py`, `_core/__init__.py`, `docs/index.md`, `docs/api.md`,
`pyproject.toml` (treedata, networkx).
**Interfaces (produces):**
- `PHYLO_KEY = "phylo"`; `TreeData` (re-exported type)
- `tree_from_edges(edges: Iterable[tuple[str, str, float]]) -> nx.DiGraph` (edge attribute `length`)
- `get_tree(tdata: TreeData) -> nx.DiGraph`
- `make_treedata(X, *, obs, var, tree: nx.DiGraph | None, x_kind: XKind, source: str) -> TreeData`
- `bt.datasets.toy() -> TreeData`: 6 samples x 8 features. `obs["group"]` A (s1-s3) / B (s4-s6);
  kingdom..genus with `f8` genus missing; phylum totals give archetypes f3, f6, f7.

- [ ] **Step 1: Failing tests**
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
- [ ] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_tree.py tests/datasets -q` -> ImportError.
- [ ] **Step 3: Implement the tree helpers**
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
- [ ] **Step 4: Implement the toy dataset**
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
- [ ] **Step 5: Docs.** `docs/guide/index.md` (toctree of guide pages) linked
  from `docs/index.md`; `docs/guide/data_model.md`: the slot table and the
  samples-as-rows rule, written for users (from [data-model-slots](/contracts/data-model-slots.md)).
  Add `datasets.toy` to `docs/api.md`.
- [ ] **Step 6: Run, expect pass** - `uv run --group test pytest -q` (includes the doctest) and
  `uv run --group doc sphinx-build -W -b html docs docs/_build/html`.
- [ ] **Step 7: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(datasets): add in-memory toy TreeData and core tree helpers"`

### Task 1.4: `pp.relative`

**Files:** create `src/biotapy/pp/__init__.py`, `src/biotapy/pp/_transform.py`,
`tests/pp/test_transform.py`, `docs/guide/transforms.md`; modify
`tests/conftest.py`, `src/biotapy/__init__.py`, `docs/api.md`, `docs/guide/index.md`.
**Interfaces:** consumes `as_csr`, `add_provenance`; produces
`bt.pp.relative(adata: AnnData) -> AnnData` adding `layers["relative"]`; test fixture `assert_unchanged`.

- [ ] **Step 1: Purity fixture** - append to `tests/conftest.py`:
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
- [ ] **Step 2: Failing tests**
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
- [ ] **Step 3: Run, expect failure** - `uv run --group test pytest tests/pp -q` -> `AttributeError: module 'biotapy' has no attribute 'pp'`.
- [ ] **Step 4: Implement**
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
- [ ] **Step 5: Docs** - `docs/guide/transforms.md` (what `relative` does, the
  zero-sample difference from phyloseq); add to guide toctree and `docs/api.md`.
- [ ] **Step 6: Run, expect pass** - tests, doctest and `sphinx-build -W`.
- [ ] **Step 7: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(pp): add relative abundance transform"`

### Task 1.5: `pp.tax_glom`

**Files:** create `src/biotapy/pp/_glom.py`, `tests/pp/test_glom.py`,
`docs/guide/aggregation.md`; modify `src/biotapy/pp/__init__.py`, `docs/api.md`, `docs/guide/index.md`.
**Interfaces:** consumes `split_ranks`, `as_csr`, `sum_by`, `argmax_by`,
`feature_subset`, `add_provenance`; produces
`bt.pp.tax_glom(adata: AnnData, rank: str, *, dropna: bool = True) -> AnnData`.
Semantics follow phyloseq exactly:[^phyloseq-glom] group by the lineage string
joined with `";_;"` (missing -> `"NA"`), archetype = most abundant member (first
on ties), ranks below `rank` set to `NaN`.

- [ ] **Step 1: Failing tests**
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
- [ ] **Step 2: Run, expect failure** - `uv run --group test pytest tests/pp/test_glom.py -q` -> `AttributeError: ... 'tax_glom'`.
- [ ] **Step 3: Implement**
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
- [ ] **Step 4: Docs** - `docs/guide/aggregation.md`: lineage grouping, the
  archetype rule, what happens to the tree and to derived slots; add to toctree and `docs/api.md`.
- [ ] **Step 5: Run, expect pass** - tests, doctests, `sphinx-build -W`.
- [ ] **Step 6: Gate and commit** - `uvx prek run --all-files`; `git add -A && git commit -m "feat(pp): add tax_glom with phyloseq archetype semantics"`

### Checkpoint A
- [ ] Review slice 1A against every contract (superpowers:requesting-code-review).
- [ ] Write `Module` concepts `.knowledge/modules/core.md` and `.knowledge/modules/pp.md`
  (codebase-map templates), replace the "modules - not yet documented" line in
  `.knowledge/index.md` with `* [modules](modules/index.md) - ...`, create `modules/index.md`, log it.
- [ ] Ask the user to review before slice 1B.

---

## Slice 1B - Readers, datasets, golden infrastructure

### Task 1.6: Spike - phyloseq objects through `rdata` (timebox: 1 day)
- **Question:** can `rdata.read_rda`/`read_rds` with `constructor_dict` turn
  phyloseq's `GlobalPatterns.RData` (S4 `phyloseq` with `otu_table`,
  `sample_data`, `taxonomyTable`, ape `phylo`) into arrays and a Newick-able tree?
  `rdata` 1.1.0 supports S4 as `SimpleNamespace`; a real phyloseq object is untested.
- **Output:** a `Decision` concept `decisions/phyloseq-import-route.md`: native
  `rdata` route, or an R export script shipped in the docs (write BIOM + Newick + TSV).
- **Done when:** the decision is written and approved; throwaway spike code is deleted.

### Task 1.7: `io.read_biom`
- **Interface:** `read_biom(path: str | Path, *, tree: str | Path | None = None) -> TreeData`.
- **Files:** `src/biotapy/io/_biom.py`, `src/biotapy/io/_taxonomy.py` (normalizer
  shared by all readers), `_core/_tree.py` gains `tree_from_newick(text: str) -> nx.DiGraph`
  via `skbio.TreeNode.read([text])`, unnamed internal nodes get unique names; tests under `tests/io/`, fixtures `tests/data/` (< 1 MB, generated by a committed script).
- **Tests:** JSON and HDF5 BIOM; samples become rows; observation `taxonomy`
  metadata -> normalized rank columns ([data-model-slots](/contracts/data-model-slots.md) convention 1);
  integer ids cast to `str`; duplicate ids raise naming them; tree tips not in the table and features not in the tree -> intersection plus one warning with both counts.
- **Done when:** tests pass; `x_kind="counts"`; provenance `io.read_biom`.

### Task 1.8: `io.read_qiime2`
- **Interface:** `read_qiime2(table: str | Path, *, taxonomy=None, tree=None, metadata=None) -> TreeData`
  where each argument is a `.qza` (zip; payload under `<uuid>/data/`) or, for `metadata`, a QIIME 2 metadata TSV.
- **Tests:** tiny `.qza` fixtures built with `zipfile`; `#q2:types` row skipped; confidence column kept as `var["confidence"]`.
- **Done when:** reuses 1.7's BIOM parsing and normalizer; no QIIME 2 install needed.

### Task 1.9: `io.read_dada2`
- **Interface:** `read_dada2(seqtab: str | Path, taxa: str | Path | None = None, *, tree=None) -> TreeData`;
  `seqtab`/`taxa` as CSV/TSV or `.rds` matrices (via `rdata`). Sequences go to
  `var["sequence"]`; feature names `ASV1..n`; no transpose (DADA2 is already samples x ASVs).
- **Done when:** a fixture from the DADA2 tutorial shape round-trips.

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
