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
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-03T08:10:00Z }
commit: 2b9fc24
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
> (recommended) or superpowers:executing-plans. Every task through 1.23 has full TDD steps.

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
| 1.10 | runtime | rdata, xarray | read phyloseq `.RData`/`.rds` natively ([phyloseq-import-route](/decisions/phyloseq-import-route.md)) - approved 2026-09-26, added when stage 2 first imports them |
| 1.7c | runtime | biom-format `>=2.1.16` | BIOM 1.0 JSON and 2.1 HDF5; no CPython 3.14 wheels yet, builds from source - approved 2026-09-26 |
| 1.7a | runtime | scikit-bio `>=0.7.4,<0.8` | Newick parsing, diversity, ordination - approved 2026-09-26 |
| 1.11 | runtime | pooch | cached dataset downloads - approved 2026-09-27 |
| 1.12b | test | pyarrow | declined 2026-09-27 - golden files are gzip CSV read by pandas |
| 1.15a | R image | picante (CRAN, the image's pinned P3M snapshot) | Faith PD golden file (`picante::pd`) - approved 2026-09-27 |
| 1.17 | runtime | scikit-learn `>=1.8` | non-metric MDS (scikit-bio has none); 1.8 renamed `dissimilarity` to `metric` - approved 2026-09-27 |
| 1.18 | runtime | matplotlib `>=3.8` | `pl`, imported inside its functions so `import biotapy` stays fast - approved 2026-09-27 |
| 1.21 | dev | asv `>=0.6.6` | benchmarks, with asv's own `uv` environment plugin - approved 2026-09-27 |
| Checkpoint D | runtime | threadpoolctl `>=3.5` | limit `tl.permanova`'s OpenMP F-statistic to one thread (24 s to 0.01 s on a busy machine); already installed through scikit-learn - approved 2026-10-03 |

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
- [x] **Step 6: Done when** the user approves the decision (status -> `stable`)
  and the scratchpad spike files are deleted. (approved 2026-09-26; spike files deleted)

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
  naming them. A single-tip tree becomes a one-node graph (fix round 1).
- `make_treedata(X, *, obs, var, tree, x_kind, source) -> TreeData`: the
  signature is unchanged. New behaviour:
  - obs/var ids are cast to `str`; duplicates raise `ValueError("duplicate <axis> ids: [...]")`;
  - with a tree, only features that are tips are kept, and tips outside the
    table are pruned (ancestors kept, as TreeData subsetting does);
  - any mismatch gives one `UserWarning` naming both counts, attributed to the
    first frame outside biotapy (`skip_file_prefixes`, Python 3.12+);
  - no shared feature raises `ValueError`;
  - the input frames are never mutated.

- [x] **Step 1: Dependency.** Add `"scikit-bio>=0.7.4,<0.8"` to `[project] dependencies`, then run
  `uv sync --group dev --group test --group doc`. On Python 3.14 this builds
  biom-format from source (see stage-1 constraints).
- [x] **Step 2: Failing tests** - append to `tests/core/test_tree.py`:
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
- [x] **Step 3: Run, expect failure** - `uv run --group test pytest tests/core/test_tree.py -q` -> `ImportError: cannot import name 'tree_from_newick'`.
- [x] **Step 4: Implement.** In `src/biotapy/_core/_tree.py`:
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
- [x] **Step 5: Run, expect pass** - `uv run --group test pytest tests/core -q`. The existing `toy()`
  tests still pass with no warning, because toy's tree tips equal its features.
  Also record `uv run python -X importtime -c "import biotapy" 2>&1 | tail -1` in
  the report: rules.md R10.1 wants a measurement before anyone makes skbio
  import lazily.
- [x] **Step 6: Knowledge and gate.**
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

- [x] **Step 1: Failing tests** - append to `tests/core/test_taxonomy.py`:
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
- [x] **Step 2: Run, expect failure** - `uv run --group test pytest tests/core/test_taxonomy.py -q` -> ImportError.
- [x] **Step 3: Implement** in `src/biotapy/_core/_taxonomy.py`. Add imports `re`,
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
- [x] **Step 4: Run, expect pass** -> 10 new tests pass.
- [x] **Step 5: Gate and commit** - `uvx prek run --all-files`; commit `feat(core): normalize taxonomy ranks and split lineage strings`.

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

- [x] **Step 1: Dependency.** Add `"biom-format>=2.1.16"` to `[project] dependencies`, then run `uv sync`.
- [x] **Step 2: Fixtures** - `tests/io/conftest.py`:
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
- [x] **Step 3: Failing tests** - `tests/io/test_biom.py`:
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
- [x] **Step 4: Run, expect failure** - `uv run --group test pytest tests/io -q` -> `AttributeError: module 'biotapy' has no attribute 'io'`.
- [x] **Step 5: Implement** `src/biotapy/io/_biom.py`:
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
- [x] **Step 6: Docs.**
  - `docs/guide/reading_data.md`, titled "Reading and writing data", with a
    "BIOM" section covering: samples become rows (one transpose); taxonomy
    dialects become rank columns; ids are strings; how a tree/table mismatch is
    handled; duplicate ids are rejected by biom-format.
  - Add it to the guide toctree.
  - `docs/api.md` gets an "Input and output" autosummary block (`.. module:: biotapy.io`) listing `io.read_biom`.
- [x] **Step 7: Run, expect pass** - the tests, the doctest and `sphinx-build -W`.
- [x] **Step 8: Knowledge and gate.**
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

- [x] **Step 1: Failing tests** - append to `tests/io/test_biom.py`:
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
- [x] **Step 2: Run, expect failure** -> `AttributeError: ... 'write_biom'`.
- [x] **Step 3: Implement.** Append to `src/biotapy/io/_biom.py`. Add imports
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
- [x] **Step 4: Docs.**
  - Add a "Writing BIOM" subsection to `docs/guide/reading_data.md`: what is not
    written, and why ranks carry prefixes.
  - Add `io.write_biom` to `docs/api.md`.
- [x] **Step 5: Run, expect pass**; gate; commit `feat(io): write BIOM tables`.

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

- [x] **Step 1: Fixture** - append to `tests/io/conftest.py`:
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
- [x] **Step 2: Failing tests** - `tests/io/test_qiime2.py`:
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
- [x] **Step 3: Run, expect failure** -> `AttributeError: ... 'read_qiime2'`.
- [x] **Step 4: Implement** `src/biotapy/io/_qiime2.py`:
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
- [x] **Step 5: Docs** - add a "QIIME 2" section to `docs/guide/reading_data.md`
  (which artifacts are read, the metadata rules, no QIIME 2 install needed), and
  add `io.read_qiime2` to `docs/api.md`.
- [x] **Step 6: Run, expect pass**; gate; commit `feat(io): read QIIME 2 artifacts and metadata`.

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
    (changed at Checkpoint B1: options keyword-only, R3.1; tree tips must be
    sequences. The code block below is left as the historical plan, not the
    current signature - see `bt.io._dada2.read_dada2`.)

**Formats** (DADA2 tutorial, R `write.table` docs):
- `seqtab` is samples x sequences, and its column names are the sequences.
- `assignTaxonomy` gives sequences x `Kingdom..Species`, with `NA` where a
  rank is unassigned.
- `write.csv` leaves the first header cell empty and writes row names.
- `.rds` input waits for the 1.6 decision (stage 2).

- [x] **Step 1: Failing tests.** Add to `tests/core/test_tree.py`:
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
- [x] **Step 2: Run, expect failure** -> ImportError for `relabel_tips`; `AttributeError` for `read_dada2`.
- [x] **Step 3: Implement.** In `_core/_tree.py` (add `from collections.abc import Mapping`):
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
- [x] **Step 4: Docs** - add a "DADA2" section to `docs/guide/reading_data.md`
  (the expected orientation, ASV naming, tree tips by sequence), and add
  `io.read_dada2` to `docs/api.md`.
- [x] **Step 5: Run, expect pass**; gate; commit `feat(io): read DADA2 sequence tables`.

### Checkpoint B1 - stage 1 review
- [x] Whole-branch review of stage 1 against every contract and the stage-1
  review focus (superpowers:requesting-code-review); fix pass.
- [x] Knowledge:
  - add a `Module` concept `.knowledge/modules/io.md`;
  - update `.knowledge/modules/core.md` (Newick, alignment, taxonomy normalizer);
  - update `modules/index.md` and the log.
- [x] The PR's CI is green, including the Python 3.14 jobs that build
  biom-format from source.
- [x] Ask the user to review stage 1 and approve the stage-2 expansion.

### Stage 2 - expanded after the 1.6 decision

Expanded 2026-09-27 from three research passes run in the scratchpad against
rdata 1.1.0, pooch and the rocker images, with a prototype run on
GlobalPatterns, enterotype and esophagus. Execution order:
**1.12a -> 1.10 -> 1.9b -> 1.11 -> 1.12b**. The R container comes first
because it writes the `.RData`/`.rds` test fixtures the readers need.

#### Stage 2 design
- **No third-party data files are committed.** phyloseq is AGPL-3 and
  biotapy is BSD-3.
  - Reader tests use a synthetic phyloseq object: biotapy's own `toy()`
    numbers, written by the R script in the container.
  - GlobalPatterns and enterotype are downloaded at run time by pooch from
    phyloseq's repository, pinned to a commit and a SHA-256.
  - Golden files hold only derived numbers.
- **Where the code goes.**
  - `io/_rdata.py` holds everything `rdata`-specific: the phyloseq
    `constructor_dict`, and reading a plain R matrix (DADA2's `.rds`).
    Both `read_phyloseq` and `read_dada2` use it.
  - Tree building from ape's `phylo` edge matrix goes in
    `_core.tree_from_phylo`, which keeps networkx inside `_core/_tree.py`
    (tree-access). Internal nodes are named with the same collision-free
    scheme as `tree_from_newick`.
- **Golden tests need GlobalPatterns, which is a download**, so they carry
  both the `golden` and `network` markers. One CI job runs
  `pytest -m "network or golden"` with a cached pooch directory.
- **Research facts the tasks rely on** (every one re-checked by the
  implementer, R2.2):
  - **Missing slots.** An unset phyloseq slot comes back as the string
    `"\x01NULL\x01"` with no class, so it never reaches a constructor. The
    `phyloseq` constructor maps it to `None`.
  - **What each slot arrives as.**
    - `otu_table` arrives already shaped as an `xarray.DataArray`; the
      constructor adds `taxa_are_rows`.
    - `taxonomyTable` arrives as a flat, column-major object array; the
      constructor reshapes it from `attrs["dim"]`/`attrs["dimnames"]`.
    - `sample_data` works with `rdata.conversion.dataframe_constructor`
      unchanged.
    - `phylo` arrives as a dict with keys `edge` (1-based, tips `1..Ntip`),
      `edge.length`, `tip.label`, `Nnode`.
  - **Plain `.rds` matrices.**
    - A plain numeric matrix (`seqtab`) reads straight into a `DataArray`.
    - A plain character matrix (`taxa`) has no class. `read_rds` flattens it
      silently, so it is read with `rdata.parser.parse_file` plus
      `convert_attrs`/`convert_char`.
  - **Typing.** mypy `--strict` needs every constructor typed
    `(obj: Any, attrs: Mapping[str, Any])`, narrowed with `cast`, because
    callable parameters are contravariant. `ConstructorDict` is not exported
    by rdata, so a local alias is used. rdata and xarray ship `py.typed`.
  - **Measurements.** GlobalPatterns reads in about 2.5 s; its OTU table has
    104,578 nonzeros out of 499,616 cells (20.9%).

#### Stage 2 global constraints (in addition to stage 1)
- **Runtime deps.**
  - `rdata>=1.1,<2` and `xarray`: approved 2026-09-26, added in 1.10.
  - `pooch`: added in 1.11, pending approval.
- **Golden format:** gzip CSV (`.csv.gz`), not parquet. The user declined
  `pyarrow` (2026-09-27), so tests read the files with plain pandas and the
  R image needs no `arrow`. zlib writes no timestamp in the gzip header, so
  reruns stay byte-identical.
- **Datasets.** Base URL
  `https://raw.githubusercontent.com/joey711/phyloseq/8a6c2350b985afb909428d396081605e9b3e2f0b/data/`,
  with these SHA-256 hashes:
  - `GlobalPatterns.RData`:
    `bea90c3c48275ea874e0c9400b133da1647e4cddd11b39f89a3d8ffd78512d2d`
    (435,652 B)
  - `enterotype.RData`:
    `0701dd010023344a917bd31680f78580c076bf039befe829830dc43bfc56b8db`
    (195,260 B)
- **The R image.**
  - `rocker/r-ver:4.5.3` (Ubuntu noble; CRAN pinned by rocker to the P3M
    snapshot of 2026-04-23) plus Bioconductor 3.22.
  - Only `phyloseq` (which brings Biostrings) is installed now; there is no
    `arrow` (golden files are CSV).
    vegan and mia are added in slice 1C with their golden files (R2.3).
  - Building and running it needs the user's approval, requested with this
    plan.
- **Size cap.** Every file under `tests/data/` and `tests/golden/` stays under
  1 MB (R6.6). A test enforces this.
- **Dataset loader examples**, and `read_phyloseq`'s example, which loads
  GlobalPatterns, are `# doctest: +SKIP` because they download (R11.4). The
  network CI job runs the same calls for real. This is an accepted exception
  to R8.2's "runnable example": the only phyloseq files biotapy may ship are
  test fixtures, and docstrings cannot reach those.

#### Stage 2 review focus
1. **Absent phyloseq slots**, which appear as the NULL sentinel: no taxonomy,
   no sample data, no tree, no refseq. The reader returns empty frames or no
   tree, never crashes. Tests in 1.10.
2. **`taxa_are_rows` either way.** Samples always become rows, with exactly
   one transpose. Tests in 1.10.
3. **Populated `refseq`.** One warning naming `Biostrings::writeXStringSet`;
   the counts, taxonomy and tree still load. Tests in 1.10.
   > Amended 2026-09-27: shipped as a `ValueError` naming the R fix, not a
   > warning - rdata 1.1.0 cannot parse a populated `refseq` at all. See
   > [phyloseq-import-route](/decisions/phyloseq-import-route.md) Amendment.
4. **`.RData` holding several objects.** `name=` selects one; without it,
   several phyloseq objects raise an error that lists their names; none
   raises. Tests in 1.10.
5. **A DADA2 `taxa.rds` character matrix.** Its shape and dimnames are
   recovered, never silently flat; R's `NA` becomes NaN. Tests in 1.9b.

### Task 1.12a: R container, golden files and R-written fixtures

**Files:**
- Create:
  - `tests/r/Dockerfile`
  - `tests/r/export_golden.R`
  - `tests/test_data_files.py`
  - `.knowledge/playbooks/regenerate-golden-files.md`
- Generate with the container:
  - `tests/golden/global_patterns/{relative,tax_glom_phylum,tax_glom_genus}.csv.gz`
  - `tests/golden/VERSIONS.txt`
  - `tests/data/phyloseq/{toy.rds,toy.RData,two_objects.RData,samples_as_rows.rds,with_refseq.rds}`
  - `tests/data/dada2/{seqtab.rds,taxa.rds}`
- Modify:
  - `.knowledge/contracts/r-golden-parity.md`: golden tests also carry
    `network`; they run in the network CI job
  - `.knowledge/playbooks/index.md`
  - `.knowledge/log.md`

**Interfaces (produces):** the files above, with the layouts below.
- `relative.csv.gz` is long format: `sample_id`, `taxon_id`, `value`
  (float64), nonzeros only, sorted by taxon then sample.
- `tax_glom_{phylum,genus}.csv.gz` is dense: a `sample_id` column, then one
  column per kept taxon, in phyloseq's order.
- The fixtures hold biotapy's `toy()` numbers:
  - `toy`: taxa as rows, taxonomy Kingdom..Genus with `f8`'s Genus `NA`,
    sample data `group`, and the toy tree.
  - `toy_b`: `toy` pruned to s1 and s2.
  - `samples_as_rows.rds`: taxa as columns, with sample data only (no
    taxonomy, no tree).
  - `with_refseq.rds`: `toy` plus a `DNAStringSet`.
  - The DADA2 `.rds` files hold the same values as the CSV strings in
    `tests/io/test_dada2.py`.

- [x] **Step 1: Dockerfile** - `tests/r/Dockerfile`:
  ```dockerfile
  # syntax=docker/dockerfile:1
  # biotapy's R golden-file image (contracts/r-golden-parity). Bumping a pin here is its own commit,
  # followed by regenerating every file (playbooks/regenerate-golden-files).
  FROM rocker/r-ver:4.5.3
  # rocker pins CRAN to a dated Posit Package Manager snapshot, so install.packages() is reproducible
  # and installs Ubuntu binaries.
  ENV DEBIAN_FRONTEND=noninteractive
  RUN Rscript -e 'install.packages("BiocManager")' \
      && Rscript -e 'BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)'
  RUN Rscript -e 'BiocManager::install("phyloseq", version = "3.22", ask = FALSE, update = FALSE)'
  WORKDIR /work
  CMD ["Rscript", "tests/r/export_golden.R"]
  ```
- [x] **Step 2: Export script** - `tests/r/export_golden.R`:
  ```r
  # Writes biotapy's R golden files and R-only test fixtures. Run only in tests/r/Dockerfile's image, from the
  # repo root: docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work biotapy-golden
  suppressPackageStartupMessages({
    library(phyloseq)
    library(Biostrings)
  })

  write_golden <- function(df, path) {
    # write.csv keeps 15 significant digits; zlib's gzip header has no timestamp, so reruns are identical.
    con <- gzfile(path, "w")
    write.csv(df, con, row.names = FALSE)
    close(con)
    if (file.size(path) >= 1e6) stop(path, " is ", file.size(path), " bytes; rules.md R6.6 caps data files at 1 MB")
  }

  samples_as_rows <- function(physeq) {
    m <- as(otu_table(physeq), "matrix")
    if (taxa_are_rows(physeq)) t(m) else m
  }

  dense_frame <- function(physeq) {
    m <- samples_as_rows(physeq)
    data.frame(sample_id = rownames(m), m, check.names = FALSE, row.names = NULL)
  }

  for (dir in c("tests/golden/global_patterns", "tests/data/phyloseq", "tests/data/dada2")) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }

  ## GlobalPatterns golden files (derived numbers only)
  data(GlobalPatterns)
  rel <- samples_as_rows(transform_sample_counts(GlobalPatterns, function(x) x / sum(x)))
  nz <- which(rel != 0, arr.ind = TRUE)
  nz <- nz[order(nz[, "col"], nz[, "row"]), , drop = FALSE]
  write_golden(
    data.frame(sample_id = rownames(rel)[nz[, "row"]], taxon_id = colnames(rel)[nz[, "col"]], value = rel[nz]),
    "tests/golden/global_patterns/relative.csv.gz"
  )
  write_golden(dense_frame(tax_glom(GlobalPatterns, "Phylum")), "tests/golden/global_patterns/tax_glom_phylum.csv.gz")
  write_golden(dense_frame(tax_glom(GlobalPatterns, "Genus")), "tests/golden/global_patterns/tax_glom_genus.csv.gz")

  ## Synthetic phyloseq fixtures: biotapy's toy() numbers, no third-party data
  counts <- rbind(
    c(10, 5, 20, 30, 0, 2, 1, 0), c(8, 7, 25, 22, 3, 0, 0, 1), c(12, 4, 18, 35, 1, 5, 2, 0),
    c(2, 1, 5, 10, 0, 40, 15, 3), c(0, 2, 3, 12, 2, 38, 20, 5), c(1, 0, 4, 8, 1, 45, 12, 2)
  )
  dimnames(counts) <- list(paste0("s", 1:6), paste0("f", 1:8))
  firm <- c("Bacteria", "Firmicutes", "Clostridia")
  bact <- c("Bacteria", "Bacteroidota", "Bacteroidia", "Bacteroidales")
  prot <- c("Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae")
  taxonomy <- rbind(
    c(firm, "Lachnospirales", "Lachnospiraceae", "Blautia"), c(firm, "Lachnospirales", "Lachnospiraceae", "Roseburia"),
    c(firm, "Oscillospirales", "Ruminococcaceae", "Faecalibacterium"), c(bact, "Bacteroidaceae", "Bacteroides"),
    c(bact, "Bacteroidaceae", "Bacteroides"), c(bact, "Prevotellaceae", "Prevotella"),
    c(prot, "Escherichia"), c(prot, NA)
  )
  dimnames(taxonomy) <- list(paste0("f", 1:8), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
  newick <- "(((f1:0.1,f2:0.12)n4:0.05,f3:0.2)n1:0.1,((f4:0.02,f5:0.03)n5:0.05,f6:0.15)n2:0.1,(f7:0.1,f8:0.12)n3:0.2)root;"
  groups <- sample_data(data.frame(group = rep(c("A", "B"), each = 3), row.names = paste0("s", 1:6)))

  toy <- phyloseq(otu_table(t(counts), taxa_are_rows = TRUE), tax_table(taxonomy), groups, phy_tree(ape::read.tree(text = newick)))
  toy_b <- prune_samples(c("s1", "s2"), toy)
  saveRDS(toy, "tests/data/phyloseq/toy.rds")
  save(toy, file = "tests/data/phyloseq/toy.RData")
  save(toy, toy_b, file = "tests/data/phyloseq/two_objects.RData")
  saveRDS(phyloseq(otu_table(counts, taxa_are_rows = FALSE), groups), "tests/data/phyloseq/samples_as_rows.rds")
  sequences <- DNAStringSet(setNames(c("ACGT", "ACGA", "ACGC", "ACGG", "TCGT", "TCGA", "TCGC", "TCGG"), paste0("f", 1:8)))
  saveRDS(merge_phyloseq(toy, sequences), "tests/data/phyloseq/with_refseq.rds")

  ## DADA2 .rds fixtures: the same values as the CSV strings in tests/io/test_dada2.py
  seqtab <- rbind(S1 = c(10L, 0L, 5L, 0L), S2 = c(0L, 0L, 0L, 0L), S3 = c(3L, 7L, 1L, 0L))
  colnames(seqtab) <- c("ACGTACGT", "TTGACCAA", "GGGCCCAA", "CCCCAAAA")
  saveRDS(seqtab, "tests/data/dada2/seqtab.rds")
  taxa <- rbind(ACGTACGT = c("Bacteria", "Firmicutes", "Blautia"), TTGACCAA = c("Bacteria", "Bacteroidota", NA))
  colnames(taxa) <- c("Kingdom", "Phylum", "Genus")
  saveRDS(taxa, "tests/data/dada2/taxa.rds")

  writeLines(c(
    R.version.string,
    paste("Bioconductor", as.character(BiocManager::version())),
    paste0("phyloseq ", packageVersion("phyloseq"))
  ), "tests/golden/VERSIONS.txt")
  ```
- [x] **Step 3: Build and run** (the user approved the Docker build with this plan). From the repo root:
  ```bash
  docker build -t biotapy-golden tests/r
  docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work biotapy-golden
  sha256sum tests/golden/global_patterns/*.csv.gz tests/data/phyloseq/* tests/data/dada2/* > /tmp/claude-1000/golden-run1.sha
  docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work biotapy-golden
  sha256sum -c /tmp/claude-1000/golden-run1.sha
  ```
  Expected:
  - The build takes about 10-20 minutes. If a step starts compiling large
    packages from source for much longer than that, stop and report (R14.2).
  - Both runs finish, and the second run's output files hash the same as the
    first run's (bit-identical, contract rule 3).
  - If `relative.csv.gz` hits the 1 MB `stop()`, report the size and stop.
    The controller then rules on a documented subset.
  - Record in the report the file sizes, the build time, and
    `docker image inspect biotapy-golden --format '{{.Size}}'`.
- [x] **Step 4: Size guard test** - `tests/test_data_files.py`:
  ```python
  from pathlib import Path

  import pytest

  TESTS = Path(__file__).parent
  DATA_FILES = sorted(p for d in ("data", "golden") for p in (TESTS / d).rglob("*") if p.is_file())


  @pytest.mark.parametrize("path", DATA_FILES, ids=lambda p: str(p.relative_to(TESTS)))
  def test_data_files_stay_under_one_megabyte(path):
      # rules.md R6.6
      assert path.stat().st_size < 1_000_000


  def test_versions_file_names_r_and_bioconductor():
      text = (TESTS / "golden" / "VERSIONS.txt").read_text()
      assert "R version 4.5.3" in text and "Bioconductor 3.22" in text
  ```
  Run `uv run --group test pytest tests/test_data_files.py -q`. Expect every
  test to pass. The test runs after the files exist; its job is to keep
  catching growth later.
- [x] **Step 5: Knowledge.**
  - Write `.knowledge/playbooks/regenerate-golden-files.md` (`type: Playbook`):
    - when to regenerate: a pin bump, or new golden functions, each in its own
      commit;
    - the build and run commands from Step 3, plus the bit-identical check;
    - the 1 MB rule;
    - never edit golden files by hand.

    Add its line to `playbooks/index.md`.
  - `contracts/r-golden-parity.md`:
    - Enforced-by: golden tests carry `golden` and `network`, because the
      input is downloaded by pooch, and they run in the network CI job;
    - Statement 1: the image installs only what current golden files need.
    - Statement 2: golden files are gzip CSV (`<function>.csv.gz`), not
      parquet; the user declined `pyarrow` on 2026-09-27.
  - Add a log line.
- [x] **Step 6: Gate and commit.**
  - Run `uvx prek run --all-files` and `uv run --group test pytest`.
  - Commit `test(r): add the R golden container, golden files and R-written fixtures`.
  - Stage the generated binaries by explicit path.

### Task 1.10: `io.read_phyloseq`

> **Amended 2026-09-27 (Checkpoint B):** the `refseq` handling below (Steps 2,
> 5, 6) plans a warned skip. The shipped behaviour is a `ValueError` naming
> the R fix instead: rdata 1.1.0 cannot parse a populated `refseq` at all, so
> there is no four-slots-good, one-slot-skipped result to warn and return.
> See [phyloseq-import-route](/decisions/phyloseq-import-route.md) Amendment
> 2026-09-27 and rules.md R7.4. The steps below are the historical plan and
> are not rewritten.

**Files:**
- Create:
  - `src/biotapy/io/_rdata.py`
  - `src/biotapy/io/_phyloseq.py`
  - `tests/io/test_phyloseq.py`
- Modify:
  - `src/biotapy/_core/_tree.py` (add `tree_from_phylo`)
  - `src/biotapy/_core/__init__.py`
  - `tests/core/test_tree.py`
  - `src/biotapy/io/__init__.py`
  - `pyproject.toml`: add `rdata>=1.1,<2` and `xarray`
  - `docs/guide/reading_data.md`, `docs/api.md`
  - `.knowledge/decisions/optional-heavy-dependencies.md`
  - `.knowledge/contracts/tree-access.md` (it gains `tree_from_phylo`)
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:**
  - `make_treedata`, `normalize_ranks`, `infer_x_kind`, `warn_user`,
    `as_csr`;
  - `io._join._join_to`;
  - the fixtures from 1.12a.
- **Produces:**
  - `_core.tree_from_phylo(edge: npt.ArrayLike, lengths: npt.ArrayLike | None, tips: Sequence[str]) -> nx.DiGraph[str]`.
  - `io._rdata`:
    - `PHYLOSEQ_SLOTS`;
    - `load_phyloseq(path: Path, *, name: str | None) -> dict[str, Any]`,
      where each slot is the converted value or `None`;
    - `read_matrix_rds(path: Path) -> pd.DataFrame`, used by 1.9b.
  - `bt.io.read_phyloseq(path: str | Path, *, name: str | None = None) -> TreeData`.

- [x] **Step 1: Dependency.**
  - Add `"rdata>=1.1,<2"` and `"xarray"` to `[project] dependencies`, then run
    `uv sync --group dev --group test --group doc`.
  - Both ship `py.typed`. If mypy disagrees, handle it per P8.
- [x] **Step 2: Failing tests.**
  - Append to `tests/core/test_tree.py`, with imports in the top block (P7):
    ```python
    from biotapy._core import tree_from_phylo


    def test_tree_from_phylo_names_tips_and_internal_nodes():
        tree = tree_from_phylo([[3, 1], [3, 2]], [0.5, 1.0], ["a", "b"])
        assert _leaves(tree) == {"a", "b"} and tree.number_of_nodes() == 3
        root = next(n for n in tree.nodes if tree.in_degree(n) == 0)
        assert tree.edges[root, "a"]["length"] == 0.5


    def test_tree_from_phylo_internal_names_skip_tip_names():
        tree = tree_from_phylo([[3, 1], [3, 2]], None, ["n0", "b"])
        assert tree.number_of_nodes() == 3 and _leaves(tree) == {"n0", "b"}


    def test_tree_from_phylo_without_lengths_uses_nan():
        tree = tree_from_phylo([[3, 1], [3, 2]], None, ["a", "b"])
        assert all(math.isnan(d["length"]) for *_, d in tree.edges(data=True))
    ```
  - `tests/io/test_phyloseq.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pytest

    import biotapy as bt

    PHYLOSEQ = Path(__file__).parents[1] / "data" / "phyloseq"


    def test_read_phyloseq_rds_matches_toy():
        tdata, toy = bt.io.read_phyloseq(PHYLOSEQ / "toy.rds"), bt.datasets.toy()
        assert list(tdata.obs_names) == list(toy.obs_names) and list(tdata.var_names) == list(toy.var_names)
        np.testing.assert_array_equal(tdata.X.toarray(), toy.X.toarray())
        assert tdata.obs["group"].astype(str).tolist() == toy.obs["group"].astype(str).tolist()
        assert tdata.var.loc["f1", "genus"] == "Blautia" and np.isnan(tdata.var.loc["f8", "genus"])
        tree = tdata.vart["phylo"]
        assert {n for n in tree.nodes if tree.out_degree(n) == 0} == {f"f{i}" for i in range(1, 9)}
        assert tdata.uns["biotapy"]["x_kind"] == "counts"


    def test_read_phyloseq_single_object_rdata():
        assert bt.io.read_phyloseq(PHYLOSEQ / "toy.RData").shape == (6, 8)


    def test_read_phyloseq_samples_as_rows_and_null_slots():
        tdata = bt.io.read_phyloseq(PHYLOSEQ / "samples_as_rows.rds")
        np.testing.assert_array_equal(tdata.X.toarray(), bt.datasets.toy().X.toarray())
        assert "phylo" not in tdata.vart and not set(tdata.var.columns) & {"kingdom", "genus"}
        assert "group" in tdata.obs.columns


    def test_read_phyloseq_skips_refseq_with_one_warning():
        with pytest.warns(UserWarning, match="writeXStringSet") as record:
            tdata = bt.io.read_phyloseq(PHYLOSEQ / "with_refseq.rds")
        assert len(record) == 1 and tdata.shape == (6, 8) and "sequence" not in tdata.var.columns


    def test_read_phyloseq_two_objects_need_a_name():
        with pytest.raises(ValueError, match=r"name=.*toy.*toy_b"):
            bt.io.read_phyloseq(PHYLOSEQ / "two_objects.RData")
        assert bt.io.read_phyloseq(PHYLOSEQ / "two_objects.RData", name="toy_b").shape == (2, 8)


    def test_read_phyloseq_unknown_name_is_named():
        with pytest.raises(KeyError, match="nope"):
            bt.io.read_phyloseq(PHYLOSEQ / "two_objects.RData", name="nope")


    def test_read_phyloseq_rejects_a_non_phyloseq_file():
        with pytest.raises(ValueError, match="path=.*phyloseq"):
            bt.io.read_phyloseq(Path(__file__).parents[1] / "data" / "dada2" / "seqtab.rds")
    ```
- [x] **Step 3: Run, expect failure** - `uv run --group test pytest tests/core/test_tree.py tests/io/test_phyloseq.py -q` -> ImportError for `tree_from_phylo`; `AttributeError` for `read_phyloseq`.
- [x] **Step 4: Implement `_core.tree_from_phylo`** in `_core/_tree.py` (add `from collections.abc import Sequence`):
  ```python
  def tree_from_phylo(edge: npt.ArrayLike, lengths: npt.ArrayLike | None, tips: Sequence[str]) -> nx.DiGraph[str]:
      """Build a tree from an ape ``phylo`` edge matrix (1-based; tips are ``1..len(tips)``).

      Internal nodes get the same collision-free ``n<i>`` names as ``tree_from_newick``;
      missing branch lengths are NaN.
      """
      edges = np.asarray(edge, dtype=np.int64).reshape(-1, 2)
      names = [str(tip) for tip in tips]
      _require_unique_names(list(names))
      taken = set(names)
      fresh = (name for name in (f"n{i}" for i in itertools.count()) if name not in taken)
      internal = {int(node): next(fresh) for node in np.unique(edges) if node > len(names)}
      label = {**{i + 1: tip for i, tip in enumerate(names)}, **internal}
      length = np.full(len(edges), np.nan) if lengths is None else np.asarray(lengths, dtype=np.float64)
      tree = tree_from_edges(
          (label[int(parent)], label[int(child)], float(value))
          for (parent, child), value in zip(edges, length, strict=True)
      )
      tree.add_nodes_from(names)
      return tree
  ```
  Export it from `_core/__init__.py`. `numpy.typing` is already available as `npt`; import it if it was dropped.
- [x] **Step 5: Implement `io/_rdata.py`.**
  ```python
  """R data files through rdata: phyloseq objects and plain matrices (.RData/.rda/.rds)."""

  from collections.abc import Callable, Mapping
  from pathlib import Path
  from types import SimpleNamespace
  from typing import Any, cast

  import numpy as np
  import pandas as pd
  import rdata
  import xarray as xr
  from rdata.conversion import DEFAULT_CLASS_MAP, convert_attrs, convert_char, dataframe_constructor
  from rdata.parser import RObjectType

  # An unset "...OrNULL" S4 slot is serialised as this marker string without a class,
  # so no constructor sees it; the phyloseq constructor maps it to None.
  _R_NULL = "\x01NULL\x01"
  PHYLOSEQ_SLOTS = ("otu_table", "tax_table", "sam_data", "phy_tree", "refseq")
  # rdata does not export its ConstructorDict alias; constructors take Any because callable
  # parameters are contravariant, and narrow with cast.
  _Constructors = Mapping[str | bytes, Callable[[Any, Mapping[str, Any]], Any]]


  def _otu_table(obj: Any, attrs: Mapping[str, Any]) -> tuple[xr.DataArray, bool]:
      return cast(xr.DataArray, obj), bool(attrs["taxa_are_rows"][0])


  def _char_matrix(obj: Any, attrs: Mapping[str, Any]) -> pd.DataFrame:
      # rdata leaves character matrices flat (column-major) and drops dim/dimnames.
      rows, columns = attrs["dimnames"]
      values = np.reshape(np.asarray(obj, dtype=object), attrs["dim"], order="F")
      return pd.DataFrame(values, index=rows, columns=columns)


  def _passthrough(obj: Any, attrs: Mapping[str, Any]) -> Any:
      return obj


  def _or_none(value: Any) -> Any:
      return None if isinstance(value, str) and value == _R_NULL else value


  def _phyloseq(obj: Any, attrs: Mapping[str, Any]) -> dict[str, Any]:
      slots = cast(SimpleNamespace, obj)
      return {name: _or_none(getattr(slots, name)) for name in PHYLOSEQ_SLOTS}


  _PHYLOSEQ: _Constructors = {
      **DEFAULT_CLASS_MAP,
      "otu_table": _otu_table,
      "taxonomyTable": _char_matrix,
      "sample_data": dataframe_constructor,
      "phylo": _passthrough,
      "phyloseq": _phyloseq,
  }


  def _is_phyloseq(value: object) -> bool:
      return isinstance(value, dict) and tuple(value) == PHYLOSEQ_SLOTS


  def load_phyloseq(path: Path, *, name: str | None) -> dict[str, Any]:
      """The phyloseq object in ``path``: an ``.rds`` file, or ``name`` (or the only one) in an ``.RData``."""
      if path.suffix.lower() == ".rds":
          objects = {"": rdata.read_rds(path, constructor_dict=_PHYLOSEQ)}
      else:
          objects = rdata.read_rda(path, constructor_dict=_PHYLOSEQ)
      found = {key: value for key, value in objects.items() if _is_phyloseq(value)}
      if name is not None:
          if name not in found:
              msg = f"name={name!r} is not a phyloseq object in path={str(path)!r}; found {sorted(found)}"
              raise KeyError(msg)
          return cast(dict[str, Any], found[name])
      if len(found) != 1:
          msg = f"path={str(path)!r} holds {len(found)} phyloseq objects; pass name= with one of {sorted(found)}"
          raise ValueError(msg)
      return cast(dict[str, Any], next(iter(found.values())))


  def read_matrix_rds(path: Path) -> pd.DataFrame:
      """A plain R matrix saved with saveRDS, rows and columns named from its dimnames."""
      parsed = rdata.parser.parse_file(path)
      if parsed.object.info.type is RObjectType.STR:
          attrs = convert_attrs(parsed.object, lambda node: rdata.conversion.convert(node))
          flat = [convert_char(cell, default_encoding=None, force_default_encoding=False) for cell in parsed.object.value]
          return _char_matrix(flat, attrs)
      return cast(xr.DataArray, rdata.conversion.convert(parsed)).to_pandas()
  ```
  - The rdata calls come from the research prototypes (`phyloseq_constructors.py`,
    `plain_matrix.py`). Confirm each signature in the installed rdata (R2.2),
    especially `convert_char`'s keyword names and whether `convert_attrs` wants
    the `RData` wrapper or the node.
  - `load_phyloseq`'s "none found" error must say "phyloseq". The test matches
    `path=.*phyloseq`, and the `len(found) != 1` message already contains it.
  - Stay within R5 limits: split `load_phyloseq`'s name handling into a helper
    if ruff asks.
- [x] **Step 6: Implement `io/_phyloseq.py`.**
  ```python
  """phyloseq objects saved from R (.RData/.rda/.rds), read natively (decisions/phyloseq-import-route)."""

  from pathlib import Path
  from typing import Any

  import numpy as np
  import pandas as pd
  import scipy.sparse as sp

  from biotapy._core import TreeData, as_csr, infer_x_kind, make_treedata, normalize_ranks, tree_from_phylo, warn_user

  from ._join import _join_to
  from ._rdata import load_phyloseq


  def read_phyloseq(path: str | Path, *, name: str | None = None) -> TreeData:
      """Read a phyloseq object saved from R, without R.

      Parameters
      ----------
      path
          ``.rds`` file (``saveRDS``) or ``.RData``/``.rda`` file (``save``).
      name
          Object to read from an ``.RData`` file holding several phyloseq objects.

      Returns
      -------
      TreeData
          ``otu_table`` in ``X`` (samples are rows); ``tax_table`` as rank columns in
          ``var``; ``sample_data`` in ``obs``; ``phy_tree`` in ``vart['phylo']``.

      Raises
      ------
      ValueError
          The file holds no phyloseq object, or several and ``name`` is not given.
      KeyError
          ``name`` is not a phyloseq object in the file.

      Warns
      -----
      UserWarning
          The ``refseq`` slot holds sequences; they are skipped (export them with
          ``Biostrings::writeXStringSet``).

      Notes
      -----
      R equivalent: ``phyloseq::phyloseq``
      Guide: :doc:`/guide/reading_data`

      The OTU table is read densely once (8 bytes x samples x taxa) and stored sparse.

      Examples
      --------
      >>> import biotapy as bt
      >>> tdata = bt.datasets.global_patterns()  # doctest: +SKIP
      >>> tdata.shape  # doctest: +SKIP
      (26, 19216)
      """
      slots = load_phyloseq(Path(path), name=name)
      X, samples, features = _counts(slots["otu_table"], path)
      var = pd.DataFrame(index=features)
      if slots["tax_table"] is not None:
          var = _join_to(normalize_ranks(slots["tax_table"]), features, argument=f"path={str(path)!r} tax_table")
      obs = pd.DataFrame(index=samples)
      if slots["sam_data"] is not None:
          obs = _join_to(slots["sam_data"], samples, argument=f"path={str(path)!r} sam_data")
      if slots["refseq"] is not None:
          warn_user(f"path={str(path)!r}: the refseq slot (sequences) is skipped; export it with Biostrings::writeXStringSet")
      return make_treedata(X, obs=obs, var=var, tree=_tree(slots["phy_tree"]), x_kind=infer_x_kind(X), source="io.read_phyloseq")


  def _counts(otu: Any, path: str | Path) -> tuple[sp.csr_matrix, list[str], list[str]]:
      if otu is None:
          msg = f"path={str(path)!r}: the phyloseq object has no otu_table"
          raise ValueError(msg)
      array, taxa_are_rows = otu
      rows, columns = ([str(v) for v in array.coords[dim].values] for dim in array.dims)
      values = np.asarray(array.values)
      # Samples are rows (R6.1): transpose once when phyloseq stored taxa as rows.
      return (as_csr(values.T), columns, rows) if taxa_are_rows else (as_csr(values), rows, columns)


  def _tree(phylo: Any) -> Any:
      if phylo is None:
          return None
      return tree_from_phylo(phylo["edge"], phylo.get("edge.length"), [str(t) for t in phylo["tip.label"]])
  ```
  - Check `_join_to`'s real signature and argument wording in `io/_join.py`
    (fix pass I2) and adapt the `argument=` text to it.
  - `_tree` returns `nx.DiGraph[str] | None`, but `io` may not import networkx
    (R4.5). Annotate it via `TYPE_CHECKING` only if ruff TID251 allows that in
    `io`. Otherwise add a `_core` type alias such as `Tree = nx.DiGraph[str]`,
    exported from `_core`.
  - Export `read_phyloseq` from `io/__init__.py`.
- [x] **Step 7: Docs.**
  - Add a "phyloseq" section to `docs/guide/reading_data.md`:
    - `.rds` vs `.RData`, and `name=`;
    - what is read and what is skipped (`refseq`);
    - no R needed.
  - Add `io.read_phyloseq` to `docs/api.md`.
- [x] **Step 8: Knowledge.**
  - Dependency decision: one sentence, "Task 1.10 added rdata (`>=1.1,<2`) and xarray".
  - `tree-access` statement 2 gains `tree_from_phylo`.
  - One log line.
  - Record in the report `uv run python -X importtime -c "import biotapy" 2>&1 | tail -1`,
    because xarray adds import time (R10.1).
- [x] **Step 9: Run, expect pass; gate; commit** `feat(io): read phyloseq objects saved from R`.

### Task 1.9b: `.rds` input for `io.read_dada2`
Deferred from 1.9 until the rdata route was approved (1.6) and R-written fixtures existed (1.12a).

**Files:** modify `src/biotapy/io/_dada2.py`, `tests/io/test_dada2.py`,
`docs/guide/reading_data.md`.
**Interfaces:** consumes `io._rdata.read_matrix_rds`. `read_dada2`'s signature
does not change; `seqtab` and `taxa` may each be a `.rds` file.

- [x] **Step 1: Failing tests** - append to `tests/io/test_dada2.py`, with imports at the top:
  ```python
  from pathlib import Path

  DADA2 = Path(__file__).parents[1] / "data" / "dada2"


  def test_read_dada2_rds_matches_csv(seqtab, taxa):
      # Both taxa tables cover 2 of the 4 sequences, so both reads warn (checked join).
      with pytest.warns(UserWarning, match="taxa="):
          from_rds = bt.io.read_dada2(DADA2 / "seqtab.rds", taxa=DADA2 / "taxa.rds")
      with pytest.warns(UserWarning, match="taxa="):
          from_csv = bt.io.read_dada2(seqtab, taxa=taxa)
      np.testing.assert_array_equal(from_rds.X.toarray(), from_csv.X.toarray())
      assert list(from_rds.obs_names) == list(from_csv.obs_names)
      assert from_rds.var["sequence"].tolist() == from_csv.var["sequence"].tolist()


  def test_read_dada2_rds_taxa_keeps_matrix_shape():
      with pytest.warns(UserWarning, match="taxa="):
          var = bt.io.read_dada2(DADA2 / "seqtab.rds", taxa=DADA2 / "taxa.rds").var
      assert var.loc["ASV1", ["kingdom", "phylum", "genus"]].tolist() == ["Bacteria", "Firmicutes", "Blautia"]
      assert np.isnan(var.loc["ASV2", "genus"])
  ```
  - Confirm that the warning text in `io/_join.py` matches `taxa=`, and adapt
    the pattern if it does not.
  - Never leave a warning unasserted: the suite must pass under `-W error::UserWarning`.
- [x] **Step 2: Run, expect failure** (the `.rds` file is parsed as text).
- [x] **Step 3: Implement.** In `_dada2.py`'s `_read_table`, dispatch `.rds` to
  `read_matrix_rds`. The `.rds` values come back typed; confirm that the counts
  are integers or whole floats (`infer_x_kind` handles both) and that `NA`
  becomes None/NaN before `normalize_ranks`. Update the docstring: "CSV, TSV
  or `.rds`".
- [x] **Step 4: Docs** - the DADA2 guide section mentions `.rds`.
- [x] **Step 5: Run, expect pass; gate; commit** `feat(io): read DADA2 .rds tables`.

### Task 1.11: `datasets.global_patterns()`, `datasets.enterotype()`

**Files:**
- Create `src/biotapy/datasets/_remote.py`, `tests/datasets/test_remote.py`.
- Modify:
  - `src/biotapy/datasets/__init__.py`
  - `pyproject.toml` (`pooch`)
  - `.github/workflows/test.yaml` (network job)
  - `docs/api.md`, `docs/guide/data_model.md` or a datasets section
  - `.knowledge/modules/datasets.md`
  - `.knowledge/decisions/optional-heavy-dependencies.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:** `bt.io.read_phyloseq`.
- **Produces:** `bt.datasets.global_patterns() -> TreeData` and
  `bt.datasets.enterotype() -> TreeData`.

- [x] **Step 1: Dependency.** Add `"pooch"` to `[project] dependencies`, pending user approval, then run `uv sync`.
- [x] **Step 2: Failing tests** - `tests/datasets/test_remote.py`:
  ```python
  from pathlib import Path

  import pytest

  import biotapy as bt
  from biotapy.datasets import _remote

  TOY_RDS = Path(__file__).parents[1] / "data" / "phyloseq" / "toy.rds"


  def test_loaders_read_the_fetched_file(monkeypatch):
      fetched = []
      monkeypatch.setattr(_remote, "_fetch", lambda name: fetched.append(name) or str(TOY_RDS))
      assert bt.datasets.global_patterns().shape == (6, 8)
      assert bt.datasets.enterotype().shape == (6, 8)
      assert fetched == ["GlobalPatterns.RData", "enterotype.RData"]


  @pytest.mark.network
  def test_global_patterns_downloads_and_loads():
      tdata = bt.datasets.global_patterns()
      tree = tdata.vart["phylo"]
      assert tdata.shape == (26, 19216) and sum(1 for n in tree.nodes if tree.out_degree(n) == 0) == 19216
      assert tdata.uns["biotapy"]["x_kind"] == "counts" and "phylum" in tdata.var.columns


  @pytest.mark.network
  def test_enterotype_downloads_as_relative_abundance():
      tdata = bt.datasets.enterotype()
      assert tdata.shape == (280, 553) and tdata.uns["biotapy"]["x_kind"] == "relative"
      assert "phylo" not in tdata.vart
  ```
  The first test touches a private helper, which R4.9 allows only for `_core`.
  That is a deliberate exception for network isolation (R11.4): it is the
  only way to test the loaders offline without committing third-party data.
  Say so in a one-line comment in the test.
- [x] **Step 3: Run, expect failure** (`ImportError` for `_remote`).
- [x] **Step 4: Implement** `src/biotapy/datasets/_remote.py`:
  ```python
  """Example datasets from phyloseq's repository, downloaded once and cached with pooch."""

  from functools import cache

  import pooch

  from biotapy._core import TreeData
  from biotapy.io import read_phyloseq

  # Pinned to one phyloseq commit so the SHA-256 hashes stay valid.
  _BASE_URL = "https://raw.githubusercontent.com/joey711/phyloseq/8a6c2350b985afb909428d396081605e9b3e2f0b/data/"
  _REGISTRY = {
      "GlobalPatterns.RData": "sha256:bea90c3c48275ea874e0c9400b133da1647e4cddd11b39f89a3d8ffd78512d2d",
      "enterotype.RData": "sha256:0701dd010023344a917bd31680f78580c076bf039befe829830dc43bfc56b8db",
  }


  @cache
  def _pooch() -> pooch.Pooch:
      # BIOTAPY_DATA_DIR overrides the per-user cache directory.
      return pooch.create(path=pooch.os_cache("biotapy"), base_url=_BASE_URL, registry=_REGISTRY, env="BIOTAPY_DATA_DIR")


  def _fetch(name: str) -> str:
      return str(_pooch().fetch(name))


  def global_patterns() -> TreeData:
      """GlobalPatterns: 26 samples from 9 environments, 19,216 OTUs, with taxonomy and a tree.

      Downloaded once (435 kB) from phyloseq's repository and cached.

      Returns
      -------
      TreeData
          Counts in ``X``; sample types in ``obs``; kingdom to species in ``var``;
          the tree in ``vart['phylo']``.

      Notes
      -----
      R equivalent: ``phyloseq::data(GlobalPatterns)``
      Guide: :doc:`/guide/reading_data`

      References
      ----------
      Caporaso JG et al. (2011) Global patterns of 16S rRNA diversity at a depth of millions
      of sequences per sample. PNAS 108:4516-4522.

      Examples
      --------
      >>> import biotapy as bt
      >>> bt.datasets.global_patterns().shape  # doctest: +SKIP
      (26, 19216)
      """
      return read_phyloseq(_fetch("GlobalPatterns.RData"))


  def enterotype() -> TreeData:
      """enterotype: 280 gut samples, 553 genera, relative abundances; no tree.

      Returns
      -------
      TreeData
          Relative abundances in ``X`` (``x_kind == 'relative'``); sample metadata in ``obs``.

      Notes
      -----
      R equivalent: ``phyloseq::data(enterotype)``
      Guide: :doc:`/guide/reading_data`

      References
      ----------
      Arumugam M et al. (2011) Enterotypes of the human gut microbiome. Nature 473:174-180.

      Examples
      --------
      >>> import biotapy as bt
      >>> bt.datasets.enterotype().shape  # doctest: +SKIP
      (280, 553)
      """
      return read_phyloseq(_fetch("enterotype.RData"))
  ```
  - Export both from `datasets/__init__.py`.
  - Confirm the `pooch.create`/`fetch` signatures in the installed pooch (R2.2).
  - `@cache` keeps construction lazy, so importing does no work (R4.7).
- [x] **Step 5: CI job.** Add a `network` job to `.github/workflows/test.yaml`, matching the file's style: pinned action SHAs, `uv`, `shell: bash`.
  - It runs on `ubuntu-latest` with Python 3.13.
  - It caches `BIOTAPY_DATA_DIR=${{ github.workspace }}/.pooch` with
    `actions/cache` pinned to `55cc8345863c7cc4c66a329aec7e433d2d1c52a9 # v6.1.0`
    (research), keyed on `hashFiles('src/biotapy/datasets/_remote.py')`.
  - It runs `uv run --group test pytest -m "network or golden"`. The CLI `-m`
    replaces `addopts`' `-m "not network and not r"`; the research confirmed
    this in pytest's source.
  - Add the job to `check.needs`.
  - Extend `tests/test_ci.py` with one test that the workflow has a job running
    `-m "network or golden"` and that it is in `check.needs`, in the existing
    style.
- [x] **Step 6: Docs and knowledge.**
  - Add both loaders to `docs/api.md` (Datasets block).
  - Add a short "Example datasets" section to the reading guide:
    - they are downloaded and cached;
    - `BIOTAPY_DATA_DIR`;
    - their licences stay with phyloseq's authors.
  - `modules/datasets.md` gains the remote loaders, and the gotcha that
    examples are `+SKIP` because they download.
  - Add the dependency sentence, and a log line.
- [x] **Step 7: Run** `uv run --group test pytest -m network tests/datasets -q` once locally. It downloads about 630 kB. Expect 2 passed; record it in the report.
- [x] **Step 8: Gate** (the normal suite excludes network); commit `feat(datasets): add GlobalPatterns and enterotype via pooch`.

### Task 1.12b: golden tests for `relative` and `tax_glom`

**Files:**
- Create `tests/pp/test_transform_golden.py` and `tests/pp/test_glom_golden.py`.
- No new dependency: pandas reads `.csv.gz` itself (the user declined `pyarrow`).

**Interfaces:** consumes the golden files from 1.12a and `bt.datasets.global_patterns()` from 1.11.

- [x] **Step 1: Tests.**
  - `tests/pp/test_transform_golden.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pytest

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    def test_relative_matches_phyloseq_transform_sample_counts():
        # OTU ids look numeric; read them as text so they match var_names.
        golden = pd.read_csv(GOLDEN / "relative.csv.gz", dtype={"sample_id": str, "taxon_id": str})
        out = bt.pp.relative(bt.datasets.global_patterns())
        rel = out.layers["relative"].tocoo()
        ours = pd.DataFrame(
            {"sample_id": out.obs_names[rel.row], "taxon_id": out.var_names[rel.col], "value": rel.data}
        )
        merged = golden.merge(ours, on=["sample_id", "taxon_id"], how="outer", suffixes=("_r", "_py"), indicator=True)
        assert (merged["_merge"] == "both").all()
        np.testing.assert_allclose(merged["value_py"], merged["value_r"], rtol=1e-7)
    ```
  - `tests/pp/test_glom_golden.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pytest

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    @pytest.mark.parametrize("rank", ["phylum", "genus"])
    def test_tax_glom_matches_phyloseq(rank):
        golden = pd.read_csv(GOLDEN / f"tax_glom_{rank}.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
        out = bt.pp.tax_glom(bt.datasets.global_patterns(), rank)
        assert list(out.obs_names) == list(golden.index)
        assert list(out.var_names) == list(golden.columns)
        np.testing.assert_allclose(out.X.toarray(), golden.to_numpy(), rtol=1e-7)
    ```
- [x] **Step 2: Run** `uv run --group test pytest -m "golden" tests/pp -q`.
  - These tests check against R output that already exists, so they may pass
    first time. Say that in the report and do not fake a RED.
  - If a comparison fails, investigate the semantics before touching code:
    NArm, lineage keys, archetype order. Report findings; never loosen a
    tolerance without a one-line reason (R11.3).
- [x] **Step 3: Gate; commit** `test(pp): compare relative and tax_glom with phyloseq golden files`.

### Checkpoint B - review slice 1B stage 2
- [x] Whole-branch review of stage 2 against every contract and the stage-2
  review focus; fix pass.
- [x] Knowledge: update `modules/io.md` (phyloseq and `.rds`), `modules/datasets.md`
  and `modules/core.md` (`tree_from_phylo`), plus the log.
- [x] The PR's CI is green, including the new network/golden job.
- [x] Ask the user to review slice 1B before slice 1C.

---

## Slice 1C - Filtering and diversity

Expanded 2026-09-27 with superpowers:writing-plans (rules.md R1.2a), from three
research passes and a full prototype. Every file below was written into a
scratch copy of the repo at 11ee8dd and passed ruff, `ruff format`, `mypy
--strict`, import-linter, pytest with doctests, and `sphinx-build -W`. The
golden tests were also run against R output from the `biotapy-golden` image,
with picante installed in a throwaway container. Implementers still re-check
every API they call (R2.2).

The branch then moved to e14e997 (PR #7, docs site). That PR adds
`docs/guide/datasets.md`, points dataset `Guide:` lines at it, and adds pandas
to intersphinx. The docs steps below are written against e14e997, and
`sphinx-build -W` passed there with all of this slice's docs edits applied.

**Goal:** phyloseq's filtering, rarefaction, alpha and beta diversity, UniFrac,
PCoA, NMDS and PERMANOVA on TreeData, each matching R within the
[r-golden-parity](/contracts/r-golden-parity.md) tolerances.

**Architecture:**
- `pp` gains two filters and `rarefy`. They go through `_core.feature_subset`,
  or through AnnData indexing for samples.
- A new `tl` package wraps scikit-bio, plus scikit-learn for NMDS. Each
  function returns its result; with `inplace=True` it writes the
  [data-model-slots](/contracts/data-model-slots.md) keys instead.
- `_core` gains the only networkx-to-`skbio.TreeNode` converter.

**Tech stack:**
- scikit-bio 0.7.4: `alpha_diversity`, `beta_diversity`, `subsample_counts`,
  `pcoa`, `permanova`.
- scikit-learn `>=1.8`: `sklearn.manifold.MDS`.
- SciPy: `procrustes`, used in tests.
- The R image: phyloseq 1.54.2, vegan 2.7.3, ape 5.8.1 and picante 1.8.2.

**Execution order:** **1.15a -> 1.13 -> 1.14 -> 1.15b -> 1.15c -> 1.15 -> 1.16 -> 1.17 -> Checkpoint C.**
- The R exports (1.15a) run first. They need no biotapy code, and every later
  task's golden test reads them.
- The esophagus loader (1.15b) and the tree converter (1.15c) are split out of
  1.15/1.16, because each is its own function with its own tests (R1.4).

### Slice 1C design
- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `pp/_filter.py` | `filter_features`, `filter_samples` |
  | `pp/_rarefy.py` | `rarefy` |
  | `tl/_alpha.py` | `alpha` |
  | `tl/_beta.py` | `beta`, `unifrac`; from 1.17 also `stored_distances`, which `_ordination.py` and `_permanova.py` import (same subpackage, module-boundaries rule 1) |
  | `tl/_ordination.py` | `pcoa`, `nmds` |
  | `tl/_permanova.py` | `permanova` |
  | `_core/_tree.py` | `get_skbio_tree`, the networkx-to-`TreeNode` converter (tree-access) |
  | `_core/_slots.py` | `feature_subset`: keeps `X` (1.13); drops the ordination summaries in `uns["biotapy"]` (1.17) |
  | `datasets/_remote.py` | `esophagus` |

- **What `tl` returns and writes.** Everything defaults to `inplace=False`.
  With `inplace=True` a function writes the keys below and returns `None`
  ([pure-by-default](/decisions/pure-by-default.md)).

  | Function | Returns | `inplace=True` writes |
  |---|---|---|
  | `tl.alpha` | `DataFrame`, samples x metrics | `obs["alpha_<metric>"]` |
  | `tl.beta` | `DataFrame`, samples x samples | `obsp["braycurtis"]` or `obsp["jaccard"]` |
  | `tl.unifrac` | `DataFrame`, samples x samples | `obsp["unweighted_unifrac"]` or `obsp["weighted_unifrac"]` |
  | `tl.pcoa` | `(coords, axes)`: `PC1..` per sample; `eigenvalue`, `proportion_explained` per axis | `obsm["X_pcoa"]`, `uns["biotapy"]["pcoa"] = {"eigenvalues", "proportion_explained"}` |
  | `tl.nmds` | `(coords, stress)`: `NMDS1..` per sample; Kruskal stress-1 | `obsm["X_nmds"]`, `uns["biotapy"]["nmds"] = {"stress"}` (new key, 1.17) |
  | `tl.permanova` | scikit-bio's result `Series` | nothing; no `inplace` (ruling in 1.17) |

- **Dense input.** scikit-bio's `_ingest_table` turns a sparse matrix into a
  0-d object array. It then raises `TypeError: '<class>' is not a supported
  table format.`, for a raw `csr_matrix` and for an AnnData or TreeData with
  sparse `X` alike.
  - `tl.alpha` densifies at most `2**20` values (8 MiB) at a time.
  - `tl.beta` and `tl.unifrac` densify once; pairwise distances need every row.
- **Trees.** scikit-bio's Faith PD and UniFrac raise `The tree must be rooted.`
  when the root has more than two children. `toy()`'s root has three.
  - `get_skbio_tree` keeps the first child and moves the rest under one new
    zero-length node, which changes no root-to-tip length.
  - Unary nodes left by TreeData pruning, and NaN lengths (counted as 0), pass
    through.
- **Randomness.** `as_generator(seed)` is called once per public call.
  - scikit-bio (`subsample_counts`, `permanova`) takes that Generator.
  - scikit-learn takes an `int` drawn from it, because its `random_state`
    rejects Generators.
- **Found while prototyping (pre-existing bug).** On anndata 0.13.4,
  `list(adata.layers.keys())` is `[None]`: `X` itself is `layers[None]`.
  - `_core.feature_subset` deletes every `layers` key, so it deletes `X` and
    returns `X=None`.
  - `pp.tax_glom` hides this by reassigning `out.X`.
  - `filter_features` would ship broken, so Task 1.13 fixes it first, with its
    own test and commit.
- **Research facts the tasks rely on.** All were measured on the prototype
  against the R reference; re-check each (R2.2).
  - **Filters.** On GlobalPatterns, `filter_taxa` keeps 12,711 taxa with
    `sum(x > 0) >= 0.1 * length(x)` and 14,185 with `sum(x) >= 5` (13,711 with
    `> 5`). `filter_features` matches both exactly.
  - **Rarefaction.** At `sample.size = 1e5` phyloseq drops only TRRsed1 (58,688
    reads). `rarefy` drops the same sample and takes 1.7 s for the other 25.
    `skbio.stats.subsample_counts(counts, n, replace=False, seed=...)` accepts a
    Generator and returns `int64`.
  - **Alpha.**
    - `estimate_richness` agrees with scikit-bio within 1.2e-15 relative, for
      Shannon (natural log, `vegan::diversity`), Gini-Simpson and Chao1.
      scikit-bio's `chao1(bias_corrected=True)` is vegan's
      `S.obs + a1 (a1 - 1) / (2 (a2 + 1))`.
    - Faith PD matches `picante::pd(include.root = TRUE)` within 2.2e-16. On
      esophagus, `include.root` TRUE and FALSE agree, because every sample
      spans both root children.
    - An all-zero row gives observed 0, chao1 0, faith_pd 0 and Shannon/Simpson
      NaN in scikit-bio; vegan gives Shannon 0 and Simpson 1.
  - **Beta.** Bray-Curtis, binary Jaccard and both UniFracs, on GlobalPatterns
    and esophagus, agree within 6.7e-16 absolute. UniFrac takes 0.4 s on
    GlobalPatterns, and both datasets' trees are rooted (`ape::is.rooted`).
  - **PCoA.**
    - GlobalPatterns Bray-Curtis has no negative eigenvalues (ape's "no
      correction was applied" branch).
    - Eigenvalues and `Relative_eig` agree within 2e-15 relative, and
      coordinates within 4e-12 after aligning each axis's sign.
    - ape's `Relative_eig` divides by the trace. scikit-bio's
      `proportion_explained` does so only when `dimensions < n`, and otherwise
      divides by the sum of the positive eigenvalues. So `pcoa` asks for at most
      `n_obs - 1` axes.
  - **NMDS.**
    - `metaMDS` reaches stress 0.171061. sklearn with `n_init=20` reaches 0.1736
      (difference 0.0026), with Procrustes r >= 0.99995 for seeds 0 to 3, in
      about 2 s.
    - With `n_init=1`, r is 0.79 to 0.88 for three of four seeds, which fails
      the contract. With `init="classical_mds"`, r is 0.72.
  - **PERMANOVA.** F = 4.425943578211 in both, and p = 1e-4 in both at 9,999
    permutations. scikit-bio takes 9-14 s.
  - **Import time** is unchanged at about 1.3 s; `sklearn.manifold` adds 0.15 s.
  - **esophagus.** sha256
    `0b06d9c35f2e694c34461308af149eb54419453fcb98763de20ab61980b87e46`,
    1,840 B. It reads through `read_phyloseq` with no warning: 3 x 58 float64
    counts (samples B, C, D), no `obs`/`var` columns, and a tree with 114 edges.

### Slice 1C global constraints (in addition to the Phase 1 list)
- **scikit-learn.** Runtime dependency `scikit-learn>=1.8`, approved
  2026-09-27 and added in 1.17. `MDS` is called with every argument explicit:
  `n_components`, `metric_mds=False`, `metric="precomputed"`, `n_init=20`,
  `init="random"`, `normalized_stress="auto"`, `random_state=<int>`.
  - Never pass `dissimilarity=`: it is deprecated in 1.8 and removed in 1.10.
  - `n_init` and `init` must be explicit because their defaults change:
    `n_init` went from 4 to 1 in 1.9, and `init` becomes `"classical_mds"` in
    1.10.
- **The R image.** CRAN `picante` is added from the image's pinned P3M snapshot
  (approved 2026-09-27), in its own commit (1.15a). vegan and ape already come
  with phyloseq; nothing else is added.
- **mypy.** scikit-learn 1.9.1 ships no `py.typed`, so it gets the same
  `follow_untyped_imports` override as scikit-bio (1.17).
  `untyped_calls_exclude` gains `"skbio.stats._subsample"` (1.14) and
  `"sklearn"` (1.17). Never use `# type: ignore` for these.
- **Dense copies.** scikit-bio never receives sparse input. Each wrapper states
  its dense copy and its memory cost in `Notes` (R6.2).
- **Trees.** networkx stays in `_core/_tree.py`; `tl` gets a `TreeNode` only
  from `_core.get_skbio_tree` (R4.5).
- **Keyword-only options.** Every defaulted parameter is keyword-only
  (function-shape statement 1). That includes `alpha(metrics=)`,
  `beta(metric=)` and `rarefy(depth=)`, which the draft interfaces listed as
  positional.
- **Mutation.** `tl` functions never mutate their input unless
  `inplace=True`, and then they write only the data-model-slots keys. No
  `key_added` (R2.3).
- **Golden tolerances** are the r-golden-parity table's, unchanged:
  - filter, alpha, beta and UniFrac: `rtol=1e-7`;
  - PCoA: `rtol=1e-6`, coordinates per axis up to sign;
  - NMDS: stress within 0.02, Procrustes r > 0.95;
  - PERMANOVA: F within `rtol=1e-7`, p within 0.02 at 9,999 permutations;
  - rarefaction: invariants only.
- **Docs.** `docs/conf.py` is `nitpicky`.
  - A public signature returning pandas types needs the pandas intersphinx
    entry. Without it `-W` failed on 11ee8dd; `docs/conf.py` has had one since
    eb7538e (PR #7), so no task adds it.
  - A docstring cross-reference to a function from a later task breaks
    `sphinx-build -W` at the earlier task.
  - Dataset loaders' `Guide:` line points at `/guide/datasets` (b4e7811).
- **Commits** stage explicit paths only. Never stage `.claude/`,
  `.superpowers/`, `.worktrees/` or `notebooks/`.
- **Every task's last commit** also stages `.knowledge/roadmap/phase-1-core.md`
  with that task's boxes ticked, and `.knowledge/log.md` with its dated line
  (R12.4). "Tick ... here" in a Knowledge step means this.

### Slice 1C review focus
1. **All-zero samples** reach every function. Expected:
   - `rarefy`'s default depth skips them, and they are dropped with one warning;
   - `alpha` gives 0 or NaN per metric and never raises;
   - two of them are NaN apart under Bray-Curtis;
   - `pcoa`, `nmds` and `permanova` then raise a `ValueError` naming the NaN,
     not scikit-bio's symmetry error.

   Tests:
   - 1.14 `test_all_zero_sample_is_dropped_by_default`;
   - 1.15 `test_all_zero_sample_is_nan_or_zero_never_raises`;
   - 1.16 `test_all_zero_sample_and_feature`;
   - 1.17 `test_pcoa_nan_distances_raise`.
2. **Stale results after a feature change.** Distances, ordinations and their
   `uns["biotapy"]` summaries must not survive `filter_features` or `rarefy`,
   while `X` must. Tests:
   - 1.13 `test_feature_subset_keeps_x` and
     `test_drops_derived_slots_and_records_provenance`;
   - 1.17 `test_feature_changes_drop_stored_ordinations` and
     `test_feature_subset_drops_ordination_metadata`.
3. **Trees scikit-bio rejects or could mis-measure**: a root with three or
   more children, unary nodes left by filtering, and NaN branch lengths. Faith
   PD and UniFrac must run, with root-to-tip lengths unchanged. Tests:
   - 1.15c `test_get_skbio_tree_*`;
   - 1.16 `test_unifrac_multifurcating_root_is_accepted` and
     `test_unifrac_after_filtering_keeps_path_lengths`.
4. **Stored zeros and decimal thresholds.** An explicitly stored CSR zero is
   not "present", and 3 of 10 samples meets `min_prevalence=0.3`. Tests: 1.13
   `test_explicit_zeros_do_not_count_as_present` and
   `test_decimal_prevalence_boundary_is_kept`.
5. **Calls in the wrong order or shape.** Each error says what to do:
   - ordinating before `tl.beta` raises a `KeyError` naming the exact
     `bt.tl...` call;
   - a bare string as `metrics`, or `faith_pd` on a plain AnnData, raises a
     `TypeError`.

   Tests:
   - 1.17 `test_pcoa_missing_distance_names_the_call` and
     `test_permanova_missing_distance_names_the_call`;
   - 1.15 `test_string_metrics_raise` and `test_faith_pd_needs_a_tree`.

### Task 1.15a: R golden exports for slice 1C

**Files:**
- Modify:
  - `tests/r/Dockerfile`: picante, in its own commit
  - `tests/r/export_golden.R`
  - `.knowledge/contracts/r-golden-parity.md`
  - `.knowledge/playbooks/regenerate-golden-files.md`
  - `.knowledge/log.md`
- Generate with the container:
  - 13 files under `tests/golden/global_patterns/`
  - 2 files under `tests/golden/esophagus/`
  - `tests/golden/VERSIONS.txt`, which gains three lines

**Interfaces:**
- **Consumes:**
  - the `biotapy-golden` image from 1.12a;
  - phyloseq's `data(GlobalPatterns)` and `data(esophagus)`, inside R.
- **Produces:** these gzip CSV files, read by 1.13 to 1.17.

  | File under `tests/golden/` | Columns | R call |
  |---|---|---|
  | `global_patterns/filter_features.csv.gz` | `taxon_id`, `prevalence`, `total` (TRUE/FALSE for all 19,216 taxa) | `filter_taxa(GP, function(x) sum(x > 0) >= 0.1 * length(x))`; `filter_taxa(GP, function(x) sum(x) >= 5)` |
  | `global_patterns/rarefy.csv.gz` | `sample_id`, `sample_sum` (the 25 kept samples; `1e+05` each) | `rarefy_even_depth(GP, sample.size = 1e5, rngseed = FALSE, replace = FALSE)` after `set.seed(20260927)` |
  | `global_patterns/alpha.csv.gz` | `sample_id`, `Observed`, `Chao1`, `se.chao1`, `Shannon`, `Simpson` | `estimate_richness(GP, measures = c("Observed", "Chao1", "Shannon", "Simpson"))` |
  | `global_patterns/alpha_faith_pd.csv.gz` | `sample_id`, `pd`, `sr` | `picante::pd(samples x taxa, phy_tree(GP), include.root = TRUE)` |
  | `global_patterns/beta_braycurtis.csv.gz` | `sample_id`, then one column per sample | `phyloseq::distance(GP, "bray")` |
  | `global_patterns/beta_jaccard.csv.gz` | as above | `phyloseq::distance(GP, "jaccard", binary = TRUE)` |
  | `global_patterns/unifrac_unweighted.csv.gz`, `unifrac_weighted.csv.gz` | as above | `UniFrac(GP, weighted = FALSE)`; `UniFrac(GP, weighted = TRUE, normalized = TRUE)` |
  | `esophagus/unifrac_unweighted.csv.gz`, `unifrac_weighted.csv.gz` | as above (samples B, C, D) | the same on `esophagus` |
  | `global_patterns/pcoa_braycurtis_vectors.csv.gz` | `sample_id`, `Axis.1` .. `Axis.10` | `ordinate(GP, "PCoA", bray)$vectors` (`ape::pcoa`, `correction = "none"`) |
  | `global_patterns/pcoa_braycurtis_values.csv.gz` | `axis`, `eigenvalue`, `relative_eig` (10 rows) | `$values$Eigenvalues`, `$values$Relative_eig` |
  | `global_patterns/nmds_braycurtis_points.csv.gz` | `sample_id`, `MDS1`, `MDS2` | `ordinate(GP, "NMDS", bray)$points` after `set.seed(20260927)` |
  | `global_patterns/nmds_braycurtis_stress.csv.gz` | `stress` (one row) | `$stress` |
  | `global_patterns/permanova_sampletype.csv.gz` | `df`, `sum_of_sqs`, `r2`, `f`, `p` (the `Model` row) | `vegan::adonis2(bray ~ SampleType, data, permutations = 9999)` after `set.seed(20260927)` |

- [x] **Step 1: Record today's hashes.** Before changing anything, from the repo root:
  ```bash
  git ls-files tests/golden tests/data | grep -v VERSIONS.txt | xargs sha256sum > /tmp/claude-1000/golden-before.sha
  ```
  `VERSIONS.txt` is the only existing file this task changes: Step 3 appends three lines.
- [x] **Step 2: Dockerfile, in its own commit** (r-golden-parity statement 1).
  - In `tests/r/Dockerfile`, replace the last `RUN Rscript -e 'stopifnot(...)'` line with:
    ```dockerfile
    # picante (CRAN, same P3M snapshot): the Faith PD golden file. vegan and ape come with phyloseq.
    RUN Rscript -e 'install.packages("picante")'
    RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE))'
    ```
  - Build: `docker build -t biotapy-golden tests/r`.
    - The earlier layers come from the build cache.
    - The image's CRAN repo is `https://p3m.dev/cran/__linux__/noble/2026-04-23`.
      From it, `install.packages("picante")` installs only picante 1.8.2,
      because ape, nlme and vegan are already there. This was verified in a
      throwaway container, not by building the image.
    - If the cache is gone and phyloseq rebuilds, Step 4's before/after check
      catches any drift. If it fails, stop and report (R14.2).
  - Check that `docker run --rm biotapy-golden Rscript -e 'cat(format(packageVersion("picante")))'`
    prints `1.8.2`.
  - Commit: `git add tests/r/Dockerfile && git commit -m "build(r): add picante to the golden image"`.
- [x] **Step 3: Export script.** In `tests/r/export_golden.R`:
  - The directory loop gains the esophagus directory:
    ```r
    for (dir in c("tests/golden/global_patterns", "tests/golden/esophagus", "tests/data/phyloseq", "tests/data/dada2")) {
    ```
  - Insert this block after the `tax_glom_genus.csv.gz` line and before `## Synthetic phyloseq fixtures`:
    ```r
    ## Slice 1C golden files: filtering, rarefaction, diversity, ordination (GlobalPatterns and esophagus)
    square_frame <- function(d) {
      m <- as.matrix(d)
      data.frame(sample_id = rownames(m), m, check.names = FALSE, row.names = NULL)
    }
    gp <- "tests/golden/global_patterns"
    keep_prevalence <- filter_taxa(GlobalPatterns, function(x) sum(x > 0) >= 0.1 * length(x))
    keep_total <- filter_taxa(GlobalPatterns, function(x) sum(x) >= 5)
    write_golden(
      data.frame(taxon_id = names(keep_prevalence), prevalence = keep_prevalence, total = keep_total, row.names = NULL),
      file.path(gp, "filter_features.csv.gz")
    )
    # rngseed = FALSE draws from the global stream seeded here; rngseed = <n> fails before R's first random draw.
    set.seed(20260927)
    rarefied <- rarefy_even_depth(GlobalPatterns, sample.size = 1e5, rngseed = FALSE, replace = FALSE, verbose = FALSE)
    write_golden(
      data.frame(sample_id = sample_names(rarefied), sample_sum = sample_sums(rarefied), row.names = NULL),
      file.path(gp, "rarefy.csv.gz")
    )
    richness <- estimate_richness(GlobalPatterns, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
    write_golden(
      data.frame(sample_id = sample_names(GlobalPatterns), richness, check.names = FALSE, row.names = NULL),
      file.path(gp, "alpha.csv.gz")
    )
    faith <- picante::pd(samples_as_rows(GlobalPatterns), phy_tree(GlobalPatterns), include.root = TRUE)
    write_golden(data.frame(sample_id = rownames(faith), pd = faith$PD, sr = faith$SR, row.names = NULL), file.path(gp, "alpha_faith_pd.csv.gz"))
    # Biostrings (attached after phyloseq) masks distance() with IRanges' generic.
    bray <- phyloseq::distance(GlobalPatterns, "bray")
    write_golden(square_frame(bray), file.path(gp, "beta_braycurtis.csv.gz"))
    # vegdist's jaccard is quantitative unless binary = TRUE; scikit-bio's is presence/absence.
    write_golden(square_frame(phyloseq::distance(GlobalPatterns, "jaccard", binary = TRUE)), file.path(gp, "beta_jaccard.csv.gz"))
    data(esophagus)
    for (name in c("global_patterns", "esophagus")) {
      physeq <- if (name == "esophagus") esophagus else GlobalPatterns
      write_golden(square_frame(UniFrac(physeq, weighted = FALSE)), file.path("tests/golden", name, "unifrac_unweighted.csv.gz"))
      write_golden(square_frame(UniFrac(physeq, weighted = TRUE, normalized = TRUE)), file.path("tests/golden", name, "unifrac_weighted.csv.gz"))
    }
    pcoa <- ordinate(GlobalPatterns, "PCoA", bray)
    axes <- 1:10
    write_golden(
      data.frame(sample_id = rownames(pcoa$vectors), pcoa$vectors[, axes], check.names = FALSE, row.names = NULL),
      file.path(gp, "pcoa_braycurtis_vectors.csv.gz")
    )
    write_golden(
      data.frame(axis = axes, eigenvalue = pcoa$values$Eigenvalues[axes], relative_eig = pcoa$values$Relative_eig[axes]),
      file.path(gp, "pcoa_braycurtis_values.csv.gz")
    )
    set.seed(20260927)
    # ordinate() runs metaMDS on the dist object (no autotransform) and prints every random start; keep the log short.
    invisible(capture.output(nmds <- ordinate(GlobalPatterns, "NMDS", bray)))
    write_golden(
      data.frame(sample_id = rownames(nmds$points), nmds$points, check.names = FALSE, row.names = NULL),
      file.path(gp, "nmds_braycurtis_points.csv.gz")
    )
    write_golden(data.frame(stress = nmds$stress), file.path(gp, "nmds_braycurtis_stress.csv.gz"))
    set.seed(20260927)
    adonis <- vegan::adonis2(bray ~ SampleType, data = data.frame(sample_data(GlobalPatterns)), permutations = 9999)
    write_golden(
      data.frame(df = adonis$Df[1], sum_of_sqs = adonis$SumOfSqs[1], r2 = adonis$R2[1], f = adonis$F[1], p = adonis[["Pr(>F)"]][1]),
      file.path(gp, "permanova_sampletype.csv.gz")
    )
    ```
  - Replace the closing `writeLines` call with:
    ```r
    writeLines(c(
      R.version.string,
      paste("Bioconductor", as.character(BiocManager::version())),
      paste0("phyloseq ", packageVersion("phyloseq")),
      paste0("vegan ", packageVersion("vegan")),
      paste0("ape ", packageVersion("ape")),
      paste0("picante ", packageVersion("picante"))
    ), "tests/golden/VERSIONS.txt")
    ```
  - Why the block looks as it does (each point verified in the container):
    - `phyloseq::distance` is namespaced because Biostrings, attached after
      phyloseq, masks `distance()` with IRanges' generic. Unqualified, it fails:
      `unable to find an inherited method for function 'distance' for signature 'x = "phyloseq", y = "character"'`.
    - Each random step has its own `set.seed(20260927)`, so no file depends on
      section order. `rarefy_even_depth(rngseed = FALSE)` draws from that
      stream; `rngseed = <n>` fails with `object '.Random.seed' not found` in a
      fresh session.
    - `ordinate(GP, "NMDS", bray)` calls `metaMDS(ps.dist)` on the `dist`
      object, so no autotransform happens. It prints all 20 random starts, which
      `capture.output` swallows.
    - `adonis2`'s default `by = NULL` gives one overall test; row 1 is `Model`.
- [x] **Step 4: Run twice and check bit-identical.** From the repo root:
  ```bash
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
  sha256sum tests/golden/*/*.csv.gz tests/golden/VERSIONS.txt tests/data/phyloseq/* tests/data/dada2/* > /tmp/claude-1000/golden-run1.sha
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
  sha256sum -c --quiet /tmp/claude-1000/golden-run1.sha
  sha256sum -c --quiet /tmp/claude-1000/golden-before.sha
  ```
  Expected (as in the prototype):
  - Both checks print nothing and exit 0: the second run is bit-identical to
    the first, and all 14 pre-existing files (3 golden CSVs, 11 fixtures) are
    unchanged.
  - One run takes about 1 minute; Faith PD on GlobalPatterns alone takes about 20 s.
  - The new files' sizes in bytes:

    | File | Bytes |
    |---|---|
    | `filter_features` | 78,594 |
    | `rarefy` | 163 |
    | `alpha` | 1,214 |
    | `alpha_faith_pd` | 376 |
    | `beta_braycurtis` | 3,894 |
    | `beta_jaccard` | 3,955 |
    | `unifrac_unweighted` (GlobalPatterns) | 3,975 |
    | `unifrac_weighted` (GlobalPatterns) | 3,989 |
    | `pcoa_braycurtis_vectors` | 2,635 |
    | `pcoa_braycurtis_values` | 261 |
    | `nmds_braycurtis_points` | 657 |
    | `nmds_braycurtis_stress` | 47 |
    | `permanova_sampletype` | 97 |
    | esophagus `unifrac_unweighted` | 107 |
    | esophagus `unifrac_weighted` | 105 |
  - `VERSIONS.txt` ends with `vegan 2.7.3`, `ape 5.8.1` and `picante 1.8.2`.
  - Spot checks:
    - `zcat tests/golden/global_patterns/permanova_sampletype.csv.gz` shows
      `8,7.77436045388947,0.675619248642622,4.42594357821139,1e-04`;
    - the NMDS stress is `0.171061337156148`;
    - `rarefy.csv.gz` has 25 rows, with no `TRRsed1`.
- [x] **Step 5: Size guard.** Run `uv run --group test pytest tests/test_data_files.py -q`.
  Every file stays under 1 MB; the largest new one is 78,594 B. There is no
  docs step: the files are test data, and the playbook covers how to
  regenerate them.
- [x] **Step 6: Knowledge.**
  - `contracts/r-golden-parity.md` statement 1: replace "(`phyloseq`, which
    brings `Biostrings`)" with "(`phyloseq`, which brings `Biostrings`,
    `vegan` and `ape`, plus CRAN `picante` for Faith PD)".
  - `playbooks/regenerate-golden-files.md`:
    - the Steps' `sha256sum` glob becomes `tests/golden/*/*.csv.gz`, because
      golden files now live in two dataset directories;
    - add two Common mistakes:
      - Biostrings masks `distance()`, so write `phyloseq::distance`;
      - every random call needs its own preceding `set.seed(20260927)`, or
        reruns stop being bit-identical.
  - Add a log line. Tick 1.15a here.
- [x] **Step 7: Gate and commit.**
  - Run `uvx prek run --all-files` and `uv run --group test pytest`.
  - Commit `test(r): export slice 1C golden files for filtering, diversity and ordination`,
    staging `tests/r/export_golden.R`, the 15 new `.csv.gz` files,
    `tests/golden/VERSIONS.txt` and the three knowledge files by explicit path.

### Task 1.13: `pp.filter_features`, `pp.filter_samples`

**Files:**
- Create:
  - `src/biotapy/pp/_filter.py`
  - `tests/pp/test_filter.py`
  - `tests/pp/test_filter_golden.py`
  - `docs/guide/filtering.md`
- Modify:
  - `src/biotapy/_core/_slots.py` and `tests/core/test_slots.py` (the `X` fix,
    its own commit)
  - `src/biotapy/pp/__init__.py`
  - `docs/api.md`, `docs/guide/index.md`
  - `.knowledge/modules/core.md`, `.knowledge/modules/pp.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:**
  - `feature_subset`, `as_csr`, `add_provenance`;
  - `tests/golden/global_patterns/filter_features.csv.gz` (1.15a).
- **Produces:**
  - `_core.feature_subset(adata, index)`: signature unchanged; it now keeps `X`.
  - `bt.pp.filter_features(adata: AnnData, *, min_prevalence: float | None = None, min_total: float | None = None) -> AnnData`.
    Prevalence is the fraction of samples with a stored non-zero value. Both
    thresholds are inclusive, and with both given a feature must pass both. It
    goes through `feature_subset` and records provenance
    `pp.filter_features`. Nothing passing raises `ValueError`.
  - `bt.pp.filter_samples(adata: AnnData, min_depth: float) -> AnnData`.
    `min_depth` is required, so it is positional (R3.1; a `None` default would
    only raise). The threshold is inclusive. It subsets with AnnData indexing,
    keeping every slot, and records provenance `pp.filter_samples`.
  - phyloseq's `prune_*`/`subset_*` map to AnnData indexing in the
    Coming-from-R table; there are no wrappers (R2.1).

- [x] **Step 1: Failing test for the `feature_subset` bug.** Append to `tests/core/test_slots.py`:
  ```python
  def test_feature_subset_keeps_x():
      np.testing.assert_array_equal(feature_subset(_adata(), np.array([0, 2])).X.toarray(), [[0, 2], [3, 5]])
  ```
  Run `uv run --group test pytest tests/core/test_slots.py -q`. It fails with
  `AttributeError: 'NoneType' object has no attribute 'toarray'`.
- [x] **Step 2: Fix.** In `src/biotapy/_core/_slots.py`'s `feature_subset`, replace the inner loop
  (`for key in list(mapping.keys()):` and its `del` line) with:
  ```python
          # anndata 0.13 lists X itself as layers[None]; deleting that key would delete X.
          for key in [key for key in mapping.keys() if key is not None]:
              del mapping[key]
  ```
  - Run the same command and expect it to pass.
  - `test_feature_subset_drops_derived_slots` still passes, because
    `bool(adata.layers)` is `False` when only `None` is listed.
  - Knowledge, in the same commit: add a Gotcha to `.knowledge/modules/core.md`:
    > anndata 0.13 exposes `X` as `layers[None]`: `list(adata.layers.keys())`
    > includes `None`, and deleting that key deletes `X`. `feature_subset`
    > skips it. Before Task 1.13 it returned `X=None`, which `pp.tax_glom` hid
    > by reassigning `out.X` (`_slots.py:feature_subset`,
    > `tests/core/test_slots.py:test_feature_subset_keeps_x`).

    Add a log line.
  - Gate: `uvx prek run --all-files` and `uv run --group test pytest`.
  - Commit `fix(core): keep X when feature_subset drops derived slots`.
- [x] **Step 3: Failing tests.**
  - `tests/pp/test_filter.py`:
    ```python
    import json

    import anndata as ad
    import numpy as np
    import pandas as pd
    import pytest
    import scipy.sparse as sp
    from hypothesis import given
    from hypothesis import strategies as st
    from hypothesis.extra.numpy import arrays

    import biotapy as bt
    from biotapy._core import get_tree


    def _adata(dense) -> ad.AnnData:
        return ad.AnnData(
            X=sp.csr_matrix(dense),
            obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
            var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
        )


    def _tips(tdata) -> set[str]:
        tree = get_tree(tdata)
        return {n for n in tree.nodes if tree.out_degree(n) == 0}


    # toy(): f1..f8 are non-zero in 5, 5, 6, 6, 4, 5, 5, 4 of 6 samples and total 33, 19, 75, 117, 7, 130, 50, 11.


    def test_min_prevalence_is_inclusive():
        # f5 and f8 are in exactly 4 of 6 samples.
        assert bt.pp.filter_features(bt.datasets.toy(), min_prevalence=4 / 6).n_vars == 8
        assert list(bt.pp.filter_features(bt.datasets.toy(), min_prevalence=5 / 6).var_names) == [
            "f1",
            "f2",
            "f3",
            "f4",
            "f6",
            "f7",
        ]


    def test_min_total_is_inclusive():
        out = bt.pp.filter_features(bt.datasets.toy(), min_total=19)
        assert list(out.var_names) == ["f1", "f2", "f3", "f4", "f6", "f7"]


    def test_both_thresholds_must_pass():
        out = bt.pp.filter_features(bt.datasets.toy(), min_prevalence=1.0, min_total=100)
        assert list(out.var_names) == ["f4"]


    def test_decimal_prevalence_boundary_is_kept():
        # 3 / 10 == 0.3 in floating point, while 0.3 * 10 is just above 3.
        dense = np.ones((10, 2), dtype=np.int64)
        dense[3:, 0] = 0
        assert list(bt.pp.filter_features(_adata(dense), min_prevalence=0.3).var_names) == ["f0", "f1"]


    def test_explicit_zeros_do_not_count_as_present():
        adata = _adata(np.array([[1, 3], [2, 4]]))
        adata.X.data[0] = 0  # a stored zero: still 4 stored entries
        assert adata.X.nnz == 4
        assert list(bt.pp.filter_features(adata, min_prevalence=1.0).var_names) == ["f1"]


    def test_all_zero_feature_and_sample():
        dense = np.array([[0, 0, 0], [1, 0, 2], [3, 0, 0]])
        out = bt.pp.filter_features(_adata(dense), min_prevalence=0.1)
        assert list(out.var_names) == ["f0", "f2"] and out.n_obs == 3


    def test_single_sample():
        assert list(bt.pp.filter_features(_adata(np.array([[0, 4, 1]])), min_total=1).var_names) == ["f1", "f2"]


    def test_keeps_missing_ranks_and_prunes_the_tree():
        out = bt.pp.filter_features(bt.datasets.toy(), min_total=11)
        assert "f8" in out.var_names and pd.isna(out.var.loc["f8", "genus"])
        assert _tips(out) == set(out.var_names) == {"f1", "f2", "f3", "f4", "f6", "f7", "f8"}


    def test_drops_derived_slots_and_records_provenance():
        tdata = bt.pp.relative(bt.datasets.toy())
        tdata.obsp["braycurtis"] = np.zeros((6, 6))
        out = bt.pp.filter_features(tdata, min_total=19)
        assert "relative" not in out.layers and "braycurtis" not in out.obsp
        entry = json.loads(out.uns["biotapy"]["provenance"][-1])
        assert entry["step"] == "pp.filter_features" and entry["params"] == {"min_prevalence": None, "min_total": 19}


    def test_input_unchanged(assert_unchanged):
        tdata = bt.datasets.toy()
        before = tdata.copy()
        bt.pp.filter_features(tdata, min_prevalence=0.9)
        assert_unchanged(before, tdata)


    def test_neither_threshold_raises():
        with pytest.raises(ValueError, match="min_prevalence=, min_total="):
            bt.pp.filter_features(bt.datasets.toy())


    def test_prevalence_outside_zero_to_one_raises():
        with pytest.raises(ValueError, match="min_prevalence"):
            bt.pp.filter_features(bt.datasets.toy(), min_prevalence=5)


    def test_nothing_passing_raises():
        with pytest.raises(ValueError, match="no feature passes"):
            bt.pp.filter_features(bt.datasets.toy(), min_total=10_000)


    @given(
        arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 8)), elements=st.integers(0, 20)),
        st.sampled_from([0.0, 0.25, 0.5, 1.0]),
        st.integers(0, 30),
    )
    def test_keeps_exactly_the_features_passing_both(dense, prevalence, total):
        passing = ((dense > 0).mean(axis=0) >= prevalence) & (dense.sum(axis=0) >= total)
        if not passing.any():
            with pytest.raises(ValueError, match="no feature passes"):
                bt.pp.filter_features(_adata(dense), min_prevalence=prevalence, min_total=total)
            return
        out = bt.pp.filter_features(_adata(dense), min_prevalence=prevalence, min_total=total)
        assert list(out.var_names) == [f"f{i}" for i in np.flatnonzero(passing)]
        np.testing.assert_array_equal(out.X.toarray(), dense[:, passing])


    # toy() sample depths: s1..s6 = 68, 66, 77, 76, 82, 73.


    def test_min_depth_is_inclusive():
        assert list(bt.pp.filter_samples(bt.datasets.toy(), 73).obs_names) == ["s3", "s4", "s5", "s6"]


    def test_filter_samples_keeps_every_slot():
        tdata = bt.pp.relative(bt.datasets.toy())
        tdata.obsp["braycurtis"] = np.arange(36.0).reshape(6, 6)
        tdata.obsm["X_pcoa"] = np.arange(12.0).reshape(6, 2)
        out = bt.pp.filter_samples(tdata, 73)
        np.testing.assert_array_equal(out.obsp["braycurtis"], np.arange(36.0).reshape(6, 6)[2:, 2:])
        np.testing.assert_array_equal(out.obsm["X_pcoa"], np.arange(12.0).reshape(6, 2)[2:])
        assert "relative" in out.layers and out.n_vars == 8 and _tips(out) == set(out.var_names)
        assert json.loads(out.uns["biotapy"]["provenance"][-1])["step"] == "pp.filter_samples"


    def test_filter_samples_drops_an_all_zero_sample_and_keeps_all_zero_features():
        out = bt.pp.filter_samples(_adata(np.array([[0, 0], [0, 3]])), 1)
        assert list(out.obs_names) == ["s1"] and out.n_vars == 2


    def test_filter_samples_single_sample():
        assert bt.pp.filter_samples(_adata(np.array([[2, 3]])), 5).n_obs == 1


    def test_filter_samples_input_unchanged(assert_unchanged):
        tdata = bt.datasets.toy()
        before = tdata.copy()
        bt.pp.filter_samples(tdata, 73)
        assert_unchanged(before, tdata)


    def test_filter_samples_nothing_passing_raises():
        with pytest.raises(ValueError, match="min_depth=1000"):
            bt.pp.filter_samples(bt.datasets.toy(), 1000)


    @given(
        arrays(np.int64, st.tuples(st.integers(1, 8), st.integers(1, 5)), elements=st.integers(0, 10)), st.integers(0, 20)
    )
    def test_filter_samples_keeps_exactly_the_deep_samples(dense, depth):
        deep = dense.sum(axis=1) >= depth
        if not deep.any():
            with pytest.raises(ValueError, match="min_depth"):
                bt.pp.filter_samples(_adata(dense), depth)
            return
        assert list(bt.pp.filter_samples(_adata(dense), depth).obs_names) == [f"s{i}" for i in np.flatnonzero(deep)]
    ```
  - `tests/pp/test_filter_golden.py`:
    ```python
    from pathlib import Path

    import pandas as pd
    import pytest

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    @pytest.mark.parametrize(
        ("column", "threshold"), [("prevalence", {"min_prevalence": 0.1}), ("total", {"min_total": 5})]
    )
    def test_filter_features_matches_phyloseq_filter_taxa(column, threshold):
        # prevalence: filter_taxa(GP, function(x) sum(x > 0) >= 0.1 * length(x)); total: function(x) sum(x) >= 5.
        golden = pd.read_csv(GOLDEN / "filter_features.csv.gz", dtype={"taxon_id": str})
        out = bt.pp.filter_features(bt.datasets.global_patterns(), **threshold)
        assert list(out.var_names) == golden.loc[golden[column], "taxon_id"].tolist()
    ```
- [x] **Step 4: Run, expect failure.** Run `uv run --group test pytest tests/pp/test_filter.py -q`.
  It fails with `AttributeError: module 'biotapy.pp' has no attribute 'filter_features'`,
  and with the same error for `'filter_samples'`.
- [x] **Step 5: Implement.** `src/biotapy/pp/_filter.py`:
  ```python
  """Filters that keep a subset of features or samples."""

  import numpy as np
  from anndata import AnnData

  from biotapy._core import add_provenance, as_csr, feature_subset


  def filter_features(adata: AnnData, *, min_prevalence: float | None = None, min_total: float | None = None) -> AnnData:
      """Keep features that are present in enough samples and have enough reads.

      Parameters
      ----------
      adata
          Samples x features.
      min_prevalence
          Keep features that are non-zero in at least this fraction of samples, from 0 to 1.
      min_total
          Keep features whose total over all samples is at least this.

      Returns
      -------
      AnnData
          Same type as ``adata`` with the kept features in their original order; a
          TreeData keeps their subtree. Both thresholds are inclusive and, when both
          are given, a feature must pass both. ``layers``, ``obsm``, ``obsp``,
          ``varm``, ``varp`` and every non-``biotapy`` ``uns`` key are dropped
          because they described the old features.

      Raises
      ------
      ValueError
          Neither threshold is given, ``min_prevalence`` is outside 0 to 1, or no
          feature passes.

      Notes
      -----
      R equivalent: ``phyloseq::filter_taxa``
      Guide: :doc:`/guide/filtering`

      In R: ``filter_taxa(physeq, function(x) sum(x > 0) >= p * length(x), prune = TRUE)``
      for ``min_prevalence=p``, and ``function(x) sum(x) >= n`` for ``min_total=n``.

      Examples
      --------
      >>> import biotapy as bt
      >>> bt.pp.filter_features(bt.datasets.toy(), min_prevalence=1.0).n_vars
      2
      """
      if min_prevalence is None and min_total is None:
          msg = "pass min_prevalence=, min_total= or both"
          raise ValueError(msg)
      X = as_csr(adata.X)
      keep = np.ones(adata.n_vars, dtype=bool)
      if min_prevalence is not None:
          if not 0 <= min_prevalence <= 1:
              msg = f"min_prevalence must be between 0 and 1, got {min_prevalence}"
              raise ValueError(msg)
          # Count stored non-zero values per column: a CSR matrix may also store explicit zeros.
          present = np.bincount(X.indices[X.data != 0], minlength=adata.n_vars)
          # Divide rather than multiply: 3 / 10 >= 0.3 holds, 3 >= 0.3 * 10 does not.
          keep &= present / adata.n_obs >= min_prevalence
      if min_total is not None:
          keep &= np.asarray(X.sum(axis=0)).ravel() >= min_total
      if not keep.any():
          msg = f"no feature passes min_prevalence={min_prevalence}, min_total={min_total}"
          raise ValueError(msg)
      out = feature_subset(adata, np.flatnonzero(keep))
      add_provenance(out, "pp.filter_features", min_prevalence=min_prevalence, min_total=min_total)
      return out


  def filter_samples(adata: AnnData, min_depth: float) -> AnnData:
      """Keep samples with at least ``min_depth`` reads.

      Parameters
      ----------
      adata
          Samples x features.
      min_depth
          Keep samples whose total over all features is at least this.

      Returns
      -------
      AnnData
          Same type as ``adata`` with the kept samples in their original order.
          Every slot is kept and subset by AnnData indexing, so ``obsp`` distances
          stay valid; features that are now all-zero are kept.

      Raises
      ------
      ValueError
          No sample has ``min_depth`` reads.

      Notes
      -----
      R equivalent: none
      Guide: :doc:`/guide/filtering`

      In R: ``prune_samples(sample_sums(physeq) >= min_depth, physeq)``.

      Examples
      --------
      >>> import biotapy as bt
      >>> bt.pp.filter_samples(bt.datasets.toy(), 70).n_obs
      4
      """
      depth = np.asarray(as_csr(adata.X).sum(axis=1)).ravel()
      keep = np.flatnonzero(depth >= min_depth)
      if keep.size == 0:
          msg = f"no sample has min_depth={min_depth} reads; the deepest has {depth.max():g}"
          raise ValueError(msg)
      out = adata[keep].copy()
      add_provenance(out, "pp.filter_samples", min_depth=min_depth)
      return out
  ```
  - `src/biotapy/pp/__init__.py` gains `from ._filter import filter_features, filter_samples`.
  - Its `__all__` becomes `["filter_features", "filter_samples", "relative", "tax_glom"]`.
  - `np.bincount` over the column indices of stored non-zero values counts
    presence. A sparse `X != 0` also works at run time, but scipy-stubs types
    it as `bool`, so mypy rejects it.
- [x] **Step 6: Docs.**
  - Create `docs/guide/filtering.md`:
    ````markdown
    # Filtering and rarefaction

    Filters keep a subset of features or samples; rarefaction subsamples every sample to the same
    depth. Each returns a new object and leaves its input alone.

    ## Features: `filter_features`

    `bt.pp.filter_features` keeps features that are present in enough samples (`min_prevalence`, a
    fraction from 0 to 1) and have enough reads in total (`min_total`). Both thresholds are
    inclusive. Give either or both; with both, a feature must pass both:

    ```python
    import biotapy as bt

    tdata = bt.datasets.toy()
    out = bt.pp.filter_features(tdata, min_prevalence=0.5, min_total=20)
    ```

    Present means a non-zero count; a zero stored explicitly in the sparse matrix counts as absent.
    In phyloseq the same filters are
    `filter_taxa(physeq, function(x) sum(x > 0) >= 0.5 * length(x), prune = TRUE)` and
    `filter_taxa(physeq, function(x) sum(x) >= 20, prune = TRUE)`.

    Filtering changes the feature set, so everything computed from the old one is dropped: every
    entry in `layers`, `obsm` and `obsp`. A TreeData keeps the subtree of the kept features, with
    branch lengths unchanged.

    ## Samples: `filter_samples`

    `bt.pp.filter_samples(tdata, min_depth)` keeps samples with at least `min_depth` reads (again
    inclusive). It only removes rows, so every slot is kept and subset, and distances in `obsp` stay
    valid for the samples that remain. Features that become all-zero are kept; follow with
    `filter_features` to drop them.
    ````
  - Add `filtering` to the `docs/guide/index.md` toctree, after `aggregation`.
  - Add `pp.filter_features` and `pp.filter_samples` to the Preprocessing
    block of `docs/api.md`, alphabetically before `pp.relative`.
- [x] **Step 7: Knowledge.**
  - `.knowledge/modules/pp.md`:
    - In the Responsibility, "Does NOT own filtering or rarefaction
      (`pp.filter_features`, `pp.filter_samples`, `pp.rarefy` - Slice 1C, not
      yet written)" becomes "Owns filtering (`filter_features`,
      `filter_samples`); does NOT own rarefaction (`pp.rarefy`, Task 1.14)".
    - Add Entry points:
      - `_filter.py:filter_features`: prevalence and total thresholds,
        through `feature_subset`;
      - `_filter.py:filter_samples`: a depth threshold, through AnnData
        indexing, keeping every slot.
    - Add Invariants:
      - thresholds are inclusive;
      - prevalence counts stored non-zero values and divides by `n_obs`,
        rather than multiplying the threshold;
      - nothing passing is a `ValueError`, never an empty object.
  - Add a log line. Tick 1.13 here.
- [x] **Step 8: Run, gate and commit.**
  - Tests: `uv run --group test pytest tests/pp tests/core -q`.
  - Golden: `uv run --group test pytest -m golden tests/pp/test_filter_golden.py -q`
    (2 passed; the GlobalPatterns download is cached by pooch).
  - Gates: `uvx prek run --all-files`, `uv run --group test pytest`, and
    `uv run --group doc sphinx-build -W -b html docs docs/_build/html`.
  - Commit `feat(pp): add filter_features and filter_samples`.

### Task 1.14: `pp.rarefy`

**Files:**
- Create:
  - `src/biotapy/pp/_rarefy.py`
  - `tests/pp/test_rarefy.py`
  - `tests/pp/test_rarefy_golden.py`
- Modify:
  - `src/biotapy/pp/__init__.py`
  - `pyproject.toml` (`untyped_calls_exclude`)
  - `docs/guide/filtering.md`, `docs/api.md`
  - `.knowledge/modules/pp.md`, `.knowledge/modules/core.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:**
  - `require_counts`, `as_csr`, `as_generator`, `feature_subset`,
    `add_provenance`, `warn_user`;
  - `skbio.stats.subsample_counts`;
  - `tests/golden/global_patterns/rarefy.csv.gz`.
- **Produces:** `bt.pp.rarefy(adata: AnnData, *, depth: int | None = None, seed: int | np.random.Generator | None = None) -> AnnData`.
  - `depth` defaults to the smallest non-zero sample depth.
  - Samples below `depth` (strict `<`, as phyloseq) are dropped with one
    `UserWarning` naming up to 5 of them.
  - Draws are without replacement, per row.
  - Features left all-zero are dropped through `feature_subset`.
  - `X` becomes `int64` counts, and provenance records `pp.rarefy` with `depth`.

- [x] **Step 1: Failing tests.**
  - `tests/pp/test_rarefy.py`:
    ```python
    import json
    import warnings

    import anndata as ad
    import numpy as np
    import pandas as pd
    import pytest
    import scipy.sparse as sp
    from hypothesis import given
    from hypothesis import strategies as st
    from hypothesis.extra.numpy import arrays

    import biotapy as bt
    from biotapy._core import get_tree


    def _adata(dense) -> ad.AnnData:
        return ad.AnnData(
            X=sp.csr_matrix(dense),
            obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
            var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
        )


    def _depths(adata) -> list[int]:
        return np.asarray(adata.X.sum(axis=1)).ravel().astype(int).tolist()


    # toy() sample depths: s1..s6 = 68, 66, 77, 76, 82, 73.


    def test_default_depth_is_the_smallest_sample():
        out = bt.pp.rarefy(bt.datasets.toy(), seed=0)
        assert out.n_obs == 6 and _depths(out) == [66] * 6


    def test_sample_at_exactly_depth_is_kept():
        with pytest.warns(UserWarning, match=r"dropped 2 sample\(s\).*\['s1', 's2'\]"):
            out = bt.pp.rarefy(bt.datasets.toy(), depth=73, seed=0)
        assert list(out.obs_names) == ["s3", "s4", "s5", "s6"] and _depths(out) == [73] * 4


    def test_no_count_exceeds_the_original():
        tdata = bt.datasets.toy()
        out = bt.pp.rarefy(tdata, depth=50, seed=1)
        assert (out.X.toarray() <= tdata[:, out.var_names].X.toarray()).all()


    def test_same_seed_same_result_and_generator_accepted():
        first = bt.pp.rarefy(bt.datasets.toy(), depth=40, seed=7)
        again = bt.pp.rarefy(bt.datasets.toy(), depth=40, seed=np.random.default_rng(7))
        assert (first.X != again.X).nnz == 0


    def test_features_left_all_zero_are_dropped():
        dense = np.array([[5, 1, 0], [6, 0, 1]])
        out = bt.pp.rarefy(_adata(dense), depth=5, seed=0)
        assert out.X.toarray().sum(axis=0).min() > 0


    def test_toy_tree_is_pruned_to_kept_features():
        out = bt.pp.rarefy(bt.datasets.toy(), depth=5, seed=0)
        tree = get_tree(out)
        assert {n for n in tree.nodes if tree.out_degree(n) == 0} == set(out.var_names)


    def test_all_zero_sample_is_dropped_by_default():
        with pytest.warns(UserWarning, match=r"\['s0'\]"):
            out = bt.pp.rarefy(_adata(np.array([[0, 0], [3, 4], [5, 5]])), seed=0)
        assert list(out.obs_names) == ["s1", "s2"] and _depths(out) == [7, 7]


    def test_single_sample():
        assert _depths(bt.pp.rarefy(_adata(np.array([[3, 9, 2]])), depth=10, seed=0)) == [10]


    def test_drops_derived_slots_keeps_counts_and_records_provenance():
        out = bt.pp.rarefy(bt.pp.relative(bt.datasets.toy()), depth=60, seed=0)
        assert "relative" not in out.layers and out.uns["biotapy"]["x_kind"] == "counts"
        entry = json.loads(out.uns["biotapy"]["provenance"][-1])
        assert entry["step"] == "pp.rarefy" and entry["params"] == {"depth": 60}


    def test_input_unchanged(assert_unchanged):
        tdata = bt.datasets.toy()
        before = tdata.copy()
        bt.pp.rarefy(tdata, depth=60, seed=0)
        assert_unchanged(before, tdata)


    def test_rejects_non_counts():
        tdata = bt.datasets.toy()
        tdata.uns["biotapy"]["x_kind"] = "relative"
        with pytest.raises(ValueError, match="pp.rarefy needs raw counts"):
            bt.pp.rarefy(tdata)


    def test_depth_below_one_raises():
        with pytest.raises(ValueError, match="depth must be at least 1"):
            bt.pp.rarefy(bt.datasets.toy(), depth=0)


    def test_depth_above_every_sample_raises():
        with pytest.raises(ValueError, match="depth=1000"):
            bt.pp.rarefy(bt.datasets.toy(), depth=1000)


    @given(
        arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 30)), st.integers(1, 40)
    )
    def test_kept_rows_sum_to_depth_and_never_exceed_the_original(dense, depth):
        deep = dense.sum(axis=1) >= depth
        if not deep.any():
            with pytest.raises(ValueError, match="depth="):
                bt.pp.rarefy(_adata(dense), depth=depth, seed=0)
            return
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # shallow samples are dropped with a warning
            out = bt.pp.rarefy(_adata(dense), depth=depth, seed=0)
        assert list(out.obs_names) == [f"s{i}" for i in np.flatnonzero(deep)]
        assert _depths(out) == [depth] * int(deep.sum())
        kept = dense[deep][:, [int(name[1:]) for name in out.var_names]]
        assert (out.X.toarray() <= kept).all()
    ```
  - `tests/pp/test_rarefy_golden.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pytest

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    def test_rarefy_drops_the_same_samples_as_phyloseq():
        # R and NumPy generators differ, so only invariants are compared (contracts/r-golden-parity).
        golden = pd.read_csv(GOLDEN / "rarefy.csv.gz", dtype={"sample_id": str})
        tdata = bt.datasets.global_patterns()
        with pytest.warns(UserWarning, match="TRRsed1"):
            out = bt.pp.rarefy(tdata, depth=100_000, seed=0)
        assert list(out.obs_names) == golden["sample_id"].tolist()
        np.testing.assert_array_equal(np.asarray(out.X.sum(axis=1)).ravel(), golden["sample_sum"].to_numpy())
        assert (out.X.toarray() <= tdata[out.obs_names, out.var_names].X.toarray()).all()
    ```
- [x] **Step 2: Run, expect failure.** Run `uv run --group test pytest tests/pp/test_rarefy.py -q`.
  It fails with `AttributeError: module 'biotapy.pp' has no attribute 'rarefy'`.
- [x] **Step 3: Implement.** `src/biotapy/pp/_rarefy.py`:
  ```python
  """Rarefaction: subsample every sample to the same depth."""

  import numpy as np
  import numpy.typing as npt
  import scipy.sparse as sp
  from anndata import AnnData
  from skbio.stats import subsample_counts

  from biotapy._core import add_provenance, as_csr, as_generator, feature_subset, require_counts, warn_user


  def rarefy(adata: AnnData, *, depth: int | None = None, seed: int | np.random.Generator | None = None) -> AnnData:
      """Subsample every sample to ``depth`` reads, without replacement.

      Parameters
      ----------
      adata
          Samples x features with raw counts in ``X``.
      depth
          Reads to keep per sample. Defaults to the smallest non-zero sample depth.
      seed
          Seed or generator for the subsampling.

      Returns
      -------
      AnnData
          Same type as ``adata``. Samples with fewer than ``depth`` reads are dropped,
          with one warning naming them. Every kept sample sums to ``depth``. Features
          left all-zero are dropped, a TreeData keeps the remaining features' subtree,
          and ``layers``, ``obsm``, ``obsp``, ``varm``, ``varp`` and every
          non-``biotapy`` ``uns`` key are dropped.

      Raises
      ------
      ValueError
          ``X`` does not hold counts, ``depth`` is below 1, or no sample has ``depth`` reads.

      Warns
      -----
      UserWarning
          Samples with fewer than ``depth`` reads were dropped.

      Notes
      -----
      R equivalent: ``phyloseq::rarefy_even_depth``
      Guide: :doc:`/guide/filtering`

      Draws are without replacement (scikit-bio ``subsample_counts``), so no count
      exceeds the original; phyloseq defaults to ``replace = TRUE``. The default depth
      skips all-zero samples, where phyloseq's ``min(sample_sums(physeq))`` would be 0.
      R and NumPy random generators differ, so the counts never match phyloseq's.

      Examples
      --------
      >>> import biotapy as bt
      >>> out = bt.pp.rarefy(bt.datasets.toy(), depth=60, seed=0)
      >>> out.X.sum(axis=1).A1.tolist()
      [60, 60, 60, 60, 60, 60]
      """
      require_counts(adata, func="pp.rarefy")
      X = as_csr(adata.X)
      sums = np.asarray(X.sum(axis=1)).ravel()
      depth = _smallest_nonzero(sums) if depth is None else depth
      if depth < 1:
          msg = f"depth must be at least 1, got {depth}"
          raise ValueError(msg)
      # phyloseq drops samples with sample_sums(physeq) < sample.size; a sample at exactly depth stays.
      kept = np.flatnonzero(sums >= depth)
      if kept.size == 0:
          msg = f"no sample has depth={depth} reads; the deepest has {sums.max():g}"
          raise ValueError(msg)
      if kept.size < adata.n_obs:
          dropped = adata.obs_names[sums < depth].tolist()
          warn_user(f"pp.rarefy dropped {len(dropped)} sample(s) with fewer than depth={depth} reads: {dropped[:5]}")
      counts = _subsample_rows(X[kept], depth, as_generator(seed))
      present = np.unique(counts.indices).astype(np.intp)
      out = feature_subset(adata[kept], present)
      out.X = counts[:, present]
      add_provenance(out, "pp.rarefy", depth=depth)
      return out


  def _smallest_nonzero(sums: npt.NDArray[np.float64]) -> int:
      # phyloseq's default is min(sample_sums(physeq)); skipping empty samples drops them instead of
      # rarefying everything to 0. With every sample empty this is 1, and the caller's no-sample check raises.
      nonzero = sums[sums > 0]
      return int(nonzero.min()) if nonzero.size else 1


  def _subsample_rows(X: sp.csr_matrix, depth: int, rng: np.random.Generator) -> sp.csr_matrix:
      # Only the stored values are drawn from; zeros cannot be drawn, so positions carry over.
      data = np.concatenate(
          [
              subsample_counts(X.data[start:end].astype(np.int64), depth, replace=False, seed=rng)
              for start, end in zip(X.indptr[:-1], X.indptr[1:], strict=True)
          ]
      )
      out = sp.csr_matrix((data, X.indices.copy(), X.indptr.copy()), shape=X.shape)
      out.eliminate_zeros()
      return out
  ```
  - `pp/__init__.py` gains `from ._rarefy import rarefy`, and `"rarefy"` joins `__all__`.
  - Why the code looks as it does:
    - Only a row's stored values are subsampled, and positions carry over,
      because a zero can never be drawn.
    - `indices`/`indptr` are copied and `data` is new. `as_csr` may share
      buffers with the input, so nothing is written back (core invariant).
    - mypy flags `Call to untyped function "subsample_counts" in typed context`,
      so in `pyproject.toml` replace the `untyped_calls_exclude` comment and line with:
      ```toml
      # biom-format ships no type annotations at all, nor does scikit-bio's subsample_counts; exempt
      # only calls into them (not our own code) from strict's disallow_untyped_calls
      # (mypy/checkexpr.py matches by callee fullname).
      untyped_calls_exclude = [ "biom", "skbio.stats._subsample" ]
      ```
- [x] **Step 4: Docs.**
  - Append to `docs/guide/filtering.md`:
    ````markdown

    ## Rarefaction: `rarefy`

    `bt.pp.rarefy` subsamples every sample to `depth` reads without replacement, so no count grows
    and every kept sample sums to exactly `depth`:

    ```python
    out = bt.pp.rarefy(tdata, depth=60, seed=0)
    ```

    - `depth` defaults to the smallest non-zero sample depth. phyloseq's `rarefy_even_depth` uses
      the smallest depth, zero included.
    - Samples with fewer than `depth` reads are dropped, with one warning naming them; a sample with
      exactly `depth` reads is kept.
    - Features left all-zero are dropped, as with phyloseq's `trimOTUs = TRUE`.
    - `X` must hold raw counts (`uns["biotapy"]["x_kind"] == "counts"`).
    - phyloseq samples with replacement by default (`replace = TRUE`). biotapy always samples
      without, like `replace = FALSE`.
    - The same `seed` gives the same result. R and NumPy use different random generators, so the
      counts never match phyloseq's draw for draw.
    ````
  - Add `pp.rarefy` to the Preprocessing block of `docs/api.md`.
- [x] **Step 5: Knowledge.**
  - `.knowledge/modules/pp.md`:
    - the Responsibility now owns rarefaction too;
    - add the `_rarefy.py:rarefy` entry point;
    - add Invariants: without replacement, the default depth, the strict `<`
      drop with a warning, all-zero features dropped, and `X` becomes `int64`.
  - `.knowledge/modules/core.md`: in the `as_csr` invariant, "matters for the
    future `pp.rarefy`" becomes "`pp.rarefy` copies `indices`/`indptr` and
    builds new `data`".
  - Add a log line. Tick 1.14 here.
- [x] **Step 6: Run, gate and commit.**
  - Tests: `uv run --group test pytest tests/pp -q`.
  - Golden: `uv run --group test pytest -m golden tests/pp/test_rarefy_golden.py -q`
    (1 passed, about 4 s).
  - The three gates from Task 1.13 Step 8.
  - Commit `feat(pp): add rarefy`.

### Task 1.15b: `datasets.esophagus()`

**Files:**
- Modify:
  - `src/biotapy/datasets/_remote.py`, `src/biotapy/datasets/__init__.py`
  - `tests/datasets/test_remote.py`
  - `docs/api.md`, `docs/guide/datasets.md`, `docs/guide/reading_data.md`
  - `.knowledge/modules/datasets.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:** `read_phyloseq`, `_fetch`.
- **Produces:** `bt.datasets.esophagus() -> TreeData`: 3 samples (B, C, D) x
  58 OTUs, with a tree and no taxonomy or sample data. 1.16's UniFrac golden
  test uses it.

- [x] **Step 1: Failing tests.** In `tests/datasets/test_remote.py`:
  - Extend `test_loaders_read_the_fetched_file`: after the enterotype
    assertion, add `assert bt.datasets.esophagus().shape == (6, 8)`. The
    `fetched` expectation becomes
    `["GlobalPatterns.RData", "enterotype.RData", "esophagus.RData"]`.
  - Append:
    ```python
    @pytest.mark.network
    def test_esophagus_downloads_with_a_tree_and_no_metadata():
        tdata = bt.datasets.esophagus()
        tree = tdata.vart["phylo"]
        assert tdata.shape == (3, 58) and list(tdata.obs_names) == ["B", "C", "D"]
        assert sum(1 for n in tree.nodes if tree.out_degree(n) == 0) == 58
        assert tdata.obs.columns.empty and tdata.var.columns.empty and tdata.uns["biotapy"]["x_kind"] == "counts"
    ```
- [x] **Step 2: Run, expect failure.** Run `uv run --group test pytest tests/datasets -q`.
  It fails with `AttributeError: module 'biotapy.datasets' has no attribute 'esophagus'`.
- [x] **Step 3: Implement.**
  - In `_remote.py`, `_REGISTRY` gains
    `"esophagus.RData": "sha256:0b06d9c35f2e694c34461308af149eb54419453fcb98763de20ab61980b87e46",`.
    It is at the same pinned commit; the file is 1,840 B.
  - Append:
    ```python
    def esophagus() -> TreeData:
        """esophagus: 3 esophageal biopsies, 58 OTUs, with a tree; no taxonomy or sample data.

        Downloaded once (2 kB) from phyloseq's repository and cached.

        Returns
        -------
        TreeData
            Counts in ``X`` (samples ``B``, ``C``, ``D``) and the tree in ``vart['phylo']``;
            ``obs`` and ``var`` have no columns.

        Notes
        -----
        R equivalent: ``utils::data``
        Guide: :doc:`/guide/datasets`

        In R: ``data(esophagus, package = "phyloseq")``.

        References
        ----------
        Pei Z et al. (2004) Bacterial biota in the human distal esophagus. PNAS 101:4250-4255.

        Examples
        --------
        >>> import biotapy as bt
        >>> bt.datasets.esophagus().shape  # doctest: +SKIP
        (3, 58)
        """
        return read_phyloseq(_fetch("esophagus.RData"))
    ```
  - `datasets/__init__.py`: `from ._remote import enterotype, esophagus, global_patterns`
    and `__all__ = ["enterotype", "esophagus", "global_patterns", "toy"]`.
- [x] **Step 4: Docs.**
  - In `docs/api.md`'s Datasets block, add `datasets.esophagus` after `datasets.enterotype`.
  - `docs/guide/datasets.md`, the example-datasets page since PR #7:
    - The intro "ships three example datasets" becomes "ships four example datasets".
    - The heading `` ## `global_patterns` and `enterotype` `` becomes
      `` ## `global_patterns`, `enterotype` and `esophagus` ``.
    - Its first sentence becomes "`bt.datasets.global_patterns()`,
      `bt.datasets.enterotype()` and `bt.datasets.esophagus()` are three
      well-known phyloseq example datasets, read through `bt.io.read_phyloseq`:".
    - After the enterotype bullet, add:
      ```markdown
      - **esophagus**: 3 esophageal biopsies (samples `B`, `C`, `D`), 58 OTUs, with a tree but no
        taxonomy or sample data; counts in `X`. biotapy's UniFrac golden tests use it.
      ```
    - In "Caching with pooch", "Each of `global_patterns` and `enterotype`"
      becomes "Each of `global_patterns`, `enterotype` and `esophagus`".
    - In "Licensing", "`global_patterns` and `enterotype` download data"
      becomes "`global_patterns`, `enterotype` and `esophagus` download data".
  - `docs/guide/reading_data.md`, under "Example datasets":
    "`bt.datasets.global_patterns()` and `bt.datasets.enterotype()`, two
    well-known phyloseq example datasets" becomes
    "`bt.datasets.global_patterns()`, `bt.datasets.enterotype()` and
    `bt.datasets.esophagus()`, three well-known phyloseq example datasets".
- [x] **Step 5: Knowledge.**
  - `.knowledge/modules/datasets.md`:
    - the Responsibility and Entry points gain `_remote.py:esophagus`;
    - the Invariants' "both are downloaded" becomes "all three are downloaded";
    - the Dependencies' pooch line names `esophagus.RData`.
  - Add a log line. Tick 1.15b here.
- [x] **Step 6: Run, gate and commit.**
  - Run `uv run --group test pytest tests/datasets -q`.
  - Run once `BIOTAPY_DATA_DIR=/tmp/claude-1000/pooch uv run --group test pytest -m network tests/datasets -q`
    and expect 4 passed.
  - The three gates.
  - Commit `feat(datasets): add esophagus via pooch`.

### Task 1.15c: `_core.get_skbio_tree`

**Files:**
- Modify:
  - `src/biotapy/_core/_tree.py`, `src/biotapy/_core/__init__.py`
  - `tests/core/test_tree.py`
  - `.knowledge/contracts/tree-access.md`, `.knowledge/modules/core.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:** `get_tree`, `PHYLO_KEY`, `skbio.TreeNode`.
- **Produces:** `_core.get_skbio_tree(adata: AnnData) -> skbio.TreeNode`.
  - A plain AnnData raises `TypeError`, "needs a TreeData with a tree in vart['phylo']".
  - No `vart['phylo']` raises `KeyError`, from `get_tree`.
  - A root with more than two children gets a zero-length split.
  - Used by `tl.alpha` (faith_pd) and `tl.unifrac`. Only `_tree.py` may import
    networkx (R4.5), so this lives in `_core` whatever the caller count.

- [x] **Step 1: Failing tests.** In `tests/core/test_tree.py`:
  - Add `import networkx as nx`, `from skbio import TreeNode` and `import biotapy as bt`
    to the import block.
  - Add `get_skbio_tree` to the `from biotapy._core import (...)` list.
  - Append:
    ```python
    def test_get_skbio_tree_keeps_names_and_lengths():
        tree = get_skbio_tree(_make(tree_from_edges([("r", "a", 1.0), ("r", "b", 2.0)])))
        assert isinstance(tree, TreeNode) and tree.name == "r"
        assert {tip.name: tip.length for tip in tree.tips()} == {"a": 1.0, "b": 2.0}


    def test_get_skbio_tree_splits_a_multifurcating_root_without_changing_paths():
        tdata = bt.datasets.toy()  # the root has three children
        tree = get_skbio_tree(tdata)
        assert len(tree.children) == 2
        expected = nx.shortest_path_length(get_tree(tdata), "root", weight="length")
        assert {tip.name: tip.distance(tree) for tip in tree.tips()} == pytest.approx(
            {f"f{i}": expected[f"f{i}"] for i in range(1, 9)}
        )


    def test_get_skbio_tree_keeps_nan_lengths():
        assert math.isnan(get_skbio_tree(_make(tree_from_newick("(a,b:2);"))).find("a").length)


    def test_get_skbio_tree_without_phylogeny_names_the_key():
        with pytest.raises(KeyError, match="phylo"):
            get_skbio_tree(_make(None))
        with pytest.raises(TypeError, match="needs a TreeData"):
            get_skbio_tree(bt.datasets.toy().to_adata())
    ```
- [x] **Step 2: Run, expect failure.** Run `uv run --group test pytest tests/core/test_tree.py -q`.
  It fails with `ImportError: cannot import name 'get_skbio_tree' from 'biotapy._core'`.
- [x] **Step 3: Implement.**
  - In `_core/_tree.py`, add `from anndata import AnnData` next to the other
    third-party imports, and this function before `tree_from_newick`:
    ```python
    def get_skbio_tree(adata: AnnData) -> TreeNode:
        """The phylogeny in ``vart['phylo']`` as a scikit-bio ``TreeNode``, rooted where it is drawn.

        scikit-bio's Faith PD and UniFrac accept a root with at most two children. A
        root with more (an unrooted Newick tree, or the toy tree) keeps its first child
        and gets the others under one new zero-length node, which changes no
        root-to-tip path length. Missing (NaN) branch lengths stay NaN; scikit-bio
        counts them as zero. A plain AnnData raises ``TypeError``; a TreeData without
        ``vart['phylo']`` raises ``KeyError``.
        """
        if not isinstance(adata, TreeData):
            msg = f"needs a TreeData with a tree in vart[{PHYLO_KEY!r}], got {type(adata).__name__}"
            raise TypeError(msg)
        tree = get_tree(adata)
        root = next(node for node in tree if tree.in_degree(node) == 0)
        nodes = {root: TreeNode(name=root)}
        for parent, child in nx.bfs_edges(tree, root):
            nodes[child] = TreeNode(name=child, length=tree.edges[parent, child]["length"])
            nodes[parent].append(nodes[child])
        top = nodes[root]
        if len(top.children) > 2:
            split = TreeNode(length=0.0)
            split.extend(top.children[1:])
            top.append(split)
        return top
    ```
  - Export `get_skbio_tree` from `_core/__init__.py`, both the import and `__all__`.
  - Checked in the prototype:
    - toy's tip distances to the root are unchanged;
    - GlobalPatterns (19,216 tips) converts in 0.1 s;
    - `TreeNode.append` and `extend` do not copy.
- [x] **Step 4: Knowledge.**
  - `contracts/tree-access.md` statement 2 gains "converting the phylogeny to
    a scikit-bio `TreeNode` (`get_skbio_tree`)". Add a Gotcha: scikit-bio's
    Faith PD and UniFrac need a root with at most two children, so
    `get_skbio_tree` splits a wider root with a zero-length node, and phyloseq
    instead roots an unrooted tree at a random tip.
  - `.knowledge/modules/core.md`:
    - add the Entry point `_tree.py:get_skbio_tree`;
    - Dependencies: scikit-bio is also used for `TreeNode` conversion.
  - Add a log line. Tick 1.15c here.
  - There is no docs step: `_core` is private and not in `docs/api.md`. Users
    meet the behaviour through the `tl.alpha` and `tl.unifrac` docstrings.
- [x] **Step 5: Run, gate and commit.**
  - Run `uv run --group test pytest tests/core -q`, then `uvx prek run --all-files`
    and `uv run --group test pytest`.
  - Commit `feat(core): convert the phylogeny to a scikit-bio TreeNode`.

### Task 1.15: `tl.alpha`

**Files:**
- Create:
  - `src/biotapy/tl/__init__.py`, `src/biotapy/tl/_alpha.py`
  - `tests/tl/test_alpha.py`, `tests/tl/test_alpha_golden.py`
  - `docs/guide/diversity.md`
- Modify:
  - `src/biotapy/__init__.py`
  - `docs/api.md`, `docs/guide/index.md`
  - `.knowledge/contracts/data-model-slots.md`, `.knowledge/modules/pp.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:**
  - `as_csr`, `require_counts`, `get_skbio_tree` (1.15c);
  - `skbio.diversity.alpha_diversity`;
  - `alpha.csv.gz` and `alpha_faith_pd.csv.gz`.
- **Produces:**
  - `bt.tl.alpha(adata: AnnData, *, metrics: Sequence[AlphaMetric] = ("observed_features", "shannon", "simpson", "chao1"), inplace: bool = False) -> pd.DataFrame | None`.
  - `AlphaMetric = Literal["observed_features", "shannon", "simpson", "chao1", "faith_pd"]`.
  - The package `bt.tl`.

- [x] **Step 1: Failing tests.**
  - `tests/tl/test_alpha.py`:
    ```python
    import tracemalloc

    import anndata as ad
    import numpy as np
    import pandas as pd
    import pytest
    import scipy.sparse as sp
    from hypothesis import given
    from hypothesis import strategies as st
    from hypothesis.extra.numpy import arrays

    import biotapy as bt

    ALL = ["observed_features", "shannon", "simpson", "chao1", "faith_pd"]


    def _adata(dense) -> ad.AnnData:
        return ad.AnnData(
            X=sp.csr_matrix(dense),
            obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
            var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
        )


    def test_default_metrics_by_hand():
        out = bt.tl.alpha(_adata(np.array([[1, 2, 3, 0], [4, 4, 0, 0]])))
        assert list(out.columns) == ["observed_features", "shannon", "simpson", "chao1"]
        p = np.array([1, 2, 3]) / 6
        np.testing.assert_allclose(out.loc["s0"].to_numpy(), [3, -(p * np.log(p)).sum(), 1 - (p**2).sum(), 3], rtol=1e-12)
        np.testing.assert_allclose(out.loc["s1"].to_numpy(), [2, np.log(2), 0.5, 2], rtol=1e-12)


    def test_chao1_is_bias_corrected():
        # 2 singletons, 1 doubleton: 4 + 2 * 1 / (2 * (1 + 1)) = 4.5
        assert bt.tl.alpha(_adata(np.array([[1, 1, 2, 5]])), metrics=["chao1"]).loc["s0", "chao1"] == 4.5


    def test_faith_pd_sums_branches_to_the_root():
        # s1 lacks f5 and f8: every branch except n5-f5 (0.03) and n3-f8 (0.12); toy's branches total 1.34.
        out = bt.tl.alpha(bt.datasets.toy(), metrics=["faith_pd"])
        assert out.loc["s1", "faith_pd"] == pytest.approx(1.34 - 0.03 - 0.12)


    def test_all_zero_sample_is_nan_or_zero_never_raises():
        tdata = bt.datasets.toy()
        dense = tdata.X.toarray()
        dense[0] = 0
        tdata.X = sp.csr_matrix(dense)
        out = bt.tl.alpha(tdata, metrics=ALL)
        assert out.loc["s1", ["observed_features", "chao1", "faith_pd"]].tolist() == [0, 0, 0]
        assert out.loc["s1", ["shannon", "simpson"]].isna().all()


    def test_all_zero_feature_changes_nothing():
        dense = np.array([[1, 2, 0], [3, 1, 0]])
        pd.testing.assert_frame_equal(bt.tl.alpha(_adata(dense)), bt.tl.alpha(_adata(dense[:, :2])))


    def test_single_sample():
        assert bt.tl.alpha(_adata(np.array([[4, 0, 1]])), metrics=["observed_features"]).shape == (1, 1)


    def test_missing_rank_is_irrelevant():
        tdata = bt.datasets.toy()
        tdata.var["genus"] = np.nan
        assert bt.tl.alpha(tdata).shape == (6, 4)


    def test_inplace_writes_obs_and_returns_none():
        tdata = bt.datasets.toy()
        expected = bt.tl.alpha(tdata, metrics=["shannon", "faith_pd"])
        assert bt.tl.alpha(tdata, metrics=["shannon", "faith_pd"], inplace=True) is None
        np.testing.assert_array_equal(tdata.obs["alpha_shannon"], expected["shannon"])
        np.testing.assert_array_equal(tdata.obs["alpha_faith_pd"], expected["faith_pd"])


    def test_input_unchanged(assert_unchanged):
        tdata = bt.datasets.toy()
        before = tdata.copy()
        bt.tl.alpha(tdata, metrics=ALL)
        assert_unchanged(before, tdata)


    def test_densifies_in_bounded_chunks():
        # 2**19 + 1 features: each chunk of at most 2**20 values holds one sample (4 MiB);
        # densifying all 8 samples at once would take 32 MiB.
        n_obs, n_vars = 8, 2**19 + 1
        rng = np.random.default_rng(0)
        rows = np.repeat(np.arange(n_obs), 20)
        cols = rng.choice(n_vars, size=rows.size)
        X = sp.csr_matrix((rng.integers(1, 50, size=rows.size).astype(np.float64), (rows, cols)), shape=(n_obs, n_vars))
        wide = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"s{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"f{i}" for i in range(n_vars)]),
        )
        tracemalloc.start()
        out = bt.tl.alpha(wide)
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
        assert peak < n_obs * n_vars * 8 / 2
        pd.testing.assert_frame_equal(out, bt.tl.alpha(wide[:, np.unique(X.indices)].copy()))


    def test_string_metrics_raise():
        with pytest.raises(TypeError, match="not the string 'shannon'"):
            bt.tl.alpha(bt.datasets.toy(), metrics="shannon")


    def test_unknown_or_empty_metrics_raise():
        with pytest.raises(ValueError, match="metrics must name"):
            bt.tl.alpha(bt.datasets.toy(), metrics=["pielou"])
        with pytest.raises(ValueError, match="metrics must name"):
            bt.tl.alpha(bt.datasets.toy(), metrics=[])


    def test_count_metrics_need_counts_others_do_not():
        rel = bt.datasets.toy()
        rel.uns["biotapy"]["x_kind"] = "relative"
        with pytest.raises(ValueError, match=r"tl.alpha with \['chao1'\] needs raw counts"):
            bt.tl.alpha(rel, metrics=["shannon", "chao1"])
        assert bt.tl.alpha(rel, metrics=["shannon", "simpson", "faith_pd"]).shape == (6, 3)


    def test_faith_pd_needs_a_tree():
        with pytest.raises(TypeError, match="needs a TreeData"):
            bt.tl.alpha(_adata(np.array([[1, 2]])), metrics=["faith_pd"])
        tdata = bt.datasets.toy()
        del tdata.vart["phylo"]
        with pytest.raises(KeyError, match="phylo"):
            bt.tl.alpha(tdata, metrics=["faith_pd"])


    @given(arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 8)), elements=st.integers(0, 20)))
    def test_observed_features_counts_nonzero_features(dense):
        out = bt.tl.alpha(_adata(dense), metrics=["observed_features", "simpson"])
        np.testing.assert_array_equal(out["observed_features"], (dense > 0).sum(axis=1))
        assert (out["observed_features"] <= dense.shape[1]).all()
        simpson = out["simpson"].dropna()
        assert ((simpson >= 0) & (simpson < 1)).all()
    ```
  - `tests/tl/test_alpha_golden.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pytest

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    def test_alpha_matches_phyloseq_estimate_richness():
        golden = pd.read_csv(GOLDEN / "alpha.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
        out = bt.tl.alpha(bt.datasets.global_patterns())
        assert list(out.index) == list(golden.index)
        for metric, column in [
            ("observed_features", "Observed"),
            ("shannon", "Shannon"),
            ("simpson", "Simpson"),
            ("chao1", "Chao1"),
        ]:
            np.testing.assert_allclose(out[metric], golden[column], rtol=1e-7, err_msg=metric)


    def test_faith_pd_matches_picante_pd_with_root():
        golden = pd.read_csv(GOLDEN / "alpha_faith_pd.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
        out = bt.tl.alpha(bt.datasets.global_patterns(), metrics=["faith_pd"])
        assert list(out.index) == list(golden.index)
        np.testing.assert_allclose(out["faith_pd"], golden["pd"], rtol=1e-7)
    ```
  - What the tests pin:
    - `test_densifies_in_bounded_chunks` asserts the R6.2 chunking through the
      public API. The tracemalloc peak was 8.1 MiB against 32 MiB for full
      densification, and the test takes 0.04 s.
    - It then checks that the chunked result equals the one-chunk result on
      the non-zero columns.
- [x] **Step 2: Run, expect failure.** Run `uv run --group test pytest tests/tl -q`.
  It fails with `AttributeError: module 'biotapy' has no attribute 'tl'`.
- [x] **Step 3: Implement.**
  - `src/biotapy/tl/_alpha.py`:
    ```python
    """Alpha diversity: one value per sample and metric."""

    import math
    from collections.abc import Sequence
    from typing import Any, Literal, get_args

    import pandas as pd
    from anndata import AnnData
    from skbio.diversity import alpha_diversity

    from biotapy._core import as_csr, get_skbio_tree, require_counts

    AlphaMetric = Literal["observed_features", "shannon", "simpson", "chao1", "faith_pd"]
    # phyloseq's estimate_richness: natural-log Shannon (vegan::diversity) and bias-corrected Chao1 (vegan::estimateR).
    _PARAMS: dict[str, dict[str, Any]] = {"shannon": {"base": math.e}, "chao1": {"bias_corrected": True}}
    # scikit-bio needs dense rows (rules.md R6.2); densify at most 2**20 values (8 MiB of float64) at a time.
    _CHUNK_VALUES = 2**20


    def alpha(
        adata: AnnData,
        *,
        metrics: Sequence[AlphaMetric] = ("observed_features", "shannon", "simpson", "chao1"),
        inplace: bool = False,
    ) -> pd.DataFrame | None:
        """Alpha diversity of every sample.

        Parameters
        ----------
        adata
            Samples x features. ``"faith_pd"`` needs a TreeData with ``vart['phylo']``.
        metrics
            Any of ``"observed_features"``, ``"shannon"`` (natural log), ``"simpson"``
            (Gini-Simpson, ``1 - sum(p**2)``), ``"chao1"`` (bias-corrected) and ``"faith_pd"``.
        inplace
            Write ``obs['alpha_<metric>']`` for every metric and return ``None``.

        Returns
        -------
        pandas.DataFrame or None
            One row per sample (index ``obs_names``) and one column per metric, in the
            order given. An all-zero sample gets 0 for ``observed_features``, ``chao1``
            and ``faith_pd`` and NaN for ``shannon`` and ``simpson``.

        Raises
        ------
        TypeError
            ``metrics`` is a string, or ``"faith_pd"`` is asked of an AnnData that is not a TreeData.
        ValueError
            ``metrics`` is empty or names an unknown metric, or ``"observed_features"``
            or ``"chao1"`` is asked of data that is not counts.
        KeyError
            ``"faith_pd"`` is asked of a TreeData without ``vart['phylo']``.

        Notes
        -----
        R equivalent: ``phyloseq::estimate_richness``, ``picante::pd``
        Guide: :doc:`/guide/diversity`

        scikit-bio needs dense input, so rows are densified in chunks of at most 2**20
        values (8 MiB of float64). Faith PD includes the root, as
        ``picante::pd(include.root = TRUE)``; a root with more than two children first
        gets a zero-length split, which changes no root-to-tip distance. For an
        all-zero sample phyloseq reports Shannon 0 and Simpson 1; biotapy returns NaN.

        Examples
        --------
        >>> import biotapy as bt
        >>> out = bt.tl.alpha(bt.datasets.toy(), metrics=["observed_features", "shannon"])
        >>> out.loc["s1"].round(3).tolist()
        [6.0, 1.361]
        """
        _check_metrics(adata, metrics)
        params = {metric: _PARAMS.get(metric, {}) for metric in metrics}
        if "faith_pd" in metrics:
            params["faith_pd"] = {"taxa": adata.var_names.tolist(), "tree": get_skbio_tree(adata)}
        X = as_csr(adata.X)
        step = max(1, _CHUNK_VALUES // max(adata.n_vars, 1))
        frames = []
        for start in range(0, adata.n_obs, step):
            dense, ids = X[start : start + step].toarray(), adata.obs_names[start : start + step].tolist()
            frames.append(pd.DataFrame({m: alpha_diversity(m, dense, ids=ids, **params[m]) for m in metrics}))
        result = pd.concat(frames)
        if not inplace:
            return result
        for metric in metrics:
            adata.obs[f"alpha_{metric}"] = result[metric].to_numpy()
        return None


    def _check_metrics(adata: AnnData, metrics: Sequence[str]) -> None:
        if isinstance(metrics, str):
            msg = f"metrics must be a list of metric names, not the string {metrics!r}"
            raise TypeError(msg)
        known = get_args(AlphaMetric)
        if not metrics or any(metric not in known for metric in metrics):
            msg = f"metrics must name one or more of {list(known)}, got {list(metrics)}"
            raise ValueError(msg)
        needs_counts = [metric for metric in metrics if metric in ("observed_features", "chao1")]
        if needs_counts:
            require_counts(adata, func=f"tl.alpha with {needs_counts}")
    ```
  - `src/biotapy/tl/__init__.py`: `from ._alpha import alpha` and
    `__all__ = ["alpha"]`. Imports only, R4.1.
  - `src/biotapy/__init__.py`: `from . import datasets, io, pp, tl` and
    `__all__ = ["__version__", "datasets", "io", "pp", "tl"]`.
  - import-linter already declares the `(tl) | (fn)` layer, so `pyproject.toml`
    does not change.
- [x] **Step 4: Docs.**
  - `docs/conf.py` needs no change. Its `intersphinx_mapping` has mapped
    pandas since eb7538e, and the `nitpicky` build needs that for
    `pandas.DataFrame` in the signature: without it `-W` failed in the
    prototype on 11ee8dd. Confirm the entry is still there.
  - `docs/api.md` gains a section after Preprocessing:
    ````markdown
    ## Tools

    ```{eval-rst}
    .. module:: biotapy.tl
    .. currentmodule:: biotapy

    .. autosummary::
        :toctree: generated

        tl.alpha
    ```
    ````
  - Create `docs/guide/diversity.md` and add `diversity` to the guide toctree after `filtering`:
    ````markdown
    # Diversity

    ## Alpha diversity: `tl.alpha`

    `bt.tl.alpha` computes one value per sample for each metric, through scikit-bio:

    | Metric | What it is | In R |
    |---|---|---|
    | `observed_features` | features with a non-zero count | `estimate_richness`: `Observed` |
    | `shannon` | Shannon entropy, natural log | `estimate_richness`: `Shannon` |
    | `simpson` | Gini-Simpson, `1 - sum(p**2)` | `estimate_richness`: `Simpson` |
    | `chao1` | bias-corrected Chao1 | `estimate_richness`: `Chao1` |
    | `faith_pd` | Faith's phylogenetic diversity, root included | `picante::pd(include.root = TRUE)` |

    ```python
    import biotapy as bt

    tdata = bt.datasets.toy()
    table = bt.tl.alpha(tdata)  # samples x metrics
    bt.tl.alpha(tdata, metrics=["shannon", "faith_pd"], inplace=True)  # obs["alpha_shannon"], obs["alpha_faith_pd"]
    ```

    - `observed_features` and `chao1` need raw counts, as `estimate_richness` does. `shannon`,
      `simpson` and `faith_pd` also run on relative abundances.
    - An all-zero sample gets 0 for `observed_features`, `chao1` and `faith_pd`, and `NaN` for
      `shannon` and `simpson`; phyloseq reports Shannon 0 and Simpson 1.
    - `faith_pd` needs a TreeData with a tree. A root with three or more children, common in a
      tree read from unrooted Newick, gets a zero-length split, which changes no root-to-tip
      distance.
    - scikit-bio needs dense input, so biotapy densifies at most 2**20 values (8 MiB) at a time.
    ````
- [x] **Step 5: Knowledge.**
  - `contracts/data-model-slots.md` Convention 2: "Functions that need raw
    counts (rarefy, chao1) call `_core.require_counts`" becomes "Functions
    that need raw counts (`pp.rarefy`, and `tl.alpha` for `observed_features`
    and `chao1`, which phyloseq's `estimate_richness` refuses on non-integers)
    call `_core.require_counts`".
  - `.knowledge/modules/pp.md`: "any diversity/ordination computation (`tl`,
    later phases)" becomes "(`tl`, Slice 1C)". The `tl` Module concept itself
    comes at Checkpoint C.
  - Add a log line. Tick 1.15 here.
- [x] **Step 6: Run, gate and commit.**
  - Tests: `uv run --group test pytest tests/tl -q`.
  - Golden: `uv run --group test pytest -m golden tests/tl/test_alpha_golden.py -q` (2 passed).
  - The three gates.
  - Commit `feat(tl): add alpha diversity`.

### Task 1.16: `tl.beta`, `tl.unifrac`

**Files:**
- Create:
  - `src/biotapy/tl/_beta.py`
  - `tests/tl/test_beta.py`, `tests/tl/test_beta_golden.py`
- Modify:
  - `src/biotapy/tl/__init__.py`
  - `docs/api.md`, `docs/guide/diversity.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:**
  - `as_csr`, `get_skbio_tree`, `TreeData`;
  - `skbio.diversity.beta_diversity`, `skbio.DistanceMatrix`;
  - `bt.datasets.esophagus` (1.15b);
  - `bt.pp.filter_features` (1.13), in a test;
  - the `beta_*` and `unifrac_*` golden files.
- **Produces:**
  - `bt.tl.beta(adata: AnnData, *, metric: BetaMetric = "braycurtis", inplace: bool = False) -> pd.DataFrame | None`,
    with `BetaMetric = Literal["braycurtis", "jaccard"]`.
  - `bt.tl.unifrac(tdata: TreeData, *, weighted: bool = False, normalized: bool = True, inplace: bool = False) -> pd.DataFrame | None`.
  - `inplace=True` writes `obsp["braycurtis" | "jaccard" | "unweighted_unifrac" | "weighted_unifrac"]`
    as a dense `ndarray`.

- [x] **Step 1: Failing tests.**
  - `tests/tl/test_beta.py`:
    ```python
    import anndata as ad
    import numpy as np
    import pandas as pd
    import pytest
    import scipy.sparse as sp
    from hypothesis import given
    from hypothesis import strategies as st
    from hypothesis.extra.numpy import arrays

    import biotapy as bt
    from biotapy._core import make_treedata, tree_from_edges


    def _adata(dense) -> ad.AnnData:
        return ad.AnnData(
            X=sp.csr_matrix(dense),
            obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
            var=pd.DataFrame(index=[f"f{i}" for i in range(dense.shape[1])]),
        )


    def test_braycurtis_by_hand():
        out = bt.tl.beta(_adata(np.array([[1, 2, 3], [3, 2, 1]])))
        assert out.loc["s0", "s1"] == pytest.approx(4 / 12)
        assert list(out.index) == list(out.columns) == ["s0", "s1"]


    def test_jaccard_is_presence_absence():
        out = bt.tl.beta(_adata(np.array([[1, 0, 30], [5, 5, 0]])), metric="jaccard")
        assert out.loc["s0", "s1"] == pytest.approx(2 / 3)


    def test_all_zero_sample_and_feature():
        out = bt.tl.beta(_adata(np.array([[0, 0, 0], [1, 2, 0], [0, 0, 0]])))
        assert out.loc["s0", "s1"] == 1.0 and np.isnan(out.loc["s0", "s2"])


    def test_single_sample():
        assert bt.tl.beta(_adata(np.array([[1, 2]]))).to_numpy().tolist() == [[0.0]]


    def test_missing_rank_is_irrelevant():
        tdata = bt.datasets.toy()
        tdata.var["genus"] = np.nan
        assert bt.tl.beta(tdata).shape == (6, 6)


    def test_beta_inplace_writes_obsp(assert_unchanged):
        tdata = bt.datasets.toy()
        expected = bt.tl.beta(tdata, metric="jaccard")
        assert bt.tl.beta(tdata, metric="jaccard", inplace=True) is None
        np.testing.assert_array_equal(tdata.obsp["jaccard"], expected.to_numpy())


    def test_beta_input_unchanged(assert_unchanged):
        tdata = bt.datasets.toy()
        before = tdata.copy()
        bt.tl.beta(tdata)
        assert_unchanged(before, tdata)


    def test_unknown_metric_raises():
        with pytest.raises(ValueError, match="metric must be one of"):
            bt.tl.beta(bt.datasets.toy(), metric="euclidean")


    @given(
        arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 20)),
        st.sampled_from(["braycurtis", "jaccard"]),
    )
    def test_beta_is_symmetric_with_zero_diagonal(dense, metric):
        out = bt.tl.beta(_adata(dense), metric=metric).to_numpy()
        np.testing.assert_array_equal(out, out.T)
        np.testing.assert_array_equal(np.diag(out), 0.0)


    def test_unifrac_by_hand():
        # root -> a (1), root -> b (3); s0 holds a, s1 holds b: unweighted 4/4, weighted normalized 1.
        tdata = make_treedata(
            np.array([[2, 0], [0, 5]]),
            obs=pd.DataFrame(index=["s0", "s1"]),
            var=pd.DataFrame(index=["a", "b"]),
            tree=tree_from_edges([("r", "a", 1.0), ("r", "b", 3.0)]),
            x_kind="counts",
            source="test",
        )
        assert bt.tl.unifrac(tdata).loc["s0", "s1"] == pytest.approx(1.0)
        assert bt.tl.unifrac(tdata, weighted=True).loc["s0", "s1"] == pytest.approx(1.0)
        assert bt.tl.unifrac(tdata, weighted=True, normalized=False).loc["s0", "s1"] == pytest.approx(4.0)


    def test_unifrac_multifurcating_root_is_accepted():
        # toy's root has three children; scikit-bio alone raises "The tree must be rooted."
        out = bt.tl.unifrac(bt.datasets.toy())
        assert out.shape == (6, 6) and out.loc["s1", "s4"] == pytest.approx(0.0916030534)


    def test_unifrac_after_filtering_keeps_path_lengths():
        # Filtering leaves unary nodes in the tree; distances between the kept features' samples must not change.
        tdata = bt.datasets.toy()
        dense = tdata.X.toarray()
        dense[:, [4, 7]] = 0  # f5, f8 absent everywhere
        tdata.X = sp.csr_matrix(dense)
        filtered = bt.pp.filter_features(tdata, min_total=1)
        assert filtered.n_vars == 6
        pd.testing.assert_frame_equal(bt.tl.unifrac(filtered), bt.tl.unifrac(tdata))
        pd.testing.assert_frame_equal(bt.tl.unifrac(filtered, weighted=True), bt.tl.unifrac(tdata, weighted=True))


    def test_unifrac_all_zero_sample_single_sample_and_missing_rank():
        tdata = bt.datasets.toy()
        tdata.var["genus"] = np.nan
        assert bt.tl.unifrac(tdata[:1].copy()).shape == (1, 1)
        dense = tdata.X.toarray()
        dense[0] = 0
        tdata.X = sp.csr_matrix(dense)
        assert bt.tl.unifrac(tdata, weighted=True).shape == (6, 6)


    def test_unifrac_inplace_writes_the_contract_keys():
        tdata = bt.datasets.toy()
        assert bt.tl.unifrac(tdata, inplace=True) is None
        assert bt.tl.unifrac(tdata, weighted=True, inplace=True) is None
        assert {"unweighted_unifrac", "weighted_unifrac"} <= set(tdata.obsp.keys())


    def test_unifrac_input_unchanged(assert_unchanged):
        tdata = bt.datasets.toy()
        before = tdata.copy()
        bt.tl.unifrac(tdata, weighted=True)
        assert_unchanged(before, tdata)


    def test_unifrac_without_a_tree_names_it():
        with pytest.raises(TypeError, match="needs a TreeData"):
            bt.tl.unifrac(_adata(np.array([[1, 2]])))
        tdata = bt.datasets.toy()
        del tdata.vart["phylo"]
        with pytest.raises(KeyError, match="phylo"):
            bt.tl.unifrac(tdata)


    @given(arrays(np.int64, st.tuples(st.integers(1, 5), st.just(8)), elements=st.integers(0, 20)), st.booleans())
    def test_unifrac_is_symmetric_with_zero_diagonal(dense, weighted):
        tdata = bt.datasets.toy()[: dense.shape[0]].copy()
        tdata.X = sp.csr_matrix(dense)
        out = bt.tl.unifrac(tdata, weighted=weighted).to_numpy()
        np.testing.assert_array_equal(out, out.T)
        np.testing.assert_array_equal(np.diag(out), 0.0)
    ```
  - `tests/tl/test_beta_golden.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pytest

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    def _square(path: Path) -> pd.DataFrame:
        return pd.read_csv(path, dtype={"sample_id": str}).set_index("sample_id")


    @pytest.mark.parametrize("metric", ["braycurtis", "jaccard"])
    def test_beta_matches_phyloseq_distance(metric):
        # jaccard: phyloseq::distance(GP, "jaccard", binary = TRUE).
        golden = _square(GOLDEN / "global_patterns" / f"beta_{metric}.csv.gz")
        out = bt.tl.beta(bt.datasets.global_patterns(), metric=metric)
        assert list(out.index) == list(golden.index) and list(out.columns) == list(golden.columns)
        np.testing.assert_allclose(out.to_numpy(), golden.to_numpy(), rtol=1e-7)


    @pytest.mark.parametrize("dataset", ["global_patterns", "esophagus"])
    @pytest.mark.parametrize("weighted", [False, True])
    def test_unifrac_matches_phyloseq_unifrac(dataset, weighted):
        name = "weighted" if weighted else "unweighted"
        golden = _square(GOLDEN / dataset / f"unifrac_{name}.csv.gz")
        out = bt.tl.unifrac(getattr(bt.datasets, dataset)(), weighted=weighted)
        assert list(out.index) == list(golden.index)
        np.testing.assert_allclose(out.to_numpy(), golden.to_numpy(), rtol=1e-7)
    ```
- [x] **Step 2: Run, expect failure.** Run `uv run --group test pytest tests/tl/test_beta.py -q`.
  It fails with `AttributeError: module 'biotapy.tl' has no attribute 'beta'`.
- [x] **Step 3: Implement.**
  - `src/biotapy/tl/_beta.py`:
    ```python
    """Beta diversity: sample x sample distance matrices in ``obsp``."""

    from typing import Literal, get_args

    import pandas as pd
    from anndata import AnnData
    from skbio import DistanceMatrix
    from skbio.diversity import beta_diversity

    from biotapy._core import TreeData, as_csr, get_skbio_tree

    BetaMetric = Literal["braycurtis", "jaccard"]


    def beta(adata: AnnData, *, metric: BetaMetric = "braycurtis", inplace: bool = False) -> pd.DataFrame | None:
        """Distances between every pair of samples.

        Parameters
        ----------
        adata
            Samples x features.
        metric
            ``"braycurtis"``, or ``"jaccard"`` on presence/absence.
        inplace
            Write the matrix to ``obsp[metric]`` and return ``None``.

        Returns
        -------
        pandas.DataFrame or None
            Symmetric samples x samples distances with a zero diagonal, indexed by
            ``obs_names``. Two all-zero samples are NaN apart under Bray-Curtis.

        Raises
        ------
        ValueError
            ``metric`` is not ``"braycurtis"`` or ``"jaccard"``.

        Notes
        -----
        R equivalent: ``phyloseq::distance``
        Guide: :doc:`/guide/diversity`

        ``"jaccard"`` matches ``phyloseq::distance(physeq, "jaccard", binary = TRUE)``:
        without ``binary = TRUE``, vegan computes a quantitative Jaccard instead.
        scikit-bio needs dense input, so ``X`` is densified once (8 bytes x samples x
        features) and the result takes 8 bytes x samples x samples.

        Examples
        --------
        >>> import biotapy as bt
        >>> round(float(bt.tl.beta(bt.datasets.toy()).loc["s1", "s4"]), 3)
        0.708
        """
        if metric not in get_args(BetaMetric):
            msg = f"metric must be one of {list(get_args(BetaMetric))}, got {metric!r}"
            raise ValueError(msg)
        # scikit-bio rejects sparse input (rules.md R6.2): one dense copy of X.
        distances = beta_diversity(metric, as_csr(adata.X).toarray(), ids=adata.obs_names.tolist())
        return _store(adata, distances, key=metric, inplace=inplace)


    def unifrac(
        tdata: TreeData, *, weighted: bool = False, normalized: bool = True, inplace: bool = False
    ) -> pd.DataFrame | None:
        """UniFrac distances between every pair of samples.

        Parameters
        ----------
        tdata
            Samples x features with the phylogeny in ``vart['phylo']``.
        weighted
            Weight branches by abundance (weighted UniFrac) instead of presence.
        normalized
            Scale weighted UniFrac to 0-1, as phyloseq does. Ignored when unweighted.
        inplace
            Write the matrix to ``obsp['unweighted_unifrac']`` or
            ``obsp['weighted_unifrac']`` and return ``None``.

        Returns
        -------
        pandas.DataFrame or None
            Symmetric samples x samples distances with a zero diagonal, indexed by ``obs_names``.

        Raises
        ------
        TypeError
            ``tdata`` is an AnnData that is not a TreeData.
        KeyError
            ``tdata`` has no ``vart['phylo']``.

        Notes
        -----
        R equivalent: ``phyloseq::UniFrac``
        Guide: :doc:`/guide/diversity`

        scikit-bio needs a root with at most two children. A root with more keeps its
        first child and gets the others under one new zero-length branch, so the tree
        is used rooted where it is drawn; phyloseq instead roots an unrooted tree at a
        random tip. scikit-bio's weighted UniFrac is unnormalized by default; biotapy
        passes ``normalized`` explicitly and defaults to phyloseq's ``TRUE``.
        ``X`` is densified once (8 bytes x samples x features).

        Examples
        --------
        >>> import biotapy as bt
        >>> round(float(bt.tl.unifrac(bt.datasets.toy()).loc["s1", "s4"]), 3)
        0.092
        """
        tree = get_skbio_tree(tdata)
        counts, ids, taxa = as_csr(tdata.X).toarray(), tdata.obs_names.tolist(), tdata.var_names.tolist()
        if weighted:
            distances = beta_diversity("weighted_unifrac", counts, ids=ids, taxa=taxa, tree=tree, normalized=normalized)
        else:
            distances = beta_diversity("unweighted_unifrac", counts, ids=ids, taxa=taxa, tree=tree)
        return _store(tdata, distances, key="weighted_unifrac" if weighted else "unweighted_unifrac", inplace=inplace)


    def _store(adata: AnnData, distances: DistanceMatrix, *, key: str, inplace: bool) -> pd.DataFrame | None:
        frame = distances.to_data_frame()
        if not inplace:
            return frame
        adata.obsp[key] = frame.to_numpy()
        return None
    ```
  - `tl/__init__.py` gains `from ._beta import beta, unifrac`.
  - `__all__` becomes `["alpha", "beta", "unifrac"]`.
  - scikit-bio auto-qualifies `jaccard` to presence/absence
    (`_qualitative_metrics`). `weighted_unifrac`'s own default is
    `normalized=False`, so `normalized` is always passed.
  - This state of the file was gated alone in the prototype: ruff, format,
    mypy, and 32 tl tests.
- [x] **Step 4: Docs.**
  - Add `tl.beta` and `tl.unifrac` to the Tools block of `docs/api.md`.
  - Append to `docs/guide/diversity.md`:
    ````markdown

    ## Beta diversity: `tl.beta` and `tl.unifrac`

    `bt.tl.beta` computes Bray-Curtis (`metric="braycurtis"`, the default) or Jaccard
    (`metric="jaccard"`) distances between every pair of samples; `bt.tl.unifrac` computes
    unweighted or weighted UniFrac along the tree:

    ```python
    bt.tl.beta(tdata, inplace=True)  # obsp["braycurtis"]
    bt.tl.unifrac(tdata, weighted=True, inplace=True)  # obsp["weighted_unifrac"]
    ```

    - Jaccard is on presence/absence, like `phyloseq::distance(physeq, "jaccard", binary = TRUE)`.
      Without `binary = TRUE` phyloseq computes vegan's quantitative Jaccard, a different number.
    - Weighted UniFrac is normalized to 0-1 by default, as in phyloseq; pass `normalized=False` for
      the raw value.
    - A tree whose root has three or more children is used rooted where it is drawn. phyloseq
      instead roots such a tree at a random tip, so its UniFrac changes from run to run.
    - Two all-zero samples are `NaN` apart under Bray-Curtis. Drop empty samples with
      `bt.pp.filter_samples(tdata, 1)` before ordinating.
    - The table is densified once: 8 bytes x samples x features, plus 8 bytes x samples x samples
      for the result.
    ````
- [x] **Step 5: Knowledge.**
  - No concept states anything this changes: `obsp` keys are already in
    data-model-slots.
  - Add a log line. Tick 1.16 here.
- [x] **Step 6: Run, gate and commit.**
  - Tests: `uv run --group test pytest tests/tl -q`.
  - Golden: `uv run --group test pytest -m golden tests/tl/test_beta_golden.py -q`
    (6 passed; the esophagus golden test also downloads 1,840 B).
  - The three gates.
  - Commit `feat(tl): add beta diversity and UniFrac`.

### Task 1.17: `tl.pcoa`, `tl.nmds`, `tl.permanova`

**Files:**
- Create:
  - `src/biotapy/tl/_ordination.py`, `src/biotapy/tl/_permanova.py`
  - `tests/tl/test_ordination.py`, `tests/tl/test_ordination_golden.py`
  - `tests/tl/test_permanova.py`, `tests/tl/test_permanova_golden.py`
  - `docs/guide/ordination.md`
- Modify:
  - `src/biotapy/tl/_beta.py` (add `stored_distances`), `src/biotapy/tl/__init__.py`
  - `src/biotapy/_core/_slots.py` and `tests/core/test_slots.py`
    (`uns["biotapy"]` propagation)
  - `pyproject.toml`, `uv.lock` (scikit-learn)
  - `docs/api.md`, `docs/guide/index.md`, `docs/guide/data_model.md`, `docs/guide/filtering.md`
  - `.knowledge/contracts/data-model-slots.md`
  - `.knowledge/decisions/optional-heavy-dependencies.md`
  - `.knowledge/decisions/pure-by-default.md` (see the ruling in Step 6)
  - `.knowledge/modules/core.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:**
  - `obsp` matrices from 1.16, and `as_generator`;
  - `skbio.stats.ordination.pcoa`, `skbio.stats.distance.permanova`,
    `sklearn.manifold.MDS`;
  - the `pcoa_*`, `nmds_*` and `permanova_*` golden files.
- **Produces:**
  - `bt.tl.pcoa(adata: AnnData, *, distance: str = "braycurtis", n_components: int = 10, inplace: bool = False) -> tuple[pd.DataFrame, pd.DataFrame] | None`.
  - `bt.tl.nmds(adata: AnnData, *, distance: str = "braycurtis", n_components: int = 2, seed: int | np.random.Generator | None = None, inplace: bool = False) -> tuple[pd.DataFrame, float] | None`.
  - `bt.tl.permanova(adata: AnnData, grouping: str, *, distance: str = "braycurtis", permutations: int = 999, seed: int | np.random.Generator | None = None) -> pd.Series`.
  - `tl._beta.stored_distances(adata: AnnData, key: str) -> DistanceMatrix`
    (private; `tl` only).
  - `_core.feature_subset` keeps only `x_kind` and `provenance` in `uns["biotapy"]`.

- [x] **Step 1: Dependency.**
  - In `[project] dependencies`, add `"scikit-learn>=1.8",` after `"scikit-bio>=0.7.4,<0.8",`.
  - Add the mypy overrides after the pooch pair:
    ```toml
      # scikit-learn 1.9.1 ships no py.typed marker either; same treatment.
      { module = "sklearn", follow_untyped_imports = true, implicit_reexport = true },
      { module = "sklearn.*", follow_untyped_imports = true, implicit_reexport = true },
    ```
  - Replace the `untyped_calls_exclude` comment and line with:
    ```toml
    # biom-format and scikit-learn ship no type annotations at all, nor does scikit-bio's
    # subsample_counts; exempt only calls into them (not our own code) from strict's
    # disallow_untyped_calls (mypy/checkexpr.py matches by callee fullname).
    untyped_calls_exclude = [ "biom", "skbio.stats._subsample", "sklearn" ]
    ```
  - Run `uv sync --group dev --group test --group doc`. It resolves scikit-learn
    1.9.1 plus joblib, threadpoolctl and cloudpickle; wheels exist for
    CPython 3.12-3.14.
  - Confirm that `uv run python -c "import sklearn; print(sklearn.__version__)"`
    prints `1.9.1` or later, and that `sklearn/py.typed` does not exist
    (R2.2). `pyproject-fmt` leaves this `pyproject.toml` unchanged, as in the
    prototype.
- [x] **Step 2: Failing tests.**
  - Append to `tests/core/test_slots.py`:
    ```python
    def test_feature_subset_drops_ordination_metadata():
        adata = _adata()
        adata.uns["biotapy"] = {
            "x_kind": "counts",
            "provenance": ["{}"],
            "pcoa": {"eigenvalues": np.ones(2)},
            "nmds": {"stress": 0.1},
        }
        assert feature_subset(adata, np.array([0])).uns["biotapy"] == {"x_kind": "counts", "provenance": ["{}"]}
    ```
  - `tests/tl/test_ordination.py`:
    ```python
    import anndata as ad
    import numpy as np
    import pandas as pd
    import pytest
    import scipy.sparse as sp
    import treedata as td
    from scipy.spatial import procrustes
    from scipy.spatial.distance import pdist, squareform

    import biotapy as bt


    def _toy_with(metric="braycurtis"):
        tdata = bt.datasets.toy()
        bt.tl.beta(tdata, metric=metric, inplace=True)
        return tdata


    def test_pcoa_reproduces_euclidean_distances():
        # PCoA of Euclidean distances recovers them exactly from all positive axes.
        tdata = bt.datasets.toy()
        dense = tdata.X.toarray().astype(float)
        tdata.obsp["euclidean"] = squareform(pdist(dense))
        coords, axes = bt.tl.pcoa(tdata, distance="euclidean")
        assert list(coords.columns) == [f"PC{i}" for i in range(1, 6)]  # at most n_obs - 1 axes
        np.testing.assert_allclose(squareform(pdist(coords.to_numpy())), tdata.obsp["euclidean"], atol=1e-9)
        assert axes["proportion_explained"].sum() == pytest.approx(1.0)


    def test_pcoa_proportion_divides_by_all_eigenvalues():
        coords, axes = bt.tl.pcoa(_toy_with(), n_components=2)
        assert coords.shape == (6, 2) and list(axes.index) == ["PC1", "PC2"]
        assert (axes["proportion_explained"] <= 1).all() and axes["eigenvalue"].is_monotonic_decreasing


    def test_pcoa_inplace_writes_obsm_and_uns():
        tdata = _toy_with()
        coords, axes = bt.tl.pcoa(tdata, n_components=3)
        assert bt.tl.pcoa(tdata, n_components=3, inplace=True) is None
        np.testing.assert_array_equal(tdata.obsm["X_pcoa"], coords.to_numpy())
        np.testing.assert_array_equal(tdata.uns["biotapy"]["pcoa"]["eigenvalues"], axes["eigenvalue"].to_numpy())
        np.testing.assert_array_equal(
            tdata.uns["biotapy"]["pcoa"]["proportion_explained"], axes["proportion_explained"].to_numpy()
        )


    def test_pcoa_input_unchanged(assert_unchanged):
        tdata = _toy_with()
        before = tdata.copy()
        bt.tl.pcoa(tdata)
        assert_unchanged(before, tdata)


    def test_pcoa_missing_distance_names_the_call():
        with pytest.raises(KeyError, match=r"bt.tl.unifrac\(tdata, weighted=True, inplace=True\)"):
            bt.tl.pcoa(bt.datasets.toy(), distance="weighted_unifrac")


    def test_pcoa_nan_distances_raise():
        tdata = bt.datasets.toy()
        dense = tdata.X.toarray()
        dense[[0, 1]] = 0  # two all-zero samples are NaN apart under Bray-Curtis
        tdata.X = sp.csr_matrix(dense)
        bt.tl.beta(tdata, inplace=True)
        with pytest.raises(ValueError, match="NaN"):
            bt.tl.pcoa(tdata)


    def test_pcoa_single_sample_and_bad_components_raise():
        with pytest.raises(ValueError, match="at least 2 samples"):
            bt.tl.pcoa(_toy_with()[:1].copy())
        with pytest.raises(ValueError, match="n_components"):
            bt.tl.pcoa(_toy_with(), n_components=0)


    def test_pcoa_after_filter_samples_uses_the_subset():
        tdata = bt.pp.filter_samples(_toy_with(), 70)
        coords, _ = bt.tl.pcoa(tdata)
        assert list(coords.index) == ["s3", "s4", "s5", "s6"] and coords.shape[1] == 3


    def test_nmds_same_seed_same_result():
        first, stress = bt.tl.nmds(_toy_with(), seed=3)
        again, stress_again = bt.tl.nmds(_toy_with(), seed=np.random.default_rng(3))
        np.testing.assert_array_equal(first.to_numpy(), again.to_numpy())
        assert stress == stress_again and 0 <= stress < 0.2
        assert list(first.columns) == ["NMDS1", "NMDS2"]


    def test_nmds_recovers_a_planar_configuration():
        points = np.random.default_rng(0).uniform(size=(10, 2))
        adata = ad.AnnData(X=sp.csr_matrix(np.ones((10, 1))), obs=pd.DataFrame(index=[f"s{i}" for i in range(10)]))
        adata.obsp["euclidean"] = squareform(pdist(points))
        coords, stress = bt.tl.nmds(adata, distance="euclidean", seed=0)
        _, _, disparity = procrustes(points, coords.to_numpy())
        assert stress < 0.05 and np.sqrt(1 - disparity) > 0.99


    def test_nmds_inplace_writes_obsm_and_stress():
        tdata = _toy_with()
        coords, stress = bt.tl.nmds(tdata, seed=0)
        assert bt.tl.nmds(tdata, seed=0, inplace=True) is None
        np.testing.assert_array_equal(tdata.obsm["X_nmds"], coords.to_numpy())
        assert tdata.uns["biotapy"]["nmds"] == {"stress": stress}


    def test_nmds_input_unchanged(assert_unchanged):
        tdata = _toy_with()
        before = tdata.copy()
        bt.tl.nmds(tdata, seed=0)
        assert_unchanged(before, tdata)


    def test_nmds_too_few_samples_raise():
        with pytest.raises(ValueError, match="more than n_components"):
            bt.tl.nmds(_toy_with()[:3].copy(), seed=0)


    def test_feature_changes_drop_stored_ordinations():
        tdata = _toy_with()
        bt.tl.pcoa(tdata, inplace=True)
        bt.tl.nmds(tdata, seed=0, inplace=True)
        out = bt.pp.filter_features(tdata, min_total=19)
        assert not {"X_pcoa", "X_nmds"} & set(out.obsm.keys())
        assert set(out.uns["biotapy"]) == {"x_kind", "provenance"}


    def test_stored_results_survive_h5td(tmp_path):
        tdata = _toy_with()
        bt.tl.alpha(tdata, inplace=True)
        bt.tl.unifrac(tdata, weighted=True, inplace=True)
        bt.tl.pcoa(tdata, inplace=True)
        bt.tl.nmds(tdata, seed=0, inplace=True)
        tdata.write_h5td(tmp_path / "toy.h5td")
        back = td.read_h5td(tmp_path / "toy.h5td")
        np.testing.assert_array_equal(back.obsp["weighted_unifrac"], tdata.obsp["weighted_unifrac"])
        np.testing.assert_array_equal(back.obsm["X_pcoa"], tdata.obsm["X_pcoa"])
        np.testing.assert_array_equal(
            back.uns["biotapy"]["pcoa"]["eigenvalues"], tdata.uns["biotapy"]["pcoa"]["eigenvalues"]
        )
        assert back.uns["biotapy"]["nmds"]["stress"] == tdata.uns["biotapy"]["nmds"]["stress"]
        np.testing.assert_array_equal(back.obs["alpha_shannon"], tdata.obs["alpha_shannon"])
    ```
  - `tests/tl/test_ordination_golden.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pytest
    from scipy.spatial import procrustes

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    def _with_braycurtis():
        tdata = bt.datasets.global_patterns()
        bt.tl.beta(tdata, inplace=True)
        return tdata


    def test_pcoa_matches_ape_pcoa():
        vectors = pd.read_csv(GOLDEN / "pcoa_braycurtis_vectors.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
        values = pd.read_csv(GOLDEN / "pcoa_braycurtis_values.csv.gz")
        coords, axes = bt.tl.pcoa(_with_braycurtis(), n_components=10)
        assert list(coords.index) == list(vectors.index)
        np.testing.assert_allclose(axes["eigenvalue"], values["eigenvalue"], rtol=1e-6)
        np.testing.assert_allclose(axes["proportion_explained"], values["relative_eig"], rtol=1e-6)
        ours, theirs = coords.to_numpy(), vectors.to_numpy()
        # Eigenvector signs are arbitrary: align each axis before comparing.
        signs = np.sign((ours * theirs).sum(axis=0))
        np.testing.assert_allclose(ours * signs, theirs, rtol=1e-6)


    def test_nmds_matches_vegan_metamds():
        points = pd.read_csv(GOLDEN / "nmds_braycurtis_points.csv.gz", dtype={"sample_id": str}).set_index("sample_id")
        r_stress = pd.read_csv(GOLDEN / "nmds_braycurtis_stress.csv.gz")["stress"].iloc[0]
        coords, stress = bt.tl.nmds(_with_braycurtis(), seed=0)
        assert list(coords.index) == list(points.index)
        assert abs(stress - r_stress) < 0.02
        _, _, disparity = procrustes(points.to_numpy(), coords.to_numpy())
        assert np.sqrt(1 - disparity) > 0.95
    ```
  - `tests/tl/test_permanova.py`:
    ```python
    import numpy as np
    import pytest

    import biotapy as bt


    def _toy_with():
        tdata = bt.datasets.toy()
        bt.tl.beta(tdata, inplace=True)
        return tdata


    def test_permanova_pseudo_f_by_hand():
        tdata = _toy_with()
        d2 = tdata.obsp["braycurtis"] ** 2
        total = d2[np.triu_indices(6, 1)].sum() / 6
        within = sum(d2[np.ix_(idx, idx)][np.triu_indices(3, 1)].sum() / 3 for idx in ([0, 1, 2], [3, 4, 5]))
        expected = (total - within) / (within / 4)
        result = bt.tl.permanova(tdata, "group", seed=0)
        assert result["test statistic"] == pytest.approx(expected)
        assert result["number of groups"] == 2 and result["sample size"] == 6


    def test_permanova_same_seed_same_p_value():
        first = bt.tl.permanova(_toy_with(), "group", permutations=99, seed=1)
        again = bt.tl.permanova(_toy_with(), "group", permutations=99, seed=np.random.default_rng(1))
        assert first["p-value"] == again["p-value"] and first["number of permutations"] == 99


    def test_permanova_input_unchanged(assert_unchanged):
        tdata = _toy_with()
        before = tdata.copy()
        bt.tl.permanova(tdata, "group", seed=0)
        assert_unchanged(before, tdata)


    def test_permanova_unknown_grouping_is_named():
        with pytest.raises(KeyError, match="grouping='site'"):
            bt.tl.permanova(_toy_with(), "site")


    def test_permanova_missing_group_raises():
        tdata = _toy_with()
        tdata.obs["group"] = tdata.obs["group"].cat.add_categories("C")
        tdata.obs.loc["s1", "group"] = np.nan
        with pytest.raises(ValueError, match=r"missing for 1 sample"):
            bt.tl.permanova(tdata, "group")


    def test_permanova_missing_distance_names_the_call():
        with pytest.raises(KeyError, match=r"bt.tl.beta\(adata, metric='jaccard', inplace=True\)"):
            bt.tl.permanova(bt.datasets.toy(), "group", distance="jaccard")
    ```
  - `tests/tl/test_permanova_golden.py`:
    ```python
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import pytest

    import biotapy as bt

    GOLDEN = Path(__file__).parents[1] / "golden" / "global_patterns"
    pytestmark = [pytest.mark.golden, pytest.mark.network]


    def test_permanova_matches_vegan_adonis2():
        golden = pd.read_csv(GOLDEN / "permanova_sampletype.csv.gz").iloc[0]
        tdata = bt.datasets.global_patterns()
        bt.tl.beta(tdata, inplace=True)
        result = bt.tl.permanova(tdata, "SampleType", permutations=9999, seed=0)
        np.testing.assert_allclose(result["test statistic"], golden["f"], rtol=1e-7)
        assert result["number of groups"] == golden["df"] + 1
        assert abs(result["p-value"] - golden["p"]) <= 0.02
    ```
  - Test choices:
    - `test_nmds_recovers_a_planar_configuration` uses 10 random planar points
      rather than toy. On toy, NMDS falls into the degenerate two-cluster
      solution (stress 0.002), where a rank-order check fails (73% agreement)
      without any bug.
    - The p-value tolerance in the PERMANOVA golden test is the contract's 0.02.
      Both sides give 1e-4, the smallest p at 9,999 permutations.
- [x] **Step 3: Run, expect failure.** Run `uv run --group test pytest tests/tl tests/core/test_slots.py -q`.
  - `AttributeError: module 'biotapy.tl' has no attribute 'pcoa'` (and `nmds`, `permanova`).
  - `test_feature_subset_drops_ordination_metadata` fails with an
    `AssertionError`: the result still holds `pcoa` and `nmds`.
- [x] **Step 4: Implement.**
  - `_core/_slots.py`:
    - after `DERIVED_SLOTS` add:
      ```python
      # The uns['biotapy'] keys that survive a feature change; the others (pcoa, nmds) described the old features.
      KEPT_META = ("x_kind", "provenance")
      ```
    - replace `feature_subset`'s last `out.uns = ...` line with:
      ```python
          meta = out.uns.get("biotapy", {"x_kind": "counts"})
          out.uns = {"biotapy": {key: meta[key] for key in KEPT_META if key in meta}}
      ```
  - `tl/_beta.py`:
    - add `import numpy as np` to the imports;
    - after `BetaMetric` add:
      ```python
      # The call that writes each obsp key, for the error raised when a key is missing.
      _WRITTEN_BY = {
          "braycurtis": "bt.tl.beta(adata, metric='braycurtis', inplace=True)",
          "jaccard": "bt.tl.beta(adata, metric='jaccard', inplace=True)",
          "unweighted_unifrac": "bt.tl.unifrac(tdata, inplace=True)",
          "weighted_unifrac": "bt.tl.unifrac(tdata, weighted=True, inplace=True)",
      }
      ```
    - append:
      ```python
      def stored_distances(adata: AnnData, key: str) -> DistanceMatrix:
          """``obsp[key]`` as a scikit-bio DistanceMatrix; a missing key names the call that writes it."""
          if key not in adata.obsp:
              call = _WRITTEN_BY.get(key, "bt.tl.beta or bt.tl.unifrac with inplace=True")
              msg = f"distance={key!r}: no obsp[{key!r}]; run {call} first"
              raise KeyError(msg)
          values = np.asarray(adata.obsp[key], dtype=np.float64)
          if np.isnan(values).any():
              msg = f"distance={key!r}: obsp[{key!r}] holds NaN, as between two all-zero samples; drop them first"
              raise ValueError(msg)
          return DistanceMatrix(values, ids=adata.obs_names.tolist())
      ```
  - `src/biotapy/tl/_ordination.py`:
    ```python
    """Ordination of a stored distance matrix: PCoA and non-metric MDS."""

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from skbio.stats.ordination import pcoa as skbio_pcoa
    from sklearn.manifold import MDS

    from biotapy._core import as_generator

    from ._beta import stored_distances

    # vegan::metaMDS's try = 20: the best of 20 random starts.
    _NMDS_STARTS = 20


    def pcoa(
        adata: AnnData, *, distance: str = "braycurtis", n_components: int = 10, inplace: bool = False
    ) -> tuple[pd.DataFrame, pd.DataFrame] | None:
        """Principal coordinates analysis of a distance matrix in ``obsp``.

        Parameters
        ----------
        adata
            Samples x features with ``obsp[distance]``, written by :func:`biotapy.tl.beta`
            or :func:`biotapy.tl.unifrac` with ``inplace=True``.
        distance
            The ``obsp`` key to ordinate.
        n_components
            Axes to keep; at most ``n_obs - 1`` are returned.
        inplace
            Write the coordinates to ``obsm['X_pcoa']`` and the axes to
            ``uns['biotapy']['pcoa']`` (``eigenvalues``, ``proportion_explained``), and return ``None``.

        Returns
        -------
        tuple of pandas.DataFrame, or None
            Coordinates (samples x ``PC1``..), and per axis its ``eigenvalue`` and
            ``proportion_explained`` (eigenvalue / sum of all eigenvalues, as ``ape::pcoa``'s ``Relative_eig``).

        Raises
        ------
        KeyError
            ``obsp[distance]`` is missing; the message names the call that writes it.
        ValueError
            ``n_components`` is below 1, there are fewer than 2 samples, or the distances hold NaN.

        Notes
        -----
        R equivalent: ``phyloseq::ordinate``, ``ape::pcoa``
        Guide: :doc:`/guide/ordination`

        Like ``ordinate(physeq, "PCoA", ...)`` (``ape::pcoa`` with ``correction = "none"``).
        Axis signs are arbitrary, as in R.

        Examples
        --------
        >>> import biotapy as bt
        >>> tdata = bt.datasets.toy()
        >>> bt.tl.beta(tdata, inplace=True)
        >>> coords, axes = bt.tl.pcoa(tdata, n_components=2)
        >>> coords.shape
        (6, 2)
        """
        distances = stored_distances(adata, distance)
        if n_components < 1 or adata.n_obs < 2:
            msg = f"tl.pcoa needs n_components >= 1 and at least 2 samples, got {n_components} and {adata.n_obs}"
            raise ValueError(msg)
        # With fewer axes than samples, scikit-bio divides by the trace of the centred matrix, as ape::pcoa does.
        dimensions = min(n_components, adata.n_obs - 1)
        result = skbio_pcoa(distances, dimensions=dimensions)
        names = [f"PC{i}" for i in range(1, dimensions + 1)]
        coords = pd.DataFrame(result.samples.to_numpy(), index=adata.obs_names, columns=names)
        axes = pd.DataFrame(
            {"eigenvalue": result.eigvals.to_numpy(), "proportion_explained": result.proportion_explained.to_numpy()},
            index=names,
        )
        if not inplace:
            return coords, axes
        adata.obsm["X_pcoa"] = coords.to_numpy()
        adata.uns.setdefault("biotapy", {})["pcoa"] = {
            "eigenvalues": axes["eigenvalue"].to_numpy(),
            "proportion_explained": axes["proportion_explained"].to_numpy(),
        }
        return None


    def nmds(
        adata: AnnData,
        *,
        distance: str = "braycurtis",
        n_components: int = 2,
        seed: int | np.random.Generator | None = None,
        inplace: bool = False,
    ) -> tuple[pd.DataFrame, float] | None:
        """Non-metric multidimensional scaling of a distance matrix in ``obsp``.

        Parameters
        ----------
        adata
            Samples x features with ``obsp[distance]``, written by :func:`biotapy.tl.beta`
            or :func:`biotapy.tl.unifrac` with ``inplace=True``.
        distance
            The ``obsp`` key to ordinate.
        n_components
            Dimensions of the configuration.
        seed
            Seed or generator for the random starts.
        inplace
            Write the configuration to ``obsm['X_nmds']`` and the stress to
            ``uns['biotapy']['nmds']['stress']``, and return ``None``.

        Returns
        -------
        tuple of pandas.DataFrame and float, or None
            The configuration (samples x ``NMDS1``..) and its Kruskal stress-1.

        Raises
        ------
        KeyError
            ``obsp[distance]`` is missing; the message names the call that writes it.
        ValueError
            ``n_components`` is below 1, there are not more samples than
            ``n_components + 1``, or the distances hold NaN.

        Notes
        -----
        R equivalent: ``phyloseq::ordinate``, ``vegan::metaMDS``
        Guide: :doc:`/guide/ordination`

        scikit-learn's SMACOF, non-metric, keeps the best of 20 random starts, as
        ``vegan::metaMDS`` does by default (``try = 20``). The configuration is not
        rotated or scaled afterwards, unlike metaMDS's ``postMDS``, so compare
        configurations up to rotation and scale (Procrustes).

        Examples
        --------
        >>> import biotapy as bt
        >>> tdata = bt.datasets.toy()
        >>> bt.tl.beta(tdata, inplace=True)
        >>> coords, stress = bt.tl.nmds(tdata, seed=0)
        >>> coords.shape
        (6, 2)
        """
        distances = stored_distances(adata, distance)
        if n_components < 1 or adata.n_obs <= n_components + 1:
            msg = f"tl.nmds needs n_components >= 1 and more than n_components + 1 samples, got {n_components} and {adata.n_obs}"
            raise ValueError(msg)
        # scikit-learn takes an int or a RandomState, not a Generator: draw its seed from ours.
        random_state = int(as_generator(seed).integers(2**31 - 1))
        mds = MDS(
            n_components=n_components,
            metric_mds=False,
            metric="precomputed",
            n_init=_NMDS_STARTS,
            init="random",
            normalized_stress="auto",
            random_state=random_state,
        )
        names = [f"NMDS{i}" for i in range(1, n_components + 1)]
        coords = pd.DataFrame(mds.fit_transform(distances.data), index=adata.obs_names, columns=names)
        stress = float(mds.stress_)
        if not inplace:
            return coords, stress
        adata.obsm["X_nmds"] = coords.to_numpy()
        adata.uns.setdefault("biotapy", {})["nmds"] = {"stress": stress}
        return None
    ```
  - `src/biotapy/tl/_permanova.py`:
    ```python
    """PERMANOVA: do groups of samples differ in a stored distance matrix?"""

    import numpy as np
    import pandas as pd
    from anndata import AnnData
    from skbio.stats.distance import permanova as skbio_permanova

    from biotapy._core import as_generator

    from ._beta import stored_distances


    def permanova(
        adata: AnnData,
        grouping: str,
        *,
        distance: str = "braycurtis",
        permutations: int = 999,
        seed: int | np.random.Generator | None = None,
    ) -> pd.Series:
        """Permutational multivariate analysis of variance of a distance matrix in ``obsp``.

        Parameters
        ----------
        adata
            Samples x features with ``obsp[distance]``, written by :func:`biotapy.tl.beta`
            or :func:`biotapy.tl.unifrac` with ``inplace=True``.
        grouping
            The ``obs`` column holding each sample's group.
        distance
            The ``obsp`` key to test.
        permutations
            Permutations for the p-value.
        seed
            Seed or generator for the permutations.

        Returns
        -------
        pandas.Series
            scikit-bio's result: ``test statistic`` (pseudo-F), ``p-value``,
            ``sample size``, ``number of groups``, ``number of permutations`` and the
            method and statistic names.

        Raises
        ------
        KeyError
            ``grouping`` is not an ``obs`` column, or ``obsp[distance]`` is missing.
        ValueError
            ``grouping`` has missing values, or the distances hold NaN.

        Notes
        -----
        R equivalent: ``vegan::adonis2``
        Guide: :doc:`/guide/ordination`

        Matches ``adonis2(distance ~ grouping, data, permutations)`` with one term,
        whose ``F`` is the test statistic. P-values agree only up to permutation noise:
        R and NumPy random generators differ.

        Examples
        --------
        >>> import biotapy as bt
        >>> tdata = bt.datasets.toy()
        >>> bt.tl.beta(tdata, inplace=True)
        >>> result = bt.tl.permanova(tdata, "group", seed=0)
        >>> int(result["number of groups"])
        2
        """
        if grouping not in adata.obs.columns:
            msg = f"grouping={grouping!r} is not a column of obs"
            raise KeyError(msg)
        groups = adata.obs[grouping]
        if groups.isna().any():
            msg = f"grouping={grouping!r} is missing for {int(groups.isna().sum())} sample(s); drop them first"
            raise ValueError(msg)
        distances = stored_distances(adata, distance)
        result: pd.Series = skbio_permanova(
            distances, groups.to_numpy(), permutations=permutations, seed=as_generator(seed)
        )
        return result
    ```
  - `tl/__init__.py`:
    ```python
    from ._alpha import alpha
    from ._beta import beta, unifrac
    from ._ordination import nmds, pcoa
    from ._permanova import permanova

    __all__ = ["alpha", "beta", "nmds", "pcoa", "permanova", "unifrac"]
    ```
  - Why the code looks as it does:
    - `pcoa` asks scikit-bio for at most `n_obs - 1` axes. With exactly `n`
      axes, `proportion_explained` would divide by the sum of positive
      eigenvalues instead of ape's trace.
    - `result: pd.Series = ...` in `permanova` is a typed binding, not a cast.
      scikit-bio's string annotation `'pd.Series'` resolves to `Any` under
      mypy, and returning it fails `no-any-return`.
- [x] **Step 5: Docs.**
  - `docs/api.md` Tools block: add `tl.nmds`, `tl.pcoa` and `tl.permanova`.
    The block is then `tl.alpha`, `tl.beta`, `tl.nmds`, `tl.pcoa`,
    `tl.permanova`, `tl.unifrac`.
  - Add `ordination` to the guide toctree after `diversity`.
  - Create `docs/guide/ordination.md`:
    ````markdown
    # Ordination and PERMANOVA

    Ordination and PERMANOVA read a distance matrix that `tl.beta` or `tl.unifrac` stored in
    `obsp`; `distance=` names it (default `"braycurtis"`). When it is missing, the error names the
    call to run.

    ```python
    import biotapy as bt

    tdata = bt.datasets.toy()
    bt.tl.beta(tdata, inplace=True)
    coords, axes = bt.tl.pcoa(tdata)  # PC1.. per sample; eigenvalue and proportion per axis
    coords, stress = bt.tl.nmds(tdata, seed=0)  # NMDS1, NMDS2 and Kruskal stress-1
    result = bt.tl.permanova(tdata, "group", seed=0)
    result["test statistic"], result["p-value"]
    ```

    ## PCoA

    `bt.tl.pcoa` matches `ordinate(physeq, "PCoA", distance)`, which runs `ape::pcoa` with no
    correction for negative eigenvalues. `proportion_explained` is each eigenvalue over the sum of
    all of them, ape's `Relative_eig`. At most `n_obs - 1` axes are returned, and axis signs are
    arbitrary, in R too. With `inplace=True` the coordinates go to `obsm["X_pcoa"]` and the axes to
    `uns["biotapy"]["pcoa"]`.

    ## NMDS

    `bt.tl.nmds` runs scikit-learn's non-metric SMACOF on the stored distances and keeps the best
    of 20 random starts, as `vegan::metaMDS` does by default. vegan then rotates and rescales its
    result (`postMDS`); biotapy does not, so compare configurations up to rotation and scale. With
    `inplace=True`: `obsm["X_nmds"]` and `uns["biotapy"]["nmds"]["stress"]`.

    ## PERMANOVA

    `bt.tl.permanova(tdata, grouping)` tests whether the groups in `obs[grouping]` differ, like
    `vegan::adonis2(distance ~ grouping, data)`. The pseudo-F statistic matches vegan exactly;
    p-values come from permutations and agree only up to that noise. Drop samples with a missing
    group first.

    ## After filtering

    `pp.filter_features` and `pp.rarefy` drop stored distances and ordinations, including
    `uns["biotapy"]["pcoa"]` and `["nmds"]`, since they described the old features.
    `pp.filter_samples` keeps them, subset to the remaining samples: distances stay valid, but an
    ordination is not recomputed, so run `tl.pcoa` again for the ordination of the subset.
    ````
  - `docs/guide/data_model.md`:
    - the `uns["biotapy"]` row becomes
      `` | `uns["biotapy"]` | biotapy's own bookkeeping, and ordination summaries | `x_kind`, `provenance`, `pcoa`, `nmds` | ``;
    - in "What survives a filter or an aggregation", "`varm`, `varp`, and every
      `uns` key other than `uns["biotapy"]`" becomes "`varm`, `varp`, the
      ordination summaries `uns["biotapy"]["pcoa"]` and `["nmds"]`, and every
      `uns` key other than `uns["biotapy"]`".
  - `docs/guide/filtering.md`: "every entry in `layers`, `obsm` and `obsp`."
    becomes "every entry in `layers`, `obsm` and `obsp`, and the ordination
    summaries in `uns["biotapy"]`."
- [x] **Step 6: Knowledge.**
  - `contracts/data-model-slots.md`:
    - The Slots row `uns["biotapy"]` keys become: `x_kind`, `provenance`,
      `pcoa` (`eigenvalues`, `proportion_explained`), `nmds` (`stress`).
    - Propagation, Feature-changing row:
      - Keeps "`uns["biotapy"]`" becomes "`uns["biotapy"]["x_kind"]` and `["provenance"]`".
      - Drops gains "`uns["biotapy"]["pcoa"]`, `["nmds"]`".
  - `.knowledge/modules/core.md`: in the `feature_subset` invariant, "every
    `uns` key except `biotapy`" becomes "every `uns` key except `biotapy`, and
    every `uns["biotapy"]` key except `x_kind` and `provenance`
    (`_slots.py:KEPT_META`)".
  - `decisions/optional-heavy-dependencies.md`:
    - "scikit-learn (Phase 1, NMDS - pending approval)" becomes
      "scikit-learn (`>=1.8`, Phase 1, NMDS)";
    - add: "Task 1.17 added scikit-learn (`>=1.8`, approved 2026-09-27):
      scikit-bio has no non-metric MDS, and 1.8 renamed `dissimilarity` to
      `metric`. It brings joblib, threadpoolctl and cloudpickle, and adds
      about 0.15 s to `import biotapy`."
  - `decisions/pure-by-default.md`, a verified decision, so edit it only
    because the user approves this plan's ruling. Add a Consequences bullet:
    "`tl.permanova` returns a test result, not per-sample or per-pair values,
    so it has no slot and no `inplace`. No `tl` function takes `key_added`
    until a use case needs one (rules.md R2.3)."
  - Add a log line. Tick 1.17 here.
- [x] **Step 7: Run, gate and commit.**
  - Tests: `uv run --group test pytest tests/tl tests/core tests/pp -q`.
  - Golden:
    `uv run --group test pytest -m golden tests/tl/test_ordination_golden.py tests/tl/test_permanova_golden.py -q`.
    Expect 3 passed. PERMANOVA takes about 10 s at 9,999 permutations, and
    NMDS about 5 s.
  - Record `uv run python -X importtime -c "import biotapy" 2>&1 | tail -1` in
    the report (R10.1; the prototype measured about 1.3 s, unchanged).
  - The three gates.
  - Commit `feat(tl): add PCoA, NMDS and PERMANOVA`, staging `uv.lock`
    together with `pyproject.toml`.

### Checkpoint C - review slice 1C
- [x] Review the whole slice (superpowers:requesting-code-review) against:
  - every contract: function-shape, data-model-slots, module-boundaries,
    tree-access and r-golden-parity;
  - the pure-by-default decision;
  - the slice 1C review focus.

  Then a fix pass, one commit per finding, each with a test.
- [x] Knowledge: write the `Module` concept `.knowledge/modules/tl.md` with the
  codebase-map template:
  - **Responsibility:** diversity, ordination and PERMANOVA over `obsp`; owns
    no reader and no transform.
  - **Entry points:** the six functions, and `_beta.py:stored_distances`.
  - **Invariants:**
    - `inplace` writes only the data-model-slots keys;
    - scikit-bio gets dense input, `alpha` in chunks of `2**20` values;
    - trees come only from `_core.get_skbio_tree`;
    - `pcoa` returns at most `n_obs - 1` axes, and its proportions divide by
      the trace;
    - `nmds` keeps the best of 20 starts;
    - `permanova` has no `inplace`.
  - **Dependencies:** `_core` (`as_csr`, `require_counts`, `get_skbio_tree`,
    `as_generator`), scikit-bio and scikit-learn.
  - **Verification:** `uv run --group test pytest tests/tl -q`, and
    `-m golden tests/tl` with a pooch cache.
  - **Gotchas:**
    - the all-zero-sample results per metric;
    - PCoA axis signs are arbitrary, and tests align them per axis;
    - NMDS on toy is degenerate (stress 0.002), so tests use planar points;
    - the PERMANOVA golden test takes about 10 s at 9,999 permutations.

  Add its line to `modules/index.md`, and a log line.
- [x] Push `phase-1c` and open the PR, after the user approves that push (R13.3).
  The PR's CI must be green, including the network/golden job and the Python
  3.14 jobs.
- [x] Ask the user to review slice 1C before slice 1D.

---

## Slice 1D - Plots, validation, docs, release

Expanded 2026-09-28 with superpowers:writing-plans (rules.md R1.2a), from three
research passes and a prototype. The prototype was a scratch copy of the repo
at 37da645 with matplotlib 3.11.2 and asv 0.6.6 added. Every source, test,
docs and benchmark file below was written there and passed ruff,
`ruff format`, `mypy --strict`, import-linter, pytest with doctests (random
order, 4 xdist workers), and `sphinx-build -W` with all three notebooks
executed against a pooch cache. The `pl` tests also passed on the floor,
matplotlib 3.8.4 on Python 3.12, and the asv suite ran end to end. Mutating
the stacking loop, the NA recoding and one `R equivalent:` line each made the
tests fail. Implementers still re-check every API they call (R2.2).

**Goal:** close Phase 1. phyloseq's plots draw the results `tl` and `pp` store;
a generated Coming-from-R table is checked against every docstring; executed
tutorials redo phyloseq's analysis vignette; asv baselines, a coverage gate and
a refreshed knowledge bundle follow; biotapy 0.1.0 goes to PyPI.

**Architecture:**
- A new top-layer package `pl` reads slots and computes nothing. Each function
  draws on the `ax` it is given, or on a new pyplot figure, and returns it.
  matplotlib is imported inside the functions, so `import biotapy` stays as
  fast as before.
- A local Sphinx extension, `docs/extensions/coming_from_r.py`, writes the
  Coming-from-R table at `builder-inited` from each public docstring's
  `R equivalent:` line plus `docs/_data/r_idioms.toml`. The same module is the
  parser the docstring contract test uses.
- Notebooks run on every docs build (`nb_execution_mode = "cache"`), and a new
  CI job builds the docs with the network job's pooch cache.
- `benchmarks/` holds an asv suite on a synthetic 5,000 x 50,000 table; its
  results are recorded in `docs/performance.md`.

**Tech stack:**
- matplotlib `>=3.8` (runtime) and asv `>=0.6.6` (dev), both approved
  2026-09-27.
- Stdlib `tomllib`; no numpydoc, no pyyaml in the doc group.
- myst-nb 1.4.0 with jupyter-cache 1.0.1: myst-nb depends on jupyter-cache
  directly, so the doc group already has it.
- Sphinx 9.1 and sphinx-autodoc-typehints 3.13.7. The latter resolves the
  `TYPE_CHECKING`-only `Axes` import, and the matplotlib intersphinx entry
  links it.

**Execution order:** **1.18 -> 1.19 -> 1.19b -> 1.20 -> 1.21 -> 1.22 -> Checkpoint D -> 1.23.**
- 1.19's contract test covers every public subpackage, `pl` included, so `pl`
  comes first.
- 1.22 runs after every code task: any later change under a concept's `paths`
  makes its refreshed `commit` stale again.
- Checkpoint D re-checks staleness after its fix pass, and 1.23 after the
  version bump.

### Slice 1D design
- **Where the code goes.**

  | File | Holds |
  |---|---|
  | `pl/_common.py` | helpers the pl topic files share: `new_axes`, `table`, `obs_groups`, `groups`, `scatter`, `label_ticks` (same subpackage, module-boundaries rule 1) |
  | `pl/_abundance.py` | `bar`, `heatmap` |
  | `pl/_richness.py` | `richness` |
  | `pl/_ordination.py` | `ordination`, `scree` |
  | `docs/extensions/coming_from_r.py` | `r_equivalents` (the parser), `public_functions`, `rows`, `render`, and the `builder-inited` hook |
  | `docs/_data/r_idioms.toml` | R calls that are plain AnnData/TreeData code, and the "not in 0.1" rows |
  | `benchmarks/asv.conf.json`, `benchmarks/benchmarks/{_data,pp,tl}.py` | the asv suite and its synthetic table |

- **What `pl` reads.** A missing slot raises `KeyError` naming the call that
  writes it.

  | Function | Reads | Missing names |
  |---|---|---|
  | `pl.bar(adata, fill, *, x=None, layer=None, ax=None)` | `X` or `layers[layer]`; `fill` in `var` or `obs`; `x` in `obs` | `bt.pp.relative(adata)` for `layer="relative"` |
  | `pl.heatmap(adata, *, layer=None, ax=None)` | `X` or `layers[layer]` | as above |
  | `pl.richness(adata, metric, *, x=None, color=None, ax=None)` | `obs["alpha_<metric>"]` | `bt.tl.alpha(adata, metrics=[<metric>], inplace=True)` |
  | `pl.ordination(adata, *, basis="pcoa", components=(1, 2), color=None, ax=None)` | `obsm["X_<basis>"]` and `uns["biotapy"][<basis>]` | `bt.tl.pcoa(adata, inplace=True)` or `bt.tl.nmds(adata, inplace=True)` |
  | `pl.scree(adata, *, ax=None)` | `uns["biotapy"]["pcoa"]["proportion_explained"]` | `bt.tl.pcoa(adata, inplace=True)` |

- **Rulings.** Each gives what was decided, why, and what it costs if wrong.
  1. **`pl.bar` sums the features of each `fill` group before drawing.**
     phyloseq stacks one black-outlined rectangle per feature and sample; on
     GlobalPatterns that is 19,216 x 26 = 499,616 rectangles. Summing gives the
     same bar and segment heights without the per-feature outlines.
     - Measured: phylum (67 groups) draws 1,742 rectangles in 0.51 s plus
       0.40 s to save. Per-feature stacking took 17.5 s plus 4.5 s for 2,000
       of the 19,216 OTUs, about 170 s plus 45 s for all of them.
     - The docstring states the difference.
     - Cost: a user cannot see individual features inside a segment; for that
       they subset to the features they want first.
  2. **`fill` may name a `var` column (a rank) or an `obs` column; `x` names an
     `obs` column that groups samples into bars.** The vignette's two
     `plot_bar` calls use both (`x = "SeqTech", fill = "Enterotype"`), which
     makes this a real use case (R2.3). A name in both or neither of `var` and
     `obs` raises `KeyError`.
     - `facet_grid` has no parameter: the vignette draws one subset per axes of
       `plt.subplots`.
     - Cost: two lookups to document; a name in both needs renaming.
  3. **Scree is its own function, `pl.scree`, mapping `phyloseq::plot_scree`,**
     rather than a mode of `pl.ordination`. One function is one verb (R3.1),
     and phyloseq keeps `plot_scree` separate too (`plot_ordination(type =
     "scree")` only calls it).
     - This adds a fifth `pl` function and a 32nd row to the table.
     - It plots the stored `proportion_explained`, which is phyloseq's
       `x / sum(x)`, for the axes `tl.pcoa` kept (10 by default). phyloseq
       plots every axis.
     - Cost: one more public name.
  4. **Parameter names follow phyloseq**: `x`, `fill`, `color`, and `metric`
     for one measure. `components=(1, 2)` stands for phyloseq's `axes = 1:2`,
     since `axes`/`ax` would collide in matplotlib code. `basis` is a
     `Literal["pcoa", "nmds"]`.
     - `metric` is `str`, like `tl.pcoa(distance=)`: it names a stored column,
       and `tl`'s `AlphaMetric` is private to `tl` (module-boundaries rule 1).
     - Cost: none known.
  5. **`pl.richness` draws one metric per call.** phyloseq facets every
     measure, but a `pl` function returns one `Axes` (pure-by-default). The
     vignette loops over `plt.subplots(1, 4)`.
     - There are no standard-error bars, because `tl.alpha` stores no SE.
     - Samples with a NaN value are left out, as `geom_point(na.rm = TRUE)`.
     - Cost: a loop in user code for a multi-metric figure.
  6. **Groups.** A categorical column keeps its category order; other values are
     sorted. Missing values form a last group, `NA`, drawn grey (ggplot2's
     default).
     - A numeric, non-bool column raises `TypeError` with an
       `.astype("category")` hint, as `tl.permanova` does. enterotype's `Age`
       would otherwise give 31 legend entries (30 ages and `NA`).
     - Colours: tab10, then tab20 up to 20 groups, then `turbo` resampled.
     - Cost: numeric colouring (a gradient) is not in 0.1.
  7. **`ax=None` draws on `plt.figure().add_subplot()`**, with pyplot imported
     inside `_common.py:new_axes`. A figure made without pyplot is not shown by
     Jupyter's inline backend, nor by `plt.show()`.
     - `import biotapy` never imports matplotlib: it uses a `TYPE_CHECKING`
       import, and `tests/pl/test_init.py` pins it.
     - Cost: the first plot pays about 0.44 s for pyplot.
  8. **`pl.heatmap` keeps `obs`/`var` order** (controller ruling).
     - It uses phyloseq's colours and log scale: `#000033` to `#66CCFF`, zeros
       black through `LogNorm` and the colormap's "bad" colour.
     - Tick labels are drawn only up to 250 names per axis, phyloseq's
       `max.label`; the same cap applies to `pl.bar` and `pl.richness`.
     - A table with no positive value raises `ValueError`. Without that check,
       matplotlib fails inside `colorbar` with `Invalid vmin or vmax`.
     - The guide and the vignette show how to sort by NMDS angle first.
     - Cost: a `pl.heatmap` alone does not reproduce phyloseq's picture.
  9. **An NMDS plot notes the stored stress in its lower right corner**
     (`uns["biotapy"]["nmds"]["stress"]`), which phyloseq does not. PCoA axis
     labels append `[xx.x%]` from `proportion_explained`, like phyloseq's
     `round(100 * Relative_eig, 1)`; NMDS labels carry none.
     - `ordination` needs both the `obsm` coordinates and the `uns` summary.
       `tl.pcoa`/`tl.nmds` write them together, and `feature_subset` drops
       them together.
     - Cost: the note must be removed by hand (`ax.texts[0].remove()`) when
       unwanted.
  10. **No label, title, colour-map or `plot_kwargs` options** (R2.3). The
      function returns the `Axes`, so matplotlib relabels, retitles and
      restyles; the vignette shows `set_xticks(..., labels=...)` in place of
      phyloseq's `sample.label`/`taxa.label`.
      - Cost: one matplotlib line per relabel.
  11. **The Coming-from-R page is committed (`docs/coming_from_r.md`) and
      `{include}`s a generated fragment, `docs/generated/coming_from_r_table.md`.**
      `docs/generated/` is already git-ignored and left out of the sdist.
      - `conf.py` excludes the fragment from the source list, so `-W` does not
        fail on "document isn't included in any toctree".
      - The page's prose stays editable; only the table is generated (R8.4).
      - The include does not depend on Sphinx 9.1's second `find_files()`
        scan (research gotcha 7), because the fragment is read when the page
        is parsed.
      - Cost: the fragment path is named in two places; a test ties them
        together.
  12. **Idioms are TOML (`docs/_data/r_idioms.toml`, stdlib `tomllib`),** not
      the roadmap's YAML: the doc group has no YAML parser, and pyyaml is not
      added.
      - Values are MyST inserted verbatim, so a "not in 0.1" row is plain text
        and an idiom is a code span.
      - A call mapped by both a docstring and the file raises `ValueError`.
      - Cost: TOML literal strings (`'...'`) are needed wherever a value holds
        double quotes.
  13. **The vignette is a MyST text notebook, `docs/tutorials/phyloseq_analysis.md`,**
      like `quick_tour.md`, not an `.ipynb`. It diffs as plain text, never
      carries committed outputs, and matches the existing tutorial.
      - The exit gate's wording is updated to match (done with this plan).
      - Cost: ruff does not lint its code cells (it lints `*.ipynb`), and a
        user who wants an `.ipynb` converts it with jupytext.
  14. **One execution setting.** `conf.py` sets `nb_execution_mode = "cache"`
      and `nb_execution_raise_on_error = True`, and `quick_tour.md` drops its
      page-level `execution_mode: force` override.
      - CI and Read the Docs start with an empty `docs/_build/.jupyter_cache`,
        so they always execute.
      - Cost: locally, an unchanged notebook keeps cached outputs after a code
        change until `uvx hatch run docs:clean`.
      - The default 30 s per-cell timeout stays. The slowest notebook ran in
        11.4 s in total.
  15. **The vignette computes everything; it reads no golden CSV** (controller
      ruling). Unweighted UniFrac on GlobalPatterns takes 0.32 s, so it does
      not need the R vignette's precomputed file.
  16. **No CI job runs the benchmarks.** Phase 1 has no speed gate, and the
      full suite took 11 minutes and a 2.73 GB peak in the prototype.
      - The `lint` job runs `asv check`, which imports and discovers the suite
        in seconds and needs no machine file. It catches a benchmark broken by
        an API change.
      - Cost: a regression in speed goes unseen until the next manual run.
  17. **asv runs with `HOME` pointed at the git-ignored `.asv/`**
      (`uv run --group dev env HOME="$PWD/../.asv" asv ...`). asv 0.6.6 keeps
      its machine file only at `~/.asv-machine.json` (`asv/machine.py`), and
      the user's rules forbid files outside the project.
      - `uv run` itself keeps the normal `HOME`, so its cache is untouched.
      - Cost: a longer command line, written out on the performance page.
  18. **Coverage gate: `report.fail_under = 90` on total line coverage.**
      coverage.py has no per-function view, so the whole package (public and
      private modules) is the proxy for R11.6.
      - The CI `test` job's `cov-report` step (hatch's `coverage report`) then
        fails below 90%.
      - Per-file numbers are reviewed at Checkpoint D.
      - Cost: one weak public module can hide under a high total, as
        `datasets/_remote.py` (89% without network) does today.
  19. **"Uncheckable" concepts stay uncheckable.** The 12 are the 10 decisions
      plus `phase-5-beyond` and `spec-review`, which have no `paths:` key, by
      design: an ADR is not tied to a code path.
      - The 1.22 gate is 0 stale: 21 current, 0 stale, 12 uncheckable.
      - Cost: a decision made wrong by code is caught only in review.
  20. **matplotlib floor `>=3.8`.** 3.8.0 is the first release with type
      information (`py.typed`), which `mypy --strict` needs for `Axes`, and
      the first with CPython 3.12 wheels.
      - Every call `pl` makes exists there: `Colormap.with_extremes` (3.4),
        the `matplotlib.colormaps` registry and `set_xticks(labels=)` (3.5).
      - 3.8.0 to 3.8.3 pin `numpy<2`, so with NumPy 2 the resolver takes
        3.8.4 or later. The tests passed on 3.8.4 with Python 3.12.
      - Cost: an older matplotlib cannot be installed alongside biotapy.
  21. **No golden tests for `pl`** (controller ruling), although R11.2 asks
      for one when an R equivalent exists: plots draw numbers that `tl`'s golden
      tests already check.
      - Task 1.18 records the exception in r-golden-parity, add-a-function and
        rules.md R11.2. It edits rules.md only because approving this plan
        approves that edit.
      - Cost: none beyond the wording.
  22. **`optional-heavy-dependencies` is corrected where the code disagrees**
      (R0.2), in Task 1.18, which edits that concept anyway. matplotlib was
      listed as added in "Phase 0-1" but was not a dependency. The concept also
      says CI runs an all-extras job, but no extra exists yet; only
      `import-without-extras` does.
- **Research facts the tasks rely on.** All were measured on the prototype;
  re-check each (R2.2).
  - **Import time** before this slice, on 37da645:
    - `uv run python -X importtime -c "import biotapy" 2>&1 | tail -1`:
      1.21-1.31 s cumulative over three runs;
    - wall clock: 1.40-1.51 s;
    - `matplotlib` is not in `sys.modules` afterwards.
  - **GlobalPatterns facts.**
    - `filter_features(min_total=1)` keeps 18,988 taxa, phyloseq's
      `taxa_sums > 0`.
    - 3 OTUs have no phylum.
    - The Crenarchaeota subset is 106 OTUs with no empty sample.
    - Unweighted UniFrac takes 0.32 s, PCoA at 25 axes 0.02 s, and NMDS 1.45 s
      (stress 0.145).
    - PC1-PC3 explain 43% of the unweighted UniFrac PCoA.
  - **enterotype.** `Enterotype` is categorical (`1`, `2`, `3`) with 9 missing
    values; `Age` is float64 with 243 missing; 2 features have no genus. The
    10 most abundant features start with `-1`, whose genus is missing.
  - **Import time after `pl`:** `biotapy.pl` adds 0.47-0.54 ms (its
    cumulative `-X importtime` line). Totals of 1.19-1.40 s over eight runs
    differ from the baseline only by machine noise, and no `matplotlib` line
    appears.
  - **Coverage** with `pl` in place: 560 tests passed; `TOTAL` 99%; each `pl`
    file 100%; `datasets/_remote.py` 89%. With `fail_under` above the total,
    `coverage report` exits 2.
  - **Docs build.** `sphinx-build -W` took 34 s from clean. The notebooks
    executed in 5.1-5.9 s (getting started), 11.4-13.5 s (vignette) and 2.3 s
    (quick tour). The generated table had 48 rows, and `Axes` linked to
    matplotlib.org. A failing cell made the build exit 2 with
    `WARNING: Executing notebook failed: CellExecutionError [mystnb.exec]`.
    Changing one notebook re-executed only that one (cache mode).
  - **asv on 5,000 x 50,000** (16-core laptop, 62 GB; `asv run --python=same`,
    660 s in all):
    - generating the table: 1.7 s, 5,000,000 non-zeros, a 100,006-node tree,
      and a smallest sample depth of 43,495;
    - `relative` 897 ms, `tax_glom` to genus 2.90 s, `rarefy` 5.44 s;
    - `alpha` (four metrics) 1.39 s, with a 478 MB process peak;
    - `faith_pd` 43.5 s. Under cProfile, 97% of its time went to scikit-bio
      re-indexing and re-validating the tree for each of the 250 chunks
      (`vectorize_counts_and_tree`, `_validate_taxa_and_tree`), against
      0.46 s for `get_skbio_tree` once. This measures the per-chunk
      re-indexing Checkpoint C deferred to this slice;
    - `beta`: 5.62 s, 16.2 s and 54.2 s at 500, 1,000 and 2,000 samples
      (quadratic, so about 6 minutes at 5,000). The 5,000-sample run peaked at
      2.73 GB: X densified once (2.0 GB) plus the distances, as the
      docstring states.
    - `asv check` took 3.6 s with an empty `HOME`, wrote nothing there, and
      exited 1 when a benchmark module failed to import.
  - **conda-forge.** `treedata` is on neither conda-forge nor bioconda
    (api.anaconda.org, 404 on both). Every other runtime dependency is on
    conda-forge.

### Slice 1D global constraints (in addition to the Phase 1 list)
- **Dependencies.**
  - Runtime `matplotlib>=3.8` and dev `asv>=0.6.6`, both approved 2026-09-27,
    are added in 1.18 and 1.21.
  - No other dependency: numpydoc, pyyaml (doc group) and grayskull are not
    added (R9.1).
- **Lazy matplotlib.** `pl` modules import matplotlib only inside functions.
  `Axes` annotations are strings, resolved by an `if TYPE_CHECKING:` import.
  `import biotapy` must not load matplotlib, and `tests/pl/test_init.py` pins
  this. matplotlib is a runtime dependency, not an extra, so it does not go
  through `import_optional`.
- **`pl` computes nothing.**
  - It reads the slots in the design table, and a missing slot is a `KeyError`
    naming the writing call.
  - It writes nothing to the AnnData, and every test asserts
    `assert_unchanged`.
  - `ax: "Axes | None" = None` is the last keyword-only parameter, and every
    function returns the `Axes`.
  - No `plot_kwargs`, no label or title options (ruling 10).
- **Tests draw headless.**
  - `tests/conftest.py` calls `matplotlib.use("Agg")`. `tests/` is collected
    before `src/biotapy`, so the doctests get it too; CI also sets
    `MPLBACKEND=agg`.
  - Tests draw on the `ax` fixture (`Figure().add_subplot()`, not tracked by
    pyplot). Only the `ax=None` test uses pyplot, and it closes its figure.
  - Read-back: `ax.patches` (bars), `ax.collections[i].get_offsets()`
    (points), `ax.images[0].get_array()` (heatmap), `ax.texts` (stress note),
    and `ax.get_legend().get_texts()`.
- **Size limits (R5) and types.** Function bodies stay within 30 statements,
  complexity 8, 8 branches, 6 arguments (3 positional) and 4 returns, with 120
  characters per line. `mypy --strict` passes on `src/biotapy`. The prototype
  passed all of these as written below; keep it that way when adapting.
- **Docs.**
  - `nitpicky` stays on, and `-W` must pass.
  - Guide pages are plain markdown with `python` blocks, not executed.
    Tutorials are MyST notebooks, executed at build time.
- **Commits.**
  - Stage explicit paths only. Never stage `.claude/`, `.superpowers/`,
    `.worktrees/` or `notebooks/`.
  - `uv.lock` is git-ignored (`/uv.lock` in `.gitignore`, "resolve fresh in
    CI"). `git add uv.lock` fails, so it is never staged; slice 1C's "stage
    `uv.lock`" instruction did not apply.
- **Every task's last commit** also stages `.knowledge/roadmap/phase-1-core.md`
  with that task's boxes ticked, and `.knowledge/log.md` with its dated line
  (R12.4).
- **Approval before every outward action** (R13.3): push, open a PR, merge,
  tag, create a GitHub release, and publish to PyPI. Each needs the user's
  explicit approval at that step.

### Slice 1D review focus
1. **Missing values in the columns plots group by.** enterotype's `Enterotype`
   (9 missing) and GlobalPatterns' `phylum` (3 OTUs) are examples. Expected:
   a last group `NA`, drawn grey; bars still add up to the sample totals; NaN
   alpha values (all-zero samples) are left out, not drawn at 0. Tests:
   - 1.18 `test_bar_missing_rank_is_its_own_group_and_bars_keep_every_read`;
   - 1.18 `test_bar_missing_obs_value_is_its_own_group`;
   - 1.18 `test_bar_heights_add_up_to_sample_totals` (Hypothesis, with `None`
     labels);
   - 1.18 `test_richness_color_draws_one_scatter_per_group`;
   - 1.18 `test_richness_nan_values_are_left_out`.
2. **Stale or half-present results.** A plot after `pp.filter_features`
   dropped the ordination, or `obsm` present without its `uns` summary, must
   raise a `KeyError` naming the `tl` call, not an `IndexError` from numpy.
   Tests:
   - 1.18 `test_ordination_is_gone_after_a_feature_change`;
   - 1.18 `test_ordination_missing_names_the_call`;
   - 1.18 `test_scree_missing_pcoa_names_the_call`;
   - 1.18 `test_bar_missing_layer_names_the_call`.
3. **Numeric columns as groups** (enterotype `Age`, `obs["alpha_*"]`) raise
   `TypeError` with the `.astype("category")` hint, as `tl.permanova` does.
   Tests: 1.18 `test_bar_numeric_x_raises` and `test_richness_bad_columns_raise`.
4. **Tables too big or too empty to draw.** Past 250 names an axis gets no
   tick labels; an all-zero table gets biotapy's `ValueError`, not
   matplotlib's `Invalid vmin or vmax` from inside `colorbar`. Tests: 1.18
   `test_heatmap_labels_at_most_250_names` and
   `test_heatmap_nothing_positive_raises`.
5. **A notebook user calling `bt.pl.bar(tdata, "phylum")` with no `ax`.**
   The figure must show (pyplot-managed), while `import biotapy` stays free of
   matplotlib. Tests: 1.18
   `test_bar_without_ax_draws_on_a_new_pyplot_figure` and
   `test_import_biotapy_does_not_load_matplotlib`; the executed tutorials in
   1.20 call `pl.ordination`, `pl.richness` and `pl.scree` without `ax`.

### Task 1.18: `pl.bar`, `pl.heatmap`, `pl.richness`, `pl.ordination`, `pl.scree`

**Files:**
- Create:
  - `src/biotapy/pl/__init__.py`, `_common.py`, `_abundance.py`, `_richness.py`, `_ordination.py`
  - `tests/pl/conftest.py`, `test_abundance.py`, `test_richness.py`, `test_ordination.py`, `test_init.py`
  - `docs/guide/plotting.md`
- Modify:
  - `pyproject.toml` (matplotlib)
  - `src/biotapy/__init__.py`
  - `tests/conftest.py` (the Agg backend)
  - `docs/conf.py` (matplotlib intersphinx), `docs/api.md`, `docs/guide/index.md`
  - `rules.md` (R11.2, ruling 21)
  - `.knowledge/decisions/optional-heavy-dependencies.md`
  - `.knowledge/contracts/r-golden-parity.md`
  - `.knowledge/playbooks/add-a-function.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:**
  - the slots written by `tl.alpha`, `tl.pcoa` and `tl.nmds` with
    `inplace=True`, and by `pp.relative`;
  - `_core.as_csr` and `_core.sum_by`.
- **Produces** (1.19 reads each docstring's `R equivalent:` line; 1.20's
  notebooks call each function):
  - `bt.pl.bar(adata: AnnData, fill: str, *, x: str | None = None, layer: str | None = None, ax: Axes | None = None) -> Axes`
  - `bt.pl.heatmap(adata: AnnData, *, layer: str | None = None, ax: Axes | None = None) -> Axes`
  - `bt.pl.richness(adata: AnnData, metric: str, *, x: str | None = None, color: str | None = None, ax: Axes | None = None) -> Axes`
  - `bt.pl.ordination(adata: AnnData, *, basis: Literal["pcoa", "nmds"] = "pcoa", components: tuple[int, int] = (1, 2), color: str | None = None, ax: Axes | None = None) -> Axes`
  - `bt.pl.scree(adata: AnnData, *, ax: Axes | None = None) -> Axes`
  - `Axes` is `matplotlib.axes.Axes`. The R equivalents are
    ``phyloseq::plot_bar``, ``plot_heatmap``, ``plot_richness``,
    ``plot_ordination`` and ``plot_scree``, and every `Guide:` line is
    `:doc:`/guide/plotting``.

- [x] **Step 1: Baseline and dependency.**
  - Before any change, record `uv run python -X importtime -c "import biotapy" 2>&1 | tail -1`
    three times. The prototype measured 1.21-1.31 s cumulative on 37da645
    (wall clock 1.40-1.51 s).
  - In `[project] dependencies`, add `"matplotlib>=3.8",` after
    `"biom-format>=2.1.16",`; pyproject-fmt keeps that order.
  - Run `uv sync --group dev --group test --group doc`. It resolves
    matplotlib 3.11.2 or later, and with it contourpy, cycler, fonttools,
    kiwisolver, pillow and pyparsing.
  - Confirm (R2.2) that
    `uv run python -c "import matplotlib, pathlib; print(matplotlib.__version__, (pathlib.Path(matplotlib.__file__).parent / 'py.typed').exists())"`
    prints the version and `True`: mypy needs the marker to type `Axes`.
  - No mypy override is needed, unlike scikit-bio and scikit-learn.
- [x] **Step 2: Failing tests.**
  - `tests/conftest.py`:
    - add `import matplotlib` to the imports, after `from collections.abc import Callable`;
    - after `from biotapy._core import as_csr`, add:
      ```python

      # Headless and identical to CI's MPLBACKEND=agg; tests/ is collected before the src/ doctests, so they get it too.
      matplotlib.use("Agg")
      ```
  - `tests/pl/conftest.py`:
    ```python
    import pytest
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


    @pytest.fixture
    def ax() -> Axes:
        """Axes on a figure pyplot does not track, so tests leave no open figures behind."""
        return Figure().add_subplot()
    ```
  - `tests/pl/test_abundance.py`:
    ```python
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import pytest
    import scipy.sparse as sp
    from hypothesis import given
    from hypothesis import strategies as st
    from hypothesis.extra.numpy import arrays
    from matplotlib.figure import Figure

    import biotapy as bt


    def _adata(dense, var=None) -> ad.AnnData:
        return ad.AnnData(
            X=sp.csr_matrix(dense),
            obs=pd.DataFrame(index=[f"s{i}" for i in range(dense.shape[0])]),
            var=pd.DataFrame(var, index=[f"f{i}" for i in range(dense.shape[1])]),
        )


    def _heights(ax, n_bars):
        # ax.bar adds one Rectangle per bar, one group after another: reshape to bars x groups.
        return np.array([patch.get_height() for patch in ax.patches]).reshape(-1, n_bars).T


    def _legend(ax):
        return [text.get_text() for text in ax.get_legend().get_texts()]


    def _by_phylum(tdata):
        frame = pd.DataFrame(tdata.X.toarray(), columns=tdata.var["phylum"].to_numpy())
        return frame.T.groupby(level=0).sum().T  # samples x phyla, phyla sorted


    def test_bar_sums_each_fill_group_per_sample(ax):
        tdata = bt.datasets.toy()
        bt.pl.bar(tdata, "phylum", ax=ax)
        expected = _by_phylum(tdata)
        assert _legend(ax) == list(expected.columns) == ["Bacteroidota", "Firmicutes", "Proteobacteria"]
        np.testing.assert_array_equal(_heights(ax, 6), expected.to_numpy())
        bottoms = np.array([patch.get_y() for patch in ax.patches]).reshape(-1, 6).T
        np.testing.assert_array_equal(bottoms[:, 1:], np.cumsum(expected.to_numpy(), axis=1)[:, :-1])
        assert [label.get_text() for label in ax.get_xticklabels()] == list(tdata.obs_names)


    def test_bar_x_sums_the_samples_of_each_group(ax):
        tdata = bt.datasets.toy()
        bt.pl.bar(tdata, "phylum", x="group", ax=ax)
        expected = _by_phylum(tdata).groupby(tdata.obs["group"].to_numpy()).sum()
        np.testing.assert_array_equal(_heights(ax, 2), expected.to_numpy())
        assert [label.get_text() for label in ax.get_xticklabels()] == ["A", "B"]


    def test_bar_fill_from_obs_colours_whole_samples(ax):
        tdata = bt.datasets.toy()
        bt.pl.bar(tdata, "group", ax=ax)
        totals = tdata.X.sum(axis=1).A1
        np.testing.assert_array_equal(_heights(ax, 6), np.c_[np.r_[totals[:3], 0, 0, 0], np.r_[0, 0, 0, totals[3:]]])
        assert _legend(ax) == ["A", "B"]


    def test_bar_missing_rank_is_its_own_group_and_bars_keep_every_read(ax):
        tdata = bt.datasets.toy()  # f8 has no genus
        bt.pl.bar(tdata, "genus", ax=ax)
        assert _legend(ax)[-1] == "NA"
        np.testing.assert_array_equal(_heights(ax, 6)[:, -1], tdata.X[:, 7].toarray().ravel())
        np.testing.assert_array_equal(_heights(ax, 6).sum(axis=1), tdata.X.sum(axis=1).A1)


    def test_bar_missing_obs_value_is_its_own_group(ax):
        tdata = bt.datasets.toy()
        tdata.obs["site"] = pd.Categorical(["u", None, "v", "u", "v", None])
        bt.pl.bar(tdata, "phylum", x="site", ax=ax)
        assert [label.get_text() for label in ax.get_xticklabels()] == ["u", "v", "NA"]
        np.testing.assert_array_equal(_heights(ax, 3).sum(), tdata.X.sum())


    def test_bar_layer_plots_relative_abundance(ax):
        bt.pl.bar(bt.pp.relative(bt.datasets.toy()), "phylum", layer="relative", ax=ax)
        np.testing.assert_allclose(_heights(ax, 6).sum(axis=1), 1.0)
        assert ax.get_ylabel() == "relative"


    def test_bar_all_zero_sample_and_feature(ax):
        dense = np.array([[0, 0, 0], [1, 0, 2], [3, 0, 0]])
        bt.pl.bar(_adata(dense, {"phylum": ["a", "b", "a"]}), "phylum", ax=ax)
        np.testing.assert_array_equal(_heights(ax, 3), [[0, 0], [3, 0], [3, 0]])


    def test_bar_single_sample(ax):
        bt.pl.bar(_adata(np.array([[4, 0, 1]]), {"phylum": ["a", "b", "a"]}), "phylum", ax=ax)
        np.testing.assert_array_equal(_heights(ax, 1), [[5, 0]])


    def test_bar_input_unchanged(assert_unchanged, ax):
        tdata = bt.pp.relative(bt.datasets.toy())
        before = tdata.copy()
        bt.pl.bar(tdata, "genus", x="group", layer="relative", ax=ax)
        assert_unchanged(before, tdata)


    def test_bar_without_ax_draws_on_a_new_pyplot_figure():
        ax = bt.pl.bar(bt.datasets.toy(), "phylum")
        assert ax.figure.number in plt.get_fignums()
        plt.close(ax.figure)


    def test_bar_missing_layer_names_the_call(ax):
        with pytest.raises(KeyError, match=r"run bt.pp.relative\(adata\) first"):
            bt.pl.bar(bt.datasets.toy(), "phylum", layer="relative", ax=ax)


    def test_bar_fill_must_be_in_exactly_one_of_var_and_obs(ax):
        with pytest.raises(KeyError, match="fill='site'"):
            bt.pl.bar(bt.datasets.toy(), "site", ax=ax)
        tdata = bt.datasets.toy()
        tdata.obs["phylum"] = "x"
        with pytest.raises(KeyError, match="exactly one of var and obs"):
            bt.pl.bar(tdata, "phylum", ax=ax)


    def test_bar_numeric_x_raises(ax):
        tdata = bt.datasets.toy()
        tdata.obs["depth"] = tdata.X.sum(axis=1).A1
        with pytest.raises(TypeError, match=r"x='depth' is a numeric column"):
            bt.pl.bar(tdata, "phylum", x="depth", ax=ax)
        with pytest.raises(KeyError, match="x='site' is not a column of obs"):
            bt.pl.bar(tdata, "phylum", x="site", ax=ax)


    @given(
        arrays(np.int64, st.tuples(st.integers(1, 6), st.integers(1, 6)), elements=st.integers(0, 20)),
        st.lists(st.sampled_from(["a", "b", None]), min_size=6, max_size=6),
    )
    def test_bar_heights_add_up_to_sample_totals(dense, labels):
        ax = Figure().add_subplot()
        bt.pl.bar(_adata(dense, {"phylum": labels[: dense.shape[1]]}), "phylum", ax=ax)
        np.testing.assert_array_equal(_heights(ax, dense.shape[0]).sum(axis=1), dense.sum(axis=1))


    def test_heatmap_draws_the_table_features_by_samples(ax):
        tdata = bt.datasets.toy()
        image = bt.pl.heatmap(tdata, ax=ax).images[0]
        np.testing.assert_array_equal(image.get_array(), tdata.X.toarray().T)
        assert [label.get_text() for label in ax.get_xticklabels()] == list(tdata.obs_names)
        assert [label.get_text() for label in ax.get_yticklabels()] == list(tdata.var_names)
        assert len(ax.figure.axes) == 2  # the heatmap and its colour bar


    def test_heatmap_zeros_are_drawn_black_on_a_log_scale(ax):
        image = bt.pl.heatmap(bt.datasets.toy(), ax=ax).images[0]
        assert np.ma.is_masked(image.norm(np.array([0.0])))  # LogNorm masks zeros; the colormap draws them "bad"
        assert image.cmap.get_bad().tolist() == [0.0, 0.0, 0.0, 1.0]


    def test_heatmap_labels_at_most_250_names(ax):
        bt.pl.heatmap(_adata(np.ones((1, 251))), ax=ax)
        assert ax.get_yticks().size == 0
        bt.pl.heatmap(_adata(np.ones((1, 250))), ax=Figure().add_subplot())


    def test_heatmap_layer_and_all_zero_sample(ax):
        tdata = bt.pp.relative(bt.datasets.toy())
        image = bt.pl.heatmap(tdata, layer="relative", ax=ax).images[0]
        np.testing.assert_allclose(image.get_array().sum(axis=0), 1.0)
        dense = np.array([[0, 0], [1, 2]])
        np.testing.assert_array_equal(
            bt.pl.heatmap(_adata(dense), ax=Figure().add_subplot()).images[0].get_array(), dense.T
        )


    def test_heatmap_nothing_positive_raises(ax):
        with pytest.raises(ValueError, match="no positive value"):
            bt.pl.heatmap(_adata(np.zeros((2, 3))), ax=ax)


    def test_heatmap_input_unchanged(assert_unchanged, ax):
        tdata = bt.datasets.toy()
        before = tdata.copy()
        bt.pl.heatmap(tdata, ax=ax)
        assert_unchanged(before, tdata)
    ```
  - `tests/pl/test_richness.py`:
    ```python
    import numpy as np
    import pandas as pd
    import pytest
    import scipy.sparse as sp

    import biotapy as bt


    def _toy_with(*metrics):
        tdata = bt.datasets.toy()
        bt.tl.alpha(tdata, metrics=list(metrics), inplace=True)
        return tdata


    def _points(ax):
        return np.concatenate([collection.get_offsets() for collection in ax.collections])


    def test_richness_one_point_per_sample(ax):
        tdata = _toy_with("shannon")
        bt.pl.richness(tdata, "shannon", ax=ax)
        np.testing.assert_array_equal(_points(ax), np.c_[np.arange(6), tdata.obs["alpha_shannon"]])
        assert [label.get_text() for label in ax.get_xticklabels()] == list(tdata.obs_names)
        assert ax.get_ylabel() == "shannon"


    def test_richness_x_places_samples_by_group(ax):
        tdata = _toy_with("observed_features")
        bt.pl.richness(tdata, "observed_features", x="group", ax=ax)
        np.testing.assert_array_equal(_points(ax)[:, 0], [0, 0, 0, 1, 1, 1])
        assert [label.get_text() for label in ax.get_xticklabels()] == ["A", "B"]


    def test_richness_color_draws_one_scatter_per_group(ax):
        tdata = _toy_with("shannon")
        tdata.obs["site"] = pd.Categorical(["u", "v", None, "u", "v", "u"])
        bt.pl.richness(tdata, "shannon", x="group", color="site", ax=ax)
        assert [text.get_text() for text in ax.get_legend().get_texts()] == ["u", "v", "NA"]
        assert [len(collection.get_offsets()) for collection in ax.collections] == [3, 2, 1]
        np.testing.assert_array_equal(np.sort(_points(ax)[:, 1]), np.sort(tdata.obs["alpha_shannon"]))


    def test_richness_nan_values_are_left_out(ax):
        tdata = bt.datasets.toy()
        dense = tdata.X.toarray()
        dense[0] = 0
        tdata.X = sp.csr_matrix(dense)
        bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)  # s1 is all-zero: Shannon NaN
        bt.pl.richness(tdata, "shannon", color="group", ax=ax)
        assert len(_points(ax)) == 5 and np.isfinite(_points(ax)).all()


    def test_richness_input_unchanged(assert_unchanged, ax):
        tdata = _toy_with("shannon")
        before = tdata.copy()
        bt.pl.richness(tdata, "shannon", x="group", color="group", ax=ax)
        assert_unchanged(before, tdata)


    def test_richness_missing_metric_names_the_call(ax):
        with pytest.raises(KeyError, match=r"bt.tl.alpha\(adata, metrics=\['chao1'\], inplace=True\)"):
            bt.pl.richness(_toy_with("shannon"), "chao1", ax=ax)


    def test_richness_bad_columns_raise(ax):
        tdata = _toy_with("shannon")
        with pytest.raises(KeyError, match="color='site' is not a column of obs"):
            bt.pl.richness(tdata, "shannon", color="site", ax=ax)
        with pytest.raises(TypeError, match="x='alpha_shannon' is a numeric column"):
            bt.pl.richness(tdata, "shannon", x="alpha_shannon", ax=ax)
    ```
  - `tests/pl/test_ordination.py`:
    ```python
    import numpy as np
    import pytest

    import biotapy as bt


    def _toy_with(*bases):
        tdata = bt.datasets.toy()
        bt.tl.beta(tdata, inplace=True)
        if "pcoa" in bases:
            bt.tl.pcoa(tdata, inplace=True)
        if "nmds" in bases:
            bt.tl.nmds(tdata, seed=0, inplace=True)
        return tdata


    def _points(ax):
        return np.concatenate([collection.get_offsets() for collection in ax.collections])


    def test_ordination_draws_stored_pcoa_axes_with_percentages(ax):
        tdata = _toy_with("pcoa")
        bt.pl.ordination(tdata, ax=ax)
        np.testing.assert_array_equal(_points(ax), tdata.obsm["X_pcoa"][:, :2])
        share = tdata.uns["biotapy"]["pcoa"]["proportion_explained"]
        assert ax.get_xlabel() == f"PC1 [{100 * share[0]:.1f}%]"
        assert ax.get_ylabel() == f"PC2 [{100 * share[1]:.1f}%]"


    def test_ordination_components_pick_the_axes(ax):
        tdata = _toy_with("pcoa")
        bt.pl.ordination(tdata, components=(1, 3), ax=ax)
        np.testing.assert_array_equal(_points(ax), tdata.obsm["X_pcoa"][:, [0, 2]])
        assert ax.get_ylabel().startswith("PC3 [")


    def test_ordination_nmds_has_no_percentages_and_notes_the_stress(ax):
        tdata = _toy_with("nmds")
        bt.pl.ordination(tdata, basis="nmds", color="group", ax=ax)
        np.testing.assert_array_equal(np.sort(_points(ax), axis=0), np.sort(tdata.obsm["X_nmds"], axis=0))
        assert (ax.get_xlabel(), ax.get_ylabel()) == ("NMDS1", "NMDS2")
        assert ax.texts[0].get_text() == f"stress {tdata.uns['biotapy']['nmds']['stress']:.3f}"
        assert [text.get_text() for text in ax.get_legend().get_texts()] == ["A", "B"]


    def test_ordination_input_unchanged(assert_unchanged, ax):
        tdata = _toy_with("pcoa", "nmds")
        before = tdata.copy()
        bt.pl.ordination(tdata, basis="nmds", color="group", ax=ax)
        bt.pl.scree(tdata, ax=ax)
        assert_unchanged(before, tdata)


    def test_ordination_missing_names_the_call(ax):
        with pytest.raises(KeyError, match=r"run bt.tl.pcoa\(adata, inplace=True\) first"):
            bt.pl.ordination(_toy_with(), ax=ax)
        with pytest.raises(KeyError, match=r"run bt.tl.nmds\(adata, inplace=True\) first"):
            bt.pl.ordination(_toy_with("pcoa"), basis="nmds", ax=ax)


    def test_ordination_is_gone_after_a_feature_change(ax):
        tdata = bt.pp.filter_features(_toy_with("pcoa"), min_total=19)
        with pytest.raises(KeyError, match="bt.tl.pcoa"):
            bt.pl.ordination(tdata, ax=ax)


    def test_ordination_bad_basis_or_components_raise(ax):
        tdata = _toy_with("pcoa")
        with pytest.raises(ValueError, match="basis must be one of"):
            bt.pl.ordination(tdata, basis="tsne", ax=ax)
        with pytest.raises(ValueError, match="components must be two axis numbers from 1 to 5"):
            bt.pl.ordination(tdata, components=(1, 6), ax=ax)
        with pytest.raises(ValueError, match="components"):
            bt.pl.ordination(tdata, components=(0, 1), ax=ax)


    def test_scree_draws_the_stored_proportions(ax):
        tdata = _toy_with("pcoa")
        bt.pl.scree(tdata, ax=ax)
        heights = [patch.get_height() for patch in ax.patches]
        np.testing.assert_array_equal(heights, tdata.uns["biotapy"]["pcoa"]["proportion_explained"])
        assert [label.get_text() for label in ax.get_xticklabels()] == ["PC1", "PC2", "PC3", "PC4", "PC5"]


    def test_scree_missing_pcoa_names_the_call(ax):
        with pytest.raises(KeyError, match="bt.tl.pcoa"):
            bt.pl.scree(_toy_with("nmds"), ax=ax)
    ```
  - `tests/pl/test_init.py`:
    ```python
    import subprocess
    import sys


    def test_import_biotapy_does_not_load_matplotlib():
        # pl imports matplotlib inside its functions, so `import biotapy` stays as fast as before (Task 1.18).
        code = "import sys, biotapy; print('matplotlib' in sys.modules)"
        assert (
            subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout.strip()
            == "False"
        )
    ```
  - What the tests pin:
    - `_heights` reads `ax.patches` back as bars x groups. `ax.bar` adds one
      `Rectangle` per bar, one group after another.
    - `test_bar_heights_add_up_to_sample_totals` is the invariant behind
      ruling 1: with `None` labels mixed in, every read stays in its bar.
    - The toy NMDS is degenerate (stress about 0.002), which does not matter
      here: the tests compare plotted points with the stored ones, sorted.
- [x] **Step 3: Run, expect failure.** Run `uv run --group test pytest tests/pl -q`.
  - Every test except `test_import_biotapy_does_not_load_matplotlib` fails
    with `AttributeError: module 'biotapy' has no attribute 'pl'`.
  - That one passes already: it guards the implementation against a
    module-level import.
- [x] **Step 4: Implement.**
  - `src/biotapy/pl/_common.py`:
    ```python
    """Helpers shared by the pl topic files: new axes, groups of a column, colours, the plotted table."""

    from typing import TYPE_CHECKING, Literal, cast

    import numpy as np
    import numpy.typing as npt
    import pandas as pd
    import scipy.sparse as sp
    from anndata import AnnData

    from biotapy._core import as_csr

    if TYPE_CHECKING:
        from matplotlib.axes import Axes

    RGBA = tuple[float, float, float, float]
    # One integer code per sample or feature, the group labels, and one colour per label.
    Groups = tuple[npt.NDArray[np.intp], list[str], list[RGBA]]

    # ggplot2's defaults, which phyloseq's plots inherit: missing values are their own group, "NA", drawn grey50 and last.
    MISSING_LABEL = "NA"
    MISSING_COLOR = "#7f7f7f"
    # phyloseq's plot_heatmap(max.label = 250): past this many names, tick labels overlap into a solid block.
    MAX_LABELS = 250
    # The call that writes each layer pl can read, for the error raised when it is missing.
    _LAYER_WRITTEN_BY = {"relative": "bt.pp.relative(adata)"}


    def new_axes(ax: "Axes | None") -> "Axes":
        """``ax``, or the axes of a new pyplot figure, which a notebook displays and ``plt.show()`` shows."""
        if ax is not None:
            return ax
        # Imported here: pyplot costs about 0.4 s, paid by the first plot and never by `import biotapy`.
        import matplotlib.pyplot as plt

        return plt.figure().add_subplot()


    def table(adata: AnnData, layer: str | None) -> sp.csr_matrix:
        """``X``, or ``layers[layer]``; a missing layer names the call that writes it."""
        if layer is None:
            return as_csr(adata.X)
        if layer not in adata.layers:
            call = _LAYER_WRITTEN_BY.get(layer, "the function that writes it")
            msg = f"layer={layer!r}: no layers[{layer!r}]; run {call} first"
            raise KeyError(msg)
        return as_csr(adata.layers[layer])


    def obs_groups(adata: AnnData, column: str, *, arg: str) -> Groups:
        """Groups of ``obs[column]``; ``arg`` names the parameter in errors."""
        if column not in adata.obs.columns:
            msg = f"{arg}={column!r} is not a column of obs"
            raise KeyError(msg)
        # anndata types obs columns as Series | DataArray (its lazy variant); the data model guarantees a Series.
        return groups(cast("pd.Series", adata.obs[column]), arg=arg)


    def groups(values: pd.Series, *, arg: str) -> Groups:
        """Codes, labels and colours: a categorical's order, else sorted values, with missing values last as NA."""
        if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
            msg = (
                f"{arg}={values.name!r} is a numeric column ({values.dtype}); plots group by category, "
                'so convert it with .astype("category") for one group per value'
            )
            raise TypeError(msg)
        categorical = pd.Categorical(values).remove_unused_categories()
        codes = categorical.codes.astype(np.intp)
        labels = [str(category) for category in categorical.categories]
        missing = bool((codes < 0).any())
        if missing:
            codes = np.where(codes < 0, len(labels), codes)
            labels.append(MISSING_LABEL)
        return codes, labels, _colors(len(labels) - missing, missing=missing)


    def _colors(n: int, *, missing: bool) -> list[RGBA]:
        from matplotlib import colormaps
        from matplotlib.colors import to_rgba

        # tab10 and tab20 as long as they last; past 20 groups, n evenly spaced colours, as ggplot2's hue scale.
        cmap = colormaps["tab10"] if n <= 10 else colormaps["tab20"] if n <= 20 else colormaps["turbo"].resampled(n)
        colors = [cmap(i) for i in range(n)]
        return [*colors, to_rgba(MISSING_COLOR)] if missing else colors


    def scatter(
        ax: "Axes", x: npt.NDArray[np.float64], y: npt.NDArray[np.float64], *, by: Groups | None, title: str | None
    ) -> None:
        """One scatter, or one per group with a legend titled ``title``."""
        if by is None:
            ax.scatter(x, y)
            return
        codes, labels, colors = by
        for code, (label, color) in enumerate(zip(labels, colors, strict=True)):
            keep = codes == code
            ax.scatter(x[keep], y[keep], color=color, label=label)
        ax.legend(title=title)


    def label_ticks(ax: "Axes", names: list[str], *, axis: Literal["x", "y"]) -> None:
        """One tick per name, or none past ``MAX_LABELS`` names."""
        shown = names if len(names) <= MAX_LABELS else []
        if axis == "x":
            ax.set_xticks(np.arange(len(shown)), labels=shown, rotation=90)
        else:
            ax.set_yticks(np.arange(len(shown)), labels=shown)
    ```
  - `src/biotapy/pl/_abundance.py`:
    ```python
    """Abundance plots: stacked bars and a heatmap of the table."""

    from typing import TYPE_CHECKING, cast

    import numpy as np
    import numpy.typing as npt
    import pandas as pd
    import scipy.sparse as sp
    from anndata import AnnData

    from biotapy._core import sum_by

    from ._common import RGBA, groups, label_ticks, new_axes, obs_groups, table

    if TYPE_CHECKING:
        from matplotlib.axes import Axes

    # phyloseq's plot_heatmap colours: low "#000033", high "#66CCFF", and na.value "black", which zeros become on its log scale.
    _LOW, _HIGH, _ZERO = "#000033", "#66CCFF", "black"


    def bar(
        adata: AnnData, fill: str, *, x: str | None = None, layer: str | None = None, ax: "Axes | None" = None
    ) -> "Axes":
        """Stacked bars of abundance, one segment per group of ``fill``.

        Parameters
        ----------
        adata
            Samples x features.
        fill
            A ``var`` column, such as a taxonomic rank, or an ``obs`` column; each of its
            values is one coloured segment.
        x
            An ``obs`` column: one bar per value, summing its samples. By default one bar per sample.
        layer
            Plot ``layers[layer]``, such as ``"relative"`` from :func:`biotapy.pp.relative`, instead of ``X``.
        ax
            Axes to draw on; by default a new figure's.

        Returns
        -------
        matplotlib.axes.Axes
            Bars in ``obs`` order (samples) or group order (``x``), segments in group order,
            and a legend titled ``fill``.

        Raises
        ------
        KeyError
            ``fill`` is in neither or both of ``var`` and ``obs``, ``x`` is not an ``obs``
            column, or ``layers[layer]`` is missing.
        TypeError
            ``fill`` or ``x`` is a numeric column.

        Notes
        -----
        R equivalent: ``phyloseq::plot_bar``
        Guide: :doc:`/guide/plotting`

        phyloseq stacks one outlined segment per feature and sample, which is 500,000
        rectangles on GlobalPatterns. biotapy sums the features of each ``fill`` group
        first: bar and segment heights are the same, without the per-feature outlines.
        Groups are categories in their order, else sorted values; missing values form
        the last group, ``NA``, so bars still sum every feature. ``facet_grid`` has no
        equivalent: draw one subset per axes of ``matplotlib.pyplot.subplots``.

        Examples
        --------
        >>> import biotapy as bt
        >>> ax = bt.pl.bar(bt.datasets.toy(), "phylum")
        >>> len(ax.patches)  # 6 samples x 3 phyla
        18
        """
        values = table(adata, layer)
        segments, labels, colors = _segments(adata, values, fill)
        heights, names = _by_x(adata, segments, x)
        ax = new_axes(ax)
        positions, bottom = np.arange(len(names)), np.zeros(len(names))
        for column, (label, color) in enumerate(zip(labels, colors, strict=True)):
            ax.bar(positions, heights[:, column], bottom=bottom, color=color, label=label)
            bottom += heights[:, column]
        label_ticks(ax, names, axis="x")
        ax.set_xlabel(x or "sample")
        ax.set_ylabel(layer or "abundance")
        ax.legend(title=fill)
        return ax


    def heatmap(adata: AnnData, *, layer: str | None = None, ax: "Axes | None" = None) -> "Axes":
        """Heatmap of the table, samples as columns and features as rows, on a log colour scale.

        Parameters
        ----------
        adata
            Samples x features, drawn in ``obs`` and ``var`` order.
        layer
            Plot ``layers[layer]``, such as ``"relative"`` from :func:`biotapy.pp.relative`, instead of ``X``.
        ax
            Axes to draw on; by default a new figure's.

        Returns
        -------
        matplotlib.axes.Axes
            The image, features x samples, with a colour bar. Zeros are drawn black.

        Raises
        ------
        KeyError
            ``layers[layer]`` is missing.
        ValueError
            The table holds no positive value, so the log scale has nothing to show.

        Notes
        -----
        R equivalent: ``phyloseq::plot_heatmap``
        Guide: :doc:`/guide/plotting`

        phyloseq orders samples and features by the angle of their scores on a
        two-axis NMDS that it computes on the fly; biotapy keeps ``obs`` and ``var``
        order and computes nothing, so sort first (see the guide). Colours are
        phyloseq's: ``#000033`` to ``#66CCFF`` on a log scale, zeros black. Tick labels
        are drawn up to 250 names per axis, as phyloseq's ``max.label``. matplotlib needs
        a dense array: the table is densified once (8 bytes x samples x features).

        Examples
        --------
        >>> import biotapy as bt
        >>> ax = bt.pl.heatmap(bt.datasets.toy())
        >>> ax.images[0].get_array().shape  # features x samples
        (8, 6)
        """
        from matplotlib.colors import LinearSegmentedColormap, LogNorm

        # matplotlib's imshow needs a dense array (rules.md R6.2): one dense copy, features x samples.
        dense = table(adata, layer).T.toarray()
        if not (dense > 0).any():
            msg = "the table holds no positive value to draw on a log scale"
            raise ValueError(msg)
        ax = new_axes(ax)
        cmap = LinearSegmentedColormap.from_list("phyloseq", [_LOW, _HIGH]).with_extremes(bad=_ZERO)
        image = ax.imshow(dense, cmap=cmap, norm=LogNorm(), aspect="auto", interpolation="nearest")
        ax.figure.colorbar(image, ax=ax, label=layer or "abundance")
        label_ticks(ax, adata.obs_names.tolist(), axis="x")
        label_ticks(ax, adata.var_names.tolist(), axis="y")
        ax.set_xlabel("sample")
        ax.set_ylabel("feature")
        return ax


    def _segments(
        adata: AnnData, values: sp.csr_matrix, fill: str
    ) -> tuple[npt.NDArray[np.float64], list[str], list[RGBA]]:
        in_var, in_obs = fill in adata.var.columns, fill in adata.obs.columns
        if in_var == in_obs:
            msg = f"fill={fill!r} must be a column of exactly one of var and obs"
            raise KeyError(msg)
        if in_var:
            codes, labels, colors = groups(cast("pd.Series", adata.var[fill]), arg="fill")
            # Samples x groups: small, unlike X, so it is safe to densify.
            return np.asarray(sum_by(values, codes, len(labels)).toarray(), dtype=np.float64), labels, colors
        codes, labels, colors = obs_groups(adata, fill, arg="fill")
        segments = np.zeros((adata.n_obs, len(labels)))
        segments[np.arange(adata.n_obs), codes] = np.asarray(values.sum(axis=1)).ravel()
        return segments, labels, colors


    def _by_x(
        adata: AnnData, segments: npt.NDArray[np.float64], x: str | None
    ) -> tuple[npt.NDArray[np.float64], list[str]]:
        if x is None:
            return segments, adata.obs_names.tolist()
        codes, names, _ = obs_groups(adata, x, arg="x")
        indicator = np.zeros((len(names), adata.n_obs))
        indicator[codes, np.arange(adata.n_obs)] = 1
        return indicator @ segments, names
    ```
  - `src/biotapy/pl/_richness.py`:
    ```python
    """Alpha diversity per sample, as stored by ``tl.alpha``."""

    from typing import TYPE_CHECKING

    import numpy as np
    from anndata import AnnData

    from ._common import label_ticks, new_axes, obs_groups, scatter

    if TYPE_CHECKING:
        from matplotlib.axes import Axes


    def richness(
        adata: AnnData, metric: str, *, x: str | None = None, color: str | None = None, ax: "Axes | None" = None
    ) -> "Axes":
        """One point per sample for a stored alpha diversity metric.

        Parameters
        ----------
        adata
            Samples x features with ``obs['alpha_<metric>']``, written by
            :func:`biotapy.tl.alpha` with ``inplace=True``.
        metric
            The metric to plot, such as ``"shannon"``.
        x
            An ``obs`` column whose values place the points. By default one position per sample.
        color
            An ``obs`` column that colours the points, with a legend.
        ax
            Axes to draw on; by default a new figure's.

        Returns
        -------
        matplotlib.axes.Axes
            The points, one per sample with a finite value.

        Raises
        ------
        KeyError
            ``obs['alpha_<metric>']`` is missing (the message names the call that writes
            it), or ``x`` or ``color`` is not an ``obs`` column.
        TypeError
            ``x`` or ``color`` is a numeric column.

        Notes
        -----
        R equivalent: ``phyloseq::plot_richness``
        Guide: :doc:`/guide/plotting`

        phyloseq computes ``estimate_richness`` and draws every measure in its own
        facet; biotapy reads what :func:`biotapy.tl.alpha` stored and draws one metric
        per axes. Samples with a NaN value, such as Shannon of an all-zero sample, are
        left out, as ``geom_point(na.rm = TRUE)`` does. There are no standard-error bars:
        ``tl.alpha`` stores none.

        Examples
        --------
        >>> import biotapy as bt
        >>> tdata = bt.datasets.toy()
        >>> bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)
        >>> ax = bt.pl.richness(tdata, "shannon", x="group")
        >>> ax.collections[0].get_offsets().shape
        (6, 2)
        """
        column = f"alpha_{metric}"
        if column not in adata.obs.columns:
            msg = f"metric={metric!r}: no obs[{column!r}]; run bt.tl.alpha(adata, metrics=[{metric!r}], inplace=True) first"
            raise KeyError(msg)
        values = np.asarray(adata.obs[column], dtype=np.float64)
        by = obs_groups(adata, color, arg="color") if color is not None else None
        if x is None:
            positions, names = np.arange(adata.n_obs, dtype=np.float64), adata.obs_names.tolist()
        else:
            codes, names, _ = obs_groups(adata, x, arg="x")
            positions = codes.astype(np.float64)
        finite = np.isfinite(values)
        if by is not None:
            by = (by[0][finite], by[1], by[2])
        ax = new_axes(ax)
        scatter(ax, positions[finite], values[finite], by=by, title=color)
        label_ticks(ax, names, axis="x")
        ax.set_xlabel(x or "sample")
        ax.set_ylabel(metric)
        return ax
    ```
  - `src/biotapy/pl/_ordination.py`:
    ```python
    """Ordination plots: samples on two axes, and the scree of a PCoA."""

    from typing import TYPE_CHECKING, Literal, get_args

    import numpy as np
    import numpy.typing as npt
    from anndata import AnnData

    from ._common import new_axes, obs_groups, scatter

    if TYPE_CHECKING:
        from matplotlib.axes import Axes

    Basis = Literal["pcoa", "nmds"]
    # The call that writes obsm["X_<basis>"] and uns["biotapy"][<basis>], for the error raised when they are missing.
    _WRITTEN_BY = {"pcoa": "bt.tl.pcoa(adata, inplace=True)", "nmds": "bt.tl.nmds(adata, inplace=True)"}


    def ordination(
        adata: AnnData,
        *,
        basis: Basis = "pcoa",
        components: tuple[int, int] = (1, 2),
        color: str | None = None,
        ax: "Axes | None" = None,
    ) -> "Axes":
        """Samples on two axes of a stored ordination.

        Parameters
        ----------
        adata
            Samples x features with ``obsm['X_pcoa']`` or ``obsm['X_nmds']``, written by
            :func:`biotapy.tl.pcoa` or :func:`biotapy.tl.nmds` with ``inplace=True``.
        basis
            ``"pcoa"`` or ``"nmds"``.
        components
            The two axes to draw, counted from 1.
        color
            An ``obs`` column that colours the samples, with a legend.
        ax
            Axes to draw on; by default a new figure's.

        Returns
        -------
        matplotlib.axes.Axes
            One point per sample. PCoA axis labels carry the percentage each axis
            explains; an NMDS plot notes the stress in its lower right corner.

        Raises
        ------
        ValueError
            ``basis`` is not ``"pcoa"`` or ``"nmds"``, or ``components`` does not name two stored axes.
        KeyError
            The ordination is missing (the message names the call that writes it), or
            ``color`` is not an ``obs`` column.
        TypeError
            ``color`` is a numeric column.

        Notes
        -----
        R equivalent: ``phyloseq::plot_ordination``
        Guide: :doc:`/guide/plotting`

        As in phyloseq, PCoA labels append ``round(100 * Relative_eig, 1)`` and NMDS
        labels carry no percentage. The stress note is biotapy's; phyloseq shows none.
        Taxa, biplot and split plots are not in 0.1.

        Examples
        --------
        >>> import biotapy as bt
        >>> tdata = bt.datasets.toy()
        >>> bt.tl.beta(tdata, inplace=True)
        >>> bt.tl.pcoa(tdata, inplace=True)
        >>> bt.pl.ordination(tdata, color="group").get_xlabel()
        'PC1 [93.8%]'
        """
        coords, summary = _stored(adata, basis)
        first, second = _components(components, coords.shape[1])
        by = obs_groups(adata, color, arg="color") if color is not None else None
        ax = new_axes(ax)
        scatter(ax, coords[:, first], coords[:, second], by=by, title=color)
        if basis == "pcoa":
            shares = np.asarray(summary["proportion_explained"], dtype=np.float64)
            ax.set_xlabel(f"PC{first + 1} [{100 * shares[first]:.1f}%]")
            ax.set_ylabel(f"PC{second + 1} [{100 * shares[second]:.1f}%]")
        else:
            ax.set_xlabel(f"NMDS{first + 1}")
            ax.set_ylabel(f"NMDS{second + 1}")
            stress = float(np.asarray(summary["stress"]))
            ax.text(0.99, 0.01, f"stress {stress:.3f}", transform=ax.transAxes, ha="right", va="bottom")
        return ax


    def scree(adata: AnnData, *, ax: "Axes | None" = None) -> "Axes":
        """Bars of the proportion of variation each stored PCoA axis explains.

        Parameters
        ----------
        adata
            Samples x features with ``uns['biotapy']['pcoa']``, written by
            :func:`biotapy.tl.pcoa` with ``inplace=True``.
        ax
            Axes to draw on; by default a new figure's.

        Returns
        -------
        matplotlib.axes.Axes
            One bar per stored axis: ``PC1``, ``PC2``, ...

        Raises
        ------
        KeyError
            The PCoA is missing; the message names the call that writes it.

        Notes
        -----
        R equivalent: ``phyloseq::plot_scree``
        Guide: :doc:`/guide/plotting`

        The bars are ``proportion_explained``, each eigenvalue over their sum, as
        phyloseq's ``x / sum(x)``. phyloseq draws every axis; biotapy draws the axes
        ``tl.pcoa`` stored, 10 by default, so pass ``n_components`` for more.

        Examples
        --------
        >>> import biotapy as bt
        >>> tdata = bt.datasets.toy()
        >>> bt.tl.beta(tdata, inplace=True)
        >>> bt.tl.pcoa(tdata, inplace=True)
        >>> len(bt.pl.scree(tdata).patches)  # at most n_obs - 1 axes
        5
        """
        _, summary = _stored(adata, "pcoa")
        shares = np.asarray(summary["proportion_explained"], dtype=np.float64)
        ax = new_axes(ax)
        positions = np.arange(shares.size)
        ax.bar(positions, shares)
        ax.set_xticks(positions, labels=[f"PC{i}" for i in range(1, shares.size + 1)])
        ax.set_xlabel("axis")
        ax.set_ylabel("proportion explained")
        return ax


    def _stored(adata: AnnData, basis: str) -> tuple[npt.NDArray[np.float64], dict[str, object]]:
        if basis not in get_args(Basis):
            msg = f"basis must be one of {list(get_args(Basis))}, got {basis!r}"
            raise ValueError(msg)
        key = f"X_{basis}"
        summary = adata.uns.get("biotapy", {}).get(basis)
        if key not in adata.obsm or summary is None:
            msg = f"basis={basis!r}: no obsm[{key!r}] or uns['biotapy'][{basis!r}]; run {_WRITTEN_BY[basis]} first"
            raise KeyError(msg)
        return np.asarray(adata.obsm[key], dtype=np.float64), dict(summary)


    def _components(components: tuple[int, int], n_axes: int) -> tuple[int, int]:
        if len(components) != 2 or not all(1 <= axis <= n_axes for axis in components):
            msg = f"components must be two axis numbers from 1 to {n_axes}, got {components!r}"
            raise ValueError(msg)
        return components[0] - 1, components[1] - 1
    ```
  - `src/biotapy/pl/__init__.py`, imports only (R4.1):
    ```python
    """Plots of what tl and pp stored; pl reads slots and computes nothing (contracts/module-boundaries)."""

    from ._abundance import bar, heatmap
    from ._ordination import ordination, scree
    from ._richness import richness

    __all__ = ["bar", "heatmap", "ordination", "richness", "scree"]
    ```
  - `src/biotapy/__init__.py`: `from . import datasets, io, pl, pp, tl` and
    `__all__ = ["__version__", "datasets", "io", "pl", "pp", "tl"]`.
  - import-linter already declares the `(pl) | (ml) | (da)` top layer, so
    `pyproject.toml` needs no other change.
  - Why the code looks as it does:
    - `cast("pd.Series", adata.obs[column])`: anndata types a column as
      `Series | DataArray`, and mypy rejected passing that union to `groups`.
      The same cast pattern is in `pp/_glom.py`.
    - `scatter(ax, x, y, *, by, title)`: ruff's `PLR0917` allows at most 3
      positional parameters, and a fourth failed the prototype.
    - `_segments` densifies only the samples x groups sums, never `X`.
      `heatmap` densifies the table once, which R6.2 allows because `imshow`
      needs a dense array; its `Notes` state the cost.
    - In `_by_x`, a missing `x` value is coded `len(labels) - 1`, the `NA` row.
- [x] **Step 5: Run.**
  - Run `uv run --group test pytest tests/pl src/biotapy/pl -q`: 37 tests and
    5 doctests pass (42 in the prototype, 5.4 s).
  - Then run
    `uv run --group test --with pytest-xdist --with pytest-randomly pytest tests/pl -q -n 4`,
    because CI's hatch-test environment runs xdist with pytest-randomly. It
    must pass in random order.
- [x] **Step 6: Docs.**
  - `docs/conf.py`: in `intersphinx_mapping`, after the networkx entry, add
    `"matplotlib": ("https://matplotlib.org/stable/", None),`. Without it the
    `nitpicky` build cannot resolve `matplotlib.axes.Axes` in the signatures.
  - `docs/api.md`: add a section after Tools, the order scanpy uses
    (`pp`, `tl`, `pl`):
    ````markdown

    ## Plots

    ```{eval-rst}
    .. module:: biotapy.pl
    .. currentmodule:: biotapy

    .. autosummary::
        :toctree: generated

        pl.bar
        pl.heatmap
        pl.ordination
        pl.richness
        pl.scree
    ```
    ````
  - `docs/guide/index.md`: add `plotting` to the toctree after `ordination`.
  - Create `docs/guide/plotting.md`:
    ````markdown
    # Plotting

    `bt.pl` draws what `tl` and `pp` stored; it computes no diversity and no ordination. Each
    function takes `ax=` to draw on existing axes, or makes a new figure, and returns the
    `matplotlib.axes.Axes`, so you finish the plot with matplotlib. When a stored result is
    missing, the error names the call to run:

    | Function | Reads | Written by |
    |---|---|---|
    | `pl.bar` | `X` or `layers[layer]`, and `var`/`obs` columns | `pp.relative` for `layer="relative"` |
    | `pl.heatmap` | `X` or `layers[layer]` | `pp.relative` for `layer="relative"` |
    | `pl.richness` | `obs["alpha_<metric>"]` | `tl.alpha(..., inplace=True)` |
    | `pl.ordination` | `obsm["X_pcoa"]` or `obsm["X_nmds"]`, and `uns["biotapy"]["pcoa"]` or `["nmds"]` | `tl.pcoa` or `tl.nmds` with `inplace=True` |
    | `pl.scree` | `uns["biotapy"]["pcoa"]["proportion_explained"]` | `tl.pcoa(..., inplace=True)` |

    ```python
    import biotapy as bt

    tdata = bt.datasets.toy()
    bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)
    bt.tl.beta(tdata, inplace=True)
    bt.tl.pcoa(tdata, inplace=True)

    bt.pl.bar(tdata, "phylum")  # one bar per sample, one segment per phylum
    bt.pl.richness(tdata, "shannon", x="group")
    ax = bt.pl.ordination(tdata, color="group")  # axis labels carry the % explained
    ax.set_title("Bray-Curtis PCoA")
    ```

    ## Groups and colours

    `fill`, `x` and `color` name a column. Its categories keep their order; other values are
    sorted. Missing values form a last group, `NA`, drawn grey, as in ggplot2. A numeric column
    raises a `TypeError`: convert it with `.astype("category")` for one group per value, or bin it
    with `pandas.cut`.

    ## Bars: `pl.bar`

    phyloseq's `plot_bar` stacks one outlined rectangle per feature and sample, about 500,000 on
    GlobalPatterns. `pl.bar` sums the features of each `fill` group first, so the bars and
    segments have the same heights without the per-feature outlines, and a phylum plot of
    GlobalPatterns draws in under a second. `fill` can be a `var` column (a rank) or an `obs`
    column; `x` groups samples into one bar per value. `layer="relative"` plots proportions from
    `bt.pp.relative`.

    A `fill` with hundreds of groups is slow and its legend unreadable; aggregate first with
    `bt.pp.tax_glom` or keep the most abundant features.

    ## Facets

    `facet_grid` becomes one subset per axes:

    ```python
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, sharey=True, figsize=(12, 4))
    for ax, phylum in zip(axes, ["Bacteroidota", "Firmicutes", "Proteobacteria"]):
        bt.pl.bar(tdata[:, tdata.var["phylum"] == phylum], "group", ax=ax)
        ax.set_title(phylum)
    ```

    ## Heatmaps: `pl.heatmap`

    `pl.heatmap` draws the table in `obs` and `var` order, features as rows, on phyloseq's log
    colour scale with zeros black. phyloseq instead orders samples and taxa by their angle on a
    two-axis NMDS it computes on the fly. To do the same, compute the NMDS and sort first:

    ```python
    import numpy as np

    bt.tl.nmds(tdata, seed=0, inplace=True)
    x, y = tdata.obsm["X_nmds"].T
    ax = bt.pl.heatmap(tdata[np.argsort(np.arctan2(y, x))])
    ```

    Tick labels name samples and features, up to 250 per axis. Relabel them with matplotlib,
    for example by sample type:
    `ax.set_xticks(range(tdata.n_obs), labels=tdata.obs["group"], rotation=90)`.
    ````
- [x] **Step 7: Knowledge and rules.**
  - `decisions/optional-heavy-dependencies.md` (ruling 22):
    - In the first Decision bullet, "scikit-bio (`>=0.7.4,<0.8`), matplotlib,
      pooch (Phase 0-1); scikit-learn (`>=1.8`, Phase 1, NMDS); mudata
      (Phase 2)." becomes "scikit-bio (`>=0.7.4,<0.8`), pooch (Phase 0-1);
      scikit-learn (`>=1.8`, Phase 1, NMDS); matplotlib (`>=3.8`, Phase 1,
      `pl`); mudata (Phase 2)."
    - At the end of that bullet, after the Task 1.17 sentence, add: "Task 1.18
      added matplotlib (`>=3.8`, approved 2026-09-27) for `pl`. 3.8 is its
      first release with type information and CPython 3.12 wheels. `pl`
      imports it inside its functions, and pyplot only when it must make a
      figure, so `import biotapy` does not load it. It brings contourpy,
      cycler, fonttools, kiwisolver, pillow and pyparsing."
    - Consequences: "CI runs one job with all extras and one job with none, so
      a lazy import leaking to module level fails fast." becomes "CI imports
      every module with no extras installed (`import-without-extras`), so a
      lazy import leaking to module level fails fast. A job with all extras
      comes with the first extra."
  - `contracts/r-golden-parity.md`: add Statement 7: "`pl` functions have an
    R equivalent but no golden test. They draw numbers that `tl` stores, and
    `tl`'s golden tests check those numbers (controller ruling 2026-09-27;
    rules.md R11.2)."
  - `playbooks/add-a-function.md`, Step 4: the bullet "golden test when an R
    equivalent exists, per r-golden-parity." becomes "golden test when an R
    equivalent exists, per r-golden-parity; not in `pl`, whose plots draw
    `tl` results that have one."
  - `rules.md` R11.2: "a golden test when an R equivalent exists
    ([r-golden-parity](...))" becomes "a golden test when an R equivalent
    exists, except in `pl`, whose plots draw `tl` results that have one
    ([r-golden-parity](...))". Make this edit only because approving this
    plan approves it (ruling 21).
  - Add a log line. Tick 1.18 here.
- [x] **Step 8: Run, gate and commit.**
  - Tests: `uv run --group test pytest tests/pl tests/tl src/biotapy/pl -q`.
  - Import time:
    - run `uv run python -X importtime -c "import biotapy" 2>&1 | tail -1`
      three times and record the times next to Step 1's. The prototype
      measured 1.19-1.40 s: the same range, since the machine's load moves
      it by about 0.1 s;
    - run `uv run python -X importtime -c "import biotapy" 2>&1 | grep -E "biotapy\.pl$"`.
      `pl`'s own cumulative time was 0.47-0.54 ms;
    - `uv run python -X importtime -c "import biotapy" 2>&1 | grep -c matplotlib`
      must print `0`.
  - The three gates: `uvx prek run --all-files`, `uv run --group test pytest`,
    and `uv run --group doc sphinx-build -W -b html docs docs/_build/html`.
  - Commit `feat(pl): add bar, heatmap, richness, ordination and scree plots`.
    Stage `pyproject.toml`, the five `src/biotapy/pl/` files,
    `src/biotapy/__init__.py`, `tests/conftest.py`, the five `tests/pl/`
    files, the three docs files and `docs/guide/plotting.md`, `rules.md`, the
    three knowledge concepts, the roadmap and the log.

### Task 1.19: Docstring contract test and Coming-from-R table

**Files:**
- Create:
  - `docs/extensions/coming_from_r.py`, `docs/_data/r_idioms.toml`, `docs/coming_from_r.md`
  - `tests/test_docstrings.py`, `tests/test_coming_from_r.py`
- Modify:
  - `docs/conf.py` (`exclude_patterns`), `docs/index.md` (toctree)
  - `.knowledge/contracts/function-shape.md`
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:** every name in every public subpackage's `__all__`: `datasets`,
  `io`, `pl` (1.18), `pp` and `tl`. `_core` is not in `biotapy.__all__`, so it
  is skipped without a special case.
- **Produces:**
  - In `docs/extensions/coming_from_r.py`, imported by the tests through
    `importlib.util.spec_from_file_location` (they never import Sphinx; the
    module imports it only under `TYPE_CHECKING`):
    - `r_equivalents(doc: str) -> list[str]`, which raises `ValueError`
      unless there is exactly one parseable `R equivalent:` line;
    - `public_functions() -> list[tuple[str, Callable[..., object]]]`;
    - `rows() -> dict[str, list[str]]`;
    - `render() -> str`;
    - `write_table(app) -> None`, the `builder-inited` hook;
    - `TABLE = Path("generated") / "coming_from_r_table.md"` and `IDIOMS`.
  - The page `docs/coming_from_r.md`, published at
    `https://biotapy.readthedocs.io/page/coming_from_r.html`, which 1.23's
    README links.

- [x] **Step 1: Failing tests.**
  - `tests/test_docstrings.py`:
    ```python
    """Every public function's docstring keeps the contracts/function-shape skeleton (rules.md R8.2)."""

    import importlib.util
    import re
    from pathlib import Path

    import pytest

    SPEC = importlib.util.spec_from_file_location(
        "coming_from_r", Path(__file__).resolve().parents[1] / "docs" / "extensions" / "coming_from_r.py"
    )
    coming_from_r = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(coming_from_r)
    PUBLIC = coming_from_r.public_functions()


    def test_every_public_subpackage_is_covered():
        assert {name.split(".")[1] for name, _ in PUBLIC} == {"datasets", "io", "pl", "pp", "tl"}


    @pytest.mark.parametrize(("name", "function"), PUBLIC, ids=[name for name, _ in PUBLIC])
    def test_docstring_has_the_contract_sections(name, function):
        doc = function.__doc__ or ""
        coming_from_r.r_equivalents(doc)  # raises unless there is exactly one parseable line
        assert re.search(r"^\s*Notes\n\s*-----\n\s*R equivalent: .+\n\s*Guide: :doc:`/guide/\w+`$", doc, re.MULTILINE), name
        assert re.search(r"^\s*Examples\n\s*--------\n\s*>>> ", doc, re.MULTILINE), name


    @pytest.mark.parametrize(
        "doc",
        [
            "Notes\n-----\nGuide: :doc:`/guide/x`\n",
            "R equivalent: ``phyloseq::a``\nR equivalent: ``phyloseq::b``\n",
            "R equivalent: phyloseq::a\n",
            "R equivalent: ``phyloseq::a``; ``vegan::b``\n",
        ],
        ids=["missing", "twice", "no backticks", "semicolon"],
    )
    def test_malformed_r_lines_raise(doc):
        with pytest.raises(ValueError, match="exactly one line"):
            coming_from_r.r_equivalents(doc)


    def test_r_line_items_and_none():
        assert coming_from_r.r_equivalents("    R equivalent: ``phyloseq::ordinate``, ``ape::pcoa``\n") == [
            "phyloseq::ordinate",
            "ape::pcoa",
        ]
        assert coming_from_r.r_equivalents("R equivalent: none\n") == []
    ```
  - `tests/test_coming_from_r.py`:
    ```python
    """The generated Coming-from-R table covers phyloseq's 31 core functions (Phase 1 exit gate)."""

    import importlib.util
    from pathlib import Path
    from types import SimpleNamespace

    import pytest

    ROOT = Path(__file__).resolve().parents[1]
    SPEC = importlib.util.spec_from_file_location("coming_from_r", ROOT / "docs" / "extensions" / "coming_from_r.py")
    coming_from_r = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(coming_from_r)

    # .knowledge/roadmap/phase-1-core.md, "The 31 phyloseq functions the table must cover".
    PHYLOSEQ_31 = """otu_table sample_data tax_table phy_tree refseq nsamples ntaxa sample_names taxa_names
    sample_sums taxa_sums rank_names sample_variables get_taxa_unique prune_taxa prune_samples subset_taxa
    subset_samples filter_taxa transform_sample_counts rarefy_even_depth tax_glom estimate_richness distance
    UniFrac ordinate plot_bar plot_richness plot_ordination plot_heatmap import_biom""".split()
    NOT_IN_0_1 = ["tip_glom", "merge_samples", "psmelt", "plot_tree", "plot_net"]


    def test_the_list_has_31_functions():
        assert len(set(PHYLOSEQ_31)) == 31


    @pytest.mark.parametrize("name", PHYLOSEQ_31)
    def test_table_maps_each_of_the_31(name):
        cells = coming_from_r.rows()[f"phyloseq::{name}"]
        assert cells and "not in 0.1" not in cells


    @pytest.mark.parametrize("name", NOT_IN_0_1)
    def test_uncovered_functions_are_marked(name):
        assert coming_from_r.rows()[f"phyloseq::{name}"] == ["not in 0.1"]


    def test_plot_functions_link_to_pl():
        assert coming_from_r.rows()["phyloseq::plot_bar"] == ["{func}`bt.pl.bar <biotapy.pl.bar>`"]


    def test_render_is_one_sorted_markdown_table():
        lines = coming_from_r.render().splitlines()
        assert lines[:2] == ["| R | biotapy |", "|---|---|"]
        names = [line.split("`")[1] for line in lines[2:]]
        assert names == sorted(names, key=str.lower) and len(names) == len(coming_from_r.rows())


    def test_an_idiom_a_docstring_already_maps_raises(tmp_path, monkeypatch):
        idioms = tmp_path / "r_idioms.toml"
        idioms.write_text('[idioms]\n"phyloseq::tax_glom" = "`x`"\n', encoding="utf-8")
        monkeypatch.setattr(coming_from_r, "IDIOMS", idioms)
        with pytest.raises(ValueError, match=r"\['phyloseq::tax_glom'\] are mapped by a docstring"):
            coming_from_r.rows()


    def test_write_table_writes_the_file_the_page_includes(tmp_path):
        coming_from_r.write_table(SimpleNamespace(srcdir=str(tmp_path)))
        assert (tmp_path / coming_from_r.TABLE).read_text(encoding="utf-8") == coming_from_r.render()
        page = (ROOT / "docs" / "coming_from_r.md").read_text(encoding="utf-8")
        assert f"```{{include}} {coming_from_r.TABLE.as_posix()}\n```" in page
    ```
  - What the tests pin:
    - The contract test uses the generator's own parser, so "parseable"
      means exactly what the table needs.
    - The `Examples` check asks only for a `>>> ` line. The `io` readers'
      examples start with `import tempfile`, and the contract allows that.
- [x] **Step 2: Run, expect failure.**
  `uv run --group test pytest tests/test_docstrings.py tests/test_coming_from_r.py -q`
  errors while collecting both files with
  `FileNotFoundError: [Errno 2] No such file or directory: '.../docs/extensions/coming_from_r.py'`.
- [x] **Step 3: Implement.**
  - `docs/extensions/coming_from_r.py`. `conf.py` already puts
    `docs/extensions` on `sys.path` and loads every `.py` there as an
    extension, so the file needs no registration:
    ```python
    """Write the Coming-from-R table from the public docstrings' ``R equivalent:`` lines (rules.md R8.4)."""

    import importlib
    import re
    import tomllib
    from collections.abc import Callable
    from pathlib import Path
    from typing import TYPE_CHECKING

    if TYPE_CHECKING:
        from sphinx.application import Sphinx

    # contracts/function-shape: exactly one such line, holding ``pkg::fn`` items separated by ", ", or "none".
    R_LINE = re.compile(r"^\s*R equivalent:(?P<rest>.*)$", re.MULTILINE)
    R_ITEMS = re.compile(r"none|``[\w.]+::[\w.]+``(?:, ``[\w.]+::[\w.]+``)*")
    R_ITEM = re.compile(r"``([\w.]+::[\w.]+)``")
    IDIOMS = Path(__file__).resolve().parents[1] / "_data" / "r_idioms.toml"
    # Relative to the docs source directory; docs/generated/ is git-ignored and left out of the sdist.
    TABLE = Path("generated") / "coming_from_r_table.md"


    def public_functions() -> list[tuple[str, Callable[..., object]]]:
        """``("bt.<module>.<name>", function)`` for every name in every public subpackage's ``__all__``."""
        biotapy = importlib.import_module("biotapy")
        modules = [importlib.import_module(f"biotapy.{name}") for name in biotapy.__all__ if name != "__version__"]
        return [
            (f"bt.{module.__name__.split('.')[-1]}.{name}", getattr(module, name))
            for module in modules
            for name in module.__all__
        ]


    def r_equivalents(doc: str) -> list[str]:
        """The ``pkg::fn`` items of the docstring's one ``R equivalent:`` line; ``[]`` when it says ``none``."""
        lines = R_LINE.findall(doc)
        if len(lines) != 1 or R_ITEMS.fullmatch(lines[0].strip()) is None:
            msg = f"need exactly one line 'R equivalent: ``pkg::fn``, ...' or 'R equivalent: none', found {lines}"
            raise ValueError(msg)
        return R_ITEM.findall(lines[0])


    def rows() -> dict[str, list[str]]:
        """R call -> the biotapy cells for it: the docstrings' functions, then the idioms file."""
        table: dict[str, list[str]] = {}
        for name, function in public_functions():
            for r_name in r_equivalents(function.__doc__ or ""):
                table.setdefault(r_name, []).append(f"{{func}}`{name} <biotapy.{name.removeprefix('bt.')}>`")
        idioms: dict[str, str] = tomllib.loads(IDIOMS.read_text(encoding="utf-8"))["idioms"]
        both = sorted(set(idioms) & set(table))
        if both:
            msg = f"{both} are mapped by a docstring and by {IDIOMS.name}; keep only the docstring"
            raise ValueError(msg)
        return table | {r_name: [python] for r_name, python in idioms.items()}


    def render() -> str:
        """The table as MyST markdown, one row per R call, sorted by package and name."""
        table = rows()
        lines = ["| R | biotapy |", "|---|---|"]
        lines += [f"| `{r_name}` | {', '.join(table[r_name])} |" for r_name in sorted(table, key=str.lower)]
        return "\n".join(lines) + "\n"


    def write_table(app: "Sphinx") -> None:
        """``builder-inited`` hook: write the table before Sphinx reads any source."""
        path = Path(app.srcdir) / TABLE
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render(), encoding="utf-8")


    def setup(app: "Sphinx") -> dict[str, bool]:
        """Register the hook."""
        app.connect("builder-inited", write_table)
        return {"parallel_read_safe": True, "parallel_write_safe": True}
    ```
  - `docs/_data/r_idioms.toml`: 19 idioms, the 18 of the 31 that no
    docstring maps plus `get_variable`, which the vignette uses; and the 5
    "not in 0.1" rows:
    ```toml
    # R calls with no biotapy function of their own, for the generated Coming-from-R table
    # (docs/extensions/coming_from_r.py). Every other row comes from a public docstring's
    # "R equivalent:" line; a call listed there must not be listed here too.
    # Values are MyST markdown, inserted as they are.

    [idioms]
    "phyloseq::otu_table" = '`tdata.X` (samples x features)'
    "phyloseq::sample_data" = '`tdata.obs`'
    "phyloseq::tax_table" = '`tdata.var`, one column per rank'
    "phyloseq::phy_tree" = '`tdata.vart["phylo"]`'
    "phyloseq::refseq" = '`tdata.var["sequence"]`'
    "phyloseq::nsamples" = '`tdata.n_obs`'
    "phyloseq::ntaxa" = '`tdata.n_vars`'
    "phyloseq::sample_names" = '`tdata.obs_names`'
    "phyloseq::taxa_names" = '`tdata.var_names`'
    "phyloseq::sample_sums" = '`tdata.X.sum(axis=1)`'
    "phyloseq::taxa_sums" = '`tdata.X.sum(axis=0)`'
    "phyloseq::rank_names" = '`tdata.var.columns`'
    "phyloseq::sample_variables" = '`tdata.obs.columns`'
    "phyloseq::get_variable" = '`tdata.obs[column]`'
    "phyloseq::get_taxa_unique" = '`tdata.var[rank].unique()`'
    "phyloseq::prune_taxa" = '`tdata[:, keep].copy()`'
    "phyloseq::prune_samples" = '`tdata[keep].copy()`'
    "phyloseq::subset_taxa" = '`tdata[:, tdata.var["phylum"] == "Chlamydiae"].copy()`'
    "phyloseq::subset_samples" = '`tdata[tdata.obs["SampleType"] == "Feces"].copy()`'
    "phyloseq::tip_glom" = "not in 0.1"
    "phyloseq::merge_samples" = "not in 0.1"
    "phyloseq::psmelt" = "not in 0.1"
    "phyloseq::plot_tree" = "not in 0.1"
    "phyloseq::plot_net" = "not in 0.1"
    ```
  - Run the Step 2 command: 72 passed in the prototype, against the 25 public
    functions of 1.18's tree.
  - Checked on the prototype: all 25 real docstrings parse.
    - 18 name one R function.
    - 5 name two: `pp.relative`, `pp.tax_glom`, `tl.alpha`, `tl.nmds`,
      `tl.pcoa`.
    - 2 say `none`: `datasets.toy`, `pp.filter_samples`.

    Removing the double backticks from `pl.richness`'s R line failed exactly
    `test_docstring_has_the_contract_sections[bt.pl.richness]`.
- [x] **Step 4: Docs.**
  - Create `docs/coming_from_r.md`. Its prose is written by hand; only the
    included table is generated (R8.4):
    ````markdown
    # Coming from R

    Every public biotapy function names its R equivalent in its docstring. This table is generated
    from those lines each time the docs are built, plus a short list of phyloseq accessors that are
    plain AnnData/TreeData code (`docs/_data/r_idioms.toml`). Rows marked "not in 0.1" have no
    biotapy equivalent yet.

    biotapy keeps samples as rows, so `tdata.X` is phyloseq's `otu_table` with
    `taxa_are_rows = FALSE`. Functions return new objects instead of changing yours; `tl` functions
    store their result in the object only with `inplace=True`.

    ```{include} generated/coming_from_r_table.md
    ```
    ````
  - `docs/conf.py`: replace the `exclude_patterns` line with:
    ```python
    # The Coming-from-R table is written by extensions/coming_from_r.py and included by coming_from_r.md.
    exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints", "generated/coming_from_r_table.md"]
    ```
  - `docs/index.md`: in the "User guide" toctree, add `coming_from_r.md`
    after `guide/index.md`.
  - Build with `uv run --group doc sphinx-build -W -b html docs docs/_build/html`,
    then check:
    - `grep -c "<tr" docs/_build/html/coming_from_r.html` prints `49` (48
      rows and the header);
    - `grep -c 'href="generated/biotapy.pl.bar.html#biotapy.pl.bar"' docs/_build/html/coming_from_r.html`
      prints `1`: the `{func}` roles resolved under `nitpicky`.
- [x] **Step 5: Knowledge.**
  - `contracts/function-shape.md`, Enforced by: "`tests/test_docstrings.py`
    (Phase 1, task 1.19) parses every public function's `R equivalent:` line."
    becomes: "`tests/test_docstrings.py` checks every public function's
    docstring. It must have exactly one `R equivalent:` line, parsed by
    `docs/extensions/coming_from_r.py`, with `Guide:` on the next line, and an
    `Examples` section. The same parser writes the Coming-from-R table at each
    docs build, with `docs/_data/r_idioms.toml` for the R calls that are plain
    AnnData code (task 1.19)."
  - Add a log line. Tick 1.19 here.
- [x] **Step 6: Gate and commit.**
  - The three gates.
  - Commit `docs: generate the Coming-from-R table and test every public docstring`,
    staging the five new files, `docs/conf.py`, `docs/index.md`, the concept,
    the roadmap and the log.

### Task 1.19b: Coverage gate

**Files:** modify `pyproject.toml` (`[tool.coverage]`), `tests/test_ci.py`,
`.knowledge/log.md`.

**Interfaces:** CI's `test` job already runs `uvx hatch run <env>:cov-report`,
which is hatch-test's `coverage report`. After this task that step exits
non-zero below 90% total line coverage (ruling 18).

- [x] **Step 1: Failing test.** In `tests/test_ci.py`:
  - Replace the header imports and `WORKFLOW` with:
    ```python
    import tomllib
    from pathlib import Path

    import yaml

    ROOT = Path(__file__).resolve().parents[1]
    WORKFLOW = yaml.safe_load((ROOT / ".github/workflows/test.yaml").read_text(encoding="utf-8"))
    ```
  - Append:
    ```python


    def test_coverage_below_90_percent_fails_the_test_job():
        coverage = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["tool"]["coverage"]
        assert coverage["report"]["fail_under"] == 90
        steps = WORKFLOW["jobs"]["test"]["steps"]
        assert any(":cov-report" in step.get("run", "") for step in steps)
    ```
  - `uv run --group test pytest tests/test_ci.py -q` fails with `KeyError: 'report'`.
- [x] **Step 2: Implement.** In `pyproject.toml`, after `run.source = [ "biotapy" ]`
  (pyproject-fmt moves the key there if it is written earlier), add:
  ```toml
  # rules.md R11.6 through a line-coverage proxy: coverage.py has no per-function view, so the whole
  # package (public and private modules) must reach 90%; CI's `coverage report` step fails below it.
  report.fail_under = 90
  ```
  Run the test again; it passes.
- [x] **Step 3: Measure.** Run
  `uv run --group test coverage run -m pytest -q && uv run --group test coverage report`
  and read the per-file table:
  - `TOTAL` must be at least 90; it was 99% before this slice;
  - every `pl` file should be at or near 100%.

  Name any file under 90% in the report. `datasets/_remote.py` is expected at
  89%, because its downloads run only in the network job.
- [x] **Step 4: Gate and commit.** Add a log line and tick 1.19b. Run the
  three gates. Commit `ci: fail the test job below 90% line coverage`,
  staging `pyproject.toml`, `tests/test_ci.py`, the roadmap and the log.

### Task 1.20: Docs, the vignette notebook and the docs CI job

**Files:**
- Create: `docs/tutorials/getting_started.md`, `docs/tutorials/phyloseq_analysis.md`
- Modify:
  - `docs/conf.py` (notebook execution)
  - `docs/tutorials/quick_tour.md` (drop the page override), `docs/tutorials/index.md`
  - `.github/workflows/test.yaml` (the `docs` job; `check` needs it)
  - `tests/test_ci.py`
  - `.knowledge/log.md`
- The exit gate's notebook path is already updated to
  `docs/tutorials/phyloseq_analysis.md` (ruling 13).

**Interfaces:**
- **Consumes:**
  - `bt.pl.*` (1.18); `tl.alpha`, `beta`, `unifrac`, `pcoa`, `nmds`;
    `pp.filter_features`;
  - `datasets.global_patterns`, `enterotype`, `esophagus`, through pooch.
- **Produces:**
  - the CI job `docs`, required by `check`;
  - the pages `tutorials/getting_started` and `tutorials/phyloseq_analysis`,
    which 1.23's README links.

- [x] **Step 1: Failing tests.** Append to `tests/test_ci.py`:
  ```python


  def test_docs_job_builds_the_docs_with_the_pooch_cache():
      steps = WORKFLOW["jobs"]["docs"]["steps"]
      build = [step for step in steps if step.get("run", "").strip() == "uvx hatch run docs:build"]
      assert build and "BIOTAPY_DATA_DIR" in build[0]["env"]


  def test_docs_job_blocks_merges():
      assert "docs" in WORKFLOW["jobs"]["check"]["needs"]
  ```
  `uv run --group test pytest tests/test_ci.py -q` fails with `KeyError: 'docs'`.
- [x] **Step 2: The CI job.** In `.github/workflows/test.yaml`:
  - Insert before the `knowledge-touched` job's comment:
    ```yaml
      # Builds the docs as Read the Docs does, executing every notebook (Phase 1 exit gate): the
      # phyloseq vignette downloads GlobalPatterns, enterotype and esophagus through the network
      # job's pooch cache.
      docs:
        runs-on: ubuntu-latest
        steps:
          - uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1 # v7.0.1
            with:
              filter: blob:none
              fetch-depth: 0
              persist-credentials: false
          - name: Install uv
            uses: astral-sh/setup-uv@bec219d24cd3e171d82865faccec33120bb574f4 # v10.1.0
            with:
              python-version: "3.13"
          - name: Cache pooch downloads
            uses: actions/cache@55cc8345863c7cc4c66a329aec7e433d2d1c52a9 # v6.1.0
            with:
              path: ${{ github.workspace }}/.pooch
              key: pooch-${{ hashFiles('src/biotapy/datasets/_remote.py') }}
          - name: Build the docs, executing every notebook
            env:
              BIOTAPY_DATA_DIR: ${{ github.workspace }}/.pooch
            run: uvx hatch run docs:build

    ```
  - In `check`'s `needs`, add `- docs` after `- network`.
  - The job runs the command Read the Docs runs (`.readthedocs.yaml`), on the
    same Python (3.13).
  - It shares the network job's cache key. When both jobs miss, the second
    save only logs "Unable to reserve cache".
  - The prototype's zizmor 1.24.1 run on this file reported no findings.
  - Run the Step 1 command: 8 passed (5 from before this slice, 1 from 1.19b).
- [x] **Step 3: Execution settings.**
  - In `docs/conf.py`, replace `nb_execution_mode = "off"` with:
    ```python
    # Every notebook runs at build time; a failing cell fails the build. The cache lives in docs/_build.
    nb_execution_mode = "cache"
    nb_execution_raise_on_error = True
    ```
  - In `docs/tutorials/quick_tour.md`, delete the three front-matter lines
    `mystnb:`, `  execution_mode: force` and `  execution_raise_on_error: true`
    (ruling 14).
  - Keep myst-nb's default 30 s cell timeout: the slowest notebook ran in
    11.4 s in total.
- [x] **Step 4: Notebooks.**
  - Create `docs/tutorials/getting_started.md`:
    ````markdown
    ---
    jupytext:
      text_representation:
        extension: .md
        format_name: myst
        format_version: 0.13
        jupytext_version: 1.16.4
    kernelspec:
      display_name: Python 3
      language: python
      name: python3
    ---

    # Getting started

    Install biotapy with `pip install biotapy` (Python 3.12 or newer). This page downloads
    phyloseq's GlobalPatterns dataset once (435 kB) and draws its first ordination:

    ```{code-cell} ipython3
    import biotapy as bt

    tdata = bt.datasets.global_patterns()  # 26 samples x 19,216 OTUs, with taxonomy and a tree
    bt.tl.beta(tdata, inplace=True)  # Bray-Curtis distances in obsp["braycurtis"]
    bt.tl.pcoa(tdata, inplace=True)  # coordinates in obsm["X_pcoa"]
    bt.pl.ordination(tdata, color="SampleType");
    ```

    Each `tl` call with `inplace=True` stores its result in the object; `pl` draws what is stored.
    Alpha diversity works the same way:

    ```{code-cell} ipython3
    bt.tl.alpha(tdata, metrics=["shannon"], inplace=True)  # obs["alpha_shannon"]
    bt.pl.richness(tdata, "shannon", x="SampleType");
    ```

    Next: the [quick tour](quick_tour.md) of the data model, and the
    [phyloseq vignette](phyloseq_analysis.md) redone with biotapy.
    ````
  - Create `docs/tutorials/phyloseq_analysis.md`:
    ````markdown
    ---
    jupytext:
      text_representation:
        extension: .md
        format_name: myst
        format_version: 0.13
        jupytext_version: 1.16.4
    kernelspec:
      display_name: Python 3
      language: python
      name: python3
    ---

    # The phyloseq analysis vignette in biotapy

    This notebook redoes the sections of phyloseq's
    [analysis vignette](https://github.com/joey711/phyloseq/blob/master/vignettes/phyloseq-analysis.Rmd)
    that biotapy 0.1 covers, on the same three datasets: GlobalPatterns, enterotype and esophagus.
    Each section names the R chunk it follows. Everything runs in biotapy; nothing is read from R.

    **Not in 0.1**, so left out:

    - `plot_tree` (exploratory tree plots) and `plot_net` (sample networks);
    - correspondence analysis (`ordinate(..., "CCA")`) and DPCoA, with their scree, species and biplot plots;
    - the `ACE` richness estimator, and `betadiver` distances such as `distance(esophagus, "g")`;
    - `hclust` dendrograms: SciPy's `scipy.cluster.hierarchy.linkage(..., method="average")` does it;
    - multiple testing and differential abundance (biotapy 0.3).

    ```{code-cell} ipython3
    import matplotlib.pyplot as plt
    import numpy as np

    import biotapy as bt
    ```

    ## Data

    `data(GlobalPatterns)`; `prune_taxa(taxa_sums(GP) > 0, GP)`; a `human` variable:

    ```{code-cell} ipython3
    global_patterns = bt.datasets.global_patterns()
    gp = bt.pp.filter_features(global_patterns, min_total=1)  # taxa_sums > 0 for counts
    gp.obs["human"] = gp.obs["SampleType"].isin(["Feces", "Mock", "Skin", "Tongue"])
    gp.shape
    ```

    ## Richness

    `plot_richness(GP, "human", "SampleType", measures = alpha_meas)`, one panel per measure, then a
    box plot per group (`+ geom_boxplot()`). InvSimpson is `1 / (1 - simpson)`.

    ```{code-cell} ipython3
    metrics = ["observed_features", "chao1", "shannon", "simpson"]
    bt.tl.alpha(gp, metrics=metrics, inplace=True)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4), layout="constrained")
    for ax, metric in zip(axes, metrics):
        bt.pl.richness(gp, metric, x="human", color="SampleType", ax=ax)
        values = gp.obs[f"alpha_{metric}"]
        ax.boxplot([values[~gp.obs["human"]], values[gp.obs["human"]]], positions=[0, 1], manage_ticks=False)
        ax.get_legend().remove()
    fig.legend(*axes[0].get_legend_handles_labels(), title="SampleType", loc="outside right center");
    ```

    ## Bar plots

    `data(enterotype)`; a rank-abundance bar plot of the 30 most abundant genera, averaged over
    samples (`barplot(sort(taxa_sums(enterotype), TRUE)[1:30] / nsamples(enterotype))`):

    ```{code-cell} ipython3
    enterotype = bt.datasets.enterotype()
    mean_abundance = enterotype.X.sum(axis=0).A1 / enterotype.n_obs
    top = np.argsort(-mean_abundance, kind="stable")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(enterotype.var_names[top[:30]], mean_abundance[top[:30]])
    ax.tick_params(axis="x", labelrotation=90)
    ```

    `rank_names(enterotype)` and, after keeping the 10 most abundant genera,
    `sample_variables(ent10)`:

    ```{code-cell} ipython3
    ent10 = enterotype[:, top[:10]].copy()  # prune_taxa(TopNOTUs, enterotype)
    list(enterotype.var.columns), list(ent10.obs.columns)
    ```

    `plot_bar(ent10, "SeqTech", fill = "Enterotype", facet_grid = ~Genus)`: one axes per genus,
    bars per sequencing technology, segments per enterotype. The unassigned genus is `NA`.

    ```{code-cell} ipython3
    fig, axes = plt.subplots(1, 10, sharey=True, figsize=(20, 4), layout="constrained")
    for ax, feature, genus in zip(axes, ent10.var_names, ent10.var["genus"].fillna("NA")):
        bt.pl.bar(ent10[:, [feature]], "Enterotype", x="SeqTech", ax=ax)
        ax.set_title(genus)
        ax.get_legend().remove()
    fig.legend(*axes[0].get_legend_handles_labels(), title="Enterotype", loc="outside right center");
    ```

    Later in the vignette, `plot_bar(GP, x = "human", fill = "SampleType", facet_grid = ~Phylum)` on
    the 200 most abundant OTUs of the five most abundant phyla among them:

    ```{code-cell} ipython3
    top200 = global_patterns[:, np.argsort(-global_patterns.X.sum(axis=0).A1, kind="stable")[:200]].copy()
    phylum_sums = top200.to_df().T.groupby(top200.var["phylum"].to_numpy()).sum().sum(axis=1)
    top5 = phylum_sums.sort_values(ascending=False).index[:5]
    gp200 = top200[:, top200.var["phylum"].isin(top5).to_numpy()].copy()
    gp200.obs["human"] = gp200.obs["SampleType"].isin(["Feces", "Mock", "Skin", "Tongue"])
    fig, axes = plt.subplots(1, 5, sharey=True, figsize=(20, 4), layout="constrained")
    for ax, phylum in zip(axes, top5):
        bt.pl.bar(gp200[:, (gp200.var["phylum"] == phylum).to_numpy()], "SampleType", x="human", ax=ax)
        ax.set_title(phylum)
        ax.get_legend().remove()
    fig.legend(*axes[0].get_legend_handles_labels(), title="SampleType", loc="outside right center");
    ```

    ## Heat map

    `plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family")` on the Crenarchaeota. phyloseq
    orders samples and taxa by their angle on an NMDS it computes; `pl.heatmap` keeps the order it
    is given, so the NMDS and the sort come first. Taxa are placed by the weighted average of the
    sample positions, as vegan's species scores. The NMDS differs from R's, so the order does too.

    ```{code-cell} ipython3
    gpac = global_patterns[:, (global_patterns.var["phylum"] == "Crenarchaeota").to_numpy()].copy()
    bt.tl.beta(gpac, inplace=True)
    bt.tl.nmds(gpac, seed=0, inplace=True)
    samples = gpac.obsm["X_nmds"]
    counts = gpac.X.toarray()
    taxa = counts.T @ samples / counts.sum(axis=0)[:, None]
    by_angle = gpac[np.argsort(np.arctan2(samples[:, 1], samples[:, 0])), np.argsort(np.arctan2(taxa[:, 1], taxa[:, 0]))]
    fig, ax = plt.subplots(figsize=(8, 10))
    bt.pl.heatmap(by_angle, ax=ax)
    ax.set_xticks(range(by_angle.n_obs), labels=by_angle.obs["SampleType"], rotation=90)
    ax.set_yticks(range(by_angle.n_vars), labels=by_angle.var["family"].fillna("NA"), fontsize=5);
    ```

    ## Ordination: PCoA on unweighted UniFrac

    The vignette loads a precomputed `UniFrac(GlobalPatterns)` because it is slow in R; here it
    takes about a second. `ordinate(GlobalPatterns, "PCoA", GPUF)`, then `plot_scree`:

    ```{code-cell} ipython3
    bt.tl.unifrac(global_patterns, inplace=True)
    bt.tl.pcoa(global_patterns, distance="unweighted_unifrac", n_components=25, inplace=True)
    share = global_patterns.uns["biotapy"]["pcoa"]["proportion_explained"]
    print(f"The first three axes explain {100 * share[:3].sum():.0f}%, the fourth another {100 * share[3]:.0f}%.")
    bt.pl.scree(global_patterns);
    ```

    Figure 5 of the Global Patterns article on two plots: axes 1 and 2, then 1 and 3.

    ```{code-cell} ipython3
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    bt.pl.ordination(global_patterns, color="SampleType", ax=axes[0]).get_legend().remove()
    bt.pl.ordination(global_patterns, components=(1, 3), color="SampleType", ax=axes[1]);
    ```

    ## Ordination: NMDS on unweighted UniFrac

    `ordinate(GlobalPatterns, "NMDS", GPUF)`:

    ```{code-cell} ipython3
    bt.tl.nmds(global_patterns, distance="unweighted_unifrac", seed=0, inplace=True)
    bt.pl.ordination(global_patterns, basis="nmds", color="SampleType");
    ```

    ## Distances

    `distance(esophagus, "bray")`, `"wunifrac"` and `"jaccard"`. phyloseq's `"jaccard"` is vegan's
    quantitative Jaccard; biotapy's is presence/absence, phyloseq's
    `distance(esophagus, "jaccard", binary = TRUE)`. The `betadiver` method `"g"` is not in 0.1.

    ```{code-cell} ipython3
    esophagus = bt.datasets.esophagus()
    bt.tl.beta(esophagus)
    ```

    ```{code-cell} ipython3
    bt.tl.unifrac(esophagus, weighted=True)
    ```

    ```{code-cell} ipython3
    bt.tl.beta(esophagus, metric="jaccard")
    ```
    ````
  - `docs/tutorials/index.md`: the toctree becomes `getting_started`,
    `quick_tour`, `phyloseq_analysis`.
  - How the vignette maps to phyloseq's chunks (research table, 39 chunks):
    - In scope, and redone: 4-7, 10-16, 21-27, the data steps of 28, 34
      and 37.
    - Listed as "not in 0.1": 8-9 (`plot_tree` and its Chlamydiae subset),
      18 (`plot_net`), 29-33 (CCA), 35-36 (DPCoA), the `"g"` method of 37,
      and 39 (`hclust`).
    - Left out as cosmetics, pointers or unevaluated code: 1-3, 17, 19, 20
      and 38.
  - Two differences from R are stated in the notebook:
    - phyloseq's `distance(esophagus, "jaccard")` is vegan's quantitative
      Jaccard, and biotapy's is its `binary = TRUE`;
    - the heat map's order follows biotapy's NMDS, not R's.
  - `manage_ticks=False` on `ax.boxplot` keeps the `False`/`True` group
    labels, which `boxplot` otherwise replaces with `0`/`1`. The prototype's
    first render showed that bug.
  - The faceted cells use `layout="constrained"` with
    `fig.legend(..., loc="outside right center")`. With plain
    `loc="center right"`, the shared legend covered the last panel.
- [x] **Step 5: Build and read the output.**
  - Build with `uv run --group doc sphinx-build -W -b html docs docs/_build/html`,
    with `BIOTAPY_DATA_DIR` set to a pooch cache (default: the per-user
    cache). CI runs the equivalent `uvx hatch run docs:build`.
  - Expect three `Executed notebook in ... seconds [mystnb]` lines (5.1 s,
    13.5 s and 2.3 s in the prototype) and `build succeeded.`
  - `grep -c '<img' docs/_build/html/tutorials/phyloseq_analysis.html` prints
    `8`.
  - Open the eight figures in `docs/_build/html/_images/` and check them
    against the vignette's:
    - the richness panels with `False`/`True` under boxes;
    - ten enterotype panels, one titled `NA`;
    - five phylum panels;
    - the heat map;
    - the scree;
    - PCoA axes 1-2 and 1-3 with `[%]` labels;
    - the NMDS with its stress note.
  - Check that the gate bites: add a cell `bt.pl.scree(bt.datasets.esophagus())`
    to the vignette and rebuild. The build must exit 2 with
    `WARNING: Executing notebook failed: CellExecutionError [mystnb.exec]`.
    The cell raises the `KeyError` naming `bt.tl.pcoa`, and myst-nb does not
    print tracebacks. Remove the cell and rebuild.
  - Kernel lines such as `[IPKernelApp] WARNING | Kernel is running over TCP
    without encryption` are the kernel's stderr, not Sphinx warnings; `-W`
    ignores them.
- [x] **Step 6: Gate and commit.**
  - Add a log line and tick 1.20. Run the three gates.
  - Commit `docs: add getting started and the phyloseq vignette, executed on every build`,
    staging the two new notebooks, `docs/conf.py`, `docs/tutorials/quick_tour.md`,
    `docs/tutorials/index.md`, `.github/workflows/test.yaml`,
    `tests/test_ci.py`, the roadmap and the log.
  - The job's first CI run is at Checkpoint D's PR.

### Task 1.21: asv benchmarks

**Files:**
- Create:
  - `benchmarks/asv.conf.json`
  - `benchmarks/benchmarks/__init__.py` (empty), `_data.py`, `pp.py`, `tl.py`
  - `docs/performance.md`
- Modify:
  - `pyproject.toml` (asv in the dev group)
  - `.github/workflows/test.yaml` (`asv check` in `lint`), `tests/test_ci.py`
  - `docs/index.md` (the Project toctree)
  - `.knowledge/contracts/engine-parity.md`, `.knowledge/modules/tl.md`
  - `.knowledge/log.md`
- `benchmarks/` needs no sdist change: `build.targets.sdist.include` is a whitelist
  (`/CHANGELOG.md`, `/docs`, `/src/biotapy`, `/tests`). On the prototype,
  `uv build --sdist` held no `benchmarks/` or `.asv` entry.

**Interfaces:**
- **Consumes:**
  - `bt.pp.relative`, `tax_glom`, `rarefy` and `bt.tl.alpha`, `beta`;
  - `biotapy._core.make_treedata` and `tree_from_edges`, the only
    TreeData/tree constructors ([tree-access](/contracts/tree-access.md)).
    Benchmarks are development code outside the package, as tests are.
- **Produces:**
  - `benchmarks/benchmarks/_data.py:synthetic(n_obs=5_000, n_vars=50_000) -> TreeData`,
    the dataset [engine-parity](/contracts/engine-parity.md) statement 5 names;
  - the baseline table in `docs/performance.md`, which ticks the exit gate's
    "asv baselines recorded".

- [x] **Step 1: Dependency.**
  - In `[dependency-groups] dev`, add `"asv>=0.6.6",` before
    `"import-linter>=2.15",`.
  - 0.6.6 is the version prototyped: it ships the `uv` environment plugin
    (`asv/plugins/uv.py`), and its wheels cover Linux, macOS and Windows for
    every CPython the hatch-test matrix installs the dev group on.
  - Run `uv sync --group dev --group test --group doc`, then confirm that
    `uv run --group dev python -c "import asv; print(asv.__version__)"`
    prints `0.6.6` or later.
- [x] **Step 2: Failing test.** Append to `tests/test_ci.py`:
  ```python


  def test_lint_job_imports_the_benchmarks():
      steps = WORKFLOW["jobs"]["lint"]["steps"]
      check = [step for step in steps if step.get("run", "").strip() == "uv run --group dev asv check --python=same"]
      assert check and check[0]["working-directory"] == "benchmarks"
  ```
  It fails: `check` is empty.
- [x] **Step 3: The suite.**
  - `benchmarks/asv.conf.json`, as plain JSON:
    - asv 0.6 would accept JSON5 comments, but the repo's biome-format hook
      rejects them; it failed the prototype's first draft.
    - `"environment_type": "uv"` is asv's own uv plugin: uv creates each
      environment and pip installs into it. `--python=same` below skips it
      and uses the current environment.
    - The three `../.asv/` directories sit under the repository's git-ignored
      `/.asv/`.
    - `"branches": ["master"]` must name an existing branch. The prototype's
      scratch repository, on `main`, failed with
      `Unknown branch master in configuration`.
    ```json
    {
        "version": 1,
        "project": "biotapy",
        "project_url": "https://github.com/pedrocr83/biotapy",
        "repo": "..",
        "branches": ["master"],
        "environment_type": "uv",
        "pythons": ["3.13"],
        "benchmark_dir": "benchmarks",
        "env_dir": "../.asv/env",
        "results_dir": "../.asv/results",
        "html_dir": "../.asv/html"
    }
    ```
  - `benchmarks/benchmarks/_data.py`:
    ```python
    """Synthetic benchmark data: a sparse count table with taxonomy and a balanced tree."""

    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    # Benchmarks build their TreeData with biotapy's own constructors, as the readers do (contracts/tree-access).
    from biotapy._core import TreeData, make_treedata, tree_from_edges

    N_OBS, N_VARS, DENSITY, SEED = 5_000, 50_000, 0.02, 0
    # Features per group at each rank: 1,000 genera, 100 families, 20 orders, 10 classes and 5 phyla at 50,000 features.
    _RANK_SIZES = {"phylum": 10_000, "class": 5_000, "order": 2_500, "family": 500, "genus": 50}


    def synthetic(n_obs: int = N_OBS, n_vars: int = N_VARS) -> TreeData:
        """``n_obs`` x ``n_vars`` counts from 1 to 99 at 2% density, seed 0, with taxonomy and a tree."""
        rng = np.random.default_rng(SEED)
        X = sp.random(
            n_obs, n_vars, density=DENSITY, format="csr", rng=rng, data_rvs=lambda size: rng.integers(1, 100, size)
        )
        features = [f"f{j}" for j in range(n_vars)]
        ranks = {rank: [f"{rank[0]}{j // size}" for j in range(n_vars)] for rank, size in _RANK_SIZES.items()}
        var = pd.DataFrame({"kingdom": "Bacteria"} | ranks, index=features)
        obs = pd.DataFrame(index=[f"s{i}" for i in range(n_obs)])
        tree = tree_from_edges(_balanced_edges(features, rng))
        return make_treedata(X, obs=obs, var=var, tree=tree, x_kind="counts", source="benchmarks")


    def _balanced_edges(tips: list[str], rng: np.random.Generator) -> list[tuple[str, str, float]]:
        """Join neighbours pairwise until one root is left; branch lengths uniform in [0.01, 0.1)."""
        level, edges = list(tips), []
        while len(level) > 1:
            pairs = [level[i : i + 2] for i in range(0, len(level), 2)]
            level = [f"n{len(edges)}_{k}" for k in range(len(pairs))]
            for parent, pair in zip(level, pairs, strict=True):
                edges += [(parent, child, float(rng.uniform(0.01, 0.1))) for child in pair]
        return edges
    ```
  - `benchmarks/benchmarks/pp.py`:
    ```python
    """Preprocessing at 5,000 samples x 50,000 features."""

    import biotapy as bt
    from biotapy._core import TreeData

    from ._data import synthetic


    class Preprocessing:
        """``relative``, ``tax_glom`` and ``rarefy`` on the full synthetic table."""

        number, repeat, rounds, timeout = 1, 3, 1, 300

        def setup_cache(self) -> TreeData:
            """Build the table once; asv pickles it for every benchmark."""
            return synthetic()

        def time_relative(self, tdata: TreeData) -> None:
            """``pp.relative``."""
            bt.pp.relative(tdata)

        def time_tax_glom(self, tdata: TreeData) -> None:
            """``pp.tax_glom`` to 1,000 genera."""
            bt.pp.tax_glom(tdata, "genus")

        def time_rarefy(self, tdata: TreeData) -> None:
            """``pp.rarefy`` to the smallest sample depth."""
            bt.pp.rarefy(tdata, seed=0)
    ```
  - `benchmarks/benchmarks/tl.py`:
    ```python
    """Diversity at 5,000 samples x 50,000 features, and Bray-Curtis on sample subsets."""

    import biotapy as bt
    from biotapy._core import TreeData

    from ._data import synthetic


    class Alpha:
        """``tl.alpha``: rows densified in chunks of 2**20 values."""

        number, repeat, rounds, timeout = 1, 1, 1, 900

        def setup_cache(self) -> TreeData:
            """Build the table once; asv pickles it for every benchmark."""
            return synthetic()

        def time_alpha(self, tdata: TreeData) -> None:
            """The four default metrics."""
            bt.tl.alpha(tdata)

        def peakmem_alpha(self, tdata: TreeData) -> None:
            """Peak memory of the four default metrics: chunking keeps it far below one dense copy (2 GB)."""
            bt.tl.alpha(tdata)

        def time_faith_pd(self, tdata: TreeData) -> None:
            """Faith PD on a 100,000-node tree, which scikit-bio re-indexes for each chunk of 20 samples."""
            bt.tl.alpha(tdata, metrics=["faith_pd"])


    class Beta:
        """``tl.beta`` Bray-Curtis: pairwise, so time grows with the square of the samples."""

        number, repeat, rounds, timeout = 1, 1, 1, 1800
        params, param_names = [500, 1_000, 2_000], ["n_obs"]

        def setup_cache(self) -> TreeData:
            """Build the table once; asv pickles it for every benchmark."""
            return synthetic()

        def time_beta(self, tdata: TreeData, n_obs: int) -> None:
            """The first ``n_obs`` samples, all 50,000 features."""
            bt.tl.beta(tdata[:n_obs].copy())


    class BetaFullSize:
        """Peak memory of ``tl.beta`` on all 5,000 samples: X densified once (2 GB) plus the distances."""

        number, repeat, rounds, timeout = 1, 1, 1, 1800

        def setup_cache(self) -> TreeData:
            """Build the table once; asv pickles it for every benchmark."""
            return synthetic()

        def peakmem_beta(self, tdata: TreeData) -> None:
            """One call, about 6 minutes."""
            bt.tl.beta(tdata)
    ```
  - Why the suite looks as it does:
    - `setup_cache` builds the table once per class, and asv passes the
      pickled result to every benchmark in it.
    - `number, repeat, rounds = 1, ...`: each call takes seconds, so asv's
      automatic calibration would only multiply the run time.
    - `timeout` is raised: asv's default is 60 s.
    - `beta` runs on the first 500, 1,000 and 2,000 samples. At all 5,000 one
      call takes about 6 minutes and needs X densified (2 GB), so
      `BetaFullSize` measures it once, for peak memory only.
    - ruff lints `benchmarks/` like any Python file. The docstrings satisfy
      its `D` rules, and asv shows them as descriptions.
  - In `.github/workflows/test.yaml`, in the `lint` job, add after the prek step:
    ```yaml
          # Imports and discovers the asv suite without running it (docs/performance.md), so an API
          # change that breaks a benchmark module fails here; asv check needs no machine file.
          - name: Import the asv benchmarks
            working-directory: benchmarks
            run: uv run --group dev asv check --python=same
    ```
  - Run `cd benchmarks && uv run --group dev asv check --python=same`. It
    prints `No problems found.`
  - Run `uv run --group test pytest tests/test_ci.py -q`: 9 passed.
  - `contracts/engine-parity.md` statement 5: "an `asv` benchmark on the
    5,000 x 50,000 sparse synthetic dataset" becomes "an `asv` benchmark on
    the 5,000 x 50,000 sparse synthetic dataset
    (`benchmarks/benchmarks/_data.py:synthetic`, baselines in
    `docs/performance.md`)". Add a log line.
  - Run the three gates. Commit `perf: add asv benchmarks on a synthetic 5,000 x 50,000 table`,
    staging `pyproject.toml`, the five `benchmarks/` files,
    `.github/workflows/test.yaml`, `tests/test_ci.py`, the concept and the log.
    The baseline is measured on this commit.
- [x] **Step 4: Run the baseline.** It takes about 11 minutes and 3 GB of free
  memory. Close other heavy work first, since the numbers are the record.
  With `git status --short` empty, from the repository root:
  ```bash
  cd benchmarks
  uv run --group dev env HOME="$PWD/../.asv" asv machine --yes
  uv run --group dev env HOME="$PWD/../.asv" asv run --python=same --set-commit-hash "$(git rev-parse HEAD)" --show-stderr
  uv run --group dev env HOME="$PWD/../.asv" asv show "$(git rev-parse HEAD)"
  ```
  - `HOME` points asv's machine file into the git-ignored `.asv/` (ruling
    17), while `uv run` keeps the normal `HOME`.
  - Expect 8 benchmarks, each `ok`. The prototype's results are in the design
    section's asv facts; a result more than twice those, or a failure, is
    reported before it is recorded (R14.1).
  - `git status --short` must show nothing new: `.asv/` is ignored.
- [x] **Step 5: The performance page.**
  - Create `docs/performance.md` from the prototype's page below. Replace the
    commit, date, machine line and every number with Step 4's `asv show`
    output: the page records this run, not the prototype's.
    ````markdown
    # Performance

    Baselines measured with [asv](https://asv.readthedocs.io/) on a synthetic table: 5,000 samples x
    50,000 features, 2% non-zero counts from 1 to 99, a taxonomy of 1,000 genera in 5 phyla and a
    balanced tree of 100,006 nodes, all from seed 0 (`benchmarks/benchmarks/_data.py`). They are a
    record, not a gate: release 0.1 sets no speed target.

    Measured on commit `6498f54`, 2026-09-28: a 16-core laptop with 62 GB of RAM, Linux, Python
    3.13, in the development environment (`asv run --python=same`).

    | Benchmark | Result |
    |---|---|
    | `pp.relative` | 0.90 s |
    | `pp.tax_glom` to genus | 2.90 s |
    | `pp.rarefy` to the smallest depth | 5.44 s |
    | `tl.alpha`, the four default metrics | 1.39 s |
    | `tl.alpha`, the four default metrics, peak memory | 478 MB |
    | `tl.alpha`, `faith_pd` | 43.5 s |
    | `tl.beta`, Bray-Curtis, 500 samples | 5.62 s |
    | `tl.beta`, Bray-Curtis, 1,000 samples | 16.2 s |
    | `tl.beta`, Bray-Curtis, 2,000 samples | 54.2 s |
    | `tl.beta`, Bray-Curtis, 5,000 samples, peak memory | 2.73 GB |

    - `tl.beta` compares every pair of samples, so its time grows with the square of the samples:
      about 6 minutes at 5,000. Its memory is one dense copy of `X` (8 bytes x 5,000 x 50,000 =
      2 GB) plus the distances, as its docstring says.
    - `tl.alpha` densifies at most 2**20 values at a time, so its peak stays far below one dense
      copy of `X`.
    - `faith_pd` spends almost all its time in scikit-bio, which re-indexes the 100,006-node tree for
      each chunk of 20 samples; converting the tree once takes 0.46 s.

    ## Running the benchmarks

    ```bash
    uv sync --group dev
    cd benchmarks
    # asv keeps its machine description in ~/.asv-machine.json; HOME points it at the git-ignored .asv/.
    uv run --group dev env HOME="$PWD/../.asv" asv machine --yes
    uv run --group dev env HOME="$PWD/../.asv" asv run --python=same --set-commit-hash "$(git rev-parse HEAD)"
    uv run --group dev env HOME="$PWD/../.asv" asv show "$(git rev-parse HEAD)"
    ```

    `--python=same` runs in the current environment, and `--set-commit-hash` keeps the results, in
    `.asv/results`. The whole suite takes about 11 minutes and needs about 3 GB of free memory.
    `uv run --group dev asv check --python=same` imports the suite without running it; CI runs it.
    ````
  - `docs/index.md`: in the "Project" toctree, add `performance.md` after
    `design.md`.
- [x] **Step 6: Knowledge.**
  - `modules/tl.md`, add a Gotcha: "`faith_pd` on 5,000 x 50,000 takes about
    45 s (asv, Task 1.21). scikit-bio re-indexes and re-validates the tree on
    every `alpha_diversity` call, once per chunk of `2**20 // n_vars` samples,
    and that is about 97% of the time; converting the tree once
    (`get_skbio_tree`) takes 0.46 s. A fix needs a profile-driven perf task
    (rules.md R10.1), not a change here."
  - Add a log line. Tick 1.21 here.
- [x] **Step 7: Gate and commit.**
  - The three gates.
  - Commit `docs: record the asv baselines for 0.1`, staging
    `docs/performance.md`, `docs/index.md`, `.knowledge/modules/tl.md`, the
    roadmap and the log.

### Task 1.22: Knowledge

**Files:**
- Create: `.knowledge/modules/pl.md`
- Modify:
  - `.knowledge/modules/index.md`
  - the 18 concepts `bash scripts/knowledge_stale.sh` reports stale (listed in Step 2)
  - `.knowledge/log.md`

**Interfaces:**
- **Consumes:** the code of 1.18-1.21, committed.
- **Produces:** `bash scripts/knowledge_stale.sh --against HEAD` prints
  `21 current, 0 stale, 12 uncheckable`. Checkpoint D and 1.23 re-run it after
  their own commits.

- [x] **Step 1: The `pl` Module concept.** Create `.knowledge/modules/pl.md`.
  - Set `generated.by` to `claude-code/` plus the implementing model's id.
  - Set `generated.at` to the output of `date -u +%FT%TZ`.
  - Set `commit` to `git rev-parse --short HEAD`, the last code commit.
  ```markdown
  ---
  type: Module
  title: pl (plots)
  description: Plots of what tl and pp stored - stacked bars, heatmap, richness, ordination and scree - drawn with matplotlib on the given or a new Axes, computing nothing.
  resource: /src/biotapy/pl/
  paths: ["src/biotapy/pl/**"]
  tags: [pl, plots, matplotlib]
  generated: { by: claude-code/<implementing model id>, at: <date -u +%FT%TZ> }
  commit: <git rev-parse --short HEAD>
  status: stable
  ---

  # Responsibility

  Owns the `bt.pl.*` verbs that draw stored results:
  - `bar` and `heatmap` draw `X` or a layer;
  - `richness` draws `obs["alpha_<metric>"]`;
  - `ordination` draws `obsm["X_pcoa" | "X_nmds"]` with its `uns["biotapy"]`
    summary;
  - `scree` draws `uns["biotapy"]["pcoa"]`.

  It computes no diversity, distance or ordination; those are
  [tl](/modules/tl.md)'s. A missing slot is a `KeyError` naming the `tl` or
  `pp` call that writes it.

  # Entry points
  - `_abundance.py:bar` - stacked bars. `_segments` sums features per `fill`
    group (a `var` or an `obs` column); `_by_x` sums samples per `x` group.
  - `_abundance.py:heatmap` - `imshow` of the table densified once, features x
    samples, on a `LogNorm` scale with phyloseq's colours.
  - `_richness.py:richness` - one point per sample; NaN values are left out.
  - `_ordination.py:ordination` and `_ordination.py:scree` - both read through
    `_ordination.py:_stored`.
  - `_common.py` - helpers the topic files share: `new_axes`, `table`,
    `obs_groups`/`groups`, `scatter` and `label_ticks`.

  # Invariants
  - `import biotapy` never imports matplotlib. Annotations use a
    `TYPE_CHECKING` import, and functions import matplotlib when they draw
    (`tests/pl/test_init.py`).
  - Every function takes `ax=None` last and keyword-only, and returns the
    `Axes`. With `ax=None` it draws on a new pyplot figure
    (`_common.py:new_axes`), which notebooks display.
  - Nothing is written to the AnnData
    ([pure-by-default](/decisions/pure-by-default.md)); every test asserts the
    input unchanged.
  - Groups (`_common.py:groups`) keep a categorical's order and sort other
    values. Missing values are a last `NA` group, drawn grey. A numeric column
    raises `TypeError`, as `tl.permanova` does.
  - `bar` heights equal the sample (or `x` group) totals of the plotted table:
    features with a missing rank are a group, not dropped (a Hypothesis test).
  - `heatmap` keeps `obs`/`var` order and densifies the table once (rules.md
    R6.2); it raises `ValueError` when nothing is positive.
  - Past 250 names an axis gets no tick labels (`_common.py:MAX_LABELS`,
    phyloseq's `max.label`).

  # Dependencies
  - [core](/modules/core.md): `as_csr`, `sum_by`.
  - matplotlib `>=3.8`, a runtime dependency
    ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).
  - The slots [tl](/modules/tl.md) and `pp.relative` write
    ([data-model-slots](/contracts/data-model-slots.md)).

  # Verification

  `uv run --group test pytest tests/pl src/biotapy/pl -q`. The tutorials draw
  every function on GlobalPatterns, enterotype and esophagus in the docs build:
  `BIOTAPY_DATA_DIR=<cache> uv run --group doc sphinx-build -W -b html docs docs/_build/html`.

  # Gotchas
  - `ax.bar` adds one `Rectangle` per bar, one group after another, so tests
    reshape the `ax.patches` heights to bars x groups.
  - `LogNorm` masks zeros, which the colormap draws in its "bad" colour
    (black). On an all-zero table matplotlib would fail inside `colorbar` with
    `Invalid vmin or vmax`, so `heatmap` raises its own `ValueError` first.
  - `ax.boxplot` on the same axes replaces the group tick labels unless it
    gets `manage_ticks=False` (the vignette's richness cell).
  - `bar` with hundreds of `fill` groups is slow: GlobalPatterns by genus is
    984 groups and about 8 s. Aggregate first.
  - matplotlib 3.8.0-3.8.3 pin `numpy<2`, so with NumPy 2 the resolver takes
    3.8.4 or later.
  - `pl` has R equivalents but no golden tests: it draws numbers that `tl`'s
    golden tests check ([r-golden-parity](/contracts/r-golden-parity.md)
    statement 7).
  ```
  - Add to `modules/index.md`, after the `tl` line: `* [pl](pl.md) - Plots of
    what tl and pp stored - stacked bars, heatmap, richness, ordination and
    scree - drawn with matplotlib on the given or a new Axes, computing
    nothing.`
- [x] **Step 2: Refresh the stale concepts.**
  - `bash scripts/knowledge_stale.sh` reported 18 stale concepts before this
    slice (research, 2026-09-27):
    - the roadmap: `phase-0-foundation`, `phase-1-core`, `phase-2-function`,
      `phase-3-stats`, `phase-4-ml-multiomics`;
    - the modules: `core`, `io`, `pp`, `tl`;
    - the contracts: `data-model-slots`, `engine-parity`, `function-shape`,
      `module-boundaries`, `r-golden-parity`, `tree-access`;
    - the playbooks: `add-a-function`, `cut-a-release`,
      `regenerate-golden-files`.
  - `datasets` and `maintain-knowledge` stay current: nothing under their
    `paths` changes in 1D.
  - For each of the 18:
    1. List what moved under its globs:
       `git diff --stat <its commit>..HEAD -- <its paths>`.
    2. Read the concept, and open the file behind every statement those
       changes could touch (R0.2).
    3. Fix a statement only if the code now disagrees (R12.1). A concept can
       be correct and only need its `commit` bumped.
    4. Set `commit:` to `git rev-parse --short HEAD`, and `generated` to the
       implementing model and `date -u +%FT%TZ`. Never add `verified`
       (R12.3).
  - Edits already known:
    - `roadmap/phase-1-core.md`, the "For agentic workers" note: every task
      now has full TDD steps, so replace its second and third sentences with
      "Every task through 1.23 has full TDD steps."
    - `contracts/data-model-slots.md`: after the Slots table, add "`pl` reads
      these slots and writes none ([pl](/modules/pl.md))."
    - `contracts/function-shape.md`: statement 1's "except a documented
      `plot_kwargs` in `pl`" stays. It is a permission, and 0.1 uses none.
      Check that the Enforced-by line matches 1.19's edit.
    - `playbooks/add-a-function.md`, step 7: append "The `R equivalent:` line
      feeds the generated Coming-from-R page; `tests/test_docstrings.py` fails
      if it cannot be parsed."
  - Report every other edit in the task report, one line per concept: what
    changed, or "commit only".
- [x] **Step 3: Verify.**
  - `uv run --group test pytest tests/test_knowledge_bundle.py -q`: every
    concept has a `type` and is listed in its `index.md`, including `pl`.
  - `bash scripts/knowledge_stale.sh --against HEAD` must print
    `21 current, 0 stale, 12 uncheckable`.
    - The default `--against` is `origin/master`. Until the branch merges,
      that ref lacks these commits, so every refreshed concept reads stale
      against it; `--against HEAD` checks the branch as it will merge.
    - The 12 uncheckable are the 10 decisions, `phase-5-beyond` and
      `spec-review`, which have no `paths:` key (ruling 19).
- [x] **Step 4: Commit.** Add a log line naming every concept refreshed, and
  tick 1.22. Commit `docs(knowledge): add the pl module and refresh every Phase 1 concept`,
  staging the concepts, `modules/index.md`, the roadmap and the log.

### Checkpoint D - review slice 1D
- [x] Review the whole slice (superpowers:requesting-code-review) against:
  - every contract: function-shape, data-model-slots, module-boundaries,
    tree-access and r-golden-parity;
  - the pure-by-default decision;
  - the slice 1D review focus and rulings;
  - rules.md R3-R8.

  Then a fix pass, one commit per finding, each with a test. Also read the
  per-file coverage table from 1.19b's command, and report any file under 90%.
- [x] Knowledge after the fix pass: run
  `bash scripts/knowledge_stale.sh --against HEAD`.
  - For each concept a fix made stale: check it against the fix, then bump its
    `commit` to the new `HEAD`, in one commit
    `docs(knowledge): refresh concepts after the Checkpoint D fixes`.
  - The script must again print `21 current, 0 stale, 12 uncheckable`.
- [ ] Push `phase-1d` and open the PR, each after the user approves it
  (R13.3).
  - The PR's CI must be green: every hatch-test job (Python 3.12-3.14 on
    three OSes, the coverage gate included), `lint` (with `asv check`),
    `import-without-extras`, `network` (every golden test) and the new
    `docs` job.
  - The `docs` job is the exit gate's "executes in CI". Read its log for the
    three `Executed notebook` lines.
- [ ] After the user approves, merge with a merge commit, not a squash
  (maintain-knowledge: a squash drops the branch SHAs that concepts cite).
  Then:
  - `git switch master && git pull --ff-only && bash scripts/knowledge_stale.sh`
    prints 0 stale against `origin/master`;
  - `curl -s -o /dev/null -w "%{http_code}\n" https://biotapy.readthedocs.io/en/latest/tutorials/phyloseq_analysis.html`
    prints `200` once Read the Docs has built `master`.
    - Also check that the page shows the eight figures. This is the first
      build to download datasets on Read the Docs; the research found its
      builders have network access, which this confirms.
    - If the page is missing, read the build log on readthedocs.org and report
      it; do not work around it.
- [ ] Ask the user to review slice 1D before the release.

### Task 1.23: Release 0.1

Follows [cut-a-release](/playbooks/cut-a-release.md). Publishing to PyPI is
irreversible: a version can never be uploaded twice. Every outward step waits
for the user's explicit approval of that step (R13.3).

**Files:**
- Modify: `pyproject.toml` (version), `CHANGELOG.md`, `README.md`
- Modify (knowledge): the concepts the version bump makes stale,
  `.knowledge/playbooks/cut-a-release.md`,
  `.knowledge/roadmap/phase-1-core.md`, `.knowledge/roadmap/index.md`,
  `.knowledge/roadmap/phase-2-function.md`, `.knowledge/log.md`
- Proposals, written in the scratchpad and never committed or submitted: a
  conda-forge recipe note and a scverse ecosystem listing.

**Interfaces:**
- **Consumes:**
  - `master` after Checkpoint D;
  - `.github/workflows/release.yaml`, which runs on `release: published`:
    `uv build`, then `pypa/gh-action-pypi-publish` v1.14.2 in environment
    `pypi` with `id-token: write`.
- **Produces:** `biotapy 0.1.0` on PyPI, Phase 1 `phase_state: done`, and
  Phase 2 active.

- [ ] **Step 1: User action - trusted publisher.** Claude cannot see PyPI's
  settings. Ask the user to confirm, and stop until they do:
  - on pypi.org, project `biotapy`, Publishing lists a trusted publisher with
    owner `pedrocr83`, repository `biotapy`, workflow `release.yaml` and
    environment `pypi`. 0.0.1 was published through it, so it should exist
    and no longer be pending;
  - on GitHub, Settings > Environments has `pypi`. If it requires a reviewer,
    the user approves the deployment in the Actions run at Step 9.
- [ ] **Step 2: Check the publish action.** Run
  `gh api repos/pypa/gh-action-pypi-publish/releases/latest --jq .tag_name`.
  - If it prints anything newer than `v1.14.2`, stop and report it. A stale
    pin failed 0.0.1 on Metadata-Version 2.5 (the playbook's Common
    mistakes), so the Dependabot update must merge first.
  - On 2026-09-27 the research found v1.14.2 to be the latest.
- [ ] **Step 3: Branch.**
  `git switch master && git pull --ff-only && git switch -c release-0.1.0`.
- [ ] **Step 4: Version and changelog.**
  - `pyproject.toml`: `version = "0.0.1"` becomes `version = "0.1.0"`.
  - `CHANGELOG.md`: append to the `### Added` list under `## [Unreleased]`:
    ```markdown
    - `bt.pl.bar`, `bt.pl.heatmap`, `bt.pl.richness`, `bt.pl.ordination` and
      `bt.pl.scree`: phyloseq's plots, drawn from the results `tl` and `pp`
      store, returning matplotlib axes.
    - A Coming-from-R page, generated from every public docstring's R
      equivalent, and a test that keeps those docstrings parseable.
    - Tutorials: getting started, and phyloseq's analysis vignette redone with
      biotapy. Every notebook now runs on each docs build.
    - asv benchmarks on a synthetic 5,000 x 50,000 table, with the 0.1
      baselines in the docs.
    - New runtime dependency: `matplotlib>=3.8`, for `pl`, imported only when a
      plot is drawn.
    ```
  - Append to its `### Changed` list:
    ```markdown
    - CI builds the docs, executing every notebook, and fails when total line
      coverage drops below 90%.
    ```
  - Rename `## [Unreleased]` to `## [0.1.0] - ` followed by `date -u +%F`.
    Insert a new empty `## [Unreleased]` above it, as Keep a Changelog does.
- [ ] **Step 5: README.** The README is PyPI's project page. It still says
  0.0.1 is a placeholder and that 0.1 is coming. The parts true before the
  release were done on 2026-10-02 (the Datasets, Preprocessing, Tools and
  Plots bullets, `esophagus` in Example datasets and licensing, the data
  model and pure-by-default paragraphs, and the Coming from R paragraph with
  its two links); check that they still hold, then change only:
  - The Status section's first line and its closing paragraph:
    - "**biotapy is pre-release.** The API can still change without notice."
      becomes "**biotapy 0.1 is an early release.** The API can still change
      between minor versions."
    - "What works today" becomes "What 0.1 does".
    - The "Release 0.1 will publish all of this on PyPI. ..." paragraph
      becomes: "Next, in 0.2: functional profiles from HUMAnN, PICRUSt2 and
      MetaPhlAn. See the [roadmap][roadmap]; no dates are promised."
    - Its link definition becomes
      `[roadmap]: https://github.com/pedrocr83/biotapy/blob/master/.knowledge/roadmap/index.md`.
  - Installation: replace everything from "biotapy is not functional on PyPI
    yet" through the closing `-->` of the commented-out block with:
    ````markdown
    ```bash
    pip install biotapy
    ```

    or, with [uv][]:

    ```bash
    uv add biotapy
    ```

    The development version installs straight from GitHub:
    `pip install git+https://github.com/pedrocr83/biotapy.git`.
    ````
  - `docs/tutorials/getting_started.md`: its install sentence names the GitHub
    install line until now; replace it with `pip install biotapy`.
- [ ] **Step 6: Build check.**
  - `rm -rf dist && uv build && uvx twine check --strict dist/*`: both files
    print `PASSED` (0.0.1's sdist and wheel did).
  - `tar tzf dist/biotapy-0.1.0.tar.gz | grep -c benchmarks` prints `0`.
  - `unzip -p dist/biotapy-0.1.0-py3-none-any.whl 'biotapy-0.1.0.dist-info/METADATA' | grep -E "^(Version|Requires-Dist): (0.1.0|matplotlib)"`
    prints `Version: 0.1.0` and `Requires-Dist: matplotlib>=3.8`.
  - `rm -rf dist`.
- [ ] **Step 7: Gate and commit.** Run the three gates. Commit
  `chore: release 0.1.0`, staging `pyproject.toml`, `CHANGELOG.md` and
  `README.md`.
- [ ] **Step 8: Knowledge, second commit.**
  - Run `bash scripts/knowledge_stale.sh --against HEAD`. The version bump
    re-stales the concepts whose `paths` hold `pyproject.toml` or
    `CHANGELOG.md`: `module-boundaries`, `tree-access`, `phase-0-foundation`
    and `cut-a-release`.
    - Check each against the bump, which changes nothing they state, and bump
      `commit` to `HEAD`.
  - `playbooks/cut-a-release.md`: add Step 2b: "Update `README.md` wherever
    it describes the previous release: it is PyPI's project page, and 0.1.0
    replaced 0.0.1's placeholder text." Add any other step this release
    needed.
  - `roadmap/phase-1-core.md`: tick the first five exit-gate items. Each was
    proven on Checkpoint D's PR:
    - the `docs` job;
    - the `network` job's golden tests;
    - `tests/test_coming_from_r.py`;
    - `docs/performance.md`;
    - the coverage gate.
  - Add a log line and tick Steps 1-8 here.
  - `bash scripts/knowledge_stale.sh --against HEAD` prints
    `21 current, 0 stale, 12 uncheckable`.
  - Commit `docs(knowledge): refresh concepts for the 0.1.0 release`.
- [ ] **Step 9: PR and merge.**
  - After the user approves: push `release-0.1.0` and open the PR.
  - CI must be green.
  - After the user approves: merge it with a merge commit.
- [ ] **Step 10: STOP AND ASK.**
  - Tagging and releasing publishes to PyPI, and that cannot be undone. Ask
    the user, in one message, for explicit approval to:
    1. tag the merged `master` commit `v0.1.0`;
    2. push only that tag;
    3. create the GitHub release, which triggers the upload.
  - Wait for a yes that names this release. An earlier approval does not
    count.
- [ ] **Step 11: Tag and release**, only after Step 10's approval:
  ```bash
  git switch master && git pull --ff-only
  git log -1 --format=%s   # the merge of release-0.1.0
  git tag v0.1.0
  git push origin v0.1.0
  awk '/^## \[0.1.0\]/{on=1;next} /^## \[/{on=0} on' CHANGELOG.md > /tmp/claude-1000/notes-0.1.0.md
  gh release create v0.1.0 --title "0.1.0" --notes-file /tmp/claude-1000/notes-0.1.0.md
  gh run list --workflow release.yaml --limit 3 --json databaseId,headBranch,status
  ```
  - Pick the run whose `headBranch` is `v0.1.0`, since the list can still
    show 0.0.1's run for a few seconds. Then run
    `gh run watch <its databaseId> --exit-status`.
  - If the environment asks for a reviewer, the user approves it in the run.
  - If the run fails before the upload, follow the playbook's Common
    mistakes: fix `master`, then delete the release and tag
    (`gh release delete v0.1.0 --cleanup-tag`) and re-tag, only with the
    user's approval.
  - If the run fails after the upload, report it. The version is spent, and
    the fix is 0.1.1.
- [ ] **Step 12: Verify** (playbook Verification). Both print `0.1.0`:
  ```bash
  curl -s https://pypi.org/pypi/biotapy/json | python3 -c "import json,sys; print(json.load(sys.stdin)['info']['version'])"
  uv run --no-project --with biotapy==0.1.0 python -c "import biotapy as bt; print(bt.__version__, bt.pl.__all__)"
  ```
  The second also prints `['bar', 'heatmap', 'ordination', 'richness', 'scree']`.
- [ ] **Step 13: Proposals, not submitted.** Write both into the scratchpad,
  show them to the user, and submit neither.
  - **conda-forge.**
    - A `recipe/recipe.yaml` draft: `noarch: python`, source from the PyPI
      sdist and its sha256, `pip install . --no-deps` over hatchling, and
      the `[project] dependencies` as run requirements, using conda-forge's
      `matplotlib-base`.
    - The blocker, first: `treedata` is on neither conda-forge nor bioconda
      (the research checked api.anaconda.org on 2026-09-27; re-check). A
      biotapy recipe cannot be merged until a `treedata` recipe is. The
      proposal is therefore two staged-recipes PRs, `treedata` first; ask the
      user whether to ask treedata's maintainers (YosefLab) about it.
    - Every other runtime dependency is on conda-forge.
  - **scverse ecosystem.** A `packages/biotapy/meta.yaml` for
    `scverse/ecosystem-packages`, in the format of its existing entries
    (`name`, `description`, `project_home`, `documentation_home`,
    `tutorials_home`, `install: {pypi: biotapy}`, `tags`,
    `primary_category`, `language: Python`, `license: BSD-3-Clause`,
    `contact: [pedrocr83]`, `test_command`).
    - Read `scripts/src/ecosystem_scripts/schema.json` in that repository
      for the controlled vocabulary first. On 2026-09-27 no tag or category
      named microbiome existed, and the README says to add one to the enum
      in the same PR.
    - Propose `microbiome` as a new tag, alongside the existing
      `preprocessing`, `visualization`, `dimensionality reduction` and
      `interoperability`.
    - Answer the README's mandatory checklist item by item, from facts:
      BSD-3 license, 0.1.0 on PyPI, tests and CI, the API docs, and
      AnnData/TreeData throughout.
- [ ] **Step 14: Close Phase 1.** On a branch `close-phase-1` from `master`:
  - `roadmap/phase-1-core.md`:
    - tick the last exit-gate item;
    - set `phase_state: done`;
    - tick this task.
  - `roadmap/phase-2-function.md`: set `phase_state: in-progress`.
  - `roadmap/index.md`:
    - move Phase 1 under "Phases" with `**phase_state: done** (0.1.0 on
      PyPI, <release date>)`;
    - put Phase 2 under "Active phase" with `**phase_state: in-progress**`.
  - Add a log line.
  - Run `uv run --group test pytest tests/test_knowledge_bundle.py -q`.
  - Commit `docs(knowledge): close Phase 1 after the 0.1.0 release`.
  - Push, open the PR and merge, each after the user approves.
  - Then `bash scripts/knowledge_stale.sh` (against `origin/master`) prints
    0 stale.
  - Phase 2 starts only when the user asks. Its tasks get expanded first
    (R1.2a).

---

# Exit gate
- [ ] `docs/tutorials/phyloseq_analysis.md` (a MyST notebook; slice 1D ruling 13) executes in CI and reproduces the
  vignette's in-scope sections on GlobalPatterns, enterotype and esophagus:
  bar plots, richness, UniFrac + PCoA with scree, NMDS, heatmap, `distance()`
  examples. Out-of-scope sections are listed in the notebook as "not in 0.1":
  `plot_tree`, `plot_net`, CCA, DPCoA.
- [ ] All `golden` tests pass per [r-golden-parity](/contracts/r-golden-parity.md).
- [ ] The generated Coming-from-R table maps all 31 functions below.
- [x] asv baselines recorded.
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
