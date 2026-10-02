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

# Quick tour

This page runs at build time, against `bt.datasets.toy()` only - no network,
no external files. It walks through biotapy's data model, its two
preprocessing functions, and reading and writing your own data.

```{code-cell} ipython3
import biotapy as bt

tdata = bt.datasets.toy()
tdata
```

## The data model

Every biotapy object is a `TreeData`: an [AnnData](https://anndata.readthedocs.io/)
that also carries a phylogeny. Samples are always rows:

```{code-cell} ipython3
tdata.shape  # (6 samples, 8 features)
```

`X` stays a sparse matrix, never a dense array:

```{code-cell} ipython3
type(tdata.X)
```

Taxonomy lives in `var`, one column per rank from `kingdom` down; `f8` has no
`genus` assignment, which biotapy represents as `NaN` rather than a string:

```{code-cell} ipython3
tdata.var
```

The phylogeny lives in `vart["phylo"]`, keyed by the same feature ids as
`var_names`:

```{code-cell} ipython3
tree = tdata.vart["phylo"]
tree.number_of_nodes(), tree.number_of_edges()
```

`uns["biotapy"]` is biotapy's own bookkeeping: `x_kind` records what `X`
currently holds, and `provenance` lists every biotapy function that produced
this object, in order:

```{code-cell} ipython3
tdata.uns["biotapy"]
```

## Relative abundance

`bt.pp.relative` adds per-sample relative abundance as a new layer and
leaves `X` untouched - every biotapy `pp` function is pure by default:

```{code-cell} ipython3
import numpy as np

rel = bt.pp.relative(tdata)
rel.layers["relative"][0].sum()  # each sample's relative abundances sum to 1
```

```{code-cell} ipython3
np.array_equal(tdata.X.toarray(), rel.X.toarray())  # X is unchanged
```

## Aggregating to a rank

`bt.pp.tax_glom` merges features that share a lineage down to a chosen rank,
each merged feature represented by its most abundant member:

```{code-cell} ipython3
glom = bt.pp.tax_glom(tdata, "phylum")
glom.n_vars, list(glom.var_names)
```

## Reading and writing BIOM

`bt.io.write_biom` and `bt.io.read_biom` round-trip `X`, taxonomy and sample
metadata through a BIOM 2.1 (HDF5) table:

```{code-cell} ipython3
import tempfile
from pathlib import Path

with tempfile.TemporaryDirectory() as tmp:
    biom_path = Path(tmp) / "toy.biom"
    bt.io.write_biom(tdata, biom_path)
    back = bt.io.read_biom(biom_path)

back.shape, np.array_equal(tdata.X.toarray(), back.X.toarray())
```

## Saving your work

A `TreeData` round-trips through treedata's own format, tree included, with
no biotapy-specific step:

```{code-cell} ipython3
import treedata as td

with tempfile.TemporaryDirectory() as tmp:
    h5_path = Path(tmp) / "toy.h5td"
    tdata.write_h5td(h5_path)
    restored = td.read_h5td(h5_path)

type(restored), restored.shape
```
