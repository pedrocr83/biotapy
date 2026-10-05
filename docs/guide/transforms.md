# Transforms

Transforms add a layer (or, for PhILR, an `obsm` table) to a copy of your data; they never touch `X` or drop
anything else, so you can always get back to what you started from.

## Relative abundance

`bt.pp.relative` scales every sample (row) so its features sum to one, and stores the result in
`layers["relative"]`:

```python
import biotapy as bt

tdata = bt.datasets.toy()
out = bt.pp.relative(tdata)
out.layers["relative"][0].sum()  # 1.0
```

`X` still holds the original counts; `out.uns["biotapy"]["provenance"]` records that
`pp.relative` ran.

### All-zero samples

A sample with no reads has nothing to divide by. Where phyloseq's
`transform_sample_counts` divides by zero and returns `NaN`, biotapy leaves an
all-zero sample as all zeros, so downstream steps do not have to special-case `NaN`.

## Centred log-ratio (CLR)

Microbiome counts are compositional: a sequencing run fixes the total, so only ratios between
features carry information. `bt.pp.clr` takes the logarithm of each value relative to the
geometric mean of its sample, and stores it in `layers["clr"]`:

```python
out = bt.pp.clr(tdata)  # pseudocount=0.5
out.layers["clr"][0].sum()  # 0.0, up to rounding
```

The logarithm of zero is undefined, so a pseudocount is added to every value first, zeros or
not, as `mia::transformAssay(method = "clr", pseudocount = 0.5)` and
`vegan::decostand(x, "clr", pseudocount = 0.5)` do. The default 0.5 suits counts. On relative
abundances or CPM it would swamp the rarest features, so biotapy warns when the pseudocount is
larger than the smallest non-zero value; pass one on the data's scale, for example half that value.

An all-zero sample gives an all-zero CLR row.

`layers["clr"]` is dense: CLR has no zeros, so it takes 8 bytes per sample and feature, and the call peaks at three to five such arrays, the most when `X` is mostly non-zero.

## PhILR

`bt.pp.philr` turns a sample into balances along its phylogeny: one per internal node with two
children, the log-ratio of the geometric means of the taxa under the node's two children, scaled so
that the balances are orthonormal coordinates. They go to `obsm["X_philr"]`, a table with one
column per balance:

```python
tdata = bt.datasets.toy()[:, :6].copy()
out = bt.pp.philr(tdata)  # pseudocount=0.5, as in pp.clr
out.obsm["X_philr"]  # columns root, n1, n4, n2, n5
```

A balance is positive when the node's first child is more abundant than its second, as in
`philr::philr` with its default uniform weights. Distances between samples' balances equal
their Aitchison distances, the Euclidean distances between their CLR vectors.

PhILR needs a rooted binary tree. Filtering features leaves nodes with one child; those define
no balance and are skipped. A node with three or more children raises: that includes the root
of an unrooted tree (`bt.datasets.toy()` has one, which is why the example keeps f1-f6). Root
the tree and resolve its polytomies before reading it, for example with `ape::multi2di` in R, or
in Python with scikit-bio, which resolves them arbitrarily, as `multi2di` does:

```python
from skbio import TreeNode

tree = TreeNode.read("tree.nwk")  # root it first (outgroup or midpoint) if it is unrooted: your choice
tree.bifurcate()  # in place; every node ends with two children
tree.write("binary.nwk")
# then read your table again with tree="binary.nwk"
```

Balance names are the node names in `vart["phylo"]`, which differ from the ones `philr` in R
makes with `makeNodeLabel`; to compare balances across tools, match them by their numerator and
denominator taxa, not by name.

Filtering features with `bt.pp.filter_features` drops `layers["clr"]` and `obsm["X_philr"]`;
plain `adata[:, ...]` slicing keeps the old values, so transform after subsetting.
