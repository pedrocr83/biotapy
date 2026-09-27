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
