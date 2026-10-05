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

Filtering changes the feature set, so these are dropped: every entry in `layers`, `obsm`,
`obsp`, `varm` and `varp`, every `uns` key other than `uns["biotapy"]`, and the ordination
summaries (`pcoa`, `nmds`) inside it. `obs` is kept whole, as mia keeps `colData`, so the
`alpha_*` columns that `bt.tl.alpha(..., inplace=True)` writes survive `filter_features`,
`rarefy` and `tax_glom` while still describing the old features: recompute them. A TreeData
keeps the subtree of the kept features, with branch lengths unchanged.

## Samples: `filter_samples`

`bt.pp.filter_samples(tdata, min_depth)` keeps samples with at least `min_depth` reads (again
inclusive). It only removes rows, so every slot is kept and subset, and distances in `obsp` stay
valid for the samples that remain. Features that become all-zero are kept; follow with
`filter_features` to drop them.

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
- `X` must hold raw counts: `uns["biotapy"]["x_kind"] == "counts"` and non-negative whole numbers.
- phyloseq samples with replacement by default (`replace = TRUE`). biotapy always samples
  without, like `replace = FALSE`.
- The same `seed` gives the same result. R and NumPy use different random generators, so the
  counts never match phyloseq's draw for draw.
