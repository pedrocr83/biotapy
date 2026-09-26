# Data model

Every biotapy dataset is one object: a [treedata](https://treedata.readthedocs.io/) `TreeData`
(an [AnnData](https://anndata.readthedocs.io/) that also carries a phylogeny). This page covers
the one layout rule every function relies on, and what each slot on that object holds.

## Samples are always rows

`X` is samples (rows) by features (columns) - never the other way round. If you load data from
a tool that stores features as rows (most do), the biotapy reader transposes it for you, once,
so everything downstream can assume the same orientation.

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.shape  # (6 samples, 8 features)
```

## Slots

| Slot | Holds | Example keys |
|---|---|---|
| `X` | Counts (or another abundance) as a sparse matrix | - |
| `layers` | Transforms of `X` with the same shape | `relative` |
| `obs` | Sample metadata, and per-sample results | `group`, `alpha_shannon` |
| `var` | Taxonomy, one column per rank, and sequences | `kingdom` .. `species`, `sequence` |
| `vart` | The phylogeny, leaves named after your features | `phylo` |
| `obsm` | Ordinations and embeddings | `X_pcoa` |
| `obsp` | Sample-by-sample distance matrices | `braycurtis`, `weighted_unifrac` |
| `uns["biotapy"]` | biotapy's own bookkeeping | `x_kind`, `provenance` |

`uns["biotapy"]["x_kind"]` records what `X` currently holds (`counts`, `relative`, `rpk`, `cpm`
or `abundance`); a missing key means `counts`. `uns["biotapy"]["provenance"]` lists every
biotapy function that has touched the object, in order, so you can always tell how it reached
its current state.

## What survives a filter or an aggregation

Functions that change which features exist (dropping rare taxa, aggregating to a rank) drop
`layers`, `obsm` and `obsp`, because a transform or a distance computed on the old features
would silently misdescribe the new ones. Functions that only add a layer, or that subset
samples, leave everything else in place.

See the [data-model-slots contract](https://github.com/pedrocr83/biotapy/blob/master/.knowledge/contracts/data-model-slots.md)
for the exact rules every biotapy function follows.
