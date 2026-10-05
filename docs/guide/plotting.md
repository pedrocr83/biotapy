# Plotting

`bt.pl` draws what `tl`, `pp` and `fn` give; it computes no diversity and no ordination. Each
function takes `ax=` to draw on existing axes, or makes a new figure, and returns the
`matplotlib.axes.Axes`, so you finish the plot with matplotlib. When a stored result is
missing, the error names the call to run:

| Function | Reads | Written by |
|---|---|---|
| `pl.bar` | `X` or `layers[layer]`, and `var`/`obs` columns | `pp.relative` for `layer="relative"` |
| `pl.heatmap` | `X` or `layers[layer]` | `pp.relative` for `layer="relative"` |
| `pl.contributions` | a stratified function table's `X` and `var["function"]`, `var["taxon"]` | `io.read_humann`, `io.read_picrust2` (the `"function_by_taxon"` modality) |
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

## Contributions: `pl.contributions`

`pl.contributions(mdata["function_by_taxon"], "PWY-5100")` draws one bar per sample for one
function, one segment per taxon carrying it: the table `bt.fn.contributions` returns. `top=8`
(the default) keeps the eight taxa with the largest total and draws the rest as one grey
segment, `other`; `top=None` draws every taxon. Bars are the strata as stored, so a pathway's
bar need not reach its community value; plot the stratified modality of
`bt.fn.renorm(mdata, "relab")` to read shares of each sample's community total.

```python
mdata = bt.datasets.toy_humann()
bt.pl.contributions(mdata["function_by_taxon"], "2.7.1.2", top=5)
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
