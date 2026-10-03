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
box plot per group (`+ geom_boxplot()`). biotapy stores Gini-Simpson as `simpson` and has no
InvSimpson metric; derive phyloseq's from it as `1 / (1 - gp.obs["alpha_simpson"])`.

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
