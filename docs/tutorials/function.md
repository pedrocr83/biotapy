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

# Function in an IBD cohort

This tutorial reads the pathway abundance of a real cohort, the Integrative Human Microbiome
Project's inflammatory bowel disease study (HMP2), renormalises it, asks which species carry a
butyrate pathway, and compares Crohn's disease (CD), ulcerative colitis (UC) and non-IBD
controls. `bt.datasets.hmp2()` holds each participant's first stool metagenome: HUMAnN 3
pathways in the `"function"` and `"function_by_taxon"` modalities and MetaPhlAn 3 species in
`"taxa"`.

:::{note}
The HMP2 tables are downloaded from the [IBDMDB](https://ibdmdb.org/) on first use (23 MB) and
cached. The IBDMDB states no licence for them; biotapy ships none of them. Cite the study when
you use them: Lloyd-Price J et al. (2019) Multi-omics of the gut microbial ecosystem in
inflammatory bowel diseases. *Nature* 569:655-662.
:::

Two parts of `bt.fn` are not shown here, and the [function guide](../guide/function.md) shows
both on small tables:

- **[Grouping along a hierarchy](../guide/function.md#aggregating-along-a-hierarchy)** (`bt.fn.func_glom`): HMP2's pathways are MetaCyc pathways,
  whose classes biotapy does not ship or download (their licence does not allow it), and HMP2
  publishes its enzyme (EC) table only per sample or as a 113 MB merged file.
- **[Functional redundancy](../guide/function.md#functional-redundancy)** (`bt.fn.functional_redundancy`): it needs each species' gene copy
  numbers, which HUMAnN does not write.

```{code-cell} ipython3
import matplotlib.pyplot as plt

import biotapy as bt
```

## The cohort

```{code-cell} ipython3
mdata = bt.datasets.hmp2()
mdata
```

```{code-cell} ipython3
mdata.obs["diagnosis"].value_counts()
```

## Renormalising

HMP2's pathway tables are in copies per million (CPM), and most of each sample (a median of 96%
here) sits in `UNMAPPED` and `UNINTEGRATED`, the reads HUMAnN could not place in a pathway.
`bt.fn.renorm(mdata, "relab", special=False)` drops those rows and divides every pathway, and
every per-species row, by the total of the pathways left, as `humann_renorm_table --special n`
does. The `"taxa"` modality is kept as it is.

```{code-cell} ipython3
relab = bt.fn.renorm(mdata, "relab", special=False)
community = relab["function"]
community.X.sum(axis=1)[:5]
```

## Who carries butyrate production

`PWY-5676`, acetyl-CoA fermentation to butanoate II, is a butyrate pathway HUMAnN finds in 113
of the 130 samples. `bt.fn.contributions` splits it by species, here keeping the five with
the largest total and summing the rest into `other`. HMP2's table has per-species rows
(`unclassified` included) for 88 of those 113 samples, so the other 25 draw an empty bar
below. The values are shares of each sample's mapped pathway total (after the `renorm` above),
not of `PWY-5676` itself, and a pathway's strata need not add up to its community value (here
they add up to a median of a third of it; the [function guide](../guide/function.md#contributions)
says why). In the mean per diagnosis, the 25 samples with the pathway but no strata count as
zero for every species:

```{code-cell} ipython3
by_taxon = relab["function_by_taxon"]
butyrate = bt.fn.contributions(by_taxon, "PWY-5676", top=5)
butyrate.groupby(by_taxon.obs["diagnosis"], observed=True).mean()
```

`bt.pl.contributions` draws the same table, one bar per sample. Sorting the samples by
diagnosis first puts each group's bars side by side; 130 sample names are too many to read, so
the x axis names the groups instead:

```{code-cell} ipython3
order = by_taxon.obs.sort_values("diagnosis", kind="stable").index
fig, ax = plt.subplots(figsize=(14, 4))
bt.pl.contributions(by_taxon[order].copy(), "PWY-5676", top=5, ax=ax)
ends = by_taxon.obs["diagnosis"].value_counts(sort=False).cumsum()
for end in ends.iloc[:-1]:
    ax.axvline(end - 0.5, color="black", linewidth=0.8)
ax.set_xticks((ends + ends.shift(fill_value=0)) / 2 - 0.5, ends.index, rotation=0)
ax.legend(title="taxon", loc="upper left", bbox_to_anchor=(1, 1));
```

## Comparing diagnoses

Bray-Curtis distances between the samples' pathway profiles, their principal coordinates, and
PERMANOVA by diagnosis. Each participant contributes one sample, so the test counts each
person once.

```{code-cell} ipython3
bt.tl.beta(community, inplace=True)
bt.tl.pcoa(community, inplace=True)
bt.pl.ordination(community, color="diagnosis");
```

```{code-cell} ipython3
bt.tl.permanova(community, "diagnosis", seed=0)
```

The species side of the same samples: Shannon diversity of the MetaPhlAn profiles per
diagnosis.

```{code-cell} ipython3
taxa = relab["taxa"]
bt.tl.alpha(taxa, metrics=["shannon"], inplace=True)
bt.pl.richness(taxa, "shannon", x="diagnosis", color="diagnosis");
```
