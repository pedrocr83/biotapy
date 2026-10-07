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

# Differential abundance in GlobalPatterns

This tutorial asks which genera differ between human-associated samples and the rest of
GlobalPatterns, runs two differential abundance methods on the same question, and reports where
they agree. It then shows the same report with the two methods that run in R.

:::{note}
The contrast is a demonstration of the workflow, not a biological finding. GlobalPatterns has 9
human-associated samples (feces, skin and tongue) and 17 others: 14 from soil, sediment,
freshwater and the ocean, and 3 mock communities. The data are downloaded from phyloseq's
repository on first use and cached; they are licensed to phyloseq's authors under AGPL-3.
:::

```{code-cell} ipython3
import numpy as np
import pandas as pd

import biotapy as bt
```

## The data

Counts are merged to genus, and genera present in fewer than 20% of the samples are dropped. This
is the one filter: the methods never filter on their own, so each of them tests the same 636
genera.

```{code-cell} ipython3
global_patterns = bt.datasets.global_patterns()
genus = bt.pp.tax_glom(global_patterns, "genus")
tdata = bt.pp.filter_features(genus, min_prevalence=0.2)
human = tdata.obs["SampleType"].isin(["Feces", "Skin", "Tongue"])
tdata.obs["host"] = pd.Categorical(np.where(human, "human", "other"))
tdata.obs["host"].value_counts()
```

## Two methods, chosen first

ANCOM-BC2 and LinDA both run in Python. They are chosen before looking at any result, and run with
the same `group` and `reference`, so their tables compare the same thing (`"human vs other"`): a
positive `effect` means more abundant in human-associated samples.

```{code-cell} ipython3
results = [
    bt.da.ancombc2(tdata, "host", reference="other"),
    bt.da.linda(tdata, "host", reference="other"),
]
{table["method"].iloc[0]: int((table["qvalue"] < 0.05).sum()) for table in results}
```

ANCOM-BC2 leaves zeros out, so it cannot fit the 36 genera with no read in one of the two groups:
their rows are NaN, and they count as not tested, not as "not significant".

```{code-cell} ipython3
int(results[0]["pvalue"].isna().sum())
```

## Where they agree

`bt.da.consensus` counts, per genus, the methods that call it at `qvalue < 0.05` and whether they
agree on its direction. A consensus genus is called by at least `min_methods` methods (2 by
default), all with the same sign.

```{code-cell} ipython3
table = bt.da.consensus(results)
table["n_significant"].value_counts().sort_index()
```

```{code-cell} ipython3
table[["consensus", "conflict"]].sum()
```

Which method calls what, genus by genus:

```{code-cell} ipython3
pd.crosstab(table["significant_ancombc2"], table["significant_linda"])
```

The dot plot shows the genera called by most methods first, then by the size of their effect.
Feature ids in GlobalPatterns are OTU numbers, so the rows are labelled with the genus name; a
name can appear twice, because GlobalPatterns files some genera under more than one family and
`tax_glom` keeps the lineages apart.

```{code-cell} ipython3
genera = table.rename(index=tdata.var["genus"])
bt.pl.consensus(genera, top=30);
```

The genera only one method calls are where the two disagree. A hollow dot is a genus the other
method tested without calling it:

```{code-cell} ipython3
bt.pl.consensus(genera[genera["n_significant"] == 1], top=30);
```

## Adding the methods that run in R

ALDEx2 and MaAsLin 3 run in R ({ref}`Methods that run in R <da-methods-in-r>`), which this
documentation is built without, so the code below is shown but not run here:

```python
results = [
    bt.da.ancombc2(tdata, "host", reference="other"),
    bt.da.linda(tdata, "host", reference="other"),
    bt.da.aldex2(tdata, "host", reference="other", seed=0),
    bt.da.maaslin3(tdata, "host", reference="other", seed=0),
]
table = bt.da.consensus(results)
```

biotapy's continuous integration runs exactly this on every pull request, in R 4.5.3 with
Bioconductor 3.22's ALDEx2 and maaslin3, and checks that it gives these numbers:

| Method | Genera called |
|---|---|
| `ancombc2` | 208 |
| `linda` | 118 |
| `aldex2` | 13 |
| `maaslin3` | 52 |

| Methods calling a genus | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| Genera | 413 | 114 | 62 | 35 | 12 |

With four methods, 109 genera are a consensus at `min_methods=2`, none is a conflict, and the 12
genera all four methods call are all more abundant in human-associated samples. The methods differ
a lot in how many genera they call; the consensus table says which calls do not depend on the
choice of method.

## Reading the result

- A genus called by one method only is a result of that method, not of the data alone.
- `min_methods` is part of the analysis: decide it, like the methods, before looking.
- The {doc}`method pages </methods/index>` say what each `effect` measures and how each method's
  numbers compare with its R package.
