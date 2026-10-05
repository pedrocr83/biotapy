# Differential abundance

Differential abundance (DA) asks which features are more abundant in one group of samples than in
another. Every `bt.da` method takes the same arguments and returns the same table, so their answers
can be put side by side.

## One question, one table

A method takes the `obs` column whose effect you want (`group`), the `covariates` to adjust for,
and the `reference` level that the other level is compared with:

```python
import biotapy as bt

tdata = bt.datasets.toy()
table = bt.da.linda(tdata, "group")  # "B vs A": A, the first category, is the reference
table = bt.da.linda(tdata, "group", reference="B")  # "A vs B": every effect changes sign
```

`group` is a categorical, string or bool column with two levels, or a numeric column, whose effect
is then a slope. Covariates are numeric, or categorical with one indicator per level against
their first. There are no formula strings: R and patsy take the alphabetically first level as the
reference and silently drop samples with a missing value. biotapy takes the reference you give, or
the first category, says which way round the effect is in the `contrast` column, and raises on a
missing value, on a constant column, on a group with more than two levels and on collinear covariates.

The table has one row per feature, in `var_names` order:

| Column | Meaning |
|---|---|
| `effect` | log2 fold change of the other level over the reference, or the slope of a numeric group |
| `se` | its standard error |
| `pvalue` | the method's p-value |
| `qvalue` | Benjamini-Hochberg adjusted p-value, over the features the method tested |
| `direction` | sign of `effect`: -1, 0 or 1 |
| `method` | the method's name |
| `contrast` | `"<level> vs <reference>"`, or the name of a numeric group |

A feature a method cannot test keeps its row, with NaN values and `direction` 0. Methods never
filter features: filter once, with `bt.pp.filter_features`, before running any method, so every
method tests the same features. Methods need raw counts in `X` and a read in every sample
(`bt.pp.filter_samples(tdata, min_depth=1)` drops empty ones).

## LinDA

`bt.da.linda` fits one linear model per feature to the log2 centred log-ratios of the counts (with
0.5 added to every count when the table has a zero). Compositionality biases every coefficient by
the same amount; LinDA estimates that bias as the mode of all features' coefficients and subtracts
it, so most features end up near zero and the ones that changed stand out:

```python
tdata.obs["age"] = [30, 41, 52, 38, 45, 60]
table = bt.da.linda(tdata, "group", covariates=["age"])
table[table["qvalue"] < 0.05]
```

It equals `MicrobiomeStat::linda(..., is.winsor = FALSE)` in R with fixed effects. Numeric columns
are scaled to unit variance first, as LinDA does, so a numeric group's effect is per standard
deviation. Not available: winsorisation (MicrobiomeStat's default), random effects such as
`(1 | subject)`, and LinDA's own prevalence filters.

## ANCOM-BC2

`bt.da.ancombc2` runs scikit-bio's ANCOM-BC2 (Lin and Peddada 2024): it estimates each sample's
sampling fraction, corrects the log counts for it and the coefficients for their shared bias, and
fits one linear model per feature. Zeros are treated as missing rather than given a pseudocount,
so a feature with no read in one of the groups cannot be fitted: its row is NaN and it is left out
of the Benjamini-Hochberg correction (R's `ANCOMBC::ancombc2` reports it with p = 1 and counts it).

```python
table = bt.da.ancombc2(tdata, "group")
```

ANCOM-BC2 reports natural logs; biotapy divides `effect` and `se` by ln 2, so they are log2 like
every other method's; a numeric `group`'s effect is per unit, where `da.linda`'s is per standard
deviation. The settings are R's `ancombc2(..., p_adj_method = "BH", prv_cut = 0,
pseudo_sens = FALSE)`: biotapy does not run R's pseudocount sensitivity analysis (its
`passed_ss` flag) or its 10% prevalence filter. On the GlobalPatterns genera that the golden tests
use, biotapy's effects are within 0.012 log2 of R's and the significant genera are the same; the
small difference comes from the bias estimate, whose iterations stop at R's cap of 100 before they
have converged on that data, in R as in scikit-bio.

## Where methods agree

`bt.da.consensus` puts the tables of several methods side by side and counts, for each feature,
the methods that call it significant (`qvalue < alpha`, strictly) and whether they agree on its
direction:

```python
results = [bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")]
table = bt.da.consensus(results)  # alpha=0.05, min_methods=2
table[table["consensus"]]
```

A feature is a consensus hit when at least `min_methods` methods call it and all of them give it
the same sign. Methods that call it in opposite directions mark a `conflict`, which is never a
consensus. A method that could not test a feature does not count against it: `n_tested` says how
many methods tested each feature. The tables must come from different methods and compare the same
`contrast`, so run every method with the same `group` and `reference`.

Methods disagree a lot on real data, and the literature is split on what to do about it: Nearing
et al. (2022) recommend a consensus of several methods, Pelto et al. (2025) one simple method. Either
way, choose the methods before you look at their results; trying methods until one finds what you
hoped for is selective reporting, and a consensus table does not protect against it. Methods that
share a model, such as ANCOM-BC and ANCOM-BC2, also agree more often for that reason alone.

## Plotting the consensus

`bt.pl.consensus` draws the consensus table as a dot matrix: one row per feature that at least one
method calls, one column per method. A filled dot is a call, red for a positive effect and blue for
a negative one; a hollow dot is a feature the method tested without calling it; no dot, a feature
the method could not test. Consensus features have bold labels. Rows are sorted by the number of
methods calling them, then by their mean absolute effect, and `top` keeps the first 30:

```python
ax = bt.pl.consensus(table, top=20)
```

Feature ids are often accession numbers; to label rows with a rank, rename the table's index first,
for example `table.rename(index=tdata.var["genus"])`.
