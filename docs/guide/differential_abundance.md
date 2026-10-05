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
missing value, on a group with more than two levels and on collinear covariates.

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
