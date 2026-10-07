# MaAsLin 3

{func}`bt.da.maaslin3 <biotapy.da.maaslin3>` runs the abundance model of `maaslin3::maaslin3` in R
through rpy2: one linear model of each feature's log2 relative abundance, on the samples where the
feature is present, tested against the median coefficient. It needs R, the R package maaslin3 and
biotapy's `r` extra ({ref}`Methods that run in R <da-methods-in-r>`).

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.obs["age"] = [30, 41, 52, 38, 45, 60]
table = bt.da.maaslin3(tdata, "group", covariates=["age"], seed=0)
```

## Model

Counts are divided by their sample's total (total sum scaling). For feature $j$, on the samples
where it is present ($y_{ij} > 0$):

$$
\log_2 \frac{y_{ij}}{\sum_k y_{ik}} = x_i^\top \beta_j + \varepsilon_{ij}
$$

with numeric columns of the design $x_i$ standardised. Compositionality shifts every coefficient
of a term by about the same amount, so MaAsLin 3 tests each coefficient against the median
coefficient $\tilde\beta$ of that term, taken over the features without a fit error whose own
p-value is below 0.95; the test draws 10,000 normal samples, which is why `seed` matters.
biotapy reports

$$
\mathrm{effect}_j = \hat\beta_j - \tilde\beta
$$

(`subtract_median = TRUE`), so the sign of `effect` is the side of the median the p-value is
about. `qvalue` is the Benjamini-Hochberg correction of the group's p-values. A feature whose fit
reports an error is NaN in the table (not tested), as MaAsLin 3 leaves it out of its own
correction.

## Units

`effect` is a log2 fold change of the relative abundance where the feature is present, minus the
median. Numeric columns are standardised, so a numeric `group`'s effect is per standard
deviation, as in `da.linda`.

## Compared with R

biotapy calls `maaslin3::maaslin3` (maaslin3 1.2.0) with:

| `maaslin3` argument | R default | biotapy |
|---|---|---|
| `formula` | `NULL` | `~ x0 + x1 + ...`, the `group` and `covariates` columns under plain names, so a column such as `"body site"` needs no quoting |
| `reference` | `NULL` | not passed: categorical columns arrive as factors with string levels, `reference` first |
| `random_effects`, `group_effects`, `ordered_effects`, `strata_effects` | `NULL` | not offered |
| `min_abundance`, `min_prevalence` | `0`, `0` | the same: no filter |
| `normalization`, `transform`, `standardize` | `"TSS"`, `"LOG"`, `TRUE` | the same |
| `median_comparison_abundance` | `TRUE` | the same |
| `subtract_median` | `FALSE` | `TRUE`: `effect` is what the p-value tests |
| `evaluate_only` | `NULL` (both models) | `"abundance"`: the prevalence model's log-odds cannot share a column with fold changes |
| `warn_prevalence` | `TRUE` | `FALSE`, which `evaluate_only` requires |
| `correction` | `"BH"` | BH over the group's p-values only; MaAsLin 3's `qval_individual` corrects them together with every covariate's |
| `plot_summary_plot`, `plot_associations` | `TRUE` | `FALSE`; the `output` folder is temporary and deleted |
| `cores`, `verbosity` | `1`, `"FINEST"` | `1`, `"ERROR"` |

With one numeric covariate on the GlobalPatterns genera below, `qval_individual` calls 20 genera
where biotapy's `qvalue` calls 45; without covariates the two agree (52).

## Agreement with R

MaAsLin 3's median test simulates, but `seed` becomes one integer for R's `set.seed`, so biotapy
and R draw the same numbers. The golden test runs `set.seed(...); maaslin3::maaslin3(...)` with
that integer on the 636 GlobalPatterns genera in at least 20% of samples, human hosts against the
rest, with and without a sequencing-depth covariate: the same 36 genera fail to fit, and `effect`,
`se` and `pvalue` agree with R's `coef`, `stderr` and `pval_individual` to a relative 7.2e-13 or
better (checked at 1e-7). With other seeds only the p-values move, by up to about 1.6e-3.

## Choosing the reference

Swapping `reference` negates `effect` exactly, but the median test's simulation draws around the
coefficients rather than their negatives, so the same `seed` moves `pvalue` by up to about 1.2e-3 and
`qvalue` by up to about 3.3e-3 on the GlobalPatterns genera (the calls are the same).

## Reference

Nickols WA, Kuntz T, Shen J, Maharjan S, Mallick H, Franzosa EA, Thompson KN, Nearing JT,
Huttenhower C (2026) MaAsLin 3: refining and extending generalized multivariable linear models
for meta-omic association discovery. *Nature Methods* 23:554-564.
[doi:10.1038/s41592-025-02923-9](https://doi.org/10.1038/s41592-025-02923-9)
