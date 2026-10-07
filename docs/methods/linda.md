# LinDA

{func}`bt.da.linda <biotapy.da.linda>` fits one linear model per feature to the log2 centred
log-ratios of the counts and removes the bias that compositionality adds to every coefficient.
It runs in Python; no R is needed.

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.obs["age"] = [30, 41, 52, 38, 45, 60]
table = bt.da.linda(tdata, "group", covariates=["age"])
table[table["qvalue"] < 0.05]
```

## Model

For counts $y_{ij}$ of feature $j$ in sample $i$ ($n$ samples, $m$ features), 0.5 is added to
every count when the table holds a zero. Each sample's log2 counts are centred:

$$
w_{ij} = \log_2 y_{ij} - \frac{1}{m} \sum_{k=1}^{m} \log_2 y_{ik}
$$

Each feature gets one least-squares fit $w_{\cdot j} = X \beta_j + \varepsilon_j$ on the design
$X$ (an intercept, `group` and the `covariates`), with $p$ columns. The group coefficient
$\tilde\beta_j$ estimates the log2 fold change plus a bias $b$ that is the same for every feature,
because the centring divides by the geometric mean of a composition. LinDA estimates the bias as
the mode of the coefficients and subtracts it:

$$
\hat b = \frac{1}{\sqrt n}\, \operatorname{mode}\{\sqrt n\, \tilde\beta_1, \dots, \sqrt n\, \tilde\beta_m\},
\qquad \mathrm{effect}_j = \tilde\beta_j - \hat b
$$

The mode is found by a Gaussian mean shift started at the mean of the shortest half of the values
(`shorth`), with `bw.nrd0`'s bandwidth, as `modeest::mlv(method = "meanshift")` does. The
p-value is a two-sided t-test of $\mathrm{effect}_j / \mathrm{se}_j$ on $n - p$ degrees of
freedom, where $\mathrm{se}_j$ is the least-squares standard error of $\tilde\beta_j$, and
`qvalue` is the Benjamini-Hochberg correction over all features. Most features end up near zero;
the ones that changed stand out.

## Units

`effect` is a log2 fold change of the other level over `reference`. Numeric columns, a numeric
`group` included, are scaled to unit variance first, as LinDA does, so a numeric group's effect
is per standard deviation.

## Compared with R

biotapy equals `MicrobiomeStat::linda(t(counts), meta, "~group + covariates", feature.dat.type =
"count", is.winsor = FALSE)` (MicrobiomeStat 1.4) with fixed effects:

| `linda` argument | R default | biotapy |
|---|---|---|
| `formula` | required | `group` and `covariates`, fixed effects only: no `(1 \| subject)` |
| `feature.dat.type` | `"count"` | counts only: `X` must hold raw counts |
| `prev.filter`, `mean.abund.filter`, `max.abund.filter` | `0` | no filter: run `bt.pp.filter_features` once, before any method |
| `is.winsor`, `outlier.pct` | `TRUE`, `0.03` | no winsorisation |
| `adaptive`, `zero.handling`, `pseudo.cnt` | `TRUE`, `"pseudo-count"`, `0.5` | 0.5 added to every count when the table has a zero, which is what R runs (below) |
| `p.adj.method` | `"BH"` | Benjamini-Hochberg |
| `alpha` | `0.05` | not an argument: `bt.da.consensus` makes the calls |

MicrobiomeStat 1.4 announces an imputation of zeros when library size depends on the model, but
its switch compares `"Imputation"` with `"imputation"` and never fires, so it always adds the
pseudocount. biotapy computes what R returns.

## Agreement with R

The golden test compares `effect`, `se`, `pvalue` and `qvalue` with `MicrobiomeStat::linda` on
the 636 GlobalPatterns genera in at least 20% of samples, human hosts against the rest, with and
without a sequencing-depth covariate: they agree to a relative 1e-7.

## Reference

Zhou H, He K, Chen J, Zhang X (2022) LinDA: linear models for differential abundance analysis of
microbiome compositional data. *Genome Biology* 23:95.
[doi:10.1186/s13059-022-02655-5](https://doi.org/10.1186/s13059-022-02655-5)
