# ANCOM-BC2

{func}`bt.da.ancombc2 <biotapy.da.ancombc2>` runs scikit-bio's ANCOM-BC2, which corrects each
sample's log counts for the share of its ecosystem that was sequenced before it fits one linear
model per feature. It runs in Python; no R is needed.

```python
import biotapy as bt

table = bt.da.ancombc2(bt.datasets.toy(), "group")
```

## Model

The log of an observed count is the log of the feature's abundance in the ecosystem plus the log
of the sample's sampling fraction $\theta_i$, which differs between samples and biases every
comparison of raw counts:

$$
\log o_{ij} = \theta_i + \log a_{ij}, \qquad \log a_{ij} = x_i^\top \beta_j + \varepsilon_{ij}
$$

ANCOM-BC2 (scikit-bio's {func}`~skbio.stats.composition.ancombc2`):

1. treats zeros as missing (no pseudocount) and fits one least-squares model per feature to the
   log counts on the design $x_i$ (an intercept, `group` and the `covariates`);
2. estimates the bias the coefficients of each term share with an E-M algorithm that models them as
   a mixture of a null, a negative and a positive normal component, at most 100 iterations;
3. estimates each sample's $\theta_i$ from the bias-corrected fit, subtracts it from the sample's
   log counts and fits the models again;
4. widens every variance by the bias estimate's own variance and then by the 5% quantile of all
   variances, so that rare features with tiny standard errors do not dominate, and tests
   $W_j = \hat\beta_j / \mathrm{se}_j$ two-sided against a t distribution whose degrees of freedom
   are the feature's observed samples minus the model's rank.

`qvalue` is the Benjamini-Hochberg correction over the features ANCOM-BC2 could test. A feature
whose zeros leave one level of `group` without an observed value, or that has no more observed
samples than model terms, cannot be tested: its row is NaN.

## Units

ANCOM-BC2 reports natural logs; biotapy divides `effect` and `se` by ln 2, so they are log2 like
every other method's. A numeric `group`'s effect is per unit, where `da.linda`'s and
`da.maaslin3`'s are per standard deviation.

## Compared with R

scikit-bio's implementation matches `ANCOMBC::ancombc2` (ANCOMBC 2.12.0) with these settings:

| `ancombc2` argument | R default | biotapy |
|---|---|---|
| `fix_formula` | required | `group` and `covariates` |
| `rand_formula` | `NULL` | not offered: no random effects |
| `p_adj_method` | `"holm"` | Benjamini-Hochberg |
| `prv_cut`, `lib_cut` | `0.10`, `0` | `0`, `0`: run `bt.pp.filter_features` once, before any method |
| `pseudo` | `0` | `0`: zeros are missing |
| `pseudo_sens` | `TRUE` | not run, so the table has no `passed_ss` flag |
| `s0_perc` | `0.05` | `0.05` |
| `struc_zero`, `neg_lb` | `FALSE` | not run |
| `em_control` | `tol = 1e-5`, `max_iter = 100` | the same |
| `global`, `pairwise`, `dunnet`, `trend` | `FALSE` | not offered: `group` has two levels or is numeric |
| `alpha` | `0.05` | not an argument: `bt.da.consensus` makes the calls |

R reports a feature it cannot test with p = 1 and counts it in the correction; biotapy reports
NaN and leaves it out, so its q-values are BH over the features actually tested.

## Agreement with R

The golden test compares biotapy with `ANCOMBC::ancombc2` on the 636 GlobalPatterns genera in at
least 20% of samples, human hosts against the rest. With a sequencing-depth covariate, `effect`,
`se` and `pvalue` agree to 1e-6. Without it, the bias E-M stops at R's cap of 100 iterations
before it has converged, on a slightly different iterate in scikit-bio: effects agree within
0.015 log2, standard errors to a relative 2e-3 and p-values within 0.02, and the significant
genera are the same. With 1,000 iterations R moves every effect by about -0.28 log2 and calls 220
genera instead of 208; biotapy keeps R's default.

## Choosing the reference

Unlike LinDA, ANCOM-BC2 is not antisymmetric in `reference`: the bias-corrected fit is made
against the reference level, in R as here, so swapping it changes more than the sign of
`effect`. On the GlobalPatterns genera the calls at q < 0.05 go from 208 to 230, and an effect
plus its swap is about -0.38 log2, not 0; `da.consensus` with LinDA goes from 104 to 112 genera.
Choose `reference` on the biology (the control or baseline level), not to change the results.

## Reference

Lin H, Peddada SD (2024) Multigroup analysis of compositions of microbiomes with covariate
adjustments and repeated measures. *Nature Methods* 21:83-91.
[doi:10.1038/s41592-023-02092-7](https://doi.org/10.1038/s41592-023-02092-7)
