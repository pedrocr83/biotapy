# ALDEx2

{func}`bt.da.aldex2 <biotapy.da.aldex2>` runs `ALDEx2::aldex` in R through rpy2: it draws
plausible relative abundances for each sample, takes their log2 centred log-ratios and compares
the two groups with a Welch t-test on every draw. It needs R, the R package ALDEx2 and biotapy's
`r` extra ({ref}`Methods that run in R <da-methods-in-r>`).

```python
import biotapy as bt

table = bt.da.aldex2(bt.datasets.toy(), "group", seed=0)
```

## Model

For each sample $i$, ALDEx2 draws `mc_samples` proportion vectors from a Dirichlet distribution
whose parameters are the sample's counts plus 0.5, so a draw reflects how uncertain proportions
are at that sequencing depth. Each draw is turned into log2 centred log-ratios over all features
(`denom = "all"`):

$$
c^{(k)}_{ij} = \log_2 p^{(k)}_{ij} - \frac{1}{m} \sum_{l=1}^{m} \log_2 p^{(k)}_{il},
\qquad p^{(k)}_{i\cdot} \sim \operatorname{Dirichlet}(y_{i\cdot} + 0.5)
$$

For each draw $k$ and feature $j$, a Welch t-test compares the two groups' values $c^{(k)}_{ij}$;
`pvalue` is the mean of the draws' two-sided p-values (ALDEx2's `we.ep`). `effect` is ALDEx2's
`diff.btw`, the median of the differences between the two groups' values over all draws, already
log2. ALDEx2
reports no standard error, so `se` is NaN for every feature, and `qvalue` is the
Benjamini-Hochberg correction of `pvalue`, as for every method. A feature with no read in any
sample is dropped by ALDEx2 and is NaN in the table (not tested).

## Units

`effect` is the median log2 difference of the centred log-ratios between the other level and
`reference`. ALDEx2 compares two groups only: a numeric `group`, covariates and a level in fewer
than two samples are refused.

## Compared with R

biotapy calls `ALDEx2::aldex` (ALDEx2 1.42.0) with its defaults:

| `aldex` argument or output | R default | biotapy |
|---|---|---|
| `conditions` | required | the labels `"0"` (`reference`) and `"1"`: ALDEx2 fails on a factor and sorts character labels in the R session's locale |
| `mc.samples` | `128` | `mc_samples=128`; R warns below 128, and biotapy re-emits it as a `UserWarning` |
| `test`, `effect`, `denom` | `"t"`, `TRUE`, `"all"` | the same |
| `paired.test`, `iterate`, `gamma` | `FALSE`, `FALSE`, `NULL` | the same: no paired test, no scale model |
| `we.ep` | output | `pvalue` |
| `we.eBH` | output | not carried: it averages the draws' corrected values, a different quantity; `qvalue` is BH of `we.ep` |
| `diff.btw` | output | `effect` |
| `effect`, `overlap`, `wi.ep`, `wi.eBH`, `rab.*`, `diff.win` | outputs | not carried: ALDEx2's `effect` is a standardised size, not a fold change |
| random state | R's global seed | `seed` sets R's seed for the call, and the R session's random state is put back afterwards |

On the GlobalPatterns genera below, BH of `we.ep` calls 11 genera and ALDEx2's `we.eBH` 19.

## Agreement with R

ALDEx2 is random, but `seed` becomes one integer for R's `set.seed`, so biotapy and R draw the same
numbers. The golden test runs `set.seed(...); ALDEx2::aldex(...)` with that integer on the 636
GlobalPatterns genera in at least 20% of samples, human hosts against the rest: `effect` and
`diff.btw` agree to a relative 4e-15 and `pvalue` and `we.ep` to 8.4e-13 (checked at 1e-7). The
golden passes `reference="human"`, the level R sorts first, so both sides draw in the same label
order. With other seeds, effects correlate
at 0.993 (Spearman) and 11 or 12 genera are called.

## Choosing the reference

Swapping `reference` does more than flip the sign: ALDEx2 takes its Monte Carlo draws in label
order, so the same `seed` gives different effects, and the two runs are not exact mirror images.
They differ from exact antisymmetry by up to 0.25 log2 on `toy()` (whose effects are about 4 log2
wide) and up to 0.65 log2 on the GlobalPatterns genera, 15 of 636 of which then do not change
direction. The p-values and the calls are the same.

## Reference

Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB (2014) Unifying the
analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing
and selective growth experiments by compositional data analysis. *Microbiome* 2:15.
[doi:10.1186/2049-2618-2-15](https://doi.org/10.1186/2049-2618-2-15)
