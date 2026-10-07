# Consensus

{func}`bt.da.consensus <biotapy.da.consensus>` puts the result tables of several methods side by
side and counts, per feature, the methods that call it and whether they agree on its direction.
It runs no method itself: you choose and run the methods, and pass their tables.
{func}`bt.pl.consensus <biotapy.pl.consensus>` draws the table it returns.

```python
import biotapy as bt

tdata = bt.datasets.toy()
results = [bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")]
table = bt.da.consensus(results, alpha=0.05, min_methods=2)
```

## Rule

For feature $f$ and method $m$, with $q_{fm}$ the method's `qvalue` and $d_{fm}$ its `direction`:

| Column | Definition |
|---|---|
| `significant_<m>` | $q_{fm} < \alpha$, strictly: $q = \alpha$ is not a call |
| `n_tested` | methods whose `pvalue` for $f$ is finite; a method that could not test $f$, or whose table lacks it, does not count |
| `n_significant` | methods that call $f$ |
| `direction` | the sign $d_{fm}$ every calling method shares; 0 when none calls $f$ or their signs differ |
| `consensus` | `n_significant` $\ge$ `min_methods` and `direction` is not 0 |
| `conflict` | calling methods give $f$ opposite signs; never a consensus |

A call whose `effect` is exactly 0 has no direction: it counts in `n_significant`, but then
`direction` is 0, with no consensus and no conflict. The table also carries each method's
`effect_<m>` and `qvalue_<m>`, in the order of `results`.

## Why it is defined this way

- **One correction.** Every method's `qvalue` is Benjamini-Hochberg over the features it tested,
  so "significant" means the same false discovery rate in every column. R's defaults would mix
  Holm (ANCOM-BC2), BH over the group's p-values (LinDA), BH pooled with the covariates'
  (MaAsLin 3) and averaged per-draw BH (ALDEx2).
- **Strict calls, recomputed.** Calls come from `qvalue` with one rule, never from a method's own
  flag: LinDA's `reject` uses $q \le \alpha$, and ANCOM-BC2's `Signif` sits on Holm by default.
- **Untested is not "not significant".** ANCOM-BC2 and MaAsLin 3 cannot fit a feature that is
  absent from one group; counting that as a negative would let one method veto a feature it never
  tested. `min_methods=len(results)` still asks for every method.
- **Same contrast.** The tables must come from different methods and compare the same `contrast`;
  `bt.da.consensus` raises otherwise, so run every method with the same `group` and `reference`.
- **The user picks the methods.** A one-call consensus would hide which methods ran and invite
  trying methods until one agrees. Methods disagree a lot on real data: Nearing et al. (2022)
  recommend a consensus of several methods, Pelto et al. (2025) one simple method. Either way,
  choose the methods before you look at their results. Methods that share a model, such as
  ANCOM-BC and ANCOM-BC2, agree more often for that reason alone.

The rule is recorded, with the options rejected, in the knowledge bundle's decision
[da-consensus-agreement](https://github.com/pedrocr83/biotapy/blob/master/.knowledge/decisions/da-consensus-agreement.md).
There is no R equivalent.

## References

Nearing JT, Douglas GM, Hayes MG, MacDonald J, Desai DK, Allward N, Jones CMA, Wright RJ, Dhanani
AS, Comeau AM, Langille MGI (2022) Microbiome differential abundance methods produce different
results across 38 datasets. *Nature Communications* 13:342.
[doi:10.1038/s41467-022-28034-z](https://doi.org/10.1038/s41467-022-28034-z)

Pelto J, Auranen K, Kujala JV, Lahti L (2025) Elementary methods provide more replicable results
in microbial differential abundance analysis. *Briefings in Bioinformatics* 26:bbaf130.
[doi:10.1093/bib/bbaf130](https://doi.org/10.1093/bib/bbaf130)
