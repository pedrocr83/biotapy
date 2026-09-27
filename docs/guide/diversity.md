# Diversity

## Alpha diversity: `tl.alpha`

`bt.tl.alpha` computes one value per sample for each metric, through scikit-bio:

| Metric | What it is | In R |
|---|---|---|
| `observed_features` | features with a non-zero count | `estimate_richness`: `Observed` |
| `shannon` | Shannon entropy, natural log | `estimate_richness`: `Shannon` |
| `simpson` | Gini-Simpson, `1 - sum(p**2)` | `estimate_richness`: `Simpson` |
| `chao1` | bias-corrected Chao1 | `estimate_richness`: `Chao1` |
| `faith_pd` | Faith's phylogenetic diversity, root included | `picante::pd(include.root = TRUE)` |

```python
import biotapy as bt

tdata = bt.datasets.toy()
table = bt.tl.alpha(tdata)  # samples x metrics
bt.tl.alpha(tdata, metrics=["shannon", "faith_pd"], inplace=True)  # obs["alpha_shannon"], obs["alpha_faith_pd"]
```

- `observed_features` and `chao1` need raw counts, as `estimate_richness` does. `shannon`,
  `simpson` and `faith_pd` also run on relative abundances.
- An all-zero sample gets 0 for `observed_features`, `chao1` and `faith_pd`, and `NaN` for
  `shannon` and `simpson`; phyloseq reports Shannon 0 and Simpson 1.
- `faith_pd` needs a TreeData with a tree. A root with three or more children, common in a
  tree read from unrooted Newick, gets a zero-length split, which changes no root-to-tip
  distance.
- scikit-bio needs dense input, so biotapy densifies at most 2**20 values (8 MiB) at a time.
