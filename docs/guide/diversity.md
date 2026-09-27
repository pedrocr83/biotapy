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

## Beta diversity: `tl.beta` and `tl.unifrac`

`bt.tl.beta` computes Bray-Curtis (`metric="braycurtis"`, the default) or Jaccard
(`metric="jaccard"`) distances between every pair of samples; `bt.tl.unifrac` computes
unweighted or weighted UniFrac along the tree:

```python
bt.tl.beta(tdata, inplace=True)  # obsp["braycurtis"]
bt.tl.unifrac(tdata, weighted=True, inplace=True)  # obsp["weighted_unifrac"]
```

- Jaccard is on presence/absence, like `phyloseq::distance(physeq, "jaccard", binary = TRUE)`.
  Without `binary = TRUE` phyloseq computes vegan's quantitative Jaccard, a different number.
- Weighted UniFrac is normalized to 0-1 by default, as in phyloseq; pass `normalized=False` for
  the raw value.
- Weighted UniFrac needs raw counts: `x_kind` `"counts"` and whole numbers. scikit-bio's tree code
  casts abundances to integers, so proportions would all become 0. Unweighted UniFrac uses
  presence only and runs on relative abundances too.
- A tree whose root has three or more children is used rooted where it is drawn. phyloseq
  instead roots such a tree at a random tip, so its UniFrac changes from run to run.
- Two all-zero samples are `NaN` apart under Bray-Curtis and 0 apart under Jaccard and UniFrac
  (scikit-bio's convention; vegan's binary Jaccard gives `NaN`). Drop empty samples with
  `bt.pp.filter_samples(tdata, 1)` before ordinating.
- The table is densified once: 8 bytes x samples x features, plus 8 bytes x samples x samples
  for the result.
