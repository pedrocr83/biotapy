# Ordination and PERMANOVA

Ordination and PERMANOVA read a distance matrix that `tl.beta` or `tl.unifrac` stored in
`obsp`; `distance=` names it (default `"braycurtis"`). When it is missing, the error names the
call to run.

```python
import biotapy as bt

tdata = bt.datasets.toy()
bt.tl.beta(tdata, inplace=True)
coords, axes = bt.tl.pcoa(tdata)  # PC1.. per sample; eigenvalue and proportion per axis
coords, stress = bt.tl.nmds(tdata, seed=0)  # NMDS1, NMDS2 and Kruskal stress-1
result = bt.tl.permanova(tdata, "group", seed=0)
result["test statistic"], result["p-value"]
```

## PCoA

`bt.tl.pcoa` matches `ordinate(physeq, "PCoA", distance)`, which runs `ape::pcoa` with no
correction for negative eigenvalues. `proportion_explained` is each eigenvalue over the sum of
all of them, ape's `Relative_eig`. At most `n_obs - 1` axes are returned, and axis signs are
arbitrary, in R too. With `inplace=True` the coordinates go to `obsm["X_pcoa"]` and the axes to
`uns["biotapy"]["pcoa"]`.

When the distances are not Euclidean, some eigenvalues can be negative. scikit-bio sets them, and
their coordinates, to 0, whereas `ape::pcoa` reports those eigenvalues as negative. The
proportions divide by the trace, the sum of all eigenvalues with the negative ones included, so
they still match ape's `Relative_eig` on the positive axes.

## NMDS

`bt.tl.nmds` runs scikit-learn's non-metric SMACOF on the stored distances and keeps the best
of 20 random starts, as `vegan::metaMDS` does by default. vegan then rotates and rescales its
result (`postMDS`); biotapy does not, so compare configurations up to rotation and scale. With
`inplace=True`: `obsm["X_nmds"]` and `uns["biotapy"]["nmds"]["stress"]`.

## PERMANOVA

`bt.tl.permanova(tdata, grouping)` tests whether the groups in `obs[grouping]` differ, like
`vegan::adonis2(distance ~ grouping, data)`. The pseudo-F statistic matches vegan exactly;
p-values come from permutations and agree only up to that noise. Drop samples with a missing
group first.

## After filtering

`pp.filter_features` and `pp.rarefy` drop stored distances and ordinations, including
`uns["biotapy"]["pcoa"]` and `["nmds"]`, since they described the old features.
`pp.filter_samples` keeps them, subset to the remaining samples: distances stay valid, but an
ordination is not recomputed, so run `tl.pcoa` again for the ordination of the subset.
