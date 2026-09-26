# Aggregation

`bt.pp.tax_glom` merges features that share the same taxonomy down to a chosen rank, the
biotapy equivalent of phyloseq's `tax_glom`.

## Grouping by lineage, not by label

Features are grouped by their full lineage up to and including `rank` - kingdom through that
rank, joined into one string - not by the value at `rank` alone. Two features both labelled
`"uncultured"` at genus but sitting under different families stay separate groups, exactly as
in phyloseq:

```python
import biotapy as bt

tdata = bt.datasets.toy()
out = bt.pp.tax_glom(tdata, "phylum")
list(out.var_names)  # ["f3", "f6", "f7"]
```

## The archetype

Each group is represented by one of its members: the most abundant one (by total count across
samples), first on a tie. That feature's row in `var` - and its place in the tree - is what the
merged feature keeps. Every rank below `rank` is set to `NaN` on the output, since the merged
feature no longer has a single value there.

## Unassigned features

By default (`dropna=True`) features with no value at `rank` are dropped before grouping. Pass
`dropna=False` to keep them as one extra group instead.

## What happens to the tree and to derived slots

`tax_glom` goes through the same feature-changing path as every other biotapy function that
drops taxa: `layers`, `obsm` and `obsp` are dropped, since they described the old features. On
a `TreeData`, the phylogeny in `vart["phylo"]` is pruned to the kept archetypes' subtree -
unary internal nodes stay, so root-to-tip path lengths (and therefore Faith PD and UniFrac) are
unaffected by the merge.
