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
`dropna=False` to keep them instead. They are not merged into one shared group: an unassigned
feature is grouped like any other, by its full lineage up to `rank` (with `rank` itself counted as
`"NA"`), so two unassigned features under different lineages land in one group *each*, not
together:

```python
import biotapy as bt

tdata = bt.datasets.toy()
tdata.var.loc["f3", "genus"] = None  # f8 is already unassigned at genus; now f3 is too
bt.pp.tax_glom(tdata, "genus", dropna=False).n_vars  # 7 groups, not 6: f3 and f8 stay apart
```

## What happens to the tree and to derived slots

`tax_glom` goes through the same feature-changing path as every other biotapy function that
drops taxa: `layers`, `obsm`, `obsp`, `varm`, `varp` and every non-`biotapy` key in `uns` are
dropped, since they described the old features (a plotted `uns["group_colors"]`, for example,
disappears along with them). On a `TreeData`, the phylogeny in `vart["phylo"]` is pruned to the
kept archetypes' subtree - unary internal nodes stay, so the *kept* tips' root-to-tip path lengths
are unchanged. Faith PD and UniFrac still change, though: removing tips removes their branches
from the tree, which is exactly what those metrics measure.
