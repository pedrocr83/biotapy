# Function

`bt.fn` works on functional profiles - gene families, EC numbers, pathways -
the way `bt.pp.tax_glom` works on taxonomy. Its rules follow HUMAnN's own
utility scripts, and biotapy's golden tests compare the two on the same files.

## Two tables per HUMAnN file

A HUMAnN table holds community rows (`PWY-5100`) and the same functions split
by taxon (`PWY-5100|g__Bacteroides.s__Bacteroides_ovatus`). For pathways the
community value is not the sum of its strata, so neither can be derived from
the other. `bt.io.read_humann` therefore returns a
[MuData](https://mudata.scverse.org/) with two modalities over the same
samples:

| Modality | Rows | `var` columns |
|---|---|---|
| `"function"` | community rows | `name`, `special` |
| `"function_by_taxon"` | stratified rows | `function`, `name`, `taxon`, `genus`, `species`, `special` |

```python
import biotapy as bt

mdata = bt.io.read_humann("pathabundance.tsv")
mdata["function"]  # samples x pathways
mdata["function_by_taxon"]  # samples x (pathway, taxon) pairs
```

`UNMAPPED`, `READS_UNMAPPED`, `UNINTEGRATED` and `UNGROUPED` stay features,
flagged in `var["special"]`, so a sample's total keeps what HUMAnN could not
assign. Read one table per call: a gene family table and a pathway table both
hold `UNMAPPED`.

## Aggregating along a hierarchy

`bt.fn.func_glom` sums functions into their parents at one level of a
hierarchy, with `humann_regroup_table`'s rules:

- **Many-to-many.** A function with two parents counts in full toward both,
  so a level's total can exceed the input's.
- **`UNGROUPED`.** Functions with no parent at that level are summed into
  `UNGROUPED`, per taxon in a stratified table.
- **Specials pass through.** `UNMAPPED`, `READS_UNMAPPED` and `UNINTEGRATED`
  keep their own rows. (HUMAnN 3.9 sums `READS_UNMAPPED` into `UNGROUPED`;
  biotapy follows HUMAnN master, which keeps it.)
- **`agg="mean"`** divides by the members present in the table, not by the
  group's size in the hierarchy.

Apply it to each modality to keep the pair together:

```python
import mudata

import biotapy as bt

mdata = bt.datasets.toy_humann()
edges = bt.datasets.enzyme()
by_class = mudata.MuData({key: bt.fn.func_glom(mod, "class", hierarchy=edges) for key, mod in mdata.mod.items()})
```

A hierarchy is an edge table with columns `child`, `parent`, `level` and,
optionally, `parent_name`, which becomes `var["name"]`. Every ancestor is a
row, not only the direct parent, so a table already grouped one level up
still reaches the levels above.

## Where hierarchies come from

- **EC numbers:** `bt.datasets.enzyme()` downloads the open ENZYME hierarchy
  (CC BY 4.0) once and caches it; see the [datasets guide](datasets.md).
- **Your own mapping files:** `bt.fn.load_hierarchy(path, level)` reads
  `humann_regroup_table --custom` files and PICRUSt2 mapping files
  (`parent<TAB>child<TAB>...`, the default `layout="parent_first"`), or
  `child<TAB>parent...` tables with `layout="child_first"`.

biotapy ships and downloads no KEGG or MetaCyc mapping: their licences do not
allow it to. If you hold a KEGG or MetaCyc licence, export the mapping you
need and read it with `load_hierarchy`.
