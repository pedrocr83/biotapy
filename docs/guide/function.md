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

## PICRUSt2 tables

`bt.io.read_picrust2` gives PICRUSt2 predictions the same two modalities. Its
`"function_by_taxon"` features are functions per ASV (`2.7.1.1|ASV1`), read
from the long contribution table. For gene families a sample's
contributions sum to its community value; for pathways they need not, as in
HUMAnN. EC numbers lose PICRUSt2's `EC:` prefix, so they match
`bt.datasets.enzyme()`:

```python
import biotapy as bt

mdata = bt.io.read_picrust2("pred_metagenome_unstrat.tsv.gz", contrib="pred_metagenome_contrib.tsv.gz")
by_class = bt.fn.func_glom(mdata["function"], "class", hierarchy=bt.datasets.enzyme())
```

PICRUSt2's own mapping files write `EC:1.1.1.1`; remove the prefix from the
edge table's `child` column before regrouping with them, or `func_glom`
raises because nothing maps.

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
- **Units.** The output keeps `uns["biotapy"]["x_kind"]` only for a sum in
  which every function has at most one parent at the level. A mean, or a
  sum in which a function counts toward several parents, is labelled
  `"abundance"`: its values are no longer counts of reads or proportions that
  sum to 1, so `bt.pp.rarefy` refuses it.

If no function of a community table has a parent at the level, `func_glom`
raises and shows a few ids from each side: this is almost always an id
format mismatch, such as `EC:1.1.1.1` in the table against `1.1.1.1` in the
hierarchy. A stratified table does not raise, because a pathway need not
have strata; its rows go to `UNGROUPED|<taxon>`, as in HUMAnN.

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
still reaches the levels above. The output's provenance records the
hierarchy's `attrs["source"]` and, when it has one, `attrs["license"]`.
`pd.concat` keeps `attrs` only when every input has the same ones, so
`pd.concat([bt.datasets.enzyme(), my_edges])` is recorded with no source or
licence: set `attrs` on the result yourself.

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

## Renormalising

`bt.fn.renorm(mdata, "relab")` (or `"cpm"`) divides every row of both
modalities by the sample's community total, as `humann_renorm_table` does by
default. Stratified rows are divided by the community total too, so a
pathway's strata keep their share of the community. `special=False` drops
the special rows before totalling, as `--special n` does.

HUMAnN's `--mode levelwise`, where each modality is scaled by its own total,
is `bt.pp.relative` applied to each modality.

A sample whose community total is zero stays all zero, and `renorm` warns
once, naming up to three such samples and how many there are, as HUMAnN does.

## Contributions

`bt.fn.contributions` answers "which taxa carry this function?". It takes the
stratified modality and one function id, and returns a samples x taxa
`DataFrame` of that function's strata, as stored:

```python
import biotapy as bt

mdata = bt.datasets.toy_humann()
table = bt.fn.contributions(mdata["function_by_taxon"], "2.7.1.2", top=5)
```

- **Raw strata.** Each row sums to the function's strata in that sample. For
  pathways that is not the community value, because a pathway's strata need
  not add up to it. Pass the stratified modality of
  `bt.fn.renorm(mdata, "relab")` to read the values as shares of each
  sample's community total.
- **Order.** Taxa are ordered by their total over all samples, largest first
  (ties by id). `top=5` keeps the five largest and sums the rest into a last
  column, `"other"`.
- **Every stratum is a taxon.** `unclassified` (HUMAnN) and `RARE`
  (PICRUSt2) are columns like any other; special functions such as
  `UNINTEGRATED` can be queried too.
- **Regrouped tables work too.** `bt.fn.func_glom`'s output for the
  stratified modality has the same `var` columns, so
  `contributions(by_class, "2.-.-.-")` shows the taxa behind an enzyme class.

An unknown id raises a `KeyError` listing up to three close ids, which
catches typos and an `EC:` prefix.

`bt.pl.contributions` draws the same table as stacked bars, one per sample;
see the [plotting guide](plotting.md).

## Functional redundancy

`bt.fn.functional_redundancy(adata, traits=traits)` measures, per sample, how
much the taxa present overlap in what their genomes can do. It follows Tian
et al. (2020, *Nature Communications* 11:6217; the equation numbers below are
theirs). For the taxa of a sample, with relative abundances $p_i$ summing
to 1, and a functional distance $d_{ij}$ between taxa $i$ and $j$:

$$
\mathrm{TD} = \sum_{i} \sum_{j \ne i} p_i p_j = 1 - \sum_i p_i^2 \qquad \text{(Gini-Simpson, Eq. 2)}
$$

$$
\mathrm{FD} = \sum_{i} \sum_{j \ne i} d_{ij} p_i p_j \qquad \text{(Rao's quadratic entropy, Eq. 3)}
$$

$$
\mathrm{FR} = \mathrm{TD} - \mathrm{FD} = \sum_{i} \sum_{j \ne i} (1 - d_{ij}) p_i p_j \qquad \text{(Eqs. 1 and 4)}
$$

$d_{ij}$ is the weighted Jaccard distance between the two taxa's gene copy
numbers $G_{ia}$ (Eq. 7):

$$
d_{ij} = 1 - \frac{\sum_a \min(G_{ia}, G_{ja})}{\sum_a \max(G_{ia}, G_{ja})}
$$

FR is the chance that two individuals drawn from the sample belong to
different taxa, weighted by how much of their genomes those taxa share. The
result also holds `normalized_redundancy`, $\mathrm{FR} / \mathrm{TD}$, which
compares samples whose taxonomic diversity differs; Tian et al. report about
0.4 for most human body sites. biotapy computes $d_{ij}$ with SciPy's Bray-Curtis
dissimilarity $\mathrm{BC}$, as $2\,\mathrm{BC} / (1 + \mathrm{BC})$, which
equals the weighted Jaccard distance for non-negative vectors. Two taxa
with no gene at all share nothing: their distance is 1.

**Inputs.** `traits` is a taxa x genes table of copy numbers, such as
PICRUSt2's per-ASV predictions read by `bt.io.read_picrust2_traits`.
`adata` holds the same taxa as features. Tian et al. use relative organism
abundances (MetaPhlAn2 profiles), so divide 16S read counts by each ASV's
predicted 16S copy number first. PICRUSt2 writes those copy numbers to
`marker_predicted_and_nsti.tsv.gz`, which the same reader reads:

```python
import biotapy as bt

asvs = bt.io.read_biom("table.biom")  # samples x ASVs, read counts
traits = bt.io.read_picrust2_traits("picrust2_out/EC_predicted.tsv.gz")
copies = bt.io.read_picrust2_traits("picrust2_out/marker_predicted_and_nsti.tsv.gz")["16S_rRNA_Count"]

cells = asvs[:, asvs.var_names.isin(copies.index)].copy()
cells.X = cells.X.multiply(1 / copies[cells.var_names].to_numpy()).tocsr()
fr = bt.fn.functional_redundancy(cells, traits=traits)
```

The corrected table is no longer read counts; use it for this calculation
only. A MetaPhlAn profile is already in organism abundances.

- **Taxa without traits** (ASVs PICRUSt2 dropped above its NSTI cutoff) are
  left out with a warning, abundance included, so each sample's $p$ sums to
  1 over the taxa that have a genome, as in the paper. PICRUSt2's `RARE`
  group exists only in contribution tables, never in an ASV table.
- **Empty samples.** A sample with no abundance on a taxon with traits gets
  NaN in every column; a sample with a single such taxon has no diversity
  (0) and an undefined `normalized_redundancy` (NaN).
- **Size.** The distances form a dense taxa x taxa matrix, and the time grows
  as taxa² x genes: 2,000 ASVs take seconds, 10,000 take minutes and over a
  gigabyte. Filter rare ASVs first with `bt.pp.filter_features`.
