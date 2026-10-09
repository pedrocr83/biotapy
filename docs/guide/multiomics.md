# Multi-omics

When several data types are measured on the same samples - taxa, a function
table, metabolites, host data - biotapy keeps them in one
[MuData](https://mudata.scverse.org/): one AnnData per data type, called a
modality, over shared samples.

## Modality names

| Modality | Holds |
|---|---|
| `"taxa"` | taxa, usually a TreeData with its phylogeny |
| `"function"` | a function table's community rows |
| `"function_by_taxon"` | a function table's per-taxon rows |
| `"metabolites"` | metabolite intensities or counts |
| `"host"` | host measurements, such as transcripts or clinical values |

biotapy functions that read a modality take its name as an argument and
default to these names. Other names work too.

## Building one

`bt.io.to_mudata` takes a mapping of name to AnnData:

```python
import biotapy as bt

table = bt.io.read_humann("pathabundance.tsv")  # already a MuData
mdata = bt.io.to_mudata({"taxa": tdata, **table.mod, "metabolites": metabolites})
```

A function table is already a MuData of two modalities, and a MuData cannot
hold another MuData, so `**table.mod` adds `function` and `function_by_taxon`
as two entries.

Every modality keeps only the samples that all of them have, in the first
modality's order, and biotapy warns with how many samples each modality lost.
Each modality is a copy, so the tables you passed are unchanged.

The global `mdata.obs` starts without columns; `mdata.pull_obs()` gathers
the modalities' sample metadata into it.

## Saving

`mdata.write_h5mu("study.h5mu")` saves every modality, but a TreeData
modality reads back as a plain AnnData without its tree. Save that modality
with its tree as well:

```python
mdata.write_h5mu("study.h5mu")
mdata["taxa"].write_h5td("taxa.h5td")
```
