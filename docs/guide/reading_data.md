# Reading and writing data

Readers turn a file format into a `TreeData` that follows the one
[data model](data_model.md) every biotapy function relies on.

## BIOM

`bt.io.read_biom` reads both BIOM dialects - JSON 1.0 and HDF5 2.1 - through
[biom-format](http://biom-format.org/), which detects the dialect from the file itself:

```python
import biotapy as bt

tdata = bt.io.read_biom("table.biom")
```

### Samples become rows

BIOM stores its matrix as features by samples. `read_biom` transposes it exactly
once, so `X` comes out samples (rows) by features (columns), like every other
biotapy object.

### Taxonomy becomes rank columns

Observation metadata's `taxonomy` entry - a list of ranks in HDF5, or whatever
string or list a producer wrote in JSON - is parsed into the canonical
`kingdom` .. `species` columns in `var`, the same way for every taxonomy
dialect (Greengenes, RESCRIPT, SILVA). A table with no `taxonomy` metadata
gets a `var` with no rank columns at all, rather than columns full of `NaN`.

### Sample metadata becomes obs columns

Whatever a table's sample metadata holds becomes columns of `obs`, one column
per key, unchanged.

### Ids are always strings

BIOM allows unquoted integer ids in its JSON dialect. `read_biom` casts every
sample and observation id to `str`, so `obs_names`/`var_names` are always
strings regardless of what the file stored.

### Duplicate ids

biom-format validates ids itself: a table with a repeated sample or
observation id raises `biom.exception.TableException` before biotapy ever
sees the table.

### Attaching a tree

Pass a Newick file's path as `tree=` to attach a phylogeny in `vart["phylo"]`.
If the tree's tips and the table's features disagree, only the shared
features are kept and one `UserWarning` names both counts - the same
tree/table alignment every biotapy reader uses.
