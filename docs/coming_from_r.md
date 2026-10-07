# Coming from R

Every public biotapy function names its R equivalent in its docstring. This table is generated
from those lines each time the docs are built, plus a short list of phyloseq accessors that are
plain AnnData/TreeData code (`docs/_data/r_idioms.toml`). Rows marked "not in 0.3" have no
biotapy equivalent yet.

biotapy keeps samples as rows, so `tdata.X` is phyloseq's `otu_table` with
`taxa_are_rows = FALSE`. Functions return new objects instead of changing yours; `tl` functions
store their result in the object only with `inplace=True`.

```{include} generated/coming_from_r_table.md
```
