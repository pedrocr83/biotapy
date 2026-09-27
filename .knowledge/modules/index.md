# Modules

* [core](core.md) - Private kernel: sparse group math, taxonomic rank order, x_kind/provenance/slot rules, and the sole gateway to TreeData and networkx.
* [io](io.md) - File readers and writer for BIOM, QIIME 2 artifacts and DADA2 sequence tables, building every TreeData through _core.make_treedata.
* [pp](pp.md) - Pure transforms over AnnData/TreeData that scale abundances per sample, aggregate features along the taxonomy, or filter samples/features.
* [datasets](datasets.md) - In-memory and (future) cached example TreeData objects for docs, doctests and tests.
