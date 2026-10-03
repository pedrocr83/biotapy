# Modules

* [core](core.md) - Private kernel: sparse group math, taxonomic rank order, x_kind/provenance/slot rules, and the sole gateway to TreeData and networkx.
* [io](io.md) - File readers and writer for BIOM, QIIME 2 artifacts, DADA2 sequence tables and phyloseq objects, building every TreeData through _core.make_treedata.
* [pp](pp.md) - Pure transforms over AnnData/TreeData that scale abundances per sample, aggregate features along the taxonomy, or filter samples/features.
* [tl](tl.md) - Diversity, ordination and PERMANOVA over AnnData/TreeData - alpha, beta, UniFrac, PCoA, NMDS and PERMANOVA through scikit-bio and scikit-learn, returning results or writing the data-model-slots keys.
* [pl](pl.md) - Plots of what tl and pp stored - stacked bars, heatmap, richness, ordination and scree - drawn with matplotlib on the given or a new Axes, computing nothing.
* [datasets](datasets.md) - In-memory and (future) cached example TreeData objects for docs, doctests and tests.
