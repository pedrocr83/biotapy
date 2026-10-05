# Modules

* [core](core.md) - Private kernel: sparse group math, taxonomic rank order, function-table construction, x_kind/provenance/slot rules, and the sole gateway to TreeData and networkx.
* [io](io.md) - File readers and writer for BIOM, QIIME 2 artifacts, DADA2 sequence tables, phyloseq objects and MetaPhlAn profiles (each a TreeData through _core.make_treedata), HUMAnN and PICRUSt2 tables (a MuData through _core.make_function_mudata) and PICRUSt2 per-ASV trait tables (a DataFrame).
* [pp](pp.md) - Pure transforms over AnnData/TreeData that scale abundances per sample, aggregate features along the taxonomy, or filter samples/features.
* [tl](tl.md) - Diversity, ordination and PERMANOVA over AnnData/TreeData - alpha, beta, UniFrac, PCoA, NMDS and PERMANOVA through scikit-bio and scikit-learn, returning results or writing the data-model-slots keys.
* [pl](pl.md) - Plots of what tl, pp and fn give - stacked bars, heatmap, a function's contributions per taxon, richness, ordination and scree - drawn with matplotlib on the given or a new Axes, computing nothing.
* [fn](fn.md) - Function hierarchies, aggregation along them, HUMAnN-style renormalisation, per-taxon contributions and functional redundancy (Tian 2020) over function tables; owns no reader and no download.
* [datasets](datasets.md) - In-memory and pooch-cached example data for docs, doctests and tests - TreeData objects, a HUMAnN-style function MuData, the HMP2 cohort as a three-modality MuData, and the ENZYME hierarchy as an edge table.
