# Decisions

* [TreeData and MuData are the only containers](treedata-as-container.md) - biotapy owns no data class; one data type is a TreeData, several are a MuData, every feature is a function over them.
* [Samples are rows](samples-as-rows.md) - Every matrix is samples x features (scverse orientation); importers transpose once, no other function checks orientation.
* [Pure by default, one inplace convention for tl](pure-by-default.md) - io/pp return new objects and never mutate input; tl returns results, and inplace=True writes them to the documented slot; pl returns Axes.
* [Python first, compiled code last](python-first-compiled-last.md) - Pure Python/NumPy by default; delegate to compiled libraries; Numba then Rust only for a benchmarked hotspot; never new C/C++.
* [Heavy dependencies are optional extras](optional-heavy-dependencies.md) - torch, rpy2, plotnine, numba and unifrac install only through extras and are imported lazily.
* [Knowledge in OKF, user docs in Sphinx](docs-okf-and-sphinx.md) - Contributor and agent knowledge is an OKF v0.2 bundle in .knowledge/; user docs are a Sphinx site in docs/.
* [R bridge before native ports](r-bridge-before-ports.md) - R-only DA methods ship first through an optional rpy2 bridge; native ports only for the most used.
* [No bundled KEGG mapping files](no-bundled-kegg.md) - Functional hierarchy mappings are downloaded on first use, never shipped in the wheel.
* [Package name biotapy](package-name-biotapy.md) - Distribution and import name is biotapy, hosted at github.com/pedrocr83/biotapy.
* [Read phyloseq via rdata, refseq warned and skipped](phyloseq-import-route.md) - A ~35-line rdata constructor_dict reads GlobalPatterns/enterotype/esophagus with zero shape mismatches and zero residual warnings; native rdata (+xarray) is the route for Tasks 1.9/1.10, an R export script is rejected as the default, and a populated refseq is warned-and-skipped rather than guessed at (no test fixture has one).
