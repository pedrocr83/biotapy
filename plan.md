# Python Microbiome Toolkit: Development Report

Sep 26, 2026 · @Pedro Ribeiro

Build a mia-style microbiome package for Python on AnnData/TreeData, reuse scikit-bio for the math, keep every feature a small pure function, and move to compiled code only where profiling proves Python is the bottleneck.

## Positioning

The package is "mia for Python": a microbiome data object and analysis toolkit that fits the scverse ecosystem, treats function as seriously as taxonomy, and is ready for machine learning out of the box.

| Existing tool | What it does well | Gap this package fills |
| --- | --- | --- |
| [pyloseq](https://pypi.org/project/pyloseq/) (v1.1.1, June 2026) | 1:1 port of phyloseq, R golden-file tests, reads .qza | Custom object, R-style API, no AnnData/scverse interop, no functional or multi-omics layer |
| [scikit-bio](https://github.com/scikit-bio/scikit-bio/releases) 0.7.x | Diversity, UniFrac, ordination, PERMANOVA, ANCOM-BC, mmvec, optional C++ binaries | No container object, no phyloseq-style verbs, no taxonomy or function hierarchy tools |
| [treedata](https://pypi.org/p/treedata) / AnnData | Annotated matrix with trees, h5ad/zarr, scverse compatible | Nothing microbiome-specific |
| QIIME 2 | End-to-end CLI pipelines | Heavy framework, not a notebook-first library |
| R mia / miaverse | Reference design (TreeSummarizedExperiment, miaTime, miaViz, OMA book) | R only |

Target users: microbiome researchers leaving R, Python ML practitioners who get handed microbiome data, and pharma teams running multi-cohort or multi-omics studies.

Three differentiators, in priority order:

1. **Function as a first-class hierarchy.** KO, module and pathway get the same glom, filter and plot verbs as kingdom to species, with HUMAnN stratified output kept as a taxa-to-function link.
2. **Differential abundance consensus.** One call runs several methods (MaAsLin 3, ANCOM-BC, LinDA, ALDEx2) and reports where they agree.
3. **ML-ready by default.** scikit-learn transformers for leak-free preprocessing, embeddings from microbiome foundation models in `obsm`, PyTorch loaders.

## Engineering rules

Every feature is the smallest pure function that solves one task, delegates to an existing library whenever one does the job, and carries its explanation in the docstring and docs site rather than in inline comments.

### Minimal code per task

- Before writing anything, check in order: does AnnData, scikit-bio, SciPy, NumPy, pandas or scikit-learn already do it? If yes, the function is a thin wrapper, often 3 to 10 lines.
- No speculative parameters, options or abstractions. Add a parameter only when a real use case or test needs it.
- No class hierarchies. The data object is TreeData/MuData; everything else is a function.
- Soft budget: a function body over \~30 lines is a review flag. Split it or find the library call you missed.
- Diffs stay minimal: a change touches only the functions the task needs, no drive-by refactors.
- Dependencies are weighed like code: each new one needs a written reason in the PR.

### Function-level implementation

- One task = one public function, in the style `verb(adata, ...) -> AnnData | DataFrame | ndarray`.
- Pure by default: return a new object; `inplace=True` only where scanpy users expect it, and then consistently.
- Functions are grouped by module (`pp`, `tl`, `pl`, `io`, `da`, `fn`) following scanpy conventions, so users already know the layout.
- Each function is independently testable and independently replaceable by a compiled version later without changing its signature.
- Private helpers are allowed only when shared by two or more public functions.

### Minimal comments, rich docs

- Inline comments only for the non-obvious *why* (a numerical trick, an R-compatibility quirk). Never for *what*: the name and types say that.
- Every public function has a NumPy-style docstring: summary, parameters, returns, an `R equivalent` line (phyloseq/mia function), a runnable example, and references.
- Longer explanation (method choice, math, caveats) lives in the docs site, linked from the docstring.
- Full type hints on every public function; they replace most comments.

Example of the target shape:

```python
def tax_glom(adata: AnnData, rank: str, *, dropna: bool = True) -> AnnData:
    """Aggregate features to a taxonomic rank.

    R equivalent: phyloseq::tax_glom, mia::agglomerateByRank.
    See docs/guide/aggregation.md for how unassigned taxa are handled.
    """
    groups = adata.var[rank]
    if dropna:
        adata = adata[:, groups.notna()]
        groups = groups.dropna()
    return ad.AnnData(
        X=sum_by(adata.X, groups.values),
        obs=adata.obs,
        var=first_by(adata.var, groups),
    )
```

### Enforcing it

- Ruff with a max-complexity rule and a docstring rule on public functions.
- A PR template with three checkboxes: reused an existing library call, no new abstraction, docstring plus docs page updated.
- Coding-agent instructions (CLAUDE.md or equivalent) repeat these rules, so AI-assisted contributions follow them too.

## Architecture

The package owns no data class: one data type is a TreeData (AnnData plus trees), several data types are a MuData, and every feature is a function over those.

### Data model

| Slot | Holds | phyloseq / mia equivalent |
| --- | --- | --- |
| `X` (sparse CSR) | Counts, samples × features | `otu_table` / `assay("counts")` |
| `layers` | Relative abundance, CLR, rarefied | extra assays |
| `obs` | Sample metadata | `sample_data` / `colData` |
| `var` | Taxonomy, one column per rank; sequence as a string column | `tax_table`, `refseq` / `rowData` |
| `vart["phylo"]` | Phylogenetic tree over features | `phy_tree` / `rowTree` |
| `obsm` | Ordinations, foundation-model embeddings | `reducedDim` |
| `obsp` | Beta-diversity distance matrices | none (stored separately in R) |
| `uns` | Provenance: importer, versions, parameters | metadata |
| MuData modalities | `taxa`, `function`, `metabolites`, `host` | `altExp` / MultiAssayExperiment |

Orientation follows scverse (samples are rows), the opposite of phyloseq's usual default. Importers handle the flip once, so no function has to check.

### Module layout

| Module | Scope | Built on |
| --- | --- | --- |
| `io` | phyloseq .rds, BIOM, .qza, DADA2, HUMAnN 4, PICRUSt2, MetaPhlAn | rdata, biom-format, zipfile |
| `datasets` | GlobalPatterns, enterotype, curatedMetagenomicData, MGnify | pooch caching |
| `pp` | Filter, prune, subset, rarefy, transforms (relative, CLR, PhILR), agglomerate | anndata, NumPy, SciPy |
| `tl` | Alpha and beta diversity, ordination, PERMANOVA, clustering (DMM) | scikit-bio |
| `fn` | Function hierarchies (KO → module → pathway), stratified contributions, functional redundancy | KEGG/MetaCyc mapping files |
| `da` | Differential abundance methods and the consensus runner | scikit-bio, statsmodels, optional rpy2 |
| `ml` | scikit-learn transformers, embedding plugins, PyTorch loaders | scikit-learn, torch (optional) |
| `pl` | Bar, heatmap, ordination, richness, tree plus heatmap | matplotlib, plotnine (optional) |

Heavy dependencies (torch, rpy2, plotnine) are optional extras, so `pip install <name>` stays light.

## Performance: Python first, compiled code last

Start in pure Python: the heaviest microbiome math already runs as compiled code in the libraries you will call, so Rust or C should enter only for a proven hotspot that no dependency covers, and only after Numba has been tried.

### What is already compiled for you

- **scikit-bio** ships Cython implementations and optional Numba kernels (`engine="numba"`) for PERMANOVA, Mantel, PERMDISP, distance-matrix centering and UniFrac. The [Numba PR](https://github.com/scikit-bio/scikit-bio/pull/2557) reports 7x to 13x speedups for PERMANOVA and up to \~9x for weighted UniFrac.
- **[scikit-bio-binaries](https://github.com/scikit-bio/scikit-bio-binaries)** is optional C++ for PCoA and PERMANOVA, with CUDA and AMD HIP GPU builds, loaded through ctypes.
- **[unifrac](https://github.com/biocore/unifrac)** (Striped UniFrac, C/C++) covers unweighted, weighted, generalized and variance-adjusted UniFrac plus Faith's PD, with automatic GPU use on Linux. It reads BIOM files and Newick trees from disk, so you need a thin adapter from TreeData.
- **NumPy, SciPy sparse, scikit-learn and polars** are already compiled. Most `pp` functions (filter, transform, glom) are vectorized sparse operations and will never need anything else.

### The escalation ladder

Each step is taken only when the previous one fails a benchmark on realistic data (for example 5,000 samples × 50,000 features).

1. **Vectorize** with NumPy/SciPy sparse. Profile with py-spy or scalene before assuming anything is slow.
2. **Delegate** to scikit-bio, unifrac or scikit-bio-binaries. Expose their engine or GPU switch; do not reimplement.
3. **Numba** for your own kernels: permutation loops, per-sample loops over CSR `indptr/indices/data` arrays, tree traversals over flattened node arrays. It stays in Python syntax, is an optional dependency, and mirrors scikit-bio's `engine=` pattern.
4. **Rust (PyO3 + maturin)** for a sustained hotspot Numba handles badly: pointer-heavy tree algorithms, parallel parsers for large HUMAnN stratified or BIOM files, string-heavy taxonomy parsing. [SnapATAC2](https://scverse.org/SnapATAC2/) is the scverse precedent: a Python/Rust package built on the [anndata Rust crate](https://crates.io/crates/anndata).
5. **New C/C++: no.** Use C/C++ only through existing libraries. For new native code, Rust is safer and maturin builds wheels far more easily than a C toolchain.
6. **GPU** later, through the array API support in scikit-bio 0.7.3 and CuPy, not as custom kernels.

### Rules for any compiled kernel

- It is a drop-in behind an existing public function: same signature, chosen by `engine="python" | "numba" | "rust"`.
- The pure-Python version stays as the reference and test oracle; both must match to a stated tolerance.
- Merge only with an asv benchmark showing at least 5x on the realistic dataset, and a note in the docs.
- Rust lives in one crate inside the repo (`rust/`), built by maturin; cibuildwheel produces Linux, macOS and Windows wheels so users never need a compiler.

### Likely hotspots and planned approach

| Task | First approach | Escalate to |
| --- | --- | --- |
| Beta diversity, UniFrac, Faith's PD | unifrac / scikit-bio | GPU flags |
| PCoA, PERMANOVA, Mantel | scikit-bio (Numba engine), scikit-bio-binaries | nothing needed |
| Rarefaction curves (many depths × iterations) | NumPy multivariate hypergeometric | Numba |
| Tax/func glom on large sparse tables | SciPy sparse indicator-matrix product | nothing needed |
| DA consensus across methods | joblib over methods | nothing needed |
| Permutation-based DA tests | NumPy | Numba |
| PhILR, tree-based transforms | Python over flattened tree arrays | Numba, then Rust |
| Parsing large HUMAnN stratified or BIOM files | polars / h5py | Rust |

## Roadmap

Ship v0.1 as a credible phyloseq replacement, then add the three differentiators in order; each release has an exit gate that must pass before the next one starts.

| Release | Scope | Exit gate |
| --- | --- | --- |
| 0.0.1 | Name reserved on PyPI and GitHub, repo skeleton, CI, docs site, rules in CONTRIBUTING and agent instructions | CI green, docs build, placeholder published |
| 0.1 Core | TreeData object conventions; `io` for phyloseq .rds, BIOM, .qza, DADA2; `pp` filter, prune, subset, rarefy, relative, glom; `tl` alpha, beta, UniFrac, PCoA, NMDS, PERMANOVA; `pl` bar, richness, ordination, heatmap; GlobalPatterns and enterotype datasets | Reproduces the phyloseq tutorial end to end; golden tests match R within tolerance; migration table covers the top 30 phyloseq functions |
| 0.2 Function | `io` for HUMAnN 4, PICRUSt2, MetaPhlAn; `fn` hierarchies (KO → module → pathway, EC → MetaCyc), `func_glom`, stratified contribution links, functional redundancy | One worked HUMAnN tutorial on a public cohort; results cross-checked against HUMAnN utility scripts |
| 0.3 Stats | `da` wrappers (ANCOM-BC, ALDEx2-like, LinDA, MaAsLin 3 via port or rpy2) and consensus report; CLR, PhILR | Consensus report on a benchmark dataset; per-method agreement with the R reference |
| 0.4 ML and multi-omics | MuData modalities, mmvec wrapper; `ml` scikit-learn transformers, embedding plugin interface (MGM, BiomeGPT), PyTorch loader | Leak-free CV example; one foundation model plugged in end to end |
| 0.5+ | Time series (stability, trajectories, DMM, community state types), single-cell and spatial bridge, MCP server, federated meta-analysis | Driven by users and issues |

Rough effort for one person working part time: 0.1 in 6 to 8 weeks, 0.2 and 0.3 in about 4 weeks each. Minimal-code discipline and library reuse are what keep these numbers realistic.

## Testing and validation

Correctness is proven against R, not asserted: every function with an R equivalent has a golden test generated from phyloseq, mia or the method's reference package.

- **Golden files from R.** A script under `tests/r/` (run in a pinned R container, not in normal CI) exports phyloseq/mia results on the bundled datasets to parquet. Python tests compare against them with explicit tolerances. Regenerate only when the R version is bumped.
- **One test file per module, small tests per function.** Tests mirror the function-level design: each public function gets a happy path, an edge case (empty sample, all-zero feature, missing taxonomy rank) and, when relevant, a golden comparison.
- **Property tests with Hypothesis** for invariants: glom preserves column sums, relative abundance sums to 1, CLR rows sum to 0, distances are symmetric with a zero diagonal.
- **Engine parity.** When a Numba or Rust engine exists, the same test runs across all engines.
- **Benchmarks with asv** on a synthetic large dataset (5,000 × 50,000, sparse), tracked across commits, as scikit-bio does.
- **Docs are tests.** Every docstring example runs under doctest and every tutorial notebook runs in CI, so the docs cannot go stale.
- **Coverage target** of 90% on public functions, measured but not gamed.

## Documentation

Because code comments are kept to a minimum, the docs site carries the explanation, and it is organized so each function's docstring can link to one page that explains it in depth.

| Section | Content | Format |
| --- | --- | --- |
| Getting started | Install, load a dataset, first ordination in 10 lines | Notebook |
| Coming from R | Function-by-function table: phyloseq and mia → this package, with orientation and naming differences | Reference page, generated from the `R equivalent` docstring lines |
| User guide | One page per concept: data model, aggregation, transforms, diversity, DA, function hierarchies, ML | Markdown with small runnable snippets |
| Tutorials | Full analyses on public data: 16S case-control, shotgun plus HUMAnN, multi-omics | Executed notebooks |
| API reference | Auto-generated from docstrings | Sphinx autodoc or mkdocstrings |
| Method notes | The math, assumptions and caveats behind each statistical method, with citations | Markdown |
| Design decisions | Short ADRs: why TreeData, why samples as rows, why no custom class | Markdown |
| Contributing | The engineering rules above, how to add a function, how to add an engine | Markdown |

- Tooling: Sphinx with the scverse cookiecutter theme (or MkDocs Material if you prefer Markdown everywhere), myst-nb to execute notebooks, hosted on Read the Docs.
- The "Coming from R" table is generated, not written by hand, so it can never drift from the code.

## Tooling, packaging, CI and release

Start from the scverse cookiecutter template so the repo looks and builds like every other scverse package, which also smooths a later application to become a scverse ecosystem package.

| Area | Choice | Why |
| --- | --- | --- |
| Template | cookiecutter-scverse | Same layout, docs and CI conventions as scanpy and friends |
| Environment | uv | Fast, lockfile, one tool for venvs and installs |
| Build backend | hatchling now; maturin once a Rust crate exists | Pure-Python wheels until compiled code is justified |
| Lint and format | ruff (lint, format, complexity, docstring rules) | One fast tool enforces the minimal-code rules |
| Types | pyright or mypy on public functions | Types replace most comments |
| Tests | pytest, hypothesis, pytest-doctestplus | Unit, property and docstring tests |
| Benchmarks | asv | Performance gate for compiled engines |
| CI | GitHub Actions: test matrix on Python 3.11 to 3.14, Linux, macOS, Windows; docs build; notebooks executed | Catches platform breakage early |
| Release | Tag triggers build and PyPI Trusted Publishing; conda-forge or bioconda recipe after 0.1 | No tokens to leak; conda users covered |
| Wheels (with Rust) | cibuildwheel + maturin for manylinux, macOS arm64/x86, Windows | Users never need a compiler |
| Versioning | SemVer from 0.x, CHANGELOG following Keep a Changelog | Clear upgrade expectations |
| License | BSD-3-Clause | Matches scikit-bio and scverse, friendly to pharma use |
| Agent instructions | CLAUDE.md or AGENTS.md with the engineering rules | AI-assisted code follows the same discipline |

## Risks and open questions

The biggest risk is scope: the feature list is large, so the minimal-code rule and strict release gates are what keep this shippable by one person.

| Risk | Mitigation |
| --- | --- |
| pyloseq gains traction first as "phyloseq for Python" | Position on mia, scverse and function, not on phyloseq parity; offer a pyloseq importer |
| Scope creep across ten feature areas | Release gates; nothing from 0.4+ starts before 0.3 ships |
| DA methods only exist in R (MaAsLin 3, ANCOM-BC2) | Optional rpy2 bridge first; native ports only for the most used, validated against R |
| TreeData is young (v0.3.x) and could change its API | Pin a range, wrap tree access in 2 or 3 helper functions so a change touches one place |
| Rust raises the contributor bar | Keep Rust optional and last; Python engine always exists |
| KEGG licensing limits redistributing mapping files | Download on first use from user-provided or open sources (MetaCyc, eggNOG); do not bundle |
| Side project overlaps with the day job (federated, pharma) | Check employment IP and side-project clauses before building the federated module in the open |

Open questions:

- [ ] Final name (shortlist: commensal, holobiont, taxonomia, microbiota) and GitHub org or personal account?
- [ ] Sphinx (scverse standard) or MkDocs Material for docs?
- [ ] Aim for official scverse ecosystem listing from 0.1, or later?
- [ ] Port MaAsLin 3 natively, or rpy2 only?
- [ ] Which public cohort anchors the tutorials (HMP2/IBDMDB is a strong candidate for taxa plus function plus metabolites)?

## Sources

- [pyloseq on PyPI](https://pypi.org/project/pyloseq/)
- [scikit-bio releases](https://github.com/scikit-bio/scikit-bio/releases)
- [scikit-bio Numba kernels PR #2557](https://github.com/scikit-bio/scikit-bio/pull/2557)
- [scikit-bio-binaries](https://github.com/scikit-bio/scikit-bio-binaries)
- [unifrac (Striped UniFrac)](https://github.com/biocore/unifrac)
- [treedata on PyPI](https://pypi.org/p/treedata)
- [SnapATAC2 (Python/Rust scverse package)](https://scverse.org/SnapATAC2/)
- [anndata Rust crate](https://crates.io/crates/anndata)
- [miaverse overview, OMA book](https://microbiome.github.io/OMA/docs/devel/pages/miaverse.html)
- [HUMAnN](https://github.com/biobakery/humann)
- [MaAsLin 3, Nature Methods](https://www.nature.com/articles/s41592-025-02923-9)
- [MGM foundation model](https://pmc.ncbi.nlm.nih.gov/articles/PMC13116254/)
- [BiomeGPT preprint](https://www.biorxiv.org/content/10.64898/2026.01.05.697599v1.full)
