# Performance

Baselines measured with [asv](https://asv.readthedocs.io/) on a synthetic table: 5,000 samples x
50,000 features, 2% non-zero counts from 1 to 99, a taxonomy of 1,000 genera in 5 phyla and a
balanced tree of 100,006 nodes, all from seed 0 (`benchmarks/benchmarks/_data.py`). They are a
record, not a gate: release 0.1 sets no speed target.

Measured on commit `df31dbd`, 2026-10-03: a 16-thread laptop (11th Gen Intel Core i7-11800H, 2.30
GHz) with 62 GB of RAM, Linux, Python 3.13, in the development environment (`asv run --python=same`).

Each figure comes from a single asv run of the whole suite on one machine, with no measured
variance. The load average (1, 5 and 15 minutes) was 2.02, 2.04, 1.82 before the run and 2.73,
2.54, 2.24 after it. The `peakmem_` figures include the memory of the loaded table itself.

| Benchmark | Result |
|---|---|
| `pp.relative` | 1.02 s |
| `pp.tax_glom` to genus | 3.97 s |
| `pp.rarefy` to the smallest depth | 6.29 s |
| `tl.alpha`, the four default metrics | 1.76 s |
| `tl.alpha`, the four default metrics, peak memory | 471 MB |
| `tl.alpha`, `faith_pd` | 46.6 s |
| `tl.beta`, Bray-Curtis, 500 samples | 4.13 s |
| `tl.beta`, Bray-Curtis, 1,000 samples | 17.8 s |
| `tl.beta`, Bray-Curtis, 2,000 samples | 66.7 s |
| `tl.beta`, Bray-Curtis, 5,000 samples, peak memory | 2.73 GB |

- `tl.beta` compares every pair of samples, so its time grows with the square of the samples:
  at 5,000 samples it was run once, for peak memory only, and not timed: asv recorded 7.05 minutes
  for that whole run, and extrapolating the 2,000-sample figure quadratically (66.7 s x 6.25) gives
  about 7 minutes, an estimate and not a measurement. Its memory is one dense copy of `X` (8 bytes x 5,000 x 50,000 =
  2 GB) plus the distances, as its docstring says.
- `tl.alpha` densifies at most 2**20 values at a time, so its peak stays far below one dense
  copy of `X`.
- `faith_pd` spends almost all its time in scikit-bio, which re-indexes the 100,006-node tree for
  each chunk of 20 samples (2**20 values, so 20 samples at 50,000 features); about 97% of its time
  is that re-indexing, and converting the tree once takes 0.46 s. The optimization is deferred until
  a profile-driven task (rules.md R10.1).

## Function tables

Measured on commit `2e9d268`, 2026-10-05, on the same laptop and environment, with
`asv run --python=same --bench "^fn\."` (the three classes took about 60 s). The load average was
2.48, 1.84, 1.95 before the run and 1.89, 1.78, 1.92 after it. Release 0.2 sets no speed target
either.

The function table is shaped like HMP2's pathway table: 1,600 samples, 500 functions with 43
strata each (22,000 rows), 50% of the community values and 7% of the stratified values non-zero.
The hierarchy puts every function in two of 50 groups. Functional redundancy runs on 100 samples x
2,000 taxa (5% non-zero) and 2,500 genes (70% zero copy numbers). All of it comes from seed 0
(`benchmarks/benchmarks/_data.py`).

| Benchmark | Result |
|---|---|
| `fn.func_glom`, 500 community rows to 50 groups | 10.5 ms |
| `fn.func_glom`, 21,500 stratified rows to 50 groups per taxon | 57.0 ms |
| `io.read_humann`, 22,000 rows x 1,600 samples | 2.04 s |
| `io.read_humann`, peak memory | 1.04 GB |
| `fn.functional_redundancy`, 2,000 taxa | 5.28 s |
| `fn.functional_redundancy`, 2,000 taxa, peak memory | 446 MB |

- `read_humann`'s peak is about 3.7 times the dense array its docstring names (8 bytes x 22,000 x
  1,600 = 282 MB). The figure also holds the Python process and pandas' parse of the text; no
  profile was taken, as nothing is being optimised (rules.md R10.1).
- `functional_redundancy` compares every pair of taxa, so its time grows with the square of the
  taxa times the genes, and its memory with the square of the taxa.

## Transforms and differential abundance

Measured on commit `51e0475`, 2026-10-07, on the same laptop and environment, with
`asv run --python=same --bench "^(da\.|pp\.Philr)"` (the two classes took about 55 s). The load
average was 1.66, 1.30, 1.21 before the run and 4.64, 2.23, 1.54 after it. Release 0.3 sets no
speed target either.

PhILR runs on the synthetic table above built with 1,000 samples (50,000 features, 2% non-zero,
the same balanced tree). The two differential abundance methods that run without R run on 2,000
samples in two groups of 1,000 and 10,000 features, 30% non-zero counts from 1 to 99, the size of
a large genus-level study. All of it comes from seed 0 (`benchmarks/benchmarks/_data.py`). ALDEx2
and MaAsLin 3 are not benchmarked: their time is R's.

| Benchmark | Result |
|---|---|
| `pp.philr`, 1,000 x 50,000 | 3.12 s |
| `pp.philr`, peak memory | 2.14 GB |
| `da.linda`, 2,000 x 10,000 | 425 ms |
| `da.linda`, peak memory | 1.15 GB |
| `da.ancombc2`, 2,000 x 10,000 | 3.62 s |
| `da.ancombc2`, peak memory | 1.07 GB |

- Each of the three densifies `X` once (rules.md R6.2): 400 MB for PhILR, 160 MB for the
  differential abundance table. The peaks also hold the Python process (418 MB with PhILR's table
  loaded, 345 MB with the differential abundance table).
- Measured outside asv on the same tables, the call itself added 1.72 GB for `pp.philr` (4.3 times
  the dense table, plus scikit-bio's sparse basis, as its docstring says), 822 MB for `da.linda`
  (5.1 times, as its docstring says) and 726 MB for `da.ancombc2` (4.5 times; its docstring gives
  4.5 times here and 7.2 to 7.5 times on the 400 x 1,000 and 200 x 2,000 tables it was first
  measured on, so the range depends on the table's shape). Nothing is being optimised
  (rules.md R10.1).

## Running the benchmarks

```bash
uv sync --group dev
cd benchmarks
# asv keeps its machine description in ~/.asv-machine.json; HOME points it at the git-ignored .asv/.
uv run --group dev env HOME="$PWD/../.asv" asv machine --yes
uv run --group dev env HOME="$PWD/../.asv" asv run --python=same --set-commit-hash "$(git rev-parse HEAD)"
uv run --group dev env HOME="$PWD/../.asv" asv show "$(git rev-parse HEAD)"
```

`--python=same` runs in the current environment, and `--set-commit-hash` keeps the results, in
`.asv/results`. Add `--bench "^fn\."` (or `"^(da\.|pp\.Philr)"`) to `asv run` and `asv show` to
run or show only the function (or transform and differential abundance) benchmarks. The whole suite
took about 13 minutes on this run and needs about 3 GB of free memory.
`uv run --group dev asv check --python=same` imports the suite without running it; CI runs it.
