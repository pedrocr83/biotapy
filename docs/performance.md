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
`.asv/results`. The whole suite took about 17 minutes on this run and needs about 3 GB of free memory.
`uv run --group dev asv check --python=same` imports the suite without running it; CI runs it.
