# Performance

Baselines measured with [asv](https://asv.readthedocs.io/) on a synthetic table: 5,000 samples x
50,000 features, 2% non-zero counts from 1 to 99, a taxonomy of 1,000 genera in 5 phyla and a
balanced tree of 100,006 nodes, all from seed 0 (`benchmarks/benchmarks/_data.py`). They are a
record, not a gate: release 0.1 sets no speed target.

Measured on commit `7bed7ef`, 2026-10-03: a 16-thread laptop (11th Gen Intel Core i7-11800H, 2.30
GHz) with 62 GB of RAM, Linux, Python 3.13, in the development environment (`asv run --python=same`).

| Benchmark | Result |
|---|---|
| `pp.relative` | 0.88 s |
| `pp.tax_glom` to genus | 3.41 s |
| `pp.rarefy` to the smallest depth | 5.51 s |
| `tl.alpha`, the four default metrics | 1.36 s |
| `tl.alpha`, the four default metrics, peak memory | 471 MB |
| `tl.alpha`, `faith_pd` | 48.1 s |
| `tl.beta`, Bray-Curtis, 500 samples | 8.34 s |
| `tl.beta`, Bray-Curtis, 1,000 samples | 31.1 s |
| `tl.beta`, Bray-Curtis, 2,000 samples | 88.4 s |
| `tl.beta`, Bray-Curtis, 5,000 samples, peak memory | 2.73 GB |

- `tl.beta` compares every pair of samples, so its time grows with the square of the samples:
  about 6 minutes at 5,000. Its memory is one dense copy of `X` (8 bytes x 5,000 x 50,000 =
  2 GB) plus the distances, as its docstring says.
- `tl.alpha` densifies at most 2**20 values at a time, so its peak stays far below one dense
  copy of `X`.
- `faith_pd` spends almost all its time in scikit-bio, which re-indexes the 100,006-node tree for
  each chunk of 20 samples; converting the tree once takes 0.46 s.

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
`.asv/results`. The whole suite takes about 11 minutes and needs about 3 GB of free memory.
`uv run --group dev asv check --python=same` imports the suite without running it; CI runs it.
