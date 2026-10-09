# Contributing guide

All contributions follow
[rules.md](https://github.com/pedrocr83/biotapy/blob/master/rules.md), the
binding development rules for every change, human or agent. Background and
reasoning for those rules live in the
[knowledge bundle](https://github.com/pedrocr83/biotapy/tree/master/.knowledge)
(decisions, contracts, roadmap, playbooks); new public functions follow the
[add-a-function playbook](https://github.com/pedrocr83/biotapy/blob/master/.knowledge/playbooks/add-a-function.md).

We assume you are already familiar with git and with making pull requests on
GitHub. For the absolute basics, see the
[pyopensci tutorials](https://www.pyopensci.org/learn.html) or the
[scientific Python tutorials](https://learn.scientific-python.org/development/tutorials/).

## Setting up a development environment

biotapy uses [uv](https://docs.astral.sh/uv/) to manage its environment. Sync
the groups you need:

```bash
uv sync --group dev --group test --group doc
```

- `dev`: ruff, mypy, import-linter, prek, asv - linting, type-checking and
  benchmarks.
- `test`: pytest, hypothesis, coverage.
- `doc`: sphinx, myst-nb, sphinx-book-theme and the other packages that build
  this site.

The `.venv` directory `uv sync` creates is typically auto-discovered by IDEs
such as VS Code.

## Linting and type-checking

```bash
uvx prek run --all-files
```

This runs the same hooks as CI's `lint` job: ruff lint and format (including
the size/complexity limits in rules.md R5), `mypy --strict` on `src/biotapy`
and `docs/extensions`, import-linter (the module-layer contract), and
pyproject-fmt.

## Running tests

```bash
uv run --group test pytest
```

Network and golden tests are excluded by default (`[tool.pytest]` in
`pyproject.toml` sets `-m "not network and not r and not torch"`). Run them explicitly:

```bash
uv run --group test pytest -m "network or golden"
```

These tests download `bt.datasets.global_patterns()`, `bt.datasets.enterotype()`
and `bt.datasets.esophagus()` through [pooch](https://www.fatiando.org/pooch/)
and compare biotapy with R on that data: `pp.relative`, `pp.clr`, `pp.tax_glom`, filtering, rarefaction
(its invariants), alpha and beta diversity, UniFrac, PCoA, NMDS and PERMANOVA against phyloseq,
vegan, ape and picante, `pp.philr` against philr, `da.linda` against MicrobiomeStat and `da.ancombc2` against ANCOMBC. Set
`BIOTAPY_DATA_DIR` to point the pooch cache somewhere other than the default
per-user cache directory - CI caches it across runs the same way:

```bash
BIOTAPY_DATA_DIR=.pooch uv run --group test pytest -m "network or golden"
```

### R bridge tests

Tests that call R through rpy2 (`bt.da.aldex2`, `bt.da.maaslin3`) carry the marker `r` and are excluded from the
runs above. They need R with the packages `tests/r/Dockerfile` installs and the `r` extra; those
that read GlobalPatterns also need the pooch cache:

```bash
BIOTAPY_DATA_DIR=.pooch uv run --group test --extra r pytest -m r
```

CI runs them in the `r-bridge` job, with R 4.5.3 and the Bioconductor 3.22 packages the golden image
pins.

### PyTorch tests

Tests that need PyTorch (`bt.ml.to_torch`, and its docstring example) carry the marker `torch` and are
excluded from the runs above. The `torch` extra installs torch's CPU wheel from PyTorch's index
(`[tool.uv]` in `pyproject.toml`):

```bash
uv run --group test --extra torch pytest -m torch
```

### Regenerating the R golden files

The golden CSVs under `tests/golden/` and the R-written fixtures under
`tests/data/phyloseq/` and `tests/data/dada2/` are produced by a pinned R
container, not by pytest, and are never hand-edited. See the
[regenerate-golden-files playbook](https://github.com/pedrocr83/biotapy/blob/master/.knowledge/playbooks/regenerate-golden-files.md)
for the steps and its bit-identical check.

## Building the docs locally

```bash
uv run --group doc sphinx-build -W -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html`. Read the Docs builds this same site
with `uvx hatch run docs:build` (see `.readthedocs.yaml` and the `docs` hatch
environment in `pyproject.toml`), which wraps the equivalent `sphinx-build`
invocation.

Notebooks run on every build (`nb_execution_mode = "cache"`), but a notebook is
re-executed only when its own content changes. After a code change, run
`uvx hatch run docs:clean` first, which deletes every git-ignored file under
`docs/` (`_build`, with the jupyter cache in it, and `generated`), or delete
`docs/_build/.jupyter_cache`, to force re-execution locally. CI and Read the
Docs always start clean.

If you refer to objects from another package, add an entry to
`intersphinx_mapping` in `docs/conf.py` so Sphinx can link to it. If the build
fails over a link outside your control, add an exception to `nitpick_ignore`
in the same file.

## Commit conventions

Commits follow [Conventional Commits](https://www.conventionalcommits.org/)
(`feat:`, `fix:`, `refactor:`, `docs:`, `test:`, `build:`, `ci:`, `chore:`),
one logical change per commit (rules.md R13.1).

## Publishing a release

Releases follow the
[cut-a-release playbook](https://github.com/pedrocr83/biotapy/blob/master/.knowledge/playbooks/cut-a-release.md):
bump the version, move the `[Unreleased]` changelog entry, tag, publish a
GitHub release, and let `release.yaml` upload to PyPI through trusted
publishing.
