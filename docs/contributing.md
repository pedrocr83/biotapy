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

- `dev`: ruff, mypy, import-linter, prek - linting and type-checking.
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
the size/complexity limits in rules.md R5), `mypy --strict` on `src/biotapy`,
import-linter (the module-layer contract), and pyproject-fmt.

## Running tests

```bash
uv run --group test pytest
```

Network and golden tests are excluded by default (`[tool.pytest]` in
`pyproject.toml` sets `-m "not network and not r"`). Run them explicitly:

```bash
uv run --group test pytest -m "network or golden"
```

These tests download `bt.datasets.global_patterns()`, `bt.datasets.enterotype()`
and `bt.datasets.esophagus()` through [pooch](https://www.fatiando.org/pooch/)
and compare `pp.relative` and `pp.tax_glom` against phyloseq's output on that data. Set
`BIOTAPY_DATA_DIR` to point the pooch cache somewhere other than the default
per-user cache directory - CI caches it across runs the same way:

```bash
BIOTAPY_DATA_DIR=.pooch uv run --group test pytest -m "network or golden"
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
