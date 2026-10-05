---
type: Phase
title: Phase 0 - Foundation (0.0.1)
description: Repo skeleton from cookiecutter-scverse, tooling that mechanically enforces rules.md, first _core helpers, OKF conformance test, CI, docs site, placeholder release on PyPI.
tags: [roadmap, tooling, ci, release]
status: stable
release: "0.0.1"
phase_state: done
effort: ~1 week part-time
depends_on: []
paths: ["pyproject.toml", ".pre-commit-config.yaml", ".github/**", "docs/**", "src/biotapy/__init__.py", "src/biotapy/_core/**", "tests/**", "scripts/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T11:55:00Z }
commit: ef2fe8b
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: cookiecutter
    resource: https://github.com/scverse/cookiecutter-scverse
    title: cookiecutter-scverse v0.8.0 (2026-07-17)
  - id: ruff
    resource: https://docs.astral.sh/ruff/settings/
    title: ruff settings (0.16.9)
  - id: import-linter
    resource: https://github.com/seddonym/import-linter/blob/main/docs/contract_types/layers.md
    title: import-linter layers contract (2.15)
---

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans. Steps use `- [ ]` checkboxes;
> tick them in this file as you go (rules.md R12.4).

**Goal:** a repository where every rule in `rules.md` that a machine can check
is checked in CI, and the name `biotapy` is reserved on PyPI.[^spec]

**Architecture:** generate from cookiecutter-scverse v0.8.0 and keep what it
generates (hatchling, hatch envs with the uv installer, PEP 735 dependency
groups, prek hooks, Sphinx book theme, trusted-publishing release).[^cookiecutter]
Then tighten it: ruff size limits, banned tree imports, mypy strict,
import-linter layers, an OKF conformance test.

**Tech stack:** Python >= 3.12 · cruft · hatchling · uv · prek · ruff 0.16 ·
mypy · import-linter 2.15 · pytest + hypothesis · Sphinx + myst-nb · GitHub
Actions · PyPI trusted publishing.

**Spec:** [plan.md](../../plan.md). Rules: [rules.md](../../rules.md).

# Global constraints
- Python floor **3.12**: treedata 0.3.1, mudata 0.4.1 and the template require it. The spec's 3.11 is dropped.
- CI matrix: Python 3.12, 3.13, 3.14 x ubuntu, macos, windows.
- License: BSD-3-Clause. Default branch: `master`. Remote: `git@github.com:pedrocr83/biotapy.git`.
- No dependency is added without user approval (R9.1). Task 0.1 asks once for this phase's list.
- Push, tag, PyPI, Read the Docs and GitHub settings each need explicit approval (R13.3).

# Dependencies to approve (task 0.1)
| Group | Package | Reason |
|---|---|---|
| runtime | numpy | `_core.as_generator` (task 0.3) |
| test | hypothesis | property tests (R11.2) |
| test | pyyaml | OKF conformance test parses frontmatter (task 0.5) |
| dev | import-linter | enforces [module-boundaries](/contracts/module-boundaries.md) (task 0.4) |

Everything else comes with the template (pytest, coverage, mypy, sphinx stack), approved by approving the template in task 0.1.

# Review focus
1. **Windows paths and encodings** in the knowledge test: use `pathlib` and `encoding="utf-8"`; runs in the windows matrix cells (task 0.5).
2. **A concept without frontmatter** must fail the test with its file name, not a generic error (task 0.5, step 3).
3. **import-linter on a package with no subpackages yet**: layers are optional `(pp)`, so the contract passes now and binds later (task 0.4, step 5).
4. **mypy strict and untyped third-party packages**: add overrides only for modules mypy reports; no stub packages without approval (task 0.4, step 3).
5. **An optional import leaking to module level**: CI imports `biotapy` without dev or extra groups (task 0.6, step 2).

---

### Task 0.1: Confirm decisions and commit the knowledge bundle

**Files:** modify `.knowledge/decisions/*.md` (status, verified), `.knowledge/log.md`.

- [x] **Step 1: Ask the user to confirm**, in one message:
  (a) name `biotapy` on PyPI and `pedrocr83/biotapy` ([package-name-biotapy](/decisions/package-name-biotapy.md));
  (b) Python >= 3.12;
  (c) [pure-by-default](/decisions/pure-by-default.md);
  (d) [docs-okf-and-sphinx](/decisions/docs-okf-and-sphinx.md);
  (e) [r-bridge-before-ports](/decisions/r-bridge-before-ports.md);
  (f) the dependency table above plus the template itself;
  (g) the author email to put in `pyproject.toml`.
- [x] **Step 2:** for each confirmed `draft` decision set `status: stable` and add
  `verified: { by: human:pedrocr83, at: <UTC now> }`. A rejected decision is
  rewritten before any other task starts.
- [x] **Step 3:** add a dated entry to `.knowledge/log.md`.
- [x] **Step 4: Commit**
  ```bash
  git add CLAUDE.md rules.md plan.md .knowledge
  git commit -m "docs: add rules, OKF knowledge bundle and phased roadmap"
  ```

### Task 0.2: Generate the skeleton

**Files:** create everything the template generates; delete its example code.
**Interfaces:** produces `src/biotapy/__init__.py` exposing only `__version__`.

- [x] **Step 1: Generate outside the repo** (cruft creates a new directory):
  ```bash
  cd "$(mktemp -d)"
  uvx --with prek cruft create https://github.com/scverse/cookiecutter-scverse --checkout v0.8.0
  ```
  Answers: project and package name `biotapy`; description "mia-style
  microbiome toolkit for Python on AnnData/TreeData"; author "Pedro Ribeiro";
  email from task 0.1; GitHub user `pedrocr83`; repo `biotapy`; license
  BSD 3-Clause. Keep the prompt names and answers for the commit body.
- [x] **Step 2: Copy into the repo**, keeping this repo's own files:
  ```bash
  rsync -a --exclude .git --exclude CLAUDE.md --exclude rules.md \
        --exclude plan.md --exclude .knowledge \
        ./biotapy/ /home/pedro/Desktop/Business/AI/biotapy/
  cd /home/pedro/Desktop/Business/AI/biotapy && git status
  ```
- [x] **Step 3: Delete the template's example code** (R4.8, no placeholders):
  ```bash
  rm -r src/biotapy/pp src/biotapy/tl src/biotapy/pl tests/test_basic.py \
        docs/notebooks/example.ipynb docs/template_usage.md
  ```
  Read `docs/index.md` and `docs/api.md` and remove their references to the
  deleted files. Replace `src/biotapy/__init__.py` with:
  ```python
  from importlib.metadata import version

  __all__ = ["__version__"]

  __version__ = version("biotapy")
  ```
- [x] **Step 4: Python matrix.** In `pyproject.toml` confirm
  `requires-python = ">=3.12"` and set the hatch-test matrix Python list to
  `["3.12", "3.13", "3.14"]`.
- [x] **Step 5: Verify**
  ```bash
  uv sync --all-groups
  uv run python -c "import biotapy; print(biotapy.__version__)"
  ```
  Expected: a version string, no traceback.
- [x] **Step 6: Record real commands.** Read the generated `pyproject.toml`,
  `.pre-commit-config.yaml` and `.github/workflows/*.yaml`. If group or env
  names differ from rules.md R14 (`test`, `doc`, `prek`), fix R14 now.
- [x] **Step 7: Commit**
  ```bash
  git add -A
  git commit -m "build: generate skeleton from cookiecutter-scverse v0.8.0"
  ```

### Task 0.3: `_core` kernel v0 (seeded RNG, optional imports)

**Files:** create `src/biotapy/_core/__init__.py`, `src/biotapy/_core/_rng.py`,
`src/biotapy/_core/_optional.py`, `tests/core/test_rng.py`,
`tests/core/test_optional.py`; modify `pyproject.toml` (add `numpy`).
**Interfaces (produces):**
- `as_generator(seed: int | np.random.Generator | None) -> np.random.Generator`
- `import_optional(name: str, *, extra: str) -> types.ModuleType`

- [x] **Step 1: Failing tests**
  ```python
  # tests/core/test_rng.py
  import numpy as np
  import pytest

  from biotapy._core import as_generator


  def test_int_seed_is_reproducible():
      assert as_generator(7).integers(10**9) == as_generator(7).integers(10**9)


  def test_generator_is_passed_through():
      rng = np.random.default_rng(1)
      assert as_generator(rng) is rng


  def test_none_gives_a_generator():
      assert isinstance(as_generator(None), np.random.Generator)


  def test_legacy_randomstate_is_rejected():
      with pytest.raises(TypeError, match="seed"):
          as_generator(np.random.RandomState(0))  # type: ignore[arg-type]
  ```
  ```python
  # tests/core/test_optional.py
  import pytest

  from biotapy._core import import_optional


  def test_imports_an_installed_module():
      assert import_optional("json", extra="none").dumps([]) == "[]"


  def test_missing_module_names_the_extra():
      with pytest.raises(ImportError, match=r"pip install 'biotapy\[torch\]'"):
          import_optional("biotapy_no_such_module", extra="torch")
  ```
- [x] **Step 2: Run, expect failure**
  `uv run --group test pytest tests/core -q` -> `ModuleNotFoundError: No module named 'biotapy._core'`.
- [x] **Step 3: Implement**
  ```python
  # src/biotapy/_core/_rng.py
  """Single entry point for randomness (rules.md R3.4)."""

  import numpy as np


  def as_generator(seed: int | np.random.Generator | None) -> np.random.Generator:
      """Return a NumPy Generator for ``seed`` without touching global state."""
      if isinstance(seed, np.random.Generator):
          return seed
      if seed is None or isinstance(seed, int | np.integer):
          return np.random.default_rng(seed)
      msg = f"seed must be an int, a numpy Generator or None, got {type(seed).__name__}"
      raise TypeError(msg)
  ```
  ```python
  # src/biotapy/_core/_optional.py
  """Lazy imports for optional extras (decisions/optional-heavy-dependencies)."""

  import importlib
  from types import ModuleType


  def import_optional(name: str, *, extra: str) -> ModuleType:
      """Import ``name`` or raise an ImportError that names the extra to install."""
      try:
          return importlib.import_module(name)
      except ImportError as err:
          msg = f"{name} is required here. Install it with: pip install 'biotapy[{extra}]'"
          raise ImportError(msg) from err
  ```
  ```python
  # src/biotapy/_core/__init__.py
  """Private kernel shared by biotapy subpackages (contracts/module-boundaries)."""

  from ._optional import import_optional
  from ._rng import as_generator

  __all__ = ["as_generator", "import_optional"]
  ```
  Add `"numpy"` to `[project] dependencies`.
- [x] **Step 4: Run, expect pass** - `uv run --group test pytest tests/core -q` -> 6 passed.
- [x] **Step 5: Commit** - `git add -A && git commit -m "feat(core): add seeded generator and optional-import helpers"`

### Task 0.4: Tooling that enforces rules.md

**Files:** modify `pyproject.toml`, `.pre-commit-config.yaml`.

- [x] **Step 1: ruff limits (R5) and banned imports (R4.5).** Merge into the
  existing `[tool.ruff.lint]` tables (keep the template's `select`):[^ruff]
  ```toml
  [tool.ruff.lint]
  extend-select = ["C90", "PLR0911", "PLR0912", "PLR0913", "PLR0915", "PLR0917", "T20"]

  [tool.ruff.lint.mccabe]
  max-complexity = 8

  [tool.ruff.lint.pylint]
  max-args = 6
  max-positional-args = 3
  max-branches = 8
  max-returns = 4
  max-statements = 30

  [tool.ruff.lint.flake8-tidy-imports.banned-api]
  "treedata".msg = "Use the TreeData helpers in biotapy._core (contracts/tree-access)."
  "networkx".msg = "Tree operations live in biotapy/_core/_tree.py (contracts/tree-access)."
  ```
  Add to the existing per-file ignores: `"src/biotapy/_core/_tree.py" = ["TID251"]`
  and `TID251` to the `tests/**` entry.
- [x] **Step 2: Prove the ban works.** Create `src/biotapy/_core/_scratch.py`
  containing `import treedata` and `from networkx import DiGraph`; run
  `uvx prek run ruff-check --all-files` (use the hook id the generated config shows).
  Expected: two `TID251` errors. Delete the file.
- [x] **Step 3: mypy strict.** In the existing mypy configuration set
  `strict = true`, scoped to `src/biotapy` (R7.2; tests stay unannotated): if
  the hook also checks `tests/`, add an override `module = ["tests.*"]` with
  `disallow_untyped_defs = false`. Run `uvx prek run mypy --all-files`. For each third-party
  module reported as missing stubs, add it to one override:
  ```toml
  [[tool.mypy.overrides]]
  module = ["<reported modules>"]
  ignore_missing_imports = true
  ```
- [x] **Step 4: import-linter.** Add `import-linter` to the `dev` group, then:[^import-linter]
  ```toml
  [tool.importlinter]
  root_package = "biotapy"

  [[tool.importlinter.contracts]]
  name = "Layers (contracts/module-boundaries)"
  type = "layers"
  containers = ["biotapy"]
  layers = [
      "(pl) | (ml) | (da)",
      "(tl) | (fn)",
      "(pp)",
      "(datasets)",
      "(io)",
      "_core",
  ]
  ```
  and a local hook in `.pre-commit-config.yaml`:
  ```yaml
    - repo: local
      hooks:
        - id: import-linter
          name: import-linter (contracts/module-boundaries)
          entry: uv run --group dev lint-imports
          language: system
          pass_filenames: false
          types: [python]
  ```
- [x] **Step 5: Prove the layers bind.** Temporarily create
  `src/biotapy/pp/__init__.py` (empty) and add `import biotapy.pp` to
  `src/biotapy/_core/_rng.py`. Run `uv run --group dev lint-imports`.
  Expected: contract BROKEN, naming `biotapy._core._rng -> biotapy.pp`.
  Revert both changes and rerun: contract KEPT.
- [x] **Step 6: pytest configuration.** Merge into the existing
  `[tool.pytest.ini_options]`:
  ```toml
  addopts = ["--import-mode=importlib", "--doctest-modules", "--strict-markers", "-m", "not network and not r"]
  testpaths = ["tests", "src/biotapy"]
  markers = [
      "golden: compares against R or HUMAnN golden files (contracts/r-golden-parity)",
      "network: downloads data; runs only in the dedicated CI job",
      "r: needs R and rpy2 (extra `r`)",
  ]
  ```
  Run `uv run --group test pytest -q` -> the 6 `_core` tests pass and doctest
  collection of `src/biotapy` reports no errors.
- [x] **Step 7: Full gate** - `uvx prek run --all-files` -> all hooks pass.
- [x] **Step 8: Commit** - `git commit -am "build: enforce size limits, strict typing, layer contracts and test config"`

### Task 0.5: OKF conformance test and staleness script

**Files:** create `tests/test_knowledge_bundle.py`, `scripts/knowledge_stale.sh`;
modify `pyproject.toml` (`hypothesis`, `pyyaml` in the `test` group).

- [x] **Step 1: Write the test**
  ```python
  # tests/test_knowledge_bundle.py
  """OKF v0.2 conformance for .knowledge/ (rules.md R12.2)."""

  from pathlib import Path

  import pytest
  import yaml

  BUNDLE = Path(__file__).resolve().parents[1] / ".knowledge"
  RESERVED = {"index.md", "log.md"}
  CONCEPTS = sorted(p for p in BUNDLE.rglob("*.md") if p.name not in RESERVED)


  def _frontmatter(path: Path) -> dict[str, object]:
      text = path.read_text(encoding="utf-8")
      if not text.startswith("---\n"):
          return {}
      block, found, _ = text[4:].partition("\n---\n")
      return yaml.safe_load(block) or {} if found else {}


  def _rel(path: Path) -> str:
      return path.relative_to(BUNDLE).as_posix()


  def test_bundle_has_concepts():
      assert CONCEPTS


  @pytest.mark.parametrize("path", CONCEPTS, ids=_rel)
  def test_concept_has_a_type(path: Path):
      assert str(_frontmatter(path).get("type") or "").strip(), f"{_rel(path)} needs frontmatter with a non-empty `type`"


  @pytest.mark.parametrize("path", CONCEPTS, ids=_rel)
  def test_concept_is_listed_in_its_index(path: Path):
      index = path.parent / "index.md"
      assert f"({path.name})" in index.read_text(encoding="utf-8"), f"{_rel(path)} missing from {_rel(index)}"


  def test_only_the_root_index_has_frontmatter():
      nested = [p for p in BUNDLE.rglob("index.md") if p.parent != BUNDLE]
      assert not [_rel(p) for p in nested if p.read_text(encoding="utf-8").startswith("---")]
  ```
- [x] **Step 2: Run, expect pass** on the current bundle:
  `uv run --group test pytest tests/test_knowledge_bundle.py -q`.
- [x] **Step 3: Prove it fails.** Create `.knowledge/decisions/scratch.md` with
  body text only. Rerun. Expected: FAIL naming `decisions/scratch.md` in both
  parametrized tests. Delete the file.
- [x] **Step 4: Staleness script.** Copy the codebase-map skill's checker:
  ```bash
  mkdir -p scripts
  cp ~/.claude/skills/codebase-map/scripts/stale.sh scripts/knowledge_stale.sh
  chmod +x scripts/knowledge_stale.sh
  bash scripts/knowledge_stale.sh; echo "exit=$?"
  ```
  Expected: a summary line `N current, 0 stale, M uncheckable`, exit 0.
- [x] **Step 5: Commit** - `git add -A && git commit -m "test: enforce OKF conformance of the knowledge bundle"`

### Task 0.6: CI matrix, no-extras import, PR template, CONTRIBUTING

**Files:** modify `.github/workflows/test.yaml`, `docs/contributing.md`;
create `.github/pull_request_template.md`.

- [x] **Step 1: OS matrix.** Read `.github/workflows/test.yaml`. In the test
  job add `os: [ubuntu-latest, macos-latest, windows-latest]` to the existing
  matrix, set `runs-on: ${{ matrix.os }}` and `fail-fast: false`.
- [x] **Step 2: Import without extras.** Add a job to the same workflow,
  copying `checkout` and `setup-uv` steps verbatim from the existing job
  (the template's zizmor hook audits action pinning):
  ```yaml
    import-without-extras:
      runs-on: ubuntu-latest
      steps:
        # <checkout step copied from the test job>
        # <setup-uv step copied from the test job>
        - run: uv run --no-dev python -c "import biotapy, biotapy._core"
  ```
- [x] **Step 3: Knowledge freshness report** (informational, never blocks):
  ```yaml
    knowledge-touched:
      if: github.event_name == 'pull_request'
      runs-on: ubuntu-latest
      continue-on-error: true
      steps:
        # <checkout step copied from the test job, with fetch-depth: 0>
        - env:
            BASE: ${{ github.base_ref }}
          run: bash scripts/knowledge_stale.sh --touched --against "origin/$BASE"
  ```
- [x] **Step 4: PR template** `.github/pull_request_template.md`:
  ```markdown
  ## Task
  Phase / task id from `.knowledge/roadmap/`:

  ## Checklist (rules.md)
  - [ ] Task is in the active phase (R1.1); nothing outside it changed (R1.4)
  - [ ] Reused a library call, or explained why none fits (R2.1): `<call>`
  - [ ] No new parameter, abstraction or dependency without a test or approval (R2.3, R9.1)
  - [ ] Docstring (R equivalent, Guide, example) and docs page updated (R8.2)
  - [ ] Knowledge concepts updated, or confirmed unaffected (R12.1)
  - [ ] Compiled engine only: asv benchmark >= 5x attached (R10.3)

  ## Verification
  Commands run and their result (R14):
  ```
- [x] **Step 5: CONTRIBUTING.** At the top of `docs/contributing.md` add:
  "All contributions follow
  [rules.md](https://github.com/pedrocr83/biotapy/blob/master/rules.md). Design
  knowledge lives in the
  [.knowledge bundle](https://github.com/pedrocr83/biotapy/tree/master/.knowledge);
  new functions follow the
  [add-a-function playbook](https://github.com/pedrocr83/biotapy/blob/master/.knowledge/playbooks/add-a-function.md)."
- [x] **Step 6: Verify locally** - `uvx prek run --all-files` (zizmor and biome included) passes.
- [x] **Step 7: Commit** - `git add -A && git commit -m "ci: OS matrix, no-extras import check, PR template"`
- [x] **Step 8: Push and open a PR** - only after the user approves the push (R13.3).
  Expected: 12 test cells (3 pre-release cells may fail), `lint`, `import-without-extras`
  and the build check green. Merge with a merge commit, not a squash, so the
  `commit` SHAs in `.knowledge/` stay reachable from `master`.

### Task 0.7: Docs site

**Files:** modify `docs/index.md`; create `docs/design.md`.

- [x] **Step 1:** create `docs/design.md`: a short page stating that design
  decisions, contracts and the roadmap live in `.knowledge/` (OKF v0.2), with a
  link to `https://github.com/pedrocr83/biotapy/tree/master/.knowledge`.
  Add it to the `docs/index.md` toctree.
- [x] **Step 2: Build with warnings as errors**
  `uv run --group doc sphinx-build -W -b html docs docs/_build/html` -> `build succeeded`.
- [x] **Step 3: Commit** - `git add -A && git commit -m "docs: add design page pointing to the knowledge bundle"`
- [x] **Step 4: Read the Docs** - the user imports the GitHub repo on
  readthedocs.org (outward action). Expected: first build green.

### Task 0.8: Placeholder release 0.0.1

**Files:** create `.knowledge/playbooks/cut-a-release.md`; modify
`.knowledge/playbooks/index.md`, `.knowledge/log.md`.

- [x] **Step 1: Trusted publisher (user action).** On PyPI -> Publishing ->
  "Add a new pending publisher": project `biotapy`, owner `pedrocr83`, repo
  `biotapy`, workflow `release.yaml`, environment as named in the generated
  `.github/workflows/release.yaml` (read it first).
- [x] **Step 2: Write the playbook** `.knowledge/playbooks/cut-a-release.md`
  (`type: Playbook`): update the changelog, tag `vX.Y.Z`, publish a GitHub
  release, watch `release.yaml`, verify on PyPI. List it in the playbooks index; log it.
- [x] **Step 3: Release (each command needs user approval)**
  ```bash
  git switch master && git pull --ff-only
  git tag v0.0.1
  git push origin v0.0.1
  gh release create v0.0.1 --title "0.0.1" --notes "Name reservation; no functionality yet."
  ```
- [x] **Step 4: Verify**
  ```bash
  curl -s https://pypi.org/pypi/biotapy/json | python3 -c "import json,sys; print(json.load(sys.stdin)['info']['version'])"
  uv run --no-project --with biotapy==0.0.1 python -c "import biotapy; print(biotapy.__version__)"
  ```
  Expected: `0.0.1` twice.

# Deviations (rulings made during execution)
Recorded here because the execution ledger is not committed.
- Skeleton generated non-interactively (`cruft create --no-input --extra-context`), keys from the template's `cookiecutter.json`.
- The template (v0.8.0) ships no `docs/template_usage.md`, no mypy and no `autofix.yaml`; its `conftest.py` only served the deleted example test and was removed. (Phase 1 added a repo-root `conftest.py` for the figure-closing fixture shared by `tests/` and the `src/biotapy` doctests.)
- ruff pre-commit hook bumped v0.15.21 -> v0.16.9: `PLR0917` is preview-only in 0.15.
- `mypy` added as a dev dependency with a local prek hook; `[tool.mypy]` strict on `src/biotapy`, Python 3.12. User confirmed. Phase 1 widened `files` to `docs/extensions` too, and the hook now runs with the `doc` group.
- pytest options live in the template's native `[tool.pytest]` table (pytest 9, `strict = true` implies strict markers), not `[tool.pytest.ini_options]`.
- Workflows, README badge and `docs/conf.py` target `master`, not the template's `main`.
- The test job runs `bash` on every OS (template steps use POSIX syntax) and includes the OS in its name.
- CI gained a `lint` job running every prek hook, required by the `check` job (final review finding: rules were only enforced locally).
- `knowledge_stale.sh --touched` exits 2 when it cannot diff; the CI report goes to the job summary and fails only on that error.
- `uv.lock` is not committed (`/uv.lock` in `.gitignore`).
- Outward actions (push, PR, Codecov, Read the Docs, PyPI) batched into one approval request.
- Codecov needed no manual setup: the first OIDC upload registered the repo.
- First release attempt failed before upload: the pinned `pypa/gh-action-pypi-publish` v1.14.0 bundled a twine that rejects
  Metadata-Version 2.5, which hatchling 1.32 writes. Fixed by merging Dependabot PR #2 (publish action v1.14.2 with twine 7,
  plus checkout, setup-uv v10, codecov and alls-green bumps; user-approved). Tag `v0.0.1` was moved to the fixed commit
  `68d8b15` (nothing had been published).
- Results: PR #1 merged as `d361ebe`, docs live at https://biotapy.readthedocs.io, `biotapy 0.0.1` on PyPI from `68d8b15`.

# Exit gate
- [x] `uvx prek run --all-files` and `uv run --group test pytest` green locally; output read.
- [x] CI green: 12 test cells (pre-release cells may fail), `lint`, `import-without-extras`, build; docs. (The Phase 0 job set; Phase 1 added `network` and `docs` jobs.)
- [x] Docs live on Read the Docs.
- [x] `biotapy 0.0.1` on PyPI.
- [x] rules.md R14 commands match the generated template.
- [x] This concept: `phase_state: done`; Phase 1 set to `in-progress` and moved
  under "Active phase" in the roadmap index; log entry written.
- [x] **Execution method for Phase 1.** Phase 0 runs Native (user's choice,
  2026-09-26). Before Phase 1 starts, remind the user and ask again:
  subagent-driven is recommended from Phase 1 because interfaces are fixed and
  a flawed `_core` helper spreads into every later task.

[^spec]: Python Microbiome Toolkit development report, sections Roadmap and Tooling
[^cookiecutter]: cookiecutter-scverse v0.8.0
[^ruff]: ruff settings (0.16.9)
[^import-linter]: import-linter layers contract (2.15)
