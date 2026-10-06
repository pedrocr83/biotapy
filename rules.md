# biotapy development rules

These rules are binding for every change, human or agent. `CLAUDE.md` imports
this file, so it is in context in every session. Rules are numbered so reviews
can cite them ("violates R4.3"). When a rule and a request conflict, stop and
ask; do not silently pick one.

Background and reasoning live in the knowledge bundle at
[.knowledge/index.md](.knowledge/index.md) (OKF v0.2). Rules state *what*; the
linked concepts state *why*.

---

## R0. Session protocol

- **R0.1** Before any code change, read in this order: `.knowledge/index.md` ->
  `.knowledge/roadmap/index.md` -> the active phase concept -> the contract and
  module concepts the task touches. Then open the source.
- **R0.2** The code is the truth, the knowledge bundle is a map. Never state a
  fact from `.knowledge/` without opening the file it points to. If they
  disagree, the code wins and the concept is fixed in the same change.
- **R0.3** Read `tasks/lessons.md` if it exists. After any correction from the
  user, add the pattern there.

## R1. Scope and phase discipline

- **R1.1** Work only on tasks listed in the active phase concept
  (`phase_state: in-progress`). Anything else: stop and ask.
- **R1.2** Never start a later phase before the current phase's exit gate is
  ticked. Scope creep is treated as a defect, not a bonus. Exactly one phase
  is `in-progress` at a time; it is listed under "Active phase" in
  `.knowledge/roadmap/index.md`.
- **R1.2a** Before starting a task that has no step-by-step section in its
  phase concept, expand it into TDD steps with superpowers:writing-plans and
  get the user's approval.
- **R1.3** Before editing, name in your plan what you will *not* touch.
- **R1.4** One task = one atomic change. No drive-by refactors, renames,
  reformatting, import reordering or "while I'm here" fixes. Found a problem
  elsewhere? Report it; do not fix it in this change.

## R2. Reuse before writing

- **R2.1** Before writing logic, check in order: AnnData/TreeData/MuData ->
  scikit-bio -> SciPy/NumPy -> pandas -> scikit-learn. If one does the job,
  the function is a thin wrapper (3-10 lines).
- **R2.2** Never use an API, argument or config key you have not confirmed
  exists in the installed version (read the source, the docs via Context7, or
  run it). Guessing an API is a defect.
- **R2.3** No speculative parameters, options, abstractions or "future-proof"
  hooks. A parameter exists only if a test or a real use case needs it.

## R3. Public functions

Full contract: [.knowledge/contracts/function-shape.md](.knowledge/contracts/function-shape.md).

- **R3.1** One task = one public function:
  `verb(data, <required>, *, <options>) -> result`. Options are keyword-only.
- **R3.2** Annotate `data` with the widest type that works: `AnnData` unless the
  tree is needed (`TreeData`) or several modalities are (`MuData`).
  `da.consensus` takes `da` result tables and `pl.consensus` its table instead.
- **R3.3** Purity ([pure-by-default](.knowledge/decisions/pure-by-default.md)):
  `io`/`pp`/`fn`/`da` return new objects and never mutate input; `tl` returns
  its result, and writes to the documented slot only with `inplace=True`;
  `pl` returns `matplotlib.axes.Axes`.
- **R3.4** Stochastic functions take
  `seed: int | np.random.Generator | None = None` and convert it once with
  `biotapy._core.as_generator`. Never touch global RNG state.
- **R3.5** Validate inputs only at the public boundary, with `_core` validators.
  Raise `ValueError`/`KeyError`/`TypeError` whose message names the argument.
- **R3.6** No classes, except where an external protocol requires one
  (scikit-learn estimators, `torch.utils.data.Dataset`) - and then one level of
  inheritance from that library's base, no biotapy base classes.
- **R3.7** No `**kwargs` pass-through except a documented `plot_kwargs` in `pl`.

## R4. Modularization

Full contract: [.knowledge/contracts/module-boundaries.md](.knowledge/contracts/module-boundaries.md).

- **R4.1** Layout: `src/biotapy/<module>/__init__.py` holds only imports and
  `__all__`; logic lives in private topic files `src/biotapy/<module>/_<topic>.py`.
- **R4.2** Layers, higher may import lower, never the reverse:
  `pl | ml | da` > `tl | fn` > `pp` > `datasets` > `io` > `_core`.
  Modules on the same layer do not import each other. `_core` imports only
  third-party packages.
- **R4.3** A private helper stays in its topic file when one public function
  uses it; it moves to `_core` when two or more modules use it. Never
  copy-paste logic between modules.
- **R4.4** A single-use private helper is allowed only to keep a function
  within R5 limits, and only after R2.1 found no library call.
- **R4.5** Only `biotapy/_core/_tree.py` imports `treedata` or `networkx`
  ([tree-access](.knowledge/contracts/tree-access.md)); everything else uses
  the `TreeData` type and tree helpers re-exported by `biotapy._core`
  (ruff `TID251`).
- **R4.6** Optional dependencies are imported inside the function via
  `biotapy._core.import_optional(name, extra=...)`. Never at module level.
- **R4.7** No module-level side effects other than constants.
- **R4.8** Create a subpackage only in the phase that fills it; no empty
  placeholders.
- **R4.9** Tests mirror source: `src/biotapy/pp/_glom.py` ->
  `tests/pp/test_glom.py`. Tests call the public API (`bt.pp.tax_glom`), not
  private files, except `_core` unit tests.

## R5. Size and complexity limits (enforced by ruff)

| Limit | Value | Enforced by |
|---|---|---|
| Statements per function | <= 30 | ruff `PLR0915` |
| Cyclomatic complexity | <= 8 | ruff `C901` |
| Branches per function | <= 8 | ruff `PLR0912` |
| Arguments per function | <= 6 | ruff `PLR0913` |
| Return statements | <= 4 | ruff `PLR0911` |
| Lines per topic file | <= 300 | review |
| Public functions per topic file | <= 6 | review |

- **R5.1** Exceeding a limit means split the function or find the library call
  you missed. `# noqa` for these rules requires a one-line reason and user approval.

## R6. Data rules

Full contract: [.knowledge/contracts/data-model-slots.md](.knowledge/contracts/data-model-slots.md).

- **R6.1** Samples are rows, always. Only `io` readers transpose, exactly once.
- **R6.2** `X` stays sparse CSR. Never call `.toarray()`/`.todense()` on a full
  matrix, except inside a wrapper whose delegated library requires dense input
  (scikit-bio does: its AnnData dispatch fails on sparse `X`), or inside a
  native method whose algorithm needs the full table (LinDA's log-ratios have
  no zeros). Then densify once, in that function, with a comment saying why,
  and state the memory cost in the docstring `Notes`.
- **R6.3** Results go only to the slot and key the data-model contract names
  (`layers`, `obsm`, `obsp`, `uns["biotapy"]`). No ad-hoc keys.
- **R6.4** Functions that change features (filter, glom, rarefy) follow the
  contract's propagation rules: drop derived slots they would invalidate.
- **R6.5** Never return an AnnData view to the user; return a real object.
- **R6.6** Never commit data files larger than 1 MB. Datasets are fetched and
  cached with pooch; tests use `bt.datasets.toy()` or tiny fixtures under
  `tests/data/`.

## R7. Types, errors, logging

- **R7.1** Full type hints on every function, public and private. No `Any` in
  public signatures; `Literal[...]` for string options.
- **R7.2** `mypy --strict` passes on `src/biotapy`.
- **R7.3** No `print` (ruff `T201`). Use `logging.getLogger(__name__)` for
  diagnostics and `warnings.warn(..., stacklevel=2)` for user-facing warnings.
- **R7.4** No bare `except`, no `except Exception: pass`, no silent fallbacks
  (including silent engine fallbacks).

## R8. Comments, docstrings, documentation

- **R8.1** Inline comments only for the non-obvious *why* (numerical trick, R
  compatibility quirk). Never for *what*.
- **R8.2** Every public function has a NumPy docstring with: summary,
  Parameters, Returns, a `Notes` section containing exactly one
  `R equivalent: ``pkg::fn``` line (or `R equivalent: none`) and a `Guide:`
  link, a runnable `Examples` section using `bt.datasets.toy()` (or
  `bt.datasets.toy_humann()` for function tables; a function that reads a file
  may write a small temp file in its example), and References where a method
  is cited.
- **R8.3** Where documentation goes
  ([docs-okf-and-sphinx](.knowledge/decisions/docs-okf-and-sphinx.md)):

  | Content | Location |
  |---|---|
  | Why a decision was made, contracts, playbooks, roadmap | `.knowledge/` (OKF v0.2) |
  | How to use a function, method math, tutorials, API | `docs/` (Sphinx) |
  | What one function does | its docstring |

- **R8.4** The "Coming from R" table is generated from docstrings; never edit
  it by hand.

## R9. Dependencies

- **R9.1** Never add, remove or bump a dependency (core, extra or dev) without
  asking the user first.
- **R9.2** Each approved dependency gets a written reason in the PR and, for
  runtime dependencies, a line in
  [optional-heavy-dependencies](.knowledge/decisions/optional-heavy-dependencies.md).
- **R9.3** Heavy dependencies (torch, rpy2, plotnine, numba, unifrac) are
  extras only.

## R10. Performance

- **R10.1** Profile (py-spy or scalene) before optimizing. No optimization
  without a measurement.
- **R10.2** Escalation ladder only: vectorize -> delegate to a library -> Numba
  -> Rust. Never new C/C++
  ([python-first-compiled-last](.knowledge/decisions/python-first-compiled-last.md)).
- **R10.3** A compiled engine merges only with an asv benchmark showing >= 5x
  and a parity test against the Python engine
  ([engine-parity](.knowledge/contracts/engine-parity.md)).

## R11. Testing

- **R11.1** TDD: write the failing test, run it, see it fail for the right
  reason, then implement.
- **R11.2** Every public function has: a happy path; edge cases (all-zero
  sample, all-zero feature, NaN/missing taxonomy rank, single sample); a purity
  assertion; a Hypothesis property test when an invariant exists; a golden
  test when an R equivalent exists, except in `pl`, whose plots draw `tl`
  results that have one
  ([r-golden-parity](.knowledge/contracts/r-golden-parity.md)).
- **R11.3** Tolerances are explicit in the test. Any tolerance looser than the
  contract default carries a one-line reason.
- **R11.4** Tests are deterministic: fixed seeds, no network (datasets tests
  use a local pooch cache fixture), no sleeps.
- **R11.5** Never weaken, skip or delete a test to make a change pass.
- **R11.6** Coverage is measured in CI; target >= 90% on public functions.
  Never write a test only to move the number.

## R12. Knowledge bundle (OKF)

Playbook: [maintain-knowledge](.knowledge/playbooks/maintain-knowledge.md).

- **R12.1** Update a concept in the same commit when the change invalidates
  something it states. Most changes do not; do not churn concepts.
- **R12.2** Every non-reserved `.md` in `.knowledge/` has YAML frontmatter
  with a non-empty `type`. `index.md` and `log.md` are reserved names.
- **R12.3** Agents write `generated: { by: claude-code/<model>, at: <UTC> }`.
  Only a human adds `verified`.
- **R12.4** Tick task checkboxes in the phase concept as tasks complete, and
  add a dated line to `.knowledge/log.md` for every concept created or changed.

## R13. Git

- **R13.1** Conventional commits (`feat:`, `fix:`, `refactor:`, `docs:`,
  `test:`, `build:`, `ci:`, `chore:`), one logical change per commit.
- **R13.2** Never commit secrets, `.env`, credentials, large data, or build
  outputs.
- **R13.3** Never push, tag, publish to PyPI, or change GitHub settings without
  explicit user approval for that specific action.

## R14. Definition of done

A task is done only when all of these ran and their output was read
(no `docker-compose.yml` here, so they run on the host through uv):

```bash
uvx prek run --all-files            # ruff lint + format, mypy --strict, import-linter, pyproject-fmt
uv run --group test pytest          # unit, property, doctest, golden and tests/test_knowledge_bundle.py
```

Plus, when docs changed: `uv run --group doc sphinx-build -W -b html docs docs/_build/html`.

These commands are fixed in Phase 0 (tasks 0.2-0.3). If the generated
template names things differently, update this section in the same commit.

- **R14.1** Report results faithfully: failures are quoted, skipped steps are named.
- **R14.2** Stop condition: if the same check fails twice in a row, stop
  changing code and report what was tried and what failed.
