# Knowledge bundle log

## 2026-09-26
* **Update**: [phase-1-core](roadmap/phase-1-core.md): ticked Task 1.5 (`pp.tax_glom`) step checkboxes.
* **Update**: [phase-1-core](roadmap/phase-1-core.md): ticked Task 1.4 (`pp.relative`) step checkboxes.
* **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) records treedata, networkx and types-networkx (Phase 1, task 1.3).
* **Update**: [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) records scipy, pandas and their type stubs (Phase 1, task 1.1).
* **Update**: [phase-1-core](roadmap/phase-1-core.md): Phase 1 runs subagent-driven (user's choice). Approved scipy and pandas (runtime) and pandas-stubs, scipy-stubs, types-networkx (dev) for `mypy --strict`. Added Task 1.7b `io.write_biom` at the user's request.
* **Update**: Phase 0 exit gate passed; [phase-0-foundation](roadmap/phase-0-foundation.md) is `done` and [phase-1-core](roadmap/phase-1-core.md) is `in-progress`. `biotapy 0.0.1` released; the release failure and fix are recorded in the Phase 0 Deviations and in [cut-a-release](playbooks/cut-a-release.md).
* **Update**: Phase 0 final-review follow-ups. Added a Deviations section to [phase-0-foundation](roadmap/phase-0-foundation.md); [cut-a-release](playbooks/cut-a-release.md) now pushes only the tag; [optional-heavy-dependencies](decisions/optional-heavy-dependencies.md) records numpy and session-info2; [maintain-knowledge](playbooks/maintain-knowledge.md) narrows its paths and warns against squash merges; [module-boundaries](contracts/module-boundaries.md) names the CI lint job. Refreshed `commit` to b77a226 on every concept with `paths` after re-checking it against the branch.
* **Creation**: Added the [cut-a-release](playbooks/cut-a-release.md) playbook (Phase 0, task 0.8).
* **Verification**: `human:pedrocr83` confirmed [package-name-biotapy](decisions/package-name-biotapy.md) (with Python >= 3.12), [pure-by-default](decisions/pure-by-default.md), [docs-okf-and-sphinx](decisions/docs-okf-and-sphinx.md) and [r-bridge-before-ports](decisions/r-bridge-before-ports.md); all now `stable`. Phase 0 dependencies and the cookiecutter-scverse v0.8.0 template approved. Phase 0 runs inline (Native); execution method to be re-asked at Phase 1.
* **Initialization**: Bundle created from the spec (`plan.md`) by `claude-code/claude-opus-5-5`: roadmap phases 0-5, six contracts, nine decisions, two playbooks. All concepts unverified; decisions marked `draft` await user confirmation.
