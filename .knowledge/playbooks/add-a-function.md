---
type: Playbook
title: Add a public function
description: The only sanctioned path from "we need X" to a merged public function - reuse check, failing test, minimal implementation, docstring, docs page, knowledge update, verification.
tags: [workflow, api, testing]
status: stable
paths: ["src/biotapy/**", "tests/**", "docs/**"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-05T09:34:05Z }
commit: f5236e8
---

# When
Any new entry in a subpackage `__all__`.

# Steps
1. **Phase check.** The function must be listed in the active phase under
   [roadmap](/roadmap/index.md). If it is not, stop and ask.
2. **Reuse check, in order**: AnnData/TreeData -> scikit-bio -> SciPy/NumPy ->
   pandas -> scikit-learn. If one does the job, the function is a 3-10 line
   wrapper. Write the library call you found in the PR description.
3. **Place it.** Pick the topic file per [module-boundaries](/contracts/module-boundaries.md)
   (`src/biotapy/<module>/_<topic>.py`). New topic file only if no existing
   topic fits.
4. **Failing tests first** in `tests/<module>/test_<topic>.py`:
   - happy path on `bt.datasets.toy()`;
   - edge cases: empty sample (all-zero row), all-zero feature, missing
     taxonomy rank / NaN rank value, single sample;
   - purity: input unchanged (`pp`, `tl` with `inplace=False`, `pl`);
   - property test with Hypothesis when an invariant exists (sums preserved,
     rows sum to 1, symmetric distances);
   - golden test when an R equivalent exists, per [r-golden-parity](/contracts/r-golden-parity.md);
     not in `pl`, whose plots draw results of `tl`, `pp` or `fn` that have one
     (or none, as `pl.contributions`).
   Run and see them fail for the right reason.
5. **Minimal implementation** matching [function-shape](/contracts/function-shape.md).
   No option without a test that uses it.
6. **Export**: add to the subpackage `__init__.py` import list and `__all__`.
7. **Docstring** with `R equivalent:` line, `Guide:` link, runnable example. The `R equivalent:` line
   feeds the generated Coming-from-R page; `tests/test_docstrings.py` fails
   if it cannot be parsed.
8. **Docs page**: add or extend the `docs/guide/<concept>.md` page the
   docstring links to; add the function to `docs/api.md`.
9. **Knowledge**: update `.knowledge/` only if a contract, decision or module
   concept is now wrong or incomplete - see [maintain-knowledge](/playbooks/maintain-knowledge.md).
10. **Tick** the task checkbox in the active phase concept.

# Verification
```bash
uvx prek run --all-files
uv run --group test pytest tests/<module> src/biotapy/<module>
```
All green, output read, before claiming done.

# Common mistakes
- Densifying `X` (`.toarray()`) on the full matrix. Operate on CSR.
- Returning a view (`adata[:, mask]`) instead of a real object; call `.copy()`
  or build a new TreeData.
- Forgetting the function in `__all__`, so the docs and the R table miss it.
