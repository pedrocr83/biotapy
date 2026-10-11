---
type: Playbook
title: Cut a release
description: Bump the version, move the changelog entry, tag, publish a GitHub release, and let release.yaml upload to PyPI through trusted publishing.
tags: [release, workflow]
status: stable
paths: ["pyproject.toml", "CHANGELOG.md", ".github/workflows/release.yaml"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-10-11T00:04:27Z }
commit: c757fa9
sources:
  - id: trusted-publishing
    resource: https://docs.pypi.org/trusted-publishers/
    title: PyPI trusted publishers
---

# When
A phase exit gate asks for a release (0.0.1 in Phase 0, 0.1 in Phase 1, ...).

# Preconditions
- On PyPI, project `biotapy` has a trusted publisher: owner `pedrocr83`, repo
  `biotapy`, workflow `release.yaml`, environment `pypi`.[^trusted-publishing]
  For the first release this is a *pending* publisher.
- `master` is green in CI.

# Steps
1. Set `version = "X.Y.Z"` in `pyproject.toml` (static; no VCS versioning).
2. In `CHANGELOG.md`, rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD`
   and open a new empty `## [Unreleased]` above it (Keep a Changelog). If the
   phase's pull requests left `## [Unreleased]` empty (0.2.0's did), write its
   entries first from `git log vPREVIOUS..master --no-merges`: `feat:` commits
   under Added, `fix:` commits that change a released function under Changed.
2b. Update `README.md` wherever it describes the previous release: it is PyPI's
   project page, and 0.1.0 replaced 0.0.1's placeholder text.
2c. Move the "not in X.Y" labels to the new version: `docs/_data/r_idioms.toml`,
   `docs/coming_from_r.md`, `tests/test_coming_from_r.py` and the phyloseq
   vignette (`docs/tutorials/phyloseq_analysis.md`), and the docstrings under `src/` (`grep -rn "not in 0\." src docs tests`;
   `pl.ordination`'s still said "not in 0.1" at 0.3.0). 0.2.0 did this in Task 2.13.
3. Commit `chore: release X.Y.Z` and merge it to `master` through a PR
   (merge commit, not squash). Before the PR, build and run `pytest` from the
   extracted sdist (`uv build --sdist`, `tar xzf`, then `pytest` inside it); it
   must pass with no errors. Check the wheel's `METADATA` too: `Version`, every
   runtime `Requires-Dist`, and a `Provides-Extra` line for each extra (0.3.0
   added `r`, slice 4B `torch`, slice 4C `mgm`), and the wheel's `entry_points.txt`
   (0.4.0 added `[biotapy.embeddings] mgm = biotapy.ml._mgm:embed`). The sdist must ship the root `conftest.py`
   (`pyproject.toml` `build.targets.sdist.include`): its tests need the marker hook. Run the sdist's tests with
   absolute cache paths (`BIOTAPY_DATA_DIR`, `XDG_CACHE_HOME`): `tests/core/test_download.py` compares
   `pooch.os_cache` paths, which a relative `XDG_CACHE_HOME` leaves unnormalised.
4. With explicit user approval for each (rules.md R13.3), tag the merged commit
   and push only the tag:
   ```bash
   git switch master && git pull --ff-only
   git tag vX.Y.Z
   git push origin vX.Y.Z
   gh release create vX.Y.Z --title "X.Y.Z" --notes-file <changelog excerpt>
   ```
   Publishing the GitHub release triggers `.github/workflows/release.yaml`
   (`uv build`, then `pypa/gh-action-pypi-publish` in environment `pypi`).

# Verification
```bash
curl -s https://pypi.org/pypi/biotapy/json | python3 -c "import json,sys; print(json.load(sys.stdin)['info']['version'])"
uv run --no-project --with biotapy==X.Y.Z python -c "import biotapy; print(biotapy.__version__)"
```
Both print `X.Y.Z`.

# Common mistakes
- Tagging before the version bump is merged: the wheel carries the old version.
- A PyPI version can never be re-uploaded; a broken release needs a new patch version.
- Pushing `master` directly: releases go through a merged PR; push only the tag.
- A stale publish action: `pypa/gh-action-pypi-publish` bundles its own twine, which can lag the metadata version
  hatchling writes (0.0.1 failed on Metadata-Version 2.5). Keep Dependabot's action updates merged. If a release
  fails before upload, fix `master`, then delete the release and tag (`gh release delete vX.Y.Z --cleanup-tag`) and
  re-tag: release workflows run the workflow file from the tagged commit.

[^trusted-publishing]: PyPI trusted publishers
