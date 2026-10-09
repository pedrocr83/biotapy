---
type: Decision
title: Embedding models are plugins found through entry points
description: ml.embed(adata, model) loads the callable a package registers under model in the entry-point group biotapy.embeddings, checks that it returned a finite 2-D float array with one row per sample, and returns it or writes obsm["X_<model>"]; biotapy's own MGM is registered the same way, its weights downloaded, never bundled.
tags: [ml, plugins, api]
status: draft
paths: ["src/biotapy/ml/_embed.py", "pyproject.toml"]
generated: { by: claude-code/claude-sonnet-5-5, at: 2026-10-09T16:38:49Z }
commit: 5b73a1b
sources:
  - id: spec
    resource: ../../plan.md
    title: Python Microbiome Toolkit development report
    author: human:pedrocr83
  - id: entry-points
    resource: https://packaging.python.org/en/latest/specifications/entry-points/
    title: Entry points specification (PyPA)
---

# Context
The spec plans foundation-model embeddings in `obsm` through plugins.[^spec]
Microbiome foundation models come and go (MGM, MGM2, BiomeGPT), each with its
own framework, weights and licence; biotapy cannot depend on all of them, and
a model's authors may want to ship it themselves. Python's packaging already
has a registry for this: entry points, which a package declares in its
`pyproject.toml` and `importlib.metadata` lists at run time.[^entry-points]

# Decision
- **Group** `biotapy.embeddings`. An entry point's name is the model's name;
  its value loads a callable `embed(adata) -> numpy.ndarray` that returns
  samples x dimensions and leaves `adata` unchanged.
- **Discovery at call time**: `ml.embed(adata, model)` reads
  `importlib.metadata.entry_points(group="biotapy.embeddings")` on every call,
  so installing a package is enough; nothing is imported until a model is
  asked for. An unknown name raises `KeyError` listing the installed names;
  a name two packages register raises `ValueError` naming both, never picks
  one (`ml/_embed.py:embed`).
- **biotapy checks the result**, not the plugin: a plain `numpy.ndarray`, not
  a masked array or a matrix (else `TypeError`), 2-D with one row per sample and at least one column, float,
  finite (else `ValueError`), each message naming the plugin
  (`ml/_embed.py:_checked`). A refused result is never stored, and one that
  shares memory with `X`, a layer, `obsm`, `varm`, `obsp`, `varp` or a
  top-level `uns` array (nothing nested deeper) is copied. A plugin's own exception
  propagates with its type and a note naming the plugin. The model name is
  letters, digits, `_`, `-` or `.`, so `obsm["X_<model>"]` is a plain key.
- **The plugin gets the caller's AnnData, not a copy**: a copy would double
  `X` in memory for every call. The protocol therefore forbids a plugin to
  modify it; biotapy does not enforce that (MGM's own test does:
  `tests/ml/test_mgm.py:test_keeps_the_input`).
- **Where it goes**: returned, or with `inplace=True` written to
  `obsm[f"X_{model}"]` with `None` returned, the `tl` convention
  ([pure-by-default](/decisions/pure-by-default.md),
  [data-model-slots](/contracts/data-model-slots.md)).
- **No options pass through**: the callable takes the AnnData alone (rules.md
  R3.7). The reference model, MGM, ran fastest one sample at a time on the CPU
  (phase-4 plan, task 4.4b), so no `batch_size` exists until a model needs it.
- **Weights are never bundled** (rules.md R6.6): a plugin downloads them at run
  time; biotapy's own go through `_core.make_pooch`, pinned by SHA-256.
- **biotapy's reference model registers itself the same way**:
  `[project.entry-points."biotapy.embeddings"] mgm = "biotapy.ml._mgm:embed"`
  in biotapy's `pyproject.toml`, behind the extra `mgm`
  ([optional-heavy-dependencies](/decisions/optional-heavy-dependencies.md)).

# Rejected
- **A registry function** (`bt.ml.register(name, fn)`): needs an import with
  a side effect before every use, and two packages would race for a name.
- **A base class plugins inherit from**: rules.md R3.6 allows classes only
  where an external protocol requires one; a callable is enough.
- **Passing `**kwargs` to the plugin**: rules.md R3.7.
- **Trusting the plugin's output**: a wrong row count would silently pair
  embeddings with the wrong samples.
- **A separate `biotapy-mgm` package**: no release cadence of its own yet
  (Phase 4 decision 9); the entry point keeps that move cheap.

# Consequences
- Third-party models need no change to biotapy; the guide shows the
  `pyproject.toml` lines (`docs/guide/machine_learning.md`, Embeddings).
- Tests install fake plugins as a distribution under `tmp_path` on
  `sys.path`, so they go through the same discovery
  (`tests/ml/test_embed.py:_installer`).
- A feature-changing step drops `obsm`, and with it `X_<model>`
  ([data-model-slots](/contracts/data-model-slots.md), Propagation).

[^spec]: Python Microbiome Toolkit development report, section Roadmap
[^entry-points]: Entry points specification (PyPA)
