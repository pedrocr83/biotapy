---
type: Decision
title: Read phyloseq via rdata, refseq deferred
description: A ~35-line rdata constructor_dict reads GlobalPatterns/enterotype/esophagus with zero shape mismatches and zero residual warnings; native rdata (+xarray) is the route for Tasks 1.9/1.10, an R export script is rejected as the default, and a populated refseq is warned-and-skipped rather than guessed at (no test fixture has one).
tags: [io, dependencies, phyloseq, spike]
status: draft
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T20:17:19Z }
commit: 9872d9a
sources:
  - id: phyloseq-classes
    resource: https://github.com/joey711/phyloseq/blob/master/R/allClasses.R
    title: phyloseq S4 class definitions (phyloseq, otu_table, taxonomyTable, sample_data)
  - id: rdata
    resource: https://github.com/vnmabus/rdata
    title: rdata 1.1.0 - Python reader for R .RData/.rds files
  - id: rdata-pypi
    resource: https://pypi.org/pypi/rdata/json
    title: rdata 1.1.0 PyPI metadata (runtime dependencies)
---

# Context
Task 1.6 spikes whether `rdata.read_rda` with a custom `constructor_dict` can
turn phyloseq's example `.RData` objects (S4 `phyloseq`, slots `otu_table`,
`tax_table`, `sam_data`, `phy_tree`, `refseq`) into counts, taxonomy, sample
data and a tree, well enough to base a native reader (Tasks 1.9 `.rds`, 1.10
`read_phyloseq`) on it - or whether biotapy should instead ship a short R
script (docs-only) that exports BIOM + Newick + TSV for users to read with
biotapy's existing importers.[^phyloseq-classes]

Spike code and raw output live in `.superpowers/sdd/phase-1-core/spike/`
(git-ignored, kept until the user approves this decision, per Step 6).

# Options

**A. Native `rdata` route.** `bt.io.read_phyloseq`/`.rds` support reads the
file directly in Python via `rdata.read_rda`/`read_rds` plus a small
`constructor_dict` (one function per phyloseq S4 class). New runtime deps:
`rdata` and `xarray`.[^rdata][^rdata-pypi] Fragility: depends on phyloseq's S4
layout staying stable (a `setClass` change would surface as a loud "Missing
constructor" `UserWarning`, not a silent wrong answer, per rules.md R7.4) and
on `rdata` continuing to expose `constructor_dict` the same way.

**B. R export script (docs only).** Ship a short R script under `docs/` that
a phyloseq/DADA2 user runs once, writing BIOM (otu_table + tax_table),
Newick (phy_tree) and a sample-data TSV; biotapy reads only formats it needs
for BIOM/QIIME2 support anyway. Zero new Python deps, but requires the user
to have R, phyloseq and biomformat installed just to get their data into
biotapy - a heavy ask for a Python-only toolkit's target user, and it forces
`.rds`/`.RData` input (Task 1.9) onto the same script.

**C. Hybrid - native for 4 slots, refseq warned-and-skipped (recommended).**
Option A for `otu_table`, `tax_table`, `sam_data` and `phy_tree`, all
verified below. A populated `refseq` is **not** built into a Biostrings
`XStringSet` (no phyloseq example object has one to test against - rule
R2.2: never implement against an unverified API), but it is also not a hard
failure: `read_phyloseq` reads the four other slots and emits one
`UserWarning` naming the `refseq` slot and the one-line R export
(`Biostrings::writeXStringSet`) a user can run to get the sequences
separately. This is a warned skip, not a silent fallback (rules.md R7.4 bars
silent fallbacks, not warned ones) - the four readable slots parse
correctly regardless of whether `refseq` happens to be populated.

# Rejected
- **Hard error (`NotImplementedError`) on a populated `refseq`.** Considered
  as part of Option C, rejected: the standard DADA2-to-phyloseq handoff in
  the DADA2 tutorial, `merge_phyloseq(ps, DNAStringSet(taxa_names(ps)))`,
  populates `refseq` on a large share of real-world objects. Hard-failing
  the whole read would block every user whose other four slots parse fine
  just because `refseq` happens to be set - the spec's and a reasonable
  user's expectation is that `read_phyloseq` succeeds on the majority of
  real DADA2/phyloseq objects, not only on the three slot-sparse example
  files this spike used.

# Evidence

Run: `uv run --no-project --with rdata python spike.py`, from
`.superpowers/sdd/phase-1-core/spike/`. Files downloaded with `curl` from
`raw.githubusercontent.com/joey711/phyloseq/master/data/`; sizes matched the
brief exactly (`GlobalPatterns.RData` 435,652 B, `enterotype.RData`
195,260 B, `esophagus.RData` 1,840 B), confirming reproducibility. No
install block: `uv` resolved and cached `rdata` 1.1.0, `numpy` 2.5.3,
`xarray` 2026.7.0, `pandas` 3.0.6, `typing_extensions` 4.16.0 in ~4s.

| | GlobalPatterns | enterotype | esophagus |
|---|---|---|---|
| Raw read (`DEFAULT_CLASS_MAP`) warnings | 5 (phylo, sample_data, taxonomyTable, otu_table, phyloseq) | 4 (no `phylo` - slot is NULL, untyped) | 3 (no taxonomyTable/sample_data - both NULL) |
| `otu_table` after constructor | `DataArray` (19216, 26), `taxa_are_rows=True` | `DataArray` (553, 280), `taxa_are_rows=True` | `DataArray` (58, 3), `taxa_are_rows=True` |
| `tax_table` after constructor | `DataFrame` (19216, 7), cols `Kingdom..Species` | `DataFrame` (553, 1), col `Genus` | `None` (NULL slot) |
| `sam_data` after constructor | `DataFrame` (26, 7) | `DataFrame` (280, 9) | `None` (NULL slot) |
| `phy_tree` after constructor | 19216 tips, 38430 edges | `None` (no tree) | 58 tips, 114 edges |
| `refseq` | `None` (NULL slot) | `None` (NULL slot) | `None` (NULL slot) |
| Constructed-read time | 2.306 s | 0.060 s | 0.003 s |
| Residual warnings | 0 | 0 | 0 |
| Known-shape mismatches (Step 4) | none | none | none |

Checks against the brief's known facts: GlobalPatterns is 26 samples x
19,216 taxa with 7 rank columns and a 19,216-tip tree - all confirmed exactly.
enterotype's `otu_table` is relative abundance (column sums 0.99986-1.00000,
values in [0, 0.897], checked directly, not just "looks fractional").
esophagus is 3 samples with a 58-tip tree - confirmed; it genuinely has no
`tax_table`/`sam_data` in phyloseq itself (both NULL), not a spike bug.

Both trees were independently checked with `networkx` in the repo's own
`.venv` (never importing `rdata` there, per the brief), reading a JSON dump
of `(parent, child, length)` triples: `nx.is_tree` is `True` for both, the
leaf set (`out_degree == 0`) equals the tip-label set exactly, each has one
root (`in_degree == 0`: `n19217` and `n59`), and edge counts match a fully
resolved binary tree (`2 * n_tip - 2`).

**Constructor code.** ~35 non-blank, non-comment lines total for all 5
classes: `otu_table` (3 lines: read `taxa_are_rows`, wrap), `taxonomyTable`
(9 lines: reshape the flat character array from `attrs["dim"]` with
`order="F"`, apply `attrs["dimnames"]`), `phylo` (13 lines: build
`(parent, child, length)` triples from `edge`/`edge.length`/`tip.label`,
naming internal nodes `n<i>`), `phyloseq` (1 line, pure passthrough), plus a
6-line `constructor_dict` registration. `sample_data` needed **zero** new
code - `rdata.conversion.dataframe_constructor` (already in
`DEFAULT_CLASS_MAP` for `"data.frame"`) was registered for `"sample_data"`
unchanged and worked immediately, because phyloseq's `sample_data` "contains"
`data.frame` at the R level and serializes with the same shape.

**The NULL-slot sentinel (the one real gotcha).** An absent `...OrNULL` slot
(`refseq` in all three files, plus `phy_tree`/`tax_table`/`sam_data` where
absent) does not deserialize to Python `None` or trigger any class-based
dispatch. It comes back as the literal 6-character string `"\x01NULL\x01"`,
carrying **no R `class` attribute at all** - confirmed empirically
(`spike_diag_null.py`) and confirmed absent from `rdata`'s own source and
test suite (`grep`, `_conversion.py`, `tests/test_read.py`): this is not an
`rdata` feature, it is R's own internal representation of an unset S4
`...OrNULL` slot, produced before serialization. Consequently: (a) a custom
`constructor_dict` entry is **never called** for a NULL slot, so the NULL
check has to live in the caller, not inside each constructor (this cost one
debugging pass in the spike - `phylo_constructor` first crashed on
enterotype's NULL `phy_tree` before this was understood); (b) because it is
R's own representation rather than an `rdata` quirk, no alternative
Python-side R deserializer would avoid it either.

**`refseq`/`XStringSet`.** `refseq` is `NULL` in all three example files, so
a Biostrings `XStringSet` constructor is completely unverified - Biostrings
stores sequences in a shared/compressed pool (`XVectorList`) indexed by
ranges, materially more complex than the other four slots, and phyloseq
ships no bundled example with a populated `refseq` to test a constructor
against.

# Recommendation
**Option C.** The native route reads 4 of 5 slots correctly with modest,
mostly-reusable code (~35 lines, one slot free), fast (worst case 2.3 s on
the largest file), no shape mismatches anywhere, and zero residual warnings
once all classes are registered. It adds two runtime dependencies
(`rdata`, `xarray`) on top of `numpy`/`pandas`, which biotapy already
depends on - `xarray` is not on the R9.3 heavy-dependency list, so nothing
here forces it into an extra. Option B trades that for requiring a working R
+ phyloseq + biomformat install from every user, which is a heavier and
less biotapy-native ask, and still would not avoid Task 1.9's `.rds`
requirement. `refseq` is left unimplemented rather than guessed at (R2.2):
a populated `refseq` is warned-and-skipped, not hard-failed, because it is
common enough in real DADA2/phyloseq objects (see Rejected) that failing
the whole read on it would be worse than skipping it with a named warning
and a documented R fallback.

New dependencies `rdata` and `xarray` still need the user's explicit
approval (rules.md R9.1) before Task 1.9/1.10 add them to `pyproject.toml`.

# Consequences
- **Task 1.9 (`.rds` input).** `rdata.read_rds` takes the same
  `constructor_dict`; reuse the 5 constructors from this spike unchanged
  (`read_rda` and `read_rds` share the same conversion machinery per the
  Task 1.6 research notes).
- **Task 1.10 (`read_phyloseq`).** Build on the 4 verified constructors.
  `refseq` populated and non-NULL emits one named `UserWarning` (naming the
  slot and the `Biostrings::writeXStringSet` export) and is skipped, rather
  than either guessing a constructor (R2.2) or hard-failing the whole read
  on a slot that a large share of real DADA2/phyloseq objects populate (see
  Rejected). A warned skip is not a silent fallback under R7.4.
- **Task 1.11 (datasets).** pooch can point directly at phyloseq's own
  `data/*.RData` URLs (sizes reproduced exactly in this spike) - no need to
  host pre-converted files on a biotapy release. All three files are
  fetched with pooch rather than committed, for one consistent code path
  and so every file is always sourced from upstream phyloseq with a
  pooch SHA-256 pin, not for size: at 435,652 B and 195,260 B,
  `GlobalPatterns.RData` and `enterotype.RData` are both under the R6.6
  1 MB commit limit, same as `esophagus.RData` (1,840 B).
- **Task 1.12 (golden files).** The shapes and values confirmed here
  (26x19,216 with 7 ranks and a 19,216-tip tree; enterotype's
  column-sum-to-1 relative abundances; esophagus's 58x3 with a 58-tip tree)
  become the golden assertions for `read_phyloseq` once it exists.

[^phyloseq-classes]: phyloseq S4 class definitions (phyloseq, otu_table, taxonomyTable, sample_data)
[^rdata]: rdata 1.1.0 - Python reader for R .RData/.rds files
[^rdata-pypi]: rdata 1.1.0 PyPI metadata (runtime dependencies)
