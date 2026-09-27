---
type: Playbook
title: Regenerate the R golden files
description: Rebuild the pinned R image and rerun it to regenerate golden CSVs and R-only test fixtures, with a bit-identical check before committing.
tags: [testing, r, validation]
status: stable
paths: ["tests/r/**", "tests/golden/**", "tests/data/phyloseq/**", "tests/data/dada2/**"]
generated: { by: claude-code/claude-sonnet-5, at: 2026-09-27T14:00:00Z }
commit: 6fd5344
sources:
  - id: r-golden-parity
    resource: ../contracts/r-golden-parity.md
    title: R golden parity
---

# When
- The pin in `tests/r/Dockerfile` (the R base image or the Bioconductor
  release) is bumped. This also changes which phyloseq version installs, but
  the Dockerfile does not pin phyloseq itself: `BiocManager::install("phyloseq",
  version = "3.22", ...)`'s `version=` is the Bioconductor release, not a
  phyloseq version. `tests/golden/VERSIONS.txt` records the phyloseq version
  that installed (currently 1.54.2). The bump is its own commit; regenerating
  every file is a separate, following commit (rules.md R13.1).
- A new golden function is added and needs its `.csv.gz` written.
- Never for any other reason: golden files and the R-written fixtures under
  `tests/data/phyloseq/` and `tests/data/dada2/` are never hand-edited.

# Steps
From the repo root:
```bash
docker build -t biotapy-golden tests/r
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
sha256sum tests/golden/global_patterns/*.csv.gz tests/data/phyloseq/* tests/data/dada2/* > /tmp/golden-run1.sha
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD":/work biotapy-golden
sha256sum -c /tmp/golden-run1.sha
```
- `-e HOME=/tmp` is required: R needs a writable `HOME` when the container
  runs as the host's non-root uid/gid (`--user "$(id -u):$(id -g)"`), which
  has no passwd entry inside the image.
- The `sha256sum -c` rerun must report every file `OK`. `tests/r/export_golden.R`
  writes gzip with `gzfile()`, and zlib's gzip header carries no timestamp, so
  two runs of the same image against the same inputs are bit-identical. A
  mismatch means something in the script or image is non-deterministic; do
  not commit until it is bit-identical.
- Every file under `tests/data/` and `tests/golden/` must stay under 1 MB
  (rules.md R6.6). `export_golden.R`'s `write_golden()` enforces this for the
  files it writes with a `stop()`; `tests/test_data_files.py` enforces it for
  every file at test time, so growth is caught even for files the script does
  not gate. If a file hits either limit, stop and report the size — do not
  change the golden layout without the controller's ruling.

# Verification
```bash
uv run --group test pytest tests/test_data_files.py -q
```

# Common mistakes
- Editing a `.csv.gz` or `.rds`/`.RData` fixture by hand: it will not match a
  rerun of the script and defeats the point of a pinned, reproducible source.
- Forgetting `-e HOME=/tmp`: R's startup fails trying to write to an
  unwritable home directory when the container runs as a non-root uid.
- Committing after only one run: always do the two-run, bit-identical check
  before staging generated files.
- Adding a package without its system headers: `BiocManager::install()`
  exits 0 even when a package fails to compile (e.g. a missing `zlib.h`
  broke `XVector` -> `Biostrings` -> `phyloseq` the first time this image was
  built), so `docker build` can report success on a broken image.
  `tests/r/Dockerfile`'s final `RUN Rscript -e 'stopifnot(requireNamespace(...))'`
  layer catches this at build time; keep it, and add any new package's
  `requireNamespace()` call to that check when a golden function needs one.
