---
type: Playbook
title: Maintain the knowledge bundle
description: When and how to write or update an OKF concept in .knowledge/, in the same commit as the code change, with frontmatter provenance and a log entry.
tags: [okf, docs, workflow]
status: stable
paths: [".knowledge/**"]
generated: { by: claude-code/claude-opus-5-5, at: 2026-09-26T08:21:10Z }
commit: 3b29ffe
sources:
  - id: okf-spec
    resource: https://github.com/GoogleCloudPlatform/open-knowledge-format/blob/main/SPEC.md
    title: Open Knowledge Format v0.2 specification
---

# When
- A change invalidates something a concept states (contract, decision,
  entry point, verification command).
- A decision is made that someone would otherwise relitigate -> new `Decision`.
- A module gets its first public function -> new `Module` concept.
- You had to read 5+ files to understand something -> capture it once.

Not when: the change is internal and breaks nothing a concept states. Most
changes need no knowledge edit.

# Steps
1. Edit the concept in the **same commit** as the code.
2. Frontmatter (OKF v0.2):[^okf-spec]
   ```yaml
   ---
   type: Module | Contract | Decision | Playbook | Phase   # required
   title: <display name>
   description: <one sentence; copied into index.md>
   resource: /src/biotapy/<path>/        # when it describes code
   paths: ["src/biotapy/<path>/**"]      # drives the staleness check
   tags: [..]
   status: draft | stable | deprecated
   generated: { by: <actor>, at: <UTC ISO 8601> }
   commit: <short sha the concept was checked against>
   verified: { by: human:pedrocr83, at: <UTC ISO 8601> }   # only when a human confirmed it
   sources: [{ id, resource, title }]
   ---
   ```
   Actors: `claude-code/<model>` for agents, `human:pedrocr83` for the maintainer.
3. Body: headings and lists; every behavioural claim anchored to
   `path:symbol`; links bundle-relative (`/contracts/x.md`).
4. Update the directory `index.md` line (title + description).
5. Add a dated entry to `.knowledge/log.md` (newest first).
6. Superseded decision: set `status: deprecated`, link to the replacement. Never delete.

# Verification
```bash
uv run --group test pytest tests/test_knowledge_bundle.py   # OKF conformance
bash scripts/knowledge_stale.sh --touched         # code changed under untouched concepts
```

# Common mistakes
- Restating code (signatures, line numbers) - it rots on the next commit.
- Leaving `commit` at an old SHA after re-checking a concept.
- Writing `verified` for agent-written content; only a human review earns it.

[^okf-spec]: Open Knowledge Format v0.2 specification
