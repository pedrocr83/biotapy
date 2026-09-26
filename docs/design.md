# Design

biotapy keeps its design knowledge next to the code, in the
[`.knowledge/` bundle](https://github.com/pedrocr83/biotapy/tree/master/.knowledge),
written in the Open Knowledge Format (OKF v0.2): plain markdown files with YAML frontmatter.

- **Decisions** record why things are the way they are: TreeData as the only container,
  samples as rows, pure functions by default, Python first and compiled code last.
- **Contracts** state what every function must keep true: its shape, the data-model slots,
  module boundaries, and parity with R.
- **Roadmap** lists each release as a phase, with its tasks and exit gate.

Contributors and coding agents read it before changing code; see the [contributing guide](contributing.md).
