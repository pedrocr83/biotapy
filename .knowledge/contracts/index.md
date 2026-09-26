# Contracts

* [Public function shape](function-shape.md) - One task = one public function `verb(data, required, *, options) -> result`, fully typed, keyword-only options, seeded randomness, NumPy docstring with a parseable R-equivalent line.
* [Data-model slots](data-model-slots.md) - Which AnnData/TreeData slot holds what, the exact result keys, and which slots feature-changing operations drop.
* [Module boundaries](module-boundaries.md) - Layered package (_core at the bottom, pl/ml/da at the top); public API only through subpackage __all__; private topic files; shared helpers only in _core.
* [Tree access](tree-access.md) - Only biotapy/_core/_tree.py touches the TreeData tree API, so a TreeData change touches one file.
* [Engine parity](engine-parity.md) - A compiled kernel is a drop-in behind `engine=`; the Python engine is the oracle and both must match within a stated tolerance.
* [R golden parity](r-golden-parity.md) - Every function with an R equivalent is tested against golden files exported from pinned R; deterministic outputs match numerically, stochastic outputs match invariants.
