"""Write the Coming-from-R table from the public docstrings' ``R equivalent:`` lines (rules.md R8.4)."""

import importlib
import re
import tomllib
from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sphinx.application import Sphinx

# contracts/function-shape: exactly one such line, holding ``pkg::fn`` items separated by ", ", or "none".
R_LINE = re.compile(r"^\s*R equivalent:(?P<rest>.*)$", re.MULTILINE)
R_ITEMS = re.compile(r"none|``[\w.]+::[\w.]+``(?:, ``[\w.]+::[\w.]+``)*")
R_ITEM = re.compile(r"``([\w.]+::[\w.]+)``")
IDIOMS = Path(__file__).resolve().parents[1] / "_data" / "r_idioms.toml"
# Relative to the docs source directory; docs/generated/ is git-ignored and left out of the sdist.
TABLE = Path("generated") / "coming_from_r_table.md"


def public_functions() -> list[tuple[str, Callable[..., object]]]:
    """``("bt.<module>.<name>", function)`` for every name in every public subpackage's ``__all__``."""
    biotapy = importlib.import_module("biotapy")
    modules = [importlib.import_module(f"biotapy.{name}") for name in biotapy.__all__ if name != "__version__"]
    return [
        (f"bt.{module.__name__.split('.')[-1]}.{name}", getattr(module, name))
        for module in modules
        for name in module.__all__
    ]


def r_equivalents(doc: str) -> list[str]:
    """The ``pkg::fn`` items of the docstring's one ``R equivalent:`` line; ``[]`` when it says ``none``."""
    lines = R_LINE.findall(doc)
    if len(lines) != 1 or R_ITEMS.fullmatch(lines[0].strip()) is None:
        msg = f"need exactly one line 'R equivalent: ``pkg::fn``, ...' or 'R equivalent: none', found {lines}"
        raise ValueError(msg)
    return R_ITEM.findall(lines[0])


def rows() -> dict[str, list[str]]:
    """R call -> the biotapy cells for it: the docstrings' functions, then the idioms file."""
    table: dict[str, list[str]] = {}
    for name, function in public_functions():
        for r_name in r_equivalents(function.__doc__ or ""):
            table.setdefault(r_name, []).append(f"{{func}}`{name} <biotapy.{name.removeprefix('bt.')}>`")
    idioms: dict[str, str] = tomllib.loads(IDIOMS.read_text(encoding="utf-8"))["idioms"]
    both = sorted(set(idioms) & set(table))
    if both:
        msg = f"{both} are mapped by a docstring and by {IDIOMS.name}; keep only the docstring"
        raise ValueError(msg)
    return table | {r_name: [python] for r_name, python in idioms.items()}


def render() -> str:
    """The table as MyST markdown, one row per R call, sorted by package and name."""
    table = rows()
    lines = ["| R | biotapy |", "|---|---|"]
    lines += [f"| `{r_name}` | {', '.join(table[r_name])} |" for r_name in sorted(table, key=str.lower)]
    return "\n".join(lines) + "\n"


def write_table(app: "Sphinx") -> None:
    """``builder-inited`` hook: write the table before Sphinx reads any source."""
    path = Path(app.srcdir) / TABLE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render(), encoding="utf-8")


def setup(app: "Sphinx") -> dict[str, bool]:
    """Register the hook."""
    app.connect("builder-inited", write_table)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
