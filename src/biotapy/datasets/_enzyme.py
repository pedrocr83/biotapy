"""The ENZYME (EC) hierarchy from the SIB Swiss Institute of Bioinformatics, CC BY 4.0."""

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ._remote import _fetch

LEVELS = ("class", "subclass", "subsubclass")
# enzclass.txt lines: "1. 1. 1.-    With NAD(+) or NADP(+) as acceptor."
_CLASS_LINE = re.compile(r"^(\d+\.\s*[\d-]+\.\s*[\d-]+\.-)\s+(.*?)\.?\s*$", re.MULTILINE)
_ID_LINE = re.compile(r"^ID   (\S+)$", re.MULTILINE)
_RELEASE = re.compile(r"^CC   Release of (.+)$", re.MULTILINE)
_CLASS_RELEASE = re.compile(r"^Release:\s+(\S+)", re.MULTILINE)
_URL = "https://enzyme.expasy.org/"


def enzyme() -> pd.DataFrame:
    """The ENZYME EC hierarchy as an edge table for ``bt.fn.func_glom``.

    Downloaded once (9.6 MB) from ``ftp.expasy.org`` and cached.

    Returns
    -------
    pandas.DataFrame
        One row per EC number and ancestor: columns ``child`` (an EC number
        such as ``1.1.1.1``, or an internal id such as ``1.1.1.-``),
        ``parent`` (its ancestor), ``level`` (the parent's level:
        ``"class"``, ``"subclass"`` or ``"subsubclass"``) and
        ``parent_name`` (from ``enzclass.txt``; NaN when ENZYME names none).
        ``attrs["source"]`` names the ENZYME release read and
        ``attrs["license"]`` is ``"CC BY 4.0"``.

    Raises
    ------
    ValueError
        If ``enzyme.dat`` has no ``Release of`` line, ``enzclass.txt`` has no
        ``Release:`` line or no class lines, or the two releases differ.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/datasets`

    Every ancestor is listed, not only the direct parent, so a table of EC
    numbers and one already grouped to sub-subclasses both reach any level.
    Deleted and transferred entries keep their number and so their place.

    ENZYME keeps only its current release online, so the first download is
    cached for good and ``attrs["source"]`` records which release it was.
    To take a newer one, delete the cached ``enzyme.dat`` and ``enzclass.txt``
    from the ``BIOTAPY_DATA_DIR`` directory if that variable is set, otherwise
    from pooch's per-user cache directory (``pooch.os_cache("biotapy")``,
    which depends on the platform); delete both, or the two files may come from
    different releases.

    ENZYME is copyrighted by the SIB Swiss Institute of Bioinformatics and
    distributed under the Creative Commons Attribution 4.0 (CC BY 4.0)
    licence; cite it when you publish results that use it.

    References
    ----------
    Bairoch A (2000) The ENZYME database in 2000. Nucleic Acids Res 28:304-305.

    Examples
    --------
    >>> import biotapy as bt
    >>> edges = bt.datasets.enzyme()  # doctest: +SKIP
    >>> edges.query("child == '1.1.1.1'")["parent"].tolist()  # doctest: +SKIP
    ['1.-.-.-', '1.1.-.-', '1.1.1.-']
    """
    entries_path = Path(_fetch("enzyme.dat"))
    entries = entries_path.read_text(encoding="utf-8")
    release = _RELEASE.search(entries)
    if release is None:
        raise ValueError(f"{entries_path.name} has no 'CC   Release of ...' line, so its release is unknown")
    classes_path = Path(_fetch("enzclass.txt"))
    classes = classes_path.read_text(encoding="utf-8")
    class_release = _CLASS_RELEASE.search(classes)
    if class_release is None:
        raise ValueError(f"{classes_path.name} has no 'Release: ...' line, so its release is unknown")
    if class_release[1] != release[1]:
        raise ValueError(
            f"{entries_path.name} is release {release[1]} but {classes_path.name} is release {class_release[1]}; "
            "delete both cached files and download them again"
        )
    pairs = _CLASS_LINE.findall(classes)
    if not pairs:
        raise ValueError(f"{classes_path.name} has no class lines such as '1. 1. 1.-    Name.'")
    names = {re.sub(r"\s", "", ec): name for ec, name in pairs}
    leaves = pd.Series(_ID_LINE.findall(entries), dtype=str)
    edges = _ancestors(leaves)
    internal = pd.Series(sorted(set(names) | set(edges["parent"])), dtype=str)
    edges = pd.concat([edges, _ancestors(internal)], ignore_index=True)
    edges["parent_name"] = edges["parent"].map(names).astype(pd.StringDtype(na_value=np.nan))
    edges.attrs["source"] = f"ENZYME release {release[1]}, SIB, {_URL}"
    edges.attrs["license"] = "CC BY 4.0"
    return edges


def _ancestors(ids: pd.Series) -> pd.DataFrame:
    """Every (id, ancestor, ancestor's level) row; ``1.1.1.1`` has three, ``1.-.-.-`` none."""
    depth = 4 - ids.str.count(r"\.-")
    parts = ids.str.split(".")
    frames = [
        pd.DataFrame(
            {"child": ids[depth > k], "parent": parts[depth > k].str[:k].str.join(".") + ".-" * (4 - k), "level": level}
        )
        for k, level in enumerate(LEVELS, start=1)
    ]
    return pd.concat(frames, ignore_index=True)
