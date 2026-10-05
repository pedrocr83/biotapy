"""Where differential abundance methods agree: one row per feature over the result tables the user computed."""

from collections.abc import Sequence

import numpy as np
import pandas as pd

from ._schema import validate_result


def consensus(results: Sequence[pd.DataFrame], *, alpha: float = 0.05, min_methods: int = 2) -> pd.DataFrame:
    """Count, per feature, the methods that call it significant and whether they agree on its direction.

    Parameters
    ----------
    results
        Result tables of different ``bt.da`` methods run on the same data, ``group``
        and ``reference``, such as ``[bt.da.ancombc2(t, "host"), bt.da.linda(t, "host")]``.
    alpha
        A method calls a feature significant when its ``qvalue`` is below ``alpha``.
    min_methods
        How many methods must call a feature, all in the same direction, for consensus.

    Returns
    -------
    pandas.DataFrame
        One row per feature in any table (in first-seen order), indexed by ``feature``:
        ``effect_<method>``, ``qvalue_<method>`` and ``significant_<method>`` for each
        table in ``results`` order, then ``n_tested`` (methods with a p-value),
        ``n_significant``, ``direction`` (the sign the calling methods share; 0 when
        none calls the feature or they disagree), ``consensus`` and ``conflict``
        (methods call it in opposite directions).

    Raises
    ------
    TypeError
        ``results`` is a single table, or holds something that is not a table;
        ``alpha`` is not a real number or ``min_methods`` is not an ``int``.
    ValueError
        A table is not a ``bt.da`` result; two tables come from the same method or
        compare different contrasts; ``alpha`` is not between 0 and 1;
        ``min_methods`` is not between 1 and the number of tables; ``results`` is empty.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/differential_abundance`

    Calls are strict: ``qvalue == alpha`` is not significant. A feature a method
    did not test (NaN, or absent from its table) counts as not tested, never as
    not significant. A consensus feature has ``n_significant >= min_methods`` and
    every calling method has the same non-zero direction; a conflict is never a
    consensus. A call whose ``effect`` is exactly 0 has no direction: it counts in
    ``n_significant``, but ``direction`` is then 0, with no conflict and no
    consensus. The rule is recorded in the decision ``da-consensus-agreement``.

    Agreement between methods is a robustness report, not a way to choose a
    method: decide which methods to run before looking at their results.
    Methods that share a model (ANCOM-BC and ANCOM-BC2, say) agree more often for
    that reason alone.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()
    >>> table = bt.da.consensus([bt.da.ancombc2(tdata, "group"), bt.da.linda(tdata, "group")])
    >>> table.index[table["consensus"]].tolist()
    ['f6', 'f7']
    """
    if isinstance(results, pd.DataFrame):
        msg = "results must be a list of da result tables, such as [table_a, table_b]; got one table"
        raise TypeError(msg)
    tables = [validate_result(table, arg=f"results[{i}]") for i, table in enumerate(results)]
    methods = _check(tables, alpha=alpha, min_methods=min_methods)
    features = tables[0].index.append([table.index for table in tables[1:]]).unique()
    parts = [table.reindex(features) for table in tables]
    called = np.column_stack([(part["qvalue"] < alpha).to_numpy() for part in parts])
    signs = np.column_stack([np.sign(part["effect"].fillna(0)).to_numpy() for part in parts])
    up, down = (called & (signs > 0)).sum(axis=1), (called & (signs < 0)).sum(axis=1)
    n_significant = called.sum(axis=1)
    direction = np.where((up == n_significant) & (up > 0), 1, np.where((down == n_significant) & (down > 0), -1, 0))
    columns: dict[str, object] = {}
    for method, part, significant in zip(methods, parts, called.T, strict=True):
        columns |= {f"effect_{method}": part["effect"], f"qvalue_{method}": part["qvalue"]}
        columns[f"significant_{method}"] = significant
    out = pd.DataFrame(columns, index=features.rename("feature"))
    out["n_tested"] = np.column_stack([part["pvalue"].notna().to_numpy() for part in parts]).sum(axis=1)
    out["n_significant"] = n_significant
    out["direction"] = direction.astype(np.int8)
    out["consensus"] = (n_significant >= min_methods) & (direction != 0)
    out["conflict"] = (up > 0) & (down > 0)
    return out


def _check_options(n_tables: int, *, alpha: float, min_methods: int) -> None:
    """Raise unless ``alpha`` is a real number in (0, 1) and ``min_methods`` an int in [1, n_tables]."""
    if n_tables == 0:
        msg = "results is empty; pass at least one method's table"
        raise ValueError(msg)
    if isinstance(alpha, bool) or not isinstance(alpha, int | float | np.integer | np.floating):
        msg = f"alpha must be a real number, got {type(alpha).__name__}"
        raise TypeError(msg)
    if isinstance(min_methods, bool) or not isinstance(min_methods, int | np.integer):
        msg = f"min_methods must be an int, got {type(min_methods).__name__}"
        raise TypeError(msg)
    if not 0 < alpha < 1:
        msg = f"alpha must be between 0 and 1, got {alpha}"
        raise ValueError(msg)
    if not 1 <= min_methods <= n_tables:
        msg = f"min_methods must be between 1 and the number of results ({n_tables}), got {min_methods}"
        raise ValueError(msg)


def _check(tables: list[pd.DataFrame], *, alpha: float, min_methods: int) -> list[str]:
    """The tables' method names, after checking they can be compared under ``alpha`` and ``min_methods``."""
    _check_options(len(tables), alpha=alpha, min_methods=min_methods)
    methods = [str(table["method"].iloc[0]) for table in tables]
    repeated = sorted({method for method in methods if methods.count(method) > 1})
    if repeated:
        msg = f"results repeat the method(s) {repeated}; pass one table per method"
        raise ValueError(msg)
    contrasts = sorted({str(table["contrast"].iloc[0]) for table in tables})
    if len(contrasts) > 1:
        msg = f"results compare different contrasts {contrasts}; run every method with the same group and reference"
        raise ValueError(msg)
    return methods
