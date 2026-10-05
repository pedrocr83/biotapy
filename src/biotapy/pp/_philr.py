"""PhILR: isometric log-ratio balances over the phylogeny."""

import pandas as pd
from skbio import TreeNode
from skbio.stats.composition import clr, tree_basis

from biotapy._core import TreeData, add_provenance, get_skbio_tree

from ._transform import pseudocounted


def philr(tdata: TreeData, *, pseudocount: float = 0.5) -> TreeData:
    """Add the phylogenetic isometric log-ratio transform (PhILR) as ``obsm['X_philr']``.

    Parameters
    ----------
    tdata
        Samples x features with a rooted binary phylogeny in ``vart['phylo']``;
        ``X`` holds counts or another non-negative abundance.
    pseudocount
        Added to every value of ``X`` before the logarithm, so zeros have one.

    Returns
    -------
    TreeData
        A copy of ``tdata`` with ``obsm['X_philr']``: a samples x balances
        ``pandas.DataFrame`` with one column per internal node of the tree, named
        after the node, in preorder. A balance is positive when its node's first
        child is more abundant than its second.

    Raises
    ------
    TypeError
        ``tdata`` is an AnnData that is not a TreeData.
    KeyError
        ``tdata`` has no ``vart['phylo']``.
    ValueError
        A node of the tree, the root included, has more than two children; a
        feature is not a tip of the tree; there are fewer than two features; or
        ``pseudocount`` or ``X`` is invalid, as in :func:`biotapy.pp.clr`.

    Warns
    -----
    UserWarning
        ``pseudocount`` is larger than the smallest non-zero value in ``X``.

    Notes
    -----
    R equivalent: ``philr::philr``, ``mia::transformAssay``
    Guide: :doc:`/guide/transforms`

    Equals ``philr::philr(x, tree, pseudocount = 0.5)`` with its default uniform part and
    ILR weights, which ``mia::transformAssay(method = "philr")`` calls; the weights are
    not offered. Each balance of node ``i`` is
    ``sqrt(r s / (r + s)) * log(g(first child) / g(second child))``, where ``r`` and
    ``s`` count the taxa under each child and ``g`` is their geometric mean.
    Subsetting a TreeData keeps one-child nodes; they define no balance and are
    skipped, as ``ape::drop.tip`` removes them in R. A wider node has no single
    balance: an unrooted tree's three-child root must be rooted, and polytomies
    resolved, before PhILR (for example with ``ape::multi2di``), as ``philr``
    requires. Balances are Euclidean coordinates: distances between samples equal
    Aitchison distances, whatever the tree.

    ``X`` is densified once (8 bytes x samples x features); the balances take
    8 bytes x samples x (features - 1).

    References
    ----------
    Silverman JD, Washburne AD, Mukherjee S, David LA (2017) A phylogenetic transform enhances
    analysis of compositional microbiota data. eLife 6:e21887.

    Examples
    --------
    >>> import biotapy as bt
    >>> tdata = bt.datasets.toy()[:, :6].copy()  # toy()'s root has three children; f1-f6 sit under two
    >>> out = bt.pp.philr(tdata)
    >>> out.obsm["X_philr"].columns.tolist()
    ['root', 'n1', 'n4', 'n2', 'n5']
    """
    if tdata.n_vars < 2:
        msg = f"pp.philr needs at least two features, got {tdata.n_vars}"
        raise ValueError(msg)
    tree = _binary_tree(get_skbio_tree(tdata, split_root=False))
    tips = pd.Index([tip.name for tip in tree.tips()])
    if len(tips) != tdata.n_vars:
        msg = f"pp.philr needs every feature to be a tip of the tree; {tdata.n_vars - len(tips)} feature(s) are not"
        raise ValueError(msg)
    values = pseudocounted(tdata, pseudocount, func="pp.philr").take(tdata.var_names.get_indexer(tips), axis=1)
    basis, nodes = tree_basis(tree)
    # tree_basis puts a node's first child in the denominator; philr::philr puts it in the numerator.
    balances = pd.DataFrame(-(clr(values) @ basis.T), index=tdata.obs_names, columns=nodes)
    out = tdata.copy()
    out.obsm["X_philr"] = balances[[node.name for node in tree.preorder() if not node.is_tip()]]
    add_provenance(out, "pp.philr", pseudocount=pseudocount)
    return out


def _binary_tree(tree: TreeNode) -> TreeNode:
    """A copy of ``tree`` without one-child nodes, children in order; raise if a node has more than two."""
    kept: dict[int, TreeNode] = {}
    wide: list[str] = []
    for node in tree.postorder(include_self=True):
        children = [kept[id(child)] for child in node.children]
        if len(children) > 2:
            wide.append(str(node.name))
        # TreeNode.prune would also drop one-child nodes, but it moves the child to the end of its parent's list.
        kept[id(node)] = children[0] if len(children) == 1 else TreeNode(node.name, children=children)
    if wide:
        msg = (
            f"pp.philr needs a rooted binary tree, but the node(s) {wide[:3]} have more than two children; "
            "root the tree and resolve its polytomies first (for example with ape::multi2di in R)"
        )
        raise ValueError(msg)
    return kept[id(tree)]
