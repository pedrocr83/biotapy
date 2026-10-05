"""Plots of what tl, pp, fn and da give, such as the da.consensus table; pl computes no statistics of its own.

Boundaries: contracts/module-boundaries.
"""

from ._abundance import bar, contributions, heatmap
from ._consensus import consensus
from ._ordination import ordination, scree
from ._richness import richness

__all__ = ["bar", "consensus", "contributions", "heatmap", "ordination", "richness", "scree"]
