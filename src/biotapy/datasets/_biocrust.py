"""A desert biocrust wetting experiment: microbes and metabolites over the same samples, from mmvec's repository."""

from mudata import MuData

from biotapy.io import read_biom, to_mudata

from ._remote import _fetch

# Fetched together, before either is parsed, so a download error comes first.
FILES = ["biocrust_microbes.biom", "biocrust_metabolites.biom"]


def biocrust() -> MuData:
    """Microbes and metabolites of a desert biological soil crust after wetting, over the samples both tables have.

    Downloaded once (135 KB) from the example of mmvec's repository
    (``examples/soils``), pinned to one commit: 466 microbes counted in 20
    samples and 85 metabolites measured in 19 of them, taken at five times
    after wetting.

    Returns
    -------
    MuData
        The 19 samples both tables have, in the microbe table's order, in
        two modalities: ``"taxa"``, the microbes' counts
        (``x_kind == "counts"``, no taxonomy columns and no tree; each
        feature name ends with its phylum), and ``"metabolites"``, the
        metabolites' intensities (``x_kind == "abundance"``), both as
        ``bt.io.read_biom`` reads them. The tables carry no sample metadata;
        a sample's name holds its time after wetting and its position.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/datasets`

    The microbe table's sample ``9hr_late`` has no metabolite profile and
    is left out without a warning; the files themselves are unchanged.

    The files are distributed in mmvec's repository under its BSD-3-Clause
    licence. biotapy ships none of them; cite the mmvec paper when you use
    them.

    References
    ----------
    Morton JT et al. (2019) Learning representations of microbe-metabolite interactions.
    Nature Methods 16:1306-1314.

    Examples
    --------
    >>> import biotapy as bt
    >>> mdata = bt.datasets.biocrust()  # doctest: +SKIP
    >>> mdata["taxa"].shape, mdata["metabolites"].shape  # doctest: +SKIP
    ((19, 466), (19, 85))
    """
    paths = [_fetch(name) for name in FILES]
    microbes, metabolites = (read_biom(path) for path in paths)
    shared = microbes.obs_names.intersection(metabolites.obs_names, sort=False)
    return to_mudata({"taxa": microbes[shared], "metabolites": metabolites[shared]})
