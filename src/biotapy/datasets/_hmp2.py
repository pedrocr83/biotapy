"""The HMP2 (IBDMDB) inflammatory bowel disease cohort: pathway abundance and taxa, one stool metagenome per person."""

from typing import cast

import pandas as pd
from anndata import AnnData
from mudata import MuData

from biotapy.io import read_humann, read_metaphlan

from ._remote import _fetch

TAXA_KEY = "taxa"
# Metadata columns kept, as hmp2_metadata_2018-08-20.csv names them; the other 483 are left out.
COLUMNS = ["Participant ID", "week_num", "diagnosis", "site_name", "sex", "consent_age", "Antibiotics"]
# A participant's metagenomes, earliest first; visit_num orders them within a week and is not kept.
ORDER = ["Participant ID", "week_num", "visit_num", "External ID"]
DIAGNOSES = ["nonIBD", "UC", "CD"]


def hmp2() -> MuData:
    """The HMP2 inflammatory bowel disease cohort: each participant's first stool metagenome.

    Downloaded once (23 MB) from the IBDMDB and cached: HUMAnN 3 pathway
    abundance, MetaPhlAn 3 profiles and the sample metadata of the HMP2
    metagenomes (1,638 samples from 130 participants), of which the first
    sample of each participant is kept.

    Returns
    -------
    MuData
        130 samples (65 CD, 38 UC, 27 nonIBD), indexed by the metadata's
        ``External ID``, in three modalities:

        - ``"function"``: community pathway abundance (``x_kind == "cpm"``),
          as ``bt.io.read_humann`` reads it;
        - ``"function_by_taxon"``: the same pathways per species;
        - ``"taxa"``: MetaPhlAn 3 species (``x_kind == "relative"``), as
          ``bt.io.read_metaphlan`` reads it, with no tree.

        The global ``obs``, copied to every modality, holds the metadata
        columns ``Participant ID``, ``week_num``, ``diagnosis`` (categorical:
        ``nonIBD``, ``UC``, ``CD``), ``site_name``, ``sex``, ``consent_age``
        and ``Antibiotics``. Every feature of the published tables is kept,
        including those absent from all 130 samples.

    Notes
    -----
    R equivalent: none
    Guide: :doc:`/guide/datasets`

    A participant's first sample is the metagenome with the lowest
    ``week_num``; equal weeks are ordered by ``visit_num`` (a missing one
    last), then by ``External ID``. One sample per
    person keeps samples independent, so group comparisons such as
    ``bt.tl.permanova`` do not count one person several times.

    Reading the pathway table builds one dense 22,113 x 1,638 ``float64``
    array (about 290 MB) before the samples are selected.

    The IBDMDB states no licence for these files. biotapy ships none of
    them; cite the study when you use them.

    References
    ----------
    Lloyd-Price J et al. (2019) Multi-omics of the gut microbial ecosystem in inflammatory
    bowel diseases. Nature 569:655-662.

    Examples
    --------
    >>> import biotapy as bt
    >>> mdata = bt.datasets.hmp2()  # doctest: +SKIP
    >>> mdata["function"].shape, mdata["taxa"].shape  # doctest: +SKIP
    ((130, 478), (130, 579))
    """
    metadata = pd.read_csv(_fetch("hmp2_metadata_2018-08-20.csv"), usecols=["data_type", *COLUMNS, *ORDER])
    metagenomes = metadata[metadata["data_type"] == "metagenomics"]
    first = metagenomes.sort_values(ORDER).drop_duplicates("Participant ID")
    obs = first.set_index("External ID")[COLUMNS].rename_axis(None)
    obs["diagnosis"] = pd.Categorical(obs["diagnosis"], categories=DIAGNOSES)
    pathways = read_humann(_fetch("pathabundances_3.tsv.gz"))
    taxa = read_metaphlan(_fetch("taxonomic_profiles_3.tsv.gz"))
    samples = obs.index
    # MuData types a modality as AnnData | MuData; read_humann's two are AnnData.
    function = cast("dict[str, AnnData]", pathways.mod)
    modalities = {key: mod[samples].copy() for key, mod in function.items()} | {TAXA_KEY: taxa[samples].copy()}
    mdata = MuData(modalities, obs=obs)
    mdata.push_obs()
    return mdata
