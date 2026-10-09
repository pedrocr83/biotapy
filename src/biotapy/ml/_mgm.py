"""MGM, the reference embedding plugin: entry point ``mgm`` in ``biotapy.embeddings``, extra ``mgm``."""

from collections.abc import Iterable
from pathlib import Path
from typing import Any, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import pooch
from anndata import AnnData

from biotapy._core import as_csr, divide_rows, finite_non_negative, import_optional, make_pooch, sum_by, warn_user

# microformer-mgm 0.5.8's wheel (MIT) carries the pretrained general model; files.pythonhosted.org never changes a file.
_WHEEL = "microformer_mgm-0.5.8-py3-none-any.whl"
_BASE_URL = (
    "https://files.pythonhosted.org/packages/4b/9c/829a1e59d5e618756ce8e57553e15d13bb16c58d37896f03233d8697515a/"
)
_SHA256 = "sha256:210891685565022ea869e88a7452769dc6fbe3f699d2a133e15338d1e30eb92b"
_MEMBERS = [
    "mgm/resources/general_model/config.json",
    "mgm/resources/general_model/pytorch_model.bin",
    "mgm/resources/phylogeny.csv",
]
# How MGM reads a genus from a column name (mgm/src/MicroCorpus.py, MicroCorpus._preprocess).
_TOKEN = r"(g__[A-Za-z0-9_]+)"
# MGM's vocabulary: <pad>, <mask>, <bos>, <eos>, then phylogeny.csv's genera in file order (mgm/resources/MicroTokenizer.pkl).
_BOS, _EOS, _FIRST_GENUS = 2, 3, 4
_MAX_TOKENS = 512


def embed(adata: AnnData) -> npt.NDArray[np.float32]:
    """MGM's embedding of each sample, the mean of its last hidden layer over the sample's tokens (see ``ml.embed``)."""
    # Checked before torch, transformers or the download, so these errors need no extra (and are tested in every job).
    if "genus" not in adata.var.columns:
        msg = "mgm reads var['genus'], which this table does not have"
        raise KeyError(msg)
    if not finite_non_negative(as_csr(adata.X)):
        msg = "mgm reads counts or relative abundances; adata.X holds negative or non-finite values"
        raise ValueError(msg)
    torch: Any = import_optional("torch", extra="mgm")
    transformers: Any = import_optional("transformers", extra="mgm")
    paths = make_pooch(_BASE_URL, {_WHEEL: _SHA256}).fetch(_WHEEL, processor=pooch.Unzip(members=_MEMBERS))
    files = {Path(path).name: Path(path) for path in paths}  # Unzip lists them in os.walk order
    # anndata types var columns as Series | DataArray (its lazy variant); the data model guarantees a Series.
    genus = cast("pd.Series[str]", adata.var["genus"])
    sentences = _sentences(adata, genus, pd.read_csv(files["phylogeny.csv"], index_col=0))
    model = transformers.GPT2Model(transformers.GPT2Config.from_json_file(files["config.json"]))
    weights = torch.load(files["pytorch_model.bin"], map_location="cpu", weights_only=True)
    # The file holds a GPT2LMHeadModel: the body under "transformer.", and the head, which an embedding does not use.
    model.load_state_dict(
        {key.removeprefix("transformer."): value for key, value in weights.items() if key != "lm_head.weight"}
    )
    model.eval()
    out = np.empty((len(sentences), model.config.n_embd), dtype=np.float32)
    # One sample at a time, unpadded: on the CPU, batches were slower and grew memory (phase-4 plan, task 4.4b).
    with torch.inference_mode():
        for row, sentence in enumerate(sentences):
            out[row] = model(input_ids=torch.from_numpy(sentence)[None]).last_hidden_state[0].mean(0).numpy()
    return out


def _sentences(adata: AnnData, genus: "pd.Series[str]", phylogeny: pd.DataFrame) -> list[npt.NDArray[np.int64]]:
    """Each sample's tokens as MGM builds them: <bos>, its genera by standardised abundance, <eos>, cut to 512."""
    tokens = ("g__" + genus.astype("string")).str.extract(_TOKEN, expand=False)
    codes = phylogeny.index.get_indexer(pd.Index(tokens))
    missing = genus.isna().to_numpy()
    unknown = (codes < 0) & ~missing
    if missing.any() or unknown.any():
        msg = (
            f"mgm leaves out {int(missing.sum() + unknown.sum())} of {codes.size} features: {int(missing.sum())} "
            f"without a genus and {int(unknown.sum())} whose genus is not one of MGM's ({_listed(genus[unknown])})"
        )
        warn_user(msg)
    # Features of one genus are summed, and relative abundance is taken over MGM's genera only, as MGM does.
    grouped = sum_by(as_csr(adata.X), codes, len(phylogeny))
    grouped.sort_indices()
    relative = divide_rows(grouped, np.asarray(grouped.sum(axis=1), dtype=np.float64).ravel())
    mean, std = phylogeny["mean"].to_numpy(), phylogeny["std"].to_numpy()
    sentences = []
    for row in range(relative.shape[0]):
        span = slice(relative.indptr[row], relative.indptr[row + 1])
        genera = relative.indices[span]
        standardised = (relative.data[span] - mean[genera]) / std[genera]
        # MGM keeps a genus whose standardised value is above that of zero abundance, and sorts with pandas.
        kept = pd.Series(standardised, index=genera)[standardised > (0 - mean[genera]) / std[genera]]
        order = kept.sort_values(ascending=False).index.to_numpy()
        sentences.append(np.r_[_BOS, order + _FIRST_GENUS, _EOS][:_MAX_TOKENS].astype(np.int64))
    empty = [name for name, sentence in zip(adata.obs_names, sentences, strict=True) if sentence.size == 2]
    if empty:
        warn_user(
            f"mgm embeds {len(empty)} sample(s) from <bos> <eos> alone, none of their genera being MGM's: {_listed(empty)}"
        )
    return sentences


def _listed(names: "Iterable[str]") -> str:
    """Up to five distinct names, sorted, for a warning."""
    distinct = sorted({str(name) for name in names})
    return ", ".join(distinct[:5]) + (", ..." if len(distinct) > 5 else "")
