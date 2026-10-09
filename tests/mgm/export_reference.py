"""Write tests/data/mgm/: a small genus table and MGM 0.5.8's own embedding of each of its samples.

microformer-mgm pins numpy 1.24, pandas 2.0, torch 2.0 and transformers 4.33, so it cannot be installed
beside biotapy. Run from the repository root, never in CI:

    uv run --no-project --python 3.11 --with microformer-mgm==0.5.8 --with torch==2.0.1+cpu \
        --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match \
        python tests/mgm/export_reference.py

The table is synthetic; its genus names are MGM's vocabulary, plus one name MGM lacks and one feature with no
genus. MGM's own code builds the tokens (`MicroCorpus`) and runs the model (`GPT2LMHeadModel`, transformers 4.33),
one sample at a time as MGM's notebook does. The embedding is the mean of the last hidden layer over the sample's
tokens (<bos>, its genera, <eos>), the "element-wise mean pooling" MGM's paper uses for the pretrained model
(Methods 4.5). MGM drops a sample with no count in its vocabulary; such a sample is embedded here from the tokens
<bos> <eos>, as biotapy does.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

# MGM, torch and transformers are imported inside the functions: pytest imports this file for doctests in every run.
OUT = Path("tests/data/mgm")
GUT = ["Bacteroides", "Prevotella", "Faecalibacterium", "Bifidobacterium", "Akkermansia", "Blautia", "Roseburia",
       "Alistipes", "Streptococcus", "Lactobacillus", "Veillonella", "Ruminococcus"]  # fmt: skip


def table() -> pd.DataFrame:
    """Features x (genus, s1..s7), counts."""
    from mgm.CLI.CLI_utils import find_pkg_resource

    vocabulary = pd.read_csv(find_pkg_resource("resources/phylogeny.csv"), index_col=0).index
    plain = [name[3:] for name in vocabulary if re.fullmatch(r"g__[A-Za-z0-9_]+", name) and name[3:] not in GUT]
    many = plain[::16][:600]
    genera = [*GUT, "Escherichia", "Escherichia-Shigella", "Notagenus", None, *many]
    counts = np.zeros((len(genera), 7), dtype=np.int64)
    gut = np.arange(1, len(GUT) + 1)
    counts[: len(GUT), 0] = gut * 10  # s1: twelve genera
    counts[: len(GUT), 1] = gut[::-1] * 7  # s2: the same, other counts, and the rest below
    counts[12:16, 1] = [3, 5, 40, 25]  # Escherichia twice (one token), a genus MGM lacks, no genus
    counts[2, 2] = 9  # s3: one genus
    counts[14:16, 3] = [8, 2]  # s4: only a genus MGM lacks and no genus; s5: all zero
    counts[16:, 5] = (np.arange(600) * 7919) % 600 + 1  # s6: 600 genera, more than the 510 tokens a sample can hold
    counts[[0, 3, 20, 400], 6] = [500, 1, 1, 30]  # s7: few genera, two of them once
    frame = pd.DataFrame(counts, columns=[f"s{i}" for i in range(1, 8)], index=[f"f{i}" for i in range(len(genera))])
    frame.insert(0, "genus", genera)
    frame.index.name = "feature"
    return frame


def embed(input_ids, attention_mask, model) -> np.ndarray:
    """MGM's notebook `cal_embed` on CPU, with the mean over the sample's tokens instead of its last token."""
    model.eval()
    hidden = model(input_ids=input_ids, attention_mask=attention_mask, output_hidden_states=True).hidden_states[-1]
    return hidden.squeeze(0)[attention_mask.squeeze(0) == 1].mean(0).detach().numpy()


def main() -> None:
    import torch
    from mgm.CLI.CLI_utils import find_pkg_resource
    from mgm.src.MicroCorpus import MicroCorpus
    from mgm.src.utils import CustomUnpickler
    from transformers import GPT2LMHeadModel

    frame = table()
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / "counts.csv")
    columns = ["g__" + genus if isinstance(genus, str) else "unclassified" for genus in frame["genus"]]
    abundance = pd.DataFrame(frame.drop(columns="genus").T.to_numpy(), index=frame.columns[1:], columns=columns)
    with open(find_pkg_resource("resources/MicroTokenizer.pkl"), "rb") as file:
        tokenizer = CustomUnpickler(file).load()
    corpus = MicroCorpus(abu=abundance.astype(float), tokenizer=tokenizer, max_len=512, preprocess=True)
    model = GPT2LMHeadModel.from_pretrained(find_pkg_resource("resources/general_model"))
    kept = list(corpus.data.index)
    empty = torch.tensor([[2, 3] + [0] * 510]), torch.tensor([[1.0, 1.0] + [0.0] * 510])
    rows = {}
    with torch.no_grad():
        for sample in abundance.index:
            if sample in kept:
                item = corpus[kept.index(sample)]
                rows[sample] = embed(item["input_ids"][None], item["attention_mask"][None], model)
            else:
                rows[sample] = embed(*empty, model)
    embeddings = pd.DataFrame(rows).T
    embeddings.columns = [f"e{i}" for i in range(embeddings.shape[1])]
    embeddings.index.name = "sample"
    embeddings.to_csv(OUT / "embeddings.csv", float_format="%.9g")


if __name__ == "__main__":
    main()
