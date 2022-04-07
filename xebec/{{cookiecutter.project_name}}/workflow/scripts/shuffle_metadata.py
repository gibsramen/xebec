import os

import numpy as np
import pandas as pd


os.makedirs(snakemake.params[0], exist_ok=True)
md = pd.read_table(snakemake.input[0], sep="\t", index_col=0)

# Generate n shuffled metadata files
n = snakemake.config["shuffle_iterations"]

for i in range(n):
    shuf_cols = []
    for col in md.columns:
        shuf = md[col].sample(frac=1, replace=False).reset_index(drop=True)
        shuf_cols.append(shuf)
    shuf_df = pd.concat(shuf_cols, axis=1)
    shuf_df.index = md.index
    shuf_df.to_csv(snakemake.output[i], sep="\t", index=True)
