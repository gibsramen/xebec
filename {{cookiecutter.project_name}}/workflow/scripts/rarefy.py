import biom
import numpy as np
import pandas as pd

from xebec import get_logger


xebec_logger = get_logger(snakemake.log[0], snakemake.rule)

table = biom.load_table(snakemake.input[0])
xebec_logger.info(f"Original shape: {table.shape}")
depths = table.sum(axis="sample")
rare_depth = round(np.quantile(
    depths,
    snakemake.config["rarefaction_depth_percentile"]/100)
)
xebec_logger.info(f"Rarefaction depth: {rare_depth}")
table_rare = table.subsample(rare_depth)
xebec_logger.info(f"Rarefied shape: {table_rare.shape}")

with biom.util.biom_open(snakemake.output[0], "w") as f:
    table_rare.to_hdf5(f, "rarefy")
