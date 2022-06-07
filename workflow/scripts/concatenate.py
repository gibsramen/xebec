from pathlib import PurePath

import pandas as pd

from xebec import get_logger
from xebec.src._concatenate import concatenate_metric_dataframes


xebec_logger = get_logger(snakemake.log[0], snakemake.rule)
xebec_logger.info("Concatenating...")
all_metrics_df = concatenate_metric_dataframes(snakemake.input)
xebec_logger.info("Finished concatenating!")
xebec_logger.info(f"\n{all_metrics_df.head()}")
all_metrics_df.to_csv(snakemake.output[0], sep="\t", index=False)
xebec_logger.info(f"Saved to {snakemake.output[0]}!")
