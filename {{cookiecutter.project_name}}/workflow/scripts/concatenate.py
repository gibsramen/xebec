from pathlib import PurePath

import pandas as pd

from helper import get_logger


def concatenate_metric_dataframes(files):
    """Concatenate results from multiple metrics."""
    def get_metric_info(f):
        """Return metric type and metric name as tuple.

        results/*_div/*/*/effect_sizes.tsv
        """
        path_parts = PurePath(f).parts
        return path_parts[2], path_parts[3]

    all_dfs = []
    all_keys = []
    for f in files:
        this_df = pd.read_table(f, sep="\t")
        this_keys = get_metric_info(f)
        all_dfs.append(this_df)
        all_keys.append(this_keys)

    total_df = pd.concat(
        all_dfs,
        keys=all_keys,
        names=["phylogenetic", "diversity_metric"]
    )
    total_df = total_df.reset_index(level=("phylogenetic", "diversity_metric"))
    return total_df


xebec_logger = get_logger(snakemake.log[0], snakemake.rule)
xebec_logger.info("Concatenating...")
all_metrics_df = concatenate_metric_dataframes(snakemake.input)
xebec_logger.info("Finished concatenating!")
xebec_logger.info(f"\n{all_metrics_df.head()}")
all_metrics_df.to_csv(snakemake.output[0], sep="\t", index=False)
xebec_logger.info(f"Saved to {snakemake.output[0]}!")
