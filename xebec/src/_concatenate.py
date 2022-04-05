from pathlib import PurePath

import pandas as pd


def concatenate_metric_dataframes(files):
    """Concatenate results from multiple metrics."""

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


def get_metric_info(f):
    """Return metric type and metric name as tuple.

    results/*_div/*/*/effect_sizes.tsv
    """
    path_parts = PurePath(f).parts
    return path_parts[2], path_parts[3]
