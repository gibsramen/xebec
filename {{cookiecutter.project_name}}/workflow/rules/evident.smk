from pathlib import PurePath

from evident import AlphaDiversityHandler, BetaDiversityHandler
from evident.effect_size import effect_size_by_category, pairwise_effect_size_by_category
import numpy as np
import pandas as pd
from skbio import DistanceMatrix


rule calculate_beta_div_effect_sizes:
    input:
        md_file = "results/filtered_metadata.tsv",
        dm_file = "results/beta_div/{is_phylo}/{beta_div_metric}/distance-matrix.tsv"
    output:
        "results/beta_div/{is_phylo}/{beta_div_metric}/effect_sizes.tsv"
    log:
        "logs/calculate_beta_div_effect_sizes.{is_phylo}.{beta_div_metric}.log"
    run:
        xebec_logger = get_logger(log[0])
        md = pd.read_table(input["md_file"], sep="\t", index_col=0)
        dm = DistanceMatrix.read(input["dm_file"])

        bdh = BetaDiversityHandler(dm, md)
        res = effect_size_by_category(bdh, md.columns).to_dataframe()
        xebec_logger.info(f"\n{res.head()}")
        res.to_csv(output[0], sep="\t", index=True)


rule calculate_beta_div_pairwise_effect_sizes:
    input:
        md_file = "results/filtered_metadata.tsv",
        dm_file = "results/beta_div/{is_phylo}/{beta_div_metric}/distance-matrix.tsv"
    output:
        "results/beta_div/{is_phylo}/{beta_div_metric}/pairwise_effect_sizes.tsv"
    log:
        "logs/calculate_beta_div_pairwise_effect_sizes.{is_phylo}.{beta_div_metric}.log"
    run:
        xebec_logger = get_logger(log[0])
        md = pd.read_table(input["md_file"], sep="\t", index_col=0)
        dm = DistanceMatrix.read(input["dm_file"])

        bdh = BetaDiversityHandler(dm, md)
        res = pairwise_effect_size_by_category(bdh, md.columns).to_dataframe()
        xebec_logger.info(f"\n{res.head()}")
        res.to_csv(output[0], sep="\t", index=True)


def concatenate_metric_dataframes(files):
    """Concatenate results from multiple metrics."""
    def get_metric_info(f):
        """Return metric type and metric name as tuple.

        results/beta_div/*/*/effect_sizes.tsv
        """
        path_parts = PurePath(f).parts
        return path_parts[2], path_parts[3]

    all_dfs = []
    all_keys = []
    for f in files:
        this_df = pd.read_table(f, sep="\t", index_col=0)
        this_keys = get_metric_info(f)
        all_dfs.append(this_df)
        all_keys.append(this_keys)

    total_df = pd.concat(all_dfs, keys=all_keys, names=["phylogenetic", "diversity_metric"])
    total_df = total_df.reset_index(level=("phylogenetic", "diversity_metric"))
    return total_df


beta_div_effect_sizes = [
    f"results/beta_div/{row['phylogenetic']}/{row['diversity_metric']}/effect_sizes.tsv"
    for i, row in beta_metrics.iterrows()
]
beta_div_pw_effect_sizes = [
    f"results/beta_div/{row['phylogenetic']}/{row['diversity_metric']}/pairwise_effect_sizes.tsv"
    for i, row in beta_metrics.iterrows()
]


# Can't use double brace syntax for Snakemake wildcards in expand because this notation is
# used for cookiecutter.
rule concatenate_beta_div_effect_sizes:
    input:
        beta_div_effect_sizes
    output:
        "results/beta_div/all_metrics_effect_sizes.tsv"
    log:
        "logs/concatenate_beta_div_effect_sizes.log"
    run:
        all_metrics_df = concatenate_metric_dataframes(input)
        all_metrics_df.to_csv(output[0], sep="\t", index=False)


rule concatenate_beta_div_pairwise_effect_sizes:
    input:
        beta_div_pw_effect_sizes
    output:
        "results/beta_div/all_metrics_pairwise_effect_sizes.tsv"
    log:
        "logs/concatenate_beta_div_pairwise_effect_sizes.log"
    run:
        all_metrics_df = concatenate_metric_dataframes(input)
        all_metrics_df.to_csv(output[0], sep="\t", index=False)
