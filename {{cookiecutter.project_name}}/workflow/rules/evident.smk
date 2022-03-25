from evident import AlphaDiversityHandler, BetaDiversityHandler
from evident.effect_size import effect_size_by_category, pairwise_effect_size_by_category
import numpy as np
import pandas as pd
from skbio import DistanceMatrix


rule filter_metadata:
    input:
        "{{cookiecutter.sample_metadata_file}}"
    output:
        "results/filtered_metadata.tsv"
    log:
        "logs/filter_metadata.log"
    run:
        metadata = pd.read_table(input[0], sep="\t", index_col=0)

        md = metadata.copy()
        all_cols = md.columns
        cols_to_drop = []
        for col in all_cols:
            # Drop columns that are not categorical
            if md[col].dtype != np.dtype("object"):
                cols_to_drop.append(col)
                continue

            # Drop columns that have an invalid number of levels
            uniq_groups = md[col].dropna().unique()
            valid_cat_count = 1 < len(uniq_groups) <= config["max_category_levels"]
            if not valid_cat_count:
                cols_to_drop.append(col)
                continue

            # Remove levels without enough samples
            level_count = md[col].value_counts()
            under_thresh = level_count[level_count < config["min_level_count"]]
            if not under_thresh.empty:
                levels_under_thresh = list(under_thresh.index)
                md[col] = md[col].replace({x: np.nan for x in levels_under_thresh})

        md = md.drop(columns=cols_to_drop)
        md.to_csv(output[0], sep="\t", index=True)


rule calculate_beta_div_effect_sizes:
    input:
        md_file = "results/filtered_metadata.tsv",
        dm_file = "results/beta_div/{is_phylo}/{beta_div_metric}/distance-matrix.tsv"
    output:
        "results/beta_div/{is_phylo}/{beta_div_metric}/effect_sizes.tsv"
    log:
        "logs/calculate_beta_div_effect_sizes.{is_phylo}.{beta_div_metric}.log"
    run:
        md = pd.read_table(input["md_file"], sep="\t", index_col=0)
        dm = DistanceMatrix.read(input["dm_file"])

        bdh = BetaDiversityHandler(dm, md)
        res = effect_size_by_category(bdh, md.columns).to_dataframe()
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
        md = pd.read_table(input["md_file"], sep="\t", index_col=0)
        dm = DistanceMatrix.read(input["dm_file"])

        bdh = BetaDiversityHandler(dm, md)
        res = pairwise_effect_size_by_category(bdh, md.columns).to_dataframe()
        res.to_csv(output[0], sep="\t", index=True)


def concatenate_metric_dataframes(files):
    """Concatenate results from multiple metrics."""
    def get_metric_info(f):
        """Return metric type and metric name as tuple."""
        path_parts = PurePath(f).parts
        return path_parts[3], path_parts[4]

    all_dfs = []
    all_keys = []
    for f in files:
        this_df = pd.read_table(f, sep="\t", index_col=0)
        this_keys = get_metric_info(f)
        all_dfs.append(this_df)
        all_keys.append(this_keys)

    total_df = pd.concat(all_dfs, keys=all_keys, names=["metric_type", "metric_name"])
    total_df = total_df.reset_index(level=("metric_type", "metric_name"))
    return total_df


beta_div_effect_sizes = [
    f"results/beta_div/{row['metric_type']}/{row['metric']}/effect_sizes.tsv"
    for i, row in beta_metrics.iterrows()
]
beta_div_pw_effect_sizes = [
    f"results/beta_div/{row['metric_type']}/{row['metric']}/pairwise_effect_sizes.tsv"
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
