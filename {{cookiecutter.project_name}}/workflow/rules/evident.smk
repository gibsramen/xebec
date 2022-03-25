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
    run:
        md = pd.read_table(input["md_file"], sep="\t", index_col=0)
        dm = DistanceMatrix.read(dm_file)

        bdh = BetaDiversityHandler(dm, md)
        res = pairwise_effect_size_by_category(bdh, md.columns).to_dataframe()
        res.to_csv(output[0], sep="\t", index=True)
