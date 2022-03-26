import biom
import numpy as np
import pandas as pd


rule filter_metadata:
    input:
        "{{cookiecutter.sample_metadata_file}}"
    output:
        "results/filtered_metadata.tsv"
    log:
        "logs/filter_metadata.log"
    run:
        xebec_logger = get_logger(log[0])
        metadata = pd.read_table(input[0], sep="\t", index_col=0)
        xebec_logger.info(f"Original metadata shape: {metadata.shape}")

        md = metadata.copy()
        all_cols = md.columns
        cols_to_drop = []
        for col in all_cols:
            # Drop columns that are not categorical
            if md[col].dtype != np.dtype("object"):
                cols_to_drop.append(col)
                xebec_logger.info(
                    f"Dropping {col} because it is not of dtype 'object'."
                )
                continue

            # Drop columns that have an invalid number of levels
            uniq_groups = md[col].dropna().unique()
            valid_cat_count = 1 < len(uniq_groups) <= config["max_category_levels"]
            if not valid_cat_count:
                cols_to_drop.append(col)
                xebec_logger.info(
                    f"Dropping {col} because it does not have a valid number "
                    f"of levels ({len(uniq_groups)})."
                )
                continue

            # Remove levels without enough samples
            level_count = md[col].value_counts()
            under_thresh = level_count[level_count < config["min_level_count"]]
            if not under_thresh.empty:
                levels_under_thresh = list(under_thresh.index)
                md[col] = md[col].replace({x: np.nan for x in levels_under_thresh})
                xebec_logger.info(
                    f"Dropping levels {levels_under_thresh} in category {col} "
                    "because they fall under the minimum level count threshold."
                )

        md = md.drop(columns=cols_to_drop)
        xebec_logger.info(f"Filtered metadata shape: {md.shape}")
        xebec_logger.info(f"\n{md.head()}")
        for col in md.columns:
            xebec_logger.info(f"\n{md[col].value_counts()}")
        md.to_csv(output[0], sep="\t", index=True)


rule rarefy:
    input:
        "{{cookiecutter.feature_table_file}}"
    output:
        "results/rarefied_table.biom"
    log:
        "logs/rarefy.log"
    run:
        xebec_logger = get_logger(log[0])
        table = biom.load_table(input[0])
        xebec_logger.info(f"Original shape: {table.shape}")
        depths = table.sum(axis="sample")
        rare_depth = round(np.quantile(depths, config["rarefaction_depth_percentile"]/100))
        xebec_logger.info(f"Rarefaction depth: {rare_depth}")
        table_rare = table.subsample(rare_depth)
        xebec_logger.info(f"Rarefied shape: {table_rare.shape}")

        with biom.util.biom_open(output[0], "w") as f:
            table_rare.to_hdf5(f, "rarefy")
