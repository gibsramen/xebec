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
    script:
        "../scripts/filter_metadata.py"


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
