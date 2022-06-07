from pathlib import PurePath

from evident import UnivariateDataHandler, MultivariateDataHandler
from evident.effect_size import effect_size_by_category, pairwise_effect_size_by_category
import numpy as np
import pandas as pd
from skbio import DistanceMatrix


rule calculate_alpha_div_effect_sizes:
    input:
        md_file = "results/filtered_metadata.tsv",
        ad_file = "results/alpha_div/{is_phylo}/{alpha_div_metric}/vector.tsv"
    output:
        "results/alpha_div/{is_phylo}/{alpha_div_metric}/effect_sizes.tsv"
    params:
        div_type = "alpha",
        pairwise = "False"
    log:
        "logs/calculate_alpha_div_effect_sizes.{is_phylo}.{alpha_div_metric}.log"
    script:
        "../scripts/run_evident.py"


rule calculate_beta_div_effect_sizes:
    input:
        md_file = "results/filtered_metadata.tsv",
        dm_file = "results/beta_div/{is_phylo}/{beta_div_metric}/distance-matrix.tsv"
    output:
        "results/beta_div/{is_phylo}/{beta_div_metric}/effect_sizes.tsv"
    params:
        div_type = "beta",
        pairwise = "False"
    log:
        "logs/calculate_beta_div_effect_sizes.{is_phylo}.{beta_div_metric}.log"
    script:
        "../scripts/run_evident.py"


rule calculate_alpha_div_pairwise_effect_sizes:
    input:
        md_file = "results/filtered_metadata.tsv",
        ad_file = "results/alpha_div/{is_phylo}/{alpha_div_metric}/vector.tsv"
    output:
        "results/alpha_div/{is_phylo}/{alpha_div_metric}/pairwise_effect_sizes.tsv"
    params:
        div_type = "alpha",
        pairwise = "True"
    log:
        "logs/calculate_alpha_div_pairwise_effect_sizes.{is_phylo}.{alpha_div_metric}.log"
    script:
        "../scripts/run_evident.py"


rule calculate_beta_div_pairwise_effect_sizes:
    input:
        md_file = "results/filtered_metadata.tsv",
        dm_file = "results/beta_div/{is_phylo}/{beta_div_metric}/distance-matrix.tsv"
    output:
        "results/beta_div/{is_phylo}/{beta_div_metric}/pairwise_effect_sizes.tsv"
    params:
        div_type = "beta",
        pairwise = "True"
    log:
        "logs/calculate_beta_div_pairwise_effect_sizes.{is_phylo}.{beta_div_metric}.log"
    script:
        "../scripts/run_evident.py"


rule concatenate_alpha_div_effect_sizes:
    input:
        alpha_div_effect_sizes
    output:
        "results/alpha_div/all_metrics_effect_sizes.tsv"
    log:
        "logs/concatenate_alpha_div_effect_sizes.log"
    script:
        "../scripts/concatenate.py"


rule concatenate_beta_div_effect_sizes:
    input:
        beta_div_effect_sizes
    output:
        "results/beta_div/all_metrics_effect_sizes.tsv"
    log:
        "logs/concatenate_beta_div_effect_sizes.log"
    script:
        "../scripts/concatenate.py"


rule concatenate_alpha_div_pairwise_effect_sizes:
    input:
        alpha_div_pw_effect_sizes
    output:
        "results/alpha_div/all_metrics_pairwise_effect_sizes.tsv"
    log:
        "logs/concatenate_alpha_div_pairwise_effect_sizes.log"
    script:
        "../scripts/concatenate.py"


rule concatenate_beta_div_pairwise_effect_sizes:
    input:
        beta_div_pw_effect_sizes
    output:
        "results/beta_div/all_metrics_pairwise_effect_sizes.tsv"
    log:
        "logs/concatenate_beta_div_pairwise_effect_sizes.log"
    script:
        "../scripts/concatenate.py"
