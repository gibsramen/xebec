import pandas as pd


alpha_metrics_path = "config/alpha_div_metrics.tsv"
alpha_metrics = pd.read_table(alpha_metrics_path, sep="\t")

beta_metrics_path = "config/beta_div_metrics.tsv"
beta_metrics = pd.read_table(beta_metrics_path, sep="\t")

non_phylo_alpha_metrics = (
    alpha_metrics
    .query("phylogenetic == 'non_phylo'")
    ["diversity_metric"]
)
phylo_alpha_metrics = (
    alpha_metrics
    .query("phylogenetic == 'phylo'")
    ["diversity_metric"]
)

non_phylo_beta_metrics = (
    beta_metrics
    .query("phylogenetic == 'non_phylo'")
    ["diversity_metric"]
)
phylo_beta_metrics = (
    beta_metrics
    .query("phylogenetic == 'phylo'")
    ["diversity_metric"]
)

alpha_div_effect_sizes = [
    f"results/alpha_div/{row['phylogenetic']}/{row['diversity_metric']}/effect_sizes.tsv"
    for i, row in alpha_metrics.iterrows()
]
alpha_div_pw_effect_sizes = [
    f"results/alpha_div/{row['phylogenetic']}/{row['diversity_metric']}/pairwise_effect_sizes.tsv"
    for i, row in alpha_metrics.iterrows()
]

beta_div_effect_sizes = [
    f"results/beta_div/{row['phylogenetic']}/{row['diversity_metric']}/effect_sizes.tsv"
    for i, row in beta_metrics.iterrows()
]
beta_div_pw_effect_sizes = [
    f"results/beta_div/{row['phylogenetic']}/{row['diversity_metric']}/pairwise_effect_sizes.tsv"
    for i, row in beta_metrics.iterrows()
]
