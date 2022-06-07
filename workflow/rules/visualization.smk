from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import output_file, save, figure
import seaborn as sns


def get_div_metrics(wildcards):
    if wildcards.diversity_type == "alpha_div":
        return alpha_metrics["diversity_metric"].tolist()
    else:
        return beta_metrics["diversity_metric"].tolist()


def get_non_phylo_div_metrics(wildcards):
    if wildcards.diversity_type == "alpha_div":
        return non_phylo_alpha_metrics
    else:
        return non_phylo_beta_metrics


def get_phylo_div_metrics(wildcards):
    if wildcards.diversity_type == "alpha_div":
        return phylo_alpha_metrics
    else:
        return phylo_beta_metrics


rule plot_effect_sizes:
    input:
        "results/{diversity_type}/all_metrics_effect_sizes.tsv"
    output:
        "results/{diversity_type}/effect_size_plot.html"
    params:
        all_div_metrics = get_div_metrics,
        non_phylo_metrics = get_non_phylo_div_metrics,
        phylo_metrics = get_phylo_div_metrics
    log:
        "logs/plot_{diversity_type}_effect_sizes.log"
    script:
        # The script path is always relative to the Snakefile
        #   containing the directive
        "../scripts/interactive_effect_sizes.py"


rule plot_pairwise_effect_sizes:
    input:
        "results/{diversity_type}/all_metrics_pairwise_effect_sizes.tsv"
    output:
        "results/{diversity_type}/pairwise_effect_size_plot.html"
    params:
        all_div_metrics = get_div_metrics,
        non_phylo_metrics = get_non_phylo_div_metrics,
        phylo_metrics = get_phylo_div_metrics
    log:
        "logs/plot_{diversity_type}_pairwise_effect_sizes.log"
    script:
        # The script path is always relative to the Snakefile
        #   containing the directive
        "../scripts/interactive_pw_effect_sizes.py"
