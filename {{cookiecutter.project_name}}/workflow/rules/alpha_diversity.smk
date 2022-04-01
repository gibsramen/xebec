import os

from bp import parse_newick, to_skbio_treenode
import pandas as pd
from skbio.diversity import alpha_diversity
from skbio import TreeNode


rule non_phylo_alpha_div:
    input:
        tbl_file = "results/rarefied_table.biom"
    output:
        "results/alpha_div/non_phylo/{alpha_div_metric}/vector.tsv"
    log:
        "logs/{alpha_div_metric}.log"
    params:
        phylo = "False",
        out_dir = "results/alpha_div/non_phylo/{alpha_div_metric}"
    script:
        "../scripts/alpha_diversity.py"


rule phylo_alpha_div:
    input:
        tbl_file = "results/rarefied_table.biom",
        tree_file = "{{cookiecutter.phylogenetic_tree_file}}"
    output:
        "results/alpha_div/phylo/{alpha_div_metric}/vector.tsv"
    log:
        "logs/{alpha_div_metric}.log"
    params:
        phylo = "True",
        out_dir = "results/alpha_div/non_phylo/{alpha_div_metric}"
    script:
        "../scripts/alpha_diversity.py"
