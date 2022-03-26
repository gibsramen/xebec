import os

from bp import parse_newick, to_skbio_treenode
import pandas as pd
from skbio.diversity import alpha_diversity
from skbio import TreeNode


rule non_phylo_alpha_div:
    input:
        "results/rarefied_table.biom",
    output:
        "results/alpha_div/non_phylo/{alpha_div_metric}/vector.tsv"
    log:
        "logs/{alpha_div_metric}.log"
    params:
        "results/alpha_div/non_phylo/{alpha_div_metric}"
    run:
        os.makedirs(params[0], exist_ok=True)
        table = biom.load_table(input[0])

        alpha_div = alpha_diversity(
            metric=wildcards.alpha_div_metric,
            counts=table.matrix_data.todense().T,
            ids=table.ids("sample")
        )
        alpha_div.to_csv(output[0], sep="\t", index=True)


rule phylo_alpha_div:
    input:
        tbl_file = "results/rarefied_table.biom",
        tree_file = "{{cookiecutter.phylogenetic_tree_file}}"
    output:
        "results/alpha_div/phylo/{alpha_div_metric}/vector.tsv"
    log:
        "logs/{alpha_div_metric}.log"
    params:
        "results/alpha_div/phylo/{alpha_div_metric}"
    run:
        os.makedirs(params[0], exist_ok=True)
        table = biom.load_table(input["tbl_file"])

        with open(input["tree_file"]) as f:
            tree = parse_newick(f.readline())
        tree = to_skbio_treenode(tree)

        alpha_div = alpha_diversity(
            metric=wildcards.alpha_div_metric,
            counts=table.matrix_data.todense().T,
            ids=table.ids("sample"),
            otu_ids=table.ids("observation"),
            tree=tree
        )
        alpha_div.to_csv(output[0], sep="\t", index=True)
