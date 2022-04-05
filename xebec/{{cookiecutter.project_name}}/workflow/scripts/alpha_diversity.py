import os

import biom
from bp import parse_newick, to_skbio_treenode
from skbio.diversity import alpha_diversity
from skbio import TreeNode

from xebec import get_logger


xebec_logger = get_logger(snakemake.log[0], snakemake.rule)
div_metric = snakemake.wildcards["alpha_div_metric"]
xebec_logger.info(f"Diversity Metric: {div_metric}")
xebec_logger.info(f"Phylogenetic: {snakemake.params['phylo']}")

os.makedirs(snakemake.params["out_dir"], exist_ok=True)
table = biom.load_table(snakemake.input["tbl_file"])

args = {
    "metric": snakemake.wildcards.alpha_div_metric,
    "counts": table.matrix_data.todense().T,
    "ids": table.ids("sample")
}
if snakemake.params["phylo"] == "True":
    with open(snakemake.input["tree_file"]) as f:
        tree = parse_newick(f.readline())
    args["tree"] = to_skbio_treenode(tree)
    args["otu_ids"] = table.ids("observation")

xebec_logger.info("Calculating alpha diversity...")
alpha_div = alpha_diversity(**args)
xebec_logger.info("Finished calculating alpha diversity!")
alpha_div.to_csv(snakemake.output[0], sep="\t", index=True)
xebec_logger.info(f"Saved to {snakemake.output[0]}")
