import os

import biom
from skbio.diversity import beta_diversity

from xebec import get_logger


xebec_logger = get_logger(snakemake.log[0], snakemake.rule)
div_metric = snakemake.wildcards["beta_div_metric"]
xebec_logger.info(f"Diversity Metric: {div_metric}")
xebec_logger.info(f"Phylogenetic: False")  # Phylo handled by ssu
os.makedirs(snakemake.params["out_dir"], exist_ok=True)
table = biom.load_table(snakemake.input[0])

xebec_logger.info("Calculating beta diversity...")
dm = beta_diversity(
    metric=snakemake.wildcards["beta_div_metric"],
    counts=table.matrix_data.todense().T,
    ids=table.ids("sample")
)
xebec_logger.info("Finished calculating beta diversity!")
dm.write(snakemake.output[0])
xebec_logger.info(f"Saved to {snakemake.output[0]}")
