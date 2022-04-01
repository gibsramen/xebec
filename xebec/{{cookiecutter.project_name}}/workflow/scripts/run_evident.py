from evident import AlphaDiversityHandler, BetaDiversityHandler
from evident.effect_size import effect_size_by_category, pairwise_effect_size_by_category
import numpy as np
import pandas as pd
from skbio import DistanceMatrix

from xebec import get_logger


xebec_logger = get_logger(snakemake.log[0], snakemake.rule)
xebec_logger.info(f"Diversity Type: {snakemake.params['div_type']}")
xebec_logger.info(f"Pairwise: {snakemake.params['pairwise']}")

md = pd.read_table(snakemake.input["md_file"], sep="\t", index_col=0)

if snakemake.params["div_type"] == "alpha":
    data = pd.read_table(snakemake.input["ad_file"], sep="\t", index_col=0).squeeze()
    dh = AlphaDiversityHandler(data, md)
    div_metric = snakemake.wildcards["alpha_div_metric"]
elif snakemake.params["div_type"] == "beta":
    data = DistanceMatrix.read(snakemake.input["dm_file"])
    dh = BetaDiversityHandler(data, md)
    div_metric = snakemake.wildcards["beta_div_metric"]
else:
    pass

if snakemake.params["pairwise"] == "True":
    func = pairwise_effect_size_by_category
elif snakemake.params["pairwise"] == "False":
    func = effect_size_by_category
else:
    pass

xebec_logger.info(f"Diversity Metric: {div_metric}")
xebec_logger.info("Calculating effect sizes...")
res = func(dh, md.columns).to_dataframe()
xebec_logger.info("Finished calculating effect sizes!")
xebec_logger.info(f"\n{res.head()}")
res.to_csv(snakemake.output[0], sep="\t", index=False)
xebec_logger.info(f"Saved to {snakemake.output[0]}!")
