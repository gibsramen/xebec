from evident import AlphaDiversityHandler, BetaDiversityHandler
from evident.effect_size import effect_size_by_category, pairwise_effect_size_by_category
import numpy as np
import pandas as pd
from skbio import DistanceMatrix

from helper import get_logger


md = pd.read_table(snakemake.input["md_file"], sep="\t", index_col=0)

if snakemake.params["div_type"] == "alpha":
    data = pd.read_table(snakemake.input["ad_file"], sep="\t", index_col=0).squeeze()
    dh = AlphaDiversityHandler(data, md)
elif snakemake.params["div_type"] == "beta":
    data = DistanceMatrix.read(snakemake.input["dm_file"])
    dh = BetaDiversityHandler(data, md)
else:
    pass

if snakemake.params["pairwise"] == "True":
    func = pairwise_effect_size_by_category
elif snakemake.params["pairwise"] == "False":
    func = effect_size_by_category
else:
    pass

xebec_logger = get_logger(snakemake.log[0])
res = func(dh, md.columns).to_dataframe()
xebec_logger.info(f"\n{res.head()}")
res.to_csv(snakemake.output[0], sep="\t", index=False)
