import os
import shutil

import biom
import numpy as np
import pandas as pd
from skbio import TreeNode

def main():
    curr_dir = os.path.dirname(__file__)
    tbl_file = os.path.join(curr_dir, "raw/127612_reference-hit.biom")
    md_file = os.path.join(curr_dir, "raw/11402_20171025-114722.txt")
    tree_file = os.path.join(curr_dir, "raw/127612_insertion_tree.relabelled.tre")

    tbl = biom.load_table(tbl_file)

    na_values = ["Missing: Not collected", "Not applicable"]
    cols_to_drop = [
        "parasites",
        "hb_iron_stats",
        "diff_btw_fecal_serum_collection",
        "collection_timestamp"
    ]
    md = pd.read_table(md_file, sep="\t", index_col=0, na_values=na_values)

    samps_in_common = set(tbl.ids()).intersection(md.index)
    min_depth = 1000
    depths = tbl.sum(axis="sample")
    samps_in_common = list(samps_in_common.intersection(
        tbl.ids()[np.where(depths >= min_depth)]
    ))

    md = md.loc[samps_in_common]
    tbl.filter(samps_in_common)

    min_prev = 10
    prev = tbl.pa(inplace=False).sum(axis="observation")
    feats_to_keep = tbl.ids("observation")[np.where(prev >= min_prev)]
    tbl.filter(feats_to_keep, "observation")

    for col in md.columns:
        col_type = md[col].dtype
        if col_type != np.dtype("object"):
            cols_to_drop.append(col)
            continue
        if not (1 < len(md[col].dropna().unique()) <= 5):
            cols_to_drop.append(col)
            continue

        level_count = md[col].value_counts()
        under_thresh = level_count[level_count < 3]
        if not under_thresh.empty:
            levels_under_thresh = list(under_thresh.index)
            md[col].replace(
                {x: np.nan for x in levels_under_thresh},
                inplace=True
            )

    md = md.drop(columns=cols_to_drop)
    print(f"Metadata shape: {md.shape}")
    print(f"Table shape: {tbl.shape}")

    with biom.util.biom_open("./processed/table.biom", "w") as f:
        tbl.to_hdf5(f, "filtered")

    md.to_csv("./processed/metadata.tsv", index=True, sep="\t")
    shutil.copyfile(tree_file, "./processed/tree.tre")

if __name__ == "__main__":
    main()
