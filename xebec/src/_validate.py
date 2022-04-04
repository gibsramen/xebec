import os

import biom
import pandas as pd
from bp import parse_newick, to_skbio_treenode


def validate_table(fpath: os.PathLike):
    biom.load_table(fpath)


def validate_metadata(fpath: os.PathLike):
    pd.read_table(fpath, sep="\t", index_col=0)


def validate_tree(fpath: os.PathLike):
    with open(fpath) as f:
        tree = parse_newick(f.readline())
    to_skbio_treenode(tree)