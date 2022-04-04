from collections import namedtuple
import os

import pytest


@pytest.fixture
def data_paths():
    curr_path = os.path.dirname(__file__)
    table_file = os.path.abspath(os.path.join(curr_path, "data/table.biom"))
    metadata_file = os.path.abspath(os.path.join(curr_path, "data/metadata.tsv"))
    tree_file = os.path.abspath(os.path.join(curr_path, "data/tree.tre"))

    xebec_paths = namedtuple(
        "xebec_paths",
        ["table_file", "metadata_file", "tree_file"]
    )
    return xebec_paths(table_file, metadata_file, tree_file)
