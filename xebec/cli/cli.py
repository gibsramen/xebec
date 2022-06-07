import logging
import os
from pathlib import Path
import subprocess

import click
from jinja2 import Template
from snakedeploy.deploy import deploy
from snakedeploy.logger import logger

from xebec import __version__
from xebec._validate import (validate_table, validate_metadata,
                             validate_tree)

CONFIG_TEMPLATE = """{
    "feature_table_file": {{ feature_table_file }},
    "sample_metadata_file": {{ sample_metadata_file }},
    "phylogenetic_tree_file" : {{ phylogenetic_tree_file }},
    "max_category_levels": {{ max_category_levels }},
    "min_level_count": {{ min_level_count }},
    "rarefaction_depth_percentile": {{ rarefaction_depth_percentile }},
}
"""


@click.command(name="xebec")
@click.option("--feature-table", "-ft", required=True, type=click.Path(),
               help="Feature table in BIOM format.")
@click.option("--metadata", "-m", required=True, type=click.Path(),
               help="Sample metadata in TSV format.")
@click.option("--tree", "-t", required=True, type=click.Path(),
               help="Phylogenetic tree in Newick format.")
@click.option("--output", "-o", required=True, type=click.Path(),
               help="Output workflow directory.")
@click.option("--max-category-levels", default=5, show_default=True,
               type=int, help="Max number of levels in a category.")
@click.option("--min-level-count", default=3, show_default=True,
               type=int, help="Min number of samples per level per category.")
@click.option("--rarefy-percentile", default=10, show_default=True,
               type=float, help="Percentile of sample depths at which to rarefy.")
@click.option("--validate-input/--no-validate-input", default=True,
              help="Whether to validate input before creating workflow.",
              show_default=True)
def xebec(
    feature_table,
    metadata,
    tree,
    output,
    max_category_levels,
    min_level_count,
    rarefy_percentile,
    validate_input,
):
    feature_table = os.path.abspath(feature_table)
    metadata = os.path.abspath(metadata)
    tree = os.path.abspath(tree)

    if validate_input:
        validate_table(feature_table)
        validate_metadata(metadata)
        validate_tree(tree)

    args={
        "feature_table_file": feature_table,
        "sample_metadata_file": metadata,
        "phylogenetic_tree_file": tree,
        "max_category_levels": max_category_levels,
        "min_level_count": min_level_count,
        "rarefaction_depth_percentile": rarefy_percentile,
    }

    cfg_text = Template(CONFIG_TEMPLATE).render(args)

    xebec_repo_url = "https://github.com/gibsramen/xebec"
    xebec_repo_tag = f"v{__version__}"

    # Suppress warning that the config file doesn't exist since we copy it over
    # afterwards with Jinja
    logger.quiet = True
    logger.logger.setLevel(logging.ERROR)

    deploy(
        xebec_repo_url, dest_path=Path(output), name="diversity-benchmark",
        tag=xebec_repo_tag, branch=None
    )

    cfg_dir = os.path.join(output, "config")
    os.makedirs(cfg_dir, exist_ok=True)
    cfg_file = os.path.join(cfg_dir, "config.yaml")
    with open(cfg_file, "w") as f:
        f.write(cfg_text)


if __name__ == "__main__":
    xebec()
