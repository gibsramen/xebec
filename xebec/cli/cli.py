import os
from pathlib import PurePath
import subprocess

import click
from cookiecutter.main import cookiecutter

from xebec import COOKIE_DIR
from xebec.src._validate import (validate_table, validate_metadata,
                                 validate_tree)


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
@click.option("--execute/--no-execute", default=False, type=bool,
               help="Whether to automatically execute the workflow.",
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
    execute
):
    feature_table = os.path.abspath(feature_table)
    metadata = os.path.abspath(metadata)
    tree = os.path.abspath(tree)

    if validate_input:
        validate_table(feature_table)
        validate_metadata(metadata)
        validate_tree(tree)

    output = PurePath(output)
    project_dir = output.parent
    project_name = os.path.basename(output)
    os.chdir(project_dir)

    args={
        "project_name": project_name,
        "feature_table_file": feature_table,
        "sample_metadata_file": metadata,
        "phylogenetic_tree_file": tree,
        "max_category_levels": max_category_levels,
        "min_level_count": min_level_count,
        "rarefaction_depth_percentile": rarefy_percentile,
    }

    cookiecutter(COOKIE_DIR, no_input=True, extra_context=args)

    if execute:
        os.chdir(project_name)
        subprocess.run(["snakemake", "--cores", "1"], check=True)


if __name__ == "__main__":
    xebec()
