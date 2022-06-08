import os
from pkg_resources import resource_filename
import shutil

import click
from jinja2 import Template

from xebec import __version__
from xebec.src._validate import (validate_table, validate_metadata,
                                 validate_tree)

help = "Create workflow to benchmark alpha and beta diversity metrics."


@click.command(name="xebec", short_help=help)
@click.version_option(__version__)
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
@click.option("--n-pcoa-components", default=3, show_default=True,
               type=int, help="Number of PCoA components to compute.")
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
    n_pcoa_components,
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
        "n_pcoa_components": n_pcoa_components
    }

    os.makedirs(output)

    # Create Snakfile
    wkflow_dir = os.path.join(output, "workflow")
    os.makedirs(wkflow_dir)
    orig_snakefile_path = resource_filename("xebec", "workflow/Snakefile")
    snkfile_template = resource_filename("xebec", "template/Snakefile.jinja2")

    with open(snkfile_template, "r") as f:
        snakefile_text = (
            Template(f.read())
            .render({"original_snakefile_path": orig_snakefile_path})
        )

    new_snakefile_path = os.path.join(wkflow_dir, "Snakefile")
    with open(new_snakefile_path, "w") as f:
        f.write(snakefile_text)

    # Create config file
    cfg_dir = os.path.join(output, "config")
    os.makedirs(cfg_dir)
    cfg_template = resource_filename("xebec", "template/config.yaml.jinja2")

    with open(cfg_template, "r") as f:
        cfg_text = Template(f.read()) .render(args)

    cfg_file_path = os.path.join(cfg_dir, "config.yaml")
    with open(cfg_file_path, "w") as f:
        f.write(cfg_text)

    # Copy diversity metric files
    orig_a_div_file = resource_filename(
        "xebec", "config/alpha_div_metrics.tsv"
    )
    orig_b_div_file = resource_filename(
        "xebec", "config/beta_div_metrics.tsv"
    )

    new_a_div_file = os.path.join(cfg_dir, "alpha_div_metrics.tsv")
    new_b_div_file = os.path.join(cfg_dir, "beta_div_metrics.tsv")

    shutil.copy(orig_a_div_file, new_a_div_file)
    shutil.copy(orig_b_div_file, new_b_div_file)


if __name__ == "__main__":
    xebec()
