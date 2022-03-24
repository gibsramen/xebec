# xebec

Snakemake pipeline for microbiome diversity effect size benchmarking

## Installation

To use xebec, you will need several dependencies:

* snakemake
* cookiecutter
* unifrac
* scikit-bio
* pandas
* evident
* gemelli

We recommend using `conda`/`mamba` to install these packages when possible.
Note that at time of writing, evident and gemelli are only available through PyPi.

From the command line, run the following command:

```
cookiecutter https://github.com/gibsramen/xebec
```

You should enter a prompt where you can input the required values to setup xebec.

* `project_name`: Name of the directory to create with the Snakemake pipeline files (defaults to `diversity-benchmark`)
* `feature_table_file`: *absolute* path to the feature table to be used in BIOM format
* `sample_metadata_file`: *absolute* path to the sample metadata file to be used in TSV format
* `phylogenetic_tree_file`: *absolute* path to the phylogenetic tree file to be used in Newick format

This will create the directory structure needed to run xebec under the project name you specified.

## Usage

Navigate inside the `<project_name>/workflow' directory.
You should see a file called `Snakefile` - this file is used to run xebec.
To perform all beta-diversity calculations, run the following command:

```
snakemake beta_diversity --cores 1
```

You should see the Snakemake pipeline start running the jobs.
The resulting distance matrices will be stored inside `<project_name>/results/beta-div`.
