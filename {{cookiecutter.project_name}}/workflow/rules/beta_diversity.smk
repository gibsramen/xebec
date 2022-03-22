import os

import biom
from skbio.diversity import beta_diversity
import gemelli

rule rpca:
    input:
        "{{cookiecutter.feature_table_file}}"
    output:
        "../results/beta_div/rpca/distance-matrix.tsv",
        "../results/beta_div/rpca/ordination.txt"
    shell:
        """
        gemelli rpca \
            --in-biom {input} \
            --output-dir ../results/beta_div/rpca \
            --n_components {config[n_components]} \
            --min-sample-count 0
        """

rule non_phylo_beta_div:
    input:
        "{{cookiecutter.feature_table_file}}"
    output:
        "../results/beta_div/{beta_div_metric}/distance-matrix.tsv"
    params:
        "../results/beta_div/{beta_div_metric}"
    run:
        os.mkdir(params[0])
        table = biom.load_table(input[0])

        dm = beta_diversity(
            metric=wildcards["beta_div_metric"],
            counts=table.matrix_data.todense().T,
            ids=table.ids("sample")
        )
        dm.write(output[0])
