import os

import biom
from skbio.diversity import beta_diversity
import unifrac

rule rpca:
    input:
        "{{cookiecutter.feature_table_file}}"
    output:
        "../results/beta_div/non_phylo/rpca/distance-matrix.tsv",
        "../results/beta_div/non_phylo/rpca/ordination.txt"
    shell:
        """
        gemelli rpca \
            --in-biom {input} \
            --output-dir ../results/beta_div/non_phylo/rpca \
            --n_components {config[n_components]} \
            --min-sample-count 0
        """

rule non_phylo_beta_div:
    input:
        "{{cookiecutter.feature_table_file}}"
    output:
        "../results/beta_div/non_phylo/{beta_div_metric}/distance-matrix.tsv"
    params:
        "../results/beta_div/non_phylo/{beta_div_metric}"
    run:
        os.makedirs(params[0], exist_ok=True)
        table = biom.load_table(input[0])

        dm = beta_diversity(
            metric=wildcards["beta_div_metric"],
            counts=table.matrix_data.todense().T,
            ids=table.ids("sample")
        )
        dm.write(output[0])


rule phylo_beta_div:
    input:
        tbl_file = "{{cookiecutter.feature_table_file}}",
        tree_file = "{{cookiecutter.phylogenetic_tree_file}}"
    output:
        "../results/beta_div/phylo/{beta_div_metric}/distance-matrix.tsv"
    params:
        "../results/beta_div/phylo/{beta_div_metric}"
    run:
        os.makedirs(params[0], exist_ok=True)

        if wildcards.beta_div_metric == "unweighted_unifrac":
            func = unifrac.unweighted_fp32
        elif wildcards.beta_div_metric == "weighted_unifrac":
            func = unifrac.weighted_normalized_fp32
        else:
            pass

        dm = func(input["tbl_file"], input["tree_file"])
        dm.write(output[0])
