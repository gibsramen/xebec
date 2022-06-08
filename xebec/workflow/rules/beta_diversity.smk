import re


## NON-PHYLOGENETIC METRICS ##
rule rpca:
    input:
        config["feature_table_file"]
    output:
        "results/beta_div/non_phylo/rpca/distance-matrix.tsv",
        "results/beta_div/non_phylo/rpca/ordination.txt"
    log:
        "logs/rpca.log"
    shell:
        """
        gemelli rpca \
            --in-biom {input} \
            --output-dir results/beta_div/non_phylo/rpca \
            --n-components {config[n_components]} \
            --min-sample-count 0 > {log} 2>&1
        """


rule non_phylo_beta_div:
    input:
        "results/rarefied_table.biom"
    output:
        "results/beta_div/non_phylo/{beta_div_metric}/distance-matrix.tsv"
    log:
        "logs/{beta_div_metric}.log"
    params:
        out_dir = "results/beta_div/non_phylo/{beta_div_metric}"
    script:
        "../scripts/beta_diversity.py"


## PHYLOGENETIC METRICS ##
rule phylo_rpca:
    input:
        tbl_file=config["feature_table_file"],
        tree_file=config["phylogenetic_tree_file"]
    output:
        "results/beta_div/phylo/phylo_rpca/distance-matrix.tsv",
        "results/beta_div/phylo/phylo_rpca/ordination.txt"
    log:
        "logs/phylo_rpca.log"
    shell:
        """
        gemelli phylogenetic-rpca \
            --in-biom {input.tbl_file} \
            --in-phylogeny {input.tree_file} \
            --output-dir results/beta_div/phylo/phylo_rpca \
            --n-components {config[n_components]} \
            --min-sample-count 0 > {log} 2>&1
        """

unifrac_regex = re.compile("(.*)_unifrac")

rule phylo_beta_div:
    input:
        tbl_file="results/rarefied_table.biom",
        tree_file=config["phylogenetic_tree_file"]
    output:
        "results/beta_div/phylo/{beta_div_metric}/distance-matrix.tsv"
    log:
        "logs/{beta_div_metric}.log"
    params:
        out_dir = "results/beta_div/phylo/{beta_div_metric}",
        method = lambda wildcards: unifrac_regex.search(wildcards.beta_div_metric).groups()[0]
    shell:
        """
        mkdir -p {params.out_dir}
        ssu -i {input.tbl_file} -t {input.tree_file} -m {params.method} -o {output} > {log} 2>&1
        """
