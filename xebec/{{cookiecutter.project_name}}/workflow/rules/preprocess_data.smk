rule filter_metadata:
    input:
        md_file = "{{cookiecutter.sample_metadata_file}}",
        tbl_file = "{{cookiecutter.feature_table_file}}"
    output:
        "results/filtered_metadata.tsv"
    log:
        "logs/filter_metadata.log"
    script:
        "../scripts/filter_metadata.py"


rule rarefy:
    input:
        "{{cookiecutter.feature_table_file}}"
    output:
        "results/rarefied_table.biom"
    log:
        "logs/rarefy.log"
    script:
        "../scripts/rarefy.py"


rule create_shuffled_metadata:
    input:
        "results/filtered_metadata.tsv"
    output:
        [f"results/shuffled/shuffled_metadata.{x+1}.tsv"
        for x in range(config["shuffle_iterations"])]
    params:
        "results/shuffled/"
    script:
        "../scripts/shuffle_metadata.py"
