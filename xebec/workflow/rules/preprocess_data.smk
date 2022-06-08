rule filter_metadata:
    input:
        md_file=config["sample_metadata_file"],
        tbl_file=config["feature_table_file"]
    output:
        "results/filtered_metadata.tsv"
    log:
        "logs/filter_metadata.log"
    script:
        "../scripts/filter_metadata.py"


rule rarefy:
    input:
        config["feature_table_file"]
    output:
        "results/rarefied_table.biom"
    log:
        "logs/rarefy.log"
    script:
        "../scripts/rarefy.py"
