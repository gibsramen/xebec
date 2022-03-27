import os

curr_path = os.path.dirname(__file__)
table_file = os.path.abspath(os.path.join(curr_path, "data/table.biom"))
metadata_file = os.path.abspath(os.path.join(curr_path, "data/metadata.tsv"))
tree_file = os.path.abspath(os.path.join(curr_path, "data/tree.tre"))


def test_bake_project(cookies):
    result = cookies.bake(extra_context={
        "project_name": "example-benchmark",
        "feature_table_file": table_file,
        "sample_metadata_file": metadata_file,
        "phylogenetic_tree_file": tree_file
    })

    assert result.exit_code == 0
    assert result.exception is None

    assert result.project_path.name == "example-benchmark"
    assert result.project_path.is_dir()

    files = os.listdir(result.project_path)
    assert set(files) == {"workflow", "config"}

    config_dir = os.path.join(result.project_path, "config")
    config_files = os.listdir(config_dir)
    assert set(config_files) == {
        "alpha_div_metrics.tsv",
        "beta_div_metrics.tsv",
        "config.yaml"
    }

    workflow_dir = os.path.join(result.project_path, "workflow")
    workflow_files = os.listdir(workflow_dir)
    assert set(workflow_files) == {"rules", "Snakefile", "scripts"}

    rules_dir = os.path.join(workflow_dir, "rules")
    rules_files = os.listdir(rules_dir)
    exp_rules = {"alpha_diversity.smk", "beta_diversity.smk",
                 "evident.smk", "visualization.smk",
                 "preprocess_data.smk", "common.smk"}
    assert set(rules_files) == exp_rules

    scripts_dir = os.path.join(workflow_dir, "scripts")
    scripts_files = os.listdir(scripts_dir)
    assert set(scripts_files) == {"interactive_effect_sizes.py"}
