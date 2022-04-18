import os


def test_bake_project(cookies, data_paths):
    result = cookies.bake(extra_context={
        "project_name": "example-benchmark",
        "feature_table_file": data_paths.table_file,
        "sample_metadata_file": data_paths.metadata_file,
        "phylogenetic_tree_file": data_paths.tree_file
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
    assert set(workflow_files) == {"rules", "Snakefile", "scripts", "report"}

    rules_dir = os.path.join(workflow_dir, "rules")
    rules_files = os.listdir(rules_dir)
    exp_rules = {
        "alpha_diversity.smk",
        "beta_diversity.smk",
        "evident.smk",
        "visualization.smk",
        "preprocess_data.smk",
    }
    assert set(rules_files) == exp_rules

    scripts_dir = os.path.join(workflow_dir, "scripts")
    scripts_files = os.listdir(scripts_dir)
    assert set(scripts_files) == {
        "interactive_effect_sizes.py",
        "interactive_pw_effect_sizes.py",
        "alpha_diversity.py",
        "beta_diversity.py",
        "concatenate.py",
        "filter_metadata.py",
        "rarefy.py",
        "run_evident.py"
    }

    report_dir = os.path.join(workflow_dir, "report")
    report_files = os.listdir(report_dir)
    assert set(report_files) == {
        "effect_size_plot.rst",
        "pw_effect_size_plot.rst",
        "workflow.rst",
    }
