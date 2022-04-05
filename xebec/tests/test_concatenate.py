import os
from pathlib import Path

import numpy as np
import pandas as pd

from xebec.src import _concatenate as conc


def test_concat_metric_dfs(tmp_path, monkeypatch):
    rng = np.random.default_rng()
    def mock_get_metric_info(div_metric):
        return rng.choice(["phylo", "non_phylo"]), Path(div_metric).stem

    monkeypatch.setattr(
        "xebec.src._concatenate.get_metric_info",
        mock_get_metric_info
    )

    n = 10

    div_metric_types = ["charmander", "cyndaquil", "torchic"]
    es_metrics = ["cohens_d", "cohens_f"]
    columns = [f"C{i+1}" for i in range(n)]
    es_metric_choices = rng.choice(es_metrics, size=n)

    files = []
    for div_metric in div_metric_types:
        div_metric_dict = {
            "effect_size": rng.gamma(0.5, size=n),
            "metric": es_metric_choices,
            "column": columns
        }
        df = pd.DataFrame.from_dict(div_metric_dict)
        df = df.sort_values(by=["metric", "effect_size"],
                            ascending=[True, False])
        out_file = os.path.join(tmp_path, f"{div_metric}.tsv")
        df.to_csv(out_file, sep="\t", index=False)
        files.append(out_file)
    total_df = conc.concatenate_metric_dataframes(files)

    assert total_df.shape == (len(div_metric_types)*n, 5)
    assert set(total_df.columns) == {
        "phylogenetic",
        "diversity_metric",
        "effect_size",
        "metric",
        "column"
    }
    assert set(total_df["diversity_metric"].unique()) == set(div_metric_types)
    assert (
        set(total_df["phylogenetic"].unique()).issubset(["non_phylo", "phylo"])
    )
    assert set(total_df["metric"].unique()).issubset(es_metrics)
    assert set(total_df["column"]) == set(columns)


def test_get_metric_info():
    x = "results/beta_div/phylo/mimikyu/effect_sizes.tsv"
    y = "results/alpha_div/non_phylo/ampharos/effect_sizes.tsv"

    assert conc.get_metric_info(x) == ("phylo", "mimikyu")
    assert conc.get_metric_info(y) == ("non_phylo", "ampharos")
