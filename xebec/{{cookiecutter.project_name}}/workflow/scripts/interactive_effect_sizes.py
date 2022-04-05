from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource
from bokeh.plotting import output_file, save, figure
import pandas as pd

import xebec.src._visualization as viz


def generate_interactive_effect_sizes(es_metric: str):
    """Generate boxplots for different diversity metrics"""
    _df = df.copy().query("metric == @es_metric")

    # Order columns by median effect size
    order = list(
        _df.groupby("column")
        .median()
        .sort_values(by="effect_size", ascending=False).index
    )

    # https://stackoverflow.com/a/27255567
    _df["column"] = _df["column"].astype("category")
    _df["column"] = _df["column"].cat.set_categories(order)
    _df = _df.sort_values(by="column")

    # Sort diversity metrics by non-phylo -> phylo
    _df["diversity_metric"] = _df["diversity_metric"].astype("category")
    _df["diversity_metric"] = (
        _df["diversity_metric"]
        .cat
        .set_categories(diversity_metric_order)
    )

    cols = _df["column"].unique()
    div_metrics = _df["diversity_metric"].unique()

    hover_points = viz.HOVER_POINTS
    hover_boxes = viz.HOVER_BOXES

    p = figure(
        tools=["pan", "reset", "box_zoom", hover_points, hover_boxes],
        y_range=order,
        width=800
    )
    p.background_fill_color = viz.BACKGROUND_FILL_COLOR

    viz.add_boxplots(p, _df, "column")

    for i, (col, col_df) in enumerate(_df.groupby("diversity_metric")):
        source = ColumnDataSource(col_df)
        phylogenetic = col_df["phylogenetic"].unique().item()
        if phylogenetic == "phylo":
            marker = "triangle"
        else:
            marker = "circle"
        scatter = p.scatter(
            source=source,
            x="effect_size",
            y="column",
            size=10,
            name="points",
            line_width=0.5,
            line_color="black",
            marker=marker,
            fill_color=palette[col],
            legend_label=col,
        )

    p.legend.location = "top_right"
    p.legend.click_policy = "hide"
    p.legend.title = "Click Entries to Toggle"
    p.legend.title_text_font_style = "bold"

    if es_metric == "cohens_d":
        p.xaxis.axis_label = "Cohen's d"
        p.title = "Binary Categories"
    else:
        p.title = "Multi-Class Categories"
        p.xaxis.axis_label = "Cohen's f"

    viz.assign_fig_parameters(p)

    return p

diversity_metric_order = snakemake.params["all_div_metrics"]
non_phylo_metrics = snakemake.params["non_phylo_metrics"]
phylo_metrics = snakemake.params["phylo_metrics"]

palette = viz.get_scatter_palette(phylo_metrics, non_phylo_metrics)

df = pd.read_table(snakemake.input[0], sep="\t")
output_file(snakemake.output[0])
p1 = generate_interactive_effect_sizes("cohens_d")
p2 = generate_interactive_effect_sizes("cohens_f")
layout = gridplot([[p1, p2]], sizing_mode="scale_width",
                  toolbar_location="right")
save(layout, title="xebec")
