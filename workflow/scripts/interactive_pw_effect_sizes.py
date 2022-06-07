from pkg_resources import resource_filename

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Select
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import output_file, save, figure
import pandas as pd

import xebec.src._visualization as viz


diversity_metric_order = snakemake.params["all_div_metrics"]
non_phylo_metrics = snakemake.params["non_phylo_metrics"]
phylo_metrics = snakemake.params["phylo_metrics"]

palette = viz.get_scatter_palette(phylo_metrics, non_phylo_metrics)

df = pd.read_table(snakemake.input[0], sep="\t")
df["comparison"] = df["group_1"] + " vs. " + df["group_2"]
df["color"] = df["diversity_metric"].map(palette)
df["marker"] = df["phylogenetic"].map(
    {"phylo": "triangle", "non_phylo": "circle"}
)
all_cols = df["column"].unique().tolist()

chosen_col = Select(options=all_cols, value=all_cols[0], title="Column")
col_df = df.copy()[df["column"] == chosen_col.value]
col_df["diversity_metric"] = col_df["diversity_metric"].astype("category")
col_df["diversity_metric"] = (
    col_df["diversity_metric"]
    .cat
    .set_categories(diversity_metric_order)
)
comparisons = col_df["comparison"].unique()
div_metrics = col_df["diversity_metric"].unique()

order = list(
    col_df
    .groupby("comparison")
    .median()
    .sort_values(by="effect_size", ascending=False)
    .index
)

hover_points = viz.HOVER_POINTS
hover_boxes = viz.HOVER_BOXES

p = figure(
    tools=["pan", "reset", "box_zoom", hover_boxes, hover_points],
    y_range=order,
    sizing_mode="stretch_both"
)

big_source = ColumnDataSource(df)
col_source = ColumnDataSource(col_df)
box_source = viz.add_boxplots(p, col_df, "comparison")

p.scatter(
    source=col_source,
    y="comparison",
    x="effect_size",
    size=10,
    name="points",
    line_width=0.5,
    line_color="black",
    fill_color="color",
    legend_field="diversity_metric",
    marker="marker"
)

callback_file = resource_filename("xebec", "js/chosen_col_callback.js")
with open(callback_file, "r") as f:
    callback_code = f.read()

callback = CustomJS(
    args=dict(
        big_source=big_source,
        box_source=box_source,
        col_source=col_source,
        yr=p.y_range
    ),
    code=callback_code
)
chosen_col.js_on_change("value", callback)

p.title = "Pairwise Comparisons"
p.xaxis.axis_label = "Cohen's d"

viz.assign_fig_parameters(p)

controls = [chosen_col]
control_panel = column(controls, width=200, height=200)
layout = row(control_panel, p)
output_file(snakemake.output[0])
save(layout, title="xebec")
