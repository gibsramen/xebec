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

callback = (
    CustomJS(
    args=dict(
        big_source=big_source,
        box_source=box_source,
        col_source=col_source,
        yr=p.y_range
    ),
    code="""
// https://stackoverflow.com/a/53660837
function median(numbers) {
    const sorted = numbers.slice().sort((a, b) => a - b);
    const middle = Math.floor(sorted.length / 2);

    if (sorted.length % 2 === 0) {
        return (sorted[middle - 1] + sorted[middle]) / 2;
    }

    return sorted[middle];
}

// https://stackoverflow.com/a/55297611
const asc = arr => arr.sort((a, b) => a - b);
function quantile(arr, q) {
    const sorted = asc(arr);
    const pos = (sorted.length - 1) * q;
    const base = Math.floor(pos);
    const rest = pos - base;
    if (sorted[base + 1] !== undefined) {
        return sorted[base] + rest * (sorted[base + 1] - sorted[base]);
    } else {
        return sorted[base];
    }
}

const data = big_source.data;
const colArray = data['column'];

const colPhyloArray = [];  // phylogenetic
const colDivMetricArray = [];  // diversity_metric
const colESArray = [];  // effect_size
const colCompArray = []; // comparison
const colColorArray = []; // color
const colMarkerArray = []; // marker

const compCompObj = {};

for (let i = 0; i < colArray.length; i++) {
    if (colArray[i] == cb_obj.value) {
        let phylo = data['phylogenetic'][i];
        let comp = data['comparison'][i];
        let divMetric = data['diversity_metric'][i];
        let effectSize = data['effect_size'][i];
        let color = data['color'][i];
        let marker = data['marker'][i];

        if (comp in compCompObj) {
            compCompObj[comp].push(effectSize)
        } else {
            compCompObj[comp] = [effectSize];
        }

        colPhyloArray.push(phylo);
        colCompArray.push(comp);
        colDivMetricArray.push(divMetric);
        colESArray.push(effectSize);
        colColorArray.push(color);
        colMarkerArray.push(marker);
    } else {
    }
}

const uniqueComps = [...new Set(colCompArray)];

const newBoxSource = {
    'comparison': [],
    'index': [],
    'upper': [],
    'lower': [],
    'q1': [],
    'q2': [],
    'q3': [],
    'qmax': [],
    'qmin': [],
};

var i = 1;
for (let [key, values] of Object.entries(compCompObj)) {
    newBoxSource['comparison'].push(key);
    let qmax = Math.max(...values);
    let qmin = Math.min(...values);
    let q1 = quantile(values, 0.25);
    let q3 = quantile(values, 0.75);
    newBoxSource['qmax'].push(qmax);
    newBoxSource['qmin'].push(qmin);
    newBoxSource['q1'].push(q1);
    newBoxSource['q2'].push(median(values));
    newBoxSource['q3'].push(q3);
    newBoxSource['upper'].push(Math.max(q3, qmax));
    newBoxSource['lower'].push(Math.min(q1, qmin));

    newBoxSource['index'].push(i);
    i = i + 1;
}

box_source.data = newBoxSource;
col_source.data = {
    'phylogenetic': colPhyloArray,
    'diversity_metric': colDivMetricArray,
    'color': colColorArray,
    'marker': colMarkerArray,
    'effect_size': colESArray,
    'comparison': colCompArray
};
// Need to update y_range and factors
yr.end = i - 1;
yr.factors = newBoxSource['comparison'];
box_source.change.emit();
col_source.change.emit();
"""
    )
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
