from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, HoverTool, Select
from bokeh.models.callbacks import CustomJS
from bokeh.plotting import output_file, save, figure
import pandas as pd
import seaborn as sns


diversity_metric_order = snakemake.params["all_div_metrics"]
non_phylo_metrics = snakemake.params["non_phylo_metrics"]
phylo_metrics = snakemake.params["phylo_metrics"]

palette = dict(zip(
    non_phylo_metrics,
    sns.color_palette("Reds", len(non_phylo_metrics)).as_hex()
))
palette.update(dict(zip(
    phylo_metrics,
    sns.color_palette("Blues", len(phylo_metrics)).as_hex()
)))

df = pd.read_table(snakemake.input[0], sep="\t")
df["comparison"] = df["group_1"] + " vs. " + df["group_2"]
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

gb = col_df.groupby("comparison")["effect_size"]
q1 = gb.quantile(q=0.25)
q2 = gb.quantile(q=0.5)
q3 = gb.quantile(q=0.75)
iqr = q3 - q1
upper = q3 + 1.5*iqr
lower = q1 - 1.5*iqr

qmin = gb.quantile(q=0)
qmax = gb.quantile(q=1)

upper = [min([x, y]) for (x, y) in zip(list(qmax), upper)]
lower = [max([x, y]) for (x, y) in zip(list(qmin), lower)]

box_df = (
    col_df.groupby("comparison").median()
    .assign(q1=q1, q2=q2, q3=q3, qmin=qmin, qmax=qmax, upper=upper, lower=lower)
    .reset_index()
)
order = list(
    col_df
    .groupby("comparison")
    .median()
    .sort_values(by="effect_size", ascending=False)
    .index
)

hover_points = HoverTool(mode="mouse", names=["points"], attachment="below")
hover_points.tooltips = [
    ("Effect Size", "@effect_size{0.000}"),
    ("Diversity Metric", "@diversity_metric")
]

hover_boxes = HoverTool(mode="mouse", names=["boxes"], attachment="above")
hover_boxes.tooltips = [
    ("25%", "@q1{0.000}"),
    ("50%", "@q2{0.000}"),
    ("75%", "@q3{0.000}")
]
hover_boxes.formatters = {"@q1": "numeral"}

p = figure(
    tools=["pan", "reset", "box_zoom", hover_boxes, hover_points],
    y_range=order,
    width=800,
)

big_source = ColumnDataSource(df)
box_source = ColumnDataSource(box_df)

callback = (
    CustomJS(
    args=dict(
        big_source=big_source,
        box_source=box_source,
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

const compPhyloArray = [];  // phylogenetic
const compDivMetricArray = [];  // diversity_metric
const compESArray = [];  // effect_size
const compCompArray = []; // comparison

const compCompObj = {};

for (let i = 0; i < colArray.length; i++) {
    if (colArray[i] == cb_obj.value) {
        let phylo = data['phylogenetic'][i];
        let comp = data['comparison'][i];
        let divMetric = data['phylogenetic'][i];
        let effectSize = data['effect_size'][i];

        if (comp in compCompObj) {
            compCompObj[comp].push(effectSize)
        } else {
            compCompObj[comp] = [effectSize];
        }
    } else {
    }
}

const uniqueComps = [...new Set(compCompArray)]; 

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

var i = 0;
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
// Need to update y_range and factors
yr.end = i - 1;
yr.factors = newBoxSource['comparison'];
box_source.change.emit();
"""
    )
)
chosen_col.js_on_change("value", callback)


p.background_fill_color = "#EEEEEE"
p.grid[0].level = "image"
p.grid[1].level = "image"

lw = 2

box_args = {
    "source": box_source,
    "y": "comparison",
    "fill_color": "white",
    "line_color": "black",
    "line_width": lw,
    "height": 0.7,
    "name": "boxes"
}
boxes_1 = p.hbar(**box_args, left="q1", right="q2")
boxes_2 = p.hbar(**box_args, left="q2", right="q3")

seg_args = {"source": box_source, "y0": "comparison", "y1": "comparison",
            "line_color": "black", "line_width": lw}
seg_1 = p.segment(**seg_args, x0="upper", x1="q3")
seg_2 = p.segment(**seg_args, x0="lower", x1="q1")

whisker_args = {"source": box_source, "y": "comparison", "height": 0.5, "width": 0.00001,
                "line_color": "black", "line_width": lw}
whisk_1 = p.rect(**whisker_args, x="lower")
whisk_2 = p.rect(**whisker_args, x="upper")

# https://stackoverflow.com/a/58620263
for patch in [boxes_1, boxes_2, seg_1, seg_2, whisk_1, whisk_2]:
    patch.level = "underlay"

p.xaxis.axis_label = "Cohen's d"
controls = [chosen_col]
control_panel = column(controls, width=200, height=200)
layout = row(control_panel, p)
output_file(snakemake.output[0])
save(layout, title="xebec")
