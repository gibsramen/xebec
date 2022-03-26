from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import output_file, save, figure
import seaborn as sns


rule plot_effect_sizes:
    input:
        "results/{diversity_type}/all_metrics_effect_sizes.tsv"
    output:
        "results/{diversity_type}/effect_size_plot.html"
    run:
        if wildcards.diversity_type == "beta_div":
            non_phylo_metrics = beta_metrics.query("phylogenetic == 'non_phylo'")["diversity_metric"]
            phylo_metrics = beta_metrics.query("phylogenetic == 'phylo'")["diversity_metric"]

            diversity_metric_order = beta_metrics["diversity_metric"].tolist()
        else:
            non_phylo_metrics = alpha_metrics.query("phylogenetic == 'non_phylo'")["diversity_metric"]
            phylo_metrics = alpha_metrics.query("phylogenetic == 'phylo'")["diversity_metric"]

            diversity_metric_order = alpha_metrics["diversity_metric"].tolist()

        palette = dict(zip(
            non_phylo_metrics,
            sns.color_palette("Reds", len(non_phylo_metrics)).as_hex()
        ))
        palette.update(dict(zip(
            phylo_metrics,
            sns.color_palette("Blues", len(phylo_metrics)).as_hex()
        )))

        df = pd.read_table(input[0], sep="\t")
        output_file(output[0])
        p1 = generate_interactive(df, palette, diversity_metric_order, "cohens_d")
        p2 = generate_interactive(df, palette, diversity_metric_order, "cohens_f")
        layout = gridplot([[p1, p2]], sizing_mode="scale_width",
                          toolbar_location="right")
        save(layout, title="xebec")


def generate_interactive(df, palette, diversity_metric_order, metric="cohens_d"):
    """Generated boxplots for different diversity metrics."""
    _df = df.copy().query("metric == @metric")
    order = list(
        _df.groupby("column")
        .median()
        .sort_values(by="effect_size", ascending=False).index
    )

    # https://stackoverflow.com/a/27255567
    _df["column"] = _df["column"].astype("category")
    _df["column"] = _df["column"].cat.set_categories(order)
    _df = _df.sort_values(by="column")

    _df["diversity_metric"] = _df["diversity_metric"].astype("category")
    _df["diversity_metric"] = (
        _df["diversity_metric"]
        .cat
        .set_categories(diversity_metric_order)
    )

    cols = _df["column"].unique()
    div_metrics = _df["diversity_metric"].unique()

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

    p = figure(
        tools=["pan", "reset", "box_zoom", hover_points, hover_boxes],
        y_range=order,
        width=800
    )
    p.background_fill_color = "#EEEEEE"

    gb = _df.groupby("column")["effect_size"]
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

    lw = 2

    box_df = (
        _df.groupby("column").median()
        .assign(q1=q1, q2=q2, q3=q3, qmin=qmin, qmax=qmax, upper=upper, lower=lower)
        .reset_index()
    )
    box_source = ColumnDataSource(box_df)
    box_args = {
        "source": box_source,
        "y": "column",
        "fill_color": "white",
        "line_color": "black",
        "line_width": lw,
        "height": 0.7,
        "name": "boxes"
    }
    boxes_1 = p.hbar(**box_args, left="q1", right="q2")
    boxes_2 = p.hbar(**box_args, left="q2", right="q3")

    seg_args = {"source": box_source, "y0": "column", "y1": "column",
                "line_color": "black", "line_width": lw}
    seg_1 = p.segment(**seg_args, x0="upper", x1="q3")
    seg_2 = p.segment(**seg_args, x0="lower", x1="q1")

    whisker_args = {"source": box_source, "y": "column", "height": 0.5, "width": 0.00001,
                    "line_color": "black", "line_width": lw}
    whisk_1 = p.rect(**whisker_args, x="lower")
    whisk_2 = p.rect(**whisker_args, x="upper")

    # https://stackoverflow.com/a/58620263
    for patch in [boxes_1, boxes_2, seg_1, seg_2, whisk_1, whisk_2]:
        patch.level = "underlay"

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

    p.grid[0].level = "image"
    p.grid[1].level = "image"

    if metric == "cohens_d":
        p.xaxis.axis_label = "Cohen's d"
        p.title = "Binary Categories"
    else:
        p.title = "Multi-Class Categories"
        p.xaxis.axis_label = "Cohen's f"

    for ax in [p.xaxis, p.yaxis]:
        ax.axis_label_text_font_size = "15pt"
        ax.axis_label_text_font_style = "normal"
        ax.major_tick_line_width = 0

    p.title.text_font_size = "20pt"
    p.yaxis.major_label_text_font_size = "12pt"
    p.xaxis.major_label_text_font_size = "10pt"

    return p
