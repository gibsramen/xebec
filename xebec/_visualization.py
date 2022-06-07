from bokeh.models import ColumnDataSource, HoverTool
import seaborn as sns


# Figure parameters
TITLE_FONT_SIZE = "20pt"
Y_AXIS_MAJOR_TICK_LABEL_FONT_SIZE = "12pt"
X_AXIS_MAJOR_TICK_LABEL_FONT_SIZE = "10pt"
AXIS_LABEL_FONT_SIZE = "15pt"
AXIS_LABEL_FONT_STYLE = "normal"
AXIS_MAJOR_TICK_LINE_WIDTH = 0
BACKGROUND_FILL_COLOR = "#EEEEEE"

# Boxplot parameters
BOXPLOT_LINE_WIDTH = 2
BOXPLOT_LINE_COLOR = "black"
BOXPLOT_WHISKER_HEIGHT = 0.5
BOXPLOT_WHISKER_WIDTH = 0.00001
BOXPLOT_BOX_HEIGHT = 0.7
BOXPLOT_BOX_FILL_COLOR = "white"

# Scatter parameters
SCATTER_LINE_WIDTH = 0.5
SCATTER_POINT_SIZE = 10
SCATTER_PHYLO_PALETTE = "Blues"
SCATTER_PHYLO_MARKER = "triangle"
SCATTER_NON_PHYLO_PALETTE = "Reds"
SCATTER_NON_PHYLO_MARKER = "circle"

HOVER_POINTS = HoverTool(mode="mouse", names=["points"], attachment="below")
HOVER_POINTS.tooltips = [
    ("Effect Size", "@effect_size{0.000}"),
    ("Diversity Metric", "@diversity_metric")
]

HOVER_BOXES = HoverTool(mode="mouse", names=["boxes"], attachment="above")
HOVER_BOXES.tooltips = [
    ("25%", "@q1{0.000}"),
    ("50%", "@q2{0.000}"),
    ("75%", "@q3{0.000}")
]


def add_boxplots(figure, all_metrics_es_df, group_var) -> None:
    """Add boxplots to figure."""
    gb = all_metrics_es_df.groupby(group_var)["effect_size"]
    q1 = gb.quantile(q=0.25)
    q2 = gb.quantile(q=0.50)
    q3 = gb.quantile(q=0.75)
    iqr = q3 - q1

    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    qmin = gb.quantile(q=0)
    qmax = gb.quantile(q=1)

    upper = [min([x, y]) for (x, y) in zip(list(qmax), upper)]
    lower = [max([x, y]) for (x, y) in zip(list(qmin), lower)]

    box_df = (
        all_metrics_es_df.groupby(group_var)
        .median()
        .assign(q1=q1, q2=q2, q3=q3,
                qmin=qmin, qmax=qmax,
                upper=upper, lower=lower)
        .reset_index()
    )
    box_source = ColumnDataSource(box_df)

    box_args = {
        "source": box_source,
        "y": group_var,
        "fill_color": BOXPLOT_BOX_FILL_COLOR,
        "line_color": BOXPLOT_LINE_COLOR,
        "line_width": BOXPLOT_LINE_WIDTH,
        "height": BOXPLOT_BOX_HEIGHT,
        "name": "boxes"
    }
    boxes_1 = figure.hbar(**box_args, left="q1", right="q2")
    boxes_2 = figure.hbar(**box_args, left="q2", right="q3")

    seg_args = {
        "source": box_source,
        "y0": group_var,
        "y1": group_var,
        "line_color": BOXPLOT_LINE_COLOR,
        "line_width": BOXPLOT_LINE_WIDTH
    }
    seg_1 = figure.segment(**seg_args, x0="upper", x1="q3")
    seg_2 = figure.segment(**seg_args, x0="lower", x1="q1")

    whisker_args = {
        "source": box_source,
        "y": group_var,
        "height": BOXPLOT_WHISKER_HEIGHT,
        "width": BOXPLOT_WHISKER_WIDTH,
        "line_color": BOXPLOT_LINE_COLOR,
        "line_width": BOXPLOT_LINE_WIDTH
    }
    whisk_1 = figure.rect(**whisker_args, x="lower")
    whisk_2 = figure.rect(**whisker_args, x="upper")

    # https://stackoverflow.com/a/58620263
    for patch in [boxes_1, boxes_2, seg_1, seg_2, whisk_1, whisk_2]:
        patch.level = "underlay"

    return box_source


def get_scatter_palette(phylo_metrics, non_phylo_metrics) -> dict:
    palette = dict(zip(
        non_phylo_metrics,
        sns.color_palette(
            SCATTER_NON_PHYLO_PALETTE,
            len(non_phylo_metrics)
        ).as_hex()
    ))
    palette.update(dict(zip(
        phylo_metrics,
        sns.color_palette(
            SCATTER_PHYLO_PALETTE,
            len(phylo_metrics)
        ).as_hex()
    )))
    return palette


def assign_fig_parameters(figure):
    figure.grid[0].level = "image"
    figure.grid[1].level = "image"
    figure.background_fill_color = BACKGROUND_FILL_COLOR

    figure.title.text_font_size = TITLE_FONT_SIZE
    figure.xaxis.major_label_text_font_size = X_AXIS_MAJOR_TICK_LABEL_FONT_SIZE
    figure.yaxis.major_label_text_font_size = Y_AXIS_MAJOR_TICK_LABEL_FONT_SIZE

    for ax in [figure.xaxis, figure.yaxis]:
        ax.axis_label_text_font_size = AXIS_LABEL_FONT_SIZE
        ax.axis_label_text_font_style = AXIS_LABEL_FONT_STYLE
        ax.major_tick_line_width = AXIS_MAJOR_TICK_LINE_WIDTH
