"""
Plots two clustered heatmaps summarizing TF cluster specificity and DE calls.

1) Cluster specificity heatmap: -log10(qval_cluster) per TF and cluster.
2) DE calls heatmap: isTF_DE per TF and cluster, ordered as in (1).

The script expects an Excel file with at least the following columns:
    - 'alias'
    - 'TF'
    - 'cluster'
    - 'qval_cluster'
    - 'isTF_DE'
    - 'hasTFrelevantGOterm'
    - 'borda_rank'
"""

import math
import sys

import natsort
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

REGULONS_FILE   = sys.argv[1]       # regulons Excel file
OUT_SPEC_PREFIX = sys.argv[2]       # output prefix (no extension) for cluster specificity heatmap
OUT_DE_PREFIX   = sys.argv[3]       # output prefix (no extension) for DE calls heatmap
TOP_N_TFS       = int(sys.argv[4])  # number of top TFs (by borda_rank) to plot


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------


def simplify_cluster_name(cluster_name):
    """If tissue name and cluster name are identical, simplify the redundant name."""
    if "_Cluster_" not in cluster_name:
        return cluster_name

    tissue, cluster = cluster_name.split("_Cluster_", 1)
    if tissue == cluster:
        return f"Cluster-{cluster}"
    return cluster_name.replace("_Cluster_", "-")


def compute_min_nonzero_qval(df):
    """Return the smallest non-zero q-value in the dataframe."""
    non_zero = df.loc[df["qval_cluster"] > 0, "qval_cluster"]
    if non_zero.empty:
        return 1e-300  # extreme fallback, should almost never happen
    return float(non_zero.min())


def build_heatmap_matrices(df, top_n):
    """
    Build matrices for the two heatmaps and the row color annotation.

    Returns
    -------
    plot_qval : DataFrame, -log10(qval_cluster) per alias x cluster
    plot_de   : DataFrame, isTF_DE per alias x cluster
    row_colors: Series, color per alias (used as row_colors in seaborn.clustermap)
    """
    # work on a copy to avoid modifying the original dataframe
    df = df.copy()

    # harmonize cluster names
    df["cluster"] = df["cluster"].astype(str).apply(simplify_cluster_name)

    # select the top TFs by borda_rank
    df_top = df.sort_values("borda_rank", ascending=True).head(top_n)
    top_aliases = sorted(df_top["alias"].unique().tolist())

    # sorted cluster list
    clusters = natsort.natsorted(df["cluster"].unique(), alg=natsort.IGNORECASE)

    # initialize matrices
    plot_qval = pd.DataFrame(0.0, index=top_aliases, columns=clusters)
    plot_de = pd.DataFrame(float("nan"), index=top_aliases, columns=clusters)

    # compute smallest non-zero q-value for replacing zeros
    min_nonzero_q = compute_min_nonzero_qval(df)

    # colors for GO/TF categories
    go_color_map = {
        "relevant_known_TF": "green",
        "known_TF": "orange",
    }
    default_color = "gray"

    # fill matrices
    top_aliases_set = set(top_aliases)
    for _, row in df.iterrows():
        alias = row["alias"]
        if alias not in top_aliases_set:
            continue

        cluster = row["cluster"]
        qval = float(row["qval_cluster"])

        # cluster specificity: -log10(qval_cluster)
        if qval > 0:
            value = -math.log10(qval)
        else:
            value = -math.log10(min_nonzero_q)
        plot_qval.at[alias, cluster] = round(value)

        # DE calls
        plot_de.at[alias, cluster] = row["isTF_DE"]

    # row colors: one color per alias
    row_color_info = (
        df.loc[df["alias"].isin(top_aliases), ["alias", "hasTFrelevantGOterm"]]
        .drop_duplicates("alias")
        .set_index("alias")
    )
    row_color_info["color"] = row_color_info["hasTFrelevantGOterm"].map(
        go_color_map
    ).fillna(default_color)

    # ensure same ordering as matrices
    row_colors = row_color_info["color"].reindex(index=plot_qval.index)

    return plot_qval, plot_de, row_colors


def compute_font_size(n_rows, n_cols, min_size=3, max_size=12):
    """
    Compute a font size that scales with the matrix dimensions.

    The bigger the matrix, the smaller the labels.
    """
    max_dim = max(n_rows, n_cols)

    if max_dim <= 20:
        return max_size
    if max_dim <= 50:
        return 9
    if max_dim <= 100:
        return 7
    if max_dim <= 200:
        return 5
    if max_dim <= 400:
        return 4
    return min_size


def compute_figsize(n_rows, n_cols, base_width=6, base_height=6, max_width=20, max_height=20):
    """
    Compute a reasonable figure size based on the matrix dimensions.
    Width scales with the number of columns, height scales with the number of rows.
    """
    width = base_width + n_cols * 0.2
    height = base_height + n_rows * 0.05

    width = min(width, max_width)
    height = min(height, max_height)

    return (width, height)


def add_go_legend_for_de(cg_de):
    """
    Add legend describing GO term categories and DE statuses to the DE heatmap.
    The legend is placed just above the heatmap, with a light gray frame.
    """
    legend_colors = {
        "Associated with GO term of interest":        "#008000ff",  # green
        "Associated with any other GO term":          "#ffa500ff",  # orange
        "No GO info known":                           "#808080ff",  # gray
        "Differentially expressed (DE) TF":           "#08306bff",  # dark blue
        "Expressed (in 10% of cells) TF, not DE":     "#f7fbffff",  # white
        "Not DE, not expressed":                      "#b3b3b3ff",  # background / NA
    }
    handles = [
        Patch(facecolor=color, linewidth=0.5, edgecolor="#303030ff")
        for color in legend_colors.values()
    ]
    leg = cg_de.ax_heatmap.legend(
        handles,
        legend_colors.keys(),
        title="Legend",
        loc="lower center",
        bbox_to_anchor=(0.5, 1.06),               # just above heatmap
        bbox_transform=cg_de.ax_heatmap.transAxes,
        borderaxespad=0.0,
        frameon=True,
    )
    leg.get_frame().set_edgecolor("gray")
    leg.get_frame().set_linewidth(0.8)
    leg.get_frame().set_facecolor("white")


# -------------------------------------------------------------------------
# Main plotting logic
# -------------------------------------------------------------------------

# load input table
regulons_df = pd.read_excel(REGULONS_FILE)

# prepare matrices for plotting
plot_qval, plot_de, row_colors = build_heatmap_matrices(regulons_df, TOP_N_TFS)
plot_qval = plot_qval.loc[:, (plot_qval != 0).any(axis=0)]  # remove columns that are all zeros in qval matrix
plot_qval[plot_qval > 20] = 20  # cap very high -log10(q) values for better color scaling

# derive final dimensions
n_rows, n_cols = plot_qval.shape
font_size = compute_font_size(n_rows, n_cols)
figsize = compute_figsize(n_rows, n_cols)


# ------------------------------------------------------------------
# Cluster specificity heatmap
# ------------------------------------------------------------------

sns.set(context="notebook", style="white")
cluster_grid_spec = sns.clustermap(
    plot_qval,
    cmap="Blues",
    xticklabels=True,
    yticklabels=True,
    row_colors=row_colors,
    col_cluster=False,
    row_cluster=True,
    mask=(plot_qval == 0),
    linewidths=0.01,
    linecolor="dimgrey",
    figsize=figsize,
)

heat_ax = cluster_grid_spec.ax_heatmap
heat_ax.set_facecolor("#b3b3b3ff")  # background color begind the heatmap cells
cluster_grid_spec.ax_cbar.set_title("-log10(qval)")

# adjust tick label sizes
heat_ax.set_xticklabels(heat_ax.get_xticklabels(), rotation=90, ha="center", fontsize=font_size)
heat_ax.set_yticklabels(heat_ax.get_yticklabels(), fontsize=font_size)

# tidy row color axis
cluster_grid_spec.ax_row_colors.tick_params(bottom=False)
cluster_grid_spec.ax_row_colors.set_xticklabels("")

# force a draw so that positions are known
fig_spec = cluster_grid_spec.fig
fig_spec.canvas.draw()
renderer = fig_spec.canvas.get_renderer()

# place legend above the heatmap (centered)
legend_colors_spec = {
    "Associated with GO term of interest": "#008000ff",
    "Associated with any other GO term":   "#ffa500ff",
    "No GO info known":                    "#808080ff",
}
# list of colored square patches
handles_spec = [
    Patch(facecolor=color, linewidth=0.5, edgecolor="#303030ff")
    for color in legend_colors_spec.values()
]
legend_spec = heat_ax.legend(
    handles_spec,
    legend_colors_spec.keys(),
    title="Legend",
    loc="lower center",
    bbox_to_anchor=(0.5, 1.06),
    bbox_transform=heat_ax.transAxes,
    frameon=True,
)
legend_spec.get_frame().set_edgecolor("gray")
legend_spec.get_frame().set_linewidth(0.8)
legend_spec.get_frame().set_facecolor("white")

# compute a colorbar with same height as the legend, to the left of it
legend_bbox_pixels = legend_spec.get_window_extent(renderer=renderer)
legend_bbox_fig = legend_bbox_pixels.transformed(fig_spec.transFigure.inverted())

cbar_width = 0.015   # narrow vertical bar
gap = 0.05           # small gap between colorbar and legend

cbar_x0 = legend_bbox_fig.x0 - cbar_width - gap
cbar_y0 = legend_bbox_fig.y0
cbar_h = legend_bbox_fig.height

cbar_ax = cluster_grid_spec.ax_cbar
cbar_ax.set_position([cbar_x0, cbar_y0, cbar_width, cbar_h])

# deduce row order from dendrogram to reuse in DE heatmap
row_order = [plot_qval.index[i] for i in cluster_grid_spec.dendrogram_row.reordered_ind]


# ------------------------------------------------------------------
# DE heatmap (same order as specificity heatmap)
# ------------------------------------------------------------------

plot_de = plot_de.reindex(row_order)
plot_de = plot_de.loc[:, plot_qval.columns]  # ensure same column ordering
n_rows_de, n_cols_de = plot_de.shape
font_size_de = compute_font_size(n_rows_de, n_cols_de)
figsize_de = compute_figsize(n_rows_de, n_cols_de)

cluster_grid_de = sns.clustermap(
    plot_de,
    cmap="Blues",
    xticklabels=True,
    yticklabels=True,
    row_colors=row_colors.reindex(row_order),
    col_cluster=False,
    row_cluster=False,
    vmax=1,
    vmin=0,
    linewidths=0.01,
    linecolor="dimgrey",
    figsize=figsize_de,
)

de_heat_ax = cluster_grid_de.ax_heatmap
de_heat_ax.set_facecolor("#b3b3b3ff")

# remove colorbar (DE is binary)
cluster_grid_de.cax.set_visible(False)

de_heat_ax.set_xticklabels(de_heat_ax.get_xticklabels(), rotation=90, ha="center", fontsize=font_size_de)
de_heat_ax.set_yticklabels(heat_ax.get_yticklabels(), fontsize=font_size_de)  # keep same order as first plot

# tidy row color axis
cluster_grid_de.ax_row_colors.tick_params(bottom=False)
cluster_grid_de.ax_row_colors.set_xticklabels("")

# legend for DE heatmap (above heatmap, centered)
add_go_legend_for_de(cluster_grid_de)


# ------------------------------------------------------------------
# Save figures
# ------------------------------------------------------------------
cluster_grid_spec.savefig(f"{OUT_SPEC_PREFIX}.svg")
cluster_grid_spec.savefig(f"{OUT_SPEC_PREFIX}.png", dpi=600)
cluster_grid_de.savefig(f"{OUT_DE_PREFIX}.svg")
cluster_grid_de.savefig(f"{OUT_DE_PREFIX}.png", dpi=600)
