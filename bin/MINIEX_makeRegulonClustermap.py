"""
Plots a clustered heatmap summarizing regulon sizes (# target genes) per cluster.

Inputs:
    1) CELLID_FILE:    tab-separated; at least:
                       - column 1: cluster ID (e.g. "1")
                       - column 2: cell type / tissue name
    2) REGULONS_FILE:  tab-separated; at least:
                       - column 1: regulon/TF name
                       - column 2: cluster label (e.g. "Cluster_1")
                       - column 3: comma-separated target genes
    3) OUT_PREFIX:     output file prefix (no extension)
"""

import sys

import natsort
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MaxNLocator


CELLID_FILE   = sys.argv[1]
REGULONS_FILE = sys.argv[2]
OUT_PREFIX    = sys.argv[3]


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------


def simplify_cluster_name(cluster_name: str) -> str:
    """If tissue name and cluster name are identical, simplify the redundant name."""
    if "_Cluster_" not in cluster_name:
        return cluster_name

    tissue, cluster = cluster_name.split("_Cluster_", 1)
    if tissue == cluster:
        return f"Cluster-{cluster}"
    return cluster_name.replace("_Cluster_", "-")


def load_cell_type_mapping(cellid_path: str) -> dict:
    """
    Load mapping from cluster ID to cell type / tissue name.

    Expects a tab-separated file with:
        - column 1: cluster ID (string)
        - column 2: cell type / tissue name
    """
    cluster_to_cell_type = {}
    with open(cellid_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            cluster_id = str(parts[0])
            cell_type = parts[1]
            cluster_to_cell_type[cluster_id] = cell_type
    return cluster_to_cell_type


def build_regulon_size_matrix(regulons_path: str, cluster_to_cell_type: dict) -> pd.DataFrame:
    """
    Build a regulon size matrix:
        rows   = regulons
        cols   = cell type + cluster label
        values = # target genes (TGs) per regulon/cluster
    """
    grouped = {}

    with open(regulons_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            regulon_name = parts[0]
            cluster_label = parts[1]        # e.g. "Cluster_1"
            targets_str = parts[2]          # comma-separated TGs

            # cluster ID used in CELLID_FILE
            cluster_id = cluster_label.replace("Cluster_", "")
            cell_type = cluster_to_cell_type[cluster_id]

            key = f"{cell_type}_{cluster_label}"
            n_targets = len(targets_str.split(",")) if targets_str else 0

            if key not in grouped:
                grouped[key] = {}
            grouped[key][regulon_name] = n_targets

    df = pd.DataFrame(grouped).fillna(0)  # rows = regulons, columns = celltype_cluster

    rename_map = {col: simplify_cluster_name(col) for col in df.columns}  # clean up cluster names
    df.rename(columns=rename_map, inplace=True)

    sorted_cols = natsort.natsorted(df.columns, alg=natsort.IGNORECASE)  # natural sort columns
    df = df.reindex(sorted_cols, axis=1)

    return df


def build_column_colors(df: pd.DataFrame, cluster_to_cell_type: dict) -> pd.DataFrame:
    """
    Build a color annotation DataFrame for columns (clusters).
    Uses a 'Spectral_r' palette over unique cell types.
    """
    unique_cell_types = natsort.natsorted(
        list(set(cluster_to_cell_type.values())),
        alg=natsort.IGNORECASE,
    )

    palette = sns.color_palette("Spectral_r", len(unique_cell_types))
    colors = palette.as_hex()
    cell_type_to_color = dict(zip(unique_cell_types, colors))

    col_colors_data = []
    for col in df.columns:
        if col.startswith("Cluster"):
            tissue_or_id = col.split("-", 1)[1]  # "Cluster-X" -> X
        else:
            tissue_or_id = col.rsplit("-", 1)[0]  # "Tissue-X" -> Tissue
        color = cell_type_to_color[tissue_or_id]
        col_colors_data.append([col, color])

    col_colors = pd.DataFrame(col_colors_data, columns=["cluster", "cell type"])
    col_colors = col_colors.set_index("cluster")
    return col_colors


def get_x_font_size_for_heatmap(df: pd.DataFrame) -> int:
    """X-axis labels: a bit larger overall, still shrink with many columns."""
    n_cols = len(df.columns)
    if n_cols > 80:
        return 8
    if n_cols > 60:
        return 9
    if n_cols > 40:
        return 10
    if n_cols > 25:
        return 12
    return 14


def compute_figure_size(df: pd.DataFrame):
    """
    Compute figure size so that:
        - for small matrices, the figure is wide and not too tall
        - cells are at most square (never tall, skinny rectangles)
        - height can grow with number of rows but is capped at 2 * width
    """
    n_rows, n_cols = df.shape

    # make width depend on number of columns, but within sane bounds
    per_col = 0.35
    width = n_cols * per_col
    width = max(10.0, min(20.0, width))  # between 10 and 20 inches

    # target: height/width ~ rows/cols -> cells ~square
    if n_cols == 0:
        height = 5.0
    else:
        height = width * (float(n_rows) / float(n_cols))

    # minimal height for tiny matrices
    if height < 4.0:
        height = 4.0

    # cap height at 2 * width
    max_height = 2.0 * width
    if height > max_height:
        height = max_height

    return width, height


def plot_clustermap(df: pd.DataFrame, column_colors: pd.DataFrame, out_prefix: str) -> None:
    """Plot clustermap with layout tweaks and save as SVG + PNG."""
    sns.set_style("white")

    fig_width, fig_height = compute_figure_size(df)

    g = sns.clustermap(
        df,
        cmap="mako_r",
        yticklabels=False,
        xticklabels=True,
        col_colors=column_colors,
        figsize=(fig_width, fig_height),
    )
    g.ax_cbar.set_title("#TGs")
    g.ax_cbar.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))

    # --- Heatmap & x labels ------------------------------------------------
    heatmap_ax = g.ax_heatmap
    heatmap_pos = heatmap_ax.get_position()

    # keep full width
    heatmap_ax.set_position([heatmap_pos.x0, heatmap_pos.y0, heatmap_pos.width, heatmap_pos.height])

    labels = heatmap_ax.set_xticklabels(
        heatmap_ax.get_xticklabels(),
        rotation=90,
        horizontalalignment="center",
        verticalalignment="top",
        fontsize=get_x_font_size_for_heatmap(df),
    )
    for label in labels:
        x, y = label.get_position()
        label.set_y(y - 0.02)  # extra spacing so labels do not overlap colored squares

    # --- Column color boxes ------------------------------------------------
    color_box = g.ax_col_colors.get_position()
    color_box.y0 = heatmap_pos.y0
    color_box.y1 = color_box.y0 - 0.013
    g.ax_col_colors.set_position(
        [color_box.x0, color_box.y0, color_box.width, color_box.height]
    )
    g.ax_col_colors.tick_params(right=False)
    g.ax_col_colors.set_yticklabels("")

    # --- Column dendrogram: limit its height -------------------------------
    dendro_box = g.ax_col_dendrogram.get_position()

    min_dendro_height_in_inches = 1
    min_dendro_height_norm = min_dendro_height_in_inches / fig_height  # height is relative to figure size
    max_dendro_height = heatmap_pos.height * 0.15  # at most 15% of heatmap height
    new_dendro_height = max(min_dendro_height_norm, max_dendro_height)

    g.ax_col_dendrogram.set_position([dendro_box.x0, heatmap_pos.y1, dendro_box.width, new_dendro_height])

    # --- Colorbar: align vertically with dendrogram, with min/max height ---
    cbar_box = g.ax_cbar.get_position()
    dendro_pos = g.ax_col_dendrogram.get_position()
    heatmap_pos = g.ax_heatmap.get_position()  # re-fetch if needed

    dendro_height = dendro_pos.height

    # min and max allowed colorbar height
    max_cbar_height = 0.4
    new_cbar_height = min(new_dendro_height * 0.9, max_cbar_height)

    # center the colorbar on the dendrogram band
    center_y = dendro_pos.y0 + dendro_height / 2.0
    y0 = center_y - new_cbar_height / 2.0

    g.ax_cbar.set_position(
        [
            cbar_box.x0 + 0.05,   # keep your horizontal tweak
            y0,
            cbar_box.width * 0.5, # keep your width scaling
            new_cbar_height,
        ]
    )

    g.savefig(f"{out_prefix}.svg")
    g.savefig(f"{out_prefix}.png", dpi=600)


cluster_to_cell_type = load_cell_type_mapping(CELLID_FILE)
regulon_df = build_regulon_size_matrix(REGULONS_FILE, cluster_to_cell_type)
column_colors = build_column_colors(regulon_df, cluster_to_cell_type)
plot_clustermap(regulon_df, column_colors, OUT_PREFIX)
