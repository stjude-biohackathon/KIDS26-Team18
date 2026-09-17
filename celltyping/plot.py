"""Plotting helpers extracted from Utils.py."""

from __future__ import annotations

import math

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import Normalize
from scipy import sparse
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist

def custom_barplot(
    adata,
    var_1,
    var_2,
    var_3=None,
    plot_type="bar",
    cluster_bars=True,
    cluster_metric="euclidean",
    cluster_method="average",
    cluster_within_var_3=True,
    gap=0.3,
    figsize=None,
    bar_width=0.90,
    label_every=None,
    max_labels=50,
    title=None,
    legend=True,
    show=True,
):
    """
    Plot stacked proportions of `var_1` across categories in `var_2`.

    Parameters
    ----------
    adata
        AnnData object.

    var_1 : str
        adata.obs column used as stacked bar fill categories.

    var_2 : str
        adata.obs column used as individual bars / x-axis categories.

    var_3 : str or None, default None
        Optional adata.obs column used to group bars.

    plot_type : {"bar", "circular"}, default "bar"
        Plot style.

    cluster_bars : bool, default True
        Order bars by similarity of their var_1 proportion profiles.

    cluster_metric : str, default "euclidean"
        Distance metric passed to scipy.spatial.distance.pdist.

    cluster_method : str, default "average"
        Linkage method passed to scipy.cluster.hierarchy.linkage.

    cluster_within_var_3 : bool, default True
        When var_3 is provided, cluster bars separately within each group.
        If False, cluster all var_2/var_3 combinations together.

    gap : float, default 1.2
        Space between var_3 groups in the standard bar plot.

    figsize : tuple or None
        Matplotlib figure size.

    bar_width : float, default 0.90
        Width of bars.

    label_every : int or None
        Display every nth x-axis label. If None, selected automatically.

    max_labels : int, default 50
        Maximum approximate number of labels when label_every is None.

    title : str or None
        Custom title.

    legend : bool, default True
        Whether to display the legend.

    show : bool, default True
        Whether to call plt.show().

    Returns
    -------
    plot_df : pandas.DataFrame
        Proportion table in plotted order.

    fig : matplotlib.figure.Figure

    ax : matplotlib.axes.Axes
    """

    # =====================================================
    # VALIDATE ARGUMENTS
    # =====================================================
    valid_plot_types = {"bar", "circular"}

    if plot_type not in valid_plot_types:
        raise ValueError(
            f"`plot_type` must be one of {valid_plot_types}. "
            f"Received: {plot_type!r}"
        )

    required_columns = [var_1, var_2]

    if var_3 is not None:
        required_columns.append(var_3)

    missing_columns = [
        column
        for column in required_columns
        if column not in adata.obs.columns
    ]

    if missing_columns:
        raise KeyError(
            f"Columns not found in adata.obs: {missing_columns}"
        )

    # Work on a copy to avoid changing adata.obs
    obs = adata.obs[required_columns].copy()

    # Remove rows missing required annotations
    obs = obs.dropna(subset=required_columns)

    if obs.empty:
        raise ValueError(
            "No observations remain after removing missing values."
        )

    # =====================================================
    # VAR_1 FILL ORDER
    # =====================================================
    if isinstance(adata.obs[var_1].dtype, pd.CategoricalDtype):
        observed_values = set(
            obs[var_1].astype(str).unique()
        )

        var_1_order = [
            str(category)
            for category in adata.obs[var_1].cat.categories
            if str(category) in observed_values
        ]
    else:
        var_1_order = (
            obs[var_1]
            .astype(str)
            .unique()
            .tolist()
        )

    obs[var_1] = pd.Categorical(
        obs[var_1].astype(str),
        categories=var_1_order,
        ordered=True,
    )

    # =====================================================
    # COUNTS, PROPORTIONS, AND PIVOT TABLE
    # =====================================================
    if var_3 is not None:
        df = (
            obs.groupby(
                [var_3, var_2, var_1],
                observed=True,
            )
            .size()
            .reset_index(name="n")
        )

        df["prop"] = (
            df.groupby(
                [var_3, var_2],
                observed=True,
            )["n"]
            .transform(lambda values: values / values.sum())
        )

        plot_df = df.pivot_table(
            index=[var_3, var_2],
            columns=var_1,
            values="prop",
            fill_value=0,
            observed=True,
        )

    else:
        df = (
            obs.groupby(
                [var_2, var_1],
                observed=True,
            )
            .size()
            .reset_index(name="n")
        )

        df["prop"] = (
            df.groupby(
                var_2,
                observed=True,
            )["n"]
            .transform(lambda values: values / values.sum())
        )

        plot_df = df.pivot_table(
            index=var_2,
            columns=var_1,
            values="prop",
            fill_value=0,
            observed=True,
        )

    plot_df = (
        plot_df
        .reindex(columns=var_1_order)
        .fillna(0)
    )

    # =====================================================
    # FUNCTION TO ORDER BARS BY SIMILARITY
    # =====================================================
    def cluster_bar_order(composition_df):
        if len(composition_df) <= 2:
            return composition_df.index.tolist()

        values = composition_df.to_numpy(dtype=float)

        distances = pdist(
            values,
            metric=cluster_metric,
        )

        if len(distances) == 0 or np.allclose(distances, 0):
            return composition_df.index.tolist()

        linkage_matrix = linkage(
            distances,
            method=cluster_method,
            optimal_ordering=True,
        )

        order = leaves_list(linkage_matrix)

        return composition_df.index[order].tolist()

    # =====================================================
    # CLUSTER AND REORDER BARS
    # =====================================================
    group_centers = []
    group_labels = []
    group_boundaries = []

    if cluster_bars:
        if var_3 is not None and cluster_within_var_3:
            groups = (
                plot_df.index
                .get_level_values(var_3)
                .unique()
                .tolist()
            )

            ordered_index = []

            for group in groups:
                group_df = plot_df.xs(
                    group,
                    level=var_3,
                )

                group_order = cluster_bar_order(group_df)

                ordered_index.extend(
                    (group, category)
                    for category in group_order
                )

            plot_df = plot_df.reindex(
                pd.MultiIndex.from_tuples(
                    ordered_index,
                    names=[var_3, var_2],
                )
            )

        else:
            clustered_order = cluster_bar_order(plot_df)
            plot_df = plot_df.reindex(clustered_order)

    # =====================================================
    # COLORS
    # =====================================================
    color_key = f"{var_1}_colors"

    if color_key in adata.uns:

        stored_colors = list(adata.uns[color_key])

        if isinstance(adata.obs[var_1].dtype, pd.CategoricalDtype):
            original_categories = [
                str(category)
                for category in adata.obs[var_1].cat.categories
            ]

            full_palette = dict(
                zip(original_categories, stored_colors)
            )

            palette = {}

            for i, category in enumerate(var_1_order):
                if category in full_palette:
                    palette[category] = full_palette[category]
                else:
                    # fallback color if palette is incomplete
                    palette[category] = plt.cm.tab20(
                        i / max(len(var_1_order) - 1, 1)
                    )

        else:
            if len(stored_colors) >= len(var_1_order):
                palette = dict(
                    zip(var_1_order, stored_colors)
                )
            else:
                palette = {
                    category: plt.cm.tab20(
                        i / max(len(var_1_order) - 1, 1)
                    )
                    for i, category in enumerate(var_1_order)
                }

    else:
        # No Scanpy palette found -> use default matplotlib colors
        cmap = plt.cm.get_cmap("tab20", len(var_1_order))

        palette = {
            category: cmap(i)
            for i, category in enumerate(var_1_order)
        }

    colors = [
        palette[category]
        for category in plot_df.columns
    ]

    # =====================================================
    # LABELS
    # =====================================================
    if var_3 is not None:
        bar_labels = [
            str(category)
            for _, category in plot_df.index
        ]

        combined_labels = [
            f"{group} | {category}"
            for group, category in plot_df.index
        ]
    else:
        bar_labels = plot_df.index.astype(str).tolist()
        combined_labels = bar_labels.copy()

    n_bars = len(plot_df)

    if n_bars == 0:
        raise ValueError("No bars are available to plot.")

    if label_every is None:
        label_every = max(
            1,
            int(np.ceil(n_bars / max_labels)),
        )

    # =====================================================
    # STANDARD STACKED BARPLOT
    # =====================================================
    if plot_type == "bar":
        x_positions = []
        x_labels = []

        if var_3 is not None:
            groups = (
                plot_df.index
                .get_level_values(var_3)
                .unique()
                .tolist()
            )

            current_x = 0

            for group_number, group in enumerate(groups):
                group_df = plot_df.xs(
                    group,
                    level=var_3,
                    drop_level=False,
                )

                start_x = current_x

                for _, category in group_df.index:
                    x_positions.append(current_x)
                    x_labels.append(str(category))
                    current_x += 1

                end_x = current_x - 1

                group_centers.append(
                    (start_x + end_x) / 2
                )
                group_labels.append(str(group))

                if group_number < len(groups) - 1:
                    group_boundaries.append(
                        end_x + (gap + 1) / 2
                    )

                current_x += gap

            x_positions = np.asarray(
                x_positions,
                dtype=float,
            )

        else:
            x_positions = np.arange(
                n_bars,
                dtype=float,
            )

            x_labels = bar_labels

        if figsize is None:
            figsize = (
                max(14, n_bars * 0.45),
                7,
            )

        fig, ax = plt.subplots(figsize=figsize)

        bottom = np.zeros(n_bars)

        for fill_category, color in zip(
            plot_df.columns,
            colors,
        ):
            values = plot_df[fill_category].to_numpy()

            ax.bar(
                x_positions,
                values,
                bottom=bottom,
                width=bar_width,
                label=fill_category,
                color=color,
                edgecolor="white",
                linewidth=0.25,
            )

            bottom += values

        visible_positions = np.arange(
            0,
            n_bars,
            label_every,
        )

        ax.set_xticks(
            x_positions[visible_positions]
        )

        ax.set_xticklabels(
            [
                x_labels[position]
                for position in visible_positions
            ],
            rotation=45,
            ha="right",
            fontsize=8,
        )

        if var_3 is not None:
            for center, group_label in zip(
                group_centers,
                group_labels,
            ):
                ax.text(
                    center,
                    -0.24,
                    group_label,
                    rotation=45,
                    ha="right",
                    va="top",
                    fontsize=10,
                    fontweight="bold",
                    transform=ax.get_xaxis_transform(),
                )

            for boundary in group_boundaries:
                ax.axvline(
                    boundary,
                    linewidth=0.6,
                    linestyle="--",
                    alpha=0.35,
                )

        ax.set_xlabel(var_2)
        ax.set_ylabel(f"Proportion of {var_1}")
        ax.set_ylim(0, 1)
        ax.margins(x=0.005)

        if title is None:
            if var_3 is not None:
                title = (
                    f"{var_1} composition across {var_2}, "
                    f"grouped by {var_3}"
                )
            else:
                title = (
                    f"{var_1} composition across {var_2}"
                )

        ax.set_title(title)

        if legend:
            ax.legend(
                title=var_1,
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                frameon=False,
            )

        plt.subplots_adjust(
            bottom=0.30 if var_3 is not None else 0.20,
            right=0.80 if legend else 0.95,
        )

    # =====================================================
    # CIRCULAR STACKED BARPLOT
    # =====================================================
    else:
        if figsize is None:
            figsize = (12, 12)

        fig, ax = plt.subplots(
            figsize=figsize,
            subplot_kw={"projection": "polar"},
        )

        angles = np.linspace(
            0,
            2 * np.pi,
            n_bars,
            endpoint=False,
        )

        circular_width = (
            2 * np.pi / n_bars
        ) * bar_width

        bottom = np.zeros(n_bars)

        for fill_category, color in zip(
            plot_df.columns,
            colors,
        ):
            values = plot_df[fill_category].to_numpy()

            ax.bar(
                angles,
                values,
                width=circular_width,
                bottom=bottom,
                label=fill_category,
                color=color,
                edgecolor="white",
                linewidth=0.15,
                align="edge",
            )

            bottom += values

        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_ylim(0, 1)

        ax.set_yticks(
            [0.25, 0.50, 0.75, 1.00]
        )

        ax.set_yticklabels(
            ["25%", "50%", "75%", "100%"]
        )

        ax.set_rlabel_position(0)

        visible_positions = np.arange(
            0,
            n_bars,
            label_every,
        )

        label_angles = (
            angles[visible_positions]
            + circular_width / 2
        )

        ax.set_xticks(label_angles)

        if var_3 is not None:
            visible_labels = [
                combined_labels[position]
                for position in visible_positions
            ]
        else:
            visible_labels = [
                bar_labels[position]
                for position in visible_positions
            ]

        ax.set_xticklabels(
            visible_labels,
            fontsize=7,
        )

        # Draw separators between var_3 groups
        if var_3 is not None:
            groups = (
                plot_df.index
                .get_level_values(var_3)
                .to_numpy()
            )

            transitions = np.where(
                groups[:-1] != groups[1:]
            )[0]

            for transition in transitions:
                separator_angle = angles[transition + 1]

                ax.plot(
                    [separator_angle, separator_angle],
                    [0, 1.03],
                    linewidth=1.2,
                    alpha=0.5,
                )

        if title is None:
            if var_3 is not None:
                title = (
                    f"Circular {var_1} composition across {var_2}, "
                    f"grouped by {var_3}"
                )
            else:
                title = (
                    f"Circular {var_1} composition across {var_2}"
                )

        ax.set_title(
            title,
            pad=30,
        )

        if legend:
            ax.legend(
                title=var_1,
                bbox_to_anchor=(1.15, 1.05),
                loc="upper left",
                frameon=False,
            )

        plt.tight_layout()

    if show:
        plt.show()

    return plot_df, fig, ax

# # usage
# plot_df, fig, ax = custom_barplot(
#     adata=adata,
#     var_1="leiden_res_0_3",  # stacked bar fill categories
#     var_2="core",            # individual bars
#     var_3=None,              # optional grouping
#     plot_type="bar",
# )
# plot_df, fig, ax = custom_barplot(
#     adata=adata,
#     var_1="leiden_res_0_3",  # stacked bar fill categories
#     var_2="core",            # individual bars
#     var_3=None,              # optional grouping
#     plot_type="circular",
#     max_labels=35,
#     figsize=(14, 14),
# )
# plot_df, fig, ax = custom_barplot(
#     adata=adata,
#     var_1="leiden_res_0_3",
#     var_2="core",
#     var_3="sample",
#     plot_type="circular",
#     cluster_within_var_3=True,
#     max_labels=35,
#     figsize=(14, 14),
# )

def dotplot(
    adata,
    marker_panels,
    groupby="cell_type",
    celltype_order=None,

    fig_width=42,
    fig_height=10,
    dpi=300,
    save_path=None,

    x_tick_fontsize=6,
    y_tick_fontsize=13,
    top_label_fontsize=10,
    top_label_rotation=90,
    top_label_y=-7.5,
    bracket_y=-1.25,
    bracket_height=0.25,
    bracket_linewidth=1.2,
    axis_label_fontsize=15,
    legend_fontsize=12,
    legend_title_fontsize=13,

    dot_min=0,
    dot_min_size=2,
    dot_max_size=45,
    size_legend_pcts=(0, 25, 50, 75, 100),

    cmap="magma",
    col_min=-2.5,
    col_max=2.5,

    expression_threshold=0,

    scale=True,
    scale_by="radius",

    use_expm1=True,

    cluster_celltypes=True,
    dendrogram_method="average",
    dendrogram_metric="correlation",
    dendrogram_width=0.10,
    keep_duplicate_genes=False,
):
    if scale_by not in ["radius", "size"]:
        raise ValueError("scale_by must be either 'radius' or 'size'")

    if celltype_order is None:
        celltype_order = (
            adata.obs[groupby]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

    present_celltypes = (
        adata.obs[groupby]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    celltype_order_present = [
        ct for ct in celltype_order
        if ct in present_celltypes
    ]

    marker_panels_present = {}

    for panel, genes in marker_panels.items():
        genes_present = [g for g in genes if g in adata.var_names]

        if len(genes_present) > 0:
            marker_panels_present[panel] = genes_present

    gene_order = []
    panel_spans = []
    start = 0

    for panel, genes in marker_panels_present.items():

        if keep_duplicate_genes:
            genes_here = genes.copy()
        else:
            genes_here = [g for g in genes if g not in gene_order]

        if len(genes_here) == 0:
            continue

        gene_order.extend(genes_here)

        end = start + len(genes_here) - 1
        panel_spans.append((panel, start, end))
        start = end + 1

    gene_order_unique = list(dict.fromkeys(gene_order))

    if len(gene_order) == 0:
        raise ValueError("None of the requested marker genes were found in adata.var_names.")

    if len(celltype_order_present) == 0:
        raise ValueError("None of the requested cell types were found in adata.obs[groupby].")

    adata_plot = adata[
        adata.obs[groupby].astype(str).isin(celltype_order_present),
        gene_order_unique
    ].copy()

    X = adata_plot.X

    if sparse.issparse(X):
        X = X.toarray()

    expr = pd.DataFrame(
        X,
        index=adata_plot.obs[groupby].astype(str).values,
        columns=gene_order_unique,
    )

    if use_expm1:
        expr_for_avg = np.expm1(expr)
    else:
        expr_for_avg = expr.copy()

    avg_exp = expr_for_avg.groupby(expr.index).mean()

    pct_exp = (
        expr.gt(expression_threshold)
        .groupby(expr.index)
        .mean()
    )

    gene_order_plot = gene_order_unique

    avg_exp = avg_exp.loc[celltype_order_present, gene_order_unique]
    pct_exp = pct_exp.loc[celltype_order_present, gene_order_unique]

    avg_exp_log = np.log1p(avg_exp)

    if scale and len(celltype_order_present) > 1:
        avg_exp_scaled = avg_exp_log.copy()

        for gene in gene_order_plot:
            vals = avg_exp_log[gene].values.astype(float)
            mean_val = np.nanmean(vals)
            std_val = np.nanstd(vals)

            if std_val > 0:
                scaled_vals = (vals - mean_val) / std_val
            else:
                scaled_vals = np.zeros_like(vals)

            scaled_vals = np.clip(scaled_vals, col_min, col_max)
            avg_exp_scaled[gene] = scaled_vals

        color_expr = avg_exp_scaled
        colorbar_label = "Scaled average expression"

    else:
        color_expr = avg_exp_log.copy()
        colorbar_label = "Average expression"

        col_min = np.nanmin(color_expr.values)
        col_max = np.nanmax(color_expr.values)

    # ------------------------------------------------------------
    # Cluster cell types
    # ------------------------------------------------------------
    Z = None

    if cluster_celltypes and len(celltype_order_present) > 1:
        cluster_matrix = color_expr.loc[celltype_order_present, gene_order_plot].values

        Z = linkage(
            cluster_matrix,
            method=dendrogram_method,
            metric=dendrogram_metric,
        )

        leaf_order = leaves_list(Z)

        celltype_order_present = [
            celltype_order_present[i] for i in leaf_order
        ]

        avg_exp = avg_exp.loc[celltype_order_present, gene_order_plot]
        avg_exp_log = avg_exp_log.loc[celltype_order_present, gene_order_plot]
        pct_exp = pct_exp.loc[celltype_order_present, gene_order_plot]
        color_expr = color_expr.loc[celltype_order_present, gene_order_plot]

    gene_order_plot = gene_order if keep_duplicate_genes else gene_order_unique

    # ------------------------------------------------------------
    # Build plotting dataframe
    # ------------------------------------------------------------
    genes_for_x = gene_order if keep_duplicate_genes else gene_order_unique
    plot_df = []

    for y_idx, celltype in enumerate(celltype_order_present):
        for x_idx, gene in enumerate(genes_for_x):
            pct_here = pct_exp.at[celltype, gene] * 100

            if pct_here < dot_min:
                pct_here = np.nan

            plot_df.append({
                "celltype": celltype,
                "gene": gene,
                "x": x_idx,
                "y": y_idx,
                "avg_exp": avg_exp.at[celltype, gene],
                "avg_exp_log": avg_exp_log.at[celltype, gene],
                "avg_exp_scaled": color_expr.at[celltype, gene],
                "pct_exp": pct_here,
            })

    plot_df = pd.DataFrame(plot_df)

    pct_fraction = plot_df["pct_exp"] / 100

    if scale_by == "radius":
        sizes = dot_min_size + (
            pct_fraction ** 2
        ) * (dot_max_size - dot_min_size)
    else:
        sizes = dot_min_size + (
            pct_fraction
        ) * (dot_max_size - dot_min_size)

    plot_df["dot_size"] = sizes

    norm = Normalize(vmin=col_min, vmax=col_max)

    # ------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------
    fig = plt.figure(figsize=(fig_width, fig_height))

    if cluster_celltypes and Z is not None:
        gs = fig.add_gridspec(
            nrows=1,
            ncols=4,
            width_ratios=[dendrogram_width, 1.0, 0.025, 0.14],
            wspace=0.03,
        )

        dend_ax = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1])
        cax = fig.add_subplot(gs[0, 2])
        size_ax = fig.add_subplot(gs[0, 3])

    else:
        gs = fig.add_gridspec(
            nrows=1,
            ncols=3,
            width_ratios=[1.0, 0.025, 0.14],
            wspace=0.08,
        )

        dend_ax = None
        ax = fig.add_subplot(gs[0, 0])
        cax = fig.add_subplot(gs[0, 1])
        size_ax = fig.add_subplot(gs[0, 2])

    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["dot_size"],
        c=plot_df["avg_exp_scaled"],
        cmap=cmap,
        norm=norm,
        edgecolor="black",
        linewidth=0.25,
    )

    # ------------------------------------------------------------
    # Axes
    # ------------------------------------------------------------
    ax.set_xticks(range(len(gene_order_plot)))
    ax.set_xticklabels(
        gene_order_plot,
        rotation=90,
        ha="center",
        va="top",
        fontsize=x_tick_fontsize,
    )

    ax.set_yticks(range(len(celltype_order_present)))
    ax.set_yticklabels(
        celltype_order_present,
        rotation=0,
        ha="right",
        va="center",
        fontsize=y_tick_fontsize,
    )

    ax.set_xlabel("Features", fontsize=axis_label_fontsize)
    ax.set_ylabel("Identity", fontsize=axis_label_fontsize)

    ax.invert_yaxis()

    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", length=0)

    ax.set_xlim(-0.75, len(gene_order_plot) - 0.25)

    top_space = min(top_label_y - 0.50, bracket_y - bracket_height - 0.50)
    ax.set_ylim(len(celltype_order_present) - 0.5, top_space)

    # ------------------------------------------------------------
    # Top marker-panel brackets
    # ------------------------------------------------------------
    for panel, start, end in panel_spans:
        center = (start + end) / 2

        x_left = start - 0.45
        x_right = end + 0.45

        ax.text(
            center,
            top_label_y,
            panel,
            ha="center",
            va="bottom",
            fontsize=top_label_fontsize,
            fontweight="bold",
            rotation=top_label_rotation,
            rotation_mode="anchor",
            clip_on=False,
        )

        ax.plot(
            [x_left, x_right],
            [bracket_y, bracket_y],
            color="black",
            linewidth=bracket_linewidth,
            clip_on=False,
        )

        ax.plot(
            [x_left, x_left],
            [bracket_y, bracket_y + bracket_height],
            color="black",
            linewidth=bracket_linewidth,
            clip_on=False,
        )

        ax.plot(
            [x_right, x_right],
            [bracket_y, bracket_y + bracket_height],
            color="black",
            linewidth=bracket_linewidth,
            clip_on=False,
        )

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    # ------------------------------------------------------------
    # Dendrogram aligned to dotplot rows
    # ------------------------------------------------------------
    if dend_ax is not None:
        dendrogram(
            Z,
            orientation="left",
            ax=dend_ax,
            no_labels=True,
            color_threshold=0,
            above_threshold_color="black",
            link_color_func=lambda k: "black",
        )

        dend_ax.set_ylim(ax.get_ylim())
        dend_ax.axis("off")

    # ------------------------------------------------------------
    # Color legend
    # ------------------------------------------------------------
    cbar = fig.colorbar(
        scatter,
        cax=cax,
    )

    cbar.set_label(
        colorbar_label,
        fontsize=legend_title_fontsize,
    )

    cbar.ax.tick_params(labelsize=legend_fontsize)

    # ------------------------------------------------------------
    # Percent expressed legend
    # ------------------------------------------------------------
    size_ax.axis("off")

    y_positions = np.linspace(0.20, 0.80, len(size_legend_pcts))

    for y, p in zip(y_positions, size_legend_pcts):
        p_fraction = p / 100

        if scale_by == "radius":
            size_here = dot_min_size + (
                p_fraction ** 2
            ) * (dot_max_size - dot_min_size)
        else:
            size_here = dot_min_size + (
                p_fraction
            ) * (dot_max_size - dot_min_size)

        size_ax.scatter(
            0.25,
            y,
            s=size_here,
            facecolor="lightgray",
            edgecolor="black",
            linewidth=0.35,
        )

        size_ax.text(
            0.55,
            y,
            f"{p}%",
            va="center",
            ha="left",
            fontsize=legend_fontsize,
        )

    size_ax.text(
        0.00,
        0.95,
        "Percent Expressed",
        va="top",
        ha="left",
        fontsize=legend_title_fontsize,
    )

    size_ax.set_xlim(0, 1)
    size_ax.set_ylim(0, 1)

    plt.subplots_adjust(
        left=0.12,
        right=0.97,
        top=0.68,
        bottom=0.30,
    )

    if save_path is not None:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")

    plt.show()

    return fig, ax, plot_df

# fig, ax, dotplot_df = dotplot(
#     adata=adata_dot,
#     marker_panels=marker_panels,
#     groupby="cell_type_auto",
#     celltype_order=celltype_order,

#     scale=True,
#     scale_by="radius",
#     col_min=-2.5,
#     col_max=2.5,
#     top_label_rotation=90,
#     top_label_y=-7.5,
#     bracket_y=-1.25,
#     bracket_height=0.25,
    
#     fig_width=26,
#     fig_height=10,
#     dot_min=0,
#     dot_min_size=0,
#     dot_max_size=300,
#     size_legend_pcts=(0, 25, 50, 75, 100),
#     cmap="magma",
#     expression_threshold=0,
#     save_path=None,
# )



# 3. plot UMAP --------------------





def plot_umap(
    adata,
    color_col="cell_type_auto",
    sample_col="sample",
    basis="umap_harmony_rsc",
    mode="split",
    downsample_frac=0.10,
    seed=66,
    ncols=4,
    figsize=None,
    size=1,
    legend_loc=None,
):
    """
    Plot UMAP using coordinates stored as two .obs columns.

    color_col and sample_col may be exact .obs column names
    or unique substrings.

    Robust to duplicated .obs column names (uses the first occurrence).
    """

    # =====================================================
    # Resolve .obs column
    # =====================================================
    def resolve_obs_col(col_name):

        if col_name is None:
            return None

        cols = np.asarray(adata.obs.columns)

        # exact match
        exact = np.where(cols == col_name)[0]

        if len(exact) > 0:
            if len(exact) > 1:
                print(
                    f"Warning: found {len(exact)} duplicated '{col_name}' columns. "
                    "Using the first."
                )
            return col_name

        # substring match
        matches = [c for c in cols if col_name in c]
        matches = list(dict.fromkeys(matches))

        if len(matches) == 0:
            raise KeyError(
                f"No .obs column matching '{col_name}' found.\n"
                f"Available columns:\n{list(cols)}"
            )

        if len(matches) > 1:
            raise ValueError(
                f"Multiple .obs columns match '{col_name}':\n"
                f"{matches}\n"
                "Please use a more specific name."
            )

        return matches[0]

    # =====================================================
    # Resolve UMAP coordinate columns
    # =====================================================
    def resolve_obs_umap_pair(basis):

        cols = list(adata.obs.columns)

        matches_1 = [
            c for c in cols
            if basis in c and c.endswith("_1")
        ]

        matches_2 = [
            c for c in cols
            if basis in c and c.endswith("_2")
        ]

        if len(matches_1) == 0 or len(matches_2) == 0:
            raise KeyError(
                f"Could not find UMAP coordinate pair for basis '{basis}'.\n"
                f"Available UMAP columns:\n"
                f"{[c for c in cols if 'umap' in c.lower()]}"
            )

        if len(matches_1) > 1 or len(matches_2) > 1:
            raise ValueError(
                f"Multiple possible UMAP coordinate columns found.\n"
                f"_1: {matches_1}\n"
                f"_2: {matches_2}"
            )

        return matches_1[0], matches_2[0]

    # =====================================================
    # Helper for duplicated column names
    # =====================================================
    def get_obs_series(col):

        x = adata.obs.loc[:, adata.obs.columns == col]

        if x.shape[1] == 0:
            raise KeyError(f"Column '{col}' not found.")

        if x.shape[1] > 1:
            print(
                f"Warning: using first of "
                f"{x.shape[1]} duplicated '{col}' columns."
            )

        return x.iloc[:, 0]

    color_col = resolve_obs_col(color_col)

    if mode == "split":
        if sample_col is None:
            sample_col = "sample"
        sample_col = resolve_obs_col(sample_col)
    elif sample_col is not None:
        sample_col = resolve_obs_col(sample_col)

    x_col, y_col = resolve_obs_umap_pair(basis)

    basis_key = f"X_{basis}"

    print(f"Using color column : {color_col}")

    if sample_col is not None:
        print(f"Using sample column: {sample_col}")

    print(f"Using UMAP X column: {x_col}")
    print(f"Using UMAP Y column: {y_col}")

    # =====================================================
    # LIGHT PLOTTING OBJECT
    # =====================================================

    obs = pd.DataFrame(index=adata.obs.index)

    obs[color_col] = get_obs_series(color_col)

    if sample_col is not None:
        obs[sample_col] = get_obs_series(sample_col)

    adata_plot = ad.AnnData(
        obs=obs,
        obsm={
            basis_key: adata.obs[[x_col, y_col]].to_numpy()
        },
    )

    adata_plot.obs[color_col] = (
        adata_plot.obs[color_col]
        .astype("category")
    )

    if sample_col is not None:
        adata_plot.obs[sample_col] = (
            adata_plot.obs[sample_col]
            .astype("category")
        )

    # =====================================================
    # STRATIFIED DOWNSAMPLE
    # =====================================================

    if downsample_frac is not None and downsample_frac < 1:

        rng = np.random.default_rng(seed)
        selected_cells = []

        if sample_col is None:

            for _, df_group in adata_plot.obs.groupby(
                color_col,
                observed=True,
            ):

                n_cells = len(df_group)
                n_take = max(
                    1,
                    int(np.floor(n_cells * downsample_frac))
                )

                selected = rng.choice(
                    df_group.index.to_numpy(),
                    size=n_take,
                    replace=False,
                )

                selected_cells.extend(selected)

        else:

            for _, df_sample in adata_plot.obs.groupby(
                sample_col,
                observed=True,
            ):

                for _, df_group in df_sample.groupby(
                    color_col,
                    observed=True,
                ):

                    n_cells = len(df_group)
                    n_take = max(
                        1,
                        int(np.floor(n_cells * downsample_frac))
                    )

                    selected = rng.choice(
                        df_group.index.to_numpy(),
                        size=n_take,
                        replace=False,
                    )

                    selected_cells.extend(selected)

        adata_plot = adata_plot[selected_cells].copy()

        print(f"Original cells   : {adata.n_obs:,}")
        print(f"Downsampled cells: {adata_plot.n_obs:,}")
        print(f"Fraction kept    : {adata_plot.n_obs / adata.n_obs:.3f}")

    else:
        print(f"Cells plotted: {adata_plot.n_obs:,}")

    # =====================================================
    # ALL
    # =====================================================

    if mode == "all":

        if figsize is None:
            figsize = (8, 7)

        fig, ax = plt.subplots(figsize=figsize)

        sc.pl.embedding(
            adata_plot,
            basis=basis_key,
            color=color_col,
            size=size,
            ax=ax,
            show=False,
            title=f"All samples - {color_col}",
            legend_loc=legend_loc,
        )

        plt.tight_layout()
        plt.show()

        return adata_plot

    # =====================================================
    # SPLIT
    # =====================================================

    elif mode == "split":

        samples = adata_plot.obs[sample_col].cat.categories.tolist()

        n_panels = len(samples)
        nrows = int(np.ceil(n_panels / ncols))

        if figsize is None:
            figsize = (5 * ncols, 5 * nrows)

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=figsize,
        )

        axes = np.atleast_1d(axes).flatten()

        for i, sample in enumerate(samples):

            adata_sub = adata_plot[
                adata_plot.obs[sample_col] == sample
            ].copy()

            sc.pl.embedding(
                adata_sub,
                basis=basis_key,
                color=color_col,
                size=size,
                ax=axes[i],
                show=False,
                title=str(sample),
                legend_loc=None,
            )

        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.show()

        return adata_plot

    else:
        raise ValueError(
            "mode must be either 'all' or 'split'"
        )
# # usage
# plot_umap(adata_separated, mode="all",
#           color_col="cell_type_auto",
#           #sample_col="sample",
#           basis="umap_harmony_rsc",
#           legend_loc="right margin")

def set_scanpy_colors(adata, column, category_order, color_map):
    kelly_colors_hex = [
    "#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020",
    "#CEA262", "#817066", "#007D34", "#F6768E", "#00538A",
    "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800",
    "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16",
    ]
    observed = set(
        adata.obs[column]
        .dropna()
        .astype(str)
        .unique()
    )

    category_order = [
        category
        for category in category_order
        if category in observed
    ]

    # Include any unexpected categories at the end
    extra_categories = sorted(
        observed.difference(category_order)
    )

    category_order.extend(extra_categories)

    adata.obs[column] = pd.Categorical(
        adata.obs[column].astype(str),
        categories=category_order,
        ordered=True,
    )

    fallback_colors = iter(kelly_colors_hex)

    colors = []
    for category in category_order:
        if category in color_map:
            colors.append(color_map[category])
        else:
            colors.append(next(fallback_colors, "#817066"))

    adata.uns[f"{column}_colors"] = colors

# # # Usage
# cell_type_order = [
#     "Tumor-Prostate cancer (luminal)",
#     "Tumor-Prostate cancer (secretory luminal)",
#     "Tumor-Prostate cancer (KLK-high)",
#     "Neuroendocrine-like prostate tumor",
#     "Luminal epithelial (club-like)",
#     "Stroma-Fibroblast (CAF)",
#     "Smooth muscle-Myofibroblast",
#     "Endothelial",
#     "Plasma-B-cell",
#     "Antigen-presenting immune cells",
#     "Activated NK-cytotoxic lymphocytes",
#     "Mast cells",
#     "Erythrocytes",
# ]

# cell_type_colors = {
#     # Tumor
#     "Tumor-Prostate cancer (luminal)": "#C10020",
#     "Tumor-Prostate cancer (secretory luminal)": "#E6194B",
#     "Tumor-Prostate cancer (KLK-high)": "#F58231",
#     "Neuroendocrine-like prostate tumor": "#F64E34",

#     # Normal epithelial
#     "Luminal epithelial (club-like)": "#FFE119",

#     # Stroma
#     "Stroma-Fibroblast (CAF)": "#3CB44B",
#     "Smooth muscle-Myofibroblast": "#008080",

#     # Vasculature
#     "Endothelial": "#42D4F4",

#     # Immune
#     "Plasma-B-cell": "#4363D8",
#     "Antigen-presenting immune cells": "#000075",
#     "Activated NK-cytotoxic lymphocytes": "#469990",
#     "Mast cells": "#9A6324",

#     # Other
#     "Erythrocytes": "#808080",
# }

# Utils.set_scanpy_colors(
#     adata=adata,
#     column="celltype_from_leiden_unintegrated_r0.7",
#     category_order=cell_type_order,
#     color_map=cell_type_colors,
# )


# 6.   Spatial Plot  --------------







def _is_continuous_variable(adata, variable):
    """Return True if variable is a gene or numeric obs column."""
    names = [variable] if isinstance(variable, str) else list(variable)
    for name in names:
        if name in adata.var_names:
            continue
        if name in adata.obs.columns and pd.api.types.is_numeric_dtype(
            adata.obs[name]
        ):
            continue
        if name in adata.obs.columns:
            return False
        raise ValueError(f"{name!r} not found in adata.var_names or adata.obs.")
    return True


def SpatialCoord_plot(
    adata,
    variable,
    subset=None,
    subset_values=None,
    x_col="coord_x",
    y_col="coord_y",
    size=8,
    alpha=1.0,
    frameon=False,
    grid=False,
    ncols=4,
    figsize=None,
    plot_order=None,
    shared_legend=None,
    legend_position="bottom",
    panel_size=(5, 5),
    colorbar_loc="right",
    **kwargs,
):
    """
    Plot spatial coordinates stored in adata.obs.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing spatial coordinates and expression data.
    variable : str or list[str]
        Gene(s) in adata.var_names or column(s) in adata.obs to plot.
        Gene expression is always taken from adata.X.
    subset : str, optional
        Name of an adata.obs column used to subset the data before plotting.
    subset_values : str or list[str], optional
        Category or categories within ``subset`` to plot.
        If None, all categories in ``subset`` are plotted individually.
    x_col : str
        Name of x-coordinate column in adata.obs.
    y_col : str
        Name of y-coordinate column in adata.obs.
    size : float
        Point size.
    alpha : float, default=1.0
        Transparency of plotted cells.
    frameon : bool
        Whether to draw axes/frame.
    grid : bool, default=False
        If True, plot all subsets in a single figure.
    ncols : int, default=4
        Number of columns when grid=True.
    figsize : tuple, optional
        Figure size when grid=True. Auto-computed from ``panel_size`` and
        legend height when ``shared_legend`` is enabled.
    plot_order : {None, "rare_first"}, default None
        Controls plotting order for categorical ``obs`` columns.
        - None: keep observation order.
        - "rare_first": plot rare categories on top (ascending abundance).
    shared_legend : bool, optional
        When ``grid=True`` and ``variable`` is a categorical ``obs`` column,
        draw one figure-level legend with all categories. Defaults to True
        in that case.
    legend_position : {"bottom", "right", "none"}, default="bottom"
        Placement of the shared legend when ``shared_legend=True``.
    panel_size : tuple of float, default (5, 5)
        Width and height in inches for each grid panel.
    colorbar_loc : str or None, default "right"
        Placement of the color scale bar for continuous variables (genes or
        numeric ``obs`` columns). Passed to ``scanpy.pl.embedding()``.
        In grid mode, each panel gets its own colorbar. Set to ``None`` to
        hide colorbars.
    **kwargs
        Additional keyword arguments passed to scanpy.pl.embedding().

    Returns
    -------
    Matplotlib figure/axes or Scanpy output.
    """

    if plot_order not in {None, "rare_first"}:
        raise ValueError("plot_order must be None or 'rare_first'.")

    if legend_position not in {"bottom", "right", "none"}:
        raise ValueError(
            "legend_position must be 'bottom', 'right', or 'none'."
        )

    sort_order = False
    color_var = variable if isinstance(variable, str) else None
    is_categorical_obs = (
        color_var is not None
        and color_var in adata.obs.columns
        and not pd.api.types.is_numeric_dtype(adata.obs[color_var])
    )

    if shared_legend is None:
        shared_legend = grid and is_categorical_obs

    def _get_full_categories_and_colors(source_adata):
        if color_var is None or color_var not in source_adata.obs.columns:
            return [], {}

        series = source_adata.obs[color_var]
        if pd.api.types.is_categorical_dtype(series):
            categories = [str(c) for c in series.cat.categories]
        else:
            categories = sorted(series.dropna().astype(str).unique())

        color_key = f"{color_var}_colors"
        color_map = {}
        if color_key in source_adata.uns:
            stored_colors = list(source_adata.uns[color_key])
            color_map = {
                category: stored_colors[i]
                for i, category in enumerate(categories)
                if i < len(stored_colors)
            }

        fallback_colors = plt.cm.tab20.colors
        for i, category in enumerate(categories):
            if category not in color_map:
                color_map[category] = fallback_colors[i % len(fallback_colors)]

        return categories, color_map

    def _apply_category_preservation(ad, full_categories, color_map):
        if not full_categories or color_var is None:
            return ad

        ad = ad.copy()
        ad.obs[color_var] = pd.Categorical(
            ad.obs[color_var].astype(str),
            categories=full_categories,
            ordered=True,
        )
        ad.uns[f"{color_var}_colors"] = [
            color_map[category] for category in full_categories
        ]
        return ad

    def _apply_plot_order(ad):
        if plot_order != "rare_first" or color_var is None:
            return ad

        if color_var not in ad.obs.columns:
            return ad

        counts = ad.obs[color_var].astype(str).value_counts()
        abundance = ad.obs[color_var].astype(str).map(counts)
        order = np.argsort(abundance.to_numpy(), kind="stable")
        return ad[order].copy()

    def _prepare_adata(ad, preserve_categories=None, color_map=None):
        ad = ad.copy()

        ad = _apply_plot_order(ad)

        if preserve_categories:
            ad = _apply_category_preservation(
                ad,
                preserve_categories,
                color_map or {},
            )

        ad.obsm["spatial"] = np.column_stack(
            (
                ad.obs[x_col].to_numpy(),
                ad.obs[y_col].to_numpy(),
            )
        )

        return ad

    # ------------------------------------------------------------------
    # No subsetting
    # ------------------------------------------------------------------
    if subset is None:
        ad = _prepare_adata(adata)

        return sc.pl.embedding(
            ad,
            basis="spatial",
            color=variable,
            size=size,
            alpha=alpha,
            frameon=frameon,
            sort_order=sort_order,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # Validate subset
    # ------------------------------------------------------------------
    if subset not in adata.obs.columns:
        raise ValueError(f"{subset!r} not found in adata.obs.")

    if subset_values is None:
        subset_values = pd.unique(adata.obs[subset])

    if isinstance(subset_values, str):
        subset_values = [subset_values]

    # ------------------------------------------------------------------
    # One figure per subset
    # ------------------------------------------------------------------
    if not grid:
        plots = []

        for value in subset_values:
            ad = _prepare_adata(adata[adata.obs[subset] == value])

            plots.append(
                sc.pl.embedding(
                    ad,
                    basis="spatial",
                    color=variable,
                    title=str(value),
                    size=size,
                    alpha=alpha,
                    frameon=frameon,
                    sort_order=sort_order,
                    **kwargs,
                )
            )

        return plots

    # ------------------------------------------------------------------
    # Grid layout
    # ------------------------------------------------------------------
    nplots = len(subset_values)
    ncols = min(ncols, nplots)
    nrows = math.ceil(nplots / ncols)

    full_categories, color_map = _get_full_categories_and_colors(adata)
    use_shared_legend = shared_legend and bool(full_categories)
    is_continuous = _is_continuous_variable(adata, variable)

    panel_w, panel_h = panel_size
    legend_ncol = min(len(full_categories), 5) if full_categories else 1
    legend_rows = (
        math.ceil(len(full_categories) / legend_ncol)
        if use_shared_legend and legend_position != "none"
        else 0
    )
    legend_h = 1.2 * legend_rows

    if figsize is None:
        total_w = panel_w * ncols
        if use_shared_legend and legend_position == "right":
            total_w += max(2.0, 0.35 * len(full_categories))
        elif is_continuous and colorbar_loc is not None:
            total_w += 0.6 * ncols
        figsize = (total_w, panel_h * nrows + legend_h)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize,
        squeeze=False,
    )

    axes = axes.ravel()

    for i, value in enumerate(subset_values):
        ad = _prepare_adata(
            adata[adata.obs[subset] == value],
            preserve_categories=full_categories if use_shared_legend else None,
            color_map=color_map,
        )

        panel_legend = None
        panel_colorbar = None
        if is_continuous:
            panel_colorbar = colorbar_loc
        elif not use_shared_legend or legend_position == "none":
            if i == nplots - 1:
                panel_legend = "right margin"
        else:
            panel_legend = None
            panel_colorbar = None

        sc.pl.embedding(
            ad,
            basis="spatial",
            color=variable,
            ax=axes[i],
            show=False,
            title=str(value),
            size=size,
            alpha=alpha,
            frameon=frameon,
            sort_order=sort_order,
            legend_loc=panel_legend,
            colorbar_loc=panel_colorbar,
            **kwargs,
        )

    # Remove empty axes
    for ax in axes[nplots:]:
        fig.delaxes(ax)

    if use_shared_legend and legend_position != "none":
        legend_handles = [
            plt.Rectangle((0, 0), 1, 1, color=color_map[category])
            for category in full_categories
        ]

        if legend_position == "bottom":
            fig.legend(
                legend_handles,
                full_categories,
                title=str(color_var),
                loc="upper center",
                bbox_to_anchor=(0.5, -0.02),
                ncol=legend_ncol,
                frameon=False,
            )
            bottom_margin = 0.06 + 0.035 * legend_rows
            fig.tight_layout(rect=[0, bottom_margin, 1, 1])
        else:
            fig.legend(
                legend_handles,
                full_categories,
                title=str(color_var),
                loc="center left",
                bbox_to_anchor=(1.01, 0.5),
                frameon=False,
            )
            fig.tight_layout(rect=[0, 0, 0.82, 1])
    elif is_continuous and colorbar_loc is not None:
        fig.tight_layout()
        fig.subplots_adjust(wspace=0.45)
    else:
        plt.tight_layout()

    return fig

# # # Usage
# # -----------------------------------------------------------------------------
# # Example 1: Plot a categorical annotation
# # -----------------------------------------------------------------------------
# SpatialCoord_plot(
#     adata,
#     variable="celltype_from_leiden_unintegrated_r0.7",
#     size=4,                  # Point size
#     alpha=0.8,               # Transparency
#     plot_order="rare_first",  # None or "rare_first"
# )

# # -----------------------------------------------------------------------------
# # Example 2: Plot a gene across all samples in a grid
# # -----------------------------------------------------------------------------
# SpatialCoord_plot(
#     adata,
#     variable="EPCAM",
#     subset="sample",
#     grid=True,
#     ncols=4,
#     cmap="viridis",
#     vmin="p1",
#     vmax="p99",
#     size=2,
#     alpha=0.6,
#     colorbar_loc="right",  # one colorbar per panel (default)
# )

# # -----------------------------------------------------------------------------
# # Example 3: Plot selected samples only
# # -----------------------------------------------------------------------------
# SpatialCoord_plot(
#     adata,
#     variable="celltype_from_leiden_unintegrated_r0.7",
#     subset="sample",
#     subset_values=["Sample_1", "Sample_2", "Sample_3"],
#     grid=True,
#     ncols=3,
#     size=3,
#     alpha=0.5,
#     plot_order="rare_first",
# )
