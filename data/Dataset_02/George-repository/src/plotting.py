"""Shared plotting helpers for the Dataset_02 CosMx FOV46 workflow.

The goal is to keep Approach I, Approach II, and reference-annotation figures
visually comparable.  UMAPs keep the Scanpy frame/UMAP1/UMAP2 axes. Spatial
annotations are rendered in two complementary ways:
  1. cell-center coordinate scatter from AnnData;
  2. cell-segmentation polygons from SpatialData.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

FIG_DPI = 250
UMAP_FIGSIZE = (8, 6.5)
SPATIAL_FIGSIZE = (8, 8)
POINT_SIZE = 12


def _ensure_parent(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _categorical_palette(adata, color):
    """Return Scanpy's stored category->color mapping when available."""
    if color not in adata.obs:
        return None
    s = adata.obs[color]
    if not isinstance(s.dtype, pd.CategoricalDtype):
        return None
    key = f"{color}_colors"
    if key not in adata.uns or len(adata.uns[key]) != len(s.cat.categories):
        # Ask Scanpy to establish its normal palette without changing theme.
        sc.pl.umap(adata, color=color, frameon=True, show=False)
        plt.close("all")
    colors = adata.uns.get(key)
    if colors is None:
        return None
    return dict(zip(map(str, s.cat.categories), colors))


def plot_umap_annotation(
    adata,
    color,
    output,
    *,
    title=None,
    legend_loc="right margin",
    size=None,
    ax=None,
    show=False,
):
    """Unified UMAP theme: visible frame and UMAP1/UMAP2 axes."""
    output = _ensure_parent(output) if output is not None else None
    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=UMAP_FIGSIZE)
    else:
        fig = ax.figure

    sc.pl.umap(
        adata,
        color=color,
        ax=ax,
        show=False,
        frameon=True,
        legend_loc=legend_loc,
        size=size,
        title=title,
    )
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")

    if own_fig and output is not None:
        fig.tight_layout()
        fig.savefig(output, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
    elif show:
        plt.show()
    return ax


def plot_umap_side_by_side(
    adata,
    left_color,
    right_color,
    output,
    *,
    left_title="Approach I — cluster-based",
    right_title="Approach II — rule-based",
):
    output = _ensure_parent(output)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    plot_umap_annotation(adata, left_color, None, title=left_title, ax=axes[0])
    plot_umap_annotation(adata, right_color, None, title=right_title, ax=axes[1])
    fig.tight_layout()
    fig.savefig(output, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    return fig


def _spatial_xy(adata):
    if "spatial" in adata.obsm:
        xy = np.asarray(adata.obsm["spatial"])
        if xy.ndim == 2 and xy.shape[1] >= 2:
            return xy[:, 0], xy[:, 1]
    for xcol, ycol in [
        ("CenterX_global_px", "CenterY_global_px"),
        ("x_global_px", "y_global_px"),
        ("coord_x", "coord_y"),
    ]:
        if xcol in adata.obs and ycol in adata.obs:
            return adata.obs[xcol].to_numpy(), adata.obs[ycol].to_numpy()
    raise KeyError("No spatial coordinates found in obsm['spatial'] or known obs columns.")


def spatial_scatter(
    adata,
    color,
    output,
    *,
    title=None,
    size=POINT_SIZE,
    alpha=0.85,
    invert_y=True,
):
    """Cell-center coordinate plot with the same categorical colors as Scanpy UMAP."""
    output = _ensure_parent(output)
    x, y = _spatial_xy(adata)
    values = adata.obs[color]
    fig, ax = plt.subplots(figsize=SPATIAL_FIGSIZE)

    if isinstance(values.dtype, pd.CategoricalDtype) or values.dtype == object:
        if not isinstance(values.dtype, pd.CategoricalDtype):
            values = values.astype("category")
        palette = _categorical_palette(adata, color)
        if palette is None:
            cats = list(map(str, values.cat.categories))
            cmap = plt.get_cmap("tab20", max(len(cats), 1))
            palette = {cat: cmap(i) for i, cat in enumerate(cats)}
        for cat in values.cat.categories:
            mask = values.astype(str).to_numpy() == str(cat)
            ax.scatter(x[mask], y[mask], s=size, alpha=alpha, label=str(cat),
                       color=palette.get(str(cat)))
        ax.legend(title=color, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    else:
        pts = ax.scatter(x, y, c=np.asarray(values), s=size, alpha=alpha)
        fig.colorbar(pts, ax=ax, label=color)

    ax.set_xlabel("Global X")
    ax.set_ylabel("Global Y")
    ax.set_title(title or color)
    ax.set_aspect("equal")
    if invert_y:
        ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(output, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    return fig


def spatial_segmentation(
    adata,
    color,
    sdata_path,
    output,
    *,
    title=None,
    shapes_key="cell_boundaries",
    table_key="table",
    coordinate_system="global",
):
    """Color SpatialData cell-segmentation polygons using an AnnData obs annotation.

    Annotation is aligned by cell IDs (obs_names), so the plotting SpatialData table
    receives exactly the same labels used for the UMAP/coordinate plot.
    """
    import spatialdata as sd
    import spatialdata_plot  # noqa: F401  # registers .pl

    output = _ensure_parent(output)
    sdata = sd.read_zarr(Path(sdata_path))
    if table_key not in sdata.tables:
        raise KeyError(f"{table_key!r} not found in SpatialData tables")
    if shapes_key not in sdata.shapes:
        raise KeyError(f"{shapes_key!r} not found in SpatialData shapes")

    table = sdata[table_key]
    common = table.obs_names.intersection(adata.obs_names)
    if len(common) == 0:
        raise ValueError("No matching cell IDs between analysis AnnData and SpatialData table.")

    # Preserve cells outside the analyzed/QC-passed object as NA; analyzed cells get labels.
    aligned = pd.Series(pd.NA, index=table.obs_names, dtype="object")
    aligned.loc[common] = adata.obs.loc[common, color].astype(str).to_numpy()
    categories = list(pd.unique(adata.obs[color].astype(str)))
    table.obs[color] = pd.Categorical(aligned, categories=categories)

    # Reuse Scanpy's category colors where possible.
    palette = _categorical_palette(adata, color)
    if palette is not None:
        table.uns[f"{color}_colors"] = [palette[str(c)] for c in table.obs[color].cat.categories]

    # Current project uses spatialdata-plot 0.3.3; element= must be explicit.
    sdata.pl.render_shapes(
        element=shapes_key,
        color=color,
        fill_alpha=0.9,
        outline_alpha=0.35,
        outline_color="black",
    ).pl.show(
        coordinate_systems=coordinate_system,
        figsize=SPATIAL_FIGSIZE,
        title=title or f"{color} — cell segmentation",
    )
    fig = plt.gcf()
    fig.savefig(output, dpi=FIG_DPI, bbox_inches="tight")
    plt.close("all")
    return fig


def plot_spatial_annotation_pair(
    adata,
    color,
    sdata_path,
    output_dir,
    *,
    prefix=None,
    title=None,
    point_size=POINT_SIZE,
):
    """Write coordinate, filled-segmentation, and outline-segmentation views."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = prefix or color
    spatial_scatter(
        adata,
        color,
        output_dir / f"{prefix}_coordinates.png",
        title=(title or color) + " — cell centers",
        size=point_size,
    )
    spatial_segmentation(
        adata,
        color,
        sdata_path,
        output_dir / f"{prefix}_cell_segments.png",
        title=(title or color) + " — cell segmentation",
    )
    spatial_segmentation_outline(
        adata,
        color,
        sdata_path,
        output_dir / f"{prefix}_cell_segment_outlines.png",
        title=(title or color) + " — cell segmentation outlines",
    )


def spatial_segmentation_outline(
    adata,
    color,
    sdata_path,
    output,
    *,
    title=None,
    shapes_key="cell_boundaries",
    coordinate_system="global",
    center_size=8,
    center_alpha=0.9,
    boundary_linewidth=0.35,
):
    """Plot transparent cell polygons as white outlines plus annotation-colored centers.

    This third spatial view emphasizes the actual segmentation geometry without filling
    cell interiors. A black background and white polygon boundaries show morphology;
    small center points retain the cell-type annotation colors used in UMAP and the
    coordinate-scatter view.
    """
    import spatialdata as sd

    output = _ensure_parent(output)
    sdata = sd.read_zarr(Path(sdata_path))
    if shapes_key not in sdata.shapes:
        raise KeyError(f"{shapes_key!r} not found in SpatialData shapes")

    shapes = sdata[shapes_key]
    common = shapes.index.intersection(adata.obs_names)
    if len(common) == 0:
        raise ValueError("No matching cell IDs between analysis AnnData and SpatialData shapes.")

    # Draw only cells present in the analyzed AnnData (e.g. QC-passed cells).
    shapes_plot = shapes.loc[common]
    values = adata.obs.loc[common, color]

    if not isinstance(values.dtype, pd.CategoricalDtype):
        values = values.astype("category")
    palette = _categorical_palette(adata, color)
    if palette is None:
        cats = list(map(str, values.cat.categories))
        cmap = plt.get_cmap("tab20", max(len(cats), 1))
        palette = {cat: cmap(i) for i, cat in enumerate(cats)}

    fig, ax = plt.subplots(figsize=SPATIAL_FIGSIZE)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    # GeoPandas draws transparent polygon interiors and white cell boundaries.
    shapes_plot.boundary.plot(
        ax=ax,
        color="white",
        linewidth=boundary_linewidth,
        alpha=0.9,
    )

    # Use the same global cell-center coordinates as the coordinate view.
    x, y = _spatial_xy(adata)
    pos = pd.DataFrame({"x": x, "y": y}, index=adata.obs_names).loc[common]
    value_strings = values.astype(str)
    for cat in values.cat.categories:
        mask = value_strings.to_numpy() == str(cat)
        ax.scatter(
            pos.loc[mask, "x"],
            pos.loc[mask, "y"],
            s=center_size,
            alpha=center_alpha,
            color=palette.get(str(cat)),
            label=str(cat),
            linewidths=0,
        )

    ax.set_xlabel("Global X", color="white")
    ax.set_ylabel("Global Y", color="white")
    ax.set_title(title or f"{color} — segmentation outlines", color="white")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("white")
    legend = ax.legend(
        title=color,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    if legend is not None:
        plt.setp(legend.get_texts(), color="white")
        plt.setp(legend.get_title(), color="white")

    fig.tight_layout()
    fig.savefig(output, dpi=FIG_DPI, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    return fig
