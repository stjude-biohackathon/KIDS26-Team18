"""Plotting helpers for BioHack spatial notebook templates."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd

from lib.qc_common import ensure_output_dir

CELL_TYPE_COLORS = {
    "Tumor state A": "#d62728",
    "Tumor state B": "#ff9896",
    "CD8 T cells": "#1f77b4",
    "CD4 T cells": "#aec7e8",
    "Macrophages": "#2ca02c",
    "Fibroblasts": "#9467bd",
    "Endothelial cells": "#8c564b",
    "Unassigned": "#bdbdbd",
}


def _color_for_label(label: str, colors: dict[str, str]) -> str:
    return colors.get(str(label), "#cccccc")


def plot_spatial_core(
    adata: ad.AnnData,
    core_id: str,
    ax: plt.Axes | None = None,
    color_col: str = "cell_type",
    cell_type_colors: dict[str, str] | None = None,
    title: str | None = None,
) -> plt.Axes:
    subset = adata[adata.obs["core_id"] == core_id]
    colors = cell_type_colors or CELL_TYPE_COLORS
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    for cell_type in subset.obs[color_col].astype(str).unique():
        mask = subset.obs[color_col].astype(str) == cell_type
        color = _color_for_label(cell_type, colors)
        if mask.any():
            ax.scatter(
                subset.obsm["spatial"][mask, 0],
                subset.obsm["spatial"][mask, 1],
                s=12,
                c=color,
                label=cell_type,
                alpha=0.85,
                linewidths=0,
            )

    ax.set_title(title or f"Synthetic TMA core: {core_id}")
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
    ax.set_aspect("equal")
    return ax


def plot_composition(
    adata: ad.AnnData,
    column: str,
    core_col: str = "core_id",
    output_path: Path | None = None,
) -> plt.Figure:
    composition = (
        adata.obs.groupby([core_col, column], observed=True)
        .size()
        .unstack(fill_value=0)
        .div(adata.obs.groupby(core_col).size(), axis=0)
    )

    fig, ax = plt.subplots(figsize=(8, 4))
    composition.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel("TMA core")
    ax.set_title(f"Cell composition by {column}")
    ax.legend(title=column, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()

    if output_path is not None:
        ensure_output_dir(str(output_path))
        fig.savefig(output_path, bbox_inches="tight")
    return fig


def plot_qc_bar(
    counts: pd.Series,
    title: str,
    xlabel: str,
    ylabel: str,
    output_path: Path | None = None,
    color: str = "#4c72b0",
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(6, 4))
    counts.plot(kind="bar", ax=ax, color=color)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    if output_path is not None:
        ensure_output_dir(str(output_path))
        fig.savefig(output_path, bbox_inches="tight")
    return fig


def plot_confusion_matrix(
    crosstab: pd.DataFrame,
    title: str,
    output_path: Path | None = None,
    figsize: tuple[float, float] = (7, 6),
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(crosstab.to_numpy(), aspect="auto", cmap="Blues")
    ax.set_xticks(range(crosstab.shape[1]))
    ax.set_yticks(range(crosstab.shape[0]))
    ax.set_xticklabels(crosstab.columns.astype(str), rotation=45, ha="right")
    ax.set_yticklabels(crosstab.index.astype(str))
    ax.set_xlabel(crosstab.columns.name or "Predicted")
    ax.set_ylabel(crosstab.index.name or "Reference")
    ax.set_title(title)

    for i in range(crosstab.shape[0]):
        for j in range(crosstab.shape[1]):
            value = int(crosstab.iloc[i, j])
            ax.text(j, i, str(value), ha="center", va="center", color="black", fontsize=8)

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    if output_path is not None:
        ensure_output_dir(str(output_path))
        fig.savefig(output_path, bbox_inches="tight")
    return fig


def plot_agreement_bars(
    summary_df: pd.DataFrame,
    output_path: Path | None = None,
    metric: str = "accuracy",
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 4))
    x = summary_df["comparison"].astype(str)
    y = summary_df[metric].astype(float)
    bars = ax.bar(x, y, color=["#4c72b0", "#55a868", "#c44e52"])
    ax.set_ylim(0, 1)
    ax.set_ylabel(metric.replace("_", " ").title())
    ax.set_title(f"Cell typing agreement ({metric.replace('_', ' ')})")
    ax.tick_params(axis="x", rotation=20)
    for bar, value in zip(bars, y, strict=False):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.02, f"{value:.2f}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    if output_path is not None:
        ensure_output_dir(str(output_path))
        fig.savefig(output_path, bbox_inches="tight")
    return fig


def plot_spatial_celltyping_comparison(
    adata: ad.AnnData,
    core_id: str,
    output_path: Path | None = None,
    truth_key: str = "cell_type",
    rule_key: str = "cell_type_rule",
    scanpy_key: str = "cell_type_scanpy",
) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    specs = [
        (truth_key, "Ground truth"),
        (rule_key, "Rule-based"),
        (scanpy_key, "Scanpy clustering"),
    ]
    for ax, (col, label) in zip(axes, specs, strict=True):
        plot_spatial_core(adata, core_id=core_id, ax=ax, color_col=col, title=f"{label}: {core_id}")
    fig.suptitle("Cell typing comparison across annotation modes", y=1.02)
    fig.tight_layout()
    if output_path is not None:
        ensure_output_dir(str(output_path))
        fig.savefig(output_path, bbox_inches="tight")
    return fig
