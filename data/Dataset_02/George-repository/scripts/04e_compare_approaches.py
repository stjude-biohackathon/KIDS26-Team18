#!/usr/bin/env python
"""Compare FOV46 Approach I, Approach II, and published/reference annotation.

The comparison is categorical. We therefore use contingency tables and
row-normalized concordance tables (not Pearson/Spearman correlation matrices).
Each comparison is visualized as a heatmap and a 100% stacked bar plot.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from src.plotting import (
    plot_umap_annotation,
    plot_umap_side_by_side,
    plot_spatial_annotation_pair,
)

IN = ROOT / "results/FOV46/GSM9046088_FOV46_approach_i_ii.h5ad"
OUT = ROOT / "results/FOV46"
FIG = OUT / "figures/comparison"
TAB = OUT / "tables"
SDATA = ROOT / "data/spatial/GSM9046088_CosMx_FOV46.zarr"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

APPROACH_I = "cell_type_scanpy"
APPROACH_II = "cell_type_rule"
REFERENCE = "final_CT"  # article/metadata annotation; treated as reference, not absolute ground truth

print("=" * 70)
print("FOV46 CELL-TYPE COMPARISON")
print("=" * 70)
print("Input:", IN)
a = sc.read_h5ad(IN)

for col in (APPROACH_I, APPROACH_II):
    if col not in a.obs:
        raise KeyError(f"Required annotation column {col!r} not found in adata.obs")


def contingency_tables(data, row_col, col_col):
    """Return raw counts and P(column label | row label)."""
    counts = pd.crosstab(data.obs[row_col], data.obs[col_col], dropna=False)
    row_fraction = pd.crosstab(
        data.obs[row_col], data.obs[col_col], normalize="index", dropna=False
    )
    return counts, row_fraction


def save_heatmap(table, output, title, xlabel, ylabel):
    if table.empty:
        return
    width = max(8, 0.9 * table.shape[1] + 4)
    height = max(6, 0.55 * table.shape[0] + 3)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        table,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        vmin=0,
        vmax=1,
        linewidths=0.25,
        ax=ax,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=250, bbox_inches="tight")
    plt.close(fig)


def save_stacked_bar(row_fraction, output, title, xlabel):
    """Plot the same row-normalized contingency data as a 100% stacked bar."""
    if row_fraction.empty:
        return
    ax = row_fraction.plot(
        kind="bar",
        stacked=True,
        figsize=(max(10, 0.8 * row_fraction.shape[0] + 5), 7),
        width=0.85,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fraction of cells")
    ax.set_ylim(0, 1)
    ax.set_title(title)
    ax.legend(
        title="Assigned cell type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    plt.tight_layout()
    plt.savefig(output, dpi=250, bbox_inches="tight")
    plt.close()


def run_comparison(row_col, col_col, prefix, row_label, col_label):
    counts, row_fraction = contingency_tables(a, row_col, col_col)
    counts.to_csv(TAB / f"{prefix}_counts.csv")
    row_fraction.to_csv(TAB / f"{prefix}_row_fraction.csv")

    save_heatmap(
        row_fraction,
        FIG / f"{prefix}_concordance_heatmap.png",
        f"{row_label} vs {col_label}",
        col_label,
        row_label,
    )
    save_stacked_bar(
        row_fraction,
        FIG / f"{prefix}_stacked_bar.png",
        f"{col_label} labels within {row_label} groups",
        row_label,
    )
    return counts, row_fraction


# -----------------------------------------------------------------------------
# Approach I vs Approach II
# -----------------------------------------------------------------------------
ct_i_ii, frac_i_ii = run_comparison(
    APPROACH_I,
    APPROACH_II,
    "approach_i_vs_ii",
    "Approach I",
    "Approach II",
)

# Exact string agreement is informative only where both methods use the same
# cell-type vocabulary. Keep the per-cell flag and report it transparently.
a.obs["typing_agreement"] = (
    a.obs[APPROACH_I].astype(str) == a.obs[APPROACH_II].astype(str)
)

summary = {
    "n_cells": int(a.n_obs),
    "approach_i_vs_ii_exact_agreement_fraction": float(a.obs["typing_agreement"].mean()),
    "approach_ii_unassigned_fraction": float(
        a.obs[APPROACH_II].astype(str).eq("Unassigned").mean()
    ),
    "approach_ii_ambiguous_fraction": float(
        a.obs[APPROACH_II].astype(str).eq("Ambiguous").mean()
    ),
}

# Unified UMAP comparison and two spatial views for both approaches.
plot_umap_side_by_side(
    a, APPROACH_I, APPROACH_II, FIG / "umap_approach_i_vs_ii.png",
)
for col, prefix, title in [
    (APPROACH_I, "approach_i", "Approach I — cluster-based"),
    (APPROACH_II, "approach_ii", "Approach II — rule-based"),
    ("typing_agreement", "approach_i_vs_ii_agreement", "Approach I vs II — exact agreement"),
]:
    plot_spatial_annotation_pair(
        a, col, SDATA, FIG / "spatial", prefix=prefix, title=title, point_size=8,
    )

# -----------------------------------------------------------------------------
# Published/metadata reference annotation vs both approaches
# -----------------------------------------------------------------------------
if REFERENCE in a.obs:
    print(f"Reference annotation found: {REFERENCE}")

    run_comparison(
        REFERENCE,
        APPROACH_I,
        "reference_vs_approach_i",
        f"Reference ({REFERENCE})",
        "Approach I",
    )
    run_comparison(
        REFERENCE,
        APPROACH_II,
        "reference_vs_approach_ii",
        f"Reference ({REFERENCE})",
        "Approach II",
    )

    # Save reference label frequencies, UMAP, and both spatial views.
    a.obs[REFERENCE].value_counts(dropna=False).to_csv(TAB / "reference_label_counts.csv")
    plot_umap_annotation(
        a, REFERENCE, FIG / f"umap_{REFERENCE}.png",
        title=f"Reference annotation — {REFERENCE}",
    )
    plot_spatial_annotation_pair(
        a, REFERENCE, SDATA, FIG / "spatial",
        prefix="reference", title=f"Reference annotation — {REFERENCE}", point_size=8,
    )

    # Raw exact-string agreement is reported separately because different naming
    # granularity (e.g. Myeloid vs Macrophage, Epithelial vs AT2) can make this
    # statistic underestimate biological concordance.
    ref_valid = a.obs[REFERENCE].notna()
    if ref_valid.any():
        summary["reference_vs_approach_i_exact_string_agreement_fraction"] = float(
            (
                a.obs.loc[ref_valid, REFERENCE].astype(str)
                == a.obs.loc[ref_valid, APPROACH_I].astype(str)
            ).mean()
        )
        summary["reference_vs_approach_ii_exact_string_agreement_fraction"] = float(
            (
                a.obs.loc[ref_valid, REFERENCE].astype(str)
                == a.obs.loc[ref_valid, APPROACH_II].astype(str)
            ).mean()
        )
else:
    print(
        f"Reference column {REFERENCE!r} is not present in this object; "
        "reference comparisons are skipped."
    )

pd.Series(summary, name="value").to_csv(TAB / "typing_summary.csv")

OUT_H5AD = OUT / "GSM9046088_FOV46_typing_compared.h5ad"
a.write_h5ad(OUT_H5AD, compression="gzip")

print("=" * 70)
print("Comparison complete")
print(pd.Series(summary))
print("Saved:", OUT_H5AD)
print("Figures:", FIG)
print("Tables:", TAB)
