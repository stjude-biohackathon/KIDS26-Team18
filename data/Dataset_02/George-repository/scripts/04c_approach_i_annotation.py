#!/usr/bin/env python
"""Approach I: manual cell-type annotation after selecting a Leiden resolution.

Inputs
------
results/FOV46/GSM9046088_FOV46_leiden_sweep.h5ad
config/fov46_annotation.json

Outputs include full/top marker tables, marker visualizations, annotated UMAP/spatial
plots, and GSM9046088_FOV46_approach_i.h5ad.
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from config.markers import MARKER_PANEL
from src.plotting import (
    plot_umap_annotation,
    plot_spatial_annotation_pair,
)

IN = ROOT / "results/FOV46/GSM9046088_FOV46_leiden_sweep.h5ad"
OUT = ROOT / "results/FOV46"
FIG = OUT / "figures/approach_i"
TAB = OUT / "tables"
CFG = ROOT / "config/fov46_annotation.json"
SDATA = ROOT / "data/spatial/GSM9046088_CosMx_FOV46.zarr"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

TOP_N_TABLE = 10
TOP_N_HEATMAP = 5

print("=" * 70)
print("FOV46 APPROACH I: SELECTED LEIDEN + MANUAL ANNOTATION")
print("=" * 70)
print("Input:", IN)
print("Config:", CFG)

a = sc.read_h5ad(IN)
cfg = json.loads(CFG.read_text())
r = float(cfg["resolution"])
key = f"leiden_r{r:g}"

if key not in a.obs:
    available = [c for c in a.obs.columns if c.startswith("leiden_r")]
    raise KeyError(f"{key!r} not found. Available Leiden columns: {available}")

# Use the manually selected resolution; do not re-run graph construction or Leiden.
a.obs["leiden"] = a.obs[key].copy().astype("category")
a.uns["selected_leiden_resolution"] = r
print(f"Selected resolution: {r:g}")
print(a.obs["leiden"].value_counts().sort_index())

# -----------------------------------------------------------------------------
# Differential markers: each selected Leiden cluster vs. the rest
# -----------------------------------------------------------------------------
rank_key = "rank_genes_groups_leiden_selected"
sc.tl.rank_genes_groups(
    a,
    groupby="leiden",
    method="wilcoxon",
    use_raw=True,
    pts=True,
    key_added=rank_key,
)

markers = sc.get.rank_genes_groups_df(a, group=None, key=rank_key)
markers.to_csv(TAB / f"leiden_markers_r{r:g}_all.csv", index=False)

# Top N genes per cluster, ordered by Scanpy score.
top10 = (
    markers.sort_values(["group", "scores"], ascending=[True, False])
    .groupby("group", observed=True, sort=False)
    .head(TOP_N_TABLE)
    .copy()
)
top10.to_csv(TAB / f"leiden_top{TOP_N_TABLE}_markers_r{r:g}.csv", index=False)

# A compact gene x cluster table is useful for manual inspection/GitHub outputs.
top_gene_table = (
    top10.assign(rank=top10.groupby("group", observed=True).cumcount() + 1)
    .pivot(index="rank", columns="group", values="names")
)
top_gene_table.to_csv(TAB / f"leiden_top{TOP_N_TABLE}_marker_names_r{r:g}.csv")

# Original-style "cluster vs rest" top-gene panels.
sc.pl.rank_genes_groups(
    a,
    key=rank_key,
    n_genes=TOP_N_TABLE,
    sharey=False,
    show=False,
)
plt.gcf().savefig(
    FIG / f"leiden_top{TOP_N_TABLE}_markers_r{r:g}.png",
    dpi=250,
    bbox_inches="tight",
)
plt.close("all")

# Top DE genes across all clusters: dotplot and heatmap.
top_genes = list(dict.fromkeys(top10["names"].dropna().astype(str).tolist()))
if top_genes:
    dp = sc.pl.dotplot(
        a,
        var_names=top_genes,
        groupby="leiden",
        use_raw=True,
        standard_scale="var",
        dendrogram=True,
        show=False,
        return_fig=True,
    )
    dp.savefig(FIG / f"leiden_top_markers_dotplot_r{r:g}.png", dpi=250)
    plt.close("all")

    sc.pl.rank_genes_groups_heatmap(
        a,
        key=rank_key,
        groupby="leiden",
        n_genes=TOP_N_HEATMAP,
        use_raw=True,
        show_gene_labels=True,
        show=False,
    )
    plt.gcf().savefig(
        FIG / f"leiden_top{TOP_N_HEATMAP}_marker_heatmap_r{r:g}.png",
        dpi=250,
        bbox_inches="tight",
    )
    plt.close("all")

# Known marker-panel dotplot for biological interpretation.
raw_var_names = a.raw.var_names if a.raw is not None else a.var_names
present = {
    cell_type: [g for g in genes if g in raw_var_names]
    for cell_type, genes in MARKER_PANEL.items()
}
present = {cell_type: genes for cell_type, genes in present.items() if genes}
if present:
    dp = sc.pl.dotplot(
        a,
        var_names=present,
        groupby="leiden",
        use_raw=(a.raw is not None),
        standard_scale="var",
        dendrogram=True,
        show=False,
        return_fig=True,
    )
    dp.savefig(FIG / f"marker_panel_dotplot_selected_r{r:g}.png", dpi=250)
    plt.close("all")

# Selected-resolution UMAP and two spatial views before annotation.
plot_umap_annotation(
    a, "leiden", FIG / f"umap_selected_r{r:g}.png",
    title=f"Leiden r={r:g}", legend_loc="on data",
)
plot_spatial_annotation_pair(
    a, "leiden", SDATA, FIG / "spatial",
    prefix=f"leiden_r{r:g}", title=f"Leiden r={r:g}", point_size=8,
)

# -----------------------------------------------------------------------------
# Manual cluster -> cell-type mapping
# -----------------------------------------------------------------------------
ann = {str(k): str(v) for k, v in cfg.get("cluster_annotations", {}).items()}
clusters = sorted(a.obs["leiden"].astype(str).unique())
missing = [c for c in clusters if c not in ann or not ann[c].strip() or ann[c] == "EDIT_ME"]
extra = sorted(set(ann) - set(clusters))

if extra:
    print("WARNING: annotation config contains cluster IDs not present here:", extra)
if missing:
    raise ValueError(
        "Missing manual annotations for clusters "
        f"{missing}. Review top markers, marker-panel dotplot, UMAP and spatial plots; "
        f"then edit {CFG} and rerun."
    )

a.obs["cell_type_scanpy"] = (
    a.obs["leiden"].astype(str).map(ann).astype("category")
)

annotation_table = pd.DataFrame(
    {"leiden": clusters, "cell_type_scanpy": [ann[c] for c in clusters]}
)
annotation_table.to_csv(TAB / f"approach_i_cluster_annotations_r{r:g}.csv", index=False)
a.obs["cell_type_scanpy"].value_counts().to_csv(TAB / "approach_i_label_counts.csv")

# Annotated UMAP plus coordinate and segmentation spatial views.
plot_umap_annotation(
    a, "cell_type_scanpy", FIG / "umap_cell_type_scanpy.png",
    title="Approach I — cluster-based cell typing",
)
plot_spatial_annotation_pair(
    a, "cell_type_scanpy", SDATA, FIG / "spatial",
    prefix="approach_i", title="Approach I — cluster-based", point_size=8,
)

OUT_H5AD = OUT / "GSM9046088_FOV46_approach_i.h5ad"
a.write_h5ad(OUT_H5AD, compression="gzip")

print("=" * 70)
print("Approach I complete")
print("Resolution:", r)
print(a.obs["cell_type_scanpy"].value_counts())
print("Saved:", OUT_H5AD)
print("Figures:", FIG)
print("Tables:", TAB)
