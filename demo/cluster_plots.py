"""Notebook-style cluster visualizations (UMAP + rank-genes dotplot)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from matplotlib.figure import Figure

import anndata as ad

# From 01_celltyping_compare.ipynb — manual Leiden → Scanpy label mapping
CLUSTER_ANNOTATIONS: dict[str, str] = {
    "0": "Epithelial",
    "1": "Myofibroblast",
    "2": "Fibroblast",
    "3": "Endothelial",
    "4": "Plasma cell",
    "5": "T cell",
    "6": "Macrophage",
    "7": "Dendritic cell",
    "8": "Alveolar epithelial",
    "9": "Monocyte",
    "10": "Mast cell",
    "11": "Lymphatic endothelial",
    "12": "B cell",
    "13": "Fibroblast",
    "14": "Basal epithelial",
    "15": "Basal epithelial",
    "16": "Basal epithelial",
    "17": "Basal epithelial",
    "18": "Basal epithelial",
    "19": "Basal epithelial",
    "20": "Basal epithelial",
    "21": "Basal epithelial",
    "22": "Basal epithelial",
    "23": "PNEC",
    "24": "Basal epithelial",
    "25": "Basal epithelial",
    "26": "Basal epithelial",
    "27": "Basal epithelial",
    "28": "Basal epithelial",
    "29": "Basal epithelial",
    "30": "Basal epithelial",
    "31": "Basal epithelial",
    "32": "Basal epithelial",
    "33": "Basal epithelial",
    "34": "Basal epithelial",
}


def cluster_annotations_table(annotations: dict[str, str] | None = None) -> pd.DataFrame:
    mapping = annotations or CLUSTER_ANNOTATIONS
    rows = [{"leiden": k, "cell_type_scanpy": v} for k, v in mapping.items()]
    out = pd.DataFrame(rows)
    out["_sort"] = out["leiden"].astype(int)
    return out.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)


def plot_umap_leiden_vs_scanpy(adata: ad.AnnData) -> Figure:
    """Side-by-side UMAP: Leiden clusters vs manual Scanpy labels (notebook style)."""
    if "X_umap" not in adata.obsm:
        raise ValueError("adata.obsm['X_umap'] missing — run clustering notebook first.")
    sc.pl.umap(
        adata,
        color=["leiden", "cell_type_scanpy"],
        wspace=0.4,
        title=["Leiden", "Manual scanpy"],
        show=False,
    )
    return plt.gcf()


def _compact_figure(fig: Figure, width: float, height: float, *, bottom: float = 0.22) -> Figure:
    fig.set_size_inches(width, height, forward=True)
    fig.subplots_adjust(bottom=bottom, top=0.92, left=0.08, right=0.98)
    return fig


def plot_ncount_violin(
    adata: ad.AnnData,
    *,
    key: str = "nCount",
    groupby: str = "cell_type_rule",
    rotation: float = 90,
    figsize: tuple[float, float] | None = None,
) -> Figure:
    """Transcript counts per cell type (notebook-style violin plot)."""
    if key not in adata.obs.columns:
        raise ValueError(f"adata.obs['{key}'] missing.")
    if groupby not in adata.obs.columns:
        raise ValueError(f"adata.obs['{groupby}'] missing.")
    sc.pl.violin(adata, keys=key, groupby=groupby, rotation=rotation, show=False)
    fig = plt.gcf()
    if figsize is None:
        n_groups = adata.obs[groupby].astype(str).nunique()
        width = min(max(7.0, n_groups * 0.45), 11.0)
        height = 3.0
        figsize = (width, height)
    return _compact_figure(fig, *figsize)


def unassigned_ncount_qc(
    adata: ad.AnnData,
    *,
    key: str = "nCount",
    rule_key: str = "cell_type_rule",
    unassigned_labels: tuple[str, ...] = ("Unassigned",),
    threshold: float = 5,
    figsize: tuple[float, float] = (3.5, 2.5),
) -> tuple[Figure | None, int, float]:
    """
    QC for rule-based Unassigned cells: violin of nCount and % below threshold.

    Returns (figure, n_unassigned, pct_below_threshold).
    """
    if key not in adata.obs.columns:
        raise ValueError(f"adata.obs['{key}'] missing.")
    if rule_key not in adata.obs.columns:
        raise ValueError(f"adata.obs['{rule_key}'] missing.")

    labels = adata.obs[rule_key].astype(str)
    unassigned = adata[labels.isin(unassigned_labels)].copy()
    n_unassigned = unassigned.n_obs
    if n_unassigned == 0:
        return None, 0, float("nan")

    ncounts = pd.to_numeric(unassigned.obs[key], errors="coerce")
    pct_below = float((ncounts < threshold).mean() * 100)

    sc.pl.violin(unassigned, keys=key, show=False)
    fig = _compact_figure(plt.gcf(), *figsize, bottom=0.15)
    return fig, n_unassigned, pct_below


def plot_rank_genes_dotplot(adata: ad.AnnData, n_genes: int = 5) -> Figure:
    """Top DE genes per Leiden cluster (notebook style)."""
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("adata.uns['rank_genes_groups'] missing — run clustering notebook first.")
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=n_genes,
        groupby="leiden",
        standard_scale="var",
        show=False,
    )
    return plt.gcf()
