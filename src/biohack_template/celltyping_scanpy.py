"""Scanpy clustering + marker-panel cluster annotation for demo template."""

from __future__ import annotations

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from biohack_template.panels import load_demo_marker_panels


DEFAULT_SCANPY_CFG: dict[str, Any] = {
    "n_top_hvg": 2000,
    "n_pcs": 20,
    "n_neighbors": 15,
    "leiden_resolution": 0.5,
    "rank_genes_top_n": 50,
}


def _rank_genes_to_marker_df(adata: ad.AnnData, cluster_key: str, top_n: int = 50) -> pd.DataFrame:
    if "rank_genes_groups" not in adata.uns:
        raise ValueError("adata.uns['rank_genes_groups'] missing; run rank_genes_groups first.")

    rg = adata.uns["rank_genes_groups"]
    groups = list(rg["names"].dtype.names)
    records: list[dict[str, Any]] = []

    for group in groups:
        names = rg["names"][group][:top_n]
        scores = rg["scores"][group][:top_n] if "scores" in rg else rg["logfoldchanges"][group][:top_n]
        for rank, (gene, score) in enumerate(zip(names, scores, strict=False)):
            if gene is None or (isinstance(gene, float) and np.isnan(gene)):
                continue
            records.append({
                "group": str(group),
                "names": str(gene),
                "scores": float(score),
                "rank": rank,
            })

    return pd.DataFrame.from_records(records)


def _score_demo_cluster(sub_df: pd.DataFrame, panels: dict[str, list[str]]) -> dict[str, Any]:
    """Score cluster against demo panels directly (notebook labels)."""
    tmp = sub_df.copy()
    tmp["names"] = tmp["names"].astype(str)
    tmp["gene_upper"] = tmp["names"].str.upper()
    tmp["score_use"] = pd.to_numeric(tmp["scores"], errors="coerce").fillna(0).clip(lower=0)

    panel_scores: dict[str, float] = {}
    panel_anchors: dict[str, list[str]] = {}
    for panel_name, genes in panels.items():
        genes_upper = {g.upper() for g in genes}
        hit_df = tmp.loc[tmp["gene_upper"].isin(genes_upper)]
        panel_scores[panel_name] = float(hit_df["score_use"].sum())
        panel_anchors[panel_name] = hit_df.sort_values("scores", ascending=False)["names"].head(4).tolist()

    sorted_panels = sorted(panel_scores.items(), key=lambda x: x[1], reverse=True)
    best_panel, best_score = sorted_panels[0]
    second_panel, second_score = sorted_panels[1]
    anchors = panel_anchors.get(best_panel, [])

    n_anchor = len(anchors)
    score_gap = best_score - second_score
    if n_anchor >= 3 and best_score >= 0.5 and score_gap >= 0.1:
        confidence = "High"
    elif n_anchor >= 2 and best_score >= 0.2:
        confidence = "Medium"
    elif n_anchor >= 1:
        confidence = "Low-Medium"
    else:
        confidence = "Low"
        best_panel = "Unassigned"
        anchors = tmp.sort_values("scores", ascending=False)["names"].head(5).tolist()

    return {
        "cell_type_auto": best_panel,
        "confidence": confidence,
        "anchor_markers": ";".join(anchors),
        "note": "",
        "best_panel_score": best_score,
        "second_panel": second_panel,
        "second_panel_score": second_score,
    }


def run_scanpy_clustering(
    adata: ad.AnnData,
    cfg: dict[str, Any] | None = None,
    seed: int = 42,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Run lightweight scanpy workflow and return annotated adata + cluster table."""
    cfg = {**DEFAULT_SCANPY_CFG, **(cfg or {})}
    panels = load_demo_marker_panels()["marker_panels"]

    adata_sc = adata.copy()
    sc.pp.normalize_total(adata_sc, target_sum=1e4)
    sc.pp.log1p(adata_sc)

    n_top = min(int(cfg["n_top_hvg"]), adata_sc.n_vars - 1)
    sc.pp.highly_variable_genes(adata_sc, n_top_genes=max(n_top, 10), subset=True)
    sc.pp.scale(adata_sc, max_value=10)
    sc.tl.pca(adata_sc, n_comps=min(int(cfg["n_pcs"]), adata_sc.n_obs - 1, adata_sc.n_vars - 1))
    sc.pp.neighbors(adata_sc, n_neighbors=int(cfg["n_neighbors"]))
    sc.tl.leiden(
        adata_sc,
        resolution=float(cfg["leiden_resolution"]),
        key_added="leiden_demo",
        random_state=seed,
    )
    sc.tl.rank_genes_groups(adata_sc, groupby="leiden_demo", method="wilcoxon", use_raw=False)

    marker_df = _rank_genes_to_marker_df(adata_sc, "leiden_demo", top_n=int(cfg["rank_genes_top_n"]))
    cluster_rows: list[dict[str, Any]] = []
    for cluster_id, sub_df in marker_df.groupby("group", observed=True):
        result = _score_demo_cluster(sub_df, panels)
        cluster_rows.append({
            "cluster": str(cluster_id),
            "cell_type_scanpy": result["cell_type_auto"],
            "confidence": result["confidence"],
            "anchor_markers": result["anchor_markers"],
            "best_panel_score": result["best_panel_score"],
            "second_panel": result["second_panel"],
            "second_panel_score": result["second_panel_score"],
        })

    cluster_anno = pd.DataFrame(cluster_rows)
    cluster_map = dict(zip(cluster_anno["cluster"].astype(str), cluster_anno["cell_type_scanpy"]))
    adata_sc.obs["cell_type_scanpy"] = adata_sc.obs["leiden_demo"].astype(str).map(cluster_map).fillna("Unassigned")
    adata_sc.obs["cell_type_scanpy"] = adata_sc.obs["cell_type_scanpy"].astype("category")

    adata_sc.uns["scanpy_demo_cluster_annotations"] = cluster_anno
    return adata_sc, cluster_anno


def assign_celltypes_scanpy(
    adata: ad.AnnData,
    cfg: dict[str, Any] | None = None,
    seed: int = 42,
) -> ad.AnnData:
    """Run scanpy clustering workflow and write obs['cell_type_scanpy']."""
    adata_sc, cluster_anno = run_scanpy_clustering(adata, cfg=cfg, seed=seed)
    adata = adata.copy()
    adata.obs["leiden_demo"] = adata_sc.obs["leiden_demo"].astype(str)
    adata.obs["cell_type_scanpy"] = adata_sc.obs["cell_type_scanpy"].astype(str)
    adata.uns["scanpy_demo_cluster_annotations"] = adata_sc.uns["scanpy_demo_cluster_annotations"]
    if "X_umap" in adata_sc.obsm:
        adata.obsm["X_umap_demo"] = adata_sc.obsm["X_umap"]
    if "rank_genes_groups" in adata_sc.uns:
        adata.uns["rank_genes_groups"] = adata_sc.uns["rank_genes_groups"]
    return adata


def build_scanpy_cache(
    adata: ad.AnnData,
    cfg: dict[str, Any] | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Build a lightweight cache payload for expensive scanpy demo steps."""
    adata_sc, cluster_anno = run_scanpy_clustering(adata, cfg=cfg, seed=seed)
    return {
        "leiden_demo": adata_sc.obs["leiden_demo"].astype(str),
        "cell_type_scanpy": adata_sc.obs["cell_type_scanpy"].astype(str),
        "cluster_annotations": cluster_anno,
    }


def apply_scanpy_cache(adata: ad.AnnData, cache: dict[str, Any]) -> ad.AnnData:
    adata = adata.copy()
    adata.obs["leiden_demo"] = cache["leiden_demo"].reindex(adata.obs_names).astype(str)
    adata.obs["cell_type_scanpy"] = cache["cell_type_scanpy"].reindex(adata.obs_names).astype(str)
    adata.uns["scanpy_demo_cluster_annotations"] = cache["cluster_annotations"]
    return adata


def run_expensive_scanpy_analysis(
    adata: ad.AnnData,
    cfg: dict[str, Any] | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Alias used by the notebook expensive-computation demo."""
    return build_scanpy_cache(adata, cfg=cfg, seed=seed)
