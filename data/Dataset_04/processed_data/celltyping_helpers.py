"""Helpers for Dataset_04 celltyping notebook — load config and run panel logic."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from lib.config import load_yaml

_DEFAULT_CONFIG = Path(__file__).resolve().parent / "celltyping_config.yaml"


def load_celltyping_config(path: str | Path | None = None) -> dict[str, Any]:
    cfg = load_yaml(path or _DEFAULT_CONFIG)
    required = ("qc", "scanpy", "cluster_annotations", "rule_panels", "rule", "obs_keys")
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"celltyping_config.yaml missing sections: {missing}")
    return cfg


def qc_mask_from_config(adata: ad.AnnData, cfg: dict[str, Any]) -> pd.Series:
    qc = cfg["qc"]
    return (
        (adata.obs["nFeature"] >= qc["min_nFeature"])
        & (adata.obs["nCount"] >= qc["min_nCount"])
    )


def apply_spatial_crop(adata: ad.AnnData, cfg: dict[str, Any]) -> ad.AnnData:
    crop = cfg.get("spatial_crop", {})
    if not crop.get("enabled", False):
        return adata
    xy = adata.obsm["global"]
    mask = (
        (xy[:, 0] >= crop["x_min"])
        & (xy[:, 0] <= crop["x_max"])
        & (xy[:, 1] >= crop["y_min"])
        & (xy[:, 1] <= crop["y_max"])
    )
    return adata[mask].copy()


def _expression_matrix(adata: ad.AnnData) -> np.ndarray:
    x = adata.X
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def detection_rate(
    adata: ad.AnnData,
    genes: list[str],
    *,
    count_threshold: float = 0.0,
) -> pd.Series:
    """Fraction of cells with count > threshold for each gene present in adata."""
    present = [g for g in genes if g in adata.var_names]
    rates: dict[str, float] = {}
    dense = _expression_matrix(adata)
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    for gene in present:
        idx = gene_to_idx[gene]
        rates[gene] = float((dense[:, idx] > count_threshold).mean())
    return pd.Series(rates).sort_values(ascending=False)


def panel_assignment_counts(
    adata: ad.AnnData,
    genes: list[str],
    *,
    count_threshold: float = 0.0,
    min_hits: int | None = None,
) -> int:
    """Count cells passing a panel rule in isolation (no cross-panel competition)."""
    if not genes:
        return 0
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
    if not idx:
        return 0
    dense = _expression_matrix(adata)
    hits = (dense[:, idx] > count_threshold).sum(axis=1)
    required = len(idx) if min_hits is None else min_hits
    return int((hits >= required).sum())


def panel_dropout_curve(
    adata: ad.AnnData,
    candidates: list[str],
    *,
    max_k: int = 6,
    max_candidates: int = 10,
    count_threshold: float = 0.0,
    panel: str = "",
) -> pd.DataFrame:
    """Sweep top-k genes (by detection) and record assignment dropouts."""
    ranked = detection_rate(
        adata,
        candidates,
        count_threshold=count_threshold,
    ).head(max_candidates)
    if ranked.empty:
        return pd.DataFrame(
            columns=[
                "panel", "k", "genes", "n_assigned", "pct_assigned",
                "dropout", "dropout_pct", "is_recommended",
            ]
        )

    genes_ordered = ranked.index.tolist()
    max_k = min(max_k, len(genes_ordered))
    n_cells = adata.n_obs
    rows: list[dict[str, Any]] = []
    prev_assigned: int | None = None

    for k in range(1, max_k + 1):
        genes_k = genes_ordered[:k]
        n_assigned = panel_assignment_counts(
            adata,
            genes_k,
            count_threshold=count_threshold,
            min_hits=k,
        )
        dropout = 0 if prev_assigned is None else prev_assigned - n_assigned
        dropout_pct = 0.0 if prev_assigned in (None, 0) else dropout / prev_assigned
        rows.append({
            "panel": panel,
            "k": k,
            "genes": genes_k,
            "n_assigned": n_assigned,
            "pct_assigned": n_assigned / n_cells if n_cells else 0.0,
            "dropout": dropout,
            "dropout_pct": dropout_pct,
            "is_recommended": False,
        })
        prev_assigned = n_assigned

    return pd.DataFrame(rows)


def recommend_panel_markers(
    curve: pd.DataFrame,
    *,
    method: str = "before_max_dropout",
    min_assigned_cells: int = 50,
) -> dict[str, Any]:
    """Recommend gene set at the elbow before the largest assignment dropout."""
    if curve.empty:
        return {
            "genes": [],
            "min_hits": 0,
            "k_recommended": 0,
            "k_cliff": 0,
            "n_assigned": 0,
            "method": method,
        }

    if method != "before_max_dropout":
        raise ValueError(f"Unknown curation method: {method}")

    cliff_rows = curve.loc[curve["k"] >= 2]
    if cliff_rows.empty:
        row = curve.iloc[0]
        return {
            "genes": list(row["genes"]),
            "min_hits": int(row["k"]),
            "k_recommended": int(row["k"]),
            "k_cliff": 1,
            "n_assigned": int(row["n_assigned"]),
            "method": method,
        }

    k_cliff = int(cliff_rows.loc[cliff_rows["dropout"].idxmax(), "k"])
    k_recommended = max(1, k_cliff - 1)
    rec_row = curve.loc[curve["k"] == k_recommended].iloc[0]

    if int(rec_row["n_assigned"]) < min_assigned_cells and k_recommended > 1:
        k_recommended = max(1, k_recommended - 1)
        rec_row = curve.loc[curve["k"] == k_recommended].iloc[0]

    return {
        "genes": list(rec_row["genes"]),
        "min_hits": k_recommended,
        "k_recommended": k_recommended,
        "k_cliff": k_cliff,
        "n_assigned": int(rec_row["n_assigned"]),
        "method": method,
    }


def curate_rule_panels(
    adata: ad.AnnData,
    rule_panels: dict[str, Any],
    rule_cfg: dict[str, Any],
) -> dict[str, Any]:
    """Curate marker sets for all panels; return recommendations and dropout curves."""
    max_k = int(rule_cfg.get("max_k", 6))
    max_candidates = int(rule_cfg.get("max_candidates", 10))
    count_threshold = float(rule_cfg.get("count_threshold", 0.0))
    method = rule_cfg.get("method", "before_max_dropout")
    min_assigned_cells = int(rule_cfg.get("min_assigned_cells", 50))

    recommendations: dict[str, dict[str, Any]] = {}
    curve_frames: list[pd.DataFrame] = []

    for label, panel in rule_panels.items():
        if panel.get("markers"):
            genes = list(panel["markers"])
            n_assigned = panel_assignment_counts(
                adata,
                genes,
                count_threshold=count_threshold,
                min_hits=int(panel.get("min_hits", len(genes))),
            )
            recommendations[label] = {
                "genes": genes,
                "min_hits": int(panel.get("min_hits", len(genes))),
                "k_recommended": len(genes),
                "k_cliff": len(genes),
                "n_assigned": n_assigned,
                "method": "fixed_markers",
                "missing_genes": [g for g in genes if g not in adata.var_names],
            }
            continue

        candidates = panel.get("candidates", [])
        missing = [g for g in candidates if g not in adata.var_names]
        curve = panel_dropout_curve(
            adata,
            candidates,
            max_k=max_k,
            max_candidates=max_candidates,
            count_threshold=count_threshold,
            panel=label,
        )
        rec = recommend_panel_markers(
            curve,
            method=method,
            min_assigned_cells=min_assigned_cells,
        )
        rec["missing_genes"] = missing
        if not curve.empty and rec["k_recommended"] > 0:
            curve.loc[curve["k"] == rec["k_recommended"], "is_recommended"] = True
        recommendations[label] = rec
        curve_frames.append(curve)

    summary_rows = []
    for label, rec in recommendations.items():
        summary_rows.append({
            "panel": label,
            "recommended_genes": rec["genes"],
            "min_hits": rec["min_hits"],
            "k_recommended": rec["k_recommended"],
            "k_cliff": rec["k_cliff"],
            "n_assigned": rec["n_assigned"],
            "method": rec["method"],
            "missing_genes": rec.get("missing_genes", []),
        })

    curves = pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame()
    return {
        "recommendations": recommendations,
        "curves": curves,
        "summary_table": pd.DataFrame(summary_rows),
    }


def format_curated_yaml_snippet(curated: dict[str, Any]) -> str:
    """Return YAML snippet to paste into rule_panels after curation."""
    lines = ["# Paste under rule_panels in celltyping_config.yaml:", "rule_panels:"]
    for label, rec in curated["recommendations"].items():
        lines.append(f"  {label}:")
        lines.append("    markers:")
        for gene in rec["genes"]:
            lines.append(f"      - {gene}")
        lines.append(f"    min_hits: {rec['min_hits']}")
    return "\n".join(lines)


def resolve_panel_markers(
    adata: ad.AnnData,
    rule_panels: dict[str, Any],
    rule_cfg: dict[str, Any],
    curation_cfg: dict[str, Any] | None = None,
) -> tuple[dict[str, list[str]], dict[str, int]]:
    """Resolve marker lists and per-panel min_hits from config."""
    curation_cfg = curation_cfg or {}
    count_threshold = float(
        curation_cfg.get("count_threshold", rule_cfg.get("count_threshold", 0.0))
    )
    top_n = int(rule_cfg.get("top_n_markers", 4))
    default_min_hits = int(rule_cfg.get("min_hits", top_n))

    markers: dict[str, list[str]] = {}
    min_hits: dict[str, int] = {}

    for label, panel in rule_panels.items():
        if panel.get("markers"):
            genes = list(panel["markers"])
        else:
            candidates = panel.get("candidates", [])
            rates = detection_rate(adata, candidates, count_threshold=count_threshold)
            genes = rates.head(top_n).index.tolist()

        markers[label] = genes
        min_hits[label] = int(panel.get("min_hits", min(default_min_hits, len(genes) or 1)))

    return markers, min_hits


def select_panel_markers(
    adata: ad.AnnData,
    rule_panels: dict[str, Any],
    *,
    top_n: int,
    count_threshold: float = 0.0,
) -> dict[str, list[str]]:
    """Return resolved marker list per panel label (fixed markers or top-N by detection)."""
    selected: dict[str, list[str]] = {}
    for label, panel in rule_panels.items():
        if panel.get("markers"):
            selected[label] = list(panel["markers"])
            continue
        candidates = panel.get("candidates", [])
        rates = detection_rate(adata, candidates, count_threshold=count_threshold)
        selected[label] = rates.head(top_n).index.tolist()
    return selected


def flatten_panel_markers(panels: dict[str, list[str]]) -> list[str]:
    genes: list[str] = []
    for panel_genes in panels.values():
        genes.extend(panel_genes)
    return list(dict.fromkeys(genes))


def panel_detection_summary(
    adata: ad.AnnData,
    rule_panels: dict[str, Any],
    *,
    top_n: int,
    count_threshold: float = 0.0,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label, panel in rule_panels.items():
        genes = panel.get("markers") or panel.get("candidates", [])
        rates = detection_rate(adata, genes, count_threshold=count_threshold)
        for rank, (gene, rate) in enumerate(rates.head(top_n).items(), start=1):
            rows.append({
                "panel": label,
                "gene": gene,
                "detection_rate": rate,
                "rank": rank,
                "selected": rank <= top_n and not panel.get("markers"),
            })
    return pd.DataFrame(rows)


def assign_rule_labels_from_panels(
    adata: ad.AnnData,
    panel_markers: dict[str, list[str]],
    *,
    min_hits: int,
    count_threshold: float = 0.0,
    panel_min_hits: dict[str, int] | None = None,
    ambiguous_label: str = "Ambiguous",
    unassigned_label: str = "Unassigned",
) -> np.ndarray:
    """Assign one label per cell from multi-panel marker rules."""
    labels = np.full(adata.n_obs, unassigned_label, dtype=object)
    if not panel_markers:
        return labels

    gene_to_idx = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    panel_names = list(panel_markers.keys())
    hit_matrix = np.zeros((adata.n_obs, len(panel_names)), dtype=int)

    dense = _expression_matrix(adata)

    for j, panel_name in enumerate(panel_names):
        idx = [gene_to_idx[g] for g in panel_markers[panel_name] if g in gene_to_idx]
        if not idx:
            continue
        hit_matrix[:, j] = (dense[:, idx] > count_threshold).sum(axis=1)

    panel_thresholds = {
        name: (panel_min_hits or {}).get(name, min_hits)
        for name in panel_names
    }

    for i in range(adata.n_obs):
        hits = hit_matrix[i]
        qualifying = [
            (panel_names[j], int(hits[j]))
            for j in range(len(panel_names))
            if hits[j] >= panel_thresholds[panel_names[j]]
        ]
        if not qualifying:
            continue
        best = max(h for _, h in qualifying)
        winners = [name for name, h in qualifying if h == best]
        labels[i] = winners[0] if len(winners) == 1 else ambiguous_label

    return labels
