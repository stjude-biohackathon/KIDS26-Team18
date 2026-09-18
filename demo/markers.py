"""Marker panel definitions for Dataset 01 rule-based typing demo."""

from __future__ import annotations

from typing import Any

import anndata as ad
import numpy as np

# Same candidate panels as 01_celltyping_compare.ipynb
CANDIDATE_MARKERS: dict[str, list[str]] = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Myofibroblast": ["ACTA2", "TAGLN", "MYLK", "TPM2", "COL1A1", "COL3A1"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "EMCN", "KDR", "ERG", "RAMP2"],
    "Plasma cell": ["JCHAIN", "MZB1", "SDC1", "DERL3", "TNFRSF17", "IGKC"],
    "T cell": ["CD3D", "CD3E", "CD3G", "TRAC", "LCK", "IL32"],
    "Macrophage": ["CD68", "C1QA", "C1QB", "C1QC", "CSF1R", "FCER1G", "CTSD"],
    "Dendritic cell": ["FCER1A", "CD1C", "CLEC10A", "HLA-DRA", "HLA-DPA1", "CST3"],
    "Alveolar epithelial": ["SFTPC", "SFTPA1", "SFTPA2", "SFTPB", "ABCA3", "NAPSA", "LPCAT1"],
    "Monocyte": ["FCN1", "S100A8", "S100A9", "LYZ", "CTSS", "VCAN", "CTSD"],
    "Mast cell": ["TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2", "HDC"],
    "Lymphatic endothelial": ["PROX1", "LYVE1", "PDPN", "FLT4", "CCL21", "CCL14"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD37", "CD74", "HLA-DRA", "CD19"],
    "Basal epithelial": ["KRT5", "KRT14", "KRT15", "KRT17", "TP63", "KRT19", "NGFR"],
    "PNEC": ["CHGA", "CHGB", "ASCL1", "INSM1", "SYP", "CALCA", "NCAM1"],
}


def _detection_rates(adata: ad.AnnData, genes: list[str]) -> dict[str, float]:
    present = [g for g in genes if g in adata.var_names]
    rates: dict[str, float] = {}
    for gene in present:
        x = adata[:, gene].X
        if hasattr(x, "toarray"):
            x = x.toarray().ravel()
        else:
            x = np.asarray(x).ravel()
        rates[gene] = float((x > 0).mean())
    return rates


def top_markers_per_type(adata: ad.AnnData, n: int = 4) -> dict[str, list[str]]:
    """Pick top-n detected markers per cell type (matches notebook logic)."""
    final_markers: dict[str, list[str]] = {}
    for cell_type, genes in CANDIDATE_MARKERS.items():
        rates = _detection_rates(adata, genes)
        if not rates:
            continue
        sorted_genes = sorted(rates, key=rates.get, reverse=True)
        final_markers[cell_type] = sorted_genes[:n]
    return final_markers


def _expression_row(adata: ad.AnnData) -> np.ndarray:
    x = adata.X
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def assign_rule_labels_tunable(
    adata: ad.AnnData,
    markers_dict: dict[str, list[str]],
    *,
    min_hits: int = 1,
    count_threshold: float = 0.0,
) -> np.ndarray:
    """
    Rule-based labels with tunable hit count and count threshold.

    For each cell, count panel genes with raw count > threshold.
    Assign the panel with the most hits (ties broken by total count, then panel order).
    """
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    x = _expression_row(adata)
    cell_types = list(markers_dict.keys())
    n_cells = adata.n_obs
    labels = np.full(n_cells, "Unassigned", dtype=object)

    for cell_idx in range(n_cells):
        row = x[cell_idx]
        best_type = "Unassigned"
        best_hits = -1
        best_score = -1.0

        for panel_idx, cell_type in enumerate(cell_types):
            idxs = [gene_to_idx[g] for g in markers_dict[cell_type] if g in gene_to_idx]
            if not idxs:
                continue
            counts = row[idxs]
            hits = int((counts > count_threshold).sum())
            score = float(counts.sum())
            if hits < min_hits:
                continue
            if hits > best_hits or (hits == best_hits and score > best_score):
                best_hits = hits
                best_score = score
                best_type = cell_type
            elif hits == best_hits and score == best_score and best_type != "Unassigned":
                if panel_idx < cell_types.index(best_type):
                    best_type = cell_type

        labels[cell_idx] = best_type

    return labels


def assign_rule_labels_notebook(
    adata: ad.AnnData,
    markers_dict: dict[str, list[str]],
) -> np.ndarray:
    """Exact cellruler.rulebased logic: highest total count among qualifying panels."""
    cell_types = list(markers_dict.keys())
    n_cells = adata.n_obs
    n_panels = len(cell_types)
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    x = _expression_row(adata)

    scores = np.zeros((n_panels, n_cells), dtype=float)
    for p, cell_type in enumerate(cell_types):
        genes = markers_dict[cell_type]
        idxs = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        if idxs:
            scores[p] = x[:, idxs].sum(axis=1)

    qualifies = scores > 0
    panel_order = np.arange(n_panels).reshape(-1, 1)
    rank = scores + panel_order * 1e-12
    rank = np.where(qualifies, rank, -np.inf)
    best_panel_idx = rank.argmax(axis=0)
    has_hit = qualifies.any(axis=0)

    labels = np.full(n_cells, "Unassigned", dtype=object)
    labels[has_hit] = np.array(cell_types, dtype=object)[best_panel_idx[has_hit]]
    return labels
