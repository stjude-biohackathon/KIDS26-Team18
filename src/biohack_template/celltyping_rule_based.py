"""Rule-based cell typing: >= min_hits panel genes with count > threshold."""

from __future__ import annotations

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from biohack_template.panels import load_demo_marker_panels


def _expression_matrix(adata: ad.AnnData) -> np.ndarray:
    x = adata.X
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def assign_celltypes_rule_based(
    adata: ad.AnnData,
    panels: dict[str, list[str]] | None = None,
    min_hits: int = 3,
    count_threshold: float = 1.0,
    obs_key: str = "cell_type_rule",
) -> ad.AnnData:
    """Assign cell types when >= min_hits panel genes exceed count_threshold."""
    if panels is None:
        panels = load_demo_marker_panels()["marker_panels"]

    x = _expression_matrix(adata)
    gene_to_idx = {g: i for i, g in enumerate(adata.var_names.astype(str))}

    assignments: list[str] = []
    hit_records: list[int] = []

    for row in x:
        best_type = "Unassigned"
        best_hits = -1
        for cell_type, genes in panels.items():
            hits = 0
            for gene in genes:
                idx = gene_to_idx.get(gene)
                if idx is not None and row[idx] > count_threshold:
                    hits += 1
            if hits > best_hits:
                best_hits = hits
                best_type = cell_type
            elif hits == best_hits and hits >= min_hits:
                # Deterministic tie-break by label order.
                if best_type == "Unassigned" or cell_type < best_type:
                    best_type = cell_type

        if best_hits < min_hits:
            best_type = "Unassigned"
        assignments.append(best_type)
        hit_records.append(max(best_hits, 0))

    adata = adata.copy()
    adata.obs[obs_key] = pd.Categorical(assignments)
    adata.obs[f"{obs_key}_panel_hits"] = hit_records
    return adata


def rule_based_summary(adata: ad.AnnData, obs_key: str = "cell_type_rule") -> pd.DataFrame:
    counts = adata.obs[obs_key].value_counts()
    return counts.rename_axis("cell_type").reset_index(name="n_cells")
