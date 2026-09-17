"""Synthetic Xenium-like TMA demo data for BioHack notebook templates."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from biohack_template.panels import demo_panel_genes, load_demo_marker_panels
from lib.qc_common import ensure_output_dir


def _simulate_expression(
    rng: np.random.Generator,
    gene_names: list[str],
    cell_type: str,
    panels: dict[str, list[str]],
    n_background: int = 20,
) -> np.ndarray:
    """Simulate counts with boosted on-panel markers and low off-panel noise."""
    n_genes = len(gene_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    expression = rng.poisson(lam=0.4, size=n_genes).astype(float)

    panel = panels[cell_type]
    for gene in panel:
        idx = gene_to_idx[gene]
        # Strong on-target signal so >=3/4 panel genes exceed count threshold.
        expression[idx] = rng.poisson(lam=8 if rng.random() < 0.85 else 1)

    # Weak spillover from related panels to create realistic ambiguity.
    other_types = [ct for ct in panels if ct != cell_type]
    spill_types = rng.choice(other_types, size=min(2, len(other_types)), replace=False)
    for spill_type in spill_types:
        for gene in rng.choice(panels[spill_type], size=2, replace=False):
            idx = gene_to_idx[gene]
            if expression[idx] <= 1:
                expression[idx] = float(rng.poisson(lam=1.2))

    # Generic background expression on non-panel genes.
    non_panel_idx = [i for i, g in enumerate(gene_names) if all(g not in pg for pg in panels.values())]
    if non_panel_idx:
        chosen = rng.choice(non_panel_idx, size=min(n_background, len(non_panel_idx)), replace=False)
        expression[chosen] = rng.poisson(lam=3, size=len(chosen)).astype(float)

    return expression


def make_synthetic_tma_adata(seed: int = 42) -> ad.AnnData:
    """Generate a small synthetic Xenium-like TMA dataset with marker gene expression."""
    rng = np.random.default_rng(seed)
    panel_cfg = load_demo_marker_panels()
    panels: dict[str, list[str]] = panel_cfg["marker_panels"]
    genes = demo_panel_genes(panel_cfg)
    core_ids = [f"core_{i}" for i in range(1, 5)]

    records: list[dict[str, Any]] = []
    cell_idx = 0

    for core_i, core_id in enumerate(core_ids):
        n_cells = rng.integers(220, 320)
        center_x = 250 + core_i * 600
        center_y = 250 + (core_i % 2) * 500

        macrophage_centers = np.column_stack([
            center_x + rng.normal(0, 40, 3),
            center_y + rng.normal(0, 40, 3),
        ])

        for _ in range(n_cells):
            cell_type = rng.choice(
                list(panels.keys()),
                p=[0.18, 0.17, 0.12, 0.10, 0.14, 0.17, 0.12],
            )

            if cell_type == "Macrophages":
                anchor = macrophage_centers[rng.integers(0, 3)]
                x, y = anchor + rng.normal(0, 25, size=2)
            elif cell_type == "Tumor state A" and core_i in (0, 2):
                anchor = macrophage_centers[rng.integers(0, 3)]
                x, y = anchor + rng.normal(0, 55, size=2)
            elif cell_type == "Tumor state B":
                x, y = [center_x, center_y] + rng.normal(0, 120, size=2)
            else:
                x, y = [center_x, center_y] + rng.normal(0, 150, size=2)

            expression = _simulate_expression(rng, genes, cell_type, panels)
            tumor_state = cell_type.replace("Tumor ", "") if cell_type.startswith("Tumor") else np.nan

            records.append({
                "cell_id": f"{core_id}_cell_{cell_idx:04d}",
                "core_id": core_id,
                "cell_type": cell_type,
                "tumor_state": tumor_state,
                "x": float(x),
                "y": float(y),
                "total_counts": float(expression.sum()),
                "expression": expression,
            })
            cell_idx += 1

    obs = pd.DataFrame(records).set_index("cell_id")
    x_matrix = np.vstack(obs.pop("expression").to_numpy())
    spatial = obs[["x", "y"]].to_numpy()
    obs = obs.drop(columns=["x", "y"])

    adata = ad.AnnData(X=x_matrix, obs=obs, var=pd.DataFrame(index=genes))
    adata.obsm["spatial"] = spatial
    adata.obs["coord_x"] = spatial[:, 0]
    adata.obs["coord_y"] = spatial[:, 1]
    adata.obs["sample"] = adata.obs["core_id"].astype(str)
    adata.uns["synthetic_demo"] = True
    return adata


def load_data(use_synthetic: bool = True, seed: int = 42) -> ad.AnnData:
    """Load analysis data. Template uses synthetic demo data by default."""
    if use_synthetic:
        return make_synthetic_tma_adata(seed=seed)
    raise NotImplementedError("Replace this branch with loading of real project data.")


def save_synthetic_h5ad(adata: ad.AnnData, path: str | Path) -> None:
    """Optional helper to persist demo data for offline reruns."""
    ensure_output_dir(str(path))
    adata.write_h5ad(path)
