"""Spatial hypothesis tests for tumor-state / macrophage association."""

from __future__ import annotations

from typing import Any

import anndata as ad
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu


def compute_tumor_macrophage_distances(adata: ad.AnnData) -> pd.DataFrame:
    """Compute nearest-macrophage distance for each tumor cell."""
    macrophages = adata.obs["cell_type"] == "Macrophages"
    if macrophages.sum() == 0:
        raise ValueError("No macrophages found in the dataset.")

    tree = cKDTree(adata.obsm["spatial"][macrophages.to_numpy()])
    tumor_mask = adata.obs["cell_type"].str.startswith("Tumor").to_numpy()
    distances, _ = tree.query(adata.obsm["spatial"][tumor_mask], k=1)

    tumor_obs = adata.obs.loc[tumor_mask, ["core_id", "cell_type", "tumor_state"]].copy()
    tumor_obs["distance_to_nearest_macrophage"] = distances
    return tumor_obs.reset_index()


def compare_tumor_macrophage_distances(distance_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series, float, float]:
    """Summarize and test whether state A is closer to macrophages than state B."""
    summary_stats = (
        distance_df.groupby("tumor_state")["distance_to_nearest_macrophage"]
        .agg(n_cells="count", median="median", mean="mean")
        .reset_index()
    )

    state_a = distance_df.loc[distance_df["tumor_state"] == "state A", "distance_to_nearest_macrophage"]
    state_b = distance_df.loc[distance_df["tumor_state"] == "state B", "distance_to_nearest_macrophage"]
    mw_stat, mw_pvalue = mannwhitneyu(state_a, state_b, alternative="less")

    test_summary = pd.DataFrame({
        "comparison": ["state A vs state B (A closer?)"],
        "test": ["Mann-Whitney U"],
        "statistic": [mw_stat],
        "p_value": [mw_pvalue],
    })
    return summary_stats, test_summary, state_a, state_b, mw_stat, mw_pvalue


def run_expensive_spatial_analysis(adata: ad.AnnData) -> dict[str, Any]:
    """Lightweight stand-in for an expensive neighborhood enrichment computation."""
    distances = compute_tumor_macrophage_distances(adata)
    enrichment = (
        distances.groupby("tumor_state")["distance_to_nearest_macrophage"]
        .agg(["count", "median", "mean"])
        .rename(columns={"count": "n_tumor_cells"})
    )
    return {"enrichment_summary": enrichment, "distances": distances}
