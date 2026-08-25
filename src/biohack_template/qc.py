"""QC helpers wrapping cell_typing_v3 clustering utilities."""

from __future__ import annotations

from typing import Any

import anndata as ad
import pandas as pd

from lib.qc_common import run_qc_flagging_and_filtering


def cells_per_core_summary(adata: ad.AnnData, core_col: str = "core_id") -> pd.Series:
    return adata.obs[core_col].value_counts().sort_index()


def basic_qc_filter(
    adata: ad.AnnData,
    min_total_counts: float = 35.0,
    min_genes: int | None = None,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Apply QC using the shared clustering_common QC pipeline."""
    cfg: dict[str, Any] = {
        "do_qc": True,
        "total_counts_min": float(min_total_counts),
        "filter_qc_pass": True,
        "compute_mt": False,
        "compute_ribo": False,
    }
    if min_genes is not None:
        cfg["n_genes_min"] = int(min_genes)

    n_before = adata.n_obs
    filtered = run_qc_flagging_and_filtering(adata.copy(), cfg)
    summary = pd.DataFrame({
        "metric": ["cells_before", "cells_after", "cells_removed", "min_total_counts"],
        "value": [n_before, filtered.n_obs, n_before - filtered.n_obs, min_total_counts],
    })
    return filtered, summary
