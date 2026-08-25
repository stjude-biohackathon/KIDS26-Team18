"""Compare rule-based vs scanpy cell typing against synthetic ground truth."""

from __future__ import annotations

from typing import Any

import anndata as ad
import numpy as np
import pandas as pd


def _cohens_kappa(y1: pd.Series, y2: pd.Series) -> float:
    labels = sorted(set(y1.astype(str)) | set(y2.astype(str)))
    if not labels:
        return float("nan")

    conf = pd.crosstab(y1.astype(str), y2.astype(str), dropna=False).reindex(index=labels, columns=labels, fill_value=0)
    n = conf.to_numpy().sum()
    if n == 0:
        return float("nan")

    po = np.trace(conf.to_numpy()) / n
    pe = (conf.sum(axis=1).to_numpy() * conf.sum(axis=0).to_numpy()).sum() / (n ** 2)
    if np.isclose(pe, 1.0):
        return 1.0 if np.isclose(po, 1.0) else 0.0
    return float((po - pe) / (1 - pe))


def _accuracy(y_true: pd.Series, y_pred: pd.Series) -> float:
    return float((y_true.astype(str).to_numpy() == y_pred.astype(str).to_numpy()).mean())


def _per_type_metrics(y_true: pd.Series, y_pred: pd.Series) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cell_type in sorted(y_true.astype(str).unique()):
        true_mask = y_true.astype(str) == cell_type
        pred_mask = y_pred.astype(str) == cell_type
        tp = int((true_mask & pred_mask).sum())
        fn = int((true_mask & ~pred_mask).sum())
        fp = int((~true_mask & pred_mask).sum())
        recall = tp / (tp + fn) if (tp + fn) else np.nan
        precision = tp / (tp + fp) if (tp + fp) else np.nan
        rows.append({
            "cell_type": cell_type,
            "true_cells": int(true_mask.sum()),
            "predicted_cells": int(pred_mask.sum()),
            "tp": tp,
            "recall": recall,
            "precision": precision,
        })
    return pd.DataFrame(rows)


def compare_celltyping_methods(
    adata: ad.AnnData,
    truth_key: str = "cell_type",
    rule_key: str = "cell_type_rule",
    scanpy_key: str = "cell_type_scanpy",
) -> dict[str, Any]:
    """Return summary metrics and cross-tabs for both annotation modes."""
    y_true = adata.obs[truth_key].astype(str)
    y_rule = adata.obs[rule_key].astype(str)
    y_scanpy = adata.obs[scanpy_key].astype(str)

    summary = pd.DataFrame([
        {
            "comparison": "rule vs ground truth",
            "accuracy": _accuracy(y_true, y_rule),
            "cohen_kappa": _cohens_kappa(y_true, y_rule),
            "n_cells": adata.n_obs,
        },
        {
            "comparison": "scanpy vs ground truth",
            "accuracy": _accuracy(y_true, y_scanpy),
            "cohen_kappa": _cohens_kappa(y_true, y_scanpy),
            "n_cells": adata.n_obs,
        },
        {
            "comparison": "scanpy vs rule",
            "accuracy": _accuracy(y_rule, y_scanpy),
            "cohen_kappa": _cohens_kappa(y_rule, y_scanpy),
            "n_cells": adata.n_obs,
        },
    ])

    crosstab = pd.crosstab(y_rule, y_scanpy, rownames=["cell_type_rule"], colnames=["cell_type_scanpy"])
    rule_vs_truth = pd.crosstab(y_true, y_rule, rownames=[truth_key], colnames=[rule_key])
    scanpy_vs_truth = pd.crosstab(y_true, y_scanpy, rownames=[truth_key], colnames=[scanpy_key])

    per_type_rule = _per_type_metrics(y_true, y_rule).assign(method="rule")
    per_type_scanpy = _per_type_metrics(y_true, y_scanpy).assign(method="scanpy")
    per_type = pd.concat([per_type_rule, per_type_scanpy], ignore_index=True)

    return {
        "summary": summary,
        "crosstab": crosstab,
        "rule_vs_truth": rule_vs_truth,
        "scanpy_vs_truth": scanpy_vs_truth,
        "per_type": per_type,
    }
