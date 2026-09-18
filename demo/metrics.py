"""Comparison metrics for scanpy vs rule-based cell typing."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def cohens_kappa(y1: pd.Series, y2: pd.Series) -> float:
    labels = sorted(set(y1.astype(str)) | set(y2.astype(str)))
    if not labels:
        return float("nan")

    conf = pd.crosstab(y1.astype(str), y2.astype(str), dropna=False).reindex(
        index=labels, columns=labels, fill_value=0
    )
    n = conf.to_numpy().sum()
    if n == 0:
        return float("nan")

    po = np.trace(conf.to_numpy()) / n
    pe = (conf.sum(axis=1).to_numpy() * conf.sum(axis=0).to_numpy()).sum() / (n**2)
    if np.isclose(pe, 1.0):
        return 1.0 if np.isclose(po, 1.0) else 0.0
    return float((po - pe) / (1 - pe))


def accuracy(y_true: pd.Series, y_pred: pd.Series) -> float:
    return float((y_true.astype(str).to_numpy() == y_pred.astype(str).to_numpy()).mean())


def compare_methods(
    obs: pd.DataFrame,
    *,
    method_a_key: str = "cell_type_scanpy",
    method_b_key: str = "cell_type_rule",
    method_a_label: str = "Scanpy",
    method_b_label: str = "Rule-based",
    truth_key: str = "final_CT",
) -> dict[str, Any]:
    y_a = obs[method_a_key].astype(str)
    y_b = obs[method_b_key].astype(str)
    has_truth = truth_key in obs.columns
    y_truth = obs[truth_key].astype(str) if has_truth else None

    summary_rows = [
        {
            "comparison": f"{method_a_label} vs {method_b_label}",
            "accuracy": accuracy(y_a, y_b),
            "cohen_kappa": cohens_kappa(y_a, y_b),
            "n_cells": len(obs),
        }
    ]
    if has_truth and y_truth is not None:
        summary_rows.extend([
            {
                "comparison": f"{method_a_label} vs Published",
                "accuracy": accuracy(y_truth, y_a),
                "cohen_kappa": cohens_kappa(y_truth, y_a),
                "n_cells": len(obs),
            },
            {
                "comparison": f"{method_b_label} vs Published",
                "accuracy": accuracy(y_truth, y_b),
                "cohen_kappa": cohens_kappa(y_truth, y_b),
                "n_cells": len(obs),
            },
        ])

    disagree = y_a != y_b
    disagree_df = obs.loc[disagree, [method_a_key, method_b_key] + ([truth_key] if has_truth else [])].copy()
    disagree_df["status"] = "Disagree"

    return {
        "summary": pd.DataFrame(summary_rows),
        "crosstab": pd.crosstab(y_a, y_b, rownames=[method_a_label], colnames=[method_b_label]),
        "n_agree": int((~disagree).sum()),
        "n_disagree": int(disagree.sum()),
        "pct_agree": float((~disagree).mean() * 100),
        "disagree_df": disagree_df,
    }
