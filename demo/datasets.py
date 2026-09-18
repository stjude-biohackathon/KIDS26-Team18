"""Dataset registry and loaders for the cell-typing comparison demo."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
D4_PROCESSED = REPO_ROOT / "data/Dataset_04/processed_data"


def _ensure_import_paths() -> None:
    for path in (REPO_ROOT, SRC_ROOT, D4_PROCESSED):
        path_str = str(path)
        if path_str not in sys.path:
            sys.path.insert(0, path_str)


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    label: str
    caption: str
    clustered_h5ad: Path
    reference_h5ad: Path | None = None
    reference_obs_cols: tuple[str, ...] = ()
    config_path: Path | None = None
    rule_engine: str = "dataset01"
    method_a_col: str = "cell_type_scanpy"
    method_b_col: str = "cell_type_rule"
    method_a_label: str = "Scanpy"
    method_b_label: str = "Rule-based"
    unassigned_labels: tuple[str, ...] = ("Unassigned",)
    unassigned_on: str = "method_b"
    count_col: str = "nCount"
    obs_only: bool = False
    label_normalization: bool = False


def unassigned_column(spec: DatasetSpec) -> str:
    return spec.method_a_col if spec.unassigned_on == "method_a" else spec.method_b_col


DATASETS: dict[str, DatasetSpec] = {
    "dataset01": DatasetSpec(
        dataset_id="dataset01",
        label="Dataset 01 — Xenium lung TMA5 (GSE250346)",
        caption=(
            "GSE250346 human lung Xenium TMA5 — compare Scanpy Leiden + manual annotation "
            "against transparent marker-panel rules, with optional published `final_CT`."
        ),
        clustered_h5ad=REPO_ROOT / "data/Dataset_01/processed_data/one_sample_clustered.h5ad",
        reference_h5ad=REPO_ROOT / "data/Dataset_01/processed_data/cell_typed_TMA5_subset.h5ad",
        reference_obs_cols=("final_CT",),
        rule_engine="dataset01",
    ),
    "dataset04": DatasetSpec(
        dataset_id="dataset04",
        label="Dataset 04 — CosMx kidney (GSE282026)",
        caption=(
            "GSE282026 human kidney CosMx — compare Scanpy Leiden + manual annotation "
            "against config-driven marker-panel rules (`celltyping_config.yaml`)."
        ),
        clustered_h5ad=D4_PROCESSED / "scanpy_cluster_cellbased.h5ad",
        config_path=D4_PROCESSED / "celltyping_config.yaml",
        rule_engine="dataset04",
    ),
    "dataset05": DatasetSpec(
        dataset_id="dataset05",
        label="Dataset 05 — Xenium skin (GSE301280)",
        caption=(
            "GSE301280 human skin Xenium — compare hierarchical cell typing (`hier_label`) "
            "against the authors' published labels (`cell_type`)."
        ),
        clustered_h5ad=REPO_ROOT / "data/Dataset_05/processed_data/adata_annotated.h5ad",
        rule_engine="none",
        method_a_col="hier_label",
        method_b_col="cell_type",
        method_a_label="Hierarchical",
        method_b_label="Author",
        unassigned_labels=("unassigned", "low_quality"),
        unassigned_on="method_a",
        count_col="transcript_counts",
        obs_only=True,
        label_normalization=True,
    ),
}


def uses_rule_engine(spec: DatasetSpec) -> bool:
    return spec.rule_engine not in ("none", "")


def get_dataset(dataset_id: str) -> DatasetSpec:
    if dataset_id not in DATASETS:
        raise KeyError(f"Unknown dataset: {dataset_id}")
    spec = DATASETS[dataset_id]
    if not spec.clustered_h5ad.is_file():
        raise FileNotFoundError(
            f"Missing `{spec.clustered_h5ad}`. Run the Dataset compare notebook and save the clustered h5ad."
        )
    return spec


def _ensure_spatial_coords(adata: ad.AnnData) -> ad.AnnData:
    if "coord_x" in adata.obs.columns and "coord_y" in adata.obs.columns:
        return adata
    if "x_centroid" in adata.obs.columns and "y_centroid" in adata.obs.columns:
        adata = adata.copy()
        adata.obs["coord_x"] = adata.obs["x_centroid"]
        adata.obs["coord_y"] = adata.obs["y_centroid"]
        return adata
    if "global" in adata.obsm:
        adata = adata.copy()
        adata.obs["coord_x"] = adata.obsm["global"][:, 0]
        adata.obs["coord_y"] = adata.obsm["global"][:, 1]
        return adata
    if "spatial" in adata.obsm:
        adata = adata.copy()
        adata.obs["coord_x"] = adata.obsm["spatial"][:, 0]
        adata.obs["coord_y"] = adata.obsm["spatial"][:, 1]
        return adata
    return adata


def _load_obs_only(spec: DatasetSpec) -> ad.AnnData:
    """Load obs + spatial coords without expression matrix (large annotation files)."""
    backed = ad.read_h5ad(spec.clustered_h5ad, backed="r")
    obs = backed.obs.copy()
    obsm: dict[str, np.ndarray] = {}
    if "spatial" in backed.obsm:
        obsm["spatial"] = np.asarray(backed.obsm["spatial"][()])
    backed.file.close()
    adata = ad.AnnData(X=sp.csr_matrix((obs.shape[0], 1)), obs=obs, obsm=obsm)
    return _ensure_spatial_coords(adata)


def load_demo_adata(dataset_id: str) -> ad.AnnData:
    spec = get_dataset(dataset_id)
    if spec.obs_only:
        return _load_obs_only(spec)

    adata = ad.read_h5ad(spec.clustered_h5ad)

    if spec.reference_h5ad and spec.reference_h5ad.is_file():
        ref = ad.read_h5ad(spec.reference_h5ad)
        join_key = "cell_id" if "cell_id" in adata.obs.columns else None
        cols_to_add = [
            c for c in spec.reference_obs_cols if c in ref.obs.columns and c not in adata.obs.columns
        ]
        if cols_to_add:
            if join_key and join_key in ref.obs.columns:
                ref_indexed = ref.obs.set_index(join_key)
                for col in cols_to_add:
                    adata.obs[col] = adata.obs[join_key].map(ref_indexed[col])
            elif adata.n_obs == ref.n_obs:
                for col in cols_to_add:
                    adata.obs[col] = ref.obs[col].to_numpy()

    return _ensure_spatial_coords(adata)


def load_celltyping_config(spec: DatasetSpec) -> dict[str, Any]:
    _ensure_import_paths()
    if spec.config_path is None:
        raise ValueError(f"{spec.dataset_id} has no celltyping config.")
    from celltyping_helpers import load_celltyping_config

    return load_celltyping_config(spec.config_path)


def cluster_annotations_for_dataset(spec: DatasetSpec) -> dict[str, str]:
    if spec.rule_engine == "dataset04":
        cfg = load_celltyping_config(spec)
        return {str(k): str(v) for k, v in cfg["cluster_annotations"].items()}
    from demo.cluster_plots import CLUSTER_ANNOTATIONS

    return dict(CLUSTER_ANNOTATIONS)


def cluster_annotations_table(spec: DatasetSpec) -> pd.DataFrame:
    annotations = cluster_annotations_for_dataset(spec)
    rows = [{"leiden": k, "cell_type_scanpy": v} for k, v in annotations.items()]
    out = pd.DataFrame(rows)
    out["_sort"] = out["leiden"].astype(int)
    return out.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)


def resolve_rule_panels(
    adata: ad.AnnData,
    spec: DatasetSpec,
) -> tuple[dict[str, list[str]], dict[str, int], dict[str, Any]]:
    """Return panel markers, per-panel min_hits, and rule config defaults."""
    if spec.rule_engine == "dataset01":
        from demo.markers import top_markers_per_type

        markers = top_markers_per_type(adata, n=4)
        min_hits = {label: 1 for label in markers}
        rule_cfg = {"min_hits": 1, "count_threshold": 0.0, "max_hits_slider": 4}
        return markers, min_hits, rule_cfg

    _ensure_import_paths()
    from celltyping_helpers import load_celltyping_config, resolve_panel_markers

    cfg = load_celltyping_config(spec.config_path)
    markers, panel_min_hits = resolve_panel_markers(
        adata,
        cfg["rule_panels"],
        cfg["rule"],
        cfg.get("curation"),
    )
    rule_cfg = {
        "min_hits": int(cfg["rule"].get("min_hits", 2)),
        "count_threshold": float(cfg["rule"].get("count_threshold", 0.0)),
        "max_hits_slider": int(cfg.get("curation", {}).get("max_k", 6)),
        "ambiguous_label": cfg["rule"].get("ambiguous_label", "Ambiguous"),
        "unassigned_label": cfg["rule"].get("unassigned_label", "Unassigned"),
    }
    return markers, panel_min_hits, rule_cfg


def assign_rule_labels(
    adata: ad.AnnData,
    spec: DatasetSpec,
    panel_markers: dict[str, list[str]],
    panel_min_hits: dict[str, int],
    *,
    min_hits: int,
    count_threshold: float,
    use_config_panels: bool = True,
) -> np.ndarray:
    if spec.rule_engine == "dataset01":
        from demo.markers import assign_rule_labels_notebook, assign_rule_labels_tunable

        if use_config_panels:
            return assign_rule_labels_notebook(adata, panel_markers)
        return assign_rule_labels_tunable(
            adata,
            panel_markers,
            min_hits=min_hits,
            count_threshold=count_threshold,
        )

    _ensure_import_paths()
    from celltyping_helpers import assign_rule_labels_from_panels

    cfg = load_celltyping_config(spec)
    rule_cfg = cfg["rule"]
    if use_config_panels:
        effective_min_hits = int(rule_cfg.get("min_hits", min_hits))
        effective_threshold = float(rule_cfg.get("count_threshold", count_threshold))
        panel_hits = panel_min_hits
    else:
        effective_min_hits = min_hits
        effective_threshold = count_threshold
        panel_hits = {label: min_hits for label in panel_markers}

    return assign_rule_labels_from_panels(
        adata,
        panel_markers,
        min_hits=effective_min_hits,
        count_threshold=effective_threshold,
        panel_min_hits=panel_hits,
        ambiguous_label=rule_cfg.get("ambiguous_label", "Ambiguous"),
        unassigned_label=rule_cfg.get("unassigned_label", "Unassigned"),
    )
