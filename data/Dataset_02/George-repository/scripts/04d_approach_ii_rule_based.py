#!/usr/bin/env python
"""Approach II: rule-based cell typing for FOV46.

A candidate cell type is called positive when at least MIN_POSITIVE_GENES
marker genes that are present in the CosMx panel have raw count > 0.
Candidates with fewer than MIN_POSITIVE_GENES available markers are disabled.

Outputs include marker diagnostics, rule labels, an Approach-II UMAP, an
Approach-I-vs-II side-by-side UMAP, and an assignment-status UMAP.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from config.markers import MARKER_PANEL
from src.typing import detection_rate
from src.plotting import (
    plot_umap_annotation,
    plot_umap_side_by_side,
    plot_spatial_annotation_pair,
)

# -----------------------------------------------------------------------------
# Paths / configuration
# -----------------------------------------------------------------------------
IN = ROOT / "results/FOV46/GSM9046088_FOV46_approach_i.h5ad"
OUT = ROOT / "results/FOV46"
TAB = OUT / "tables"
FIG = OUT / "figures" / "approach_ii"
SDATA = ROOT / "data/spatial/GSM9046088_CosMx_FOV46.zarr"

TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

LAYER = "counts"
MIN_POSITIVE_GENES = 3

# -----------------------------------------------------------------------------
# Load Approach-I object so Approach I and II share the exact same cells,
# preprocessing, PCA and UMAP coordinates.
# -----------------------------------------------------------------------------
print("=" * 70)
print("FOV46 APPROACH II — RULE-BASED CELL TYPING")
print("=" * 70)
print(f"Input: {IN}")
print(f"Rule: >= {MIN_POSITIVE_GENES} detected marker genes with raw count > 0")

a = sc.read_h5ad(IN)

if LAYER not in a.layers:
    raise KeyError(
        f"Required raw-count layer {LAYER!r} is missing. "
        "Approach II must use raw biological counts."
    )
if "X_umap" not in a.obsm:
    raise KeyError(
        "X_umap is missing. Run 04b Leiden sweep and 04c Approach I first."
    )

Xall = a.layers[LAYER]

# -----------------------------------------------------------------------------
# Marker detection rates, panel availability, and rule calls
# -----------------------------------------------------------------------------
rates = []
availability = []
positives = {}

for cell_type, genes in MARKER_PANEL.items():
    # Per-gene detection rate for interpretation/QC.
    rates_for_type = detection_rate(a, genes, layer=LAYER)
    rates.extend(
        {
            "cell_type": cell_type,
            "gene": gene,
            "detection_rate": value,
        }
        for gene, value in rates_for_type.items()
    )

    present = [gene for gene in genes if gene in a.var_names]
    missing = [gene for gene in genes if gene not in a.var_names]
    available = len(present) >= MIN_POSITIVE_GENES

    availability.append(
        {
            "cell_type": cell_type,
            "n_requested": len(genes),
            "n_present": len(present),
            "min_positive_genes": MIN_POSITIVE_GENES,
            "available": available,
            "present": ",".join(present),
            "missing": ",".join(missing),
        }
    )

    if not available:
        positives[cell_type] = np.zeros(a.n_obs, dtype=bool)
        a.obs[f"rule_{cell_type}_n_positive_genes"] = 0
        a.obs[f"rule_{cell_type}_positive"] = False
        continue

    idx = [a.var_names.get_loc(gene) for gene in present]
    x = Xall[:, idx]
    x = x.toarray() if sparse.issparse(x) else np.asarray(x)

    n_positive = np.asarray((x > 0).sum(axis=1)).ravel()
    is_positive = n_positive >= MIN_POSITIVE_GENES

    positives[cell_type] = is_positive
    a.obs[f"rule_{cell_type}_n_positive_genes"] = n_positive
    a.obs[f"rule_{cell_type}_positive"] = is_positive

# Save marker diagnostics.
pd.DataFrame(rates).to_csv(
    TAB / "rule_marker_detection_rates.csv", index=False
)
availability_df = pd.DataFrame(availability)
availability_df.to_csv(
    TAB / "rule_marker_availability.csv", index=False
)

# -----------------------------------------------------------------------------
# Final Approach-II label
# -----------------------------------------------------------------------------
cell_types = list(MARKER_PANEL.keys())
positive_matrix = np.column_stack([positives[ct] for ct in cell_types])
n_positive_panels = positive_matrix.sum(axis=1)
cell_type_names = np.asarray(cell_types, dtype=object)

labels = np.full(a.n_obs, "Unassigned", dtype=object)
exactly_one = n_positive_panels == 1
labels[exactly_one] = cell_type_names[
    positive_matrix[exactly_one].argmax(axis=1)
]
labels[n_positive_panels > 1] = "Ambiguous"

a.obs["cell_type_rule"] = pd.Categorical(labels)
a.obs["rule_n_positive_panels"] = n_positive_panels

label_counts = a.obs["cell_type_rule"].value_counts(dropna=False)
label_counts.to_csv(TAB / "rule_label_counts.csv")

# -----------------------------------------------------------------------------
# Assignment-status category
# -----------------------------------------------------------------------------
status = np.full(a.n_obs, "Assigned", dtype=object)
status[labels == "Unassigned"] = "Unassigned"
status[labels == "Ambiguous"] = "Ambiguous"
a.obs["rule_assignment_status"] = pd.Categorical(
    status,
    categories=["Assigned", "Ambiguous", "Unassigned"],
)

# -----------------------------------------------------------------------------
# Unified UMAP and spatial visualization
# -----------------------------------------------------------------------------
plot_umap_annotation(
    a,
    "cell_type_rule",
    FIG / "umap_cell_type_rule.png",
    title="Approach II — rule-based cell typing",
)

if "cell_type_scanpy" in a.obs.columns:
    plot_umap_side_by_side(
        a,
        "cell_type_scanpy",
        "cell_type_rule",
        FIG / "umap_approach_i_vs_ii.png",
    )
else:
    print("WARNING: cell_type_scanpy not found; skipping side-by-side UMAP.")

plot_umap_annotation(
    a,
    "rule_assignment_status",
    FIG / "umap_rule_assignment_status.png",
    title=f"Approach II — assignment status (>= {MIN_POSITIVE_GENES} markers > 0)",
)

# Two complementary spatial views for Approach II:
# (1) cell-center coordinates; (2) actual segmentation polygons.
plot_spatial_annotation_pair(
    a,
    "cell_type_rule",
    SDATA,
    FIG / "spatial",
    prefix="approach_ii",
    title="Approach II — rule-based",
    point_size=8,
)

plot_spatial_annotation_pair(
    a,
    "rule_assignment_status",
    SDATA,
    FIG / "spatial",
    prefix="approach_ii_assignment_status",
    title="Approach II — assignment status",
    point_size=8,
)

# -----------------------------------------------------------------------------
# Save Approach I + II object
# -----------------------------------------------------------------------------
OUTPUT_H5AD = OUT / "GSM9046088_FOV46_approach_i_ii.h5ad"
a.write_h5ad(OUTPUT_H5AD, compression="gzip")

print()
print("=" * 70)
print("APPROACH II SUMMARY")
print("=" * 70)
print(label_counts.to_string())
print()
print("Marker-panel availability:")
print(availability_df.to_string(index=False))
print()
print(f"Saved h5ad: {OUTPUT_H5AD}")
print(f"Figures: {FIG}")
print("Approach II complete.")
