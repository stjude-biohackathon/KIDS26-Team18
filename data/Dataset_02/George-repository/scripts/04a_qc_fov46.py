#!/usr/bin/env python

from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from src.qc import (
    add_qc_metrics,
    recommend_thresholds,
    apply_qc,
    threshold_sweep,
)


# ============================================================
# Paths
# ============================================================

INPUT = (
    ROOT
    / "data"
    / "processed"
    / "GSM9046088_CosMx_raw_FOV46.h5ad"
)

OUT = ROOT / "results" / "FOV46"

TAB = OUT / "tables"
FIG = OUT / "figures" / "qc"

OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)


# ============================================================
# Configuration
# ============================================================

N_MADS = 3.0

print("=" * 70)
print("FOV46 QC")
print("=" * 70)

print(f"Input:\n{INPUT}")
print(f"Output:\n{OUT}")
print(f"MAD threshold: {N_MADS}")
print()


# ============================================================
# 1. Load RAW FOV46
# ============================================================

adata = sc.read_h5ad(INPUT)

print("Raw AnnData:")
print(adata)
print()

n_raw = adata.n_obs


# ============================================================
# 2. Calculate QC metrics
# ============================================================

add_qc_metrics(adata)

print("Available QC-related columns:")

for col in [
    "total_counts",
    "n_genes_by_counts",
    "control_counts",
    "control_fraction",
    "Area",
]:
    if col in adata.obs:
        print(f"  {col}")

print()


# ============================================================
# 3. Threshold sensitivity sweep
# ============================================================

print("Running QC threshold sweep...")

sweep = threshold_sweep(adata)

sweep_path = TAB / "qc_threshold_sweep.csv"

sweep.to_csv(
    sweep_path,
    index=False,
)

print(f"Saved:\n{sweep_path}")
print()


# ============================================================
# 4. Data-driven MAD thresholds
# ============================================================

thresholds = recommend_thresholds(
    adata,
    nmads=N_MADS,
)

threshold_path = (
    TAB
    / "recommended_qc_thresholds.json"
)

threshold_path.write_text(
    json.dumps(
        thresholds,
        indent=2,
    )
)

print("=" * 70)
print("QC THRESHOLDS")
print("=" * 70)

for key, value in thresholds.items():

    if value is None:
        print(f"{key:25s}: None")

    elif isinstance(value, (float, np.floating)):
        print(
            f"{key:25s}: "
            f"{float(value):.6g}"
        )

    else:
        print(
            f"{key:25s}: "
            f"{value}"
        )

print()


# ============================================================
# 5. Apply QC
# ============================================================

apply_qc(
    adata,
    thresholds,
    use_per_fov_flags=False,
)

if "qc_pass" not in adata.obs:
    raise RuntimeError(
        "apply_qc() did not create "
        "adata.obs['qc_pass']."
    )

adata.obs["qc_pass"] = (
    adata.obs["qc_pass"]
    .astype(bool)
)


# ============================================================
# 6. Explicit QC failure flags
#
# These make the QC decision interpretable.
# ============================================================

obs = adata.obs

# Start all flags as False
failure_flags = {}


# ------------------------------------------------------------
# Low total counts
# ------------------------------------------------------------

if (
    "total_counts" in obs
    and thresholds.get("min_counts")
    is not None
):

    failure_flags["fail_low_counts"] = (
        obs["total_counts"]
        < thresholds["min_counts"]
    )


# ------------------------------------------------------------
# High total counts
# ------------------------------------------------------------

if (
    "total_counts" in obs
    and thresholds.get("max_counts")
    is not None
):

    failure_flags["fail_high_counts"] = (
        obs["total_counts"]
        > thresholds["max_counts"]
    )


# ------------------------------------------------------------
# Low detected genes
# ------------------------------------------------------------

if (
    "n_genes_by_counts" in obs
    and thresholds.get("min_genes")
    is not None
):

    failure_flags["fail_low_genes"] = (
        obs["n_genes_by_counts"]
        < thresholds["min_genes"]
    )


# ------------------------------------------------------------
# High control fraction
# ------------------------------------------------------------

if (
    "control_fraction" in obs
    and thresholds.get(
        "max_control_fraction"
    )
    is not None
):

    failure_flags[
        "fail_high_control_fraction"
    ] = (
        obs["control_fraction"]
        > thresholds[
            "max_control_fraction"
        ]
    )


# ------------------------------------------------------------
# Small cell area
# ------------------------------------------------------------

if (
    "Area" in obs
    and thresholds.get(
        "min_cell_area"
    )
    is not None
):

    failure_flags["fail_small_area"] = (
        obs["Area"]
        < thresholds["min_cell_area"]
    )


# ------------------------------------------------------------
# Large cell area
# ------------------------------------------------------------

if (
    "Area" in obs
    and thresholds.get(
        "max_cell_area"
    )
    is not None
):

    failure_flags["fail_large_area"] = (
        obs["Area"]
        > thresholds["max_cell_area"]
    )


# Add flags to obs
for name, flag in failure_flags.items():

    adata.obs[name] = (
        np.asarray(flag)
        .astype(bool)
    )


# ============================================================
# 7. QC summary
# ============================================================

n_pass = int(
    adata.obs["qc_pass"].sum()
)

n_fail = n_raw - n_pass

summary_rows = [
    {
        "criterion": "all_cells",
        "n_cells": n_raw,
        "fraction": 1.0,
    },
    {
        "criterion": "qc_pass",
        "n_cells": n_pass,
        "fraction": n_pass / n_raw,
    },
    {
        "criterion": "qc_fail",
        "n_cells": n_fail,
        "fraction": n_fail / n_raw,
    },
]


for flag_name in failure_flags:

    n = int(
        adata.obs[flag_name].sum()
    )

    summary_rows.append(
        {
            "criterion": flag_name,
            "n_cells": n,
            "fraction": n / n_raw,
        }
    )


qc_summary = pd.DataFrame(
    summary_rows
)

qc_summary.to_csv(
    TAB / "qc_filter_summary.csv",
    index=False,
)


print("=" * 70)
print("QC RETENTION")
print("=" * 70)

print(
    f"Raw cells : {n_raw:,}"
)

print(
    f"QC pass   : "
    f"{n_pass:,} "
    f"({n_pass / n_raw:.1%})"
)

print(
    f"QC fail   : "
    f"{n_fail:,} "
    f"({n_fail / n_raw:.1%})"
)

print()

print("Individual failure criteria:")

for flag_name in failure_flags:

    n = int(
        adata.obs[flag_name].sum()
    )

    print(
        f"{flag_name:30s} "
        f"{n:6,d} "
        f"({n / n_raw:6.2%})"
    )

print()

print(
    "NOTE: failure criteria overlap; "
    "their counts should NOT be summed."
)


# ============================================================
# 8. Number of failed criteria per cell
# ============================================================

flag_columns = list(
    failure_flags.keys()
)

if flag_columns:

    adata.obs[
        "n_qc_failures"
    ] = (
        adata.obs[
            flag_columns
        ]
        .astype(int)
        .sum(axis=1)
    )

else:

    adata.obs[
        "n_qc_failures"
    ] = 0


failure_count_table = (
    adata.obs[
        "n_qc_failures"
    ]
    .value_counts()
    .sort_index()
    .rename_axis(
        "n_qc_failures"
    )
    .reset_index(
        name="n_cells"
    )
)

failure_count_table[
    "fraction"
] = (
    failure_count_table[
        "n_cells"
    ]
    / n_raw
)

failure_count_table.to_csv(
    TAB
    / "qc_number_of_failures.csv",
    index=False,
)


# ============================================================
# 9. Failure combinations
# ============================================================

if flag_columns:

    combination_df = (
        adata.obs[
            flag_columns
        ]
        .astype(int)
        .groupby(
            flag_columns,
            observed=True,
        )
        .size()
        .reset_index(
            name="n_cells"
        )
        .sort_values(
            "n_cells",
            ascending=False,
        )
    )

    combination_df[
        "fraction"
    ] = (
        combination_df[
            "n_cells"
        ]
        / n_raw
    )

    combination_df.to_csv(
        TAB
        / "qc_failure_combinations.csv",
        index=False,
    )


# ============================================================
# 10. Before vs after numerical summary
# ============================================================

metrics = [
    x
    for x in [
        "total_counts",
        "n_genes_by_counts",
        "control_fraction",
        "Area",
    ]
    if x in adata.obs
]


summary = []

for metric in metrics:

    before = (
        adata.obs[
            metric
        ]
        .dropna()
        .astype(float)
    )

    after = (
        adata.obs.loc[
            adata.obs[
                "qc_pass"
            ],
            metric,
        ]
        .dropna()
        .astype(float)
    )

    for stage, values in [
        ("before_qc", before),
        ("after_qc", after),
    ]:

        if len(values) == 0:
            continue

        summary.append(
            {
                "metric": metric,
                "stage": stage,
                "n": len(values),
                "mean": values.mean(),
                "median": values.median(),
                "std": values.std(),
                "q01": values.quantile(
                    0.01
                ),
                "q05": values.quantile(
                    0.05
                ),
                "q25": values.quantile(
                    0.25
                ),
                "q75": values.quantile(
                    0.75
                ),
                "q95": values.quantile(
                    0.95
                ),
                "q99": values.quantile(
                    0.99
                ),
                "min": values.min(),
                "max": values.max(),
            }
        )


pd.DataFrame(
    summary
).to_csv(
    TAB
    / "qc_before_after_summary.csv",
    index=False,
)


# ============================================================
# Plot helper
# ============================================================

def save_histogram(
    values,
    output,
    title,
    xlabel,
    lower=None,
    upper=None,
):

    values = (
        pd.Series(values)
        .replace(
            [np.inf, -np.inf],
            np.nan,
        )
        .dropna()
        .astype(float)
    )

    if len(values) == 0:
        return

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.hist(
        values,
        bins=60,
        alpha=0.8,
    )

    if lower is not None:

        ax.axvline(
            lower,
            linestyle="--",
            linewidth=2,
            label=f"lower = {lower:.4g}",
        )

    if upper is not None:

        ax.axvline(
            upper,
            linestyle="--",
            linewidth=2,
            label=f"upper = {upper:.4g}",
        )

    if (
        lower is not None
        or upper is not None
    ):

        ax.legend(
            frameon=False
        )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of cells")

    fig.tight_layout()

    fig.savefig(
        output,
        dpi=250,
        bbox_inches="tight",
    )

    plt.close(fig)


# ============================================================
# 11. Raw distributions with thresholds
# ============================================================

if "total_counts" in adata.obs:

    save_histogram(
        adata.obs[
            "total_counts"
        ],
        FIG
        / "01_total_counts_raw.png",
        "FOV46 total counts before QC",
        "Total transcript counts",
        lower=thresholds.get(
            "min_counts"
        ),
        upper=thresholds.get(
            "max_counts"
        ),
    )


if "n_genes_by_counts" in adata.obs:

    save_histogram(
        adata.obs[
            "n_genes_by_counts"
        ],
        FIG
        / "02_n_genes_raw.png",
        "FOV46 detected genes before QC",
        "Number of detected genes",
        lower=thresholds.get(
            "min_genes"
        ),
    )


if "control_fraction" in adata.obs:

    save_histogram(
        adata.obs[
            "control_fraction"
        ],
        FIG
        / "03_control_fraction_raw.png",
        "FOV46 control fraction before QC",
        "Control fraction",
        upper=thresholds.get(
            "max_control_fraction"
        ),
    )


if "Area" in adata.obs:

    save_histogram(
        adata.obs[
            "Area"
        ],
        FIG
        / "04_cell_area_raw.png",
        "FOV46 cell area before QC",
        "Cell area",
        lower=thresholds.get(
            "min_cell_area"
        ),
        upper=thresholds.get(
            "max_cell_area"
        ),
    )


# ============================================================
# 12. Before vs after histograms
# ============================================================

for i, metric in enumerate(
    metrics,
    start=5,
):

    before = (
        adata.obs[
            metric
        ]
        .replace(
            [np.inf, -np.inf],
            np.nan,
        )
        .dropna()
        .astype(float)
    )

    after = (
        adata.obs.loc[
            adata.obs[
                "qc_pass"
            ],
            metric,
        ]
        .replace(
            [np.inf, -np.inf],
            np.nan,
        )
        .dropna()
        .astype(float)
    )

    if (
        len(before) == 0
        or len(after) == 0
    ):
        continue

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    bins = np.histogram_bin_edges(
        before,
        bins=60,
    )

    ax.hist(
        before,
        bins=bins,
        density=True,
        alpha=0.45,
        label="Before QC",
    )

    ax.hist(
        after,
        bins=bins,
        density=True,
        alpha=0.45,
        label="After QC",
    )

    ax.set_xlabel(metric)
    ax.set_ylabel("Density")

    ax.set_title(
        f"{metric}: before vs after QC"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIG
        / f"{i:02d}_{metric}_before_after.png",
        dpi=250,
        bbox_inches="tight",
    )

    plt.close(fig)


# ============================================================
# 13. Counts vs genes scatter
# ============================================================

if (
    "total_counts" in adata.obs
    and
    "n_genes_by_counts"
    in adata.obs
):

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    passed = (
        adata.obs[
            "qc_pass"
        ]
        .to_numpy()
    )

    ax.scatter(
        adata.obs.loc[
            ~passed,
            "total_counts",
        ],
        adata.obs.loc[
            ~passed,
            "n_genes_by_counts",
        ],
        s=8,
        alpha=0.5,
        label="QC fail",
    )

    ax.scatter(
        adata.obs.loc[
            passed,
            "total_counts",
        ],
        adata.obs.loc[
            passed,
            "n_genes_by_counts",
        ],
        s=8,
        alpha=0.5,
        label="QC pass",
    )

    ax.axvline(
        thresholds[
            "min_counts"
        ],
        linestyle="--",
        linewidth=1,
    )

    if (
        thresholds.get(
            "max_counts"
        )
        is not None
    ):

        ax.axvline(
            thresholds[
                "max_counts"
            ],
            linestyle="--",
            linewidth=1,
        )

    ax.axhline(
        thresholds[
            "min_genes"
        ],
        linestyle="--",
        linewidth=1,
    )

    ax.set_xlabel(
        "Total counts"
    )

    ax.set_ylabel(
        "Detected genes"
    )

    ax.set_title(
        "FOV46 RNA complexity and QC status"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIG
        / "09_counts_vs_genes_qc.png",
        dpi=250,
        bbox_inches="tight",
    )

    plt.close(fig)


# ============================================================
# 14. Spatial QC visualization
# ============================================================

def find_spatial_coordinates(
    adata,
):

    if "spatial" in adata.obsm:

        coords = np.asarray(
            adata.obsm[
                "spatial"
            ]
        )

        if (
            coords.ndim == 2
            and coords.shape[1] >= 2
        ):
            return (
                coords[:, 0],
                coords[:, 1],
            )

    candidates = [
        (
            "CenterX_global_px",
            "CenterY_global_px",
        ),
        (
            "coord_x",
            "coord_y",
        ),
        (
            "x_global_px",
            "y_global_px",
        ),
    ]

    for x_col, y_col in candidates:

        if (
            x_col in adata.obs
            and y_col in adata.obs
        ):

            return (
                adata.obs[
                    x_col
                ].to_numpy(),
                adata.obs[
                    y_col
                ].to_numpy(),
            )

    return None


coords = find_spatial_coordinates(
    adata
)

if coords is not None:

    x, y = coords

    # --------------------------------------------------------
    # QC pass / fail spatial map
    # --------------------------------------------------------

    passed = (
        adata.obs[
            "qc_pass"
        ]
        .to_numpy()
    )

    fig, ax = plt.subplots(
        figsize=(8, 8)
    )

    ax.scatter(
        x[passed],
        y[passed],
        s=8,
        alpha=0.5,
        label="QC pass",
    )

    ax.scatter(
        x[~passed],
        y[~passed],
        s=14,
        alpha=0.8,
        label="QC fail",
    )

    ax.set_aspect("equal")

    ax.invert_yaxis()

    ax.set_title(
        "FOV46 spatial QC"
    )

    ax.set_xlabel(
        "X"
    )

    ax.set_ylabel(
        "Y"
    )

    ax.legend(
        frameon=False
    )

    fig.tight_layout()

    fig.savefig(
        FIG
        / "10_spatial_qc_pass_fail.png",
        dpi=250,
        bbox_inches="tight",
    )

    plt.close(fig)


    # --------------------------------------------------------
    # Individual failure reasons
    # --------------------------------------------------------

    for i, flag_name in enumerate(
        flag_columns,
        start=11,
    ):

        failed = (
            adata.obs[
                flag_name
            ]
            .to_numpy()
        )

        if failed.sum() == 0:
            continue

        fig, ax = plt.subplots(
            figsize=(8, 8)
        )

        ax.scatter(
            x[~failed],
            y[~failed],
            s=6,
            alpha=0.25,
            label="Other cells",
        )

        ax.scatter(
            x[failed],
            y[failed],
            s=16,
            alpha=0.85,
            label=flag_name,
        )

        ax.set_aspect(
            "equal"
        )

        ax.invert_yaxis()

        ax.set_title(
            f"FOV46: {flag_name}"
        )

        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        ax.legend(
            frameon=False
        )

        fig.tight_layout()

        fig.savefig(
            FIG
            / f"{i:02d}_spatial_{flag_name}.png",
            dpi=250,
            bbox_inches="tight",
        )

        plt.close(fig)


# ============================================================
# 15. Save RAW object with QC annotations
# ============================================================

raw_qc_path = (
    OUT
    / "GSM9046088_FOV46_raw_with_qc_flags.h5ad"
)

adata.write_h5ad(
    raw_qc_path,
    compression="gzip",
)


# ============================================================
# 16. Create QC-filtered biological AnnData
# ============================================================

if "is_control" in adata.var:

    biological_mask = (
        ~adata.var[
            "is_control"
        ]
        .to_numpy()
    )

else:

    biological_mask = np.ones(
        adata.n_vars,
        dtype=bool,
    )


qc_mask = (
    adata.obs[
        "qc_pass"
    ]
    .to_numpy()
)


q = adata[
    qc_mask,
    biological_mask,
].copy()


# Remove genes absent after QC
sc.pp.filter_genes(
    q,
    min_cells=1,
)


# Preserve raw biological counts
q.layers[
    "counts"
] = q.X.copy()


filtered_path = (
    OUT
    / "GSM9046088_FOV46_qc_filtered.h5ad"
)

q.write_h5ad(
    filtered_path,
    compression="gzip",
)


# ============================================================
# Final report
# ============================================================

print()
print("=" * 70)
print("FINAL QC OUTPUT")
print("=" * 70)

print(
    f"Raw FOV46: "
    f"{adata.n_obs:,} cells x "
    f"{adata.n_vars:,} features"
)

print(
    f"QC-passed biological object: "
    f"{q.n_obs:,} cells x "
    f"{q.n_vars:,} genes"
)

print(
    f"Cell retention: "
    f"{q.n_obs / adata.n_obs:.2%}"
)

print()

print(
    "Raw + QC flags:\n"
    f"{raw_qc_path}"
)

print()

print(
    "QC-filtered AnnData:\n"
    f"{filtered_path}"
)

print()

print(
    "QC tables:\n"
    f"{TAB}"
)

print()

print(
    "QC figures:\n"
    f"{FIG}"
)

print()
print(
    "QC stage complete."
)
print(
    "Inspect QC figures and tables "
    "before interpreting Leiden results."
)