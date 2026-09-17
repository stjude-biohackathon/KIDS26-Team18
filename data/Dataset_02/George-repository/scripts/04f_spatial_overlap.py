#!/usr/bin/env python

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot  # noqa: F401

from src.spatial_plot_legacy import (
    plot_gene_seg_transcripts,
)


# ============================================================
# Paths
# ============================================================

SDATA_PATH = (
    ROOT
    / "data"
    / "spatial"
    / "GSM9046088_CosMx_FOV46.zarr"
)

OUT = (
    ROOT
    / "results"
    / "FOV46"
    / "figures"
    / "spatial_overlap"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Genes
# ============================================================

GENES = [
    "EPCAM",   # epithelial
    "SFTPC",   # AT2
    "CD3D",    # T cell
    "MS4A1",   # B cell
    "C1QA",    # macrophage
    "PECAM1",  # endothelial
    "COL1A1",  # fibroblast
    "RGS5",    # mural
]


# ============================================================
# Load SpatialData
# ============================================================

print("=" * 70)
print("FOV46 SPATIAL OVERLAY")
print("=" * 70)

print(
    f"Reading:\n{SDATA_PATH}"
)

S = sd.read_zarr(
    SDATA_PATH
)

print()
print(S)
print()


# ============================================================
# Validate structure
# ============================================================

required_shapes = "cell_boundaries"
required_points = "transcripts"
required_table = "table"

if required_shapes not in S.shapes:
    raise KeyError(
        f"Missing shape element: "
        f"{required_shapes}"
    )

if required_points not in S.points:
    raise KeyError(
        f"Missing points element: "
        f"{required_points}"
    )

if required_table not in S.tables:
    raise KeyError(
        f"Missing table element: "
        f"{required_table}"
    )


table = S[required_table]

print(
    "SpatialData attrs:"
)

print(
    table.uns.get(
        "spatialdata_attrs"
    )
)

print()


# ============================================================
# Plot
# ============================================================

succeeded = []
failed = []
missing = []


for gene in GENES:

    print()
    print("-" * 70)
    print(f"Gene: {gene}")
    print("-" * 70)

    if gene not in table.var_names:

        print(
            f"SKIP: {gene} is not "
            "present in CosMx panel."
        )

        missing.append(gene)

        continue

    try:

        plot_gene_seg_transcripts(
            S,
            gene,
            shapes_key="cell_boundaries",
            points_key="transcripts",
            transcript_feature_col="target",
            coordinate_system="global",
            point_size=1.2,
            point_alpha=0.65,
            palette="orange",
            figsize=(8, 8),
            points_render_method="matplotlib",
        )

        output = (
            OUT
            / f"sdata_overlay_{gene}.png"
        )

        plt.gcf().savefig(
            output,
            dpi=250,
            bbox_inches="tight",
        )

        plt.close("all")

        succeeded.append(gene)

        print(
            f"SUCCESS:\n{output}"
        )

    except Exception as exc:

        plt.close("all")

        failed.append(
            (
                gene,
                repr(exc),
            )
        )

        print(
            f"FAILED: {gene}\n"
            f"{exc}"
        )


# ============================================================
# Summary
# ============================================================

print()
print("=" * 70)
print("SPATIAL OVERLAY SUMMARY")
print("=" * 70)

print(
    f"Succeeded ({len(succeeded)}): "
    f"{succeeded}"
)

print(
    f"Missing ({len(missing)}): "
    f"{missing}"
)

print(
    f"Failed ({len(failed)}): "
    f"{failed}"
)


# ============================================================
# Fail job if actual plotting errors occurred
# ============================================================

if failed:

    raise RuntimeError(
        f"{len(failed)} spatial overlay(s) "
        "failed. See log above."
    )


print()
print(
    "Spatial overlay stage completed successfully."
)