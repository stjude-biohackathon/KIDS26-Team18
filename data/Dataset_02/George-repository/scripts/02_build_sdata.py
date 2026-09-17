#!/usr/bin/env python
"""Build SpatialData from the existing raw AnnData and GEO CosMx spatial files.

Inputs
------
- data/processed/GSM9046088_CosMx_raw.h5ad
- polygon CSV from GEO
- transcript CSV from GEO

Output
------
- data/spatial/GSM9046088_CosMx.zarr

The five GEO files do not contain a morphology image, so this SpatialData
contains the cell table, cell-boundary polygons, and transcript points.
Global pixel coordinates are used throughout.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import dask.dataframe as dd
import geopandas as gpd
import pandas as pd
import scanpy as sc
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import PointsModel, ShapesModel, TableModel
from spatialdata.transformations import Identity

from src.cosmx_io import FILES, unique_cell_id

RAW_DIR = ROOT / "data" / "raw"
ADATA_PATH = ROOT / "data" / "processed" / "GSM9046088_CosMx_raw.h5ad"
OUT_ZARR = ROOT / "data" / "spatial" / "GSM9046088_CosMx.zarr"


def build_shapes():
    polygon_path = RAW_DIR / FILES["polygons"]
    print(f"\nReading polygons:\n{polygon_path}")

    poly = pd.read_csv(polygon_path)
    cell_col = "cell_ID" if "cell_ID" in poly.columns else "cellID"

    required = {"fov", cell_col, "x_global_px", "y_global_px"}
    missing = required.difference(poly.columns)
    if missing:
        raise ValueError(f"Polygon file missing required columns: {sorted(missing)}")

    poly["unique_cell_id"] = unique_cell_id(
        pd.to_numeric(poly["fov"]).astype(int),
        pd.to_numeric(poly[cell_col]).astype(int),
    )

    records = []
    invalid_fixed = 0
    skipped = 0

    for uid, group in poly.groupby("unique_cell_id", sort=False):
        xy = group[["x_global_px", "y_global_px"]].to_numpy(dtype=float)

        if len(xy) < 3:
            skipped += 1
            continue

        geom = Polygon(xy)

        if not geom.is_valid:
            geom = geom.buffer(0)
            invalid_fixed += 1

        if geom.is_empty:
            skipped += 1
            continue

        records.append((uid, geom))

    if not records:
        raise ValueError("No valid cell polygons could be constructed.")

    shapes = gpd.GeoDataFrame(
        {"geometry": [geom for _, geom in records]},
        index=pd.Index([uid for uid, _ in records], name="instance_id"),
        crs=None,
    )

    if not shapes.index.is_unique:
        raise ValueError("Cell-boundary shape index is not unique.")

    print(f"Valid polygons : {len(shapes):,}")
    print(f"Invalid fixed  : {invalid_fixed:,}")
    print(f"Skipped        : {skipped:,}")

    shapes_model = ShapesModel.parse(
        shapes,
        transformations={"global": Identity()},
    )

    return shapes, shapes_model


def build_table(adata, shapes):
    adata = adata.copy()

    adata.obs["region"] = "cell_boundaries"
    adata.obs["region"] = adata.obs["region"].astype("category")
    adata.obs["instance_id"] = adata.obs_names.astype(str)

    keep = adata.obs_names.intersection(shapes.index)

    print(f"\nCells in raw AnnData : {adata.n_obs:,}")
    print(f"Cells with polygons   : {len(keep):,}")
    print(f"Cells removed         : {adata.n_obs - len(keep):,}")

    if len(keep) == 0:
        raise ValueError("No overlapping cells between AnnData and polygons.")

    adata = adata[keep].copy()

    table = TableModel.parse(
        adata,
        region="cell_boundaries",
        region_key="region",
        instance_key="instance_id",
    )

    return table


def build_transcript_points():
    transcript_path = RAW_DIR / FILES["transcripts"]
    print(f"\nReading transcript points:\n{transcript_path}")

    preview = pd.read_csv(transcript_path, nrows=100)
    cols = set(preview.columns)

    print("Transcript columns:")
    print(preview.columns.tolist())

    x_col = "x_global_px"
    y_col = "y_global_px"
    feature_col = "target"

    required = {x_col, y_col, feature_col}
    missing = required.difference(cols)
    if missing:
        raise ValueError(
            f"Transcript file missing required columns: {sorted(missing)}"
        )

    # GEO supplies this table as gzip. A single gzip stream cannot be split into
    # independent Dask blocks, so blocksize=None is intentional.
    # Explicit string dtypes prevent inference failures such as CellComp being
    # inferred as float from early NA rows and later containing 'Nuclear'.
    dtype_map = {feature_col: "object"}
    if "CellComp" in cols:
        dtype_map["CellComp"] = "object"

    tx = dd.read_csv(
        transcript_path,
        blocksize=None,
        assume_missing=True,
        dtype=dtype_map,
    )

    print("\nTranscript Dask dtypes:")
    print(tx.dtypes)

    # Do not assign an instance_key to the transcript PointsModel. Some CosMx
    # transcripts may be unassigned to a segmented cell but remain valid spatial
    # molecules and should be retained for transcript-level visualization.
    points = PointsModel.parse(
        tx,
        coordinates={"x": x_col, "y": y_col},
        feature_key=feature_col,
        transformations={"global": Identity()},
    )

    return points


def main() -> None:
    print("=" * 70)
    print("Building CosMx SpatialData")
    print(f"Project root : {ROOT}")
    print(f"Raw AnnData  : {ADATA_PATH}")
    print(f"Output Zarr  : {OUT_ZARR}")
    print("=" * 70)

    if not ADATA_PATH.is_file():
        raise FileNotFoundError(
            f"Raw AnnData not found: {ADATA_PATH}\n"
            "Run scripts/01_build_anndata.py first."
        )

    print(f"\nLoading existing raw AnnData:\n{ADATA_PATH}")
    adata = sc.read_h5ad(ADATA_PATH)
    print(adata)

    shapes, shapes_model = build_shapes()
    table = build_table(adata, shapes)
    points = build_transcript_points()

    sdata = SpatialData(
        points={"transcripts": points},
        shapes={"cell_boundaries": shapes_model},
        tables={"table": table},
    )

    print("\nSpatialData constructed successfully:")
    print(sdata)

    OUT_ZARR.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nWriting SpatialData:\n{OUT_ZARR}")
    sdata.write(OUT_ZARR, overwrite=True)

    print("\nSpatialData saved successfully.")
    print(OUT_ZARR)


if __name__ == "__main__":
    main()
