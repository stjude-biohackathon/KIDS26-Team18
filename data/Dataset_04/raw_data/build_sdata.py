#!/usr/bin/env python3
"""Build SpatialData + AnnData from CosMx GEO flat files under raw_data/."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from anndata import AnnData
from scipy.sparse import csr_matrix
from spatialdata import SpatialData
from spatialdata.models import PointsModel, ShapesModel, TableModel
from spatialdata.transformations import Identity

RAW_DIR = Path(__file__).resolve().parent
FLAT_DIR = RAW_DIR / "GSE282026_RAW"
SDATA_PATH = RAW_DIR / "sdata.zarr"
H5AD_PATH = RAW_DIR / "adata.h5ad"

COUNT_INDEX_KEY = "cell_ID"
INSTANCE_KEY = "cell_uid"
REGION_KEY = "region"
GLOBAL_REGION = "cell_boundaries"
X_GLOBAL = "CenterX_global_px"
Y_GLOBAL = "CenterY_global_px"
X_LOCAL = "CenterX_local_px"
Y_LOCAL = "CenterY_local_px"
SHAPES_KEY = "cell_boundaries"
POINTS_KEY = "transcripts"
TX_FEATURE_KEY = "target"

_FLAT_PATTERNS = {
    "exprMat": "*_exprMat_file.csv.gz",
    "metadata": "*_metadata_file.csv.gz",
    "polygons": "*-polygons.csv.gz",
    "tx": "*_tx_file.csv.gz",
}


def _find_flat(file_type: str) -> Path:
    pattern = _FLAT_PATTERNS[file_type]
    matches = sorted(FLAT_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No {file_type!r} file matching {pattern!r} in {FLAT_DIR}")
    return matches[0]


def build_table() -> AnnData:
    counts_path = _find_flat("exprMat")
    meta_path = _find_flat("metadata")
    print(f"Reading counts: {counts_path.name}")
    counts = pd.read_csv(counts_path, header=0, index_col=COUNT_INDEX_KEY)
    fov = counts.pop("fov").astype(str)
    counts.index = counts.index.astype(str).str.cat(fov.values, sep="_")

    print(f"Reading metadata: {meta_path.name}")
    obs = pd.read_csv(meta_path, header=0, index_col=COUNT_INDEX_KEY)
    obs["fov"] = obs["fov"].astype(str)
    obs[REGION_KEY] = GLOBAL_REGION
    obs[INSTANCE_KEY] = obs.index.astype(str).str.cat(obs["fov"].values, sep="_")
    obs.index = obs[INSTANCE_KEY]

    if X_GLOBAL not in obs.columns or Y_GLOBAL not in obs.columns:
        raise ValueError(f"Metadata must contain {X_GLOBAL} and {Y_GLOBAL}")

    common = obs.index.intersection(counts.index)
    if len(common) < len(obs):
        print(f"Warning: {len(obs) - len(common)} metadata rows missing from counts.")

    adata = AnnData(
        csr_matrix(counts.loc[common, :].values),
        dtype=counts.values.dtype,
        obs=obs.loc[common, :],
    )
    adata.var_names = counts.columns
    adata.obsm["global"] = adata.obs[[X_GLOBAL, Y_GLOBAL]].to_numpy(dtype=np.float64)
    adata.obsm["spatial"] = adata.obs[[X_LOCAL, Y_LOCAL]].to_numpy(dtype=np.float64)
    adata.uns["dataset_id"] = counts_path.name.split("_exprMat_file")[0]
    print(f"  table: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


def build_polygons() -> ShapesModel:
    import geopandas as gpd
    from shapely.geometry import Polygon

    poly_path = _find_flat("polygons")
    print(f"Reading polygons: {poly_path.name}")
    df = pd.read_csv(poly_path)
    df.columns = [str(c).strip() for c in df.columns]
    cid_col = "cellID" if "cellID" in df.columns else "cell_ID"
    if cid_col not in df.columns or "fov" not in df.columns:
        raise ValueError(f"Unexpected polygon columns: {list(df.columns)}")

    geometries: list[Polygon] = []
    index: list[str] = []
    for (fov, cell_id), grp in df.groupby(["fov", cid_col], sort=False):
        xy = grp[["x_global_px", "y_global_px"]].to_numpy(dtype=np.float64)
        if len(xy) < 3:
            continue
        poly = Polygon(xy)
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty:
            continue
        geometries.append(poly)
        index.append(f"{cell_id}_{int(fov)}")

    print(f"  built {len(geometries):,} cell polygons")
    gdf = gpd.GeoDataFrame(geometry=geometries, index=pd.Index(index, name=INSTANCE_KEY))
    return ShapesModel.parse(gdf, transformations={"global": Identity()})


def build_transcripts(cache_parquet: Path, chunksize: int = 2_000_000):
    import dask.dataframe as dd

    tx_path = _find_flat("tx")
    if not cache_parquet.exists():
        print(f"Reading transcripts (chunked) from {tx_path.name} ...")
        writer = None
        n_rows = 0
        for chunk in pd.read_csv(tx_path, chunksize=chunksize):
            chunk.columns = [str(c).strip() for c in chunk.columns]
            out = pd.DataFrame(
                {
                    "x": chunk["x_global_px"].astype(np.float32),
                    "y": chunk["y_global_px"].astype(np.float32),
                    TX_FEATURE_KEY: chunk[TX_FEATURE_KEY].astype(str),
                    INSTANCE_KEY: chunk["cell_ID"].astype(str).str.cat(
                        chunk["fov"].astype(str), sep="_"
                    ),
                    "fov": chunk["fov"].astype(np.int32),
                }
            )
            table = pa.Table.from_pandas(out, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(cache_parquet, table.schema, compression="zstd")
            writer.write_table(table)
            n_rows += len(out)
            print(f"  ... {n_rows:,} transcripts", flush=True)
        if writer is not None:
            writer.close()
        print(f"  cached {n_rows:,} transcripts -> {cache_parquet.name}")
    else:
        print(f"Using cached transcripts: {cache_parquet.name}")

    ddf = dd.read_parquet(cache_parquet)
    return PointsModel.parse(
        ddf,
        coordinates={"x": "x", "y": "y"},
        feature_key=TX_FEATURE_KEY,
        instance_key=INSTANCE_KEY,
        transformations={"global": Identity()},
    )


def build_sdata(include_transcripts: bool = True, force: bool = False) -> None:
    if SDATA_PATH.exists() and not force:
        print(f"Skipping build — {SDATA_PATH} already exists (use --force to rebuild)")
        return

    adata = build_table()
    table = TableModel.parse(
        adata,
        region=GLOBAL_REGION,
        region_key=REGION_KEY,
        instance_key=INSTANCE_KEY,
    )
    shapes = {SHAPES_KEY: build_polygons()}
    points = {}
    if include_transcripts:
        cache = FLAT_DIR / "transcripts.global_coords.parquet"
        points[POINTS_KEY] = build_transcripts(cache)

    sdata = SpatialData(tables={"table": table}, shapes=shapes, points=points or None)
    sdata.attrs["dataset_id"] = adata.uns.get("dataset_id", "")
    sdata.attrs["coordinate_system"] = "global"

    if SDATA_PATH.exists():
        shutil.rmtree(SDATA_PATH)
    print(f"Writing {SDATA_PATH} ...")
    sdata.write(SDATA_PATH)

    print(f"Writing {H5AD_PATH} ({adata.n_obs:,} cells x {adata.n_vars:,} genes) ...")
    adata.write_h5ad(H5AD_PATH)


def validate() -> None:
    import spatialdata as sd
    from anndata import read_h5ad

    sdata = sd.read_zarr(SDATA_PATH)
    assert "table" in sdata.tables, "Missing table in SpatialData"
    adata = sdata["table"]
    h5ad = read_h5ad(H5AD_PATH)
    print(f"SpatialData elements: tables={list(sdata.tables)}, shapes={list(sdata.shapes)}, points={list(sdata.points)}")
    print(f"AnnData (zarr): {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"AnnData (h5ad): {h5ad.n_obs:,} cells x {h5ad.n_vars:,} genes")
    print(f"obsm keys: {list(adata.obsm.keys())}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Rebuild outputs")
    parser.add_argument("--skip-transcripts", action="store_true", help="Skip transcript points")
    parser.add_argument("--validate-only", action="store_true", help="Only run validation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.validate_only:
        validate()
        return

    build_sdata(include_transcripts=not args.skip_transcripts, force=args.force)
    validate()
    print("Done.")


if __name__ == "__main__":
    main()
