#!/usr/bin/env python3
"""Extract TMA5_subset.zarr.tar.gz and load the SpatialData store.

Run with the spatialdata conda env:

    /mnt/scratch1/miniconda3/envs/spatialdata/bin/python load_tma5_subset_zarr.py

Optional flags:

    --force          Re-extract even if TMA5_subset.zarr/ already exists
    --tar PATH       Path to the .tar.gz (default: TMA5_subset.zarr.tar.gz)
    --out PATH       Extraction directory (default: TMA5_subset.zarr)
"""

from __future__ import annotations

import argparse
import shutil
import tarfile
from pathlib import Path

import numpy as np
import spatialdata as sd

PROC_DIR = Path(__file__).resolve().parent
DEFAULT_TAR = PROC_DIR / "TMA5_subset.zarr.tar.gz"
DEFAULT_OUT = PROC_DIR / "TMA5_subset.zarr"


def extract_zarr(tar_path: Path, out_dir: Path, *, force: bool = False) -> Path:
    if not tar_path.is_file():
        raise FileNotFoundError(f"Archive not found: {tar_path}")

    if out_dir.exists():
        if force:
            shutil.rmtree(out_dir)
        else:
            print(f"Using existing zarr store: {out_dir}")
            return out_dir

    print(f"Extracting {tar_path.name} -> {out_dir}")
    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(out_dir.parent)
    if not out_dir.is_dir():
        raise RuntimeError(f"Expected zarr store after extract: {out_dir}")
    return out_dir


def summarize_sdata(sdata: sd.SpatialData) -> None:
    print(f"Images:  {list(sdata.images.keys())}")
    print(f"Labels:  {list(sdata.labels.keys())}")
    print(f"Shapes:  {list(sdata.shapes.keys())}")
    print(f"Points:  {list(sdata.points.keys())}")
    print(f"Tables:  {list(sdata.tables.keys())}")
    print(f"Coordinate systems: {list(sdata.coordinate_systems)}")
    if "table" in sdata.tables:
        t = sdata.tables["table"]
        print(f"\nTable: {t.n_obs:,} cells x {t.n_vars} genes")
        region = t.uns.get("spatialdata_attrs", {}).get("region")
        print(f"  region: {region}")
        if "spatial" in t.obsm:
            coords = np.asarray(t.obsm["spatial"])
            print(
                f"  obsm['spatial']: x=[{coords[:, 0].min():.0f}, {coords[:, 0].max():.0f}], "
                f"y=[{coords[:, 1].min():.0f}, {coords[:, 1].max():.0f}]"
            )
    for key in sdata.shapes:
        print(f"  shapes[{key!r}]: {len(sdata.shapes[key]):,}")
    if "transcripts" in sdata.points:
        print(f"  points['transcripts']: {len(sdata.points['transcripts']):,}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract and load TMA5_subset.zarr")
    parser.add_argument("--tar", type=Path, default=DEFAULT_TAR, help="Input .tar.gz")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output zarr directory")
    parser.add_argument("--force", action="store_true", help="Re-extract even if output exists")
    args = parser.parse_args()

    zarr_path = extract_zarr(args.tar, args.out, force=args.force)
    print(f"\nLoading SpatialData from {zarr_path} ...")
    sdata = sd.read_zarr(zarr_path)
    summarize_sdata(sdata)


if __name__ == "__main__":
    main()
