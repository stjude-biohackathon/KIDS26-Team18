#!/usr/bin/env python3
"""Subset TMA5.zarr to a ~1 GB dense-region SpatialData store.

Single-pass crop of a cell-dense ROI, keeping all elements (images, labels,
shapes, transcripts, table). Run with:

    /mnt/scratch1/miniconda3/envs/spatialdata/bin/python subset_tma5_zarr.py --half 2200

Key fixes vs. the original (broken) version
-------------------------------------------
1. Coordinate systems: ``obsm['spatial']`` and shape geometries are stored in RAW
   pixel coords; the ``global`` coordinate system used by ``query.bounding_box`` is
   RAW * ~4.706 (a ``Scale`` transform on shapes/points; images/labels are Identity).
   The ROI center from ``cell_density_roi_center`` (raw coords) is therefore scaled to
   global before cropping. The old code passed raw coords into a global query and
   dropped the entire table.
2. Raster write: multiscale bbox crops cannot be written back directly in this env.
   Each image/label is collapsed to a coarse level and re-parsed into a small, uniformly
   chunked multiscale pyramid (keeps ``scale0`` + correct transform/alignment).
3. ome_zarr <-> zarr v3 chunk bug: ``_retuple`` is patched to flatten dask's
   tuple-of-tuples ``.chunks`` into the flat integer chunk shape zarr v3 requires.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import spatialdata as sd
from spatialdata.models import Image2DModel, Labels2DModel
from spatialdata.transformations import get_transformation

# --- Patch ome_zarr <-> zarr v3 chunk incompatibility (see module docstring, fix #3) ---
import ome_zarr.writer as _ozw


def _retuple_flat(chunks, shape):
    if isinstance(chunks, int):
        return tuple([chunks] * len(shape))
    flat = tuple(max(c) if isinstance(c, (tuple, list)) else c for c in chunks)
    dims_to_add = len(shape) - len(flat)
    return (*shape[:dims_to_add], *flat)


_ozw._retuple = _retuple_flat

PROC_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROC_DIR.parents[2]
SRC_DIR = REPO_ROOT / "src"
RAW_ZARR = PROC_DIR.parent / "raw_data" / "TMA5.zarr"
DEFAULT_OUT = PROC_DIR / "TMA5_subset.zarr"

sys.path.insert(0, str(SRC_DIR))
from spatial_plot import (  # noqa: E402
    _get_xy_scale,
    build_plot_crop,
    cell_density_roi_center,
)


def dir_size_gb(path: Path) -> float:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / 1e9


def collapse_rasters_to_single_scale(
    sdata: sd.SpatialData,
    image_scale: str = "scale2",
    label_scale: str = "scale2",
    chunk: int = 1024,
    n_pyramid_levels: int = 2,
) -> sd.SpatialData:
    """Replace each multiscale image/label with a compact multiscale pyramid built from
    one coarse level, preserving its transformation to ``global`` (keeps alignment).

    A full-resolution (scale0) bbox crop is huge and also fails to write; starting from a
    coarse level keeps the output small. We rebuild a short pyramid (``scale_factors``)
    with uniform chunks so the store writes cleanly and stays multiscale.
    """
    scale_factors = [2] * n_pyramid_levels

    for key in list(sdata.images.keys()):
        tree = sdata.images[key]
        if not hasattr(tree, "children"):
            continue
        avail = sorted(tree.children.keys(), key=lambda s: int(s.replace("scale", "")))
        scale = image_scale if image_scale in avail else avail[-1]
        arr = tree[scale].ds[list(tree[scale].ds.data_vars)[0]]
        transform = get_transformation(arr, "global")
        data = arr.data.rechunk((arr.sizes["c"], chunk, chunk))
        sdata.images[key] = Image2DModel.parse(
            data,
            dims=tuple(arr.dims),
            c_coords=list(arr.coords["c"].values),
            transformations={"global": transform},
            scale_factors=scale_factors,
            chunks=(arr.sizes["c"], chunk, chunk),
        )

    for key in list(sdata.labels.keys()):
        tree = sdata.labels[key]
        if not hasattr(tree, "children"):
            continue
        avail = sorted(tree.children.keys(), key=lambda s: int(s.replace("scale", "")))
        scale = label_scale if label_scale in avail else avail[-1]
        arr = tree[scale].ds[list(tree[scale].ds.data_vars)[0]]
        transform = get_transformation(arr, "global")
        data = arr.data.rechunk((chunk, chunk))
        sdata.labels[key] = Labels2DModel.parse(
            data,
            dims=tuple(arr.dims),
            transformations={"global": transform},
            scale_factors=scale_factors,
            chunks=(chunk, chunk),
        )
    return sdata


def summarize_sdata(sdata: sd.SpatialData) -> None:
    print(f"  Images:  {list(sdata.images.keys())}")
    print(f"  Labels:  {list(sdata.labels.keys())}")
    print(f"  Shapes:  {', '.join(f'{k}={len(sdata.shapes[k]):,}' for k in sdata.shapes)}")
    print(f"  Points:  {', '.join(f'{k}={len(sdata.points[k]):,}' for k in sdata.points)}")
    if "table" in sdata.tables:
        t = sdata.tables["table"]
        print(f"  Table:   {t.n_obs:,} cells x {t.n_vars} genes")


def main() -> None:
    parser = argparse.ArgumentParser(description="Subset TMA5.zarr to a dense ~1 GB SpatialData store.")
    parser.add_argument("--half", type=float, default=2650.0,
                        help="ROI half-width in RAW pixel units (default: 2650 ~= 0.9 GB)")
    parser.add_argument("--image-scale", type=str, default="scale2",
                        help="Source pyramid level used for images (default: scale2)")
    parser.add_argument("--label-scale", type=str, default="scale2",
                        help="Source pyramid level used for labels (default: scale2)")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output zarr path")
    parser.add_argument("--input", type=Path, default=RAW_ZARR, help="Input TMA5.zarr path")
    args = parser.parse_args()

    t0 = time.time()
    print("=" * 60, flush=True)
    print("subset_tma5_zarr.py — started", flush=True)
    print("=" * 60, flush=True)
    print(f"Input:  {args.input}  ({dir_size_gb(args.input):.2f} GB)", flush=True)
    print(f"Output: {args.out}", flush=True)

    print("\n[1/4] Loading SpatialData ...", flush=True)
    sdata = sd.read_zarr(args.input)
    coords = np.asarray(sdata.tables["table"].obsm["spatial"])
    region = sdata.tables["table"].uns["spatialdata_attrs"]["region"]
    scale = _get_xy_scale(sdata, region)  # raw -> global factor (~4.706)
    print(f"  {coords.shape[0]:,} cells | raw->global scale = {scale:.4f}", flush=True)

    print(f"\n[2/4] Cropping dense ROI (raw half={args.half:.0f}) ...", flush=True)
    cx_r, cy_r = cell_density_roi_center(coords, half=args.half)
    cx_g, cy_g, half_g = cx_r * scale, cy_r * scale, args.half * scale
    print(f"  center raw=({cx_r:.0f}, {cy_r:.0f}) -> global=({cx_g:.0f}, {cy_g:.0f}) half_g={half_g:.0f}", flush=True)

    sdata_sub = build_plot_crop(sdata, cx_g, cy_g, half_g, filter_table=True)
    sdata_sub = collapse_rasters_to_single_scale(
        sdata_sub, image_scale=args.image_scale, label_scale=args.label_scale
    )
    summarize_sdata(sdata_sub)

    print("\n[3/4] Writing ...", flush=True)
    tmp_dir = args.out.with_suffix(".zarr.tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    t_write = time.time()
    sdata_sub.write(tmp_dir)
    size_gb = dir_size_gb(tmp_dir)
    print(f"  write done in {time.time() - t_write:.1f}s  size={size_gb:.2f} GB", flush=True)

    if args.out.exists():
        shutil.rmtree(args.out)
    tmp_dir.rename(args.out)

    print(f"\n[4/4] Done in {(time.time() - t0) / 60:.1f} min", flush=True)
    print(f"Output: {args.out}  ({size_gb:.2f} GB)", flush=True)
    print(f"Cells:  {sdata_sub.tables['table'].n_obs:,}", flush=True)


if __name__ == "__main__":
    main()
