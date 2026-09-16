#!/usr/bin/env python3
"""Subset TMA5.zarr to a ~1 GB dense-region SpatialData store.

Run with:
    /mnt/scratch1/miniconda3/envs/spatialdata/bin/python subset_tma5_zarr.py --target-gb 1.0
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import spatialdata as sd

PROC_DIR = Path(__file__).resolve().parent
REPO_ROOT = PROC_DIR.parents[2]
SRC_DIR = REPO_ROOT / "src"
RAW_ZARR = PROC_DIR.parent / "raw_data" / "TMA5.zarr"
DEFAULT_OUT = PROC_DIR / "TMA5_subset.zarr"

sys.path.insert(0, str(SRC_DIR))
from spatial_plot import build_plot_crop, cell_density_roi_center  # noqa: E402


def dir_size_gb(path: Path) -> float:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / 1e9


def keep_coarsest_scales(sdata: sd.SpatialData, n_keep: int = 2) -> sd.SpatialData:
    """Keep only the coarsest pyramid levels (avoids zarr chunk errors after bbox crop)."""
    for container in (sdata.images, sdata.labels):
        for key in list(container.keys()):
            tree = container[key]
            scales = sorted(tree.children.keys(), key=lambda s: int(s.replace("scale", "")))
            for scale in scales[:-n_keep]:
                if scale in tree.children:
                    del tree[scale]
    return sdata


def summarize_sdata(sdata: sd.SpatialData) -> None:
    print(f"  Images:  {list(sdata.images.keys())}")
    print(f"  Labels:  {list(sdata.labels.keys())}")
    print(f"  Shapes:  {', '.join(f'{k}={len(sdata.shapes[k]):,}' for k in sdata.shapes)}")
    print(f"  Points:  {', '.join(f'{k}={len(sdata.points[k]):,}' for k in sdata.points)}")
    if "table" in sdata.tables:
        print(f"  Table:   {sdata.tables['table'].n_obs:,} cells x {sdata.tables['table'].n_vars} genes")


def build_subset(
    sdata: sd.SpatialData,
    cx: float,
    cy: float,
    half: float,
) -> sd.SpatialData:
    sdata_sub = build_plot_crop(sdata, cx, cy, half, filter_table=True)
    # Bbox crop can leave multiscale images with inconsistent chunk shapes; keep coarse scales.
    return keep_coarsest_scales(sdata_sub, n_keep=2)


def write_and_measure(sdata_sub: sd.SpatialData, path: Path) -> float:
    if path.exists():
        shutil.rmtree(path)
    try:
        sdata_sub.write(path)
    except TypeError:
        sdata_sub = keep_coarsest_scales(sdata_sub, n_keep=1)
        sdata_sub.write(path)
    return dir_size_gb(path)


def tune_half(
    sdata: sd.SpatialData,
    coords: np.ndarray,
    target_gb: float,
    tmp_dir: Path,
    max_iters: int = 4,
) -> tuple[float, float, float, sd.SpatialData]:
    """Binary-search half-width; measure actual written zarr size each iteration."""
    cx, cy = cell_density_roi_center(coords, half=3000)
    lo, hi = 2000.0, 7000.0
    best: tuple[float, float, sd.SpatialData] | None = None

    for i in range(max_iters):
        half = (lo + hi) / 2
        print(f"\n[tune {i + 1}/{max_iters}] half={half:.0f}  center=({cx:.0f}, {cy:.0f})", flush=True)

        sdata_sub = build_subset(sdata, cx, cy, half)
        size_gb = write_and_measure(sdata_sub, tmp_dir)
        n_cells = sdata_sub.tables["table"].n_obs
        print(f"  cells={n_cells:,}  size={size_gb:.2f} GB", flush=True)

        if best is None or abs(size_gb - target_gb) < abs(best[0] - target_gb):
            best = (size_gb, half, sdata_sub)

        if size_gb < target_gb * 0.8:
            lo = half
        elif size_gb > target_gb * 1.2:
            hi = half
        else:
            return half, cx, cy, sdata_sub

    assert best is not None
    _, half, sdata_sub = best
    return half, cx, cy, sdata_sub


def main() -> None:
    parser = argparse.ArgumentParser(description="Subset TMA5.zarr to a dense ~1 GB SpatialData store.")
    parser.add_argument("--target-gb", type=float, default=1.0, help="Target output size in GB (default: 1.0)")
    parser.add_argument("--half", type=float, default=None, help="Fixed ROI half-width in px (skip auto-tuning)")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output zarr path")
    parser.add_argument("--input", type=Path, default=RAW_ZARR, help="Input TMA5.zarr path")
    args = parser.parse_args()

    t0 = time.time()
    print("=" * 60, flush=True)
    print("subset_tma5_zarr.py — started", flush=True)
    print("=" * 60, flush=True)
    print(f"Input:  {args.input}  ({dir_size_gb(args.input):.2f} GB)", flush=True)
    print(f"Output: {args.out}", flush=True)
    print(f"Target: {args.target_gb:.2f} GB", flush=True)

    print("\n[1/4] Loading SpatialData ...", flush=True)
    sdata = sd.read_zarr(args.input)
    coords = np.asarray(sdata.tables["table"].obsm["spatial"])
    print(f"  {coords.shape[0]:,} cells", flush=True)

    tmp_dir = args.out.with_suffix(".zarr.tmp")

    if args.half is not None:
        print(f"\n[2/4] Building subset with fixed half={args.half:.0f} ...", flush=True)
        cx, cy = cell_density_roi_center(coords, half=args.half)
        sdata_sub = build_subset(sdata, cx, cy, args.half)
        half = args.half
    else:
        print(f"\n[2/4] Tuning half to ~{args.target_gb:.2f} GB ...", flush=True)
        half, cx, cy, sdata_sub = tune_half(sdata, coords, args.target_gb, tmp_dir)

    print(f"\n  ROI: center=({cx:.0f}, {cy:.0f})  half={half:.0f}", flush=True)
    print(f"  bbox: x=[{cx - half:.0f}, {cx + half:.0f}]  y=[{cy - half:.0f}, {cy + half:.0f}]", flush=True)
    summarize_sdata(sdata_sub)

    print(f"\n[3/4] Writing {args.out} ...", flush=True)
    t_write = time.time()
    if not tmp_dir.exists() or args.half is not None:
        size_gb = write_and_measure(sdata_sub, tmp_dir)
    else:
        size_gb = dir_size_gb(tmp_dir)
    print(f"  write done in {time.time() - t_write:.1f}s  size={size_gb:.2f} GB", flush=True)

    if args.out.exists():
        shutil.rmtree(args.out)
    tmp_dir.rename(args.out)

    n_cells = sdata_sub.tables["table"].n_obs
    print(f"\n[4/4] Done in {(time.time() - t0) / 60:.1f} min", flush=True)
    print(f"Output: {args.out}  ({size_gb:.2f} GB)", flush=True)
    print(f"Cells:  {n_cells:,}", flush=True)


if __name__ == "__main__":
    main()
