#!/usr/bin/env python3
"""Stage GEO Xenium files and build SpatialData + AnnData under raw_data/."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import tifffile

RAW_DIR = Path(__file__).resolve().parent
GEO_RAW = RAW_DIR / "GSE301280_RAW"
BUNDLE_DIR = RAW_DIR / "xenium_bundle"
SDATA_PATH = RAW_DIR / "sdata.zarr"
H5AD_PATH = RAW_DIR / "adata.h5ad"

STAGE_MAP = {
    "cell_boundaries.parquet": "*_cell_boundaries.parquet.gz",
    "cells.parquet": "*_cells.parquet.gz",
    "nucleus_boundaries.parquet": "*_nucleus_boundaries.parquet.gz",
    "transcripts.parquet": "*_transcripts.parquet.gz",
    "morphology.ome.tif": "*_morphology.ome.tif.gz",
    "cell_feature_matrix.h5": "*_cell_feature_matrix.h5",
}

EXPERIMENT_SPECS = {
    "analysis_sw_version": "xenium-1.9.0.0",
    "pixel_size": 0.2125,
    "region": "cell_circles",
}


def _find_one(pattern: str) -> Path:
    matches = sorted(GEO_RAW.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching {pattern!r} in {GEO_RAW}")
    return matches[0]


def _write_morphology_mip(morphology_path: Path, mip_path: Path) -> None:
    """Write a single-channel 2D MIP from the Xenium morphology Z-stack."""
    print(f"Computing morphology MIP from {morphology_path.name} ...")
    with tifffile.TiffFile(morphology_path) as tif:
        stack = tif.series[0].asarray()
    if stack.ndim == 2:
        mip = stack
    elif stack.ndim == 3:
        mip = stack.max(axis=0)
    else:
        raise ValueError(f"Unexpected morphology shape: {stack.shape}")
    tifffile.imwrite(mip_path, mip[None, ...].astype(np.uint16))
    print(f"Wrote {mip_path.name} with shape {mip.shape}")


def stage_bundle(force: bool = False) -> None:
    """Decompress GEO files into a spatialdata-io-compatible Xenium bundle."""
    if BUNDLE_DIR.exists() and not force:
        required = list(STAGE_MAP) + ["experiment.xenium"]
        mip_path = BUNDLE_DIR / "morphology_mip.ome.tif"
        if all((BUNDLE_DIR / name).exists() for name in required):
            if mip_path.is_symlink() or not mip_path.exists():
                _write_morphology_mip(BUNDLE_DIR / "morphology.ome.tif", mip_path)
            print(f"Staging skipped — bundle already present: {BUNDLE_DIR}")
            return

    if BUNDLE_DIR.exists() and force:
        shutil.rmtree(BUNDLE_DIR)
    BUNDLE_DIR.mkdir(parents=True, exist_ok=True)

    for target_name, pattern in STAGE_MAP.items():
        src = _find_one(pattern)
        dst = BUNDLE_DIR / target_name
        print(f"Staging {src.name} -> {target_name}")
        if pattern.endswith(".gz"):
            with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            shutil.copy2(src, dst)

    mip_path = BUNDLE_DIR / "morphology_mip.ome.tif"
    if mip_path.exists() or mip_path.is_symlink():
        mip_path.unlink()
    _write_morphology_mip(BUNDLE_DIR / "morphology.ome.tif", mip_path)

    specs_path = BUNDLE_DIR / "experiment.xenium"
    specs_path.write_text(json.dumps(EXPERIMENT_SPECS, indent=2) + "\n")
    print(f"Wrote {specs_path}")


def build_sdata(n_jobs: int = 4, force: bool = False) -> None:
    if SDATA_PATH.exists() and not force:
        print(f"Skipping build — {SDATA_PATH} already exists (use --force to rebuild)")
        return

    from spatialdata_io import xenium

    print("Reading Xenium bundle with spatialdata_io.xenium() ...")
    sdata = xenium(
        BUNDLE_DIR,
        cells_as_circles=True,
        cells_labels=False,
        nucleus_labels=False,
        nucleus_boundaries=False,
        morphology_focus=False,
        morphology_mip=False,  # MIP staged separately; zarr image write fails on this stack size
        aligned_images=False,
        n_jobs=n_jobs,
    )

    if SDATA_PATH.exists():
        shutil.rmtree(SDATA_PATH)
    print(f"Writing {SDATA_PATH} ...")
    sdata.write(SDATA_PATH)

    adata = sdata["table"]
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
    mip_path = BUNDLE_DIR / "morphology_mip.ome.tif"
    if mip_path.exists() and not mip_path.is_symlink():
        print(f"Morphology MIP available at: {mip_path}")
    else:
        print("Warning: morphology_mip.ome.tif missing — run staging to generate it.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Rebuild staging and outputs")
    parser.add_argument("--skip-stage", action="store_true", help="Skip GEO decompression")
    parser.add_argument("--n-jobs", type=int, default=4, help="Parallel jobs for xenium reader")
    parser.add_argument("--validate-only", action="store_true", help="Only run validation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.validate_only:
        validate()
        return

    if not args.skip_stage:
        stage_bundle(force=args.force)
    build_sdata(n_jobs=args.n_jobs, force=args.force)
    validate()
    print("Done.")


if __name__ == "__main__":
    main()
