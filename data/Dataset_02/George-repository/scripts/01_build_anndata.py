#!/usr/bin/env python
"""Build and save the raw CosMx AnnData object from downloaded GEO files.

This step is intentionally independent of SpatialData construction.
The output is the single source raw AnnData used by downstream analyses and
by ``02_build_sdata.py``.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.cosmx_io import build_raw_anndata

RAW_DIR = ROOT / "data" / "raw"
OUT_DIR = ROOT / "data" / "processed"
OUT_H5AD = OUT_DIR / "GSM9046088_CosMx_raw.h5ad"


def main() -> None:
    print("=" * 70)
    print("Building raw CosMx AnnData")
    print(f"Project root : {ROOT}")
    print(f"Raw data dir : {RAW_DIR}")
    print(f"Output       : {OUT_H5AD}")
    print("=" * 70)

    if not RAW_DIR.is_dir():
        raise FileNotFoundError(
            f"Raw data directory not found: {RAW_DIR}\n"
            "Run the GEO download step first."
        )

    adata = build_raw_anndata(RAW_DIR)

    print("\nRaw AnnData constructed:")
    print(adata)
    print(f"Cells    : {adata.n_obs:,}")
    print(f"Features : {adata.n_vars:,}")
    print(f"Has counts layer: {'counts' in adata.layers}")
    print(f"Has spatial coordinates: {'spatial' in adata.obsm}")

    if not adata.obs_names.is_unique:
        raise ValueError("AnnData obs_names are not unique.")

    if "counts" not in adata.layers:
        raise ValueError("Expected raw counts in adata.layers['counts'].")

    if "spatial" not in adata.obsm:
        raise ValueError("Expected global coordinates in adata.obsm['spatial'].")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT_H5AD)

    print(f"\nSaved raw AnnData:\n{OUT_H5AD}")
    print("AnnData build completed successfully.")


if __name__ == "__main__":
    main()
