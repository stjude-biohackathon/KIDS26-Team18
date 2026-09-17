"""Parquet I/O helpers extracted from Utils.py."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

def _backup_parquet_before_overwrite(parquet_path, verbose=True):
    parquet_path = Path(parquet_path)
    backup_path = parquet_path.parent / "bk.parquet"

    if parquet_path.is_file():
        shutil.copy2(parquet_path, backup_path)
        if verbose:
            print(f"Backup saved: {backup_path}")
    elif verbose:
        print(f"No existing parquet to back up: {parquet_path}")


def save_obs_to_parquet(adata, out_parquet, verbose=True):
    """
    Write `adata.obs` to parquet, backing up the existing file to `bk.parquet`
    before overwrite.

    Safety check:
      - If the output parquet already exists and the number of rows differs
        from the new object, abort the save.
    """
    out_parquet = Path(out_parquet)
    out_parquet.parent.mkdir(parents=True, exist_ok=True)

    obs_df = adata.obs.copy()
    if "cell_id" in obs_df.columns:
        obs_df = obs_df.drop(columns=["cell_id"])
    obs_df.insert(0, "cell_id", obs_df.index.astype(str))

    # Safety check
    if out_parquet.exists():
        old_nrows = pd.read_parquet(out_parquet, columns=["cell_id"]).shape[0]
        new_nrows = obs_df.shape[0]

        if old_nrows != new_nrows:
            print(
                "ERROR: You are trying to overwrite an object with a different "
                "number of rows (cells).\n"
                "Create a new output file if that's what you wish.\n"
                f"Existing: {old_nrows:,} cells\n"
                f"New:      {new_nrows:,} cells"
            )
            return None

    _backup_parquet_before_overwrite(out_parquet, verbose=verbose)

    obs_df.to_parquet(out_parquet, index=False)

    if verbose:
        print(f"Saved: {out_parquet}")
        print(obs_df.shape)

    return obs_df

# # Usage
# out_parquet = "results/sample1/obs.parquet"

# save_obs_to_parquet(
#     adata,
#     out_parquet,
#     verbose=True
# )
