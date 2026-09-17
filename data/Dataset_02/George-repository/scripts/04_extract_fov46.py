#!/usr/bin/env python
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import numpy as np
import scanpy as sc
import spatialdata as sd
FOV="46"
# Always subset AnnData from the untouched raw object, never the processed whole-slide object.
raw=sc.read_h5ad(ROOT/"data/processed/GSM9046088_CosMx_raw.h5ad")
mask=raw.obs["fov"].astype(str).eq(FOV)
raw_fov=raw[mask].copy(); out_h5=ROOT/f"data/processed/GSM9046088_CosMx_raw_FOV{FOV}.h5ad"; raw_fov.write_h5ad(out_h5,compression="gzip")
print(f"raw FOV{FOV}: {raw_fov.n_obs:,} cells -> {out_h5}")
# SpatialData subset: first use a bounding box around FOV46, then explicitly check whether
# padding introduced neighboring FOV cells. If so, retry without padding. If that still
# contains neighbors, use the original sdata and filter table/shapes/points by FOV IDs.
sdata=sd.read_zarr(ROOT/"data/spatial/GSM9046088_CosMx.zarr")
table=sdata["table"]; tm=table.obs["fov"].astype(str).eq(FOV); coords=table.obsm["spatial"][tm]
if len(coords)==0: raise ValueError("FOV46 absent from SpatialData table")

def bbox_subset(padding):
    lo=coords.min(0)-padding; hi=coords.max(0)+padding
    return sdata.query.bounding_box(axes=("x","y"),min_coordinate=lo,max_coordinate=hi,target_coordinate_system="global",filter_table=True)

def neighbor_fovs(obj):
    vals=set(obj["table"].obs["fov"].astype(str).unique())
    return sorted(vals-{FOV})
sub=bbox_subset(50)
neighbors=neighbor_fovs(sub)
if neighbors:
    print("Padding introduced neighboring FOVs:",neighbors,"; retrying padding=0")
    sub=bbox_subset(0); neighbors=neighbor_fovs(sub)
if neighbors:
    print("BBox still overlaps neighboring FOVs:",neighbors,"; applying explicit instance filtering")
    # Strict fallback: crop to FOV46 cell extent, then retain only FOV46 table rows and linked shapes/transcripts.
    wanted=set(table.obs_names[tm].astype(str)); sub=bbox_subset(0)
    st=sub["table"]; keep=st.obs_names.astype(str).isin(wanted); sub.tables["table"]=st[keep].copy()
    if "cell_boundaries" in sub.shapes:
        sub.shapes["cell_boundaries"]=sub.shapes["cell_boundaries"].loc[sub.shapes["cell_boundaries"].index.astype(str).isin(wanted)].copy()
    if "transcripts" in sub.points and "unique_cell_id" in sub.points["transcripts"].columns:
        tx=sub.points["transcripts"]; sub.points["transcripts"]=tx[tx["unique_cell_id"].isin(list(wanted))]
    neighbors=neighbor_fovs(sub)
if neighbors: raise RuntimeError(f"Strict FOV46 filtering failed; remaining neighbors: {neighbors}")
assert sub["table"].obs["fov"].astype(str).eq(FOV).all()
out_z=ROOT/f"data/spatial/GSM9046088_CosMx_FOV{FOV}.zarr"; sub.write(out_z,overwrite=True)
print(sub); print("saved",out_z)
