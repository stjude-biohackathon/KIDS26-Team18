# Plan: Subset TMA5.zarr → ~1 GB dense-region SpatialData

**Status:** in progress (write step failing)  
**Source:** `../raw_data/TMA5.zarr` (27 GB, 628,860 cells, 121M transcripts)  
**Target:** `TMA5_subset.zarr` (~1 GB, dense cell region, all elements)

---

## Goal

Read full zarr → subset to a cell-dense bounding box → save a shareable SpatialData zarr.

---

## What already exists

| Item | Path |
|------|------|
| Source zarr | `raw_data/TMA5.zarr` |
| AnnData-only crop (28k cells) | `one_sample.h5ad` via `00_load_and_explore.ipynb` |
| ROI helpers | `src/spatial_plot.py` — `cell_density_roi_center`, `build_plot_crop` |
| Subset script (broken) | `subset_tma5_zarr.py` |

---

## Simple pipeline (target)

```python
sdata = sd.read_zarr(TMA5.zarr)
cx, cy = cell_density_roi_center(coords, half=HALF)
sdata_sub = build_plot_crop(sdata, cx, cy, HALF, filter_table=True)
sdata_sub = collapse_rasters_to_single_scale(sdata_sub)  # fix write
sdata_sub.write(TMA5_subset.zarr)
```

No binary-search tuning on first pass — pick `half` manually (~4500 px) and adjust if size is off.

---

## Known issues & fixes

### 1. Transcript under-counting
`query.bounding_box` alone drops Xenium transcript points.  
**Fix:** use `build_plot_crop()` (re-filters points on x/y in raw coords).

### 2. Multiscale image write fails after bbox crop
Errors seen:
- `TypeError: Expected an iterable of integers` (chunk shape mismatch)
- `KeyError: 'Could not find node at scale0'` (after deleting pyramid levels)

**Fix:** after crop, collapse each image/label to **single-scale** from coarsest level (`scale4`):

```python
ds = tree["scale4"].dataset
sdata.images[key] = Image2DModel.parse(
    ds["image"], dims=tuple(ds["image"].dims),
    c_coords=list(ds.coords["c"].values),
)
sdata.labels[key] = Labels2DModel.parse(ds["image"])
```

Do **not** delete `scale0` from a DataTree — spatialdata writer requires it.

### 3. Coordinate-system mismatch (CONFIRMED)
`obsm["spatial"]` range: x 47–11483, y 128–19492.  
Notebook bbox (200–4000, 16300–20000): **27,791 cells** by centroid mask.  
Same bbox via `query.bounding_box(global)`: **2 cells only**.

`cell_density_roi_center` and `build_plot_crop` use incompatible coords — explains empty crops.

**Fix options:**
- A) Use known-good notebook bbox, or convert obsm → global before bbox query
- B) Subset table by cell mask on `obsm["spatial"]`, then `sdata.subset(cell_ids)` for shapes/points/images
- C) Inspect transforms on shapes/table to derive correct global bbox

### 4. Script over-engineering
Current `subset_tma5_zarr.py` binary-searches `half` and writes up to 4× (slow, 121M transcript scan each time).  
**Fix:** single subset pass, one write, report size; optional `--half` override.

---

## Size estimates (dense ROI, single-scale images)

| half | cells (approx) | notes |
|------|----------------|-------|
| 3000 | ~174k | likely >1 GB (points dominate) |
| 4500 | ~323k | too large; try with coarse images only |
| 2000–2500 | ~80–120k | sweet spot for ~1 GB |

Points (~3.5 GB full) scale with transcript count; images small after single-scale collapse.

---

## Implementation steps

1. [x] Add `cell_density_roi_center` to `src/spatial_plot.py`
2. [ ] Fix `subset_tma5_zarr.py`:
   - remove `tune_half` loop (or make opt-in)
   - add `collapse_rasters_to_single_scale()`
   - default `half=2500`, `--target-gb` optional
3. [ ] Run script, verify output loads and size ~0.8–1.2 GB
4. [x] Update `processed_data/README.md`

---

## Run command

```bash
/mnt/scratch1/miniconda3/envs/spatialdata/bin/python \
  data/Dataset_01/processed_data/subset_tma5_zarr.py --half 2500
```

---

## Delete this file when done

Temporary working plan — remove after `TMA5_subset.zarr` is built and verified.
