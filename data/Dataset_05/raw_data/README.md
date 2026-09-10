Dump in BioHackathon_SJ_2026/Dataset_05/raw_data/ all the data from GEODataSets or whatever you have downloaded the data from


Human skin Xenium data https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280

## Raw download

GEO files are in `GSE301280_RAW/` (GSM9079610, sample `11027_A2_Healthy`).

## QC notebook

[`00_load_and_explore.ipynb`](00_load_and_explore.ipynb) — load `sdata.zarr`, verify raw counts in `.X`, and plot centroids, cell boundaries, and transcripts.

## Build SpatialData + AnnData

Requires the `spatialdata` conda env (Python ≥3.11, spatialdata ≥0.7.2):

```bash
conda activate spatialdata
python build_sdata.py
```

Options: `--force` (rebuild), `--skip-stage` (reuse `xenium_bundle/`), `--validate-only`.

## Outputs (this folder)

| File / dir | Description |
|------------|-------------|
| `xenium_bundle/` | Staged Xenium outs (decompressed GEO files + `experiment.xenium`) |
| `sdata.zarr` | SpatialData: cell table, boundaries, circles, transcript points |
| `adata.h5ad` | AnnData extracted from `sdata.tables["table"]` (9,960 cells × 100 genes) |
| `build_sdata.py` | Staging + `spatialdata_io.xenium()` build script |

Note: morphology image is omitted from `sdata.zarr` (GEO `morphology.ome.tif` has channel layout incompatible with zarr write). Nucleus boundaries are omitted (incomplete vs cell table in GEO export).
