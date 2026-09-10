Dump in BioHackathon_SJ_2026/Dataset_06/raw_data/ all the data from GEODataSets or whatever you have downloaded the data from


Human skin CosMx data https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314158

## Raw download

Nanostring flat files are in `GSE314158_RAW/` (6k panel, `Cx009_S1_TMA2_CJ_6k_CSPI_Sk_250606`).

## QC notebook

[`00_load_and_explore.ipynb`](00_load_and_explore.ipynb) — load `sdata.zarr`, verify raw counts in `.X`, and plot centroids, cell boundaries, and transcripts.

## Build SpatialData + AnnData

Requires the `spatialdata` conda env:

```bash
conda activate spatialdata
python build_sdata.py
```

Options: `--force` (rebuild), `--skip-transcripts` (table + polygons only), `--validate-only`.

First run caches transcripts to `GSE314158_RAW/transcripts.global_coords.parquet` (~4 min for ~97M detections).

## Outputs (this folder)

| File / dir | Description |
|------------|-------------|
| `sdata.zarr` | SpatialData: cell table, boundary polygons, transcript points |
| `adata.h5ad` | AnnData (71,138 cells × 6,537 genes); `obsm["global"]` and `obsm["spatial"]` |
| `build_sdata.py` | Flat-file → SpatialData builder (no morphology images in GEO export) |
