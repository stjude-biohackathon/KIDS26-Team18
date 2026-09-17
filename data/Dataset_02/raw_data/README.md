Dump in BioHackathon_SJ_2026/Dataset_02/raw_data/ all the data from GEODataSets or whatever you have downloaded the data from


Human lung CosMx data https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299786

## Raw download

Nanostring flat files are in `GSE299786_RAW/` (1k panel, GSM9046088 `Lung_Adenocarcinoma_TMA1`).

## QC notebook

[`00_load_and_explore.ipynb`](00_load_and_explore.ipynb) — load `sdata.zarr`, verify raw counts in `.X`, and plot centroids, cell boundaries, and transcripts.

## Build SpatialData + AnnData

Requires the `spatialdata` conda env:

```bash
conda activate spatialdata
python build_sdata.py
```

Options: `--force` (rebuild), `--skip-transcripts` (table + polygons only), `--validate-only`.

First run caches transcripts to `GSE299786_RAW/transcripts.global_coords.parquet` (~18M detections).

## Outputs (this folder)

| File / dir | Description |
|------------|-------------|
| `sdata.zarr` | SpatialData: cell table, boundary polygons, transcript points |
| `adata.h5ad` | AnnData extracted from `sdata.tables["table"]`; `obsm["global"]` and `obsm["spatial"]` |
| `build_sdata.py` | Flat-file → SpatialData builder (no morphology images in GEO export) |
