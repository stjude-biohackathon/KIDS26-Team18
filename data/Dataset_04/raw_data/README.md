Dump in BioHackathon_SJ_2026/Dataset_04/raw_data/ all the data from GEODataSets or whatever you have downloaded the data from

Human kidney CosMx WTX data https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282026

## Raw download

Nanostring flat files are in `GSE282026_RAW/` (WTX panel, sample `KidneyRun2Slide1` / `KidneyRun2S1`).

Download and extract:

```bash
wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE282026&format=file&file=GSE282026%5FKidneyRun2Slide1%2Etar%2Egz" \
  -O GSE282026_KidneyRun2Slide1.tar.gz
mkdir -p GSE282026_RAW
tar -xzf GSE282026_KidneyRun2Slide1.tar.gz -C GSE282026_RAW --strip-components=1
```

## QC notebook

[`00_load_and_explore.ipynb`](00_load_and_explore.ipynb) — load `sdata.zarr`, verify raw counts in `.X`, and plot centroids, cell boundaries, and transcripts.

## Build SpatialData + AnnData

Requires the `spatialdata` conda env:

```bash
conda activate spatialdata
python build_sdata.py
```

Options: `--force` (rebuild), `--skip-transcripts` (table + polygons only), `--validate-only`.

First run caches transcripts to `GSE282026_RAW/transcripts.global_coords.parquet` (large tx file; expect longer runtime than Dataset 06).

## Outputs (this folder)

| File / dir | Description |
|------------|-------------|
| `sdata.zarr` | SpatialData: cell table, boundary polygons, transcript points |
| `adata.h5ad` | AnnData (~261k cells × ~6.5k features); `obsm["global"]` and `obsm["spatial"]` |
| `build_sdata.py` | Flat-file → SpatialData builder (no morphology images in GEO export) |
