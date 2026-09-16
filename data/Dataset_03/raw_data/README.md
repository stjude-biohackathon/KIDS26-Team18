Dump in BioHackathon_SJ_2026/Dataset_03/raw_data/ all the data from GEODataSets or whatever you have downloaded the data from


Human kidney Xenium data (~377 gene panel + 27 protein targets): [10x Xenium Protein FFPE Human Renal Carcinoma](https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma)

## Raw download

Original zip: `/mnt/scratch2/Chloe/Biohackathon/Dataset 3/Xenium_V1_Human_Kidney_FFPE_Protein_updated_outs.zip`

Native Xenium outs are staged from `Unzipped_Dataset_3/` into `xenium_bundle/` (hardlink when permitted, otherwise symlink — no extra disk). Xenium v4 outs also require `cells.zarr.zip`.

## QC notebook

[`00_load_and_explore.ipynb`](00_load_and_explore.ipynb) — load `sdata.zarr`, verify raw counts in `.X`, and plot centroids, cell boundaries, and transcripts.

## Build SpatialData + AnnData

Requires the `spatialdata` conda env (Python ≥3.11, spatialdata ≥0.7.2):

```bash
conda activate spatialdata
python build_sdata.py
```

Options: `--force` (rebuild), `--skip-stage` (reuse `xenium_bundle/`), `--source <path>` (override Xenium outs dir), `--validate-only`.

## Outputs (this folder)

| File / dir | Description |
|------------|-------------|
| `xenium_bundle/` | Staged Xenium outs (linked from source + `morphology_mip.ome.tif`) |
| `sdata.zarr` | SpatialData: cell table, boundaries, circles, transcript points |
| `adata.h5ad` | AnnData extracted from `sdata.tables["table"]` (465,534 cells × 405 genes) |
| `build_sdata.py` | Staging + `spatialdata_io.xenium()` build script |

Note: morphology image is omitted from `sdata.zarr` (full `morphology.ome.tif` stack incompatible with zarr write). Nucleus boundaries are omitted (incomplete vs cell table in export). Use `xenium_bundle/morphology_mip.ome.tif` for plotting.
