# Dataset_01 — processed_data

Source: [GSE250346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346) — Human lung Xenium TMA (~340 gene panel)

## Files

| File | Description |
| --- | --- |
| `one_sample.h5ad` | Cropped TMA5 subset (~27,791 cells × 343 genes): integer counts in `.X`, spatial coords in `obsm['spatial']`, full published `.obs` metadata merged from `GSE250346_slim.h5ad`. |
| `../raw_data/TMA5.zarr` | SpatialData store (Zarr v3, spatialdata 0.7.2): morphology images, cell/nucleus labels, shapes, transcript points, cell table (628,860 cells). |
| `00_load_and_explore.ipynb` | Load `TMA5.zarr`, crop a spatial region, join h5ad metadata, build `one_sample.h5ad`. **Kernel:** `spatialdata`. |
| `01_celltypeing_compare.ipynb` | Load `one_sample.h5ad`; compare scanpy cluster-based manual annotation vs rule-based epithelial/B-cell typing. **Kernel:** `spatialdata`. |
| `outdated/` | Archived notebooks, scripts, and earlier h5ad versions — not part of the current workflow. |

## Workflow

1. **`00_load_and_explore.ipynb`** — read SpatialData from `TMA5.zarr`, crop cells by bounding box, merge all `.obs` columns from `outdated/GSE250346_slim.h5ad` on `cell_id`, optionally write `one_sample.h5ad`.
2. **`01_celltypeing_compare.ipynb`** — load `one_sample.h5ad` and run two cell-typing approaches (scanpy Leiden + manual markers; rule-based epithelial/B-cell panels), then compare against published `final_CT`.

## Key obs columns (in `one_sample.h5ad`)

`sample_id`, `donor_id`, `cell_id`, `sample_type`, `disease_status`, `TMA_core`, `final_CT`, `final_lineage`, `CNiche`, `TNiche`, `coord_x`, `coord_y`, `nCount`, `nFeature`, `cell_type_scanpy`, `cell_type_rule`, ...

## Kernel

Use **`spatialdata`** (Python ≥3.11, spatialdata ≥0.7.2) for both notebooks.
