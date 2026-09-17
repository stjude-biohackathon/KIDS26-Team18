# Dataset 01 — Lung Xenium TMA (GSE250346)

GEO human lung Xenium TMA data ([GSE250346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346)), ~340 gene panel.

Project-wide links, dataset catalog, and collaboration standards → [**`../README.md`**](../README.md)

---

## Contents

| Path | Description |
|------|-------------|
| `raw_data/TMA5.zarr` | SpatialData store (Zarr v3): morphology, labels, shapes, transcripts, cell table (628,860 cells) |
| `processed_data/one_sample.h5ad` | Cropped TMA5 subset (~27,791 cells × 343 genes) with published metadata |
| `processed_data/00_load_and_explore.ipynb` | Load zarr, crop region, build `one_sample.h5ad` |
| `processed_data/01_celltypeing_compare.ipynb` | Compare scanpy cluster-based vs rule-based cell typing |
| `metadata/` | Sample tables and variable dictionary |
| `analysis/` | Per-analyst workspaces |

Details on processed files and workflow → [`processed_data/README.md`](processed_data/README.md)

---

## Notebooks

| Notebook | Kernel | Purpose |
|----------|--------|---------|
| `processed_data/00_load_and_explore.ipynb` | `spatialdata` | SpatialData exploration, crop, h5ad export |
| `processed_data/01_celltypeing_compare.ipynb` | `spatialdata` | Cluster vs rule-based typing comparison |

Requires Python ≥3.11 and spatialdata ≥0.7.2. See [`../../env/README.md`](../../env/README.md).
