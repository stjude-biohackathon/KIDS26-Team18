# BioHackathon SJ 2026 — Collaboration Guide

This document explains how to work in this repository: how datasets are organized, what processed outputs we expect, and how analysts should collaborate without stepping on each other's work.

---

## Repository layout

Each dataset lives in its own top-level folder:

```text
BioHackathon_SJ_2026/
├── Dataset_01/
├── Dataset_02/
├── Dataset_03/
└── Dataset_04/
```

Inside every `Dataset_XX/` folder, use the same structure:

```text
Dataset_XX/
├── README.md                 # Dataset-specific notes (optional)
├── raw_data/                 # Original downloads (GEO, vendor files, etc.)
├── metadata/                 # Sample tables and variable definitions
├── processed_data/           # Shared, standardized outputs for the team
└── analysis/                 # Per-analyst workspaces
    ├── README.md
    ├── Analyst_Name01/
    │   ├── env/
    │   ├── notebooks/
    │   ├── results/
    │   ├── scripts/
    │   └── README.md
    └── Analyst_Name02/
        └── ...
```

| Folder | Purpose |
|--------|---------|
| `raw_data/` | Immutable source files. Do not edit in place; document the source in `README.md`. |
| `metadata/` | Curated sample-level metadata (`sample_metadata.csv`) and a data dictionary (`variable_dictionary.xlsx`). |
| `processed_data/` | **Team-facing** standardized objects (`adata.h5ad`, `sdata.zarr`) used for downstream analyses. |
| `analysis/<your_name>/` | **Your** notebooks, scripts, environment files, and intermediate results. |

---

## How to collaborate

### 1. Pick a dataset and claim your analyst folder

- Work inside one `Dataset_XX/` at a time unless you are explicitly coordinating across datasets.
- Rename `Analyst_Name01` / `Analyst_Name02` to your name, or add a new folder under `analysis/`:

  ```bash
  mkdir -p Dataset_02/analysis/Your_Name/{env,notebooks,results,scripts}
  touch Dataset_02/analysis/Your_Name/README.md
  ```

- Keep all personal exploration in **your** `analysis/<your_name>/` directory.
- Only write to `processed_data/` when outputs are reviewed and ready for the rest of the team.

### 2. Document raw data sources

In `raw_data/README.md`, record:

- Where the data came from (GEO accession, vendor portal, internal path, etc.)
- Download date and version
- Any preprocessing applied before placing files here

Example:

```markdown
Source: GEO GSE123456
Downloaded: 2026-07-23
Files: count matrix, spatial coordinates, H&E image
```

### 3. Share metadata early

- Fill in `metadata/sample_metadata.csv` with one row per sample (or per slide/FOV, as appropriate).
- Maintain `metadata/variable_dictionary.xlsx` so every column name and value is defined.
- Update `metadata/README.md` with cohort-specific notes (batch effects, exclusions, etc.).

### 4. Coordinate on processed outputs

`processed_data/` is the **single source of truth** for standardized objects. Before overwriting `adata.h5ad` or `sdata.zarr`:

1. Announce the change to collaborators (chat, issue, or PR description).
2. Confirm naming conventions below are met.
3. Leave a short note in `processed_data/README.md` (date, author, what changed).

Do **not** commit large intermediate `.h5ad` files from personal notebooks into `processed_data/` unless they follow the standards below.

### 5. Respect boundaries

| Do | Don't |
|----|-------|
| Work in `analysis/<your_name>/` | Edit another analyst's notebooks or results without permission |
| Read from `raw_data/` and `processed_data/` | Delete or rename files in `raw_data/` |
| Document assumptions in your `README.md` | Overwrite shared processed files silently |
| Use reproducible scripts in `scripts/` | Rely only on notebook state with hard-coded paths |

---

## Processed data standards

### AnnData (`.h5ad`)

Store the AnnData object as `processed_data/adata.h5ad`.

| Component | Requirement |
|-----------|-------------|
| Raw counts | `adata.layers["counts"]` |
| Cell metadata | `adata.obs` — final cell type column: `celltype_{your_name}` |
| Gene metadata | `adata.var` |
| Sample ID | `adata.obs`: `sample_id`, `TMA_core`, `donor_id`, `tissue_section_id` |
| Spatial coordinates | `adata.obs`: `coord_x`, `coord_y` |
| SpatialData link | `adata.obs`: `cell_id_to_sdata` (derived from `adata.obs_names`) |
| Embeddings | `adata.obs`: `UMAP1`, `UMAP2` |

### SpatialData (`.zarr`)

When cell segmentation polygons and/or transcript coordinates are available, also create `processed_data/sdata.zarr`.

| Element | Content |
|---------|---------|
| **Points** | Transcript-level coordinates (when available) |
| **Shapes** | Cell segmentation polygons |
| **Images** | Histology or microscopy images |
| **Tables** | The standardized AnnData object (`sdata.tables["table"]`) |

### Keeping AnnData and SpatialData in sync

1. Set `adata.obs["cell_id_to_sdata"]` from your cell identifiers.
2. Before updating `sdata.tables["table"]`, assign:

   ```python
   adata.obs_names = adata.obs["cell_id_to_sdata"].values
   ```

3. `adata.obs_names` must match the observation names of the AnnData stored in `sdata.tables["table"]`.

---

## Create a new dataset folder

Use this script when adding `Dataset_05`, `Dataset_06`, etc. Change `DATASET_NUM` accordingly.

```bash
DATASET_NUM=05
DATASET="Dataset_${DATASET_NUM}"

mkdir -p "$DATASET"/analysis/Analyst_Name01/{env,notebooks,results,scripts} \
         "$DATASET"/analysis/Analyst_Name02/{env,notebooks,results,scripts} \
         "$DATASET"/metadata \
         "$DATASET"/processed_data \
         "$DATASET"/raw_data

touch "$DATASET"/analysis/README.md \
      "$DATASET"/analysis/Analyst_Name01/README.md \
      "$DATASET"/analysis/Analyst_Name02/README.md \
      "$DATASET"/metadata/README.md \
      "$DATASET"/metadata/sample_metadata.csv \
      "$DATASET"/metadata/variable_dictionary.xlsx \
      "$DATASET"/processed_data/README.md \
      "$DATASET"/processed_data/adata.h5ad \
      "$DATASET"/processed_data/sdata.zarr \
      "$DATASET"/raw_data/README.md
```

Then add a short `Dataset_XX/README.md` and fill in `raw_data/README.md` with the data source.

---

## Quick checklist before sharing processed data

- [ ] Raw counts in `adata.layers["counts"]`
- [ ] Spatial coords in `adata.obs['coord_x']`, `adata.obs['coord_y']`
- [ ] Cell type column named `adata.obs['celltype_{your_name}']`
- [ ] `adata.obs['cell_id_to_sdata']` populated from `adata.obs_names`
- [ ] UMAP in `UMAP1`, `UMAP2` (if computed)
- [ ] SpatialData `.zarr` present when polygons/transcripts/images exist
- [ ] `metadata/sample_metadata.csv` and dictionary updated
- [ ] `processed_data/README.md` describes the current version

---

## Questions?

If conventions conflict with your platform (Xenium, Visium, MERFISH, etc.), document the deviation in your dataset's `README.md` and discuss with the team before merging into `processed_data/`.
