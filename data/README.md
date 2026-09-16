# KIDS26-Team18 — Data & Collaboration Guide

**Cell Identity in Spatial OMICs Imaging Platforms: Cluster-Based vs Cell-Based Annotation**

This document is the main entry point for the repository: how datasets are organized, where to download shared outputs, and how the team collaborates. Shared Python utilities live in [`../src/README.md`](../src/README.md); analysis outputs in [`../result/`](../result/).

**Team lead:** Maycon Marção ([@Mmaycon](https://github.com/Mmaycon)) — see [`../project-management/team.md`](../project-management/team.md).

---

## About

Cell type annotation in imaging-based transcriptomics sits at the core of how we interpret biology, yet we still rely on two fundamentally different strategies without fully understanding how they relate to each other. **Cluster-based** annotation assigns identity at the population level; **cell-based** (rule-based) annotation applies explicit gene rules to classify individual cells. Imaging platforms detect transcripts as spatial spots, making rule-based annotation intuitive — but thresholds, marker specificity, and agreement with clustering remain open questions.

This repo tests both approaches across public Xenium and CosMx datasets, with shared QC notebooks and utilities so the team can compare methods on equal footing.

---

## Repository layout

```text
KIDS26-Team18/
├── data/          # datasets, raw builds, QC notebooks (you are here)
├── src/           # shared Python utilities
├── result/        # analysis outputs and templates
├── env/           # conda environment
└── project-management/
```

Inside every `data/Dataset_XX/` folder:

```text
Dataset_XX/
├── raw_data/                 # Original downloads + build scripts + sdata.zarr
├── metadata/                 # Sample tables and variable definitions
├── processed_data/           # Shared, standardized outputs (when applicable)
└── analysis/                 # Per-analyst workspaces
    └── <your_name>/
        ├── notebooks/
        ├── scripts/
        └── results/
```

| Folder | Purpose |
|--------|---------|
| `raw_data/` | Immutable source files. Document the source in `raw_data/README.md`. |
| `metadata/` | `sample_metadata.csv`, `variable_dictionary.xlsx`, cohort notes. |
| `processed_data/` | Team-facing `adata.h5ad`, `sdata.zarr` for downstream analyses. |
| `analysis/<your_name>/` | Personal notebooks, scripts, and intermediate results. |

---

## Links

### Download processed SpatialData (BOX)

> **Status:** https://stjude.app.box.com/folder/416940767803

| Archive | BOX link | Contents | Extract to |
|---------|----------|----------|------------|
| `zarr.tar.gz` | *TBD* | Pre-built `sdata.zarr` stores | `data/Dataset_XX/raw_data/` |

### External data sources

| Dataset | Accession / source |
|---------|-------------------|
| 01 — Lung Xenium TMA | [GSE250346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346) |
| 02 — Lung CosMx | [GSE299786](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299786) |
| 03 — Kidney Xenium | [10x Xenium Protein FFPE Human Renal Carcinoma](https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma) |
| 04 — TBD | — |
| 05 — Skin Xenium | [GSE301280](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301280) |
| 06 — Skin CosMx | [GSE314158](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314158) |

### In-repo documentation

| Doc | Description |
|-----|-------------|
| [`../src/README.md`](../src/README.md) | Shared helpers: `spatial_plot`, `load_explore`, scanpy plotting |
| [`../env/README.md`](../env/README.md) | Conda environment (`spatialdata`, Python 3.12) |
| [`../result/`](../result/) | Analysis outputs and notebook templates |
| [`Dataset_03/raw_data/DEBUGGING.md`](Dataset_03/raw_data/DEBUGGING.md) | Known Xenium transcript overlay issues and fixes |
| [`Dataset_01/README.md`](Dataset_01/README.md) | Dataset 01 (lung TMA) specifics |

---

## Dataset catalog

| Dataset | Platform | Tissue | Folder | QC notebook | Build |
|---------|----------|--------|--------|-------------|-------|
| 01 | Xenium TMA | Lung | [`Dataset_01/`](Dataset_01/) | [`processed_data/00_load_and_explore.ipynb`](Dataset_01/processed_data/00_load_and_explore.ipynb) | — |
| 02 | CosMx | Lung | [`Dataset_02/`](Dataset_02/) | [`raw_data/00_load_and_explore.ipynb`](Dataset_02/raw_data/00_load_and_explore.ipynb) | [`build_sdata.py`](Dataset_02/raw_data/build_sdata.py) |
| 03 | Xenium | Kidney | [`Dataset_03/`](Dataset_03/) | [`raw_data/00_load_and_explore.ipynb`](Dataset_03/raw_data/00_load_and_explore.ipynb) | [`build_sdata.py`](Dataset_03/raw_data/build_sdata.py) |
| 04 | — | — | [`Dataset_04/`](Dataset_04/) | — | — |
| 05 | Xenium | Skin | [`Dataset_05/`](Dataset_05/) | [`raw_data/00_load_and_explore.ipynb`](Dataset_05/raw_data/00_load_and_explore.ipynb) | [`build_sdata.py`](Dataset_05/raw_data/build_sdata.py) |
| 06 | CosMx | Skin | [`Dataset_06/`](Dataset_06/) | [`raw_data/00_load_and_explore.ipynb`](Dataset_06/raw_data/00_load_and_explore.ipynb) | [`build_sdata.py`](Dataset_06/raw_data/build_sdata.py) |

Source paths, cell counts, and build instructions: see each dataset's [`raw_data/README.md`](Dataset_02/raw_data/README.md).

**Kernel:** `spatialdata` conda env (Python ≥3.11, spatialdata ≥0.7.2). See [`../env/README.md`](../env/README.md).

---

## How to collaborate

### 1. Pick a dataset and claim your analyst folder

- Work inside one `Dataset_XX/` at a time unless coordinating across datasets.
- Add a folder under `analysis/`:

  ```bash
  mkdir -p data/Dataset_02/analysis/Your_Name/{env,notebooks,results,scripts}
  touch data/Dataset_02/analysis/Your_Name/README.md
  ```

- Keep personal exploration in **your** `analysis/<your_name>/` directory.
- Only write to `processed_data/` when outputs are reviewed and ready for the team.

### 2. Document raw data sources

In `raw_data/README.md`, record:

- Where the data came from (GEO accession, vendor portal, internal path)
- Download date and version
- Any preprocessing applied before placing files here

### 3. Share metadata early

- Fill in `metadata/sample_metadata.csv` with one row per sample (or per slide/FOV).
- Maintain `metadata/variable_dictionary.xlsx` so every column is defined.
- Update `metadata/README.md` with cohort-specific notes.

### 4. Coordinate on processed outputs

`processed_data/` is the **single source of truth**. Before overwriting `adata.h5ad` or `sdata.zarr`:

1. Announce the change to collaborators.
2. Confirm naming conventions below are met.
3. Leave a short note in `processed_data/README.md` (date, author, what changed).

### 5. Respect boundaries

| Do | Don't |
|----|-------|
| Work in `analysis/<your_name>/` | Edit another analyst's notebooks without permission |
| Read from `raw_data/` and `processed_data/` | Delete or rename files in `raw_data/` |
| Document assumptions in your README | Overwrite shared processed files silently |
| Use reproducible scripts in `scripts/` | Rely only on hard-coded notebook paths |

---

## Processed data standards

### AnnData (`.h5ad`)

| Component | Requirement |
|-----------|-------------|
| Raw counts | `adata.layers["counts"]` |
| Cell metadata | `adata.obs` — cell type column: `celltype_{your_name}` |
| Gene metadata | `adata.var` |
| Sample ID | `adata.obs`: `sample_id`, `TMA_core`, `donor_id`, `tissue_section_id` |
| Spatial coordinates | `adata.obs`: `coord_x`, `coord_y` |
| SpatialData link | `adata.obs`: `cell_id_to_sdata` (from `adata.obs_names`) |
| Embeddings | `adata.obs`: `UMAP1`, `UMAP2` (if computed) |

### SpatialData (`.zarr`)

| Element | Content |
|---------|---------|
| **Points** | Transcript-level coordinates (when available) |
| **Shapes** | Cell segmentation polygons or circles |
| **Images** | Histology or microscopy images |
| **Tables** | Standardized AnnData (`sdata.tables["table"]`) |

### Keeping AnnData and SpatialData in sync

1. Set `adata.obs["cell_id_to_sdata"]` from your cell identifiers.
2. Before updating `sdata.tables["table"]`:

   ```python
   adata.obs_names = adata.obs["cell_id_to_sdata"].values
   ```

3. `adata.obs_names` must match the AnnData stored in `sdata.tables["table"]`.

---

## Create a new dataset folder

```bash
DATASET_NUM=07
DATASET="data/Dataset_${DATASET_NUM}"

mkdir -p "$DATASET"/analysis/Analyst_Name01/{env,notebooks,results,scripts} \
         "$DATASET"/metadata \
         "$DATASET"/processed_data \
         "$DATASET"/raw_data

touch "$DATASET"/metadata/README.md \
      "$DATASET"/metadata/sample_metadata.csv \
      "$DATASET"/processed_data/README.md \
      "$DATASET"/raw_data/README.md
```

Fill in `raw_data/README.md` with the data source and add a row to the [dataset catalog](#dataset-catalog) above.

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

If conventions conflict with your platform (Xenium, CosMx, Visium, etc.), document the deviation in the dataset's `raw_data/README.md` and discuss with the team before merging into `processed_data/`.
