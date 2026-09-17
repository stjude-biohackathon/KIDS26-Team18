# KIDS26-Team18 — shared analysis utilities

Python helpers for notebooks and scripts in this repo. Add `src/` to `sys.path`
(or call `setup_notebook_paths`) before importing.

Dataset catalog, BOX downloads, and collaboration guide → [`../data/README.md`](../data/README.md)

**Kernel:** `spatialdata` conda env (`/mnt/scratch1/miniconda3/envs/spatialdata`)

---

## Quick start (notebooks)

```python
from pathlib import Path
import sys

RAW_DIR = Path(".../data/Dataset_XX/raw_data")
sys.path.insert(0, str(RAW_DIR.parent.parent.parent / "src"))

from load_explore import setup_notebook_paths
from spatial_plot import plot_gene_seg_transcripts, build_plot_crop

setup_notebook_paths(raw_dir=RAW_DIR)
```

For AnnData / scanpy notebooks under `processed_data/`:

```python
from pathlib import Path
import sys

PROC_DIR = Path(".../data/Dataset_01/processed_data")
sys.path.insert(0, str(PROC_DIR.parent.parent.parent / "src"))

from plot import SpatialCoord_plot, custom_barplot
```

---

## Module map

### SpatialData QC (`spatial_plot.py`, `load_explore.py`)

For `00_load_and_explore.ipynb` on Xenium / CosMx zarr stores.

| Function | Purpose |
|----------|---------|
| `setup_notebook_paths` | Insert `src/` on `sys.path`, resolve repo root |
| `find_repo_root` | Walk parents to locate repo root (`src/` + `data/`) |
| `attach_morphology_mip` | Attach 2D DAPI MIP when missing from zarr |
| `transcript_roi_center` | Find ROI center from transcript density |
| `build_plot_crop` | Bbox crop with correct transcript filtering |
| `expression_roi_center` | ROI from AnnData expression peaks |
| `plot_gene_seg_transcripts` | Segmentation + transcript overlay (matplotlib) |

Works across modalities: infers `feature_name` vs `target`, `cell_id` vs `cell_uid`, optional `z`.

### AnnData / Scanpy (`load_data.py`, `plot.py`, `io.py`)

Extracted from the legacy `Utils.py` codebase. For h5ad / metadata workflows.

| Module | Key exports |
|--------|-------------|
| `load_data` | `read_adata`, `attach_metadata`, `load_expression_matrix`, `get_noncoding_genes` |
| `plot` | `custom_barplot`, `dotplot`, `plot_umap`, `SpatialCoord_plot`, `set_scanpy_colors` |
| `io` | `save_obs_to_parquet` |

### BioHack template (`biohack_template/`)

Self-contained cell-typing / QC pipeline for TMA-style workflows. Uses `lib/` internally.

| Module | Purpose |
|--------|---------|
| `celltyping_scanpy` | Scanpy clustering + marker assignment |
| `celltyping_rule_based` | Rule-based typing from marker panels |
| `celltyping_compare` | Compare typing methods |
| `qc` | Basic QC filtering, cells-per-core summary |
| `plotting` | QC bars, composition, confusion matrix, spatial core plots |
| `synthetic_tma` | Synthetic demo AnnData |
| `spatial_hypothesis` | Tumor–macrophage distance tests |

Config: `biohack_template/demo_marker_panels.yaml`

### Low-level (`lib/`)

Internal helpers for `biohack_template` — `qc_common`, `config`, `paths`.

---

## Which notebook uses what?

| Notebook | Modules |
|----------|---------|
| `data/Dataset_02–06/raw_data/00_load_and_explore.ipynb` | `load_explore`, `spatial_plot` |
| `data/Dataset_01/processed_data/01_celltypeing_compare.ipynb` | `plot` |
| `data/Dataset_01/processed_data/00_load_and_explore.ipynb` | raw spatialdata/scanpy (not wired to `src` yet) |

---

## Design notes

- **Single source of truth:** all shared code lives under `src/`. Do not add copies under `data/`.
- **Transcript overlays:** `plot_gene_seg_transcripts` forces `method="matplotlib"` because datashader misaligns large Xenium crops (see `data/Dataset_03/raw_data/DEBUGGING.md`).
- **`build_plot_crop`:** materializes filtered transcripts via `PointsModel.parse()` to avoid dask partition errors on large datasets.
