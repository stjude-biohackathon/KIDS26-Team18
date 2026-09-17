# spatialdata environment

Conda environment for **Dataset_02–06** QC notebooks (`00_load_and_explore.ipynb`), `build_sdata.py` scripts, and shared [`src/spatial_plot.py`](../src/spatial_plot.py) helpers.

**Stack:** Python 3.12, `spatialdata` 0.7.2, Zarr v3.

Dataset catalog, BOX downloads, and collaboration guide → [`../data/README.md`](../data/README.md)

---

## Use the existing cluster env (recommended)

On this machine the working env is already installed:

```bash
conda activate spatialdata
# path: /mnt/scratch1/miniconda3/envs/spatialdata
```

Jupyter notebooks for Datasets 02–06 use the **`spatialdata`** kernel (`display_name: spatialdata`).

---

## Create from this repo (new machine)

From the repository root:

```bash
conda env create -f env/environment.yml
conda activate spatialdata
```

Or install pip packages into an existing Python 3.12 env:

```bash
pip install -r env/requirements.txt
```

## Register Jupyter kernel (optional)

```bash
conda activate spatialdata
python -m ipykernel install --user --name spatialdata --display-name "spatialdata"
```

---

## Verify imports

```bash
conda activate spatialdata
python -c "
import spatialdata as sd
import spatialdata_plot
import spatialdata_io
import scanpy as sc
import zarr
import dask
import tifffile
print('spatialdata', sd.__version__)
print('zarr', zarr.__version__)
print('scanpy', sc.__version__)
print('dask', dask.__version__)
print('OK')
"
```

---

## Notebooks using this env

| Dataset | Notebook | Build script |
|---------|----------|--------------|
| 02 — Lung CosMx | `data/Dataset_02/raw_data/00_load_and_explore.ipynb` | `build_sdata.py` |
| 03 — Kidney Xenium | `data/Dataset_03/raw_data/00_load_and_explore.ipynb` | `build_sdata.py` |
| 05 — Skin Xenium | `data/Dataset_05/raw_data/00_load_and_explore.ipynb` | `build_sdata.py` |
| 06 — Skin CosMx | `data/Dataset_06/raw_data/00_load_and_explore.ipynb` | `build_sdata.py` |

Dataset_01 notebooks also run on this env (or an equivalent install with spatialdata ≥0.7.2).

Requires **Python ≥3.11** (this spec pins 3.12) and **spatialdata ≥0.7.2** for Zarr v3 stores.
