# biohack-spatial environment

Conda environment for Dataset_01 spatial analysis notebooks (`spatialdata` 0.7.2, Zarr v3).

Derived from the working `spatialdata` conda env at `/mnt/scratch1/miniconda3/envs/spatialdata`.

## Create and activate

From the repository root:

```bash
conda env create -f env/environment.yml
conda activate biohack-spatial
```

## Register Jupyter kernel (optional)

```bash
python -m ipykernel install --user --name biohack-spatial --display-name "biohack-spatial"
```

## Verify imports

```bash
python -c "
import spatialdata as sd
import spatialdata_plot
import scanpy as sc
import zarr
print('spatialdata', sd.__version__)
print('zarr', zarr.__version__)
print('scanpy', sc.__version__)
print('OK')
"
```

## Run Dataset_01 notebooks

Open in Jupyter with the `biohack-spatial` (or `spatialdata`) kernel:

- `data/Dataset_01/processed_data/00_load_and_explore.ipynb`
- `data/Dataset_01/processed_data/01_celltypeing_compare.ipynb`

Requires Python ≥3.11 (this env uses 3.12) and `spatialdata` ≥0.7.2 to open `TMA5.zarr`.
