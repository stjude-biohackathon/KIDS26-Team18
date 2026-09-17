# Jupyter workflow

These notebooks are an interactive companion to `scripts/` + `hpc/`. They do not replace the production scripts; they expose the same stages in notebook-sized checkpoints so parameters and figures can be inspected before re-running HPC jobs.

Recommended order:

1. `00_setup_and_objects.ipynb`
2. `01_whole_slide_exploration.ipynb`
3. `02_fov46_qc.ipynb`
4. `03_fov46_leiden_sweep.ipynb` — **mandatory manual checkpoint**
5. `04_approach_i.ipynb`
6. `05_approach_ii.ipynb`
7. `06_comparison_and_spatial.ipynb`

Start Jupyter/OnDemand with the `spatialdata` environment and open the notebooks from `<PROJECT_ROOT>/jupyter/`. The notebooks assume the current repository layout (`scripts/`, `src/`, `config/`, `data/`, `results/`).

The heavy build steps are intentionally still delegated to the tested Python scripts via `%run`. QC and Leiden notebooks add interactive inspection cells so thresholds/resolutions can be explored without changing the production code immediately. Once a parameter choice is accepted, update the corresponding config/script and re-run the production/HPC stage for reproducibility.
