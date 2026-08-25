# biohack-spatial environment

Conda environment for the KIDS26-Team18 spatial analysis template notebook.

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
import sys
from pathlib import Path
sys.path.insert(0, str(Path('.').resolve() / 'src'))
from biohack_template import load_data, assign_celltypes_rule_based
from plot import SpatialCoord_plot
print('OK')
"
```

## Run the template notebook

```bash
cd result/template
jupyter nbconvert --to notebook --execute TEMPLATE_spatial_analysis.ipynb --output TEMPLATE_spatial_analysis.ipynb
```

Or open `result/template/TEMPLATE_spatial_analysis.ipynb` in Jupyter and run all cells.
