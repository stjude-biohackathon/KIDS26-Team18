import anndata as ad
import spatialdata as sd
from pathlib import Path
import numpy as np

def write_adata_to_sdata(zarr_path, h5ad_path, output_path):
    # load data
    sdata = sd.read_zarr(zarr_path)
    adata = ad.read_h5ad(h5ad_path)

    # check for columns present in  adata.obs that aren't in sdata
    table = sdata.tables["table"]
    missing_cols =  list(set(adata.obs.columns.to_list()) - set(table.obs.columns.to_list()))

    # assign values from anndata to sdata table
    for col in missing_cols:
        print(f"Adding column {col}")
        sdata.tables["table"].obs[col] = adata.obs[col].values

    # Check dimensions
    sdata_dim = sdata.tables["table"].shape[0]
    adata_dim = adata.obs.shape[0]
    if sdata_dim != adata_dim:
        print(f"WARNING: sdata dimensions ({sdata_dim}) do not match adata dimensions ({adata_dim})")

    # Check that cell_ids match between sdata and adata
    sdata_labels = np.array(sdata.tables["table"].obs["cell_id"])
    adata_labels = np.array(adata.obs["cell_id"])

    not_matching_count = np.sum(sdata_labels != adata_labels)
    if not_matching_count > 0:
        print(f"WARNING: {not_matching_count} labels do not match between the sdata and adata objects!")

    # write new zarr
    storage_path = Path(output_path)
    sdata.write(storage_path, overwrite=True)
    print(f"Wrote sdata to {output_path}")
