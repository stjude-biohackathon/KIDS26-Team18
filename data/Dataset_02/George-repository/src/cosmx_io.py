from __future__ import annotations
import gzip
import re
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

DATASET_PREFIX = "GSM9046088_Lung_Adenocarcinoma_TMA1_CosMx"
FILES = {
    "expression": f"{DATASET_PREFIX}_exprMat_file.csv.gz",
    "metadata": f"{DATASET_PREFIX}_metadata_file.csv.gz",
    "fov": f"{DATASET_PREFIX}_fov_positions_file.csv.gz",
    "polygons": f"{DATASET_PREFIX}_polygons.csv.gz",
    "transcripts": f"{DATASET_PREFIX}_tx_file.csv.gz",
}
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9046nnn/GSM9046088/suppl/"

def unique_cell_id(fov, cell_id):
    return "fov_" + fov.astype(str) + "_cell_" + cell_id.astype(str)

def validate_gzip_csv(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().strip()
    if not header or "," not in header:
        raise ValueError(f"Unexpected gzip/CSV header in {path}: {header[:100]!r}")

def build_raw_anndata(raw_dir: Path) -> ad.AnnData:
    expr = pd.read_csv(raw_dir / FILES["expression"])
    meta = pd.read_csv(raw_dir / FILES["metadata"])
    fov = pd.read_csv(raw_dir / FILES["fov"])
    required = {"fov", "cell_ID"}
    for name, df in [("expression", expr), ("metadata", meta)]:
        if not required.issubset(df.columns):
            raise ValueError(f"{name} missing {required - set(df.columns)}")
        df["fov"] = pd.to_numeric(df["fov"], errors="raise").astype(int)
        df["cell_ID"] = pd.to_numeric(df["cell_ID"], errors="raise").astype(int)
        df["unique_cell_id"] = unique_cell_id(df["fov"], df["cell_ID"])
        if df["unique_cell_id"].duplicated().any():
            raise ValueError(f"Duplicated unique_cell_id in {name}")
        df.set_index("unique_cell_id", inplace=True)
    common = expr.index.intersection(meta.index, sort=False)
    if len(common) == 0:
        raise ValueError("No overlapping cells between expression and metadata")
    expr, meta = expr.loc[common].copy(), meta.loc[common].copy()
    non_gene = {"fov","FOV","cell_ID","cellID","cell","slide_ID","Slide","Run_Tissue_name"}
    candidate = [c for c in expr.columns if c not in non_gene]
    counts = expr[candidate].apply(pd.to_numeric, errors="coerce")
    counts = counts.loc[:, ~counts.isna().all(axis=0)].fillna(0)
    if (counts.to_numpy() < 0).any():
        raise ValueError("Negative expression counts detected")
    adata = ad.AnnData(
        X=sparse.csr_matrix(counts.to_numpy(dtype=np.float32)),
        obs=meta.copy(),
        var=pd.DataFrame(index=pd.Index(counts.columns.astype(str), name="gene")),
    )
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    coordinate_options = [
        ("CenterX_global_px", "CenterY_global_px"),
        ("CenterX_global_um", "CenterY_global_um"),
        ("CenterX_local_px", "CenterY_local_px"),
        ("x_global_px", "y_global_px"),
    ]
    for xcol, ycol in coordinate_options:
        if xcol in adata.obs and ycol in adata.obs:
            adata.obsm["spatial"] = adata.obs[[xcol,ycol]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
            adata.uns["spatial_coordinate_columns"] = {"x": xcol, "y": ycol}
            break
    adata.uns["fov_positions"] = fov.to_dict(orient="list")
    adata.uns["GEO_accession"] = "GSM9046088"
    adata.uns["platform"] = "NanoString CosMx SMI"
    adata.uns["sample"] = "Lung Adenocarcinoma TMA1"
    return adata

def control_mask(var_names):
    pat = re.compile(r"^(Negative|NegPrb|SystemControl|FalseCode|Control)", re.I)
    return np.array([bool(pat.search(str(g))) for g in var_names])
