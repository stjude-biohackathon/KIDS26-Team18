from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
from .cosmx_io import control_mask

def _row_sum(x): return np.asarray(x.sum(axis=1)).ravel()

def add_qc_metrics(adata):
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    adata.var["is_control"] = control_mask(adata.var_names)
    cm = adata.var["is_control"].to_numpy()
    adata.obs["control_counts"] = _row_sum(adata[:, cm].X) if cm.any() else 0
    adata.obs["control_fraction"] = (adata.obs["control_counts"] / adata.obs["total_counts"].replace(0,np.nan)).fillna(0)
    # Prefer CosMx metadata cell area if present.
    for c in ["Area.um2", "Area", "cell_area", "CellArea"]:
        if c in adata.obs:
            adata.obs["qc_cell_area"] = pd.to_numeric(adata.obs[c], errors="coerce")
            break
    return adata

def robust_bounds(x, lower_nmads=3.0, upper_nmads=4.0, log1p=False):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
    z = np.log1p(x) if log1p else x.copy()
    med = np.nanmedian(z); mad = np.nanmedian(np.abs(z-med))
    if not np.isfinite(mad) or mad == 0:
        lo, hi = np.nanquantile(z, [0.01, 0.99])
    else:
        lo, hi = med-lower_nmads*mad, med+upper_nmads*mad
    if log1p: lo, hi = np.expm1(lo), np.expm1(hi)
    return float(max(lo,0)), float(hi)

def recommend_thresholds(adata, fov_col="fov", nmads=3.0):
    # Hard floors prevent extremely low-information segmented objects from surviving
    # when an FOV itself is poor. Robust bounds then adapt to this dataset/FOV.
    global_count_lo, global_count_hi = robust_bounds(adata.obs["total_counts"], nmads, 4, True)
    global_gene_lo, global_gene_hi = robust_bounds(adata.obs["n_genes_by_counts"], nmads, 4, True)
    _, control_hi = robust_bounds(adata.obs["control_fraction"], nmads, nmads, False)
    rec = {
        "min_counts": max(20, int(np.floor(global_count_lo))),
        "min_genes": max(10, int(np.floor(global_gene_lo))),
        "max_counts": int(np.ceil(global_count_hi)),
        "max_control_fraction": min(1.0, control_hi),
    }
    if "qc_cell_area" in adata.obs:
        area_lo, area_hi = robust_bounds(adata.obs["qc_cell_area"], nmads, nmads, True)
        rec.update(min_cell_area=area_lo, max_cell_area=area_hi)
    # Per-FOV lower outlier flags: important for TMA/CosMx where FOV distributions differ.
    if fov_col in adata.obs:
        low_count = pd.Series(False,index=adata.obs_names); low_gene = low_count.copy()
        for _, idx in adata.obs.groupby(fov_col, observed=True).groups.items():
            sub=adata.obs.loc[idx]
            c_lo,_=robust_bounds(sub["total_counts"], nmads, 4, True)
            g_lo,_=robust_bounds(sub["n_genes_by_counts"], nmads, 4, True)
            low_count.loc[idx]=sub["total_counts"] < max(20,c_lo)
            low_gene.loc[idx]=sub["n_genes_by_counts"] < max(10,g_lo)
        adata.obs["qc_low_count_within_fov"] = low_count
        adata.obs["qc_low_gene_within_fov"] = low_gene
    return rec

def apply_qc(adata, thresholds, use_per_fov_flags=True):
    o=adata.obs
    keep=(o.total_counts>=thresholds["min_counts"]) & (o.n_genes_by_counts>=thresholds["min_genes"])
    if thresholds.get("max_counts") is not None: keep &= o.total_counts<=thresholds["max_counts"]
    if thresholds.get("max_control_fraction") is not None: keep &= o.control_fraction<=thresholds["max_control_fraction"]
    if "qc_cell_area" in o:
        if thresholds.get("min_cell_area") is not None: keep &= o.qc_cell_area>=thresholds["min_cell_area"]
        if thresholds.get("max_cell_area") is not None: keep &= o.qc_cell_area<=thresholds["max_cell_area"]
    if use_per_fov_flags:
        if "qc_low_count_within_fov" in o: keep &= ~o.qc_low_count_within_fov
        if "qc_low_gene_within_fov" in o: keep &= ~o.qc_low_gene_within_fov
    adata.obs["qc_pass"] = keep
    return keep

def threshold_sweep(adata, count_values=(10,20,30,50,75,100), gene_values=(5,10,15,20,30), control_values=(None,0.05,0.10,0.20)):
    rows=[]; o=adata.obs
    for mc in count_values:
      for mg in gene_values:
       for cf in control_values:
        keep=(o.total_counts>=mc)&(o.n_genes_by_counts>=mg)
        if cf is not None: keep &= o.control_fraction<=cf
        row={"min_counts":mc,"min_genes":mg,"max_control_fraction":cf,"n_cells":int(keep.sum()),"retained_fraction":float(keep.mean())}
        if "fov" in o:
            rates=keep.groupby(o["fov"]).mean()
            row.update(min_fov_retention=float(rates.min()), median_fov_retention=float(rates.median()), fov_retention_iqr=float(rates.quantile(.75)-rates.quantile(.25)))
        rows.append(row)
    return pd.DataFrame(rows)
