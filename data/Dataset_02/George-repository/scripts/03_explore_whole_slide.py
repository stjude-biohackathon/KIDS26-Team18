#!/usr/bin/env python
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from src.qc import add_qc_metrics,recommend_thresholds,apply_qc,threshold_sweep
from src.clustering import preprocess,resolution_sweep,suggest_resolution
from src.plotting import spatial_scatter,qc_by_fov_plots
IN=ROOT/"data/processed/GSM9046088_CosMx_raw.h5ad"; OUT=ROOT/"results/whole_slide"; FIG=OUT/"figures"; TAB=OUT/"tables"; FIG.mkdir(parents=True,exist_ok=True); TAB.mkdir(parents=True,exist_ok=True)
a=sc.read_h5ad(IN); add_qc_metrics(a)
# QC visualizations
for c in ["total_counts","n_genes_by_counts","control_fraction"]:
    fig,ax=plt.subplots(figsize=(7,4)); sns.histplot(a.obs[c],bins=60,ax=ax); fig.tight_layout(); fig.savefig(FIG/f"qc_{c}.png",dpi=250); plt.close(fig)
if "fov" in a.obs: qc_by_fov_plots(a,FIG)
spatial_scatter(a,"fov",FIG/"spatial_FOV.png",size=3,title="Whole slide: cells colored by FOV"); plt.close("all")
# Candidate threshold grid + robust recommendation
sweep=threshold_sweep(a); sweep.to_csv(TAB/"qc_threshold_sweep.csv",index=False)
thr=recommend_thresholds(a,nmads=3.0)
(TAB/"recommended_qc_thresholds.json").write_text(json.dumps(thr,indent=2))
print("Recommended QC thresholds:",thr); apply_qc(a,thr,use_per_fov_flags=True)
if "fov" in a.obs:
    ret=a.obs.groupby("fov",observed=True).qc_pass.agg(["count","sum","mean"]); ret.to_csv(TAB/"qc_retention_by_fov.csv")
# controls retained in raw object but excluded from biological clustering
bio=~a.var.is_control.to_numpy(); adata=a[a.obs.qc_pass,bio].copy(); sc.pp.filter_genes(adata,min_cells=1); adata.layers["counts"]=adata.X.copy()
preprocess(adata,n_pcs=30,n_neighbors=15,seed=42)
res=resolution_sweep(adata); res.to_csv(TAB/"leiden_resolution_sweep.csv",index=False); selected=suggest_resolution(res)
print("Heuristic resolution shortlist selection:",selected,"(validate with markers/spatial plots)")
selkey=f"leiden_r{selected:g}"; adata.obs["leiden"]=adata.obs[selkey].copy(); adata.uns["selected_leiden_resolution"]=selected
# UMAP multi-resolution
keys=res.key.tolist(); sc.pl.umap(adata,color=keys,ncols=3,show=False); plt.gcf().savefig(FIG/"umap_resolution_sweep.png",dpi=250,bbox_inches="tight"); plt.close("all")
spatial_scatter(adata,"leiden",FIG/"spatial_leiden.png",size=3); plt.close("all")
# quantify FOV domination of clusters
if "fov" in adata.obs:
    ct=pd.crosstab(adata.obs.leiden,adata.obs.fov,normalize="index"); ct.to_csv(TAB/"cluster_by_FOV_fraction.csv")
    fig,ax=plt.subplots(figsize=(12,6)); sns.heatmap(ct,cmap="viridis",ax=ax); ax.set_title("Fraction of each Leiden cluster contributed by each FOV"); fig.tight_layout(); fig.savefig(FIG/"cluster_by_FOV_heatmap.png",dpi=250); plt.close(fig)
sc.tl.rank_genes_groups(adata,"leiden",method="wilcoxon",use_raw=True,pts=True)
markers=sc.get.rank_genes_groups_df(adata,group=None); markers.to_csv(TAB/"leiden_markers.csv",index=False)
out=OUT/"GSM9046088_whole_slide_processed.h5ad"; adata.write_h5ad(out,compression="gzip"); a.write_h5ad(OUT/"GSM9046088_raw_with_qc_flags.h5ad",compression="gzip")
print("saved",out)
