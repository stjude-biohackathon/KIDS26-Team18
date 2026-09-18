from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

def preprocess(adata, n_pcs=30, n_neighbors=15, seed=42):
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata); adata.raw=adata.copy()
    sc.pp.scale(adata, max_value=10)
    n_pcs=min(n_pcs, adata.n_obs-1, adata.n_vars-1)
    sc.tl.pca(adata,n_comps=n_pcs,svd_solver="arpack",random_state=seed)
    sc.pp.neighbors(adata,n_neighbors=n_neighbors,n_pcs=n_pcs,random_state=seed)
    sc.tl.umap(adata,random_state=seed)
    return n_pcs

def resolution_sweep(adata, resolutions=(0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2), seed=42):
    rows=[]
    for r in resolutions:
        key=f"leiden_r{r:g}"
        sc.tl.leiden(adata,resolution=r,random_state=seed,key_added=key)
        vc=adata.obs[key].value_counts()
        rows.append({"resolution":r,"key":key,"n_clusters":len(vc),"min_cluster_size":int(vc.min()),"median_cluster_size":float(vc.median()),"tiny_cluster_fraction":float(vc[vc<max(20,int(.005*adata.n_obs))].sum()/adata.n_obs)})
    df=pd.DataFrame(rows)
    ari=[np.nan]; nmi=[np.nan]
    for a,b in zip(df.key[:-1],df.key[1:]):
        ari.append(adjusted_rand_score(adata.obs[a],adata.obs[b])); nmi.append(normalized_mutual_info_score(adata.obs[a],adata.obs[b]))
    df["ARI_vs_previous"]=ari; df["NMI_vs_previous"]=nmi
    return df

def suggest_resolution(summary, min_ari=.80, max_tiny_fraction=.05):
    # Heuristic shortlist, not a biological ground truth: prefer stable adjacent solutions
    # without a large tiny-cluster burden; choose the lowest resolution in the stable plateau.
    x=summary[(summary.ARI_vs_previous>=min_ari)&(summary.tiny_cluster_fraction<=max_tiny_fraction)].copy()
    if len(x): return float(x.iloc[0].resolution)
    y=summary.sort_values(["tiny_cluster_fraction","ARI_vs_previous"],ascending=[True,False])
    return float(y.iloc[0].resolution)
