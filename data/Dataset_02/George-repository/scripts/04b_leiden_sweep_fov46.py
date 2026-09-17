#!/usr/bin/env python
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from src.clustering import preprocess
from src.plotting import spatial_scatter
from config.markers import MARKER_PANEL

IN=ROOT/'results/FOV46/GSM9046088_FOV46_qc_filtered.h5ad'
OUT=ROOT/'results/FOV46'; FIG=OUT/'figures/leiden_sweep'; TAB=OUT/'tables'; FIG.mkdir(parents=True,exist_ok=True); TAB.mkdir(parents=True,exist_ok=True)
a=sc.read_h5ad(IN)
preprocess(a)
resolutions=[round(x/10,1) for x in range(1,11)]
rows=[]
prev=None
for r in resolutions:
    key=f'leiden_r{r:g}'
    sc.tl.leiden(a,resolution=r,random_state=42,key_added=key)
    vc=a.obs[key].value_counts()
    row={'resolution':r,'key':key,'n_clusters':len(vc),'min_cluster_size':int(vc.min()),'median_cluster_size':float(vc.median()),'tiny_cluster_fraction':float(vc[vc<max(20,int(.005*a.n_obs))].sum()/a.n_obs)}
    row['ARI_vs_previous']=float('nan') if prev is None else adjusted_rand_score(a.obs[prev],a.obs[key])
    row['NMI_vs_previous']=float('nan') if prev is None else normalized_mutual_info_score(a.obs[prev],a.obs[key])
    rows.append(row); prev=key
summary=pd.DataFrame(rows); summary.to_csv(TAB/'leiden_resolution_sweep.csv',index=False)
keys=summary['key'].tolist()
sc.pl.umap(a,color=keys,ncols=3,show=False)
plt.gcf().savefig(FIG/'umap_resolution_sweep_0.1_to_1.0.png',dpi=250,bbox_inches='tight'); plt.close('all')

present={ct:[g for g in genes if g in a.raw.var_names] for ct,genes in MARKER_PANEL.items()}
present={ct:g for ct,g in present.items() if g}
for r,key in zip(resolutions,keys):
    spatial_scatter(a,key,FIG/f'spatial_{key}.png',size=8); plt.close('all')
    if present:
        dp=sc.pl.dotplot(a,var_names=present,groupby=key,use_raw=True,standard_scale='var',dendrogram=True,show=False,return_fig=True)
        dp.savefig(FIG/f'marker_dotplot_{key}.png',dpi=220); plt.close('all')
    nxt=None
    if r < 1.0: nxt=f'leiden_r{round(r+0.1,1):g}'
    if nxt and nxt in a.obs:
        pd.crosstab(a.obs[key],a.obs[nxt],normalize='index').to_csv(TAB/f'transition_{key}_to_{nxt}.csv')
a.write_h5ad(OUT/'GSM9046088_FOV46_leiden_sweep.h5ad',compression='gzip')
print(summary.to_string(index=False))
print('\nPAUSE HERE. Inspect UMAP, marker dotplots, spatial maps, cluster sizes, ARI/NMI and transition tables.')
print('Then edit config/fov46_annotation.json and run step 04c.')
