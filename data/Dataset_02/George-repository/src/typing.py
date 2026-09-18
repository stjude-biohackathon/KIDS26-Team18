from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import sparse

def detection_rate(adata, genes, layer="counts"):
    present=[g for g in genes if g in adata.var_names]; out={}
    Xall=adata.layers[layer] if layer in adata.layers else adata.X
    for g in present:
        j=adata.var_names.get_loc(g); x=Xall[:,j]
        x=x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()
        out[g]=float((x>0).mean())
    return pd.Series(out).sort_values(ascending=False)

def rule_based_all_markers(adata, marker_panel, layer="counts"):
    """Preserve user's rule: a candidate is positive only when ALL available panel markers >0.
    If any requested marker is absent from the assay, that candidate is marked unavailable rather
    than silently relaxing the rule. Multiple positive candidates => Ambiguous.
    """
    Xall=adata.layers[layer] if layer in adata.layers else adata.X
    positives={}; availability=[]
    for ct, genes in marker_panel.items():
        present=[g for g in genes if g in adata.var_names]; missing=[g for g in genes if g not in adata.var_names]
        availability.append({"cell_type":ct,"n_requested":len(genes),"n_present":len(present),"present":",".join(present),"missing":",".join(missing)})
        if missing or not present:
            positives[ct]=np.zeros(adata.n_obs,dtype=bool); continue
        idx=[adata.var_names.get_loc(g) for g in present]; x=Xall[:,idx]
        x=x.toarray() if sparse.issparse(x) else np.asarray(x)
        positives[ct]=(x>0).all(axis=1)
        adata.obs[f"rule_{ct}_all_positive"]=positives[ct]
    mat=np.column_stack([positives[k] for k in marker_panel]); n=mat.sum(axis=1)
    names=np.array(list(marker_panel),dtype=object); labels=np.full(adata.n_obs,"Unassigned",dtype=object)
    one=n==1; labels[one]=names[mat[one].argmax(axis=1)]; labels[n>1]="Ambiguous"
    adata.obs["cell_type_rule"]=pd.Categorical(labels)
    adata.obs["rule_n_positive_panels"]=n
    return pd.DataFrame(availability)
