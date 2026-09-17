import scanpy as sc
import pandas as pd

SCANPY_CFG = {
    "n_pcs": 30,
    "n_neighbors": 15,
    "leiden_resolution": 0.5,
    "rank_genes_top_n": 20,
    "seed": 42,
}

def leiden_and_rank(adata, scanpy_config=SCANPY_CFG):
    '''
    Perform leiden clustering on a given anndata object and
    rank genes groups.
    Takes adata object and optional scanpy config override.
    '''
    adata_sc = adata.copy()
    
    # Targeted 343-gene panel: use all genes (no HVG subset)
    sc.pp.normalize_total(adata_sc, target_sum=1e4)
    sc.pp.log1p(adata_sc)
    sc.pp.scale(adata_sc, max_value=10)
    
    n_comps = min(scanpy_config["n_pcs"], adata_sc.n_obs - 1, adata_sc.n_vars - 1)
    sc.tl.pca(adata_sc, n_comps=n_comps, random_state=scanpy_config["seed"])
    sc.pp.neighbors(adata_sc, n_neighbors=scanpy_config["n_neighbors"])
    sc.tl.leiden(
        adata_sc,
        resolution=scanpy_config["leiden_resolution"],
        key_added="leiden",
        random_state=scanpy_config["seed"],
    )
    sc.tl.umap(adata_sc, random_state=scanpy_config["seed"])
    
    # Rank Genes Groups
    sc.tl.rank_genes_groups(adata_sc, groupby="leiden", method="wilcoxon", use_raw=False)

    # Add coord_x and coord_y cols
    adata.obs["coord_x"] = adata.obsm["spatial"][:, 0]
    adata.obs["coord_y"] = adata.obsm["spatial"][:, 1]
    
    return adata_sc

def get_topn_markers(adata, top_n=10):
    '''
    Return top n markers from culstered anndata.
    '''
    clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))

    marker_tables = []
    for cluster in clusters:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df = df.head(top_n).copy()
        df.insert(0, "cluster", cluster)
        marker_tables.append(df)
    
    return pd.concat(marker_tables, ignore_index=True)


def annotate_clusters(adata, cluster_defs):
    adata.obs["cell_type_scanpy"] = (
        adata.obs["leiden"]
        .astype(str)
        .map(cluster_defs)
        .fillna("Unassigned")
    )

    adata.obs["cell_type_scanpy"] = adata.obs["cell_type_scanpy"].astype("category")
    return adata