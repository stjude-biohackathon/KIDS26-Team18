'''
Utility code for clustering approach.
'''
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
    print("Leiden...")
    sc.tl.leiden(
        adata_sc,
        resolution=scanpy_config["leiden_resolution"],
        key_added="leiden",
        random_state=scanpy_config["seed"],
    )
    sc.tl.umap(adata_sc, random_state=scanpy_config["seed"])
    
    # Rank Genes Groups
    print("Rank genes group")
    sc.tl.rank_genes_groups(adata_sc, groupby="leiden", method="wilcoxon", use_raw=False)

    # Add coord_x and coord_y cols
    adata_sc.obs["coord_x"] = adata_sc.obsm["spatial"][:, 0]
    adata_sc.obs["coord_y"] = adata_sc.obsm["spatial"][:, 1]
    
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

def get_cluster_annotations(marker_dict):
    '''
    Currently this just prints an LLM prompt for assigning cell types to each cluster.
    To be replaced by a standardized method (most likely).
    '''
    # Display LLM prompt for annotating clusters.
    llm_prompt = "Paste the following prompt into your LLM of choice for updating the cluster annotations dictionary:\nCreate a python dictionary that defines the cell type for each given cluster based on the available markers. Do not label clusters as Doublets, Ambiguous, or Low quality unless explicitly instructed. Assign the most likely biological cell type based on dominant lineage markers. The markers are: " + str(marker_dict) +". Use the following format for the python dictionary so that it can be directly pasted into a Jupyter notebook: CLUSTER_ANNOTATIONS = {\n\"0\": \"Epithelial\",\n..."
    
    print(llm_prompt)

def annotate_clusters(adata, cluster_defs):
    '''
    Update cluster cell types in anndata from given cluster definitions.
    '''
    adata.obs["cell_type_scanpy"] = (
        adata.obs["leiden"]
        .astype(str)
        .map(cluster_defs)
        .fillna("Unassigned")
    )

    adata.obs["cell_type_scanpy"] = adata.obs["cell_type_scanpy"].astype("category")
    return adata

