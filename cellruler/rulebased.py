 '''
Utility code for rule-based approach.
'''
import pandas as pd
import numpy as np

def detection_rate(adata, genes):
    '''
    Fraction of cells with raw count > 0 for each gene in adata.
    '''
    present = [g for g in genes if g in adata.var_names]
    rates = {}
    for gene in present:
        x = adata[:, gene].X
        if hasattr(x, "toarray"):
            x = x.toarray().ravel()
        else:
            x = np.asarray(x).ravel()
        rates[gene] = float((x > 0).mean())
    return pd.Series(rates).sort_values(ascending=False)

def get_top_n_markers(adata, candidate_markers, n=4):
    '''
    Get top n markers for each cell type based on candidate markers.
    Gets max n markers based on what actually appears in the data.
    '''
    final_markers = {}
    
    for cell_type in candidate_markers.keys():
        rates = detection_rate(adata, candidate_markers[cell_type])
        markers = rates.head(n).index.tolist()
    
        print(f"\nSelected {cell_type} markers:")
        display(rates.head(n))
    
        final_markers[cell_type] = markers

    # Remove empty entries.
    final_markers_tmp = {}
    for key in final_markers.keys():
        if len(final_markers[key]) > 0:
            final_markers_tmp[key] = final_markers[key]
    
    final_markers = final_markers_tmp
    return final_markers

def _get_gene_counts(adata, genes):
    '''
    Return (n_cells, n_genes) array of integer transcript counts from .X.
    Pull raw transcript counts for a 4-gene panel from adata.
    '''
    x = adata[:, genes].X
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def assign_rule_labels(adata, markers_dict):
    '''
    Assign labels when all 4 genes in a panel have count > 0.
    1. For each panel, sum the 4 marker-gene counts per cell.
    2. Keep only panels that qualify.
    3. If two panels ties, pcik the one that appears first in the markers_dict.
    4. If no panel qualifies, label as "Unassigned."
    '''
    
    cell_types = list(markers_dict.keys())
    n_cells = adata.n_obs
    n_panels = len(cell_types)

    # Precompute a score matrix for all panels at once.
    # Score matrix: rows = panels, columns = cells (scores[p,i] = total transcript count for panel p in cell i)
    scores = np.zeros((n_panels, n_cells), dtype=float)

    # Fill one row of score matrix per panel
    for p, cell_type in enumerate(cell_types):
        genes = markers_dict[cell_type] # 4-gene panel for cell type
        x = _get_gene_counts(adata, genes) 
        scores[p] = x.sum(axis=1)

    qualifies = scores > 0

    panel_order = np.arange(n_panels).reshape(-1, 1)
    rank = scores + panel_order * 1e-12 # Just breaking small floating point ties

    # Disqualify panels that failed the rule by setting them to inf
    rank = np.where(qualifies, rank, -np.inf)

    # For each cell, pick panel row with highest rank
    best_panel_idx = rank.argmax(axis=0)

    has_hit = qualifies.any(axis=0) # Whether or not cells qualify

    labels = np.full(n_cells, "Unassigned", dtype=object)
    labels[has_hit] = np.array(cell_types, dtype=object)[best_panel_idx[has_hit]] # Set cell type if it has a hit

    return labels