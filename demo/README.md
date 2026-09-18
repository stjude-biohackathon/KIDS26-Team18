# Cell typing comparison demo

Interactive Streamlit app comparing **Scanpy Leiden + manual annotation** vs **rule-based marker panels**.

Supported datasets:

| ID | Tissue / platform | Rule engine |
| --- | --- | --- |
| Dataset 01 | GSE250346 lung Xenium TMA5 | `demo/markers.py` |
| Dataset 04 | GSE282026 kidney CosMx | `celltyping_config.yaml` + `celltyping_helpers.py` |
| Dataset 05 | GSE301280 skin Xenium | Precomputed `hier_label` vs author `cell_type` |

## Prerequisites

- Conda env `spatialdata` (Python ≥3.11, spatialdata ≥0.7.3)
- Processed h5ad files:
  - **Dataset 01:** `data/Dataset_01/processed_data/one_sample_clustered.h5ad` (required)
  - **Dataset 01 optional:** `cell_typed_TMA5_subset.h5ad` — published `final_CT`
  - **Dataset 04:** `data/Dataset_04/processed_data/scanpy_cluster_cellbased.h5ad` (required)
  - **Dataset 05:** `data/Dataset_05/processed_data/adata_annotated.h5ad` (required)

## Run

From the repo root:

```bash
conda activate spatialdata
streamlit run demo/app.py
```

Opens at [http://localhost:8501](http://localhost:8501).

## Views

### Main spatial comparison

- **Side-by-side** — Scanpy | Rule-based by default
- **Dataset 01 only:** optional published (`final_CT`) panel via sidebar checkbox
- **Agreement** — where the two demo methods disagree
- Respects spatial crop and filters from the sidebar
- **Dataset 04:** spatial crop defaults to the notebook crop from `celltyping_config.yaml`

### Clustering tab (notebook-style)

Separate from the spatial comparison — uses the **full sample**, not the spatial crop:

- **UMAP** — Leiden clusters vs manual Scanpy labels (`sc.pl.umap`, same as notebook)
- **Rank-genes dotplot** — top DE genes per Leiden cluster
- **Cluster annotations** table — manual leiden → label mapping from notebook (D01) or config (D04)

## Demo flow (2–3 min)

1. **Clustering tab** — show Leiden → manual annotation workflow (UMAP + dotplot)
2. **Side-by-side** — same tissue, annotation views on spatial coordinates
3. **Agreement** — red = methods disagree
4. **Rule tuning** — enable “Recompute rule labels live”, drag min hits / threshold
5. **Type counts tab** — bar charts, `nCount` violin by cell type, and Unassigned QC
6. **Metrics tab** — Cohen's κ and agreement stats for the current crop

## Files

| File | Role |
| --- | --- |
| `app.py` | Streamlit UI |
| `datasets.py` | Dataset registry, loaders, rule assignment dispatch |
| `cluster_plots.py` | Notebook-style UMAP + rank-genes dotplot |
| `markers.py` | Dataset 01 marker panels + tunable rule assignment |
| `metrics.py` | Accuracy / κ / crosstabs |
| `plots.py` | Plotly spatial and summary charts |
