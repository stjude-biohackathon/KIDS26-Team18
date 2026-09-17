# BioHackathon 2026 — Dataset_02 CosMx revised workflow

Reproducible analysis workflow for **GSM9046088 — Lung Adenocarcinoma TMA1 (NanoString/Bruker CosMx SMI)** from GEO series **GSE299786**. The project reconstructs both cell-level `AnnData` and molecule/segmentation-aware `SpatialData`, explores the whole TMA, focuses on the region labeled **FOV46**, compares two cell-typing strategies, and validates the annotations against the published/metadata reference annotation when available.

> **Terminology note.** In the source dataset, labels called `FOV` can correspond to TMA-core-level regions in the reconstructed layout. This can be confusing because a physical TMA core itself contains connected instrument fields of view. In this repository, `FOV46` follows the identifier used by the supplied metadata/code.

## Data source and biological context

The source study compares imaging-based single-cell-resolution spatial transcriptomics platforms on FFPE tumor tissue microarrays. This repository analyzes the CosMx lung adenocarcinoma TMA1 sample (`GSM9046088`).

- GEO: GSE299786 / GSM9046088
- Platform: NanoString/Bruker CosMx SMI
- Tissue: lung adenocarcinoma, FFPE TMA
- Publication: *Nature Communications* (2025), DOI page linked in the project presentation

The pipeline downloads the five CosMx flat files programmatically:

| File | Role |
|---|---|
| `*_exprMat_file.csv.gz` | Cell × feature transcript-count matrix |
| `*_metadata_file.csv.gz` | Cell metadata, QC-related fields, coordinates, and supplied annotations |
| `*_fov_positions_file.csv.gz` | Spatial placement of regions/FOVs on the slide |
| `*_polygons.csv.gz` | Cell-segmentation polygon vertices |
| `*_tx_file.csv.gz` | Individual transcript molecules with spatial coordinates and target identity |

Manual GEO downloading is therefore not required.

## Analysis overview

```text
GEO GSM9046088
      |
      v
00_download_cosmx.py
      |
      +---------------------------+
      |                           |
      v                           v
01_build_anndata.py        02_build_sdata.py
      |                           |
raw AnnData                 SpatialData
      |                    table + polygons + transcripts
      |                           |
      +-------------+-------------+
                    |
                    v
          03 whole-TMA exploration
                    |
          tissue/core heterogeneity
                    |
                    v
          extract raw FOV46 subset
                    |
                    v
              04a QC filtering
                    |
                    v
        04b Leiden sweep 0.1–1.0
                    |
                 STOP
                    |
             manual inspection
       UMAP + markers + spatial maps
                    |
                    v
        choose Leiden resolution
       + manual cluster annotation
                    |
          +---------+---------+
          |                   |
          v                   v
  04c Approach I       04d Approach II
  cluster-based        rule-based typing
          |                   |
          +---------+---------+
                    |
                    v
             04e comparison
        I vs II vs reference
                    |
                    v
          04f transcript overlap
```

A key design principle is that **raw AnnData, processed AnnData, and SpatialData remain separate**. FOV46 is extracted from the raw object rather than from an already normalized/clustering-processed object.

## Environment

Create the environment once:

```bash
module load conda3/202402
conda env create -f env/environment.yml
conda activate spatialdata
```

On the St. Jude HPC, batch jobs deliberately call the desired environment Python explicitly because non-interactive LSF shells may not initialize `conda activate` consistently:

```bash
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"
```

Update `ROOT` and `PYTHON` in `hpc/*.sh` if the repository or Conda environment is installed elsewhere.

## Run from scratch on HPC

From the repository root:

```bash
mkdir -p logs

bsub < hpc/00_download_raw.sh
bsub < hpc/01_build_anndata_objects.sh
bsub < hpc/02_build_sdata_objects.sh
bsub < hpc/03_explore_whole_slide.sh
```

After whole-TMA exploration, extract FOV46 if needed and run the staged FOV46 workflow. The exact extraction script name may depend on the checked-out revision; the important point is that the subset is created from the **raw AnnData** and the corresponding SpatialData region is also retained.

Then run QC and the Leiden sweep:

```bash
bsub < hpc/04a_fov46_qc.sh
bsub < hpc/04b_fov46_leiden_sweep.sh
```

### Mandatory manual checkpoint after 04b

**Do not automatically continue to Approach I.** Inspect the outputs from resolutions `0.1` through `1.0`, including:

- multi-resolution UMAPs;
- number and size of clusters;
- ARI/NMI stability between adjacent resolutions;
- cluster-transition tables (clustree-like inspection);
- top differential genes for each resolution;
- known-marker dotplots;
- spatial distribution of each Leiden solution.

The current analysis presentation uses **resolution 0.6** for FOV46, but this value should remain an explicit, inspectable analysis choice rather than an automatic universal optimum.

Record the selected resolution and cluster-to-cell-type mapping in:

```text
config/fov46_annotation.json
```

Then continue:

```bash
bsub < hpc/04c_fov46_approach_i.sh
bsub < hpc/04d_fov46_approach_ii.sh
bsub < hpc/04e_fov46_compare.sh
bsub < hpc/04f_fov46_spatial_overlap.sh
```

## QC strategy

CosMx is an imaging-based targeted spatial-transcriptomics assay, so the workflow does **not** blindly reuse conventional droplet-scRNA-seq cutoffs. QC is distribution-aware and considers multiple aspects of each segmented cell:

- `total_counts`: total detected transcript molecules;
- `n_genes_by_counts`: detected-feature complexity;
- control/negative-probe fraction;
- cell segmentation area when available;
- spatial distribution of QC failures.

The QC stage uses robust MAD-based thresholds together with conservative information floors and writes a threshold-sensitivity grid for inspection. The exact thresholds are data-derived and saved to:

```text
results/FOV46/tables/recommended_qc_thresholds.json
results/FOV46/tables/qc_threshold_sweep.csv
```

Additional diagnostics include counts removed by each criterion, failure combinations, before/after summaries, metric distributions with threshold lines, counts-vs-genes plots, and spatial maps of QC-pass/QC-fail cells. The goal is to remove low-information or segmentation-outlier cells without erasing a coherent biological region.

## Whole-TMA exploration

Whole-TMA Scanpy analysis showed that expression-space clustering is strongly associated with separated tissue/core regions. Importantly, the Scanpy neighbor graph used here is expression-based; physical X/Y distance is not supplied to Leiden clustering. Therefore, agreement between Leiden clusters and tissue regions indicates expression differences correlated with tissue/core origin rather than clustering caused directly by spatial distance.

This motivates a hierarchical strategy:

1. characterize whole-TMA/core heterogeneity;
2. focus on a selected region (`FOV46`);
3. resolve within-region cell identities.

## Leiden resolution selection

`04b` evaluates Leiden resolutions from **0.1 to 1.0** on the same PCA/kNN graph. Resolution is not selected solely because adjacent solutions have a high ARI/NMI: a very coarse two-cluster solution can be perfectly stable yet biologically under-resolved.

Resolution selection therefore combines:

1. **Structural stability** — cluster number/size, tiny-cluster burden, ARI/NMI, transition tables;
2. **Biological separation** — top differential genes and known marker programs;
3. **Spatial plausibility** — coordinate and segmentation-based spatial localization.

This is analogous in spirit to inspecting multiple Seurat resolutions with `clustree`, followed by marker-based biological interpretation.

## Approach I — cluster-based annotation

`04c_approach_i_annotation.py` uses the manually selected Leiden solution and produces:

- all cluster-vs-rest differential genes;
- top-10 marker genes per Leiden cluster;
- original-style `rank_genes_groups` top-gene panels;
- top-marker dotplot and heatmap;
- predefined marker-panel dotplot;
- selected-resolution UMAP and spatial maps;
- manually assigned `cell_type_scanpy` labels.

The cluster-to-cell-type mapping is intentionally kept in `config/fov46_annotation.json` so that annotation remains reviewable rather than hidden inside an automatic heuristic.

## Approach II — rule-based annotation

`04d_approach_ii_rule_based.py` applies marker-panel rules directly to the raw-count layer. A candidate type is positive when **at least three available marker genes have raw count > 0**. A candidate panel with fewer than three assayed markers is disabled rather than silently weakening the rule.

Candidate programs include epithelial, AT1/AT2, tumor epithelial, T cell, B cell, plasma cell, myeloid/macrophage, endothelial, fibroblast, and mural programs, subject to marker availability in the CosMx panel.

Cells matching:

- exactly one candidate panel → assigned that cell type;
- more than one candidate panel → `Ambiguous`;
- no candidate panel → `Unassigned`.

Approach II outputs marker detection/availability tables and UMAPs for cell type and assignment status. Approach I and II are also shown side-by-side on the **same UMAP coordinates**.

## Comparison with reference annotation

When `final_CT` is present in metadata, it is treated as a **reference annotation**, not an error-free ground truth.

`04e_compare_approaches.py` performs three categorical comparisons:

```text
Approach I  × Approach II
final_CT    × Approach I
final_CT    × Approach II
```

For each comparison the workflow saves:

- raw contingency table;
- row-normalized concordance table;
- concordance heatmap;
- 100% stacked bar plot.

The row-normalized matrix estimates the distribution of column labels within each row label. This is a categorical concordance analysis, not a Pearson/Spearman correlation matrix. Exact-string agreement is reported separately because different annotation granularities (for example, `Myeloid` vs `Macrophage`) can be biologically related without having identical labels.

## Unified UMAP and spatial visualization

Approach I, Approach II, and the reference annotation use a common plotting theme. Each annotation can be visualized in three complementary spatial representations:

1. **Coordinate map** — colored cell centers in global spatial coordinates;
2. **Filled segmentation map** — actual cell polygons filled by annotation color;
3. **Outline segmentation map** — black background, white cell boundaries, transparent polygon interiors, with annotation-colored cell-center points.

This makes it possible to distinguish expression-space organization (UMAP) from true tissue geometry and from actual cell-segmentation morphology.

## Transcript-level SpatialData overlays

`04f_fov46_spatial_overlap.py` uses the FOV46 SpatialData Zarr store to overlay:

```text
cell-level gene expression
        +
cell segmentation polygons
        +
individual transcript molecules
```

Representative markers include genes such as `EPCAM`, `CD3D`, `MS4A1`, `C1QA`, `PECAM1`, `COL1A1`, and `RGS5` when present in the panel. The plotting code uses explicit SpatialData element names (`element="cell_boundaries"` and `element="transcripts"`) for compatibility with the installed `spatialdata-plot` API.

The five GEO flat files do not provide the full morphology-image stack used by the instrument, so this repository does not fabricate DAPI/PanCK/etc. image backgrounds when those images are unavailable.

## Project structure

```text
.
├── config/
│   ├── markers.py
│   └── fov46_annotation.json
├── data/
│   ├── raw/
│   ├── processed/
│   └── spatial/
├── env/
│   └── environment.yml
├── hpc/
│   ├── 00_download_raw.sh
│   ├── 01_build_anndata_objects.sh
│   ├── 02_build_sdata_objects.sh
│   ├── 03_explore_whole_slide.sh
│   ├── 04a_fov46_qc.sh
│   ├── 04b_fov46_leiden_sweep.sh
│   ├── 04c_fov46_approach_i.sh
│   ├── 04d_fov46_approach_ii.sh
│   ├── 04e_fov46_compare.sh
│   └── 04f_fov46_spatial_overlap.sh
├── scripts/
│   ├── 00_download_cosmx.py
│   ├── 01_build_anndata.py
│   ├── 02_build_sdata.py
│   ├── 03_explore_whole_slide.py
│   ├── 04a_qc_fov46.py
│   ├── 04b_leiden_sweep_fov46.py
│   ├── 04c_approach_i_annotation.py
│   ├── 04d_approach_ii_rule_based.py
│   ├── 04e_compare_approaches.py
│   └── 04f_spatial_overlap.py
├── src/
│   ├── qc.py
│   ├── clustering.py
│   ├── typing.py
│   ├── plotting.py
│   └── spatial_plot_legacy.py
├── results/
│   ├── whole_slide/
│   └── FOV46/
├── logs/
└── README.md
```

The exact checked-out repository may contain additional helper scripts or notebooks; the staged workflow above describes the current analysis design.

## Key outputs

Important outputs to inspect include:

```text
results/FOV46/tables/
├── recommended_qc_thresholds.json
├── qc_threshold_sweep.csv
├── leiden_resolution_sweep.csv
├── leiden_top10_markers_*.csv
├── rule_marker_detection_rates.csv
├── rule_marker_availability.csv
├── approach_i_vs_ii_counts.csv
├── approach_i_vs_ii_row_fraction.csv
├── reference_vs_approach_i_*.csv
├── reference_vs_approach_ii_*.csv
└── typing_summary.csv
```

and the corresponding figures under:

```text
results/FOV46/figures/
├── qc/
├── leiden_sweep/
├── approach_i/
├── approach_ii/
├── comparison/
└── spatial_overlap/
```

## Notes on reproducibility

- Raw GEO files are downloaded programmatically and are not manually edited.
- Raw AnnData is created once and retained as the provenance object.
- SpatialData construction is a separate step from AnnData construction.
- FOV46 is extracted from raw data before FOV-specific QC and preprocessing.
- Leiden resolution selection contains an explicit manual checkpoint.
- Approach II uses the raw-count layer, not normalized/log-transformed values.
- Spatial annotation figures use the same annotation categories across UMAP, coordinate, and segmentation views whenever possible.
- Batch scripts explicitly invoke the intended Conda-environment Python on the St. Jude HPC.

## References

- GEO series: GSE299786 — Comparison of imaging-based single-cell-resolution spatial transcriptomics profiling platforms using FFPE tumor samples.
- Sample: GSM9046088 — Lung Adenocarcinoma TMA1, CosMx.
- Publication associated with the dataset: *Nature Communications* (2025), PubMed PMID 41006245.

