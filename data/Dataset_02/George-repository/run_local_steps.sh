#!/bin/bash
# Dataset_02 — GSM9046088 CosMx revised workflow
#
# This convenience runner mirrors the current staged workflow.
# IMPORTANT: execution intentionally STOPS after the FOV46 Leiden sweep so
# that the resolution and Approach-I cluster annotations can be inspected
# manually before downstream typing.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-/home/jqu/.conda/envs/spatialdata/bin/python}"

cd "$ROOT"
mkdir -p logs

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: Python executable not found or not executable: $PYTHON"
    exit 1
fi

banner() {
    echo
    echo "======================================================================"
    echo "$1"
    echo "======================================================================"
}

run_python() {
    local script="$1"
    if [[ ! -f "scripts/$script" ]]; then
        echo "ERROR: Missing scripts/$script"
        exit 1
    fi
    banner "RUNNING: scripts/$script"
    "$PYTHON" "scripts/$script"
}

manual_checkpoint() {
    banner "MANUAL CHECKPOINT — STOP AFTER 04b"
    cat <<'EOF'
Inspect the Leiden sweep before continuing:

  results/FOV46/figures/leiden_sweep/
  results/FOV46/tables/leiden_resolution_sweep.csv
  results/FOV46/tables/transition_*.csv

Compare resolutions 0.1–1.0 using:
  - UMAP structure
  - cluster sizes
  - ARI / NMI stability
  - adjacent-resolution transitions
  - marker dotplots / top markers
  - spatial localization

Then edit:
  config/fov46_annotation.json

Set the selected resolution and complete cluster_annotations.

When ready:
  bash run_local_steps.sh fov46-post
EOF
}

download_raw()       { run_python "00_download_cosmx.py"; }
build_anndata()      { run_python "01_build_anndata.py"; }
build_sdata()        { run_python "02_build_sdata.py"; }
whole_slide()        { run_python "03_explore_whole_slide.py"; }
extract_fov46()      { run_python "04_extract_fov46.py"; }
qc_fov46()           { run_python "04a_qc_fov46.py"; }
leiden_sweep_fov46() { run_python "04b_leiden_sweep_fov46.py"; }
approach_i()         { run_python "04c_approach_i_annotation.py"; }
approach_ii()        { run_python "04d_approach_ii_rule_based.py"; }
compare_methods()    { run_python "04e_compare_approaches.py"; }
spatial_overlap()    { run_python "04f_spatial_overlap.py"; }

build_workflow() {
    banner "00 — PROGRAMMATIC GEO DOWNLOAD"
    download_raw
    banner "01 — BUILD RAW ANNDATA"
    build_anndata
    banner "02 — BUILD SPATIALDATA"
    build_sdata
    banner "03 — WHOLE-SLIDE EXPLORATION"
    whole_slide
}

fov46_pre_workflow() {
    banner "04 — EXTRACT FOV46 FROM RAW OBJECTS"
    extract_fov46
    banner "04a — FOV46 QC"
    qc_fov46
    banner "04b — LEIDEN SWEEP 0.1–1.0"
    leiden_sweep_fov46
    manual_checkpoint
}

fov46_post_workflow() {
    banner "04c — APPROACH I: CLUSTER-BASED MANUAL ANNOTATION"
    approach_i

    banner "04d — APPROACH II: RULE-BASED ANNOTATION"
    echo "Rule: >= 3 available marker genes have raw count > 0."
    approach_ii

    banner "04e — APPROACH I / II / REFERENCE COMPARISON"
    compare_methods

    banner "04f — SPATIAL GENE / SEGMENTATION / TRANSCRIPT OVERLAYS"
    spatial_overlap

    banner "FOV46 POST-INSPECTION WORKFLOW COMPLETE"
    cat <<'EOF'
Current visualization outputs include:

Approach I / Approach II / reference where available:
  - UMAP with a common theme and visible axes
  - coordinate-based spatial annotation
  - filled cell-segmentation polygons
  - black-background spatial plot with white cell boundaries and
    annotation-colored cell-center points

Approach I:
  - selected-resolution top DE genes
  - cluster-vs-rest top-marker panels
  - top-marker heatmap / dotplot
  - known marker-panel dotplot

Approach II:
  - rule-based UMAP
  - Approach-I-vs-II side-by-side UMAP
  - Assigned / Ambiguous / Unassigned UMAP
  - rule uses >=3 detected marker genes with raw count >0

Comparison:
  - Approach I × Approach II contingency counts
  - row-normalized concordance tables
  - concordance heatmaps
  - 100% stacked bars
  - final_CT/reference × Approach I
  - final_CT/reference × Approach II

SpatialData overlays:
  - cell-level gene expression
  - cell segmentation
  - individual transcript molecules
EOF
}

usage() {
    cat <<EOF
Usage:
  bash $(basename "$0") <mode>

Recommended from-scratch workflow:

  bash $(basename "$0") build
  bash $(basename "$0") fov46-pre

  # STOP and inspect 04b outputs.
  # Edit config/fov46_annotation.json.

  bash $(basename "$0") fov46-post

Modes:
  build          Run 00–03
  fov46-pre      Run 04, 04a, 04b, then STOP for manual inspection
  fov46-post     Run 04c–04f after manual annotation
  download       Run 00 only
  objects        Run 01 + 02
  whole-slide    Run 03 only
  extract-fov46  Run 04 only
  fov46-qc       Run 04a only
  fov46-sweep    Run 04b only, then print manual checkpoint
  approach-i     Run 04c only
  approach-ii    Run 04d only
  compare        Run 04e only
  spatial        Run 04f only
  help           Show this message

Python defaults to:
  /home/jqu/.conda/envs/spatialdata/bin/python

Override if needed:
  PYTHON=/path/to/env/bin/python bash $(basename "$0") <mode>
EOF
}

banner "DATASET_02 COSMX WORKFLOW"
echo "Project root : $ROOT"
echo "Python       : $PYTHON"
"$PYTHON" --version

mode="${1:-help}"

case "$mode" in
    build)
        build_workflow
        ;;
    fov46-pre)
        fov46_pre_workflow
        ;;
    fov46-post)
        fov46_post_workflow
        ;;
    download)
        download_raw
        ;;
    objects)
        build_anndata
        build_sdata
        ;;
    whole-slide)
        whole_slide
        ;;
    extract-fov46)
        extract_fov46
        ;;
    fov46-qc)
        qc_fov46
        ;;
    fov46-sweep)
        leiden_sweep_fov46
        manual_checkpoint
        ;;
    approach-i)
        approach_i
        ;;
    approach-ii)
        approach_ii
        ;;
    compare)
        compare_methods
        ;;
    spatial)
        spatial_overlap
        ;;
    help|-h|--help)
        usage
        ;;
    *)
        echo "ERROR: Unknown mode: $mode"
        echo
        usage
        exit 2
        ;;
esac
