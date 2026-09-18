#!/bin/bash
#BSUB -P spatial
#BSUB -J cosmx_adata
#BSUB -q superdome
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6000]"
#BSUB -oo logs/anndata_%J.log
#BSUB -eo logs/anndata_%J.err

set -euo pipefail

ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"

cd "$ROOT"
mkdir -p logs

echo "============================================================"
echo "Build raw AnnData"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "PWD: $(pwd)"
echo "Python: $PYTHON"
"$PYTHON" --version

echo "Running scripts/01_build_anndata.py"
"$PYTHON" scripts/01_build_anndata.py

echo "Finished: $(date)"
echo "============================================================"
