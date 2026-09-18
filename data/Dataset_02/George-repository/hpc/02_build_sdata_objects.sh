#!/bin/bash
#BSUB -P spatial
#BSUB -J cosmx_sdata
#BSUB -q superdome
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -oo logs/sdata_%J.log
#BSUB -eo logs/sdata_%J.err

set -euo pipefail

ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"

cd "$ROOT"
mkdir -p logs

echo "============================================================"
echo "Build SpatialData from existing raw AnnData + GEO spatial files"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "PWD: $(pwd)"
echo "Python: $PYTHON"
"$PYTHON" --version

echo "Running scripts/02_build_sdata.py"
"$PYTHON" scripts/02_build_sdata.py

echo "Finished: $(date)"
echo "============================================================"
