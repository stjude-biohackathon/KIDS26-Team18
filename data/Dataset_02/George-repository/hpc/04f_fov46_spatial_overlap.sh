#!/bin/bash
#BSUB -P spatial
#BSUB -J f46_spatial
#BSUB -q superdome
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -oo logs/f46_spatial_%J.log
#BSUB -eo logs/f46_spatial_%J.err


set -euo pipefail

ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"

PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"

cd "$ROOT"

echo "============================================================"
echo "FOV46 spatial overlap"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Python: $PYTHON"
echo "============================================================"

"$PYTHON" scripts/04f_spatial_overlap.py

echo "============================================================"
echo "Completed successfully"
echo "End: $(date)"
echo "============================================================"