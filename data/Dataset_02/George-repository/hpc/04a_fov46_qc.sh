#!/bin/bash

#BSUB -P spatial
#BSUB -J f46_qc
#BSUB -q superdome
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6000]"
#BSUB -oo logs/f46_qc_%J.log
#BSUB -eo logs/f46_qc_%J.err

set -euo pipefail

ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"

PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"

cd "$ROOT"

echo "============================================================"
echo "FOV46 QC"
echo "Start:   $(date)"
echo "Host:    $(hostname)"
echo "Python:  $PYTHON"
echo "Root:    $ROOT"
echo "============================================================"

"$PYTHON" --version

"$PYTHON" scripts/04a_qc_fov46.py

echo "============================================================"
echo "FOV46 QC completed"
echo "End: $(date)"
echo "============================================================"