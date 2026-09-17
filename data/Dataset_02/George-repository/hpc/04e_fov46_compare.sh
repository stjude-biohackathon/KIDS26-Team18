#!/bin/bash
#BSUB -P spatial
#BSUB -J f46_cmp
#BSUB -q superdome
#BSUB -n 2
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4000]"
#BSUB -oo logs/f46_cmp_%J.log
#BSUB -eo logs/f46_cmp_%J.err
set -euo pipefail
module load conda3/202402
ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"
cd "$ROOT"
echo "Start: $(date) Host: $(hostname) Python: $PYTHON"
"$PYTHON" "scripts/04e_compare_approaches.py"
echo "End: $(date)"
