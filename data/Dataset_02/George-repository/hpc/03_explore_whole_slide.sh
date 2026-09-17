#!/bin/bash
#BSUB -P spatial
#BSUB -J cosmx_qc
#BSUB -q superdome
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6000]"
#BSUB -oo logs/explore_%J.log
#BSUB -eo logs/explore_%J.err
set -euo pipefail
module load conda3/202402
ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"
cd "$ROOT"
"$PYTHON" scripts/03_explore_whole_slide.py
