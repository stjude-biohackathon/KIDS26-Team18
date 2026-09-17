#!/bin/bash
#BSUB -P spatial
#BSUB -J cosmx_dl
#BSUB -q superdome
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2000]"
#BSUB -oo logs/download_%J.log
#BSUB -eo logs/download_%J.err
set -euo pipefail
module load conda3/202402
ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"
cd "$ROOT"
"$PYTHON" scripts/00_download_cosmx.py
