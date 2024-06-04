#!/bin/sh
#SBATCH --gres=gpu:1
#SBATCH -n 8
#SBATCH -e /home/hk-project-p0021769/hgf_ysb1444/logging/log_%j.err
#SBATCH --output /home/hk-project-p0021769/hgf_ysb1444/logging/log_%j.out
#SBATCH --time 0-02:00:00
#SBATCH --partition=accelerated-h100

echo $SLURM_JOB_ID
PYTHONPATH=/hkfs/work/workspace/scratch/hgf_ysb1444-renderingHE/rendered_HE python render.py "$@"

retVal=$?
if [ $retVal -ne 0 ]; then
    echo "Error"
    exit 100
fi
exit 0
