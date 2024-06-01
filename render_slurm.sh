#!/bin/sh
#SBATCH --gres=gpu:full:1
#SBATCH -n 8
#SBATCH --cpus-per-task=4
#SBATCH -e /home/hgf_mdc/hgf_lrx5378/logging/log_%j.err
#SBATCH --output /home/hgf_mdc/hgf_lrx5378/logging/log_%j.out
#SBATCH --time 0-08:00:00


N_GPUS=1
echo $SLURM_JOB_ID
PYTHONPATH=/hkfs/work/workspace_haic/scratch/hgf_lrx5378-darmstadt/rendered_HE python render.py "$@"

retVal=$?
if [ $retVal -ne 0 ]; then
    echo "Error"
    exit 100
fi
exit 0