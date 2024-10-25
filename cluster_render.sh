#!/bin/bash
### chdir in the correct folder and activate the correct python environment. Then run on login node with: qsub ./train.sh <args>

#$ -l gpu=1
#$ -l m_mem_free=40G
#$ -l cuda_name="A40-PCIE-45G"
#$ -cwd
#$ -V
#$ -e error_log_$JOB_ID
#$ -o out_log_$JOB_ID
#$ -l h_rt=7:00:00
#$ -A kainmueller
#$ -pe mpi 2

CUDA_VISIBLE_DEVICES=0

python "$@"
retVal=$?
if [ $retVal -ne 0 ]; then
    echo "Error"
    exit 100
fi
