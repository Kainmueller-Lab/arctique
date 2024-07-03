#!/bin/bash
### chdir in the correct folder and activate the correct python environment. Then run on login node with: qsub ./train.sh <args>
### some explanation:
### -V: currently set env variables (e.g. conda) are transferred
### -cwd: job is started in current directory
### -l h=maxg03,h=maxg04,h=maxg05,h=maxg06,h=maxg07,h=maxg08,gpu=1: nodes with the V100 GPUs installed (not T4) - will change soon!
### cuda_name = "A40-PCIE-45G"
### $ -l gpu=1 -l cuda_name="Tesla-V100-SXM2-32GB"
### $ -l gpu=1 -l cuda_name="Tesla-V100-PCIE-32GB"
### $ -l gpu=1 -l cuda_name="Tesla-V100-SXM2-16GB"
### $ -l h=maxg03,h=maxg04,h=maxg05,h=maxg06,h=maxg07,h=maxg08,gpu=1
###$ -l cuda_name = "A40-PCIE-45G"
### $ -pe mpi 1

#$ -l gpu=1
#$ -l m_mem_free=40G
#$ -l gpu_mem=16000
#$ -cwd
#$ -V
#$ -e error_log_$JOB_ID
#$ -o out_log_$JOB_ID
#$ -l h_rt=0:29:00
#$ -A kainmueller
#$ -pe mpi 2


CUDA_VISIBLE_DEVICES=0

python "$@"
retVal=$?
if [ $retVal -ne 0 ]; then
    echo "Error"
    exit 100
fi
