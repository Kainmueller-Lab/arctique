#!/bin/bash
START_BASE=0
N_SAMPLES=50
BASE_16BIT=257

# add adjustable start shift as entry from shell
if [ -z "$1" ]; then
    START_BASE=0
else
    START_BASE=$1
    fi

for i in $(seq 0 99)
do
    echo "Starting job $i"
    echo "qsub render.sh render.py --start-idx $((START_BASE + i * N_SAMPLES)) --n-samples $N_SAMPLES --base-16bit $BASE_16BIT"
    qsub render.sh render.py --start-idx $((START_BASE + i * N_SAMPLES)) --n-samples $N_SAMPLES --base-16bit $BASE_16BIT
done