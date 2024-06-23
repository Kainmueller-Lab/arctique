#!/bin/bash
start_base = 0
n_samples = 200

# add adjustable start shift as entry from shell
if [ $# -eq 1 ]
then
    start_base = $1
fi

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
do
    # start a python script
    echo "Starting job $i"
    echo "qsub render.sh render.py --start-idx $((start_base + i * n_samples)) --n-samples $n_samples"
    #qsub render.sh render.py --start-idx $((start_base + i * n_samples)) --n-samples $n_samples
done