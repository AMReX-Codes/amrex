#!/bin/bash

# Copyright 2019 Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

#BSUB -P <allocation ID>
#BSUB -W 00:10
#BSUB -nnodes 2
#BSUB -J WarpX
#BSUB -o WarpXo.%J
#BSUB -e WarpXe.%J

module load pgi
module load cuda

omp=1
export OMP_NUM_THREADS=${omp}

num_nodes=$(( $(printf '%s\n' ${LSB_HOSTS} | sort -u | wc -l) - 1 ))
jsrun -n ${num_nodes} -a 6 -g 6 -c 6 --bind=packed:${omp} <path/to/executable> <input file> > output.txt
