#!/bin/bash
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

jsrun -n 2 -a 6 -g 6 -c 6 --bind=packed:${omp} <path/to/executable> <input file> > output.txt
