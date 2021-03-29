#!/bin/bash
#BSUB -P CSC308
#BSUB -W 0:20
#BSUB -nnodes 2048
#BSUB -J amrex
#BSUB -o amrexo.%J
#BSUB -e amrexe.%J

module load gcc
module load cuda
module list
set -x

omp=1
export PAMI_DISABLE_IPC=1
export OMP_NUM_THREADS=${omp}

EXE="./main3d.gnu.TPROF.MPI.CUDA.ex"
SMPIARGS= --smpiargs="-x PAMI_DISABLE_CUDA_HOOK=1 -disable_gpu_hooks"

NUMNODES=12288
NUMCELLS=2048

jsrun -n ${NUMNODES} -a 1 -g 1 -c 1 --bind=packed:${omp} ${SMPIARGS} ${EXE} inputs.test n_cell=${NUMCELLS} verbose=0 max_fmg_iter=0 > output_2048_${LSB_JOBID}.txt
