#!/bin/bash

#BSUB -P csc190port
#BSUB -W 01:00
#BSUB -nnodes 1
#BSUB -J HeatEquation_EX1_C_OMP45
#BSUB -o HeatEquation_EX1_C_OMP45.%J.out
#BSUB -e HeatEquation_EX1_C_OMP45.%J.out

exe=main3d.ibm.ex
inputs=inputs_3d
repo=csc190

# This is not yet defined on Summit, so set it here.
if [ ! -v MEMBERWORK ]; then
  MEMBERWORK=/gpfs/alpinetds/scratch/${USER}/${repo}/
fi

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/summit/cuda/9.1.85/lib64

cp ${exe} ${inputs} ${MEMBERWORK}

cd ${MEMBERWORK}

jsrun -n 1 -r 1 -a 1 -g 1 ./${exe} ${inputs}
