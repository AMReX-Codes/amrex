#!/bin/bash

#PBS -q regular
#PBS -l nodes=1
#PBS -l walltime=00:29:59
#PBS -N Tiling_Heat_C_16Cx1T_1O_balanced

set -x

cd $PBS_O_WORKDIR

# KMP_AFFINITY sets the overall style of thread placement on the MIC card.
# "balanced" is the one which gives the best performance for the heat equation
# solver. There are other options available as well which may work better for
# certain problems. See the NERSC Babbage page for details.

export KMP_AFFINITY=balanced,verbose

# The MIC cards on Babbage have 60 cores. NERSC recommends avoid the OS core.
# However, in the native mode, it is safe to use all cores. 

export KMP_PLACE_THREADS=16Cx1T # 16 cores, 1 threads/core

WORKDIR=$SCRATCH/$PBS_JOBNAME.$PBS_JOBID
mkdir -p $WORKDIR

EXE=Heat3d.Linux.Intel.Intel.OMP.ex
INPUT=inputs_3d

cp $EXE $WORKDIR
cp $INPUT $WORKDIR

cd $WORKDIR

get_micfile

mpirun.mic -n 1 -hostfile micfile.$PBS_JOBID -ppn 1 ./$EXE $INPUT
