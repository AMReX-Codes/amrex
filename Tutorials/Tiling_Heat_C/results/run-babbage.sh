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

# KMP_PLACE_THREADS provides more fine-grained control of thread placement. A
# few things to keep in mind:
#
# 1.) KMP_PLACE_THREADS overrides OMP_NUM_THREADS. So do NOT use both of those;
# pick one or the other. The offset core should avoid the 1 core reserved for
# the OS.
#
# 2.) The MIC cards on Babbage have 60 cores, not 61. So the most you should
# ever use is 59 cores.
#
# 3.) The "1O" tells the MIC card to "offset" the threads by 1 core, i.e. start
# placing threads on core 1 instead of core 0. The default is "0O", i.e., start
# placing threads on core 0. The NERSC documentation says that the OS core is
# logical core 0 but physical core 59. I don't know if this offset parameter
# describes a logical or a physical offset, though.

export KMP_PLACE_THREADS=16Cx1T,1O # 16 cores, 1 threads/core, 1 core offset

WORKDIR=$SCRATCH/$PBS_JOBNAME.$PBS_JOBID
mkdir -p $WORKDIR

EXE=Heat3d.Linux.Intel.Intel.OMP.ex
INPUT=inputs_3d

cp $EXE $WORKDIR
cp $INPUT $WORKDIR

cd $WORKDIR

get_micfile

mpirun.mic -n 1 -hostfile micfile.$PBS_JOBID -ppn 1 ./$EXE $INPUT
