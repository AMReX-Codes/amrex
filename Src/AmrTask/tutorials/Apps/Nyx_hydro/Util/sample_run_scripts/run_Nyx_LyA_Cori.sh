#!/bin/bash

# New SLURM syntax for Cori:
#
# -p = partition/queue ("regular", "debug", etc.)
# -N = number of nodes (32 cores/node on Cori Haswell nodes)
# --ntasks-per-node = number of tasks per node (if you don't want to pack the node)
# -t = time
# -J = job name
# -o = STDOUT file
# -e = STDERR file (merges with STDOUT if -e is not given explicitly)
# -A = account to charge
# --mail-type = events to e-mail user about (ALL is short for BEGIN,END,FAIL,REQUEUE,STAGE_OUT)
# --mail-user = NIM username of user to notify
#
# SLURM by default will cd to your working directory when you submit the job,
# so while you can do "cd $SLURM_SUBMIT_DIR" if you want to (like "cd
# $PBS_O_WORKDIR" with Torque), it will do it automatically even if you don't.

#SBATCH -p debug
#SBATCH -N 2
#SBATCH --ntasks-per-node=32
#SBATCH -t 00:29:59
#SBATCH -J Nyx_LyA_512_Cori_test
#SBATCH -o Nyx_LyA_512_Cori_test.%j.out
#SBATCH -e Nyx_LyA_512_Cori_test.%j.err
#SBATCH -A m2010
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_90,TIME_LIMIT
#SBATCH --mail-user=friesen

WORKDIR=$SCRATCH/$SLURM_JOB_NAME.$SLURM_JOB_ID

[[ ! -d "$WORKDIR" ]] && mkdir -p "$WORKDIR"

EXE_DIR=$HOME/Nyx_runs/LyA/20Mpc/512
EXE=Nyx3d.Linux.Cray.Cray.MPI.Cori.ex
INPUTS=inputs_no_analysis_test_cori
PROBIN=probin
MISC_CP=/global/homes/f/friesen/Nyx/Source/HeatCool/TREECOOL_zl_dec14
MISC_LN=/project/projectdirs/m1796/Runs/20Mpc/512/512ss_20mpc.nyx

cp $INPUTS $PROBIN $MISC_CP $WORKDIR
ln -s $MISC_LN $WORKDIR
cp $EXE_DIR/$EXE $WORKDIR

cd $WORKDIR

# srun syntax (replacement for "aprun")
#
# -n = number of MPI tasks
# --ntasks-per-node = number of tasks per node (if you don't want to pack the node)
# --ntasks-per-socket = number of tasks per socket (16 cores/socket on Haswell)
# -c = number of CPUs per task (greater than 1 if you're using OpenMP)
#
# Some of these commands are redundant with SBATCH headers. I don't know which
# ones take precedence so I just put them in both places.

srun -n 64 --ntasks-per-node=32 --ntasks-per-socket=16 -c 1 ./$EXE $INPUTS
