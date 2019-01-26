#!/bin/bash -l

#SBATCH -p debug
#SBATCH -N 2
#SBATCH -t 00:29:00
#SBATCH -J Nyx_LyA_10Mpc_256_Cori_test_burst_buffer
#SBATCH -o Nyx_LyA_10Mpc_256_Cori_test_burst_buffer.%j.out
#SBATCH -e Nyx_LyA_10Mpc_256_Cori_test_burst_buffer.%j.err
#SBATCH -A cosmosim
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_90,TIME_LIMIT
#SBATCH --mail-user=friesen
#DW jobdw capacity=200GB access_mode=striped type=scratch
#DW stage_in source=/global/cscratch1/sd/friesen/bb_stage_in_256 destination=$DW_JOB_STRIPED/work type=directory
#DW stage_out source=$DW_JOB_STRIPED/work destination=/global/cscratch1/sd/friesen/bb_stage_out_256 type=directory

set -x

# We will stage in all necessary data to run the job (including the Nyx binary
# itself) in this folder.

SCRATCH_STAGE_DIR=$SCRATCH/bb_stage_in_256

# All of these files must already be in $SCRATCH_STAGE_DIR *BEFORE* you submit
# this job. DataWarp does NOT honor symlinks currently, so you must copy these
# data by hand.

EXE=Nyx3d.Linux.Intel.Intel.MPI.OMP.ex
INPUTS=inputs_cori_256_bb
PROBIN=probin
MISC_CP=TREECOOL_zl_dec14
MISC_LN=256.nyx

# This is where Nyx will dump its plot files/checkpoint files. It will be
# staged out to the location specified by the "stage_out" stanza in the SLURM
# preable above. If there is already a file/directory in the place where you
# tell DataWarp to stage this out, it will be overwritten, so be careful!

# The $DW_JOB_STRIPED variable is set by DataWarp itself, and is a temporary
# file system location on the burst buffer. It is created after job submission
# and torn down after job completion. You won't be able to access it before or
# after the job.

BB_WORK_DIR=$DW_JOB_STRIPED/work

[[ ! -d $BB_WORK_DIR ]] && mkdir -p $BB_WORK_DIR

cd $BB_WORK_DIR

# For some reason, DataWarp changes the permissions on some files during
# stage-in (usually removing executable permissions), so we have to change them
# back. Hopefully Cray will fix this soon.

chmod +x $EXE

# Run like normal. DataWarp will stage the data out after the job completes, so
# there may be a delay after job completion before you see the data in your
# stage_out location.

export OMP_NUM_THREADS=16
srun -n 4 -c 16 ./$EXE $INPUTS
