#PBS -N test_sidecar
#PBS -l mppwidth=24
#PBS -q debug
#PBS -l walltime=00:29:59
#PBS -j oe

cd $PBS_O_WORKDIR

WORKDIR=$SCRATCH/$PBS_JOBNAME.$PBS_JOBID

mkdir -p $WORKDIR

EXE=gridmovetest3d.Linux.g++.gfortran.MPI.ex
INPUTS=inputs_3d

cp $EXE $INPUTS $WORKDIR

cd $WORKDIR

aprun -n 24 ./$EXE $INPUTS
