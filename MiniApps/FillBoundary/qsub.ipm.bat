#PBS -q debug
#PBS -l walltime=00:30:00
#PBS -A accountnumber
#PBS -j eo

### total cores
#PBS -l mppwidth=64

cd $PBS_O_WORKDIR

module load ipm
setenv IPM_LOGDIR $PBS_O_WORKDIR
setenv IPM_REPORT full
setenv IPM_LOG full
setenv IPM_HPM PAPI_L1_DCM,PAPI_TOT_INS,PAPI_L1_DCA,PAPI_FP_OPS

aprun -n 64 ./fbtest3d.Linux.g++.gfortran.MPI.IPM.ex

