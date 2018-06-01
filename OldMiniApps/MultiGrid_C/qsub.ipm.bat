#PBS -q debug
#PBS -l walltime=00:30:00
#PBS -A accountnumber
#PBS -j eo

### total cores
#PBS -l mppwidth=64

###PBS -l mppwidth=1728
###PBS -l mppwidth=13824
###PBS -l mppwidth=110592

cd $PBS_O_WORKDIR

module load ipm
setenv IPM_LOGDIR $PBS_O_WORKDIR
setenv IPM_REPORT full
setenv IPM_LOG full
setenv IPM_HPM PAPI_L1_DCM,PAPI_TOT_INS,PAPI_L1_DCA,PAPI_FP_OPS

aprun -n 64 ./main3d.Linux.gcc.gfortran.MPI.IPM.ex
#aprun -n 1728 ./main3d.Linux.gcc.gfortran.MPI.IPM.ex
#aprun -n 13824 ./main3d.Linux.gcc.gfortran.MPI.IPM.ex
#aprun -n 110592 ./main3d.Linux.gcc.gfortran.MPI.IPM.ex

