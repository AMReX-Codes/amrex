#PBS -q debug
#PBS -l walltime=00:30:00
#PBS -A accountnumber
#PBS -j eo

### total cores
#PBS -l mppwidth=125

###PBS -l mppwidth=1000
###PBS -l mppwidth=10648
###PBS -l mppwidth=110592

cd $PBS_O_WORKDIR

module load ipm
setenv IPM_LOGDIR $PBS_O_WORKDIR
setenv IPM_REPORT full
setenv IPM_LOG full
setenv IPM_HPM PAPI_L1_DCM,PAPI_TOT_INS,PAPI_L1_DCA,PAPI_FP_OPS

aprun -n 125 ./fbtest3d.Linux.g++.gfortran.MPI.IPM.ex
#aprun -n 1000 ./fbtest3d.Linux.g++.gfortran.MPI.IPM.ex
#aprun -n 10648 ./fbtest3d.Linux.g++.gfortran.MPI.IPM.ex
#aprun -n 110592 ./fbtest3d.Linux.g++.gfortran.MPI.IPM.ex

