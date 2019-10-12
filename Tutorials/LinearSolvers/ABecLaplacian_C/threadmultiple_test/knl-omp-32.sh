#!/bin/bash
#SBATCH -N 32
#SBATCH -C knl
#SBATCH -q debug 
#SBATCH -t 00:10:00
#SBATCH -o 32-OMP.out

#OpenMP settings:
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

NCELLS=1024

export MPICH_MAX_THREAD_SAFETY=multiple

#run the application:
srun -n 256 -c 32 --cpu_bind=cores main3d.gnu.mic-knl.TPROF.MPI.OMP.ex inputs.test n_cell=${NCELLS}  
