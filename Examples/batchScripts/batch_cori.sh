#!/bin/bash -l

#SBATCH -N 2
#SBATCH -t 01:00:00
#SBATCH -q regular
#SBATCH -C knl
#SBATCH -J <job name>
#SBATCH -A <allocation ID>
#SBATCH -e error.txt
#SBATCH -o output.txt

export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Typically use OMP_NUM_THREADS=8 and 8 MPI ranks per node,
# without hyperthreading (HT=1)
srun --cpu_bind=cores -n 16 -c 32 <path/to/executable> <input file>
# 4 OpenMP threads, 16 MPI ranks per node, HT=1
# srun --cpu_bind=cores -n 32 -c 16 <path/to/executable> <inputs file>
# 8 OpenMP threads, 16 MPI ranks per node, HT=2
# srun --cpu_bind=cores -n 32 -c 16 <path/to/executable> <inputs file>
