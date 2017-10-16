#!/bin/bash

#SBATCH -N 16
#SBATCH -C knl,quad,cache
#SBATCH -p regular
#SBATCH -J WarpX
#SBATCH -t 00:30:00
#SBATCH --mail-user=atmyers@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --export=ALL

cd $SLURM_SUBMIT_DIR

date
srun -n 512 -c 8 --cpu_bind=cores ./main3d.intel.TPROF.MPI.ex inputs &> out
date
