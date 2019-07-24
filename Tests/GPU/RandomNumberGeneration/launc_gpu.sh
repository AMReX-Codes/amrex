#!/bin/bash

#SBATCH -J Random2
#SBATCH -A m1759
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node 1
#SBATCH -c 1
#SBATCH --gres=gpu:1



#SBATCH -C gpu
#SBATCH -p debug
#SBATCH -t 00:10:00

nvprof -o profile.nvvp --analysis-metrics ./main3d.pgi.MPI.CUDA.ex
