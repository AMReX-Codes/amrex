#!/bin/bash

#SBATCH -J Random
#SBATCH -A m1759
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH --gres=gpu:1
#SBATCH -c 10
#SBATCH --mem=30GB
#SBATCH -p debug
#SBATCH -t 00:05:00

cuda-memcheck ./main3d.pgi.CUDA.ex
