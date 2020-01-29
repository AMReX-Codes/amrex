#!/bin/bash

# Copyright 2019 Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

#SBATCH --job-name=postproc
#SBATCH --time=00:20:00
#SBATCH -C haswell
#SBATCH -N 8
#SBATCH -q regular
#SBATCH -e postproce.txt
#SBATCH -o postproco.txt
#SBATCH --mail-type=end
#SBATCH --account=m2852

export OMP_NUM_THREADS=1

# Requires python3 and yt > 3.5
srun -n 32 -c 16 python plot_parallel.py --path <path/to/plotfiles> --plotlib=yt --parallel
