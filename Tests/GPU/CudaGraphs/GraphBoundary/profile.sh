#!/bin/bash
nvprof -o nvprof_$$.nvvp ./main3d.gnu.TPROF.MPI.CUDA.ex inputs_3d
#nv-nsight-cu-cli -o nsight_$$.gui ./main3d.gnu.TPROF.MPI.CUDA.ex inputs_3d
