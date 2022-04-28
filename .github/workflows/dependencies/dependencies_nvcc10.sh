#!/usr/bin/env bash
#
# Copyright 2020 Axel Huebl
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential     \
    g++-6               \
    gfortran-6          \
    libopenmpi-dev      \
    openmpi-bin

sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list
sudo apt-get update
sudo apt-get install -y \
    cuda-command-line-tools-10-2 \
    cuda-compiler-10-2           \
    cuda-cupti-dev-10-2          \
    cuda-minimal-build-10-2      \
    cuda-nvml-dev-10-2           \
    cuda-nvtx-10-2               \
    cuda-curand-dev-10-2
sudo ln -s cuda-10.2 /usr/local/cuda
