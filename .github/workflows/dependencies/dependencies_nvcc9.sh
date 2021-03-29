#!/usr/bin/env bash
#
# Copyright 2020 Axel Huebl
#
# License: BSD-3-Clause-LBNL

# search recursive inside a folder if a file contains tabs
#
# @result 0 if no files are found, else 1
#

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential     \
    g++-6               \
    gfortran-6          \
    libopenmpi-dev      \
    openmpi-bin

sudo wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/7fa2af80.pub
sudo apt-key add 7fa2af80.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64 /" \
    | sudo tee /etc/apt/sources.list.d/cuda.list
sudo apt-get update
sudo apt-get install -y \
    cuda-command-line-tools-9-2 \
    cuda-compiler-9-2           \
    cuda-minimal-build-9-2      \
    cuda-nvml-dev-9-2           \
    cuda-nvtx-9-2               \
    cuda-curand-dev-9.2

sudo ln -s cuda-9.2 /usr/local/cuda


