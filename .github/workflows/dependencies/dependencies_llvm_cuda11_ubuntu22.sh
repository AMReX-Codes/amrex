#!/usr/bin/env bash
#
# Copyright 2020-2022 Axel Huebl
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    clang-15            \
    libc++-15-dev       \
    libc++abi-15-dev    \
    libc++1-15          \
    libc++abi1-15       \
    cmake               \
    gnupg               \
    libopenmpi-dev      \
    openmpi-bin         \
    pkg-config          \
    wget

curl -O https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.0-1_all.deb
sudo dpkg -i cuda-keyring_1.0-1_all.deb
sudo apt-get update
sudo apt-get install -y \
    cuda-command-line-tools-11-2 \
    cuda-compiler-11-2           \
    cuda-cupti-dev-11-2          \
    cuda-minimal-build-11-2      \
    cuda-nvml-dev-11-2           \
    cuda-nvtx-11-2               \
    libcurand-dev-11-2
sudo ln -s cuda-11.2 /usr/local/cuda
