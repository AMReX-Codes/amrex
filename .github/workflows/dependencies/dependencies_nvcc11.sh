#!/usr/bin/env bash
#
# Copyright 2020-2022 Axel Huebl
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

sudo apt-get -qqq update
sudo apt-get install -y \
    build-essential     \
    ca-certificates     \
    cmake               \
    g++                 \
    gfortran            \
    gnupg               \
    libopenmpi-dev      \
    openmpi-bin         \
    pkg-config          \
    wget

curl -O https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb
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

