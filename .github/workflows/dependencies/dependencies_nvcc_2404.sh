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

VERSION_DOTTED=${1-12.0} && VERSION_DASHED=$(sed 's/\./-/' <<< $VERSION_DOTTED)
curl -O https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get install -y \
    cuda-command-line-tools-$VERSION_DASHED \
    cuda-compiler-$VERSION_DASHED           \
    cuda-cupti-dev-$VERSION_DASHED          \
    cuda-minimal-build-$VERSION_DASHED      \
    cuda-nvml-dev-$VERSION_DASHED           \
    cuda-nvtx-$VERSION_DASHED               \
    libcufft-dev-$VERSION_DASHED            \
    libcurand-dev-$VERSION_DASHED           \
    libcusparse-dev-$VERSION_DASHED

sudo apt-get install -y --no-install-recommends libnvjitlink-dev-$VERSION_DASHED || true

sudo ln -s cuda-$VERSION_DOTTED /usr/local/cuda
