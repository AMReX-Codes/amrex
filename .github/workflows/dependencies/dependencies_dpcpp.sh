#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

# Ref.: https://github.com/rscohn2/oneapi-ci
# intel-basekit intel-hpckit are too large in size
wget -q -O - https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB \
  | sudo apt-key add -
echo "deb https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt-get update

# try apt install up to five times, to avoid connection splits
status=1
for itry in {1..5}
do
    sudo apt-get install -y --no-install-recommends \
        build-essential \
        intel-oneapi-dpcpp-cpp-compiler \
        intel-oneapi-compiler-fortran \
        intel-oneapi-mkl-devel \
        intel-oneapi-mpi-devel \
        && { sudo apt-get clean; status=0; break; }  \
        || { sleep 10; }
done
if [[ ${status} -ne 0 ]]; then exit 1; fi
