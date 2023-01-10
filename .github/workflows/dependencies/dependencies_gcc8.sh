#!/usr/bin/env bash
#
# Copyright 2020-2022 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update

sudo apt-get install -y --no-install-recommends \
    build-essential    \
    g++-8 gfortran-8   \
    libopenmpi-dev     \
    openmpi-bin
