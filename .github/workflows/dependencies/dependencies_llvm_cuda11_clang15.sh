#!/usr/bin/env bash
#
# Copyright 2020-2022 Axel Huebl
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

sudo apt-get -qqq update
sudo apt-get install -y \
    libunwind-15        \
    libunwind-15-dev    \
    clang-15            \
    libc++-15-dev       \
    libc++abi-15-dev    \
    libc++1-15          \
    libc++abi1-15

.github/workflows/dependencies/dependencies_nvcc.sh 11.7
