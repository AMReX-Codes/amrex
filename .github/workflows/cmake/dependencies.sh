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

sudo apt-get install -y --no-install-recommends\
    build-essential \
    cmake           \
    g++ gfortran    \
    libopenmpi-dev  \
    openmpi-bin
