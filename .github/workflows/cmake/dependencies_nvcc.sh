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
    build-essential     \
    cmake               \
    g++-5               \
    gfortran-5          \
    libopenmpi-dev      \
    openmpi-bin         \
    nvidia-cuda-dev     \
    nvidia-cuda-toolkit

# Patch broken GCC 5.5 libs in <algorithm>
#   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=76731
#   https://stackoverflow.com/a/50815334/2719194
for f in avx512fintrin.h avx512pfintrin.h avx512vlintrin.h
do
   curl \
     -H "User-Agent:Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/62.0.3202.94 Safari/537.36" \
     -o ${f} "https://gcc.gnu.org/viewcvs/gcc/branches/gcc-5-branch/gcc/config/i386/${f}?view=co&revision=245536&content-type=text%2Fplain&pathrev=245536"
done
sudo mv avx512*intrin.h /usr/lib/gcc/x86_64-linux-gnu/5/include/

