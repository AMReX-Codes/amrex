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
    libunwind-$1        \
    libunwind-$1-dev    \
    clang-$1            \
    libc++-$1-dev       \
    libc++abi-$1-dev    \
    libc++1-$1          \
    libc++abi1-$1
