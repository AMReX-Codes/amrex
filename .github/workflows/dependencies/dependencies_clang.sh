#!/usr/bin/env bash
#
# Copyright 2020-2022 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Axel Huebl

set -eu -o pipefail

clang_version=$1
install_libcxx=${2:-false}

packages=(
    build-essential
    libfftw3-dev
    gfortran
    clang-"${clang_version}"
)

case "${install_libcxx}" in
    1|true|TRUE|yes|YES|on|ON)
        packages+=(
            libc++-"${clang_version}"-dev
            libc++abi-"${clang_version}"-dev
            libc++1-"${clang_version}"
            libc++abi1-"${clang_version}"
        )
        ;;
esac

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

sudo apt-get update

sudo apt-get install -y --no-install-recommends "${packages[@]}"
