#!/usr/bin/env bash

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

if [[ ! -f /etc/apt/trusted.gpg.d/apt.llvm.org.asc ]]; then
    wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | sudo tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc
fi

source /etc/os-release # set UBUNTU_CODENAME

sudo add-apt-repository "deb http://apt.llvm.org/${UBUNTU_CODENAME}/ llvm-toolchain-${UBUNTU_CODENAME} main"
sudo add-apt-repository "deb http://apt.llvm.org/${UBUNTU_CODENAME}/ llvm-toolchain-${UBUNTU_CODENAME}-$1 main"

sudo apt-get update

sudo apt-get install -y --no-install-recommends \
    clang-tidy-$1 libomp-$1-dev
