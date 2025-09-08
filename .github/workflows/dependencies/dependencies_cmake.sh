#!/usr/bin/env bash

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

test -f /usr/share/doc/kitware-archive-keyring/copyright ||
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null

if [[ ! -f /etc/apt/trusted.gpg.d/apt.llvm.org.asc ]]; then
    wget -qO- https://apt.llvm.org/llvm-snapshot.gpg.key | sudo tee /etc/apt/trusted.gpg.d/apt.llvm.org.asc
fi

source /etc/os-release # set UBUNTU_CODENAME

echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ ${UBUNTU_CODENAME} main" | sudo tee /etc/apt/sources.list.d/kitware.list >/dev/null

sudo apt-get update

test -f /usr/share/doc/kitware-archive-keyring/copyright ||
sudo rm /usr/share/keyrings/kitware-archive-keyring.gpg

sudo apt-get install kitware-archive-keyring

sudo apt-get install -y --no-install-recommends cmake

sudo rm -f /usr/local/bin/cmake
cmake --version
