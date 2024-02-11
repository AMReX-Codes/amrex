#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

# `man apt.conf`:
#   Number of retries to perform. If this is non-zero APT will retry
#   failed files the given number of times.
echo 'Acquire::Retries "3";' | sudo tee /etc/apt/apt.conf.d/80-retries

# https://developer.codeplay.com/apt/index.html
sudo wget -qO - https://developer.codeplay.com/apt/public.key | gpg --dearmor | sudo tee /usr/share/keyrings/codeplay-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/codeplay-keyring.gpg] https://developer.codeplay.com/apt all main" | sudo tee /etc/apt/sources.list.d/codeplay.list

sudo apt-get clean
sudo apt-get update

# try apt install up to five times, to avoid connection splits
status=1
for itry in {1..5}
do
    sudo apt-get install -y --no-install-recommends \
        $1 \
        && { sudo apt-get clean; sudo apt-get update; status=0; break; }  \
        || { sleep 10; }
done
if [[ ${status} -ne 0 ]]; then exit 1; fi
