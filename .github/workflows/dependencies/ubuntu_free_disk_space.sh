#!/usr/bin/env bash
#
# Copyright 2023 The AMReX Community
#
# License: BSD-3-Clause-LBNL

# Don't want to use the following line because apt-get remove may fail if
# the package specfied does not exist.
# set -eu -o pipefail

# Large packages
dpkg-query -Wf '${Installed-Size}\t${Package}\n' | sort -n | tail -n 100

echo 'Removing some packages we do not need'

df -h

apt list --installed

sudo apt-get remove -y '^apache.*'
sudo apt-get remove -y '^aspnetcore.*'
sudo apt-get remove -y '^azure.*'
sudo apt-get remove -y '^dotnet.*'
sudo apt-get remove -y '^firebird.*'
sudo apt-get remove -y '^firefox.*'
sudo apt-get remove -y '^google.*'
sudo apt-get remove -y '^hhvm.*'
sudo apt-get remove -y '^microsoft.*'
sudo apt-get remove -y '^mongodb.*'
sudo apt-get remove -y '^mono-.*'
sudo apt-get remove -y '^monodoc-.*'
sudo apt-get remove -y '^mysql.*'
sudo apt-get remove -y '^php.*'
sudo apt-get remove -y '^powershell.*'
sudo apt-get remove -y '^snapd.*'
sudo apt-get remove -y '^temurin.*'

sudo apt-get autoremove -y

df -h
