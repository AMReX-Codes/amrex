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
sudo apt-get remove -y '^buildah$'
sudo apt-get remove -y '^containerd.*'
sudo apt-get remove -y '^docker-ce$'
sudo apt-get remove -y '^docker-ce-cli$'
sudo apt-get remove -y '^dotnet.*'
sudo apt-get remove -y '^firebird.*'
sudo apt-get remove -y '^firefox.*'
sudo apt-get remove -y '^google.*'
sudo apt-get remove -y '^hhvm.*'
sudo apt-get remove -y '^humanity-icon-theme$'
sudo apt-get remove -y '^kubectl$'
sudo apt-get remove -y '^libruby3\.2$'
sudo apt-get remove -y '^mecab-ipadic$'
sudo apt-get remove -y '^mercurial.*'
sudo apt-get remove -y '^mesa-libgallium$'
sudo apt-get remove -y '^microsoft.*'
sudo apt-get remove -y '^mongodb.*'
sudo apt-get remove -y '^mono-.*'
sudo apt-get remove -y '^monodoc-.*'
sudo apt-get remove -y '^mysql.*'
sudo apt-get remove -y '^php.*'
sudo apt-get remove -y '^podman$'
sudo apt-get remove -y '^postgresql-16$'
sudo apt-get remove -y '^powershell.*'
sudo apt-get remove -y '^python3-botocore$'
sudo apt-get remove -y '^ruby3\.2-doc$'
sudo apt-get remove -y '^shellcheck$'
sudo apt-get remove -y '^snapd.*'
sudo apt-get remove -y '^skopeo$'
sudo apt-get remove -y '^temurin.*'

sudo apt-get autoremove -y

df -h
