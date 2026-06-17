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

echo 'Removing large SDKs and cached toolchains we do not need'

sudo rm -rf \
    /home/linuxbrew \
    /opt/az \
    /opt/ghc \
    /opt/hostedtoolcache \
    /opt/microsoft \
    /usr/local/.ghcup \
    /usr/local/aws-cli \
    /usr/local/aws-sam-cli \
    /usr/local/julia* \
    /usr/local/lib/android \
    /usr/local/share/boost \
    /usr/local/share/chromedriver-linux64 \
    /usr/local/share/edge_driver \
    /usr/local/share/gecko_driver \
    /usr/local/share/powershell \
    /usr/local/share/vcpkg \
    /usr/share/dotnet \
    /usr/share/miniconda \
    /usr/share/swift

if [[ -n "${AGENT_TOOLSDIRECTORY:-}" ]]; then
    sudo rm -rf "${AGENT_TOOLSDIRECTORY}"
fi

echo 'Removing Docker images and swap storage'

sudo docker image prune --all --force || true
sudo swapoff -a || true
sudo rm -f /mnt/swapfile

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
sudo apt-get remove -y '^postgresql-[0-9]+$'
sudo apt-get remove -y '^powershell.*'
sudo apt-get remove -y '^python3-botocore$'
sudo apt-get remove -y '^libruby[0-9].*'
sudo apt-get remove -y '^ruby[0-9].*-doc$'
sudo apt-get remove -y '^shellcheck$'
sudo apt-get remove -y '^snapd.*'
sudo apt-get remove -y '^skopeo$'
sudo apt-get remove -y '^temurin.*'

sudo apt-get autoremove -y
sudo apt-get clean
sudo rm -rf /var/lib/apt/lists/* /var/lib/docker /var/lib/containerd /var/lib/containers

df -h
