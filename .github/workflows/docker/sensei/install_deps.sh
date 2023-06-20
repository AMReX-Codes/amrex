#!/bin/bash

# DNF install dependencies
dnf update -y

# Development tools
dnf install -y --setopt=install_weak_deps=False \
    libasan libtsan libubsan clang clang-tools-extra \
    ninja-build cmake git which findutils patch

# MPI dependencies
dnf install -y --setopt=install_weak_deps=False \
    openmpi-devel

# External dependencies
dnf install -y --setopt=install_weak_deps=False \
    pugixml-devel

# Python dependencies
dnf install -y --setopt=install_weak_deps=False \
    python3 python3-devel python3-numpy \
    python3-mpi4py-openmpi swig

# Clean dnf cache
dnf clean all
