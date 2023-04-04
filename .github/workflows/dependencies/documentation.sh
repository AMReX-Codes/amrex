#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL
# Authors: Andrew Myers

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential \
    pandoc \
    doxygen \
    texlive \
    texlive-latex-extra \
    texlive-lang-cjk \
    tex-gyre \
    latexmk

