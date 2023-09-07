#!/usr/bin/env bash
#
# Copyright 2022 The AMReX Community
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL
#

set -eu -o pipefail

cd Docs/Doxygen

# treat all warnings as errors (TODO)
#echo "WARN_AS_ERROR = YES" >> doxygen.conf

doxygen doxygen.conf
