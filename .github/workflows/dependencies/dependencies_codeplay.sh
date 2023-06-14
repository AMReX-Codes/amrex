#!/usr/bin/env bash
#
# Copyright 2020 The AMReX Community
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

curl -o oneapi_nvidia.sh -L "https://developer.codeplay.com/api/v1/products/download?product=oneapi&variant=nvidia&filters[]=linux&aat=$1"
chmod +x oneapi_nvidia.sh
sudo ./oneapi_nvidia.sh --yes
