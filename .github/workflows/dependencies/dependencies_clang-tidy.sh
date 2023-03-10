#!/usr/bin/env bash

set -eu -o pipefail

sudo apt-get install -y --no-install-recommends \
    clang-tidy-$1 libomp-$1-dev
