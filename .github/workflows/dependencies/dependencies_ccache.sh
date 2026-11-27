#!/usr/bin/env bash

if [[ $# -eq 2 ]]; then
  CVER=$1
else
  CVER=4.12.2 # If we update this, we might also need to update Tools/C_scripts/mmclt.py.
fi

# ccache publishes Linux release binaries for x86_64 and aarch64, named after
# `uname -m`, which covers both the x64 and the arm64 GitHub-hosted runners.
CARCH=$(uname -m)
case "${CARCH}" in
  x86_64|aarch64) ;;
  *) echo "No ccache ${CVER} release binary for ${CARCH}" >&2; exit 1 ;;
esac

wget https://github.com/ccache/ccache/releases/download/v${CVER}/ccache-${CVER}-linux-${CARCH}.tar.xz
tar xvf ccache-${CVER}-linux-${CARCH}.tar.xz
sudo cp -f ccache-${CVER}-linux-${CARCH}/ccache /usr/local/bin/
