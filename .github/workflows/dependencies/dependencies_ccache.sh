#!/usr/bin/env bash

if [[ $# -eq 2 ]]; then
  CVER=$1
else
  CVER=4.12.2 # If we update this, we might also need to update Tools/C_scripts/mmclt.py.
fi

wget https://github.com/ccache/ccache/releases/download/v${CVER}/ccache-${CVER}-linux-x86_64.tar.xz
tar xvf ccache-${CVER}-linux-x86_64.tar.xz
sudo cp -f ccache-${CVER}-linux-x86_64/ccache /usr/local/bin/
