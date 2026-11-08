#!/usr/bin/env bash
# Build Loom's C++ importer, loom-compile and libhrx from $HRX_SRC.
# The importer's embedded cxx comes from $CXX_SRC, the same checkout as the host
# compiler, so host and device passes always use one cxx revision. Only HRX's
# native AMD GPU driver (libamdf) is built; the legacy ROCr driver is off.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cmake -S "$HRX_SRC" -B "$HRX_BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DLOOM_IMPORT_CXX=ON -DLIBHRX_BUILD=ON \
  -DIREE_HAL_DRIVER_AMD_GPU=ON -DIREE_HAL_DRIVER_AMDGPU=OFF -DIREE_HAL_DRIVER_AMD_XDNA=OFF \
  -DFETCHCONTENT_SOURCE_DIR_CXX="$CXX_SRC" \
  -DCMAKE_C_COMPILER="$CLANG" -DCMAKE_CXX_COMPILER="$CLANGXX" \
  -DCMAKE_C_FLAGS="--gcc-install-dir=$GCC_INSTALL_DIR" -DCMAKE_CXX_FLAGS="--gcc-install-dir=$GCC_INSTALL_DIR"
cmake --build "$HRX_BUILD" --target loom-import-cxx loom-compile hrx -j "$BUILD_JOBS"
mkdir -p "$AMREX_LOOM_LOOMCXX_DIR"
ln -sfn "$HRX_SRC/loom/src/loom/import/cxx/include/loomcxx" "$AMREX_LOOM_LOOMCXX_DIR/loomcxx"
ls -la "$AMREX_LOOM_IMPORT" "$AMREX_LOOM_COMPILE"
