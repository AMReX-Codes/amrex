#!/usr/bin/env bash
# Build Loom's single-source driver loom-cxx (with cxx and its MLIR back end),
# the C++ importer, loom-compile and libhrx from $HRX_SRC. cxx comes from
# $CXX_SRC through cxx's own CMake project (LOOM_CXX_HOST_CODEGEN), so the
# driver links one copy of cxx for host code generation and device import.
# Only HRX's native AMD GPU driver (libamdf) is built; the legacy ROCr driver is
# off. HRX's CMake accepts pixi's MLIR package as installed (its version file
# reports no version; HRX takes LLVM's).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cmake -S "$HRX_SRC" -B "$HRX_BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DLOOM_IMPORT_CXX=ON -DLOOM_CXX_HOST_CODEGEN=ON -DLIBHRX_BUILD=ON \
  -DIREE_HAL_DRIVER_AMD_GPU=ON -DIREE_HAL_DRIVER_AMDGPU=OFF -DIREE_HAL_DRIVER_AMD_XDNA=OFF \
  -DFETCHCONTENT_SOURCE_DIR_CXX="$CXX_SRC" \
  -DMLIR_DIR="$TOOLS_ENV/lib/cmake/mlir" -DLLVM_DIR="$TOOLS_ENV/lib/cmake/llvm" \
  -DCMAKE_C_COMPILER="$CLANG" -DCMAKE_CXX_COMPILER="$CLANGXX" \
  -DCMAKE_C_FLAGS="--gcc-install-dir=$GCC_INSTALL_DIR" -DCMAKE_CXX_FLAGS="--gcc-install-dir=$GCC_INSTALL_DIR" \
  "-DCMAKE_EXE_LINKER_FLAGS=-Wl,-rpath,$TOOLS_ENV/lib"
cmake --build "$HRX_BUILD" --target loom-cxx loom-import-cxx loom-compile hrx -j "$BUILD_JOBS"
"$AMREX_LOOM_CXX_DRIVER" --help | head -1
