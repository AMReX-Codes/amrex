#!/usr/bin/env bash
# Build the cxx host compiler (with its MLIR/LLVM back end) from $CXX_SRC.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
# conda-forge's mlir package leaves PACKAGE_VERSION empty in MLIRConfigVersion.cmake,
# which fails cxx's "MLIR >= 21" check; record the installed version.
v=$TOOLS_ENV/lib/cmake/mlir/MLIRConfigVersion.cmake
if grep -q 'set(PACKAGE_VERSION "")' "$v"; then
  ver=$(sed -n 's/.*LLVM_PACKAGE_VERSION \([0-9.]*\).*/\1/p' "$TOOLS_ENV/lib/cmake/llvm/LLVMConfig.cmake" | head -1)
  sed -i "s/set(PACKAGE_VERSION \"\")/set(PACKAGE_VERSION \"${ver:-22.1.8}\")/" "$v"
fi
cmake -S "$CXX_SRC" -B "$CXX_BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCXX_ENABLE_MLIR=ON -DCXX_BUILD_TESTS=OFF \
  -DCMAKE_PREFIX_PATH="$TOOLS_ENV" -DMLIR_DIR="$TOOLS_ENV/lib/cmake/mlir" -DLLVM_DIR="$TOOLS_ENV/lib/cmake/llvm" \
  -DCMAKE_C_COMPILER="$CLANG" -DCMAKE_CXX_COMPILER="$CLANGXX" \
  -DCMAKE_C_FLAGS="--gcc-install-dir=$GCC_INSTALL_DIR" -DCMAKE_CXX_FLAGS="--gcc-install-dir=$GCC_INSTALL_DIR" \
  "-DCMAKE_EXE_LINKER_FLAGS=-Wl,-rpath,$TOOLS_ENV/lib"
cmake --build "$CXX_BUILD" --target cxx -j "$BUILD_JOBS"
"$CXX_BUILD/src/frontend/cxx" --help | head -1
