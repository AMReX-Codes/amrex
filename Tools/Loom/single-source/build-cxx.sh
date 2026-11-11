#!/usr/bin/env bash
# Prepare pixi's MLIR package for cxx's MLIR/LLVM back end. cxx itself is built
# by build-hrx.sh, inside HRX's build (LOOM_CXX_HOST_CODEGEN), together with
# Loom's single-source driver loom-cxx.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
# conda-forge's mlir package leaves PACKAGE_VERSION empty in MLIRConfigVersion.cmake,
# which fails cxx's "MLIR >= 21" check; record the installed version.
v=$TOOLS_ENV/lib/cmake/mlir/MLIRConfigVersion.cmake
if grep -q 'set(PACKAGE_VERSION "")' "$v"; then
  ver=$(sed -n 's/.*LLVM_PACKAGE_VERSION \([0-9.]*\).*/\1/p' "$TOOLS_ENV/lib/cmake/llvm/LLVMConfig.cmake" | head -1)
  sed -i "s/set(PACKAGE_VERSION \"\")/set(PACKAGE_VERSION \"${ver:-22.1.8}\")/" "$v"
fi
grep 'set(PACKAGE_VERSION' "$v"
