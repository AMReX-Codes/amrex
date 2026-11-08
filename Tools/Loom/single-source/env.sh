# Environment for single-source LOOM builds of AMReX. Source it before each
# step (build-cxx.sh, build-hrx.sh, build-amrex.sh, run-tests.sh do this).
#
# Edit the first block for your machine; everything else derives from it.

# --- Machine-specific -------------------------------------------------------
# ROCm 10.1 (provides clang/clang++ for building and linking, and libhsa).
: "${ROCM_PATH:=/opt/rocm/core-10.1}"
# GCC 14 install directory (C++23 libstdc++ for building the tools; compiler
# headers such as stddef.h for the single-source compiles).
: "${GCC_INSTALL_DIR:=/usr/lib64/gcc/x86_64-suse-linux/14}"
# glibc headers seen by the device pass (and by the host via the default sysroot).
: "${GLIBC_INCLUDE:=/usr/include}"
# Where the three clones live, and where builds go.
: "${LOOM_WORK:=$PWD}"
: "${CXX_SRC:=$LOOM_WORK/cplusplus}"      # branch loom-single-source
: "${HRX_SRC:=$LOOM_WORK/hrx}"            # branch loom-single-source
: "${AMREX_SRC:=$LOOM_WORK/amrex}"        # branch loom-single-source
: "${BUILD_JOBS:=8}"

# --- Derived ----------------------------------------------------------------
export ROCM_PATH GCC_INSTALL_DIR GLIBC_INCLUDE LOOM_WORK CXX_SRC HRX_SRC AMREX_SRC BUILD_JOBS
export TOOLS_ENV=$AMREX_SRC/Tools/Loom/single-source/.pixi/envs/default
export CXX_BUILD=$LOOM_WORK/build/cxx
export HRX_BUILD=$LOOM_WORK/build/hrx
export AMREX_BUILD=$LOOM_WORK/build/amrex
export PATH=$TOOLS_ENV/bin:$PATH
export TMPDIR=${TMPDIR:-$LOOM_WORK/tmp}
mkdir -p "$TMPDIR"
export ROCM_CLANG=$ROCM_PATH/llvm/bin/clang
export ROCM_CLANGXX=$ROCM_PATH/llvm/bin/clang++
[[ -x $ROCM_CLANGXX ]] || { export ROCM_CLANG=$ROCM_PATH/lib/llvm/bin/clang ROCM_CLANGXX=$ROCM_PATH/lib/llvm/bin/clang++; }
export IREE_HAL_AMDGPU_LIBHSA_PATH=$ROCM_PATH/lib/libhsa-runtime64.so.1

# amrex-loom-cxx (the C++ compiler AMReX is configured with).
export AMREX_LOOM_HOST_CXX=$AMREX_SRC/Tools/Loom/amrex-loom-host-cxx
export AMREX_LOOM_IMPORT=$HRX_BUILD/loom/src/loom/tools/loom-import-cxx/loom-import-cxx
export AMREX_LOOM_COMPILE=$HRX_BUILD/loom/src/loom/tools/loom-compile/loom-compile
export AMREX_LOOM_TARGET=gfx90a
export AMREX_LOOM_LOOMCXX_DIR=$LOOM_WORK/build/loomcxx-include
export AMREX_LOOM_SYSTEM_INCLUDES=$TOOLS_ENV/include/c++/v1:$GCC_INSTALL_DIR/include:$GLIBC_INCLUDE
export AMREX_LOOM_DEVICE_PRELUDE=$AMREX_SRC/Tools/Loom/amrex_loom_device_prelude_libcxx.H
export AMREX_LOOM_LAYOUT_CHECK=0
export AMREX_LOOM_HRX_LIB_DIR=$HRX_BUILD/libhrx/src/libhrx
export AMREX_LOOM_HRX_INCLUDE_DIR=$HRX_SRC/libhrx/include

# amrex-loom-host-cxx (cxx for compiles, ROCm clang++ for links, libc++ for both).
export AMREX_LOOM_CXX=$CXX_BUILD/src/frontend/cxx
export AMREX_LOOM_LIBCXX_ROOT=$TOOLS_ENV
export AMREX_LOOM_GCC_INSTALL_DIR=$GCC_INSTALL_DIR
export AMREX_LOOM_LINK_CXX=$ROCM_CLANGXX
