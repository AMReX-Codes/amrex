# Environment for single-source LOOM builds of AMReX. Source it before each
# step (build-hrx.sh, build-amrex.sh, run-tests.sh do this).
#
# Edit the first block for your machine; everything else derives from it.

# --- Machine-specific -------------------------------------------------------
# GCC 14 install directory (C++23 libstdc++ for building the tools, glibc crt
# files for linking, compiler headers such as stddef.h for the single-source
# compiles).
: "${GCC_INSTALL_DIR:=/usr/lib64/gcc/x86_64-suse-linux/14}"
# glibc headers seen by the single-source compiles.
: "${GLIBC_INCLUDE:=/usr/include}"
# Sysroot for linking. It must hold the same glibc as GLIBC_INCLUDE: glibc
# 2.38+ headers redirect strtol and friends to __isoc23_* symbols that older
# glibc (such as pixi clang's default sysroot) does not export.
: "${AMREX_LOOM_SYSROOT:=/}"
# Where the three clones live, and where builds go.
: "${LOOM_WORK:=$PWD}"
: "${CXX_SRC:=$LOOM_WORK/cplusplus}"      # branch loom-single-source, pinned at 0caf64b6ee
: "${HRX_SRC:=$LOOM_WORK/hrx}"            # branch loom-single-source, pinned at a74773957c
: "${AMREX_SRC:=$LOOM_WORK/amrex}"        # branch loom-single-source
: "${BUILD_JOBS:=8}"

# --- Derived ----------------------------------------------------------------
export GCC_INSTALL_DIR GLIBC_INCLUDE AMREX_LOOM_SYSROOT LOOM_WORK CXX_SRC HRX_SRC AMREX_SRC BUILD_JOBS
export TOOLS_ENV=$AMREX_SRC/Tools/Loom/single-source/.pixi/envs/default
export HRX_BUILD=$LOOM_WORK/build/hrx
export AMREX_BUILD=$LOOM_WORK/build/amrex
export PATH=$TOOLS_ENV/bin:$PATH
export TMPDIR=${TMPDIR:-$LOOM_WORK/tmp}
mkdir -p "$TMPDIR"
# clang from the pixi environment builds the tools and links executables.
export CLANG=$TOOLS_ENV/bin/clang
export CLANGXX=$TOOLS_ENV/bin/clang++
# HRX runs on its native AMD GPU driver over libamdf (no ROCm runtime).
export HRX_GPU_DRIVER=amdgpu

# AMReX's LOOM toolchain (amrex/Tools/Loom/LoomToolchain.cmake): loom-cxx
# compiles, clang++ links, libc++ for both.
export AMREX_LOOM_CXX_DRIVER=$HRX_BUILD/loom/src/loom/tools/loom-cxx/loom-cxx
export AMREX_LOOM_LINK_CXX=$CLANGXX
export AMREX_LOOM_LIBCXX_ROOT=$TOOLS_ENV
export AMREX_LOOM_GCC_INSTALL_DIR=$GCC_INSTALL_DIR
export AMREX_LOOM_GLIBC_INCLUDE=$GLIBC_INCLUDE
# Device target and libhrx, for AMReX's configure. loom-cxx finds Loom's
# loomcxx/ headers itself.
export AMREX_LOOM_TARGET=gfx90a
export AMREX_LOOM_HRX_LIB_DIR=$HRX_BUILD/libhrx/src/libhrx
export AMREX_LOOM_HRX_INCLUDE_DIR=$HRX_SRC/libhrx/include
