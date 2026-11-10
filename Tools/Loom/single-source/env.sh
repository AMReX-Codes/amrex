# Environment for single-source LOOM builds of AMReX. Source it before each
# step (build-cxx.sh, build-hrx.sh, build-amrex.sh, run-tests.sh do this).
#
# Edit the first block for your machine; everything else derives from it.

# --- Machine-specific -------------------------------------------------------
# GCC 14 install directory (C++23 libstdc++ for building the tools, glibc crt
# files for linking, compiler headers such as stddef.h for the single-source
# compiles).
: "${GCC_INSTALL_DIR:=/usr/lib64/gcc/x86_64-suse-linux/14}"
# glibc headers seen by the device pass (and by the host via the default sysroot).
: "${GLIBC_INCLUDE:=/usr/include}"
# Sysroot for compiling and linking host code. It must hold the same glibc as
# GLIBC_INCLUDE: glibc 2.38+ headers redirect strtol and friends to
# __isoc23_* symbols that older glibc (such as pixi clang's default sysroot)
# does not export.
: "${AMREX_LOOM_SYSROOT:=/}"
# Where the three clones live, and where builds go.
: "${LOOM_WORK:=$PWD}"
: "${CXX_SRC:=$LOOM_WORK/cplusplus}"      # branch loom-single-source, pinned at d6b18ef6c5
: "${HRX_SRC:=$LOOM_WORK/hrx}"            # branch loom-single-source, pinned at 4dcf500941
: "${AMREX_SRC:=$LOOM_WORK/amrex}"        # branch loom-single-source
: "${BUILD_JOBS:=8}"

# --- Derived ----------------------------------------------------------------
export GCC_INSTALL_DIR GLIBC_INCLUDE AMREX_LOOM_SYSROOT LOOM_WORK CXX_SRC HRX_SRC AMREX_SRC BUILD_JOBS
export TOOLS_ENV=$AMREX_SRC/Tools/Loom/single-source/.pixi/envs/default
export CXX_BUILD=$LOOM_WORK/build/cxx
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

# amrex-loom-cxx (the C++ compiler AMReX is configured with).
export AMREX_LOOM_HOST_CXX=$AMREX_SRC/Tools/Loom/amrex-loom-host-cxx
export AMREX_LOOM_IMPORT=$HRX_BUILD/loom/src/loom/tools/loom-import-cxx/loom-import-cxx
export AMREX_LOOM_COMPILE=$HRX_BUILD/loom/src/loom/tools/loom-compile/loom-compile
export AMREX_LOOM_TARGET=gfx90a
export AMREX_LOOM_LOOMCXX_DIR=$LOOM_WORK/build/loomcxx-include
export AMREX_LOOM_SYSTEM_INCLUDES=$TOOLS_ENV/include/c++/v1:$GCC_INSTALL_DIR/include:$GLIBC_INCLUDE
export AMREX_LOOM_HRX_LIB_DIR=$HRX_BUILD/libhrx/src/libhrx
export AMREX_LOOM_HRX_INCLUDE_DIR=$HRX_SRC/libhrx/include

# amrex-loom-host-cxx (cxx for compiles, clang++ for links, libc++ for both).
export AMREX_LOOM_CXX=$CXX_BUILD/src/frontend/cxx
export AMREX_LOOM_LIBCXX_ROOT=$TOOLS_ENV
export AMREX_LOOM_GCC_INSTALL_DIR=$GCC_INSTALL_DIR
export AMREX_LOOM_LINK_CXX=$CLANGXX
