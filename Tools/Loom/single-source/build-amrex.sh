#!/usr/bin/env bash
# Configure AMReX with the LOOM backend (single-source) and build the four
# single-GPU tests.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
# The toolchain file makes loom-cxx the C++ compiler and links with clang++
# and libc++ (settings from env.sh).
cmake -S "$AMREX_SRC" -B "$AMREX_BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE="$AMREX_SRC/Tools/Loom/LoomToolchain.cmake" -DCMAKE_C_COMPILER="$CLANG" \
  -DAMReX_GPU_BACKEND=LOOM -DAMReX_AMD_ARCH="$AMREX_LOOM_TARGET" -DAMReX_SPACEDIM=3 -DAMReX_PRECISION=DOUBLE \
  -DAMReX_MPI=OFF -DAMReX_OMP=OFF -DAMReX_FORTRAN=OFF -DAMReX_EB=OFF -DAMReX_PARTICLES=OFF \
  -DAMReX_FFT=OFF -DAMReX_LINEAR_SOLVERS=OFF -DAMReX_AMRLEVEL=ON \
  -DAMReX_ENABLE_TESTS=ON -DAMReX_TEST_TYPE=All \
  -DAMReX_LOOM_HRX_INCLUDE_DIR="$AMREX_LOOM_HRX_INCLUDE_DIR" -DAMReX_LOOM_HRX_LIB_DIR="$AMREX_LOOM_HRX_LIB_DIR"
cmake --build "$AMREX_BUILD" -j "$BUILD_JOBS" --target \
  Test_GPU_CNS_Sod_3d Test_GPU_ArenaLifetime_3d Test_GPU_ExternalStream_3d \
  Test_GPU_LoomCodeObjects_3d
