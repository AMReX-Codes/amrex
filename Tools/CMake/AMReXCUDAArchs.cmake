include_guard(GLOBAL)

#
# Resolve the target CUDA architecture(s) into the standard CMAKE_CUDA_ARCHITECTURES
# cache variable.
#
# This must be included BEFORE enable_language(CUDA) so that the CMP0104 NEW behavior
# does not initialize CMAKE_CUDA_ARCHITECTURES to a compiler-default architecture, which
# would defeat both our "native" default and the precedence of the legacy AMReX hints.
#
# We standardize on the CMake-native interface:
#   - CMAKE_CUDA_ARCHITECTURES cache variable (e.g. -DCMAKE_CUDA_ARCHITECTURES=80)
#   - CUDAARCHS environment variable (CMake's initializer for the above)
#
# The legacy AMReX hints are still honored (with a deprecation warning), keeping their
# historical precedence:
#   1. AMReX_CUDA_ARCH cache variable  (deprecated -> CMAKE_CUDA_ARCHITECTURES)
#   2. AMREX_CUDA_ARCH  environment variable (deprecated -> CUDAARCHS)
#   3. CMAKE_CUDA_ARCHITECTURES cache variable
#   4. CUDAARCHS environment variable
#   5. "native" (build for the locally installed GPU; requires a GPU at configure time)
#

if (DEFINED AMReX_CUDA_ARCH)
   message(WARNING
      "AMReX_CUDA_ARCH is deprecated; set the standard CMAKE_CUDA_ARCHITECTURES instead "
      "(e.g. -DCMAKE_CUDA_ARCHITECTURES=80, or 'native').")
   set(_amrex_cuda_archs "${AMReX_CUDA_ARCH}")
elseif (DEFINED ENV{AMREX_CUDA_ARCH})
   message(WARNING
      "The AMREX_CUDA_ARCH environment variable is deprecated; use the standard CUDAARCHS "
      "environment variable instead (or -DCMAKE_CUDA_ARCHITECTURES=...).")
   set(_amrex_cuda_archs "$ENV{AMREX_CUDA_ARCH}")
elseif (DEFINED CMAKE_CUDA_ARCHITECTURES)
   set(_amrex_cuda_archs "${CMAKE_CUDA_ARCHITECTURES}")
elseif (DEFINED ENV{CUDAARCHS})
   set(_amrex_cuda_archs "$ENV{CUDAARCHS}")
else ()
   set(_amrex_cuda_archs "native")
endif ()

#
# Back-compatibility normalization for the legacy hints:
#   - "Auto" (the historical default) -> "native"
#   - strip the decimal dot per entry ("8.0" -> "80", "7.5" -> "75")
# Everything the modern CUDA_ARCHITECTURES property understands is passed through unchanged:
#   native / all / all-major / integers / <NN>a / -real / -virtual suffixes.
#
set(_amrex_cuda_archs_norm)
foreach (_arch IN LISTS _amrex_cuda_archs)
   if (_arch STREQUAL "Auto")
      set(_arch "native")
   else ()
      string(REPLACE "." "" _arch "${_arch}")
   endif ()
   list(APPEND _amrex_cuda_archs_norm "${_arch}")
endforeach ()

set(CMAKE_CUDA_ARCHITECTURES "${_amrex_cuda_archs_norm}" CACHE STRING
   "CUDA architectures: 'native', 'all-major', or e.g. 80;90a (see CMake CUDA_ARCHITECTURES)" FORCE)

unset(_amrex_cuda_archs)
unset(_amrex_cuda_archs_norm)
