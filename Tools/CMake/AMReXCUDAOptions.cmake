include_guard(GLOBAL)

#
# Define a macro to print active options
#
macro (cuda_print_option _var)
   if (${_var})
      message( STATUS "   ${_var}")
   endif ()
endmacro ()


#
#  CUDA-related options
#
message(STATUS "Enabled CUDA options:")

# Resolve and report the target CUDA architecture(s). CMAKE_CUDA_ARCHITECTURES was set in
# AMReXCUDAArchs.cmake (before enable_language, honoring the user's hints); the aliases
# "native", "all" and "all-major" are resolved by CMake during compiler detection and
# exposed via CMAKE_CUDA_ARCHITECTURES_{NATIVE,ALL,ALL_MAJOR}. We expand them here because
# AMReX needs the concrete architectures: to export them to downstream projects, to check
# them against the minimum compute capability AMReX supports, and because CUDA device LTO
# is not applied to the alias forms (they compile to ordinary SASS, not to LTO objects).
set(_amrex_cuda_archs "${CMAKE_CUDA_ARCHITECTURES}")
set(_amrex_cuda_archs_alias FALSE)
if (CMAKE_CUDA_ARCHITECTURES STREQUAL "native")
   if (CMAKE_CUDA_ARCHITECTURES_NATIVE MATCHES "^[0-9]")
      # record the concrete architecture(s) instead of "native", so that the value we
      # compile with is also the value we report and export to downstream projects.
      # CMake reports them as "<NN>-real", i.e. SASS only; keep the plain integer form
      # so that PTX is embedded as well and the library can still JIT onto a newer GPU,
      # which is what both "-arch=native" and the historical "Auto" default did.
      string(REPLACE "-real" "" _amrex_cuda_archs "${CMAKE_CUDA_ARCHITECTURES_NATIVE}")
      set(_amrex_cuda_archs_alias TRUE)
      message(STATUS "   CUDA architectures: native -> ${_amrex_cuda_archs}")
   else ()
      # no GPU was visible at configure time: CMAKE_CUDA_ARCHITECTURES_NATIVE holds an
      # error message rather than an architecture, or nothing at all. Keep "native" and
      # let the CUDA compiler resolve it, which only works if a GPU is visible when
      # compiling (a Clang CUDA compiler fails right away instead).
      message(STATUS "   CUDA architectures: native (unresolved)")
      message(WARNING
         "CUDA architecture 'native' was requested, but no CUDA device was found at "
         "configure time. The architecture selection is left to the CUDA compiler and the "
         "build will fail if no GPU is visible then either. On machines without a GPU, "
         "e.g. HPC login nodes or CI runners, pass an explicit architecture instead, such "
         "as -DCMAKE_CUDA_ARCHITECTURES=80.")
   endif ()
elseif (CMAKE_CUDA_ARCHITECTURES STREQUAL "all" AND CMAKE_CUDA_ARCHITECTURES_ALL)
   # unlike "native", the alias lists keep CMake's "<NN>-real" entries: embedding PTX for
   # the newest architecture only is what "all"/"all-major" mean
   set(_amrex_cuda_archs "${CMAKE_CUDA_ARCHITECTURES_ALL}")
   set(_amrex_cuda_archs_alias TRUE)
   message(STATUS "   CUDA architectures: all -> ${_amrex_cuda_archs}")
elseif (CMAKE_CUDA_ARCHITECTURES STREQUAL "all-major" AND CMAKE_CUDA_ARCHITECTURES_ALL_MAJOR)
   set(_amrex_cuda_archs "${CMAKE_CUDA_ARCHITECTURES_ALL_MAJOR}")
   set(_amrex_cuda_archs_alias TRUE)
   message(STATUS "   CUDA architectures: all-major -> ${_amrex_cuda_archs}")
else ()
   message(STATUS "   CUDA architectures: ${_amrex_cuda_archs}")
endif ()

# expanding an alias on an older CUDA toolkit can bring in architectures below the minimum
# compute capability AMReX supports (explicit values were already checked in AMReXCUDAArchs)
amrex_filter_cuda_archs(_amrex_cuda_archs ${_amrex_cuda_archs_alias})

# The architecture(s) AMReX is compiled for: used for the AMReX targets themselves, for
# dependent test/tutorial targets (setup_target_for_cuda_compilation) and exported to
# downstream projects (AMReXConfig.cmake.in).
set(AMREX_CUDA_ARCHS "${_amrex_cuda_archs}" CACHE INTERNAL "CUDA archs AMReX is built for")
unset(_amrex_cuda_archs)
unset(_amrex_cuda_archs_alias)

option(AMReX_CUDA_FASTMATH "Enable CUDA fastmath" ON)  # Note: inconsistent with AMReX_FASTMATH defaults
cuda_print_option( AMReX_CUDA_FASTMATH )

set(AMReX_CUDA_MAXREGCOUNT "255" CACHE STRING
   "Limit the maximum number of registers available" )
message( STATUS "   AMReX_CUDA_MAXREGCOUNT = ${AMReX_CUDA_MAXREGCOUNT}")

# if this works well and does not add too much compile-time we should enable it by default
option(AMReX_CUDA_LTO "Enable CUDA link-time-optimization (requires AMReX_GPU_RDC)" OFF)
cuda_print_option(AMReX_CUDA_LTO)

# device link-time optimization needs relocatable device code: without it every
# translation unit is already device linked on its own and there is nothing left to
# optimize across. Do not silently produce a non-LTO build.
if (AMReX_CUDA_LTO AND NOT AMReX_GPU_RDC)
   message(FATAL_ERROR
      "AMReX_CUDA_LTO requires relocatable device code. Configure with -DAMReX_GPU_RDC=ON "
      "(the default) or with -DAMReX_CUDA_LTO=OFF.")
endif ()

# device LTO is applied per architecture, so it needs the concrete list: an architecture
# alias that could not be expanded above would compile to ordinary device code without any
# hint that link-time optimization was dropped
if (AMReX_CUDA_LTO AND AMREX_CUDA_ARCHS MATCHES "^(native|all|all-major)$")
   message(WARNING
      "AMReX_CUDA_LTO is ON but the CUDA architectures could not be resolved from "
      "'${AMREX_CUDA_ARCHS}', so no device link-time optimization is applied. Pass explicit "
      "architectures instead, e.g. -DCMAKE_CUDA_ARCHITECTURES=80.")
endif ()

# this warns on a typical user bug when developing on (forgiving) Power9 machines (e.g. Summit)
option(AMReX_CUDA_WARN_CAPTURE_THIS "Warn if a CUDA lambda captures a class' this" ON)
# no code should ever ship -Werror, but one can turn this on manually in CI if one likes
option(AMReX_CUDA_ERROR_CAPTURE_THIS "Error if a CUDA lambda captures a class' this" OFF)
cuda_print_option(AMReX_CUDA_WARN_CAPTURE_THIS)
cuda_print_option(AMReX_CUDA_ERROR_CAPTURE_THIS)

option(AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL
       "Error if a CUDA host function is called from a host device function" OFF)
cuda_print_option(AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL)

option(AMReX_CUDA_PTX_VERBOSE "Verbose code generation statistics in ptxas" OFF)
cuda_print_option(AMReX_CUDA_PTX_VERBOSE)

option(AMReX_CUDA_COMPILATION_TIMER "Generate CSV table with time for each compilation phase" OFF)
cuda_print_option(AMReX_CUDA_COMPILATION_TIMER)

# Default on for CMAKE_BUILD_TYPE "Debug":
# In the past, this often did not compile at all, was very sensitive to further set options, or
# compiled super slowly;
# in some cases, such as recursive function usage, apps need to increase
# `cudaLimitStackSize` in order to not stack overflow with device debug symbols
# (this costs some extra DRAM).
# Nonetheless, for CUDA approx. 11.0+, we see the opposite: we have very slow Debug builds with CUDA if
# we do not activate -G for some LinearSolvers/MLMG objects. Thus, we default-on now to -G in
# Debug builds.
if ( "${CMAKE_BUILD_TYPE}" MATCHES "Debug" )
    set(AMReX_CUDA_DEBUG_DEFAULT ON)
else ()
    set(AMReX_CUDA_DEBUG_DEFAULT OFF)
endif ()
option(AMReX_CUDA_DEBUG "Generate debug information for device code (optimizations: off)" ${AMReX_CUDA_DEBUG_DEFAULT})
cuda_print_option(AMReX_CUDA_DEBUG)

# both are performance-neutral debug symbols
option(AMReX_CUDA_SHOW_LINENUMBERS "Generate line-number information (optimizations: on)" ON)
cuda_print_option(AMReX_CUDA_SHOW_LINENUMBERS)

# off by default for the CUDA versions we support (12.2+)
# https://github.com/AMReX-Codes/amrex/issues/3215
# Nvidia Bug ID: 4088095
option(AMReX_CUDA_SHOW_CODELINES "Generate source information in PTX (optimizations: on)" OFF)
cuda_print_option(AMReX_CUDA_SHOW_CODELINES)

option(AMReX_CUDA_BACKTRACE "Generate host function symbol names (better cuda-memcheck)" ${AMReX_CUDA_DEBUG})
cuda_print_option(AMReX_CUDA_BACKTRACE)

option(AMReX_CUDA_KEEP_FILES "Keep intermediately generated files (folder: nvcc_tmp)" OFF)
cuda_print_option(AMReX_CUDA_KEEP_FILES)

option(AMReX_CUDA_OBJDIR_AS_TEMPDIR "Place intermediate files in object file folder" OFF)
cuda_print_option(AMReX_CUDA_OBJDIR_AS_TEMPDIR)
