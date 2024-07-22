include_guard(GLOBAL)

include(CMakeDependentOption)

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

set(AMReX_CUDA_ARCH_DEFAULT "Auto")
if(DEFINED CMAKE_CUDA_ARCHITECTURES)
   set(AMReX_CUDA_ARCH_DEFAULT "${CMAKE_CUDA_ARCHITECTURES}")
endif ()
if(DEFINED ENV{AMREX_CUDA_ARCH})
   set(AMReX_CUDA_ARCH_DEFAULT "$ENV{AMREX_CUDA_ARCH}")
endif()
set(AMReX_CUDA_ARCH ${AMReX_CUDA_ARCH_DEFAULT} CACHE STRING "CUDA architecture (Use 'Auto' for automatic detection)")

option(AMReX_CUDA_FASTMATH "Enable CUDA fastmath" ON)
cuda_print_option( AMReX_CUDA_FASTMATH )

set(AMReX_CUDA_MAXREGCOUNT "255" CACHE STRING
   "Limit the maximum number of registers available" )
message( STATUS "   AMReX_CUDA_MAXREGCOUNT = ${AMReX_CUDA_MAXREGCOUNT}")

# if this works well and does not add too much compile-time we should enable it by default
cmake_dependent_option(AMReX_CUDA_LTO "Enable CUDA link-time-optimization" OFF
   "CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0.167" OFF)
cuda_print_option(AMReX_CUDA_LTO)

# this warns on a typical user bug when developing on (forgiving) Power9 machines (e.g. Summit)
cmake_dependent_option(AMReX_CUDA_WARN_CAPTURE_THIS
   "Warn if a CUDA lambda captures a class' this" ON
   "CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0.194" OFF)
# no code should ever ship -Werror, but one can turn this on manually in CI if one likes
cmake_dependent_option(AMReX_CUDA_ERROR_CAPTURE_THIS
   "Error if a CUDA lambda captures a class' this" OFF
   "CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0.194" OFF)
cuda_print_option(AMReX_CUDA_WARN_CAPTURE_THIS)
cuda_print_option(AMReX_CUDA_ERROR_CAPTURE_THIS)

option(AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL
       "Error if a CUDA host function is called from a host device function" OFF)
cuda_print_option(AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL)

# makes things more robust for -Xcompiler pre-fixing unknown nvcc flags
# note: available with NVCC 10.2.89+; default in CMake 3.17.0+ for supporting NVCCs
#       https://gitlab.kitware.com/cmake/cmake/-/blob/v3.17.0/Modules/Compiler/NVIDIA-CUDA.cmake
cmake_dependent_option(CUDA_FORWARD_UNKNOWN_FLAGS_HOST
   "Forward unknown NVCC flags to the host compiler" ON
   "CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 10.2.89;CMAKE_VERSION VERSION_LESS 3.17" OFF)

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

# https://github.com/AMReX-Codes/amrex/issues/3215
# Nvidia Bug ID: 4088095
if (CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 12.1)
    set(AMReX_CUDA_SHOW_CODELINES_DEFAULT OFF)
else()
    set(AMReX_CUDA_SHOW_CODELINES_DEFAULT ON)
endif()
option(AMReX_CUDA_SHOW_CODELINES "Generate source information in PTX (optimizations: on)" AMReX_CUDA_SHOW_CODELINES_DEFAULT)
cuda_print_option(AMReX_CUDA_SHOW_CODELINES)

option(AMReX_CUDA_BACKTRACE "Generate host function symbol names (better cuda-memcheck)" ${AMReX_CUDA_DEBUG})
cuda_print_option(AMReX_CUDA_BACKTRACE)

option(AMReX_CUDA_KEEP_FILES "Keep intermediately generated files (folder: nvcc_tmp)" OFF)
cuda_print_option(AMReX_CUDA_KEEP_FILES)
