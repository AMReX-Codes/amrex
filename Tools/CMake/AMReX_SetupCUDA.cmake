#
# Collection of macros to setup the CUDA enviroment
#
#  Author: Michele Rosso
#  Date  : February 14, 2019
#
# 
include_guard(GLOBAL)

if (${CMAKE_VERSION} VERSION_LESS "3.14.0") 
   message(FATAL_ERROR "CUDA support requires CMake >= 3.14.0")
endif ()

#
# Makes sure the CUDA host compiler and CXX compiler are the same.
# CMake let you decide which host compiler to use via the env variable
# CUDAHOSTCXX and the CMake variable CMAKE_CUDA_HOST_COMPILER.
# For the time being we force the CUDA host compiler to be the C++ compiler.
#
# DO THIS BEFORE CALLING enable_language(CUDA)
#
if ( ( CMAKE_CUDA_HOST_COMPILER AND NOT ("${CMAKE_CUDA_HOST_COMPILER}" STREQUAL "${CMAKE_CXX_COMPILER}") )
      OR  ( NOT ("$ENV{CUDAHOSTCXX}" STREQUAL "") ) )
   message(WARNING
      "User-defined CUDA host compiler does not match C++ compiler: overwriting user setting.")
endif ()
unset(ENV{CUDAHOSTCXX})
set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "" FORCE)

# Enable CUDA language here
enable_language(CUDA)

# Find cuda flags for this platform
include(FindCUDA/select_compute_arch)
cuda_select_nvcc_arch_flags(_nvcc_arch_flags "Auto")
string(REPLACE ";" " " _nvcc_arch_flags "${_nvcc_arch_flags}")
set(NVCC_ARCH_FLAGS ${_nvcc_arch_flags} CACHE INTERNAL "Architecture-dependent flags")

# CUDA compiler is in the form CUDA_HOME/bin/compiler-name
# Remove bin/compiler-name to get CUDA HOME
get_filename_component(_cuda_home ${CMAKE_CUDA_COMPILER} DIRECTORY) # remove compiler from path
get_filename_component(_cuda_home ${_cuda_home} DIRECTORY) # remove bin/ from path
set( CUDA_HOME ${_cuda_home} CACHE INTERNAL "Path to CUDA library")

# Extrapolate CUDA arch number (compute capabilities number)
string(REGEX MATCHALL "compute_[0-9]+"  _cuda_arch_tmp "${_nvcc_arch_flags}")
list(REMOVE_DUPLICATES _cuda_arch_tmp)
set(_cuda_arch "")
foreach (_item ${_cuda_arch_tmp})
   string(REPLACE "compute_" "" _arch ${_item} )
   list(APPEND _cuda_arch ${_arch})      
endforeach ()

set(CUDA_ARCH "${_cuda_arch}" CACHE STRING "CUDA architecture version")

# Set CUDA versioning number
if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "8.0")
   message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION}."
      "This is unsupported. Please use CUDA toolkit version 8.0 or newer.")
endif ()

string( REPLACE "." ";" _cuda_compiler_version ${CMAKE_CUDA_COMPILER_VERSION})
list( GET _cuda_compiler_version 0 _nvcc_version_major )
list( GET _cuda_compiler_version 1 _nvcc_version_minor )

set(NVCC_VERSION_MAJOR "${_nvcc_version_major}" CACHE INTERNAL "CUDA compiler version (major)")
set(NVCC_VERSION_MINOR "${_nvcc_version_minor}" CACHE INTERNAL "CUDA compiler version (minor)")


set(_cuda_version_truncated "${NVCC_VERSION_MAJOR}.${NVCC_VERSION_MINOR}")
if (ENABLE_EB AND ( _cuda_version_truncated VERSION_EQUAL "9.2") )
   message(FATAL_ERROR "EB component of AMReX is not compatible with CUDA 9.2")
endif ()

# We gotta set CUDA flags globally since there is no other way at this time to pass CUDA flags to
# device linking stage
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --expt-extended-lambda")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets -m64 ${NVCC_ARCH_FLAGS}")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -maxrregcount=${CUDA_MAXREGCOUNT} --use_fast_math")

