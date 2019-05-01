#
# Setup the CUDA enviroment
#
#  Author: Michele Rosso
#  Date  : April 4, 2019
#
# 
include_guard(GLOBAL)

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

# 
# Enable CUDA language here
# 
enable_language(CUDA)

#
#  CUDA-related options
#
message(STATUS "CUDA options:")

set(CUDA_ARCH "Auto" CACHE STRING "CUDA architecture (Use 'Auto' for automatic detection)")

option( ENABLE_CUDA_FASTMATH "Enable CUDA fastmath" ON )
message(STATUS "   ENABLE_CUDA_FASTMATH = ${ENABLE_CUDA_FASTMATH}")

set(CUDA_MAX_THREADS "256" CACHE STRING
   "Maximum number of CUDA threads per block" )
message(STATUS "   CUDA_MAX_THREADS = ${CUDA_MAX_THREADS}")

set(CUDA_MAXREGCOUNT "255" CACHE STRING
   "Limit the maximum number of registers available" )
message(STATUS "   CUDA_MAXREGCOUNT = ${CUDA_MAXREGCOUNT}")

# 
# Error if NVCC is too old
#
If (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "8.0")
   message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION}."
      "This is unsupported. Please use CUDA toolkit version 8.0 or newer.")
endif ()

# 
# Find cuda flags for target architecture. If CUDA_ARCH is not set by the user,
# autodetection is enabled
#
include(FindCUDA/select_compute_arch)
cuda_select_nvcc_arch_flags(_nvcc_arch_flags ${CUDA_ARCH})

#
# Remove unsupported architecture: anything less the 6.0 must go
# 
string(REPLACE "-gencode;" "-gencode " _nvcc_arch_flags "${_nvcc_arch_flags}")

foreach (_item IN LISTS _nvcc_arch_flags)
   # Match one time the regex [0-9]+.
   # [0-9]+ means any number between 0 and 9 will be matched one or more times (option +)
   string(REGEX MATCH "[0-9]+" _cuda_compute_capability "${_item}")
   
   if (_cuda_compute_capability LESS 60)
      message(STATUS "Ignoring unsupported CUDA architecture ${_cuda_compute_capability}")
      list(REMOVE_ITEM _nvcc_arch_flags ${_item})
   endif ()
   
endforeach ()

if (NOT _nvcc_arch_flags)
   message(FATAL_ERROR "the given target CUDA architectures are not supported by AMReX")
endif ()

#
# Set architecture-dependent flags
# 
string(REPLACE ";" " " _nvcc_arch_flags "${_nvcc_arch_flags}")
set(NVCC_ARCH_FLAGS ${_nvcc_arch_flags} CACHE INTERNAL "CUDA architecture-dependent flags")
unset(_nvcc_arch_flags)


# CUDA compiler is in the form CUDA_HOME/bin/compiler-name
# Remove bin/compiler-name to get CUDA HOME
get_filename_component(_cuda_home ${CMAKE_CUDA_COMPILER} DIRECTORY) # remove compiler from path
get_filename_component(_cuda_home ${_cuda_home} DIRECTORY) # remove bin/ from path
set( CUDA_HOME ${_cuda_home} CACHE INTERNAL "Path to CUDA library")
unset(_cuda_home)

#
# Set CUDA version variables
#
string( REPLACE "." ";" _cuda_compiler_version ${CMAKE_CUDA_COMPILER_VERSION})
list( GET _cuda_compiler_version 0 _nvcc_version_major )
list( GET _cuda_compiler_version 1 _nvcc_version_minor )

set(NVCC_VERSION_MAJOR "${_nvcc_version_major}" CACHE INTERNAL "CUDA compiler version (major)")
set(NVCC_VERSION_MINOR "${_nvcc_version_minor}" CACHE INTERNAL "CUDA compiler version (minor)")

#
# Error if EB is used with CUDA 9.2
#
set(_cuda_version_truncated "${NVCC_VERSION_MAJOR}.${NVCC_VERSION_MINOR}")
if (ENABLE_EB AND ( _cuda_version_truncated VERSION_EQUAL "9.2") )
   message(FATAL_ERROR "EB component of AMReX is not compatible with CUDA 9.2")
endif ()
unset(_cuda_version_truncated)

# We gotta set CUDA flags globally since there is no other way at this time to pass CUDA flags to
# device linking stage
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --expt-extended-lambda")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets -m64 ${NVCC_ARCH_FLAGS}")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -maxrregcount=${CUDA_MAXREGCOUNT}")

if (ENABLE_CUDA_FASTMATH)
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")
endif ()

