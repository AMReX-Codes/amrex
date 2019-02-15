#
# Collection of macros to setup the CUDA enviroment
#
#  Author: Michele Rosso
#  Date  : February 14, 2019
#
# 


#
# MACRO: setup_cuda_host_compiler
#
# Makes sure the CUDA host compiler and CXX compiler are the same.
# CMake let you decide which host compiler to use via the env variable
# CUDAHOSTCXX and the CMake variable CMAKE_CUDA_HOST_COMPILER.
# For the time being we force the CUDA host compiler to be the C++ compiler.
#
# THIS MACRO MUST BE CALLED BEFORE ENABLE_LANGUAGE(CUDA)
#
#
macro (setup_cuda_host_compiler)
   if ( ( CMAKE_CUDA_HOST_COMPILER AND NOT ("${CMAKE_CUDA_HOST_COMPILER}" STREQUAL "${CMAKE_CXX_COMPILER}") )
         OR  ( NOT ("$ENV{CUDAHOSTCXX}" STREQUAL "") ) )
      message(WARNING
         "User-defined CUDA host compiler does not match C++ compiler: overwriting user setting.")
   endif ()
   unset(ENV{CUDAHOSTCXX})
   set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "" FORCE)
endmacro()

#
# MACRO: setup_cuda_environment
#
# It defines the following internal variables
#
#   NVCC_ARCH_FLAGS:  Architecture-dependent flags
#   CUDA_HOME:        Path to CUDA library      
# 
macro (setup_cuda_environment)

   include(FindCUDA/select_compute_arch)
   cuda_select_nvcc_arch_flags(_nvcc_arch_flags "Auto")
   set(NVCC_ARCH_FLAGS ${_nvcc_arch_flags} CACHE INTERNAL "Architecture-dependent flags")

   # CUDA compiler is in the form CUDA_HOME/bin/compiler-name
   # Remove bin/compiler-name to get CUDA HOME
   get_filename_component(_cuda_home ${CMAKE_CUDA_COMPILER} DIRECTORY) # remove compiler from path
   get_filename_component(_cuda_home ${_cuda_home} DIRECTORY) # remove bin/ from path
   set(CUDA_HOME ${_cuda_home} CACHE INTERNAL "Path to CUDA library")
   
   # Extrapolate CUDA arch number (compute capabilities number)
   string(REGEX MATCHALL "compute_[0-9]+"  _cuda_arch_tmp "${_nvcc_arch_flags}")
   list(REMOVE_DUPLICATES _cuda_arch_tmp)
   set(_cuda_arch "")
   foreach (_item ${_cuda_arch_tmp})
      string(REPLACE "compute_" "" _arch ${_item} )
      list(APPEND _cuda_arch ${_arch})      
   endforeach ()
   
   set(CUDA_ARCH "${_cuda_arch}" CACHE STRING "CUDA architecture version")
   print_option(CUDA_ARCH)
   
   string(REPLACE ";" " " tmp "${NVCC_ARCH_FLAGS}")
   set(CMAKE_CUDA_FLAGS "${tmp}")
   print(CMAKE_CUDA_FLAGS)
endmacro ()
