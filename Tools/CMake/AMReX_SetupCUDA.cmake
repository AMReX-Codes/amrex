#
# Setup the CUDA enviroment
#
#  Authors: Michele Rosso, Axel Huebl
#  Date   : April 4, 2019
#
#
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


get_property(_lang GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT ("CUDA" IN_LIST _lang ))
    message(WARNING "AMReX_SetupCUDA will not be processed because CUDA language has not been enabled.")
    return()
endif ()

#
# Makes sure the CUDA host compiler and CXX compiler are the same.
# CMake lets you decide which host compiler to use via the env variable
# CUDAHOSTCXX and the CMake variable CMAKE_CUDA_HOST_COMPILER.
# For the time being we force the CUDA host compiler to be the C++ compiler.
#
# Note: just comparing the CMAKE_..._COMPILER vars is not sufficient and raises
#       false negatives on e.g. /usr/bin/g++-8 and /usr/bin/c++
# Note: blocked by https://gitlab.kitware.com/cmake/cmake/-/issues/20901
#
if (CMAKE_CUDA_HOST_COMPILER)
  if (NOT "${CMAKE_CXX_COMPILER}" STREQUAL "${CMAKE_CUDA_HOST_COMPILER}")
    if (NOT "$ENV{CUDAHOSTCXX}" STREQUAL "" OR NOT "$ENV{CXX}" STREQUAL "")
      message(WARNING "CUDA host compiler "
                      "(${CMAKE_CUDA_HOST_COMPILER}) "
                      "does not match the C++ compiler "
                      "(${CMAKE_CXX_COMPILER})! "
                      "Consider setting the CXX and CUDAHOSTCXX environment "
                      "variables.")
    endif ()
  endif ()
endif ()

#
#  CUDA-related options
#
include(AMReXCUDAOptions)

#
# Error if NVCC is too old
#
if (CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "9.0")
   message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION}."
      "This is unsupported. Please use CUDA toolkit version 9.0 or newer.")
endif ()

#
# Find cuda flags for target architecture. If CUDA_ARCH is not set by the user,
# autodetection is enabled
#
include(FindCUDA/select_compute_arch)
cuda_select_nvcc_arch_flags(_nvcc_arch_flags ${AMReX_CUDA_ARCH})

#
# Remove unsupported architecture: anything less the 3.5 must go
#
string(REPLACE "-gencode;" "-gencode=" _nvcc_arch_flags "${_nvcc_arch_flags}")

foreach (_item IN LISTS _nvcc_arch_flags)
   # Match one time the regex [0-9]+.
   # [0-9]+ means any number between 0 and 9 will be matched one or more times (option +)
   string(REGEX MATCH "[0-9]+" _cuda_compute_capability "${_item}")

   if (_cuda_compute_capability LESS 35)
      message(STATUS "Ignoring unsupported CUDA architecture ${_cuda_compute_capability}")
      list(REMOVE_ITEM _nvcc_arch_flags ${_item})
   endif ()

endforeach ()


if (AMReX_CUDA_LTO)
    # we replace
    #   -gencode=arch=compute_NN,code=sm_NN
    # with
    #   -gencode=arch=compute_NN,code=lto_NN
    set(_nvcc_arch_flags_org ${_nvcc_arch_flags})
    foreach (_item IN LISTS _nvcc_arch_flags_org)
       string(REGEX MATCH "[0-9]+" _cuda_compute_capability "${_item}")
       string(REPLACE "code=sm_${_cuda_compute_capability}"
                      "code=lto_${_cuda_compute_capability}"
              _nvcc_arch_flags "${_nvcc_arch_flags}")
    endforeach ()
endif ()

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

return ()

# We gotta set CUDA flags globally since there is no other way at this time to pass CUDA flags to
# device linking stage
if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Windows" ) )
   # CUDA only supports 64bit builds on windows ( 32bit builds are deprecated ).
   # Thus the option "--machine 64" is being set by the msbuild configuration.
   # For Linux and MAC, we need to enforce that manually
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -m64")
endif ()
string(APPEND CMAKE_CUDA_FLAGS " --expt-relaxed-constexpr --expt-extended-lambda")
string(APPEND CMAKE_CUDA_FLAGS " -Wno-deprecated-gpu-targets ${NVCC_ARCH_FLAGS}")
string(APPEND CMAKE_CUDA_FLAGS " -maxrregcount=${AMReX_CUDA_MAXREGCOUNT}")

# This is to work around a bug with nvcc, see: https://github.com/kokkos/kokkos/issues/1473
string(APPEND CMAKE_CUDA_FLAGS " -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")

if (AMReX_CUDA_FASTMATH)
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")
endif ()

#
# Print numbers for warnings and errors
#
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --display_error_number")


#
# CUDA specific warnings
#
if (AMReX_CUDA_WARN_CAPTURE_THIS)
    string(APPEND CMAKE_CUDA_FLAGS " --Wext-lambda-captures-this")
endif()
if (AMReX_CUDA_ERROR_CAPTURE_THIS)
    # note: prefer double-dash --Werror!
    # https://github.com/ccache/ccache/issues/598
    string(APPEND CMAKE_CUDA_FLAGS " --Werror ext-lambda-captures-this")
endif()

if (AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL)
    string(APPEND CMAKE_CUDA_FLAGS " --Werror cross-execution-space-call")
endif()

#
# Forward unknown NVCC flags to the host compiler
#
if (CUDA_FORWARD_UNKNOWN_FLAGS_HOST)
    string(APPEND CMAKE_CUDA_FLAGS " --forward-unknown-to-host-compiler")
endif()

#
# Code generation
#
if (AMReX_CUDA_PTX_VERBOSE)
    string(APPEND CMAKE_CUDA_FLAGS " --ptxas-options=-v")
endif()

# keep intermediately generated files
if (AMReX_CUDA_KEEP_FILES)
    make_directory("${PROJECT_BINARY_DIR}/nvcc_tmp")
    string(APPEND CMAKE_CUDA_FLAGS " --keep --keep-dir ${PROJECT_BINARY_DIR}/nvcc_tmp")
endif ()

# compilation timings
if (AMReX_CUDA_COMPILATION_TIMER)
    file(REMOVE "${PROJECT_BINARY_DIR}/nvcc_timings.csv")
    string(APPEND CMAKE_CUDA_FLAGS " --time ${PROJECT_BINARY_DIR}/nvcc_timings.csv")
endif ()

#
# Debugging
#
if (AMReX_CUDA_DEBUG)
    # is this unsupported with MSVC?
    string(APPEND CMAKE_CUDA_FLAGS " -G")
endif()

if (AMReX_CUDA_SHOW_LINENUMBERS AND NOT AMReX_CUDA_DEBUG)
    # nvcc warning : '--device-debug (-G)' overrides '--generate-line-info (-lineinfo)'
    string(APPEND CMAKE_CUDA_FLAGS " --generate-line-info")
endif ()
if (AMReX_CUDA_SHOW_CODELINES)
    string(APPEND CMAKE_CUDA_FLAGS " --source-in-ptx")
endif ()

if (AMReX_CUDA_BACKTRACE)
    if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
        string(APPEND CMAKE_CUDA_FLAGS " -Xcompiler /Zi") # comes with Debug & RelWithDebInfo
    else ()
        string(APPEND CMAKE_CUDA_FLAGS " -Xcompiler -rdynamic")
    endif ()
endif ()
