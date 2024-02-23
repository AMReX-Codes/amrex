include_guard(GLOBAL)

include(CMakeDependentOption)

if (CMAKE_VERSION  VERSION_GREATER_EQUAL 3.20)
   message(WARNING "\nAMReX_SetupCUDA is deprecated for CMake >= 3.20: it will not be processed!\n")
   return()
endif ()

get_property(_lang GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT ("CUDA" IN_LIST _lang ))
   message(WARNING "AMReX_SetupCUDA will not be processed because CUDA language has not been enabled.")
   return()
endif ()

#
# Check CUDA compiler and host compiler
#
include(AMReXUtils)
set_mininum_compiler_version(CUDA NVIDIA 10.0)
check_cuda_host_compiler()

#
#  CUDA-related options
#
include(AMReXCUDAOptions)

#
# Find cuda flags for target architecture.
#
set_nvcc_arch_flags(AMReX_CUDA_ARCH AMReX_CUDA_LTO)


# CUDA compiler is in the form CUDA_HOME/bin/compiler-name
# Remove bin/compiler-name to get CUDA HOME
get_filename_component(_cuda_home ${CMAKE_CUDA_COMPILER} DIRECTORY) # remove compiler from path
get_filename_component(_cuda_home ${_cuda_home} DIRECTORY) # remove bin/ from path
set( CUDA_HOME ${_cuda_home} CACHE INTERNAL "Path to CUDA library")
unset(_cuda_home)


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

# and another bug related to implicit returns with if constexpr, see: https://stackoverflow.com/questions/64523302/cuda-missing-return-statement-at-end-of-non-void-function-in-constexpr-if-fun
string(APPEND CMAKE_CUDA_FLAGS " -Xcudafe --diag_suppress=implicit_return_from_non_void_function")

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

if (CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.2)
   string(APPEND CMAKE_CUDA_FLAGS " --display-error-number --diag-error 20092")
endif ()
