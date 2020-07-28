#
# This module defines the INTERFACE target SYCL and its alias AMReX::SYCL.
# These targets provides build/link requirements for the SYCL language.
# For the time being, only dpc++  is supported
#

# Provide a cache variable for the dpc++ root directory installation by probing
# the compiler
execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version  OUTPUT_VARIABLE _tmp)
string(REGEX MATCH "InstalledDir: (.*)" _tmp "${_tmp}")
unset(_tmp)

get_filename_component(DPCPP_ROOT ${CMAKE_MATCH_1} DIRECTORY CACHE)

find_file(LIBSYCL_GLIBC_OBJ libsycl-glibc.o
   PATHS ${DPCPP_ROOT} ENV LD_LIBRARY_PATH
   PATH_SUFFIXES lib
   DOC "Full path to libsycl-glibc.o")

set(_cxx_clang "$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>>") # Only Clang for now

#
# SYCL and AMReX::SYCL targets
#
add_library(SYCL INTERFACE)
add_library(AMReX::SYCL ALIAS SYCL)

target_compile_features(SYCL INTERFACE cxx_std_17)

target_link_libraries(SYCL INTERFACE ${LIBSYCL_GLIBC_OBJ})

target_compile_options( SYCL
   INTERFACE
   $<${_cxx_clang}:-Wno-error=sycl-strict -fsycl -fsycl-unnamed-lambda>
   $<${_cxx_clang}:$<$<BOOL:${ENABLE_DPCPP_SPLIT_KERNEL}>:-fsycl-device-code-split=per_kernel>>)

# TODO: use $<LINK_LANG_AND_ID:> genex for CMake >=3.17
target_link_options( SYCL
   INTERFACE
   $<${_cxx_clang}:-fsycl -device-math-lib=fp32,fp64>
   $<${_cxx_clang}:$<$<BOOL:${ENABLE_DPCPP_SPLIT_KERNEL}>:-fsycl-device-code-split=per_kernel>>)

# temporary work-around for DPC++ beta08 bug
#   define "long double" as 64bit for C++ user-defined literals
#   https://github.com/intel/llvm/issues/2187
target_compile_options( SYCL
   INTERFACE
   $<${_cxx_clang}:-mlong-double-64>)

if (ENABLE_DPCPP_AOT)
   message(FATAL_ERROR "\nAhead-of-time (AOT) compilation support not available yet.\nRe-configure with ENABLE_DPCPP_AOT=OFF.")

   #
   # TODO: remove comments to enable AOT support when the time comes
   #       (main blocker: missing math library)
   #
   # if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
   #    ## TODO: use file(READ)
   #    execute_process( COMMAND cat /sys/devices/cpu/caps/pmu_name OUTPUT_VARIABLE _cpu_long_name )
   # else ()
   #    message(FATAL_ERROR "\nENABLE_DPCPP_AOT is not supported on ${CMAKE_SYSTEM_NAME}\n")
   # endif ()

   # string(STRIP "${_cpu_long_name}" _cpu_long_name)
   # if (_cpu_long_name STREQUAL "skylake")
   #    set(_cpu_short_name "skl")
   # elseif (_cpu_long_name STREQUAL "kabylake")
   #    set(_cpu_short_name "kbl")
   # elseif (_cpu_long_name STREQUAL "cascadelake")
   #    set(_cpu_short_name "cfl")
   # else ()
   #    message(FATAL_ERROR "\n Ahead-of-time compilation for CPU ${_cpu_long_name} is not yet supported\n"
   #       "Maybe set ENABLE_DPCPP_AOT to OFF?\n")
   # endif ()

   # target_compile_options( amrex
   #    PUBLIC
   #    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>,$<NOT:$<CONFIG:Debug>>>:
   #    -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xsycl-target-backend "-device ${_cpu_short_name}" >)
   # target_link_options(amrex PUBLIC -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xsycl-target-backend "-device ${_cpu_short_name}" )
   # unset(_cpu_long_name)
   # unset(_cpu_short_name)
endif ()


unset(_cxx_clang)
