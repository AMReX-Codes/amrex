#
# This module defines the INTERFACE target SYCL and its alias AMReX::SYCL.
# These targets provides build/link requirements for the SYCL language.
# For the time being, only dpc++ is supported
#

# Provide a cache variable for the dpc++ root directory installation by probing
# the compiler
execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version  OUTPUT_VARIABLE _tmp)
string(REGEX MATCH "InstalledDir: (.*)" _tmp "${_tmp}")
unset(_tmp)

get_filename_component(DPCPP_ROOT ${CMAKE_MATCH_1} DIRECTORY CACHE)
message(STATUS "dpc++ root directory: ${DPCPP_ROOT}")

# Provide cache variable to identify the dpc++ version, including the "beta"
string(REGEX MATCH "[^//]*beta(.[^//])" DPCPP_VERSION "${DPCPP_ROOT}")
set(DPCPP_VERSION ${DPCPP_VERSION} CACHE INTERNAL "dpc++ version")
set(DPCPP_BETA_VERSION ${CMAKE_MATCH_1} CACHE INTERNAL "dpc++ beta version")
message(STATUS "dpc++ version: ${DPCPP_VERSION}")

# We do not support anything lower than beta09
if (DPCPP_BETA_VERSION LESS "09")
   message(FATAL_ERROR
      "\nUnsupported dpc++ compiler version."
      "\nAMReX requires dpc++ \"beta\" version >= 08. "
      "The current compiler \"beta\" version is ${DPCPP_BETA_VERSION}.\n")
endif ()


set(_cxx_dpcpp "$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:IntelClang>,$<CXX_COMPILER_ID:IntelDPCPP>>")
set(_cxx_dpcpp "$<AND:$<COMPILE_LANGUAGE:CXX>,${_cxx_dpcpp}>")

#
# SYCL and AMReX::SYCL targets
#
add_library(SYCL INTERFACE)
add_library(AMReX::SYCL ALIAS SYCL)

target_compile_features(SYCL INTERFACE cxx_std_17)


#
# Compiler options
#
target_compile_options( SYCL
   INTERFACE
   $<${_cxx_dpcpp}:-Wno-error=sycl-strict -Wno-pass-failed -fsycl>
   $<${_cxx_dpcpp}:$<$<BOOL:${AMReX_DPCPP_SPLIT_KERNEL}>:-fsycl-device-code-split=per_kernel>>)

# temporary work-around for DPC++ beta08 bug
#   define "long double" as 64bit for C++ user-defined literals
#   https://github.com/intel/llvm/issues/2187
target_compile_options( SYCL
   INTERFACE
   $<${_cxx_dpcpp}:-mlong-double-64 "SHELL:-Xclang -mlong-double-64">)

# Beta09 has enabled eary optimizations by default.  But this causes many
# tests to crash.  So we disable it.
target_compile_options( SYCL
   INTERFACE
   $<${_cxx_dpcpp}:-fno-sycl-early-optimizations>)

# Need this option to compile with mpiicpc
if (AMReX_MPI)
  target_compile_options( SYCL
    INTERFACE
    $<${_cxx_dpcpp}:-fsycl-unnamed-lambda>)
endif ()

#
# Link options
#
if (DPCPP_BETA_VERSION LESS "10")   # If beta < 10

   find_file(LIBSYCL_GLIBC_OBJ libsycl-glibc.o
      PATHS ${DPCPP_ROOT} ENV LD_LIBRARY_PATH
      PATH_SUFFIXES lib
      DOC "Full path to libsycl-glibc.o")

   target_link_libraries(SYCL INTERFACE ${LIBSYCL_GLIBC_OBJ})

   target_link_options( SYCL
      INTERFACE
      $<${_cxx_dpcpp}:-fsycl -device-math-lib=fp32,fp64> )

else ()  # for beta >= 10

   target_link_options( SYCL
      INTERFACE
      $<${_cxx_dpcpp}:-fsycl -fsycl-device-lib=libc,libm-fp32,libm-fp64> )

endif ()



# TODO: use $<LINK_LANG_AND_ID:> genex for CMake >=3.17
target_link_options( SYCL
   INTERFACE
   $<${_cxx_dpcpp}:$<$<BOOL:${AMReX_DPCPP_SPLIT_KERNEL}>:-fsycl-device-code-split=per_kernel>>)


if (AMReX_DPCPP_AOT)
   message(FATAL_ERROR "\nAhead-of-time (AOT) compilation support not available yet.\nRe-configure with AMReX_DPCPP_AOT=OFF.")

   #
   # TODO: remove comments to enable AOT support when the time comes
   #       (main blocker: missing math library)
   #
   # if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
   #    ## TODO: use file(READ)
   #    execute_process( COMMAND cat /sys/devices/cpu/caps/pmu_name OUTPUT_VARIABLE _cpu_long_name )
   # else ()
   #    message(FATAL_ERROR "\nAMReX_DPCPP_AOT is not supported on ${CMAKE_SYSTEM_NAME}\n")
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
   #       "Maybe set AMReX_DPCPP_AOT to OFF?\n")
   # endif ()

   # target_compile_options( amrex
   #    PUBLIC
   #    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:Clang>,$<NOT:$<CONFIG:Debug>>>:
   #    -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xsycl-target-backend "-device ${_cpu_short_name}" >)
   # target_link_options(amrex PUBLIC -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xsycl-target-backend "-device ${_cpu_short_name}" )
   # unset(_cpu_long_name)
   # unset(_cpu_short_name)
endif ()


unset(_cxx_dpcpp)
