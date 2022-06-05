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


set(_cxx_dpcpp "$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:IntelClang>,$<CXX_COMPILER_ID:IntelDPCPP>,$<CXX_COMPILER_ID:IntelLLVM>>")
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
     "$<${_cxx_dpcpp}:-mlong-double-64>"
     "$<${_cxx_dpcpp}:SHELL:-Xclang -mlong-double-64>")

# disable warning: comparison with infinity always evaluates to false in fast floating point modes [-Wtautological-constant-compare]
#                  return std::isinf(m);
# appeared since 2021.4.0
target_compile_options( SYCL
   INTERFACE
   $<${_cxx_dpcpp}:-Wno-tautological-constant-compare>)

# Need this option to compile with mpiicpc
if (AMReX_MPI)
  target_compile_options( SYCL
    INTERFACE
    $<${_cxx_dpcpp}:-fsycl-unnamed-lambda>)
endif ()

if(AMReX_DPCPP_ONEDPL)
    # TBB and PSTL are broken in oneAPI 2021.3.0
    # https://software.intel.com/content/www/us/en/develop/articles/intel-oneapi-dpcpp-library-release-notes.html#inpage-nav-2-3
    # at least since 2021.1.1 and probably won't be fixed until glibc version 10 is gone
    target_compile_definitions( SYCL
        INTERFACE
        $<${_cxx_dpcpp}:_GLIBCXX_USE_TBB_PAR_BACKEND=0 PSTL_USE_PARALLEL_POLICIES=0>)
endif()

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
   target_compile_options( SYCL
      INTERFACE
      "$<${_cxx_dpcpp}:-fsycl-targets=spir64_gen>"
      "$<${_cxx_dpcpp}:SHELL:-Xsycl-target-backend \"-device ${AMReX_INTEL_ARCH}\">" )

   target_link_options( SYCL
      INTERFACE
      "$<${_cxx_dpcpp}:-fsycl-targets=spir64_gen>"
      "$<${_cxx_dpcpp}:SHELL:-Xsycl-target-backend \"-device ${AMReX_INTEL_ARCH}\">" )
endif ()


unset(_cxx_dpcpp)
