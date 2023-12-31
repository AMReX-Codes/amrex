#
# This module defines the INTERFACE target SYCL and its alias AMReX::SYCL.
# These targets provides build/link requirements for the SYCL language.
# For the time being, only oneAPI is supported
#

set(_cxx_sycl "$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:IntelClang>,$<CXX_COMPILER_ID:IntelDPCPP>,$<CXX_COMPILER_ID:IntelLLVM>>")
set(_cxx_sycl "$<AND:$<COMPILE_LANGUAGE:CXX>,${_cxx_sycl}>")

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
   $<${_cxx_sycl}:-fsycl>
   $<${_cxx_sycl}:$<$<BOOL:${AMReX_SYCL_SPLIT_KERNEL}>:-fsycl-device-code-split=per_kernel>>)

# temporary work-around for oneAPI beta08 bug
#   define "long double" as 64bit for C++ user-defined literals
#   https://github.com/intel/llvm/issues/2187
target_compile_options( SYCL
   INTERFACE
     "$<${_cxx_sycl}:-mlong-double-64>"
     "$<${_cxx_sycl}:SHELL:-Xclang -mlong-double-64>")

# disable warning: comparison with infinity always evaluates to false in fast floating point modes [-Wtautological-constant-compare]
#                  return std::isinf(m);
# appeared since 2021.4.0
target_compile_options( SYCL
   INTERFACE
   $<${_cxx_sycl}:-Wno-tautological-constant-compare>)

if(AMReX_SYCL_ONEDPL)
    # TBB and PSTL are broken in oneAPI 2021.3.0
    # https://software.intel.com/content/www/us/en/develop/articles/intel-oneapi-dpcpp-library-release-notes.html#inpage-nav-2-3
    # at least since 2021.1.1 and probably won't be fixed until glibc version 10 is gone
    target_compile_definitions( SYCL
        INTERFACE
        $<${_cxx_sycl}:_GLIBCXX_USE_TBB_PAR_BACKEND=0 PSTL_USE_PARALLEL_POLICIES=0>)
endif()

#
# Link options
#
target_link_options( SYCL
   INTERFACE
   $<${_cxx_sycl}:-qmkl=sequential -fsycl -fsycl-device-lib=libc,libm-fp32,libm-fp64> )


# TODO: use $<LINK_LANG_AND_ID:> genex for CMake >=3.17
target_link_options( SYCL
   INTERFACE
   $<${_cxx_sycl}:$<$<BOOL:${AMReX_SYCL_SPLIT_KERNEL}>:-fsycl-device-code-split=per_kernel>>)


if (AMReX_SYCL_AOT)
   target_compile_options( SYCL
      INTERFACE
      "$<${_cxx_sycl}:-fsycl-targets=spir64_gen>" )

   set(_sycl_backend_flags "-device ${AMReX_INTEL_ARCH}")
   if (AMReX_SYCL_AOT_GRF_MODE STREQUAL "Large")
      set(_sycl_backend_flags "${_sycl_backend_flags} -internal_options -ze-opt-large-register-file")
   elseif (AMReX_SYCL_AOT_GRF_MODE STREQUAL "AutoLarge")
      set(_sycl_backend_flags "${_sycl_backend_flags} -options -ze-intel-enable-auto-large-GRF-mode")
   endif()

   target_link_options( SYCL
      INTERFACE
      "$<${_cxx_sycl}:-fsycl-targets=spir64_gen>"
      "$<${_cxx_sycl}:SHELL:-Xsycl-target-backend \"${_sycl_backend_flags}\">" )

   unset(_sycl_backend_flags)
endif ()

if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND "${CMAKE_BUILD_TYPE}" MATCHES "Debug")
   target_link_options( SYCL
      INTERFACE
      "$<${_cxx_sycl}:-fsycl-link-huge-device-code>" )
endif ()

if (AMReX_PARALLEL_LINK_JOBS GREATER 1)
   target_link_options( SYCL
      INTERFACE
      $<${_cxx_sycl}:-fsycl-max-parallel-link-jobs=${AMReX_PARALLEL_LINK_JOBS}>)
endif()

unset(_cxx_sycl)
