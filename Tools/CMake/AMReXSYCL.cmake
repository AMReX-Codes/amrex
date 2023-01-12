#
# This module defines the INTERFACE target SYCL and its alias AMReX::SYCL.
# These targets provides build/link requirements for the SYCL language.
# For the time being, only dpc++ is supported
#

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
   $<${_cxx_dpcpp}:-fsycl>
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
target_link_options( SYCL
   INTERFACE
   $<${_cxx_dpcpp}:-fsycl -fsycl-device-lib=libc,libm-fp32,libm-fp64> )


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
