#
#
#  Pthreads -- always required
#
#
set( THREADS_PREFER_PTHREAD_FLAG on )
find_package( Threads REQUIRED )
target_link_libraries( amrex PUBLIC Threads::Threads )


#
#
#  MPI
#
#
if (AMReX_MPI)
   set(_mpi_comps C CXX)  # Do we need MPI_C ?
   if (AMReX_FORTRAN_INTERFACES)
      list(APPEND _mpi_comps Fortran)
   endif ()
   find_package(MPI REQUIRED ${_mpi_comps})
   list(TRANSFORM _mpi_comps PREPEND "MPI::MPI_")
   target_link_libraries(amrex PUBLIC ${_mpi_comps})
   unset(_mpi_comps)
endif ()

#
#
#  OpenMP
#
#
if (AMReX_OMP)
   set(_omp_comps CXX)
   if (AMReX_FORTRAN)
      list(APPEND _omp_comps Fortran)
   endif ()
   find_package(OpenMP REQUIRED ${_omp_comps})
   list(TRANSFORM _omp_comps PREPEND "OpenMP::OpenMP_")
   target_link_libraries(amrex PUBLIC ${_omp_comps})
   unset(_omp_comps)
else ()
   target_compile_options( amrex
      PUBLIC
      $<$<CXX_COMPILER_ID:Cray>:-h;noomp> )
endif ()

#
#
# CUDA
#
#
#
if (  AMReX_GPU_BACKEND STREQUAL "CUDA"
      AND
      CMAKE_VERSION VERSION_GREATER_EQUAL 3.20 )

   # Check cuda compiler and host compiler
   set_mininum_compiler_version(CUDA NVIDIA 9.0)
   check_cuda_host_compiler()

   # Required CUDA flags
   set(_genex "$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>")
   target_compile_options( amrex
      PUBLIC
      $<${_genex}:
      --expt-relaxed-constexpr --expt-extended-lambda
      "SHELL:-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored"
      -maxrregcount=${AMReX_CUDA_MAXREGCOUNT}
      "SHELL:-Xcudafe --display_error_number"
      $<$<STREQUAL:$<PLATFORM_ID>,Windows>:-m64> >
      )

   # Take care of cuda archs
   set_cuda_architectures(AMReX_CUDA_ARCH)
   set_target_properties( amrex
      PROPERTIES
      CUDA_ARCHITECTURES "${AMREX_CUDA_ARCHS}"
      )

   #
   # CUDA specific warnings
   #
   set(_cuda_flags)
   if (AMReX_CUDA_WARN_CAPTURE_THIS)
      list(APPEND _cuda_flags --Wext-lambda-captures-this)
   endif()
   if (AMReX_CUDA_ERROR_CAPTURE_THIS)
      # note: prefer double-dash --Werror!
      # https://github.com/ccache/ccache/issues/598
      list(APPEND _cuda_flags "SHELL:--Werror ext-lambda-captures-this")
   endif()
   if (AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL)
      list(APPEND _cuda_flags "SHELL:--Werror cross-execution-space-call")
   endif()

   #
   # Forward unknown NVCC flags to the host compiler
   #
   if (CUDA_FORWARD_UNKNOWN_FLAGS_HOST)
      list(APPEND _cuda_flags --forward-unknown-to-host-compiler)
   endif()

   # fast math
   if (AMReX_CUDA_FASTMATH)
      list(APPEND _cuda_flags --use_fast_math)
   endif ()

   #
   # Code generation
   #
   if (AMReX_CUDA_PTX_VERBOSE)
      list(APPEND _cuda_flags --ptxas-options=-v)
   endif()

   # keep intermediately generated files
   if (AMReX_CUDA_KEEP_FILES)
      make_directory("${PROJECT_BINARY_DIR}/nvcc_tmp")
      list(APPEND _cuda_flags --keep "SHELL:--keep-dir ${PROJECT_BINARY_DIR}/nvcc_tmp")
   endif ()

   # compilation timings
   if (AMReX_CUDA_COMPILATION_TIMER)
      file(REMOVE "${PROJECT_BINARY_DIR}/nvcc_timings.csv")
      list(APPEND _cuda_flags "SHELL:--time ${PROJECT_BINARY_DIR}/nvcc_timings.csv")
   endif ()

   #
   # Debugging
   #
   if (AMReX_CUDA_DEBUG)
      # is this unsupported with MSVC?
      list(APPEND _cuda_flags -G)
   endif()

   if (AMReX_CUDA_SHOW_LINENUMBERS AND NOT AMReX_CUDA_DEBUG)
      # nvcc warning : '--device-debug (-G)' overrides '--generate-line-info (-lineinfo)'
      list(APPEND _cuda_flags --generate-line-info)
   endif ()
   if (AMReX_CUDA_SHOW_CODELINES)
      list(APPEND _cuda_flags --source-in-ptx)
   endif ()

   if (AMReX_CUDA_BACKTRACE)
      if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
         list(APPEND _cuda_flags "SHELL:-Xcompiler /Zi") # comes with Debug & RelWithDebInfo
      else ()
         list(APPEND _cuda_flags "SHELL:Xcompiler -rdynamic")
      endif ()
   endif ()

   # Flags to make it an error to write a device variable in
   # a host function.
   if (CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.2)
      list(APPEND _cuda_flag --display-error-number "SHELL:--diag-error 20092")
   endif ()

   target_compile_options( amrex PUBLIC $<${_genex}:${_cuda_flags}> )

endif ()

#
#
#  SYCL/DPCPP
#
#
if (AMReX_DPCPP)
   include(AMReXSYCL)
   target_link_libraries(amrex PUBLIC SYCL)
endif ()


#
#
# HIP
#
#
if (AMReX_HIP)

   set(_valid_hip_compilers clang++ hipcc nvcc CC)
   get_filename_component(_this_comp ${CMAKE_CXX_COMPILER} NAME)

   if (NOT (_this_comp IN_LIST _valid_hip_compilers) )
      message(WARNING "\nCMAKE_CXX_COMPILER (${_this_comp}) is likely "
         "incompatible with HIP.\n"
         "Set CMAKE_CXX_COMPILER to either AMD's clang++ (preferred) or "
         "hipcc or nvcc for HIP builds.\n")
   endif ()

   unset(_hip_compiler)
   unset(_valid_hip_compilers)

   if(NOT DEFINED HIP_PATH)
      if(NOT DEFINED ENV{HIP_PATH})
         set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
      else()
         set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
      endif()
   endif()

   set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})


   if(DEFINED AMReX_AMD_ARCH)
      # Set the GPU to compile for: semicolon-separated list
      set(GPU_TARGETS "${AMReX_AMD_ARCH}" CACHE STRING "GPU targets to compile for" FORCE)
      set(AMDGPU_TARGETS "${AMReX_AMD_ARCH}" CACHE STRING "GPU targets to compile for" FORCE)
      mark_as_advanced(AMDGPU_TARGETS)
      mark_as_advanced(GPU_TARGETS)
   endif()

   find_package(hip)

   if("${HIP_COMPILER}" STREQUAL "hcc")
      message(FATAL_ERROR "Using (deprecated) HCC compiler: please update ROCm")
   endif()

   if(hip_FOUND)
      message(STATUS "Found HIP: ${HIP_VERSION}")
      message(STATUS "HIP: Runtime=${HIP_RUNTIME} Compiler=${HIP_COMPILER} Path=${HIP_PATH}")
   else()
      message(FATAL_ERROR "Could not find HIP."
         " Ensure that HIP is either installed in /opt/rocm/hip or the variable HIP_PATH is set to point to the right location.")
   endif()

   if(${_this_comp} STREQUAL hipcc AND NOT AMReX_FORTRAN)
       message(WARNING "You are using the legacy wrapper 'hipcc' as the HIP compiler.\n"
           "This is only needed when building with Fortran support and with ROCm/HIP <=4.2.0. "
           "Use AMD's 'clang++' compiler instead.")
   endif()
   # AMD's or mainline clang++ with support for "-x hip"
   # Cray's CC wrapper that points to AMD's clang++ underneath
   if(NOT ${_this_comp} STREQUAL hipcc)
       target_link_libraries(amrex PUBLIC hip::device)

       # work-around for https://github.com/ROCm-Developer-Tools/HIP/issues/2278
       # CXX_STANDARD always adds -std=c++XX, even if the compiler default fulfills it
       #set_property(TARGET amrex PROPERTY CXX_STANDARD 17)
       # note: already bumped to C++17 (cxx_std_17) or newer in AMReX_Config.cmake

       # work-around for ROCm <=4.2
       # https://github.com/ROCm-Developer-Tools/HIP/pull/2190
       target_compile_options(amrex PUBLIC
          "$<$<COMPILE_LANGUAGE:CXX>:SHELL:-mllvm;-amdgpu-early-inline-all=true;-mllvm;-amdgpu-function-calls=false>"
       )
       target_compile_options(amrex PUBLIC
          "$<$<COMPILE_LANGUAGE:CXX>:SHELL:-x hip>"
       )
   endif()

   # Link to hiprand -- must include rocrand too
   find_package(rocrand REQUIRED CONFIG)
   find_package(rocprim REQUIRED CONFIG)
   find_package(hiprand REQUIRED CONFIG)
   if(AMReX_ROCTX)
       # To be modernized in the future, please see:
       # https://github.com/ROCm-Developer-Tools/roctracer/issues/56
       target_include_directories(amrex PUBLIC ${HIP_PATH}/../roctracer/include ${HIP_PATH}/../rocprofiler/include)
       target_link_libraries(amrex PUBLIC "-L${HIP_PATH}/../roctracer/lib/ -lroctracer64" "-L${HIP_PATH}/../roctracer/lib -lroctx64")
   endif ()
   target_link_libraries(amrex PUBLIC hip::hiprand roc::rocrand roc::rocprim)

   # avoid forcing the rocm LLVM flags on a gfortran
   # https://github.com/ROCm-Developer-Tools/HIP/issues/2275
   if(AMReX_FORTRAN)
       message(WARNING "As of ROCm/HIP <= 4.2.0, Fortran support might be flaky.\n"
                       "Especially, we cannot yet support reloctable device code (RDC)."
                       "See https://github.com/ROCm-Developer-Tools/HIP/issues/2275 "
                       "and https://github.com/AMReX-Codes/amrex/pull/2031 "
                       "for details.")
   elseif(${_this_comp} STREQUAL hipcc)
       # hipcc expects a comma-separated list
       string(REPLACE ";" "," AMReX_AMD_ARCH_HIPCC "${AMReX_AMD_ARCH}")

       target_link_libraries(amrex PUBLIC ${HIP_LIBRARIES})
       # ARCH flags -- these must be PUBLIC for all downstream targets to use,
       # else there will be a runtime issue (cannot find
       # missing gpu devices)
       target_compile_options(amrex PUBLIC
          $<$<COMPILE_LANGUAGE:CXX>:--amdgpu-target=${AMReX_AMD_ARCH_HIPCC} -Wno-pass-failed>)
   endif()

   target_compile_options(amrex PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-m64>)

   # Equivalently, relocatable-device-code (RDC) flags are needed for `extern`
   # device variable support (for codes that use global variables on device)
   # as well as our kernel fusion in AMReX, e.g. happening likely in amr regrid
   # As of ROCm 4.1, we cannot enable this with hipcc, as it looks...
   if(AMReX_GPU_RDC)
       target_compile_options(amrex PUBLIC
          $<$<COMPILE_LANGUAGE:CXX>:-fgpu-rdc> )
       if(CMAKE_VERSION VERSION_LESS 3.18)
           target_link_options(amrex PUBLIC
              -fgpu-rdc)
       else()
           target_link_options(amrex PUBLIC
              "$<$<LINK_LANGUAGE:CXX>:-fgpu-rdc>")
       endif()
   endif()

endif ()
