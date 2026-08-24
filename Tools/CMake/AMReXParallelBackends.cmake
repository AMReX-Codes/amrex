#
#
#  Pthreads -- always required
#
#
set( THREADS_PREFER_PTHREAD_FLAG on )
find_package( Threads REQUIRED )
foreach(D IN LISTS AMReX_SPACEDIM)
    target_link_libraries( amrex_${D}d PUBLIC Threads::Threads )
endforeach()


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
   foreach(D IN LISTS AMReX_SPACEDIM)
       target_link_libraries(amrex_${D}d PUBLIC ${_mpi_comps})
   endforeach()
   unset(_mpi_comps)
endif ()

#
#
#  VIR-SIMD
#
#
if (AMReX_SIMD)
   find_package(vir-simd REQUIRED)
   foreach(D IN LISTS AMReX_SPACEDIM)
       target_link_libraries(amrex_${D}d PUBLIC vir-simd::vir-simd)
   endforeach()
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
   foreach(D IN LISTS AMReX_SPACEDIM)
       target_link_libraries(amrex_${D}d PUBLIC ${_omp_comps})
   endforeach()
   unset(_omp_comps)
else ()
   foreach(D IN LISTS AMReX_SPACEDIM)
       target_compile_options(amrex_${D}d
          PUBLIC
          $<$<CXX_COMPILER_ID:Cray>:-h;noomp> )
   endforeach()
endif ()

#
#
# CUDA
#
#
#
if ( AMReX_GPU_BACKEND STREQUAL "CUDA" )

   find_package(CUDAToolkit REQUIRED)
   foreach(D IN LISTS AMReX_SPACEDIM)
       target_link_libraries(amrex_${D}d PUBLIC CUDA::curand)

       if (AMReX_LINEAR_SOLVERS)
           target_link_libraries(amrex_${D}d PUBLIC CUDA::cusparse)
       endif ()
   endforeach()

   # Check cuda compiler and host compiler
   set_mininum_compiler_version(CUDA NVIDIA 12.2)
   check_cuda_host_compiler()

   # Required CUDA flags
   set(_genex "$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>")
   foreach(D IN LISTS AMReX_SPACEDIM)
       target_compile_options(amrex_${D}d
          PUBLIC
          $<${_genex}:
          --expt-relaxed-constexpr --expt-extended-lambda
          "SHELL:-Xcudafe --diag_suppress=esa_on_defaulted_function_ignored"
          "SHELL:-Xcudafe --diag_suppress=implicit_return_from_non_void_function"
          -maxrregcount=${AMReX_CUDA_MAXREGCOUNT}
          "SHELL:-Xcudafe --display_error_number"
          $<$<STREQUAL:$<PLATFORM_ID>,Windows>:-m64> >
          )
  endforeach()

   foreach(D IN LISTS AMReX_SPACEDIM)
       set_target_properties(amrex_${D}d
          PROPERTIES
          CUDA_ARCHITECTURES "${AMREX_CUDA_ARCHS}"
       )
       if (AMREX_CUDA_IPO)
          # relocatable device code is required for device LTO and enforced in
          # AMReXCUDAOptions; CUDA_SEPARABLE_COMPILATION itself is set from AMReX_GPU_RDC
          # in setup_target_for_cuda_compilation
          set_target_properties(amrex_${D}d
             PROPERTIES
             INTERPROCEDURAL_OPTIMIZATION ON
          )
          # For a static AMReX, the CUDA device link happens in the user's target,
          # so export the device-LTO flag as an interface requirement. -dlto is
          # nvcc's flag, currently `clang -x cuda` or NVHPC do not take it
          # and CMake check_ipo_supported reports device LTO unsupported for them.
          # A consumer that enables INTERPROCEDURAL_OPTIMIZATION itself - every target
          # that goes through setup_target_for_cuda_compilation does - already gets -dlto
          # from CMake, so only add it for the ones that do not.
          target_link_options(amrex_${D}d INTERFACE
             "$<DEVICE_LINK:$<$<AND:$<CUDA_COMPILER_ID:NVIDIA>,$<NOT:$<BOOL:$<TARGET_PROPERTY:INTERPROCEDURAL_OPTIMIZATION>>>>:-dlto>>")
       endif ()
   endforeach()

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

   # place intermediate files in object file folder
   if (AMReX_CUDA_OBJDIR_AS_TEMPDIR)
      list(APPEND _cuda_flags --objdir-as-tempdir)
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
         list(APPEND _cuda_flags "SHELL:-Xcompiler -rdynamic")
      endif ()
   endif ()

   # Flags to make it an error to write a device variable in
   # a host function.
   if (CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.2)
      list(APPEND _cuda_flags --display-error-number "SHELL:--diag-error 20092")
   endif ()

   foreach(D IN LISTS AMReX_SPACEDIM)
      target_compile_options(amrex_${D}d PUBLIC $<${_genex}:${_cuda_flags}> )
   endforeach()

   unset(_genex)
   # _cuda_flags will be used later in AMReX_Config.cmake
endif ()

#
#
#  SYCL
#
#
if (AMReX_SYCL)
   include(AMReXSYCL)
   foreach(D IN LISTS AMReX_SPACEDIM)
      target_link_libraries(amrex_${D}d PUBLIC SYCL)
   endforeach()
endif ()


#
#
# HIP
#
#
if (AMReX_HIP)

   set(_valid_hip_compilers clang++ amdclang++ hipcc nvcc CC)
   get_filename_component(_this_comp ${CMAKE_CXX_COMPILER} NAME)

   if (NOT (_this_comp IN_LIST _valid_hip_compilers) )
      message(WARNING "\nCMAKE_CXX_COMPILER (${_this_comp}) is likely "
         "incompatible with HIP.\n"
         "Set CMAKE_CXX_COMPILER to either Cray's CC or AMD's clang++/amdclang++ "
         "or hipcc or nvcc for HIP builds.\n")
   endif ()

   unset(_hip_compiler)
   unset(_valid_hip_compilers)

   if(NOT DEFINED HIP_PATH)
      if(DEFINED ENV{HIP_PATH})
         set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
      elseif(DEFINED ENV{ROCM_PATH})
         set(HIP_PATH "$ENV{ROCM_PATH}/hip" CACHE PATH "Path to which HIP has been installed")
      else()
         set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
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

   # AMD's or mainline clang++ with support for "-x hip"
   # Cray's CC wrapper that points to AMD's clang++ underneath
   if(NOT ${_this_comp} STREQUAL hipcc)
       foreach(D IN LISTS AMReX_SPACEDIM)
          target_link_libraries(amrex_${D}d PUBLIC hip::device)

          # work-around for
          #   https://github.com/ROCm-Developer-Tools/hipamd/issues/12
          # not being added for Cray CC
          target_compile_options(amrex_${D}d PUBLIC
             "$<$<COMPILE_LANGUAGE:CXX>:SHELL:-mllvm;-amdgpu-early-inline-all=true;-mllvm;-amdgpu-function-calls=false>"
          )
       endforeach()
   endif()

   # Debug issues with -O0: internal compiler errors
   # work-around for
   #   https://github.com/AMReX-Codes/amrex/pull/3311
   foreach(D IN LISTS AMReX_SPACEDIM)
      target_compile_options(amrex_${D}d PUBLIC
         "$<$<CONFIG:Debug>:-O1>"
      )
   endforeach()

   # Link to hiprand -- must include rocrand too
   find_package(rocrand REQUIRED CONFIG)
   find_package(rocprim REQUIRED CONFIG)
   find_package(hiprand REQUIRED CONFIG)
   if (AMReX_LINEAR_SOLVERS)
      find_package(rocsparse REQUIRED CONFIG)
   endif()

   if(AMReX_ROCTX)
       find_package(rocprofiler-sdk-roctx REQUIRED CONFIG)
       foreach(D IN LISTS AMReX_SPACEDIM)
          target_link_libraries(amrex_${D}d PUBLIC rocprofiler-sdk-roctx::rocprofiler-sdk-roctx)
      endforeach()
   endif()
   foreach(D IN LISTS AMReX_SPACEDIM)
      target_link_libraries(amrex_${D}d PUBLIC hip::hiprand roc::rocrand roc::rocprim)
   endforeach()
   if (AMReX_LINEAR_SOLVERS)
      foreach(D IN LISTS AMReX_SPACEDIM)
         target_link_libraries(amrex_${D}d PUBLIC roc::rocsparse)
      endforeach()
   endif()

   # avoid forcing the rocm LLVM flags on a gfortran
   # https://github.com/ROCm-Developer-Tools/HIP/issues/2275
   if(${_this_comp} STREQUAL hipcc)
       # hipcc expects a comma-separated list
       string(REPLACE ";" "," AMReX_AMD_ARCH_HIPCC "${AMReX_AMD_ARCH}")

       foreach(D IN LISTS AMReX_SPACEDIM)
           target_link_libraries(amrex_${D}d PUBLIC ${HIP_LIBRARIES})
           # ARCH flags -- these must be PUBLIC for all downstream targets to use,
           # else there will be a runtime issue (cannot find
           # missing gpu devices)
           target_compile_options(amrex_${D}d PUBLIC
              $<$<COMPILE_LANGUAGE:CXX>:--offload-arch=${AMReX_AMD_ARCH_HIPCC}>)
       endforeach()
   endif()

   # ROCm 5.5: hipcc now relies on clang to offload code objects from (.a) archive files,
   # so we need to tell the offload-linker to include all code objects in archives.
   include(CheckLinkerFlag)
   check_linker_flag(
       CXX
       "SHELL:-Xoffload-linker --whole-archive"
       LINKER_HAS_WHOLE_ARCHIVE_OFFLOAD)
   if(LINKER_HAS_WHOLE_ARCHIVE_OFFLOAD)
       foreach(D IN LISTS AMReX_SPACEDIM)
           target_link_options(amrex_${D}d PUBLIC
               "$<$<LINK_LANGUAGE:HIP>:SHELL:-Xoffload-linker --whole-archive>"
               "$<$<LINK_LANGUAGE:CXX>:SHELL:-Xoffload-linker --whole-archive>")
       endforeach()
   endif()

   foreach(D IN LISTS AMReX_SPACEDIM)
       target_compile_options(amrex_${D}d PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-m64>)

       # ROCm 4.5: use unsafe floating point atomics, otherwise atomicAdd is much slower
       # 
       target_compile_options(amrex_${D}d PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-munsafe-fp-atomics>)

       # Ensure ROCm builds enable at least C++20 without overriding higher standards
       target_compile_features(amrex_${D}d PUBLIC cxx_std_20)
   endforeach()

   # Equivalently, relocatable-device-code (RDC) flags are needed for `extern`
   # device variable support (for codes that use global variables on device)
   # as well as our kernel fusion in AMReX, e.g. happening likely in amr regrid
   # As of ROCm 4.1, we cannot enable this with hipcc, as it looks...
   if(AMReX_GPU_RDC)
       foreach(D IN LISTS AMReX_SPACEDIM)
           target_compile_options(amrex_${D}d PUBLIC
              $<$<COMPILE_LANGUAGE:CXX>:-fgpu-rdc> )
           target_link_options(amrex_${D}d PUBLIC
              "$<$<LINK_LANGUAGE:HIP>:-fgpu-rdc>"
              "$<$<LINK_LANGUAGE:CXX>:-fgpu-rdc>")
       endforeach()
   endif()
endif ()
