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
   # Debug issues with -O0: internal compiler errors
   # work-around for
   #   https://github.com/AMReX-Codes/amrex/pull/3311
   foreach(D IN LISTS AMReX_SPACEDIM)
      target_compile_options(amrex_${D}d PUBLIC
         "$<$<AND:$<COMPILE_LANGUAGE:HIP>,$<CONFIG:Debug>>:-O1>"
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

   # ROCm relies on clang to offload code objects from (.a) archive files,
   # so we need to tell the offload-linker to include all code objects in archives.
   include(CheckLinkerFlag)
   check_linker_flag(
       HIP
       "SHELL:-Xoffload-linker --whole-archive"
       LINKER_HAS_WHOLE_ARCHIVE_OFFLOAD)
   if(LINKER_HAS_WHOLE_ARCHIVE_OFFLOAD)
       foreach(D IN LISTS AMReX_SPACEDIM)
           target_link_options(amrex_${D}d PUBLIC
               "$<$<LINK_LANGUAGE:HIP>:SHELL:-Xoffload-linker --whole-archive>")
       endforeach()
   endif()

   foreach(D IN LISTS AMReX_SPACEDIM)
       target_compile_options(amrex_${D}d PUBLIC $<$<COMPILE_LANGUAGE:HIP>:-m64>)

       # ROCm 4.5: use unsafe floating point atomics, otherwise atomicAdd is much slower
       # 
       target_compile_options(amrex_${D}d PUBLIC $<$<COMPILE_LANGUAGE:HIP>:-munsafe-fp-atomics>)
   endforeach()

   # Equivalently, relocatable-device-code (RDC) flags are needed for `extern`
   # device variable support (for codes that use global variables on device)
   # as well as our kernel fusion in AMReX, e.g. happening likely in amr regrid
   if(AMReX_GPU_RDC)
       foreach(D IN LISTS AMReX_SPACEDIM)
           target_compile_options(amrex_${D}d PUBLIC
              $<$<COMPILE_LANGUAGE:HIP>:-fgpu-rdc> )
           target_link_options(amrex_${D}d PUBLIC
              "$<$<LINK_LANGUAGE:HIP>:-fgpu-rdc>")
       endforeach()
   endif()
endif ()
