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
# For now this is a place holder.
# CUDA stuff will go here after we get rid of AMReXSetupCUDA
#

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

   set(_valid_hip_compilers hipcc nvcc)
   get_filename_component(_this_comp ${CMAKE_CXX_COMPILER} NAME)

   if (NOT (_this_comp IN_LIST _valid_hip_compilers) )
      message(FATAL_ERROR "\nCMAKE_CXX_COMPILER is incompatible with HIP.\n"
         "Set CMAKE_CXX_COMPILER to either hipcc or nvcc for HIP builds.\n")
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

   find_package(HIP REQUIRED)

   if("${HIP_COMPILER}" STREQUAL "hcc")
      message(FATAL_ERROR "Using (deprecated) HCC compiler: please update ROCm")
   endif()

   if(HIP_FOUND)
      message(STATUS "Found HIP: ${HIP_VERSION}")
      message(STATUS "HIP: Platform=${HIP_PLATFORM} Compiler=${HIP_COMPILER}")
   else()
      message(FATAL_ERROR "Could not find HIP."
         " Ensure that HIP is either installed in /opt/rocm/hip or the variable HIP_PATH is set to point to the right location.")
   endif()

   # Link to hiprand -- must include rocrand too
   find_package(rocrand REQUIRED CONFIG)
   find_package(hiprand REQUIRED CONFIG)
   target_link_libraries(amrex PUBLIC hip::hiprand roc::rocrand)

   # ARCH flags -- these must be PUBLIC for all downstream targets to use,
   # else there will be a runtime issue (cannot find
   # missing gpu devices)
   target_compile_options(amrex
      PUBLIC
      $<$<COMPILE_LANGUAGE:CXX>:-m64 --amdgpu-target=${AMReX_AMD_ARCH}> )

endif ()
