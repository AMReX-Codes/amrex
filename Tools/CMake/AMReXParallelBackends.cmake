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
if (ENABLE_MPI)
   set(_mpi_comps C CXX)  # Do we need MPI_C ?
   if (ENABLE_FORTRAN_INTERFACES)
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
if (ENABLE_OMP)
   set(_omp_comps CXX)
   if (ENABLE_FORTRAN_INTERFACES)
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
if (ENABLE_DPCPP)
   include(AMReXSYCL)
   target_link_libraries(amrex PUBLIC SYCL)
endif ()


#
#
# HIP
#
#
if (ENABLE_HIP)
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

   # Let's put the defines here for the time being
   target_compile_definitions( amrex PUBLIC AMREX_USE_HIP AMREX_HIP_PLATFORM=${HIP_PLATFORM} )

endif ()
