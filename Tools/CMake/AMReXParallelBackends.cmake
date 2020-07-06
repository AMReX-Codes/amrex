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
