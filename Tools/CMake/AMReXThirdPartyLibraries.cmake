#
# HDF5 -- here it would be best to create an imported target
#
if (AMReX_HDF5)
    if (AMReX_MPI)
       set(HDF5_PREFER_PARALLEL TRUE)
    endif ()
    find_package(HDF5 1.10.4 REQUIRED)

    if (AMReX_MPI AND (NOT HDF5_IS_PARALLEL))
       if (CMAKE_VERSION VERSION_LESS 3.27)
	      # The detection in earlier versions of cmake may not be reliable.
          # So we will try to do it ourselves. Work-around for:
          # https://gitlab.kitware.com/cmake/cmake/-/merge_requests/8234
          execute_process(
             COMMAND ${HDF5_C_COMPILER_EXECUTABLE} -showconfig
             OUTPUT_VARIABLE amrex_hdf5_config_output
             ERROR_VARIABLE amrex_hdf5_config_output
             OUTPUT_STRIP_TRAILING_WHITESPACE
             )
          if (amrex_hdf5_config_output MATCHES "Parallel HDF5: ([A-Za-z0-9]+)")
             if (${CMAKE_MATCH_1})
                set(HDF5_IS_PARALLEL TRUE)
             endif ()
          endif()
          unset(amrex_hdf5_config_output)
       endif ()
       if (NOT HDF5_IS_PARALLEL)
          message(FATAL_ERROR "\nHDF5 library does not support parallel I/O")
       endif ()
    endif ()

    if (HDF5_IS_PARALLEL AND (NOT AMReX_MPI))
       message(FATAL_ERROR "\nMPI enabled in HDF5 but not in AMReX, which will likely fail to build")
    endif ()

    if (TARGET hdf5::hdf5)  # CMake >= 3.19
       target_link_libraries(amrex PUBLIC hdf5::hdf5)
    else ()  # CMake < 3.19 -- Remove when minimum cmake version is bumped up
       target_include_directories(amrex PUBLIC ${HDF5_INCLUDE_DIRS})
       target_compile_definitions(amrex PUBLIC ${HDF5_DEFINITIONS})
       target_link_libraries(amrex PUBLIC ${HDF5_LIBRARIES})
    endif ()

endif ()

#
# H5Z-ZFP
#
if (AMReX_HDF5_ZFP)
   set(H5Z_ZFP_USE_STATIC_LIBS ON) # Static ON means using as a library, or OFF as an HDF5 plugin
   find_package(H5Z_ZFP 1.0.1 CONFIG)
   if (NOT AMReX_HDF5)
      message(FATAL_ERROR "\nHDF5 must be enabled for ZFP support in HDF5")
   endif ()

   if (TARGET h5z_zfp::h5z_zfp)  # CMake >= 3.19
      target_link_libraries(amrex PUBLIC h5z_zfp::h5z_zfp)
   else ()  # CMake < 3.19 -- Remove when minimum cmake version is bumped up
      target_include_directories(amrex PUBLIC ${H5Z_ZFP_INCLUDE_DIR})
      target_link_libraries(amrex PUBLIC ${H5Z_ZFP_LIBRARY})
   endif ()
endif ()

#
# Sensei
#
if (AMReX_SENSEI)
    find_package( SENSEI 4.0.0 REQUIRED )
    target_link_libraries( amrex PUBLIC sensei )
endif ()

#
#  Ascent
#
if (AMReX_ASCENT) # Ascent will find conduit, so check for Ascent first
    find_package(Ascent REQUIRED)
    if (AMReX_MPI)
        target_link_libraries( amrex PUBLIC ascent::ascent_mpi )
    else ()
        target_link_libraries( amrex PUBLIC ascent::ascent )
    endif ()
endif ()


#
# Conduit
#
if (AMReX_CONDUIT)
    find_package(Conduit REQUIRED)
    if (AMReX_MPI)
        target_link_libraries( amrex PUBLIC conduit::conduit_mpi )
    else ()
        target_link_libraries( amrex PUBLIC conduit::conduit )
    endif ()
endif ()


#
# HYPRE
#
if (AMReX_HYPRE)
    find_package(HYPRE 2.20.0 REQUIRED)
    if(AMReX_CUDA)
        find_package(CUDAToolkit REQUIRED)

        # mandatory CUDA dependencies: cuSPARSE, cuRAND
        target_link_libraries(amrex PUBLIC CUDA::cusparse CUDA::curand)
    endif()
    target_link_libraries( amrex PUBLIC HYPRE )
endif ()


#
# PETSc
#
if (AMReX_PETSC)
    find_package(PETSc 2.13 REQUIRED)
    target_link_libraries( amrex PUBLIC PETSC )
endif ()

#
# SUNDIALS
#
if (AMReX_SUNDIALS)
    if (SUNDIALS_FOUND)
        message(STATUS "SUNDIALS_FOUND is true, assuming nvecserial or gpu-specific vector found for version 6.0.0 or higher")
    else ()
       set(SUNDIALS_MINIMUM_VERSION 6.0.0 CACHE INTERNAL "Minimum required SUNDIALS version")
       find_package(SUNDIALS ${SUNDIALS_MINIMUM_VERSION} CONFIG QUIET )
    endif ()
    if (AMReX_GPU_BACKEND STREQUAL "CUDA")
       target_link_libraries( amrex PUBLIC SUNDIALS::nveccuda)
    elseif (AMReX_GPU_BACKEND STREQUAL "HIP")
       target_link_libraries( amrex PUBLIC SUNDIALS::nvechip)
    elseif (AMReX_GPU_BACKEND STREQUAL "SYCL")
       target_link_libraries( amrex PUBLIC SUNDIALS::nvecsycl)
    else ()
       target_link_libraries( amrex PUBLIC SUNDIALS::nvecserial)
    endif ()
endif ()
