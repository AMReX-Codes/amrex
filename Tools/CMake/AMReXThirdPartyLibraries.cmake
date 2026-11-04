#
# CCCL: provides cuda::fp32mp2 (<cuda/fpmp>) for AMReX_PRECISION=FP32MP2
#
if (AMReX_PRECISION STREQUAL "FP32MP2")
    # In CUDA builds, CUB and Thrust must come from the same CCCL release as
    # libcudacxx, ahead of the (older) copy bundled with the CUDA toolkit.
    # CCCL's package config uses -I (not -isystem) so that it takes priority.
    if (AMReX_CUDA)
        set(_amrex_cccl_components libcudacxx CUB Thrust)
    else ()
        set(_amrex_cccl_components libcudacxx)
    endif ()
    find_package(CCCL 3.6 CONFIG REQUIRED COMPONENTS ${_amrex_cccl_components})
    unset(_amrex_cccl_components)

    include(CheckCXXSourceCompiles)
    get_target_property(_amrex_cccl_includes CCCL::libcudacxx INTERFACE_INCLUDE_DIRECTORIES)
    set(CMAKE_REQUIRED_INCLUDES ${_amrex_cccl_includes})
    # AMReX uses the implicit conversion mode of <cuda/fpmp> (e.g., Real x = 0.5;).
    # It changes declarations, so every translation unit needs the same value.
    set(CMAKE_REQUIRED_DEFINITIONS -DCCCL_FPMP_EXPLICIT_CASTS=0)
    unset(AMReX_CCCL_HAS_FPMP CACHE)  # re-check if CCCL_DIR changed
    check_cxx_source_compiles("
        #include <cuda/fpmp>
        int main () { cuda::fp32mp2 x = 1.0; x += x; return static_cast<int>(x) - 2; }"
        AMReX_CCCL_HAS_FPMP)
    unset(CMAKE_REQUIRED_DEFINITIONS)
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(_amrex_cccl_includes)
    if (NOT AMReX_CCCL_HAS_FPMP)
        message(FATAL_ERROR "AMReX_PRECISION=FP32MP2: CCCL ${CCCL_VERSION} found in "
            "${CCCL_DIR} does not provide <cuda/fpmp> with cuda::fp32mp2 and "
            "CCCL_FPMP_EXPLICIT_CASTS=0. "
            "Point CCCL_DIR to a CCCL >= 3.6 checkout (<cccl>/lib/cmake/cccl).")
    endif ()

    foreach(D IN LISTS AMReX_SPACEDIM)
        target_link_libraries(amrex_${D}d PUBLIC CCCL::CCCL)
        target_compile_definitions(amrex_${D}d PUBLIC CCCL_FPMP_EXPLICIT_CASTS=0)
    endforeach()
endif ()

#
# FFT
#
if (AMReX_FFT)
    if (AMReX_CUDA)
        find_package(CUDAToolkit REQUIRED)
        foreach(D IN LISTS AMReX_SPACEDIM)
            target_link_libraries(amrex_${D}d PUBLIC CUDA::cufft)
        endforeach()
    elseif (AMReX_HIP)
        find_package(rocfft REQUIRED)
        foreach(D IN LISTS AMReX_SPACEDIM)
            target_link_libraries(amrex_${D}d PUBLIC roc::rocfft)
        endforeach()
    elseif (AMReX_SYCL)
        # nothing to do
    else()
        find_package(AMReXFFTW REQUIRED)
        foreach(D IN LISTS AMReX_SPACEDIM)
            target_link_libraries(amrex_${D}d PUBLIC AMReX::FFTW)
        endforeach()
    endif()
endif()

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

    foreach(D IN LISTS AMReX_SPACEDIM)
        if (TARGET hdf5::hdf5)  # CMake >= 3.19
           target_link_libraries(amrex_${D}d PUBLIC hdf5::hdf5)
        else ()  # CMake < 3.19 -- Remove when minimum cmake version is bumped up
           target_include_directories(amrex_${D}d PUBLIC ${HDF5_INCLUDE_DIRS})
           target_compile_definitions(amrex_${D}d PUBLIC ${HDF5_DEFINITIONS})
           target_link_libraries(amrex_${D}d PUBLIC ${HDF5_LIBRARIES})
        endif ()
    endforeach()

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

   foreach(D IN LISTS AMReX_SPACEDIM)
       if (TARGET h5z_zfp::h5z_zfp)  # CMake >= 3.19
          target_link_libraries(amrex_${D}d PUBLIC h5z_zfp::h5z_zfp)
       else ()  # CMake < 3.19 -- Remove when minimum cmake version is bumped up
          target_include_directories(amrex_${D}d PUBLIC ${H5Z_ZFP_INCLUDE_DIR})
          target_link_libraries(amrex_${D}d PUBLIC ${H5Z_ZFP_LIBRARY})
       endif ()
   endforeach()
endif ()

#
# Sensei
#
if (AMReX_SENSEI)
    find_package( SENSEI 4.0.0 REQUIRED )
    foreach(D IN LISTS AMReX_SPACEDIM)
        target_link_libraries(amrex_${D}d PUBLIC sensei)
    endforeach()
endif ()

#
#  Ascent
#
if (AMReX_ASCENT) # Ascent will find conduit, so check for Ascent first
    find_package(Ascent REQUIRED)
    foreach(D IN LISTS AMReX_SPACEDIM)
        if (AMReX_MPI)
            target_link_libraries(amrex_${D}d PUBLIC ascent::ascent_mpi)
        else ()
            target_link_libraries(amrex_${D}d PUBLIC ascent::ascent)
        endif ()
    endforeach()
endif ()

#
#  Catalyst
#
if (AMReX_CATALYST)
    find_package(Catalyst REQUIRED PATHS "$ENV{CATALYST_IMPLEMENTATION_PATHS}")
    foreach(D IN LISTS AMReX_SPACEDIM)
        target_link_libraries(amrex_${D}d PUBLIC catalyst::catalyst)
    endforeach()
endif ()

#
# Conduit
#
if (AMReX_CONDUIT)
    if(NOT TARGET conduit::conduit)
        find_package(Conduit REQUIRED)
    endif()
    foreach(D IN LISTS AMReX_SPACEDIM)
        if (AMReX_MPI)
            target_link_libraries(amrex_${D}d PUBLIC conduit::conduit_mpi)
        else ()
            target_link_libraries(amrex_${D}d PUBLIC conduit::conduit)
        endif ()
    endforeach()
endif ()


#
# HYPRE
#
if (AMReX_HYPRE)
    find_package(HYPRE 2.20.0 REQUIRED)
    if(AMReX_CUDA)
        find_package(CUDAToolkit REQUIRED)

        # mandatory CUDA dependencies: cuSPARSE, cuRAND
        foreach(D IN LISTS AMReX_SPACEDIM)
            target_link_libraries(amrex_${D}d PUBLIC CUDA::cusparse CUDA::curand)
        endforeach()
    endif()
    foreach(D IN LISTS AMReX_SPACEDIM)
        target_link_libraries(amrex_${D}d PUBLIC HYPRE)
    endforeach()
endif ()


#
# PETSc
#
if (AMReX_PETSC)
    find_package(PETSc 2.13 REQUIRED)
    foreach(D IN LISTS AMReX_SPACEDIM)
        target_link_libraries(amrex_${D}d PUBLIC PETSC)
    endforeach()
endif ()

#
# SUNDIALS
#
if (AMReX_SUNDIALS)
    if (SUNDIALS_FOUND)
        message(STATUS "SUNDIALS_FOUND is true, assuming nvecserial or gpu-specific vector found for version 6.0.0 or higher")
    else ()
       set(SUNDIALS_MINIMUM_VERSION 6.0.0 CACHE INTERNAL "Minimum required SUNDIALS version")
       set(SUNDIALS_COMPONENTS arkode cvode sunlinsolspgmr sunlinsolspfgmr
                               nvecserial nvecmanyvector sunnonlinsolfixedpoint)
       find_package(SUNDIALS CONFIG REQUIRED
                    COMPONENTS ${SUNDIALS_COMPONENTS}
                    OPTIONAL_COMPONENTS core) # core only available for >= 7
       if (SUNDIALS_VERSION VERSION_LESS ${SUNDIALS_MINIMUM_VERSION})
          message(FATAL_ERROR "SUNDIALS_VERSION ${SUNDIALS_MINIMUM_VERSION} or newer is required. Found version ${SUNDIALS_VERSION}.")
       endif ()
    endif ()
    foreach(D IN LISTS AMReX_SPACEDIM)
        if (SUNDIALS_VERSION VERSION_GREATER_EQUAL 7)
           target_link_libraries(amrex_${D}d PUBLIC SUNDIALS::core)
        endif ()
        foreach(comp IN LISTS SUNDIALS_COMPONENTS)
           target_link_libraries(amrex_${D}d PUBLIC SUNDIALS::${comp})
        endforeach()
        if (AMReX_GPU_BACKEND STREQUAL "CUDA")
           target_link_libraries(amrex_${D}d PUBLIC SUNDIALS::nveccuda)
        elseif (AMReX_GPU_BACKEND STREQUAL "HIP")
           target_link_libraries(amrex_${D}d PUBLIC SUNDIALS::nvechip)
        elseif (AMReX_GPU_BACKEND STREQUAL "SYCL")
           target_link_libraries(amrex_${D}d PUBLIC SUNDIALS::nvecsycl)
        else ()
           target_link_libraries(amrex_${D}d PUBLIC SUNDIALS::nvecserial)
        endif ()
    endforeach()
endif ()
