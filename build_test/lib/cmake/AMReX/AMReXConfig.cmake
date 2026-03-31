# ############################################################################ #
#
#  AMReX Configuration File
#  To import into other CMake projects
#
# ############################################################################ #

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was AMReXConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################

# Custom version of check_required_components
# Set NO_CHECK_REQUIRED_COMPONENTS_MACRO when calling
# configure_package_config_file to avoid the CMake version
# of this macro to be generated
macro(check_required_components _NAME)
   foreach(comp ${${_NAME}_FIND_COMPONENTS})
      if(NOT ${_NAME}_${comp}_FOUND)
         if(${_NAME}_FIND_REQUIRED_${comp})
            message(STATUS "Requested AMReX component '${comp}' was not found.")
            set(${_NAME}_FOUND FALSE)
         endif()
      endif()
   endforeach()
endmacro()

# Set the minimum CMake version required -- This must be the version
# of CMake used to build the library
cmake_minimum_required(VERSION 3.31.6)

# Provides find_dependency
include(CMakeFindDependencyMacro)

#
# Build type
#
set(AMReX_BUILD_TYPE  Release)

#
# Versioning
#
set(AMReX_GIT_VERSION \"26.03\")

#
# Release number
#
set(AMReX_RELEASE_NUMBER 260300)

#
# AMReX CMake modules PATH
#
set_and_check(AMReX_MODULE_PATH ${PACKAGE_PREFIX_DIR}/lib/cmake/AMReX/AMReXCMakeModules)

#
# Add AMReX modules to app code CMake
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${AMReX_MODULE_PATH})

#
# Configuration options
# Each option is treated like a "component" so that find_package can be easily
# used to check weather the option is enabled
#

# General options
set(AMReX_MPI_FOUND                 OFF)
set(AMReX_MPI_THREAD_MULTIPLE_FOUND OFF)
set(AMReX_SIMD_FOUND                OFF)
set(AMReX_OMP_FOUND                 OFF)
set(AMReX_CUDA_FOUND                OFF)
set(AMReX_SYCL_FOUND                OFF)
set(AMReX_HIP_FOUND                 OFF)
set(AMReX_DOUBLE_FOUND   ON)
set(AMReX_FORTRAN_FOUND             OFF)

# Actual components selection
set(AMReX_AMRLEVEL_FOUND            ON)
set(AMReX_EB_FOUND                  OFF)
set(AMReX_FINTERFACES_FOUND         OFF)
set(AMReX_LSOLVERS_FOUND            ON)
set(AMReX_LSOLVERS_INCFLO_FOUND     ON)
set(AMReX_LSOLVERS_EM_FOUND         ON)
set(AMReX_FFT_FOUND                 OFF)
set(AMReX_AMRDATA_FOUND             OFF)
set(AMReX_PARTICLES_FOUND           ON)
set(AMReX_PDOUBLE_FOUND ON)
set(AMReX_SENSEI_FOUND              OFF)
set(AMReX_CONDUIT_FOUND             OFF)
set(AMReX_CATALYST_FOUND            OFF)
set(AMReX_ASCENT_FOUND              OFF)
set(AMReX_HYPRE_FOUND               OFF)
set(AMReX_PETSC_FOUND               OFF)
set(AMReX_SUNDIALS_FOUND            OFF)
set(AMReX_HDF5_FOUND                OFF)
set(AMReX_HDF5_ZFP_FOUND            OFF)

# Compilation options
set(AMReX_FPE_FOUND                 OFF)
set(AMReX_PIC_FOUND                 OFF)
set(AMReX_ASSERTIONS_FOUND          OFF)
set(AMReX_FLATTEN_FOR_FOUND         OFF)
set(AMReX_COMPILER_DEFAULT_INLINE_FOUND   OFF)
set(AMReX_INLINE_LIMIT_FOUND              43210)

# Profiling options
set(AMReX_BASEP_FOUND               OFF)
set(AMReX_TINYP_FOUND               OFF)
set(AMReX_TRACEP_FOUND              OFF)
set(AMReX_MEMP_FOUND                OFF)
set(AMReX_COMMP_FOUND               OFF)
set(AMReX_PROFPARSER_FOUND          OFF)

# Plotfile tools
set(AMReX_PFTOOLS_FOUND             OFF)


# export the actual values as well.
# General options
set(AMReX_SPACEDIM                  3)
foreach(D IN LISTS AMReX_SPACEDIM)
    set(AMReX_${D}D_FOUND   ON)
endforeach()
set(AMReX_MPI                       OFF)
set(AMReX_MPI_THREAD_MULTIPLE       OFF)
set(AMReX_SIMD                      OFF)
set(AMReX_OMP                       OFF)
set(AMReX_CUDA                      OFF)
set(AMReX_SYCL                      OFF)
set(AMReX_HIP                       OFF)
set(AMReX_GPU_BACKEND               NONE)
set(AMReX_GPU_RDC                   OFF)
set(AMReX_PRECISION                 DOUBLE)
set(AMReX_FASTMATH                  OFF)
set(AMReX_FORTRAN                   OFF)

# Actual components selection
set(AMReX_AMRLEVEL                  ON)
set(AMReX_EB                        OFF)
set(AMReX_FINTERFACES               OFF)
set(AMReX_LSOLVERS                  ON)
set(AMReX_LSOLVERS_INCFLO           ON)
set(AMReX_LSOLVERS_EM               ON)
set(AMReX_FFT                       OFF)
set(AMReX_AMRDATA                   OFF)
set(AMReX_PARTICLES                 ON)
set(AMReX_PARTICLES_PRECISION       DOUBLE)
set(AMReX_SENSEI                    OFF)
set(AMReX_NO_SENSEI_AMR_INST        FALSE)
set(AMReX_CONDUIT                   OFF)
set(AMReX_CATALYST                  OFF)
set(AMReX_ASCENT                    OFF)
set(AMReX_HYPRE                     OFF)
set(AMReX_PETSC                     OFF)
set(AMReX_SUNDIALS                  OFF)
set(AMReX_HDF5                      OFF)
set(AMReX_HDF5_ZFP                  OFF)

# Compilation options
set(AMReX_FPE                       OFF)
set(AMReX_PIC                       OFF)
set(AMReX_ASSERTIONS                OFF)
set(AMReX_FLATTEN_FOR               OFF)

# Profiling options
set(AMReX_BASE_PROFILE              OFF)
set(AMReX_TINY_PROFILE              OFF)
set(AMReX_TRACE_PROFILE             OFF)
set(AMReX_MEM_PROFILE               OFF)
set(AMReX_COMM_PROFILE              OFF)
set(AMReX_PROFPARSER                OFF)
set(AMReX_PROFILE_FTOOLS            OFF)

#
# If Fortran is enabled, downstream project
# must have Fortran enabled as well
#
if (OFF AND NOT CMAKE_Fortran_COMPILER_LOADED )
   message(FATAL_ERROR
      "\nAMReX was build with AMReX_FORTRAN=ON but Fortran is not enabled for this project. "
      "Either set enable_language(Fortran) before importing AMReX or re-build AMReX with "
      "AMReX_FORTRAN=OFF.\n")
endif ()

#
# Parallel backends
#
set( THREADS_PREFER_PTHREAD_FLAG on)
find_dependency(Threads REQUIRED)

if (OFF)
   set( _mpi_components C CXX )
   if (OFF)
      list(APPEND _mpi_components Fortran)
   endif ()
   find_dependency(MPI REQUIRED ${_mpi_components})
   unset(_mpi_components)
endif()

if (OFF)
   set( _omp_components CXX )
   if (OFF)
      list(APPEND _omp_components Fortran)
   endif ()
   find_dependency(OpenMP REQUIRED ${_omp_components})
endif ()

#
# Third party libraries
#
if (OFF)
   find_dependency(vir-simd REQUIRED)
endif ()

if (OFF)
    set(SENSEI_DIR )
    find_dependency(SENSEI REQUIRED)
endif ()

if (OFF)
    find_dependency(Ascent REQUIRED)
endif ()

if (OFF)
    find_dependency(Catalyst REQUIRED)
endif ()

if (OFF)
   find_dependency(Conduit REQUIRED)
endif ()

if (OFF)
    if (NONE STREQUAL NONE)
        find_dependency(AMReXFFTW REQUIRED)
    endif()
endif()

if (OFF)
    find_dependency(HDF5 REQUIRED)
endif ()

if (OFF)
    find_dependency(H5Z_ZFP REQUIRED)
endif ()

if (OFF)
   find_dependency(HYPRE 2.20.0 REQUIRED)
endif ()

if (OFF)
   find_dependency(PETSc 2.13 REQUIRED)
endif ()

if (OFF)
   find_dependency(SUNDIALS 6.0.0 REQUIRED)
endif ()

#
# CUDA
#
# AMReX 21.06+ supports CUDA_ARCHITECTURES
if (OFF)
    if (CMAKE_VERSION VERSION_LESS 3.20)
        include(AMReX_SetupCUDA)
    else ()
        find_dependency(CUDAToolkit REQUIRED)
    endif ()
endif ()

# CMake targets
include( "${CMAKE_CURRENT_LIST_DIR}/AMReXTargets.cmake" )

# CMake targets aliases: last dimension built will be our legacy target
if (NOT TARGET AMReX::amrex)  # protection in case of multiple inclusions
    list(LENGTH AMReX_SPACEDIM list_len)
    math(EXPR list_last "${list_len} - 1")
    list(GET AMReX_SPACEDIM ${list_last} AMReX_SPACEDIM_LAST)
    add_library(AMReX::amrex ALIAS AMReX::amrex_${AMReX_SPACEDIM_LAST}d)
endif()

# More Modern CUDA CMake
if (3.31.6 VERSION_GREATER_EQUAL 3.20 AND OFF)
    foreach(D IN LISTS AMReX_SPACEDIM)
        # CUDA architectures amrex was built for -- should we make
        set(AMREX_CUDA_ARCHS  CACHE INTERNAL "CUDA archs AMReX is built for")
        set_target_properties(AMReX::amrex_${D}d
          PROPERTIES
            CUDA_ARCHITECTURES ${AMREX_CUDA_ARCHS})
    endforeach()
endif ()

#
# Check components
#
check_required_components("AMReX")
