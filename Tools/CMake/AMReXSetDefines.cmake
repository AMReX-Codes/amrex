# Populate the AMREX_DEFINES property of target "amrex".
# Target "amrex" must exist before calling this function.
# This module sets ONLY the general compile definitions.
# That means it won't add defines associated to specific compilers or
# required by optional AMReX components, like the EB or Particle code for
# example.
#
# As per xSDK requirements, if a user set the env variable CPPFLAGS,
# CPPFLAGS should overwrite AMReX_DEFINES. Since this is not possible without
# breaking the build (most of the defines here are necessary for AMReX to compile),
# for the time being we will not follow this requirement.
#
include_guard(GLOBAL)
include(AMReXGenerateConfigHeader)  # provides add_amrex_defines


# Git version
add_amrex_define( "AMREX_GIT_VERSION=\"${AMREX_GIT_VERSION}\"" NO_LEGACY )

# Release number
add_amrex_define( "AMREX_RELEASE_NUMBER=${AMREX_RELEASE_NUMBER}" NO_LEGACY )

# XSDK mode
add_amrex_define( AMREX_XSDK NO_LEGACY IF USE_XSDK_DEFAULTS )

# Debug flag
if ( "${CMAKE_BUILD_TYPE}" MATCHES "Debug" )
   add_amrex_define( AMREX_DEBUG NO_LEGACY )
endif ()

# Base profiling options
add_amrex_define( AMREX_PROFILING       IF AMReX_BASE_PROFILE )
add_amrex_define( AMREX_TRACE_PROFILING IF AMReX_TRACE_PROFILE )
add_amrex_define( AMREX_COMM_PROFILING  IF AMReX_COMM_PROFILE )

# Tiny profiler
add_amrex_define( AMREX_TINY_PROFILING NO_LEGACY IF AMReX_TINY_PROFILE )
add_amrex_define( AMREX_USE_ROCTX NO_LEGACY IF AMReX_ROCTX )

# Mem profiler
add_amrex_define( AMREX_MEM_PROFILING NO_LEGACY IF AMReX_MEM_PROFILE )

# Testing
add_amrex_define( AMREX_TESTING NO_LEGACY IF AMReX_TESTING )

# MPI
add_amrex_define( AMREX_USE_MPI IF AMReX_MPI )
add_amrex_define( AMREX_MPI_THREAD_MULTIPLE NO_LEGACY IF AMReX_MPI_THREAD_MULTIPLE)

# OpenMP -- This one has legacy definition only in Base/AMReX_omp_mod.F90
add_amrex_define( AMREX_USE_OMP IF AMReX_OMP )

# SYCL
if (AMReX_SYCL)
  add_amrex_define( AMREX_USE_SYCL NO_LEGACY )
  add_amrex_define( AMREX_USE_DPCPP NO_LEGACY )
  add_amrex_define( AMREX_SYCL_SUB_GROUP_SIZE=${AMReX_SYCL_SUB_GROUP_SIZE} NO_LEGACY )
  add_amrex_define( AMREX_USE_ONEDPL NO_LEGACY IF AMReX_SYCL_ONEDPL )
endif()

# HIP
add_amrex_define( AMREX_USE_HIP NO_LEGACY IF AMReX_HIP )
add_amrex_define( NDEBUG IF AMReX_HIP)  # This address a bug that causes slow build times

# Precision
if (AMReX_PRECISION STREQUAL "SINGLE")
   add_amrex_define(AMREX_USE_FLOAT)
endif ()

# Dimensionality
add_amrex_define( AMREX_SPACEDIM=${AMReX_SPACEDIM} )
foreach(D IN LISTS AMReX_SPACEDIM)
    target_compile_definitions(amrex_${D}d PUBLIC AMREX_SPACEDIM=${D})
endforeach()

# System -- not used anywhere in the source code
add_amrex_define( AMREX_${CMAKE_SYSTEM_NAME} )

#  Assertions
add_amrex_define( AMREX_USE_ASSERTION NO_LEGACY IF AMReX_ASSERTIONS )

# Bound checking
add_amrex_define( AMREX_BOUND_CHECK NO_LEGACY IF AMReX_BOUND_CHECK )

# Backtraces on macOS
add_amrex_define( AMREX_EXPORT_DYNAMIC NO_LEGACY IF AMReX_EXPORT_DYNAMIC )

if (AMReX_FORTRAN)

   # Fortran-specific defines, BL_LANG_FORT and AMREX_LANG_FORT do not get
   # stored in AMReX_Config.H
   foreach(D IN LISTS AMReX_SPACEDIM)
       target_compile_definitions(amrex_${D}d PRIVATE
          $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT AMREX_LANG_FORT>
          )
   endforeach()

   #
   # Fortran/C mangling scheme
   #
   include( FortranCInterface )
   if(NOT FortranCInterface_GLOBAL_FOUND)
      message(FATAL_ERROR "Failed to find the Fortran C Interface -- check the CMakeError.log")
   endif()
   include( ${FortranCInterface_BINARY_DIR}/Output.cmake )

   set( FORTLINK "" )

   if ( (FortranCInterface_GLOBAL_SUFFIX STREQUAL "" ) AND NOT
         (FortranCInterface_GLOBAL_CASE STREQUAL "") )
      set(FORTLINK "${FortranCInterface_GLOBAL_CASE}CASE" )
      message(STATUS "Fortran name mangling scheme: ${FORTLINK} (no append underscore)")
   elseif ( (FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")  AND
	 ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER" ) )
      set(FORTLINK "UNDERSCORE")
      message(STATUS "Fortran name mangling scheme: ${FORTLINK} (lower case, append underscore)")
   else ()
      # now we have to guess
      if (CMAKE_Fortran_COMPILER_ID MATCHES XL) # old IBM prior to XLClang
         set(FORTLINK "LOWERCASE")
      else ()
         set(FORTLINK "UNDERSCORE")
      endif()
      message(WARNING "Fortran to C mangling not compatible with AMReX code, assuming '${FORTLINK}'")
   endif ()

   add_amrex_define( BL_FORT_USE_${FORTLINK} )  # Only legacy form
else ()
   add_amrex_define(BL_NO_FORT)
endif ()

#
# SENSEI Insitu
#
add_amrex_define( AMREX_USE_SENSEI_INSITU NO_LEGACY IF AMReX_SENSEI )
add_amrex_define( AMREX_NO_SENSEI_AMR_INST NO_LEGACY IF AMReX_NO_SENSEI_AMR_INST )

# Conduit Support
add_amrex_define( AMREX_USE_CONDUIT NO_LEGACY IF AMReX_CONDUIT )

# Ascent Support
add_amrex_define( AMREX_USE_ASCENT NO_LEGACY IF AMReX_ASCENT )

#
# CUDA
#
add_amrex_define( AMREX_USE_CUDA NO_LEGACY IF AMReX_CUDA )
add_amrex_define( AMREX_USE_NVML NO_LEGACY IF AMReX_CUDA )

#
# General setup for any GPUs
#
if (NOT AMReX_GPU_BACKEND STREQUAL NONE)
   add_amrex_define( AMREX_USE_GPU NO_LEGACY )
   add_amrex_define( AMREX_GPU_MAX_THREADS=${AMReX_GPU_MAX_THREADS} NO_LEGACY )
   add_amrex_define( BL_COALESCE_FABS )
endif()

if (AMReX_CUDA OR AMReX_HIP)
   add_amrex_define( AMREX_USE_GPU_RDC NO_LEGACY IF AMReX_GPU_RDC )
endif ()

#
# HDF5
#
add_amrex_define(AMREX_USE_HDF5 NO_LEGACY IF AMReX_HDF5)
add_amrex_define(AMREX_USE_HDF5_ASYNC NO_LEGACY IF AMReX_HDF5_ASYNC)
add_amrex_define(AMREX_USE_HDF5_ZFP NO_LEGACY IF AMReX_HDF5_ZFP)


#
# SUNDIALS
#
add_amrex_define(AMREX_USE_SUNDIALS NO_LEGACY IF AMReX_SUNDIALS)

#
# Miscellaneous
#
add_amrex_define( AMREX_NO_PROBINIT NO_LEGACY IF_NOT AMReX_PROBINIT)

#
# Windows DLLs and Global Symbols
# https://stackoverflow.com/questions/54560832/cmake-windows-export-all-symbols-does-not-cover-global-variables/54568678#54568678
#
if(WIN32 AND AMReX_BUILD_SHARED_LIBS)
    add_amrex_define(AMREX_IS_DLL NO_LEGACY)
    foreach(D IN LISTS AMReX_SPACEDIM)
        target_compile_definitions(amrex_${D}d PRIVATE AMREX_IS_DLL_BUILDING)
    endforeach()
endif()
