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

# Mem profiler
add_amrex_define( AMREX_MEM_PROFILING NO_LEGACY IF AMReX_MEM_PROFILE )

# MPI
add_amrex_define( AMREX_USE_MPI IF AMReX_MPI )
add_amrex_define( AMREX_MPI_THREAD_MULTIPLE NO_LEGACY IF AMReX_MPI_THREAD_MULTIPLE)

# OpenMP -- This one has legacy definition only in Base/AMReX_omp_mod.F90
add_amrex_define( AMREX_USE_OMP IF AMReX_OMP )

# DPCPP
add_amrex_define( AMREX_USE_DPCPP NO_LEGACY IF AMReX_DPCPP )
add_amrex_define( AMREX_USE_GPU NO_LEGACY IF AMReX_DPCPP )

# HIP
add_amrex_define( AMREX_USE_HIP NO_LEGACY IF AMReX_HIP )
add_amrex_define( NDEBUG IF AMReX_HIP)  # This address a bug that causes slow build times

# Precision
if (AMReX_PRECISION STREQUAL "SINGLE")
   add_amrex_define(AMREX_USE_FLOAT)
endif ()

# Dimensionality
add_amrex_define( AMREX_SPACEDIM=${AMReX_SPACEDIM} )

# System -- not used anywhere in the source code
add_amrex_define( AMREX_${CMAKE_SYSTEM_NAME} )

#  Assertions
add_amrex_define( AMREX_USE_ASSERTION NO_LEGACY IF AMReX_ASSERTIONS )

# Bound checking
add_amrex_define( AMREX_BOUND_CHECK NO_LEGACY IF AMReX_BOUND_CHECK )

#
# Fortran-specific defines: BL_LANG_FORT and AMREX_LANG_FORT
#
if (AMReX_FORTRAN)

   # These defines are needed only by AMReX source files
   target_compile_definitions( amrex PRIVATE
      $<$<COMPILE_LANGUAGE:Fortran>:AMREX_SPACEDIM=${AMReX_SPACEDIM}>
      $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT AMREX_LANG_FORT>
      )

   #
   # Fortran/C mangling scheme
   #
   include( FortranCInterface )
   if(NOT FortranCInterface_GLOBAL_FOUND)
      message(FATAL_ERROR "Failed to find the Fortan C Interface -- check the CMakeError.log")
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

# SENSEI Insitu -- only legacy
add_amrex_define( BL_USE_SENSEI_INSITU IF AMReX_SENSEI )

# Conduit Support
add_amrex_define( AMREX_USE_CONDUIT NO_LEGACY IF AMReX_CONDUIT )

# Ascent Support
add_amrex_define( AMREX_USE_ASCENT NO_LEGACY IF AMReX_ASCENT )

# EB
add_amrex_define( AMREX_USE_EB NO_LEGACY IF AMReX_EB )

#
# CUDA
#
add_amrex_define( AMREX_USE_CUDA NO_LEGACY IF AMReX_CUDA )
add_amrex_define( AMREX_USE_NVML NO_LEGACY IF AMReX_CUDA )
add_amrex_define( AMREX_GPU_MAX_THREADS=${AMReX_CUDA_MAX_THREADS} NO_LEGACY
   IF AMReX_CUDA )

#
# General setup for any GPUs
#
if (AMReX_CUDA OR AMReX_HIP)
   add_amrex_define( AMREX_USE_GPU  NO_LEGACY )
   add_amrex_define( BL_COALESCE_FABS )

   add_amrex_define( AMREX_GPUS_PER_SOCKET=${GPUS_PER_SOCKET}
      NO_LEGACY IF GPUS_PER_SOCKET)

   add_amrex_define( AMREX_GPUS_PER_NODE=${GPUS_PER_NODE}
      NO_LEGACY IF GPUS_PER_NODE)
endif ()

#
# HDF5
#
add_amrex_define(AMREX_USE_HDF5 NO_LEGACY IF AMReX_HDF5)
add_amrex_define(AMREX_USE_HDF5_ASYNC NO_LEGACY IF AMReX_HDF5_ASYNC)

#
# Miscellaneous
#
add_amrex_define( AMREX_NO_PROBINIT NO_LEGACY IF_NOT AMReX_PROBINIT)

#
# Windows DLLs and Global Symbols
# https://stackoverflow.com/questions/54560832/cmake-windows-export-all-symbols-does-not-cover-global-variables/54568678#54568678
#
if(WIN32 AND BUILD_SHARED_LIBS)
  add_amrex_define(AMREX_IS_DLL NO_LEGACY)
  target_compile_definitions( amrex PRIVATE AMREX_IS_DLL_BUILDING)
endif()
