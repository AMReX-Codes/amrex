#
#
#  This file provides:
#
#    AMREX_DEFINES:          list of the cpp flags to be used on ALL sources
#    AMREX_Fortran_DEFINES:  list of Fortran-specific cpp flags
#    add_define():           function to add definitions to AMREX_DEFINES 
#
#  Once this file is included, AMREX_DEFINES and AMREX_Fortran_DEFINES will be
#  populated with the preprocessor directives needed for a succesfull build.
#  Further CPP flags can be appended manually or via add_define() after this
#  file has been included.
# 
#  As per xSDK requirements, if a user set the env variable CPPFLAGS,
#  CPPFLAGS should overwrite AMReX_DEFINES. Since this is not possible without
#  breaking the build (most of the defines here are necessary for AMReX to compile),
#  for the time being we will not follow this requirement.
#  
# 


#
#  Check if AMReX_Options.cmake has been already processed
#
if ( NOT ( DEFINED __AMREX_OPTIONS__ ) )
   message ( FATAL_ERROR "AMReX_Options.cmake must be\
included before AMReX_Configure.cmake" )
endif ()

#
# AMREX_DEFINES will contain the CPP flags for all the sources
# For Fortran sources, use AMREX_Fortran_DEFINITIONS in addition to
# AMREX_DEFINES
# 
set ( AMREX_DEFINES )
set ( AMREX_Fortran_DEFINES  "-DBL_LANG_FORT -DAMREX_LANG_FORT" )

#
# Function to accumulate preprocessor directives
#
function ( add_define new_define )

   cmake_parse_arguments ( DEFINE "" "IF" ""  ${ARGN} )

   set ( condition  1 )

   if (DEFINE_IF)
      set ( condition ${${DEFINE_IF}} )
   endif()

   # Return if flags does not need to be included
   if (NOT condition)
      return ()
   endif ()
   
   set ( definition -D${new_define} )
   
   list ( FIND AMREX_DEFINES ${definition} out )
   
   if ( ${out} EQUAL -1 )
      list ( APPEND AMREX_DEFINES ${definition} )
      #  Add legacy definition
      string (FIND ${definition} "AMREX_" out )
      if (${out} GREATER -1 )
	 string (REPLACE "AMREX_" "BL_" legacy_definition ${definition})
	 list ( APPEND AMREX_DEFINES ${legacy_definition} )
      endif ()
      set ( AMREX_DEFINES ${AMREX_DEFINES} PARENT_SCOPE )
   endif ()
   
endfunction ()

# 
# Set preprocessor flags (trying to mimic GNUMake setup) 
#

# Git version
list ( APPEND AMREX_DEFINES -DAMREX_GIT_VERSION=\"${AMREX_GIT_VERSION}\" )

# XSDK mode
add_define ( AMREX_XSDK IF USE_XSDK_DEFAULTS )

# Debug flag
if (DEBUG)
   add_define (AMREX_DEBUG)
else ()
   add_define (NDEBUG)
endif ()

# Base profiling options
add_define ( AMREX_PROFILING IF ENABLE_BASE_PROFILE )

add_define ( AMREX_TRACE_PROFILING IF ENABLE_TRACE_PROFILE )

add_define ( AMREX_COMM_PROFILING  IF ENABLE_COMM_PROFILE )

# Tiny profiler
add_define ( AMREX_TINY_PROFILING IF ENABLE_TINY_PROFILE )

# Mem profiler 
add_define ( AMREX_MEM_PROFILING IF ENABLE_MEM_PROFILE )

# Profparser 
add_define ( AMREX_USE_PROFPARSER IF ENABLE_PROFPARSER )

# Backtrace
if (ENABLE_BACKTRACE)
   add_define ( AMREX_BACKTRACING )
   add_define ( AMREX_TESTING )
endif ()

# Third party profiling
if (${TP_PROFILE} MATCHES "CRAYPAT")
   add_define ( AMREX_CRAYPAT )
elseif (${TP_PROFILE} MATCHES "FORGE")
   add_define ( AMREX_FORGE )
elseif (${TP_PROFILE} MATCHES "VTUNE")
   add_define ( AMREX_VTUNE )
endif ()

# MPI
add_define ( AMREX_USE_MPI IF ENABLE_MPI )

# OpenMP
add_define ( AMREX_USE_OMP IF ENABLE_OMP )

# Precision
if (NOT ENABLE_DP)
   add_define (AMREX_USE_FLOAT)
endif ()

# Fortran/C mangling scheme
add_define ( AMREX_FORT_USE_${FORTLINK} )

# Dimensionality
add_define ( AMREX_SPACEDIM=${DIM} )

# System
add_define ( AMREX_${CMAKE_SYSTEM_NAME} )

# Particles
if (ENABLE_PARTICLES)
   add_define ( AMREX_PARTICLES )

   if (NOT ENABLE_DP_PARTICLES)
      add_define ( AMREX_SINGLE_PRECISION_PARTICLES )
   endif ()
endif ()

#  Assertions
add_define ( AMREX_USE_ASSERTION IF ENABLE_ASSERTION )

if (ENABLE_FBASELIB)
   add_define ( AMREX_USE_FORTRAN_MPI IF ENABLE_MPI)
   add_define ( AMREX_USE_F_BASELIB )
endif ()

add_define ( AMREX_USE_EB IF ENABLE_EB )

add_define ( AMREX_USE_F_INTERFACES IF ENABLE_FORTRAN_INTERFACES )

add_define ( AMREX_USE_ASSERTION IF ENABLE_ASSERTIONS ) 

add_define ( AMREX_NO_STRICT_PREFIX )


