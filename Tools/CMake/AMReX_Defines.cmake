###############################################

# Here we set the AMReX preprocessor flags    #
# This file provides:                         #
#      AMREX_DEFINES = string with macro defs #

###############################################

#
#  Check if AMReX_Options.cmake has been already processed
#
if ( NOT AMREX_OPTIONS_SET )
   message ( FATAL_ERROR "AMReX_Options.cmake must be\
included before AMReX_Configure.cmake" )
endif ()

#
# Defines variable
#
set ( AMREX_DEFINES "" )

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
# Detect Fortran name mangling scheme for C/Fortran interface 
#
include ( FortranCInterface )
include ( ${FortranCInterface_BINARY_DIR}/Output.cmake )

set (FORTLINK "")

if ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "" )

   set (FORTLINK "${FortranCInterface_GLOBAL_CASE}CASE" )
   message (STATUS "Fortran name mangling scheme: ${FORTLINK} (no append underscore)")

elseif ( (FortranCInterface_GLOBAL_SUFFIX STREQUAL "_")  AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER" ) )

   set (FORTLINK "UNDERSCORE")
   message (STATUS "Fortran name mangling scheme: ${FORTLINK} (lower case, append underscore)")

else ()
   message (AUTHOR_WARNING "Fortran to C mangling not compatible with AMReX code")
endif ()


# 
# Set preprocessor flags (trying to mimic GNUMake setup) 
#

# Git version
list ( APPEND AMREX_DEFINES -DAMREX_GIT_VERSION=\"${AMREX_GIT_VERSION}\" )

# Debug flag
if (DEBUG)
   add_define (DEBUG)
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
   add_define ( USE_PARTICLES )
   add_define ( PARTICLES )

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

add_define ( AMREX_USE_F_INTERFACES IF ENABLE_FORTRAN_INTERFACES )

add_define ( AMREX_USE_ASSERTION IF ENABLE_ASSERTIONS ) 

#
# Add all preprocessor definitions to compile string
# 
add_definitions ( ${AMREX_DEFINES} )


