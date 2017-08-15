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

if ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "" )  AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "UPPER") )
   message(STATUS "Fortran name mangling scheme to UPPERCASE \
(upper case, no append underscore)")
   add_define ( BL_FORT_USE_UPPERCASE AMREX_DEFINES )
elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "") AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
   message(STATUS "Fortran name mangling scheme to LOWERCASE \
(lower case, no append underscore)")
   add_define ( BL_FORT_USE_LOWERCASE AMREX_DEFINES )
elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "_" ) AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
   message(STATUS "Fortran name mangling scheme to UNDERSCORE \
(lower case, append underscore)")
   add_define ( BL_FORT_USE_UNDERSCORE AMREX_DEFINES )
else ()
   message(AUTHOR_WARNING "Fortran to C mangling not backward\
 compatible with older style BoxLib code") 
endif ()


# 
# Set preprocessor flags 
# 
add_define ( AMREX_NOLINEVALUES )
add_define ( AMREX_PARALLEL_IO )
add_define ( AMREX_SPACEDIM=${BL_SPACEDIM} )
add_define ( AMREX_${CMAKE_SYSTEM_NAME} )

if ( ENABLE_DP )
   add_define ( AMREX_USE_DOUBLE )
else ()
   add_define ( AMREX_USE_FLOAT )
endif ()

add_define ( USE_PARTICLES IF ENABLE_PARTICLES )

# Profiling flags 
add_define ( AMREX_PROFILING       IF AMREX_ANY_PROFILING ) 
add_define ( AMREX_TINY_PROFILING  IF ENABLE_TINY_PROFILING )
add_define ( AMREX_COMM_PROFILING  IF ENABLE_COMM_PROFILING )
add_define ( AMREX_TRACE_PROFILING IF ENABLE_TRACE_PROFILING )
add_define ( AMREX_MEM_PROFILING   IF ENABLE_MEM_PROFILING )

if ( ENABLE_PARTICLES AND ( NOT ENABLE_DP_PARTICLES ) ) 
   add_define ( AMREX_SINGLE_PRECISION_PARTICLES )
endif ()

if (ENABLE_FBASELIB)
   add_define ( AMREX_USE_FORTRAN_MPI IF ENABLE_MPI)
   add_define ( AMREX_USE_F_BASELIB )
endif ()


add_define ( AMREX_USE_MPI IF ENABLE_MPI )

add_define ( AMREX_USE_OMP IF ENABLE_OMP )

add_define ( AMREX_USE_F_INTERFACES IF ENABLE_FORTRAN_INTERFACES )

add_define ( AMREX_USE_ASSERTION IF ENABLE_ASSERTIONS ) 

# Here I use explicitly a list because BL_GIT_VERSION (legacy define) does not
# exist
list ( APPEND AMREX_DEFINES -DAMREX_GIT_VERSION=\"${AMREX_GIT_VERSION}\" )

#
# Add all preprocessor definitions to compile string
# 
add_definitions ( ${AMREX_DEFINES} )


