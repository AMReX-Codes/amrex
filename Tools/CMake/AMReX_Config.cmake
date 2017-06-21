###############################################
                                            
# Here we define the preprocessor definition  #

###############################################
set ( AMREX_DEFINES "" )

#
#  Check if AMReX_Options.cmake has been already processed
#
if ( NOT AMREX_OPTIONS_SET )
   message ( FATAL_ERROR "AMReX_Options.cmake must be\
included before AMReX_Configure.cmake" )
endif ()


#
# This is a macro for accumulating defines into AMREX_DEFINES
#
macro ( add_define NEW_DEFINE )
   message ("${ARGV}")
   if ( ${ARGC} EQUAL 1 )
      set ( AMREX_DEFINES  ${AMREX_DEFINES} " -D${NEW_DEFINE}" )
  #    message ("FROM ADD_DEFINES IF !${ARGV}|")
   else ()
      if ( ${ARGV1} ) 
#	 message ("FROM ADD_DEFINES !${ARGV}|")
	 set ( AMREX_DEFINES  ${AMREX_DEFINES} " -D${NEW_DEFINE}" )	 
      endif ()
   endif ()
endmacro ()

#
# Decide whether or not to use PIC 
#
if (ENABLE_PIC)
  set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()

#
# Set the C/Fortran naming scheme
#
include(FortranCInterface)
include(${FortranCInterface_BINARY_DIR}/Output.cmake)

if ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "" )  AND
     ( FortranCInterface_GLOBAL_CASE STREQUAL "UPPER") )
  message(STATUS "Fortran name mangling scheme to UPPERCASE \
(upper case, no append underscore)")
  add_define ("BL_FORT_USE_UPPERCASE" )
elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "") AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
   message(STATUS "Fortran name mangling scheme to LOWERCASE \
(lower case, no append underscore)")
   add_define ("BL_FORT_USE_LOWERCASE" )
elseif ( ( FortranCInterface_GLOBAL_SUFFIX STREQUAL "_" ) AND
      ( FortranCInterface_GLOBAL_CASE STREQUAL "LOWER") )
   message(STATUS "Fortran name mangling scheme to UNDERSCORE \
(lower case, append underscore)")
   add_define ("BL_FORT_USE_UNDERSCORE" )
else ()
   message(AUTHOR_WARNING "Fortran to C mangling not backward\
 compatible with older style BoxLib code") 
endif ()

#
# This defines were always on
#
add_define (BL_NOLINEVALUES)
add_define (BL_PARALLEL_IO)
add_define (BL_SPACEDIM=${BL_SPACEDIM})
add_define (BL_${CMAKE_SYSTEM_NAME})

if ( ENABLE_DP )
   add_define (BL_USE_DOUBLE)
else ()
   add_define (BL_USE_FLOAT)
endif ()


add_define (USE_PARTICLES ENABLE_PARTICLES)
add_define (BL_PROFILING ENABLE_PROFILING) 
add_define (BL_COMM_PROFILING ENABLE_COMM_PROFILING)

# if (ENABLE_PROFILING)
#   add_define ("BL_PROFILING")
# endif ()

# if (ENABLE_COMM_PROFILING)
#   add_define ("BL_COMM_PROFILING")
# endif ()

set(AMREX_EXTRA_LIBRARIES)
set(AMREX_EXTRA_LIBRARY_PATH)
set(AMREX_EXTRA_C_INCLUDE_PATH)
set(AMREX_EXTRA_CXX_INCLUDE_PATH)
set(AMREX_EXTRA_Fortran_INCLUDE_PATH)

if (ENABLE_MPI)
   add_define ("BL_USE_MPI")
   #  list(APPEND BL_DEFINES BL_USE_MPI)
  find_package(MPI REQUIRED)
  list(APPEND AMREX_EXTRA_Fortran_INCLUDE_PATH "${MPI_Fortran_INCLUDE_PATH}")
  list(APPEND AMREX_EXTRA_C_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
  list(APPEND AMREX_EXTRA_CXX_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
  set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} ${MPI_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_FLAGS}")
endif()

if (ENABLE_OpenMP)
  set(ENABLE_OMP TRUE)
  list(APPEND BL_DEFINES BL_USE_OMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (DEFINED EXTRA_DEFS_FOR_CCSE)
  list(APPEND BL_DEFINES "${EXTRA_DEFS_FOR_CCSE}")
endif (DEFINED EXTRA_DEFS_FOR_CCSE)

set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${BL_DEFINES}")


if ( NOT ENABLE_DP_PARTICLES ) 
   add_define ( "BL_SINGLE_PRECISION_PARTICLES" )
endif ()

if (ENABLE_FORTRAN_MPI AND ENABLE_MPI)
   add_define ( "BL_USE_FORTRAN_MPI=1" )
endif ()

if (ENABLE_MG_FBOXLIB)
   add_define ("MG_USE_FBOXIB=1")
endif ()

if (ENABLE_FBASELIB)
   add_define ("BL_USE_F_BASELIB=1")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth-64 -Wno-deprecated")
endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")


message ( STATUS "   AMREX_DEFINES = ${AMREX_DEFINES}" )

#
# Add all preprocessor definitions to compile string
# 
add_definitions ( ${AMREX_DEFINES} )
