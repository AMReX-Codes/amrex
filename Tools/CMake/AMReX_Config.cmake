###############################################

# Here we configure the build                 #

###############################################

#
#  Check if AMReX_Options.cmake has been already processed
#
if ( NOT AMREX_OPTIONS_SET )
   message ( FATAL_ERROR "AMReX_Options.cmake must be\
included before AMReX_Configure.cmake" )
endif ()

#
# Find AMReX Git version
#
find_git_version ( AMREX_GIT_VERSION )

#
# Decide whether or not to use PIC 
#
if (ENABLE_PIC)
   set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()

#
# Detect Fortran name mangling scheme for C/Fortran interface 
#
include(FortranCInterface)
include(${FortranCInterface_BINARY_DIR}/Output.cmake)

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


# ------------------------------------------------------------- #
#    Set preprocessor flags 
# ------------------------------------------------------------- #

#
# This defines were always on in older version
# of AMReX/CMake. Need more details on these???
#
add_define (BL_NOLINEVALUES AMREX_DEFINES)
add_define (AMREX_NOLINEVALUES AMREX_DEFINES)
add_define (BL_PARALLEL_IO AMREX_DEFINES)
add_define (AMREX_PARALLEL_IO AMREX_DEFINES)
add_define (BL_SPACEDIM=${BL_SPACEDIM} AMREX_DEFINES)
add_define (AMREX_SPACEDIM=${BL_SPACEDIM} AMREX_DEFINES)
add_define (BL_${CMAKE_SYSTEM_NAME} AMREX_DEFINES)
add_define (AMREX_${CMAKE_SYSTEM_NAME} AMREX_DEFINES)

if ( ENABLE_DP )
   add_define (BL_USE_DOUBLE AMREX_DEFINES)
   add_define (AMREX_USE_DOUBLE AMREX_DEFINES)
else ()
   add_define (BL_USE_FLOAT AMREX_DEFINES)
   add_define (AMREX_USE_FLOAT AMREX_DEFINES)
endif ()

add_define (USE_PARTICLES AMREX_DEFINES ENABLE_PARTICLES)

add_define (BL_PROFILING AMREX_DEFINES ENABLE_PROFILING)
add_define (AMREX_PROFILING AMREX_DEFINES ENABLE_PROFILING) 

add_define (BL_TINY_PROFILING AMREX_DEFINES ENABLE_TINY_PROFILING)
add_define (AMREX_TINY_PROFILING AMREX_DEFINES ENABLE_TINY_PROFILING)

add_define (BL_COMM_PROFILING AMREX_DEFINES ENABLE_COMM_PROFILING)
add_define (AMREX_COMM_PROFILING AMREX_DEFINES ENABLE_COMM_PROFILING)
# Comm profiling needs base profiling as well
add_define (BL_PROFILING AMREX_DEFINES ENABLE_COMM_PROFILING)
add_define (AMREX_PROFILING AMREX_DEFINES ENABLE_COMM_PROFILING)

add_define (BL_TRACE_PROFILING AMREX_DEFINES ENABLE_TRACE_PROFILING)
add_define (AMREX_TRACE_PROFILING AMREX_DEFINES ENABLE_TRACE_PROFILING)


if ( ENABLE_PARTICLES AND ( NOT ENABLE_DP_PARTICLES ) ) 
   add_define ( BL_SINGLE_PRECISION_PARTICLES AMREX_DEFINES )
   add_define ( AMREX_SINGLE_PRECISION_PARTICLES AMREX_DEFINES )
endif ()

if (ENABLE_FBASELIB)
   add_define ( BL_USE_FORTRAN_MPI AMREX_DEFINES ENABLE_MPI)
   add_define ( AMREX_USE_FORTRAN_MPI AMREX_DEFINES ENABLE_MPI)
   add_define ( BL_USE_F_BASELIB  AMREX_DEFINES )
   add_define ( AMREX_USE_F_BASELIB  AMREX_DEFINES )
endif ()

add_define (BL_USE_MPI AMREX_DEFINES ENABLE_MPI)
add_define (AMREX_USE_MPI AMREX_DEFINES ENABLE_MPI)

add_define (BL_USE_OMP AMREX_DEFINES ENABLE_OMP)
add_define (AMREX_USE_OMP AMREX_DEFINES ENABLE_OMP)

add_define (BL_USE_F_INTERFACES AMREX_DEFINES ENABLE_FORTRAN_INTERFACES)

add_define (BL_USE_ASSERTION AMREX_DEFINES ENABLE_ASSERTIONS) 

add_define (AMREX_GIT_VERSION=\\"${AMREX_GIT_VERSION}\\" AMREX_DEFINES)


#
# Add all preprocessor definitions to compile string
# 
add_definitions ( ${AMREX_DEFINES} )

# ------------------------------------------------------------- #
#    Setup third party packages 
# ------------------------------------------------------------- #

if (ENABLE_MPI)
   find_package (MPI REQUIRED)
   # Includes
   list (APPEND AMREX_EXTRA_Fortran_INCLUDE_PATH "${MPI_Fortran_INCLUDE_PATH}")
   list (APPEND AMREX_EXTRA_C_INCLUDE_PATH "${MPI_C_INCLUDE_PATH}")
   list (APPEND AMREX_EXTRA_CXX_INCLUDE_PATH "${MPI_CXX_INCLUDE_PATH}")
   # Compile flags
   append ( MPI_Fortran_COMPILE_FLAGS AMREX_EXTRA_Fortran_FLAGS ) 
   append ( MPI_C_COMPILE_FLAGS AMREX_EXTRA_C_FLAGS )
   append ( MPI_CXX_COMPILE_FLAGS AMREX_EXTRA_CXX_FLAGS )
   # Libraries
   list (APPEND AMREX_EXTRA_Fortran_LIBRARIES "${MPI_Fortran_LIBRARIES}")
   list (APPEND AMREX_EXTRA_C_LIBRARIES "${MPI_C_LIBRARIES}")
   list (APPEND AMREX_EXTRA_CXX_LIBRARIES "${MPI_CXX_LIBRARIES}")
   # Link flags
   list (APPEND AMREX_EXTRA_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
   list (APPEND AMREX_EXTRA_C_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
   list (APPEND AMREX_EXTRA_CXX_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
   # Link Line
   append_to_link_line ( AMREX_EXTRA_Fortran_LIBRARIES
      AMREX_EXTRA_Fortran_LINK_LINE AMREX_EXTRA_Fortran_LINK_FLAGS )
   append_to_link_line ( AMREX_EXTRA_C_LIBRARIES
      AMREX_EXTRA_C_LINK_LINE AMREX_EXTRA_C_LINK_FLAGS )
   append_to_link_line ( AMREX_EXTRA_CXX_LIBRARIES
      AMREX_EXTRA_CXX_LINK_LINE AMREX_EXTRA_CXX_LINK_FLAGS )
endif ()

if (ENABLE_OMP)
   find_package (OpenMP REQUIRED)
   # Compile flags
   append ( OpenMP_Fortran_FLAGS AMREX_EXTRA_Fortran_FLAGS ) 
   append ( OpenMP_C_FLAGS AMREX_EXTRA_C_FLAGS )
   append ( OpenMP_CXX_FLAGS AMREX_EXTRA_CXX_FLAGS )
endif()


# ------------------------------------------------------------- #
#    Setup compiler flags 
# ------------------------------------------------------------- #
if ( AMREX_FFLAGS_OVERRIDES )
   set ( AMREX_Fortran_FLAGS ${AMREX_FFLAGS_OVERRIDES} )
else ()
   append ( AMREX_${FC_ID}_FFLAGS_${AMREX_BUILD_TYPE}
      AMREX_Fortran_FLAGS )
endif ()

if ( AMREX_CXXFLAGS_OVERRIDES )
   set ( AMREX_CXX_FLAGS ${AMREX_CXXFLAGS_OVERRIDES} )
else ()
   append ( AMREX_${CXX_ID}_CXXFLAGS_${AMREX_BUILD_TYPE}
      AMREX_CXX_FLAGS )
endif ()

append ( AMREX_EXTRA_Fortran_FLAGS AMREX_Fortran_FLAGS )
append ( AMREX_EXTRA_CXX_FLAGS AMREX_CXX_FLAGS )

# Add required flags
append ( AMREX_${FC_ID}_FFLAGS_REQUIRED AMREX_Fortran_FLAGS )
append ( AMREX_${CXX_ID}_CXXFLAGS_REQUIRED AMREX_CXX_FLAGS )

# Add FPE flags if required 
if (ENABLE_FPE)
   append ( AMREX_${FC_ID}_FFLAGS_FPE AMREX_Fortran_FLAGS )
   append ( AMREX_${CXX_ID}_CXXFLAGS_FPE AMREX_CXX_FLAGS )
endif ()

# Set CMake compiler flags
set ( CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}
   "${AMREX_Fortran_FLAGS}  ${AMREX_Fortran_DEFINITIONS}" ) 
set ( CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE} "${AMREX_CXX_FLAGS}" )


#
# Config summary
#
message( STATUS "AMReX configuration summary: ")
message( STATUS "   Build type               = ${CMAKE_BUILD_TYPE}")
message( STATUS "   Install directory        = ${CMAKE_INSTALL_PREFIX}")
message( STATUS "   Preprocessor flags       = ${AMREX_DEFINES}")
message( STATUS "   C++ compiler             = ${CMAKE_CXX_COMPILER}")
message( STATUS "   Fortran compiler         = ${CMAKE_Fortran_COMPILER}")
message( STATUS "   C++ flags                = ${CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}}")
message( STATUS "   Fortran flags            = ${CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}}")
message( STATUS "   C++ include paths        = ${AMREX_EXTRA_CXX_INCLUDE_PATH}") 
message( STATUS "   Fortran include paths    = ${AMREX_EXTRA_Fortran_INCLUDE_PATH}")
message( STATUS "   C++ external libs        = ${AMREX_EXTRA_CXX_LIBRARIES}") 
message( STATUS "   Fortran external libs    = ${AMREX_EXTRA_Fortran_LIBRARIES}")
message( STATUS "   C++ link flags           = ${AMREX_EXTRA_CXX_LINK_FLAGS}") 
message( STATUS "   Fortran link flags       = ${AMREX_EXTRA_Fortran_LINK_FLAGS}")
message( STATUS "   C++ extra link line      = ${AMREX_EXTRA_CXX_LINK_LINE}") 
message( STATUS "   Fortran extra link line  = ${AMREX_EXTRA_Fortran_LINK_LINE}")
