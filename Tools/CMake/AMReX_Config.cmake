###############################################

# Here we configure the build                 #

###############################################
if ( DEFINED __AMREX_CONFIG__ )
   return ()
endif ()
set ( __AMREX_CONFIG__ )

#
# Check if AMReX_Options.cmake  and AMREX_CMakeVariables.cmake
# have been already processed
#
if ( NOT (DEFINED __AMREX_OPTIONS__ ) )
   message ( FATAL_ERROR "AMReX_Options.cmake must be\
included before AMReX_Configure.cmake" )
endif ()

if ( NOT (DEFINED __AMREX_CMAKEVARIABLES__ ) )
   message ( FATAL_ERROR "AMReX_CMakeVariables.cmake must be\
included before AMReX_Configure.cmake" )
endif ()

#
# Decide whether or not to use PIC 
#
if (ENABLE_PIC)
   set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()

# 
#  Setup MPI 
# 
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

#
# Setup OpenMP
# 
if (ENABLE_OMP)
   find_package (OpenMP REQUIRED)
   # Compile flags
   append ( OpenMP_Fortran_FLAGS AMREX_EXTRA_Fortran_FLAGS ) 
   append ( OpenMP_C_FLAGS AMREX_EXTRA_C_FLAGS )
   append ( OpenMP_CXX_FLAGS AMREX_EXTRA_CXX_FLAGS )
else ()
   if ( ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray" ) # Cray has OMP on by default
      list ( APPEND AMREX_EXTRA_Fortran_FLAGS  "-h noomp")
   endif ()
   if ( ${CMAKE_C_COMPILER_ID} STREQUAL "Cray" ) # Cray has OMP on by default
      list ( APPEND AMREX_EXTRA_C_FLAGS  "-h noomp")
   endif ()
   if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" ) # Cray has OMP on by default
      list ( APPEND AMREX_EXTRA_CXX_FLAGS  "-h noomp")
   endif ()
endif()

#
# Setup third-party profilers
#
include ( AMReX_ThirdPartyProfilers )

# Includes
list (APPEND AMREX_EXTRA_Fortran_INCLUDE_PATH "${TPP_Fortran_INCLUDE_PATH}")
list (APPEND AMREX_EXTRA_C_INCLUDE_PATH "${TPP_C_INCLUDE_PATH}")
list (APPEND AMREX_EXTRA_CXX_INCLUDE_PATH "${TPP_CXX_INCLUDE_PATH}")

# Compile flags
append ( TPP_FFLAGS AMREX_EXTRA_Fortran_FLAGS ) 
append ( TPP_CFLAGS AMREX_EXTRA_C_FLAGS )
append ( TPP_CXXFLAGS AMREX_EXTRA_CXX_FLAGS )

# Link Line
append_to_link_line ( TPP_Fortran_LINK_LINE AMREX_EXTRA_Fortran_LINK_LINE )
append_to_link_line ( TPP_C_LINK_LINE AMREX_EXTRA_C_LINK_LINE )
append_to_link_line ( TPP_CXX_LINK_LINE AMREX_EXTRA_CXX_LINK_LINE )

#
# Set preprocessor flags
#
include ( AMReX_Defines )

#
# Setup compiler flags 
#
include ( AMReX_Compilers )

# If CMAKE_<LANG>_FLAGS is not defined by user, load defaults
if ( NOT CMAKE_CXX_FLAGS )
   set ( CMAKE_CXX_FLAGS "${AMREX_CXXFLAGS_${AMREX_BUILD_TYPE}}" )
endif ()

if ( NOT CMAKE_Fortran_FLAGS )
   set ( CMAKE_Fortran_FLAGS "${AMREX_FFLAGS_${AMREX_BUILD_TYPE}}" )
endif ()

# Use CMAKE_<LANG>_FLAGS_<BUILD_TYPE> to store required flags
# NOTE: this is not standard CMake way of doing things
set ( CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}  "${AMREX_CXXFLAGS_REQUIRED}" )
set ( CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}  "${AMREX_FFLAGS_REQUIRED}" )

if ( ENABLE_FPE )
   append ( AMREX_CXXFLAGS_FPE CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE} )
   append ( AMREX_FFLAGS_FPE   CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE} )
endif ()

      
# Add all preprocessor definitions to compile string
# For fortran, this goes into CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE} since
# it is a required flag
append ( AMREX_Fortran_DEFINES CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE} )

#  Add definition related to specific compiler ( only GNU for now )
if ( AMREX_COMPILER_DEFINES )
   foreach ( item ${AMREX_COMPILER_DEFINES} )
      add_define (${item})
   endforeach ()
endif ()

add_definitions ( ${AMREX_DEFINES} )

# Add extra flags
append ( AMREX_EXTRA_Fortran_FLAGS CMAKE_Fortran_FLAGS )
append ( AMREX_EXTRA_CXX_FLAGS CMAKE_CXX_FLAGS )

# Accumulate all the flags into AMREX_<LANG>_FLAGS: this variables will be exported
set ( AMREX_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}}" )
set ( AMREX_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}}" )

#
# Config summary
#
message( STATUS "AMReX configuration summary: ")
message( STATUS "   Build type               = ${CMAKE_BUILD_TYPE}")
message( STATUS "   Install directory        = ${CMAKE_INSTALL_PREFIX}")
message( STATUS "   Preprocessor flags       = ${AMREX_DEFINES}")
message( STATUS "   C++ compiler             = ${CMAKE_CXX_COMPILER}")
message( STATUS "   Fortran compiler         = ${CMAKE_Fortran_COMPILER}")
message( STATUS "   C++ flags                = ${AMREX_CXX_FLAGS}")
message( STATUS "   Fortran flags            = ${AMREX_Fortran_FLAGS}")
message( STATUS "   C++ include paths        = ${AMREX_EXTRA_CXX_INCLUDE_PATH}") 
message( STATUS "   Fortran include paths    = ${AMREX_EXTRA_Fortran_INCLUDE_PATH}")
message( STATUS "   C++ external libs        = ${AMREX_EXTRA_CXX_LIBRARIES}") 
message( STATUS "   Fortran external libs    = ${AMREX_EXTRA_Fortran_LIBRARIES}")
message( STATUS "   C++ link flags           = ${AMREX_EXTRA_CXX_LINK_FLAGS}") 
message( STATUS "   Fortran link flags       = ${AMREX_EXTRA_Fortran_LINK_FLAGS}")
message( STATUS "   C++ extra link line      = ${AMREX_EXTRA_CXX_LINK_LINE}") 
message( STATUS "   Fortran extra link line  = ${AMREX_EXTRA_Fortran_LINK_LINE}")
