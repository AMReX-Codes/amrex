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
# This is a macro for accumulating defines into AMREX_DEFINES
#
macro ( add_define NEW_DEFINE )
   if ( ${ARGC} EQUAL 1 )
      set ( AMREX_DEFINES  ${AMREX_DEFINES} " -D${NEW_DEFINE}" )
   else ()
      if ( ${ARGV1} ) 
	 set ( AMREX_DEFINES  ${AMREX_DEFINES} " -D${NEW_DEFINE}" )	 
      endif ()
   endif ()
endmacro ()

#
# This is a function for accumulating compiler flags
#
function ( add_compiler_flags NEW_FLAGS ALL_FLAGS )
   if ( ${NEW_FLAGS} )
      set ( ${ALL_FLAGS}  "${${ALL_FLAGS}} ${${NEW_FLAGS}}" PARENT_SCOPE )
   endif ()
endfunction ()

#
# Decide whether or not to use PIC 
#
if (ENABLE_PIC)
   set (CMAKE_POSITION_INDEPENDENT_CODE TRUE)
endif ()


# No idea why we need this.
# I think it was required for Franklin build. -- lpritch
# if(PREFER_STATIC_LIBRARIES)
#   # Prefer static libraries, but don't require that everything must be static. 
#   # This appears to be necessary on Franklin at NERSC right now.  --RTM
#   set(CMAKE_FIND_LIBRARY_SUFFIXES .a .lib)
# endif(PREFER_STATIC_LIBRARIES)



#
# Detect Fortran name mangling scheme for C/Fortran interface 
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


# ------------------------------------------------------------- #
#    Set preprocessor flags 
# ------------------------------------------------------------- #

#
# This defines were always on in older version
# of AMReX/CMake. Need more details on these???
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

if ( ENABLE_PARTICLES AND ( NOT ENABLE_DP_PARTICLES ) ) 
   add_define ( BL_SINGLE_PRECISION_PARTICLES )
endif ()

if (ENABLE_FORTRAN_MPI AND ENABLE_MPI)
   add_define ( BL_USE_FORTRAN_MPI=1 )
endif ()

add_define (MG_USE_FBOXIB=1 ENABLE_MG_FBOXLIB)
add_define (BL_USE_F_BASELIB=1  ENABLE_FBASELIB)
add_define (BL_USE_MPI ENABLE_MPI)
add_define (BL_USE_OMP ENABLE_OMP)



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
   list (APPEND AMREX_EXTRA_Fortran_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
   list (APPEND AMREX_EXTRA_C_FLAGS "${MPI_C_COMPILE_FLAGS}")
   list (APPEND AMREX_EXTRA_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
   # Libraries
   list (APPEND AMREX_EXTRA_Fortran_LIBRARIES "${MPI_Fortran_LIBRARIES}")
   list (APPEND AMREX_EXTRA_C_LIBRARIES "${MPI_C_LIBRARIES}")
   list (APPEND AMREX_EXTRA_CXX_LIBRARIES "${MPI_CXX_LIBRARIES}")
   # Link flags
   list (APPEND AMREX_EXTRA_Fortran_LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
   list (APPEND AMREX_EXTRA_C_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
   list (APPEND AMREX_EXTRA_CXX_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif ()

if (ENABLE_OMP)
   find_package (OpenMP REQUIRED)
   # Compile flags
   list (APPEND AMREX_EXTRA_Fortran_FLAGS "${OpenMP_Fortran_FLAGS}")
   list (APPEND AMREX_EXTRA_C_FLAGS "${OpenMP_C_FLAGS}")
   list (APPEND AMREX_EXTRA_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
endif()


# ------------------------------------------------------------- #
#    Setup compiler flags 
# ------------------------------------------------------------- #
if ( AMREX_FFLAGS_OVERRIDE )
   set ( AMREX_Fortran_FLAGS ${AMREX_FFLAGS_OVERRIDE} )
else ()
   add_compiler_flags ( AMREX_${FC_ID}_FFLAGS_${AMREX_BUILD_TYPE}
      AMREX_Fortran_FLAGS )
endif ()

if ( AMREX_CXXLAGS_OVERRIDE )
   set ( AMREX_CXX_FLAGS ${AMREX_CXXFLAGS_OVERRIDE} )
else ()
   add_compiler_flags ( AMREX_${CXX_ID}_CXXFLAGS_${AMREX_BUILD_TYPE}
      AMREX_CXX_FLAGS )
endif ()

# Add extra flags
add_compiler_flags ( AMREX_EXTRA_Fortran_FLAGS AMREX_Fortran_FLAGS )
add_compiler_flags ( AMREX_EXTRA_CXX_FLAGS AMREX_CXX_FLAGS )

# Add required flags
add_compiler_flags ( AMREX_${FC_ID}_FFLAGS_REQUIRED AMREX_Fortran_FLAGS )
add_compiler_flags ( AMREX_${CXX_ID}_CXXLAGS_REQUIRED AMREX_CXX_FLAGS )

# Set CMake compiler flags
set ( CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE} "${AMREX_Fortran_FLAGS}" ) 
set ( CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE} "${AMREX_CXX_FLAGS}" )
#
# Config summary
#
message( STATUS "AMReX configuration summary: ")
message( STATUS "   Build type            = ${CMAKE_BUILD_TYPE}")
message( STATUS "   Preprocessor flags    = ${AMREX_DEFINES}")
message( STATUS "   C++ compiler          = ${CMAKE_CXX_COMPILER}")
message( STATUS "   Fortran compiler      = ${CMAKE_Fortran_COMPILER}")
message( STATUS "   C++ flags             = ${CMAKE_CXX_FLAGS_${AMREX_BUILD_TYPE}}")
message( STATUS "   Fortran flags         = ${CMAKE_Fortran_FLAGS_${AMREX_BUILD_TYPE}}")
message( STATUS "   C++ include paths     = ${AMREX_EXTRA_CXX_INCLUDE_PATH}") 
message( STATUS "   Fortran include paths = ${AMREX_EXTRA_Fortran_INCLUDE_PATH}")
message( STATUS "   C++ external libs     = ${AMREX_EXTRA_CXX_LIBRARIES}") 
message( STATUS "   Fortran external libs = ${AMREX_EXTRA_Fortran_LIBRARIES}")
message( STATUS "   C++ link flags        = ${AMREX_EXTRA_CXX_LINK_FLAGS}") 
message( STATUS "   Fortran link flags    = ${AMREX_EXTRA_LINK_FLAGS}")
