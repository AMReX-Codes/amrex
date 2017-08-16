###############################################

# Here we define the default config options   #
# that can be overwritten by the user         #

###############################################

#
# Include module 
# 
include (CMakeDependentOption)

#
# Check weather the AMReX_CMakeVariables.cmake
# has been loaded; abort if not
#
if ( NOT AMREX_VARIABLES_LOADED )
   message ( FATAL_ERROR "AMReX_Options.cmake must be included\
after including AMReX_CMakeVariables.cmake" )
endif ()

#
# Define a macro to check the value of the inputs integer options 
# 
macro (print_option var)
   message ( STATUS "   ${var} = ${${var}}")
endmacro ()

#
# Populate the cache and check the value of the user-definable options 
#
message (STATUS "Configuring AMReX with the following options: ")

#
# General options
# 
option ( DEBUG "Build in debug mode" OFF )
print_option (DEBUG)

if (DEBUG)
   set (CMAKE_BUILD_TYPE "Debug")
else ()
   set (CMAKE_BUILD_TYPE "Release")
endif ()

string ( TOUPPER ${CMAKE_BUILD_TYPE} AMREX_BUILD_TYPE ) 


message ( STATUS "   CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )


set (DIM 3 CACHE INT "Dimension of AMReX build")
if ( (${DIM} GREATER 3) OR (${DIM} LESS 1) )
   message ( FATAL_ERROR "DIM must be either 1, 2 or 3.")
endif ()
set (BL_SPACEDIM ${DIM})
print_option ( DIM )


option ( ENABLE_PIC "Build position-independent code" OFF)
print_option ( ENABLE_PIC )

option ( ENABLE_MPI  "Enable MPI"  ON)
print_option ( ENABLE_MPI )

option ( ENABLE_OMP  "Enable OpenMP" OFF)
print_option ( ENABLE_OMP )

option (ENABLE_DP "Enable double precision" ON)
print_option ( ENABLE_DP )


#
# AMReX components selection
# 
option ( ENABLE_FORTRAN_INTERFACES "Build Fortran API" ON )
print_option (ENABLE_FORTRAN_INTERFACES)

option ( ENABLE_LINEAR_SOLVERS  "Build AMReX Linear solvers" ON )
print_option ( ENABLE_LINEAR_SOLVERS )

option ( ENABLE_FBASELIB "Build Fortran kernel (deprecated)" ON )
print_option ( ENABLE_FBASELIB )

option ( ENABLE_PARTICLES "Build particle classes" OFF)
print_option ( ENABLE_PARTICLES )

option ( ENABLE_AMRDATA "Build data services" OFF)
print_option ( ENABLE_AMRDATA )


# option ( ENABLE_DP_PARTICLES "Enable double-precision for particles data" ON) 
cmake_dependent_option ( ENABLE_DP_PARTICLES "Enable double-precision particle data" ON
                          "ENABLE_PARTICLES" OFF)
print_option ( ENABLE_DP_PARTICLES )


#
# Compilation options
#  
option (ENABLE_FPE "Enable Floating Point Exceptions checks" OFF)
print_option ( ENABLE_FPE )

cmake_dependent_option ( ENABLE_ASSERTION "Enable assertions" ON "DEBUG" OFF)
print_option ( ENABLE_ASSERTION )

set (AMREX_FFLAGS_OVERRIDES "" CACHE STRING "User-defined Fortran compiler flags" )

set (AMREX_CXXFLAGS_OVERRIDES "" CACHE STRING "User-defined C++ compiler flags" )

#
# Profiling options
# 
option ( ENABLE_BASE_PROFILE "Enable basic profiling" OFF )
print_option ( ENABLE_BASE_PROFILE )

option ( ENABLE_TINY_PROFILE "Enable 'tiny'-profiling" OFF )
print_option ( ENABLE_TINY_PROFILE ) 

option ( ENABLE_TRACE_PROFILE "Enable trace-profiling" OFF )
print_option ( ENABLE_TRACE_PROFILE )

option ( ENABLE_MEM_PROFILE   "Enable memory profiling" OFF )
print_option ( ENABLE_MEM_PROFILE )

option ( ENABLE_COMM_PROFILE  "Enable communicator-profiling" OFF )
print_option ( ENABLE_COMM_PROFILING )

option ( ENABLE_BACKTRACE "Enable backtracing" OFF)
print_option ( ENABLE_BACKTRACE )

option ( ENABLE_PROFPARSER "Enable profile parser" OFF)
print_option ( ENABLE_PROFPARSER )



#
# If any profiling options is activated, set AMREX_PROFILING to 1
# 
set ( AMREX_ANY_PROFILING 0 )
if ( ENABLE_PROFILING OR ENABLE_COMM_PROFILING OR ENABLE_TINY_PROFILING OR
      ENABLE_TRACE_PROFILING OR ENABLE_MEM_PROFILING)
   set ( AMREX_ANY_PROFILING 1 )
endif ()



# After the options are set, define the following variable
# so that other included file can check if this file has been
# run already
set ( AMREX_OPTIONS_SET  "TRUE" )  

