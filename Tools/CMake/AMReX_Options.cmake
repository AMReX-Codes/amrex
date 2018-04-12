###############################################

# Here we define the default config options   #
# that can be overwritten by the user         #

###############################################

#
# Include module 
# 
include (CMakeDependentOption)

if (DEFINED __AMREX_OPTIONS__)
   return ()
endif ()

# Define the following variable
# so that other included file can check if this file has been
# run already
set (__AMREX_OPTIONS__ "")

#
# Check weather the AMReX_CMakeVariables.cmake
# has been loaded; abort if not
#
if ( NOT ( DEFINED __AMREX_CMAKEVARIABLES__ ) )
   message ( FATAL_ERROR "AMReX_CMakeVariables.cmake must be included\
before including AMReX_Options.cmake" )
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
option ( ENABLE_EB "Build EB Code" OFF )
print_option (ENABLE_EB)

option ( ENABLE_FORTRAN_INTERFACES "Build Fortran API" ON )
print_option (ENABLE_FORTRAN_INTERFACES)

option ( ENABLE_LINEAR_SOLVERS  "Build AMReX Linear solvers" ON )
print_option ( ENABLE_LINEAR_SOLVERS )

option ( ENABLE_FBASELIB "Build Fortran kernel (deprecated)" ON )
print_option ( ENABLE_FBASELIB )

option ( ENABLE_AMRDATA "Build data services" OFF)
print_option ( ENABLE_AMRDATA )

option ( ENABLE_PARTICLES "Build particle classes" OFF)
print_option ( ENABLE_PARTICLES )

if ( ENABLE_PARTICLES )
   option ( ENABLE_DP_PARTICLES "Enable double-precision particle data" ON )
   print_option ( ENABLE_DP_PARTICLES )
endif ()


#
# Compilation options
#  
option (ENABLE_FPE "Enable Floating Point Exceptions checks" OFF)
print_option ( ENABLE_FPE )

if (DEBUG)
   option ( ENABLE_ASSERTION "Enable assertions" ON)
else ()
   option ( ENABLE_ASSERTION "Enable assertions" OFF)
endif ()
print_option ( ENABLE_ASSERTION )

set (AMREX_FFLAGS_OVERRIDES "" CACHE STRING "User-defined Fortran compiler flags" )

set (AMREX_CXXFLAGS_OVERRIDES "" CACHE STRING "User-defined C++ compiler flags" )


#
# Profiling options
# 
option ( ENABLE_BASE_PROFILE "Enable basic profiling" OFF )
print_option ( ENABLE_BASE_PROFILE )

option ( ENABLE_TINY_PROFILE "NOT ENABLE_BASE_DEBUG" OFF)
print_option ( ENABLE_TINY_PROFILE ) 

option ( ENABLE_TRACE_PROFILE "Enable trace-profiling" OFF )
print_option ( ENABLE_TRACE_PROFILE )

option ( ENABLE_MEM_PROFILE   "Enable memory profiling" OFF )
print_option ( ENABLE_MEM_PROFILE )

option ( ENABLE_COMM_PROFILE  "Enable communicator-profiling" OFF )
print_option ( ENABLE_COMM_PROFILE )

option ( ENABLE_BACKTRACE "Enable backtracing" OFF)
print_option ( ENABLE_BACKTRACE )

option ( ENABLE_PROFPARSER "Enable profile parser" OFF)
print_option ( ENABLE_PROFPARSER )

option ( ENABLE_VTUNE "Enable VTune" OFF)
print_option ( ENABLE_VTUNE )

option ( ENABLE_CRAYPAT "Enable CrayPat" OFF)
print_option ( ENABLE_CRAYPAT )

option ( ENABLE_ALLINEA "Enable Allinea MAP" OFF)
print_option ( ENABLE_ALLINEA )



# Modify the profiling options if needed ( because of dependencies between
# the options )
if (ENABLE_PROFPARSER)
   set (ENABLE_BASE_PROFILE ON)
   set (ENABLE_TRACE_PROFILE ON)
   set (ENABLE_AMRDATA ON)
endif ()

if (ENABLE_TRACE_PROFILE OR ENABLE_COMM_PROFILE)
   set (ENABLE_BASE_PROFILE ON)
endif ()

if (ENABLE_BASE_PROFILE)
   set (ENABLE_TINY_PROFILE OFF)
endif()

# Check consistency of options
if (  (ENABLE_VTUNE   AND ENABLE_CRAYPAT) OR
      (ENABLE_VTUNE   AND ENABLE_ALLINEA) OR
      (ENABLE_ALLINEA AND ENABLE_CRAYPAT) )
   message ( FATAL_ERROR
      "Please select ONLY one of ENABLE_VTUNE, ENABLE_CRAYPAT, ENABLE_ALLINEA and re-configure" )
endif ()

if ( (ENABLE_VTUNE OR ENABLE_CRAYPAT OR ENABLE_ALLINEA) AND
      (ENABLE_BASE_PROFILE OR ENABLE_TINY_PROFILE) )
   message (WARNING "This configuration should only be used to profile BL_PROFILE!")
endif()

# After the options are set, define the following variable
# so that other included file can check if this file has been
# run already
set ( AMREX_OPTIONS_SET  "TRUE" )  

