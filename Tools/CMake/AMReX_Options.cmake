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
# This is the option to enable/disable xSDK mode
#
option ( USE_XSDK_DEFAULTS "Enable xSDK mode"  OFF )
print_option ( USE_XSDK_DEFAULTS )

#
# Option to control the build type.
# If xSDK mode is enabled, the build type is determined by
# CMAKE_BUILD_TYPE and DEBUG is a regular variable.
# If xSDK mode is disabled, the build type is determined by DEBUG
# and CMAKE_BUILD_TYPE is set accordingly.
# 
if ( USE_XSDK_DEFAULTS )
   set ( DEBUG OFF )
   if ( ( "${CMAKE_BUILD_TYPE}" MATCHES "Debug" ) OR
	( NOT CMAKE_BUILD_TYPE ) )
      set ( DEBUG ON )
   endif ()
else ()
   option ( DEBUG "Build in debug mode" OFF )
endif ()

if (DEBUG)
   set (CMAKE_BUILD_TYPE "Debug")
else ()
   set (CMAKE_BUILD_TYPE "Release")
endif ()

if ( USE_XSDK_DEFAULTS )
   print_option (CMAKE_BUILD_TYPE)
else ()
   print_option (DEBUG)
endif ()

string ( TOUPPER ${CMAKE_BUILD_TYPE} AMREX_BUILD_TYPE ) 

#
# Print out info on install path
# 
print_option ( CMAKE_INSTALL_PREFIX )


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

set ( TP_PROFILE "None" CACHE STRING "Third-party profiling options:<CRAYPAT,FORGE,VTUNE>")
print_option ( TP_PROFILE )


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

# Check profile options
if ( NOT ( ${CMAKE_C_COMPILER_ID} MATCHES "Intel" ) AND
      ( ${TP_PROFILE} MATCHES "VTUNE") )
   message ( WARNING "VTUNE cannot be used with ${CMAKE_C_COMPILER_ID} compiler: ignoring TP_PROFILE" )
   set ( TP_PROFILE "")
endif ()

if (  ( ( ${TP_PROFILE} MATCHES "CRAYPAT" ) OR
        ( ${TP_PROFILE} MATCHES "FORGE"   ) OR
        ( ${TP_PROFILE} MATCHES "VTUNE"   )   ) AND
     (ENABLE_BASE_PROFILE OR ENABLE_TINY_PROFILE) )
   message (WARNING "This configuration should only be used to profile BL_PROFILE!")
endif()

