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
# so that other included files can check if this file has been
# processed already
set (__AMREX_OPTIONS__ "")

#
# Define a macro to check the value of the inputs integer options 
# 
macro (print_option var)
   message ( STATUS "   ${var} = ${${var}}")
endmacro ()

#
# Check if CMAKE_BUILD_TYPE is given. If not, use default
#
if ( NOT CMAKE_BUILD_TYPE )
   message(STATUS "Setting build type to Release as none was specified.")
   set( CMAKE_BUILD_TYPE Release )
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()

#
# Populate the cache and check the value of the user-definable options 
#
message (STATUS "Configuring AMReX with the following options: ")


#
# This is the option to enable/disable xSDK mode
#
option( USE_XSDK_DEFAULTS "Enable xSDK mode"  OFF )
print_option ( USE_XSDK_DEFAULTS )


#
# Option to control the type of library build: static vs shared
#
if ( USE_XSDK_DEFAULTS )
   option ( BUILD_SHARED_LIBS "Build AMReX shared library" ON )
else ()
   option ( BUILD_SHARED_LIBS "Build AMReX shared library" OFF )
endif ()
print_option ( BUILD_SHARED_LIBS ) 

#
# Print out info on install path
# 
print_option ( CMAKE_INSTALL_PREFIX )


set (DIM 3 CACHE STRING "Dimension of AMReX build")
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

option( ENABLE_CUDA  "Enable CUDA" OFF)
print_option( ENABLE_CUDA )


if ( USE_XSDK_DEFAULTS )
   set ( XSDK_PRECISION "DOUBLE" CACHE STRING "Precision:<SINGLE,DOUBLE>" )
   if ( "${XSDK_PRECISION}" STREQUAL "DOUBLE" )
      set ( ENABLE_DP ON )
   elseif ( "${XSDK_PRECISION}" STREQUAL "SINGLE" )
      set ( ENABLE_DP OFF )
   else ()
      message ( WARNING "Unsupported precision ${XSDK_PRECISION}: defaulting to DOUBLE" )
      set ( XSDK_PRECISION "DOUBLE" CACHE STRING "Precision:<SINGLE,DOUBLE>" )
      set ( ENABLE_DP ON )
   endif ()
   print_option ( XSDK_PRECISION )
else ()
   option ( ENABLE_DP "Enable double precision" ON )
   print_option ( ENABLE_DP )
endif ()


#
# AMReX components selection
# 
option ( ENABLE_EB "Build EB Code" OFF )
print_option (ENABLE_EB)

if ( USE_XSDK_DEFAULTS )
   option ( XSDK_ENABLE_Fortran "Build Fortran API" OFF )
   print_option (XSDK_ENABLE_Fortran)
   set ( ENABLE_FORTRAN_INTERFACES ${XSDK_ENABLE_Fortran} )
   print_option (ENABLE_FORTRAN_INTERFACES)
else() 
   option ( ENABLE_FORTRAN_INTERFACES "Build Fortran API" ON )
   print_option (ENABLE_FORTRAN_INTERFACES)
endif ()

option ( ENABLE_LINEAR_SOLVERS  "Build AMReX Linear solvers" ON )
print_option ( ENABLE_LINEAR_SOLVERS )

option ( ENABLE_AMRDATA "Build data services" OFF)
print_option ( ENABLE_AMRDATA )

option ( ENABLE_PARTICLES "Build particle classes" OFF)
print_option ( ENABLE_PARTICLES )

if ( ENABLE_PARTICLES )
   if ( USE_XSDK_DEFAULTS )
      if ( "${XSDK_PRECISION}" STREQUAL "DOUBLE" )
	 set ( ENABLE_DP_PARTICLES ON )
      elseif ( "${XSDK_PRECISION}" STREQUAL "SINGLE" )
	 set ( ENABLE_DP_PARTICLES OFF )
      endif ()
   else ()
      option ( ENABLE_DP_PARTICLES "Enable double-precision particle data" ON )
      print_option ( ENABLE_DP_PARTICLES )
   endif ()
endif ()

option ( ENABLE_SENSEI_INSITU "Enable SENSEI in situ infrastructure" OFF )
print_option ( ENABLE_SENSEI_INSITU )

if (ENABLE_LINEAR_SOLVERS AND (DIM EQUAL 3) AND (NOT USE_XSDK_DEFAULTS) )
   option(ENABLE_3D_NODAL_MLMG "Enable 3D nodal MLMG" OFF)
   print_option(ENABLE_3D_NODAL_MLMG)
else ()
   set(ENABLE_3D_NODAL_MLMG OFF)
endif ()

#
# This options are paths to external libraries installation directories
#
if (USE_XSDK_DEFAULTS)
   set( ALGOIM_INSTALL_DIR "" CACHE PATH
      "Path to Algoim installation directory")
   set(  BLITZ_INSTALL_DIR "" CACHE PATH
      "Path to Blitz installation directory")
endif ()

#
# Compilation options
#  
option (ENABLE_FPE "Enable Floating Point Exceptions checks" OFF)
print_option ( ENABLE_FPE )

if ( "${CMAKE_BUILD_TYPE}" MATCHES "Debug" )
   option ( ENABLE_ASSERTIONS "Enable assertions" ON)
else ()
   option ( ENABLE_ASSERTIONS "Enable assertions" OFF)
endif ()

print_option ( ENABLE_ASSERTIONS )


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

