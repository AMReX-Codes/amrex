###############################################

# Here we define the default config options   #
# that can be overwritten by the user         #

###############################################

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
macro (check_option_value NAME VALUE RANGE_STARTS RANGE_ENDS )
   
   if ( ( ${VALUE} GREATER ${RANGE_ENDS} ) OR ( ${VALUE} LESS ${RANGE_STARTS} ) )
      message ( FATAL_ERROR "Variable ${NAME} has value ${VALUE}. \
Allowed range is [${RANGE_STARTS}:${RANGE_ENDS}]." )
   endif ()

   message ( STATUS "   ${NAME} = ${VALUE} (INT: ${RANGE_STARTS},${RANGE_ENDS})" )
endmacro ()


#
# Populate the cache and check the value of the user-definable options 
#
message (STATUS "Configuring AMReX with the following options: ")

if ( NOT CMAKE_BUILD_TYPE )
   # Default to Release if no other build type specified
   set ( CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug, Release."
      FORCE )
endif ()

set ( AMREX_BUILD_TYPE ${CMAKE_BUILD_TYPE} )

# Need to be uppercase so it can be used to refer to CMAKE variable names
string ( TOUPPER ${AMREX_BUILD_TYPE} AMREX_BUILD_TYPE) 


message ( STATUS "   CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE} (STRING:\
 Debug|Release|RelWithDebInfo|MinSizeRel)" )

message ( STATUS "   CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX} (STRING: <path to install dir>)" )


set (ENABLE_PIC 0 CACHE INT
   "Compile with position-independent code enabled")
check_option_value ( "ENABLE_PIC" ${ENABLE_PIC} 0 1 )

set (BL_SPACEDIM 3 CACHE INT "Dimension of AMReX build")
check_option_value ( "BL_SPACEDIM" ${BL_SPACEDIM} 2 3 )

set (ENABLE_MPI 1 CACHE INT "Enable build with MPI")
check_option_value ( "ENABLE_MPI" ${ENABLE_MPI} 0 1 )

set (ENABLE_OMP 0 CACHE INT "Enable build with OpenMP")
check_option_value ( "ENABLE_OMP" ${ENABLE_OMP} 0 1 )

set (ENABLE_DP 1 CACHE INT "Enable double precision build")
check_option_value ( "ENABLE_DP" ${ENABLE_DP} 0 1 )

set (ENABLE_PARTICLES 0 CACHE INT "Include Particles classes in AMReX build")
check_option_value ( "ENABLE_PARTICLES" ${ENABLE_PARTICLES} 0 1 )

set (ENABLE_DP_PARTICLES 1 CACHE INT "Enable double-precision for particles data") 
check_option_value ( "ENABLE_DP_PARTICLES" ${ENABLE_DP_PARTICLES} 0 1 )

set (ENABLE_PROFILING 0 CACHE INT "Include profiling information in AMReX build")
check_option_value ( "ENABLE_PROFILING" ${ENABLE_PROFILING} 0 1 )

set (ENABLE_TINY_PROFILING 0 CACHE INT "Include 'tiny'-profiling information in AMReX build")
check_option_value ( "ENABLE_TINY_PROFILING" ${ENABLE_TINY_PROFILING} 0 1 )

set (ENABLE_TRACE_PROFILING 0 CACHE INT  "Include trace-profiling information in AMReX build" )
check_option_value ( "ENABLE_TRACE_PROFILING" ${ENABLE_TRACE_PROFILING} 0 1 )

set (ENABLE_COMM_PROFILING 0 CACHE INT  "Include comm-profiling information in AMReX build" )
check_option_value ( "ENABLE_COMM_PROFILING" ${ENABLE_COMM_PROFILING} 0 1 )

set (ENABLE_BACKTRACE 0 CACHE INT "Include backtrace information in AMReX build")
check_option_value ( "ENABLE_BACKTRACE" ${ENABLE_BACKTRACE} 0 1 )

set (ENABLE_FPE 0 CACHE INT "Enable Floating Point Exceptions checks")
check_option_value ( "ENABLE_FPE" ${ENABLE_FPE} 0 1 )

set (ENABLE_ASSERTIONS 0 CACHE INT "Enable assertions")
check_option_value ( "ENABLE_ASSERTIONS" ${ENABLE_ASSERTIONS} 0 1 )

set (ENABLE_FORTRAN_INTERFACES 1 CACHE INT "Include Fortran interfaces in AMReX build")
check_option_value ( "ENABLE_FORTRAN_INTERFACES" ${ENABLE_FORTRAN_INTERFACES} 0 1 )

set (ENABLE_LINEAR_SOLVERS 1 CACHE INT "Include Linear solvers in AMReX build")
check_option_value ( "ENABLE_LINEAR_SOLVERS" ${ENABLE_LINEAR_SOLVERS} 0 1 )

set ( ENABLE_FBASELIB 1 CACHE INT  "Enable Fortran BaseLib" )
check_option_value ( "ENABLE_FBASELIB" ${ENABLE_FBASELIB} 0 1 )


set (AMREX_FFLAGS_OVERRIDES "" CACHE STRING "User-defined Fortran compiler flags" )

set (AMREX_CXXFLAGS_OVERRIDES "" CACHE STRING "User-defined C++ compiler flags" )

# After the options are set, define the following variable
# so that other included file can check if this file has been
# run already
set ( AMREX_OPTIONS_SET  "TRUE" )  

